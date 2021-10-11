//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate electrical resistivity of a network read from an Abaqus output file (.odb)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "App_Network_From_Abaqus.h"

//Generate 3D nanoparticle network, turn it into a resitor network and find its electrical conductivity
int App_Network_From_Abaqus::Nanoparticle_resistor_network_from_odb(Input* Init)const
{
    //Time markers for total simulation
    time_t ct0, ct1;

    //Variables for CNTs
    //CNTs points
    vector<Point_3D> points_cnt;
    //CNTs radii
    vector<double> radii;
    //CNT structure, each cnts_structure[i] referes to the points in CNT_i
    vector<vector<long int> > structure_cnt;

    //Variables for GNPs
    //GNPs
    vector<GNP> gnps;
    //GNP points (only those needed are stored)
    vector<Point_3D> points_gnp;
    //GNP structure, each structure_gnp[i] referes to the points in GNP_i
    vector<vector<long int> > structure_gnp;

    //Shell vectors (used to remove nanoparticles when reducing observation window size)
    vector<vector<int> > shells_cnts;

    //----------------------------------------------------------------------
    //Read CNTs and GNPs from csv file (generated by Generate_Network) and 
    //fill the CNT and GNP variables
    //Also read the sample geometry
    hout << "Generating nanoparticle network from file......" << endl;
    ct0 = time(NULL);
    Read_Network* Reader = new Read_Network;
    if (!Reader->Generate_nanoparticle_network_from_file(Init->simu_para, Init->vis_flags, Init->geom_sample, points_cnt, radii, structure_cnt, gnps)) {
        hout << "Error in Generate_nanoparticle_network_from_file." << endl;
        return 0;
    }
    delete Reader;
    ct1 = time(NULL);
    hout << "Network generation from file time: " << (int)(ct1 - ct0) << " secs." << endl;

    //Generate vector with GNPs partially outside the sample
    vector<int> gnps_outside;
    if (!Get_gnps_partially_outside_sample(Init->geom_sample, gnps, gnps_outside))
    {
        hout << "Error in Get_gnps_partially_outside_sample." << endl;
        return 0;
    }

    //----------------------------------------------------------------------
    //Vector for GNP shells
    vector<Shell> shell_gnps(gnps.size());
    ct0 = time(NULL);
    Shells* SH = new Shells;
    if (!SH->Generate_shells(Init->simu_para.particle_type, Init->geom_sample, points_cnt, gnps, shells_cnts, shell_gnps)) {
        hout << "Error when generating shells" << endl;
        return 0;
    }
    delete SH;
    ct1 = time(NULL);
    hout << "Generate shells and structure time: " << (int)(ct1 - ct0) << " secs." << endl;
    /*for (size_t i = 0; i < shells_cnts.size(); i++)
    {
        hout << "Shell " << i << ": " << shells_cnts[i].size() << endl;
    }*/

    //----------------------------------------------------------------------
    //Determine the local networks in cutoff windows
    Cutoff_Wins* Cutwins = new Cutoff_Wins;
    //From this function I get the internal variables cnts_inside and boundary_cnt
    ct0 = time(NULL);
    //hout<<"Extract_observation_window"<<endl;
    //For the case of reading data from an Abaqus database, window is always 0
    if (!Cutwins->Extract_observation_window(0, Init->simu_para.particle_type, Init->geom_sample, Init->geom_sample.sample, Init->nanotube_geo, gnps, structure_cnt, radii, points_cnt, shells_cnts, shell_gnps, structure_gnp, points_gnp)) {
        hout << "Error when extracting observation window " << endl;
        return 0;
    }
    ct1 = time(NULL);
    hout << "Extract observation window time: " << (int)(ct1 - ct0) << " secs." << endl;

    //Initialize Abaqus C++ API
    odb_initializeAPI();

    //Open Abaqus database using the name/path intidated in the input file
    //Use a C string, since that seems to be equivalent to (or able to be cast as) an odb_String
    odb_Odb& odb = openOdb(Init->simu_para.odb_file.c_str());

    //Access the root assembly
    odb_Assembly& root_assy = odb.rootAssembly();

    //Get all frames from the steps and save them in a (pointer) variable
    odb_SequenceFrame& allFramesInStep = odb.steps()[Init->simu_para.step_name.c_str()].frames();
    //Get the number of frames in the database
    int n_frames = allFramesInStep.size();
    hout << endl << "There are " << n_frames << " frames in the Abaqus database" << endl;

    //Iterate over the number of frames
    for (int i = 0; i < n_frames; i++)
    {
        hout << "============================================================================" << endl;
        hout << "============================================================================" << endl;
        hout << "Frame " << i << endl;
        time_t it0, it1;
        it0 = time(NULL);

        //----------------------------------------------------------------------
        //Apply displacement to CNTs, GNPs and sample, except for frame 0 (which has no displacement)
        if (i)
        {
            //Read Abaqus database and apply displacements
            ct0 = time(NULL);
            //hout<<"Apply_displacements_from_Abaqus"<<endl;
            if (!Apply_displacements_from_Abaqus(Init->simu_para.particle_type, (int)structure_cnt.size(), structure_cnt, gnps_outside, root_assy, allFramesInStep[i-1], allFramesInStep[i], Init->geom_sample, points_cnt, gnps))
            {
                hout << "Error when applying displacement to nanoparticles." << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Apply displacements for current frame time: " << (int)(ct1 - ct0) << " secs." << endl;
        }

        //Window geometry is the same as that of the sample
        cuboid window_geo = Init->geom_sample.sample;
        //hout<<"window_geo = "<<window_geo.str()<<endl;

        //Export the window geometry if needed
        if (Init->vis_flags.window_domain) {
            string str = "window_" + to_string(i) + ".vtk";
            VTK_Export VTK_E;
            VTK_E.Export_cuboid(window_geo, str);
        }

        //----------------------------------------------------------------------
        //Determine the local networks inside the cutoff windows
        Contact_grid* Contacts = new Contact_grid;
        ct0 = time(NULL);
        //For the case of reading data from an Abaqus database, window is always 0
        if (!Contacts->Generate_contact_grid(0, Init->simu_para.particle_type, Init->geom_sample, window_geo, Cutwins->cnts_inside, points_cnt, structure_cnt, Cutwins->gnps_inside, gnps)) {
            hout << "Error when generating contact grid" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Generate contact grid time: " << (int)(ct1 - ct0) << " secs." << endl;

        //----------------------------------------------------------------------
        //Hoshen-Kopelman algorithm
        Hoshen_Kopelman* HoKo = new Hoshen_Kopelman;
        ct0 = time(NULL);
        if (!HoKo->Determine_clusters_and_percolation(i, Init->simu_para, Init->cutoff_dist, Init->vis_flags, Cutwins->cnts_inside, Contacts->sectioned_domain_cnts, structure_cnt, points_cnt, radii, Cutwins->boundary_cnt, Cutwins->gnps_inside, Contacts->sectioned_domain_gnps, gnps, Cutwins->boundary_gnp, structure_gnp, points_gnp)) {
            hout << "Error when finding clusters and determining percolation" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Find clusters and determine percolation: " << (int)(ct1 - ct0) << " secs." << endl;

        //Contacts are not needed anymore, so delete the object
        delete Contacts;

        //Loop over the different clusters so that the direct electrifying algorithm is aplied on each cluster
        Electrical_analysis* EA = new Electrical_analysis;
        if (!EA->Perform_analysis_on_clusters(i, window_geo, Init->simu_para, Init->electric_para, Init->cutoff_dist, Init->vis_flags, Init->out_flags, HoKo, Cutwins, structure_cnt, points_cnt, radii, points_gnp, structure_gnp, gnps)) {
            hout << "Error when performing electrical analysis" << endl;
            return 0;
        }

        //Delete objects to free memory
        delete EA;
        delete HoKo;

        it1 = time(NULL);
        hout << "Frame " << i << " time: " << (int)(it1 - it0) << " secs." << endl;
    }

    //Delete objects to free memory
    delete Cutwins;

    //Close Abaqus database
    odb.close();

    //Finalize usage of Abaqus C++ API
    odb_finalizeAPI();

    return 1;
}
//This function finds the GNPs that are partially outside the sample and stores the 
//indices of those GNPs in a vector
int App_Network_From_Abaqus::Get_gnps_partially_outside_sample(const Geom_sample& geom_sample, const vector<GNP>& gnps, vector<int>& gnps_outside)const
{
    //Iterate over the GNPs
    for (int i = 0; i < (int)gnps.size(); i++)
    {
        //Check if all eight vertices are inside the sample
        for (int j = 0; j < 8; j++)
        {
            //Check if vertex j of GNP i is outside the sample
            if (gnps[i].vertices[j].is_outside_cuboid(geom_sample.sample))
            {
                //Vertex j is outsied the sample, so add GNP i to the vector
                //of GNPs partially outside the sample
                gnps_outside.push_back(i);

                //Break the loop over j as there is no need to check the rest of vertices
                //to determine that the GNP is partially outside the sample
                break;
            }
        }
    }

    return 1;
}
//This functions adds the displacements to the CNTs, GNPs and sample
int App_Network_From_Abaqus::Apply_displacements_from_Abaqus(const string& particle_type, const int& n_cnts, const vector<vector<long int> >& structure, const vector<int>& gnps_outside, odb_Assembly& root_assy, odb_Frame& previous_frame, odb_Frame& current_frame, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<GNP>& gnps)const
{
    //Access displacement field ("U") in the current frame
    //hout << "current_fieldU" << endl;
    odb_FieldOutput& current_fieldU = current_frame.fieldOutputs()["U"];
    //Access displacement field ("U") in the current frame
    //hout << "previous_fieldU" << endl;
    odb_FieldOutput& previous_fieldU = previous_frame.fieldOutputs()["U"];

    //Apply displacements to sample
    //hout << "Apply_displacements_to_sample" << endl;
    if (!Apply_displacements_to_sample(root_assy, previous_fieldU, current_fieldU, geom_sample))
    {
        hout << "Error in Apply_displacements_from_Abaqus when calling Apply_displacements_to_sample" << endl;
        return 0;
    }

    //Check the type of nanoparticle 
    if (particle_type == "CNT_wires" || particle_type == "CNT_deposit" || particle_type == "GNP_CNT_mix") {

        //Apply displacements to CNT points
        //hout << "Apply_displacements_to_cnts" << endl;
        if (!Apply_displacements_to_cnts(structure, root_assy, previous_fieldU, current_fieldU, n_cnts, points_cnt))
        {
            hout << "Error in Apply_displacements_from_Abaqus when calling Apply_displacements_to_cnts" << endl;
            return 0;
        }
    }
    else if (particle_type == "GNP_cuboids" || particle_type == "GNP_CNT_mix") {

        //Apply displacements to GNP vertices

    }
    else if (particle_type == "Hybrid_particles") {
        hout << "Hybrid particles not yet implemented" << endl;
        return 0;
    }
    else {
        hout << "Error in Apply_displacements_from_Abaqus: the type of particles should be one of the following: CNT_wires, CNT_deposit, GNP_cuboids, or GNP_CNT_mix. Input value was: " << particle_type << endl;
        return 0;
    }

    return 1;
}
//This function gets the displacements of two of the sample nodes, applies this displacemets
//to the corresponding vertices and calculates the geometry of the deformed sample
int App_Network_From_Abaqus::Apply_displacements_to_sample(odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, Geom_sample& geom_sample)const
{
    //Name of the set for the lower left corner (it is hard coded in the python scritp too)
    string set0 = "MATRIX0";

    //Name of the set for the opposite corner (it is hard coded in the python scritp too)
    string set1 = "MATRIX1";

    //Access set0 from root assembly
    odb_Set& matrix0 = root_assy.nodeSets()[set0.c_str()];

    //Get the displacement objects of the set
    odb_FieldOutput matrix0_disp = current_fieldU.getSubset(matrix0);
    odb_FieldOutput matrix0_disp_prev = previous_fieldU.getSubset(matrix0);

    //Get the sequence of values of the displacement object for matrix0
    const odb_SequenceFieldValue& vals0 = matrix0_disp.values();
    const odb_SequenceFieldValue& vals0_prev = matrix0_disp_prev.values();

    //Check the size of vals is 1 (since it contains one node)
    if (vals0.size() != 1)
    {
        hout << "Error in Apply_displacements_to_sample. The number of node set " << set0 << " is not 1. Size is " << vals0.size() << endl;
        return 0;
    }

    //Get the acutal values
    const odb_FieldValue val0 = vals0[0];
    const odb_FieldValue val0_prev = vals0_prev[0];
    //Output node label
    //hout << "  Node: " << val0.nodeLabel() << endl;
    //Get the current data of the displacements
    int numComp = 0; //This integer is needed to call data() in the line below
    const float* const data0 = val0.data(numComp);
    //Get the previous data of the displacements
    int numComp_prev = 0; //This integer is needed to call data() in the line below
    const float* const data0_prev = val0_prev.data(numComp_prev);
    //hout << "vals0.size=" << vals0.size() << " data0.size=numComp=" << numComp << endl;

    //Update lower left corner of sample
    // 
    //Note data[i] is the displacement for frame i, no the displacement
    //with respect to the previous frame
    //Thus, calculate the increment of displacement with respect to the previous frame
    //Otherwise I would need the initial geometry of the sample for each frame
    geom_sample.sample.poi_min = geom_sample.sample.poi_min + Point_3D((double)data0[0] - (double)data0_prev[0], (double)data0[1] - (double)data0_prev[1], (double)data0[2] - (double)data0_prev[2]);
    //hout << "disp0=" << data0[0] << " " << data0[1] << " " << data0[2] << endl;

    //Access set1 from root assembly
    odb_Set& matrix1 = root_assy.nodeSets()[set1.c_str()];

    //Get the displacement objects of the set
    odb_FieldOutput matrix1_disp = current_fieldU.getSubset(matrix1);
    odb_FieldOutput matrix1_disp_prev = previous_fieldU.getSubset(matrix1);

    //Get the sequence of values of the displacement object for matrix1
    const odb_SequenceFieldValue& vals1 = matrix1_disp.values();
    const odb_SequenceFieldValue& vals1_prev = matrix1_disp_prev.values();

    //Check the size of vals is 1 (since it contains one node)
    if (vals1.size() != 1)
    {
        hout << "Error in Apply_displacements_to_sample. The number of node set " << set1 << " is not 1. Size is " << vals1.size() << endl;
        return 0;
    }

    //Get the acutal values
    const odb_FieldValue val1 = vals1[0];
    const odb_FieldValue val1_prev = vals1_prev[0];
    //Output node label
    //cout << "  Node: " << val1.nodeLabel() << endl;
    //Get the data of the displacements
    numComp = 0; //This integer is needed to call data() in the line below
    const float* const data1 = val1.data(numComp);
    numComp_prev = 0; //This integer is needed to call data() in the line below
    const float* const data1_prev = val1_prev.data(numComp_prev);
    //hout << "disp1=" << data1[0] << " " << data1[1] << " " << data1[2] << endl;

    //Update the maximum coordinates of the sample, which are the maximum coordinates of the sample
    // 
    //Note data[i] is the displacement for frame i, no the displacement
    //with respect to the previous frame
    //Thus, calculate the increment of displacement with respect to the previous frame
    //Otherwise I would need the initial geometry of the sample for each frame
    geom_sample.sample.max_x = geom_sample.sample.max_x + (double)data1[0] - (double)data1_prev[0];
    geom_sample.sample.max_y = geom_sample.sample.max_y + (double)data1[1] - (double)data1_prev[1];
    geom_sample.sample.max_z = geom_sample.sample.max_z + (double)data1[2] - (double)data1_prev[2];

    //Update the dimensions of the sample along each direction
    geom_sample.sample.len_x = geom_sample.sample.max_x - geom_sample.sample.poi_min.x;
    geom_sample.sample.wid_y = geom_sample.sample.max_y - geom_sample.sample.poi_min.y;
    geom_sample.sample.hei_z = geom_sample.sample.max_z - geom_sample.sample.poi_min.z;

    //Update the sample's volume
    geom_sample.volume = geom_sample.sample.len_x * geom_sample.sample.wid_y * geom_sample.sample.hei_z;

    /* /Access set1 from root assembly
    odb_Set& matrix2 = root_assy.nodeSets()["MATRIX2"];
    //Get the displacement object of the set
    odb_FieldOutput matrix2_disp = current_fieldU.getSubset(matrix2);
    //Get the sequence of values of the displacement object for matrix1
    const odb_SequenceFieldValue& vals2 = matrix2_disp.values();
    //Get the acutal values
    const odb_FieldValue val2 = vals2[0];
    //Get the data of the displacements
    numComp = 0; //This integer is needed to call data() in the line below
    const float* const data2 = val2.data(numComp);
    hout << "disp2=" << data2[0] << " " << data2[1] << " " << data2[2] << endl;*/

    return 1;
}
//This function applies displacements to CNT points
int App_Network_From_Abaqus::Apply_displacements_to_cnts(const vector<vector<long int> >& structure, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, const int& n_cnts, vector<Point_3D>& points_cnt)const
{

    //Iterate over the CNTs in the sample
    for (size_t i = 0; i < structure.size(); i++)
    {
        //Get the name of the set
        //Numbering of sets starts at 1
        string set_name = Get_cnt_set_name((int)i+1);

        //Access set from root assembly
        odb_Set& cnt_set = root_assy.nodeSets()[set_name.c_str()];

        //Get the displacement objects of the set
        odb_FieldOutput cnt_disp = current_fieldU.getSubset(cnt_set);
        odb_FieldOutput cnt_disp_prev = previous_fieldU.getSubset(cnt_set);

        //Get the values of the displacement object
        const odb_SequenceFieldValue& vals = cnt_disp.values();
        const odb_SequenceFieldValue& vals_prev = cnt_disp_prev.values();
        //Output size of values
        //cout << "values.size=" << vals.size() << endl;

        //Check there is the same number of points in the CNT and in the set
        if (structure[i].size() != vals.size())
        {
            hout << "Error in Apply_displacements_to_cnts. The number of points in CNT " << i << " (" << structure[i].size() << " points) is different from the number of points in set " << set_name << " ("<< vals.size() << " points)." << endl;
            return 0;
        }

        /* /Get the nodes in the set
        const odb_SequenceNode& nodeList = root_assy.nodeSets()[set_name.c_str()].nodes();
        //Compare the coordinates of the first node in the set with the
        //first and last points in the CNT
        const float* const coord = nodeList[0].coordinates();
        hout << "Node[" << nodeList[0].label() <<"]="<< coord[0]<<" "<<coord[1]<<" "<<coord[2] << endl;
        hout << "P[first]=" << points_cnt[structure[i][0]].str() << endl;
        hout<<"P[last]="<< points_cnt[structure[i].back()].str() << endl;
        if (abs(coord[0] - points_cnt[structure[i].back()].x) < 1e-6 && abs(coord[1] - points_cnt[structure[i].back()].y) < 1e-6 && abs(coord[2] - points_cnt[structure[i].back()].z) < 1e-6){
            hout << "FIRST NODE IS THE LAST POINT" << endl;
        }
        else {
            hout << "DIFS=" << abs(coord[0] - points_cnt[structure[i].back()].x) << " " << abs(coord[1] - points_cnt[structure[i].back()].y) << " " << abs(coord[2] - points_cnt[structure[i].back()].z) << endl;
        }
        hout << endl;*/

        //Iterate over the values, i.e., the points in the current CNT
        //Number used in the loop
        int v_size1 = vals.size() - 1;
        for (int j = 0; j < vals.size(); j++)
        {
            //Note that nodes are actually in reverse order given the way
            //CNTs are generated (and meshed) in Abaqus
            //Thus, I need to go in reverse order when scanning the structure vcetor
            //With this definition of idx, I can start at the last element of 
            //structure[i] and stop at its first element
            int idx = v_size1 - j;

            //Get current point number
            long int P = structure[i][idx];
            
            //Get current values
            const odb_FieldValue val = vals[j];
            //Get previous values
            const odb_FieldValue val_prev = vals_prev[j];
            //Output node label
            //cout << "  Node: " << val.nodeLabel() << endl;
            //Get the data of the displacements
            int numComp = 0; //This integer is needed to call data() in the line below
            const float* const data = val.data(numComp);
            //Get the data of the previous displacements
            int numComp_prev = 0; //This integer is needed to call data() in the line below
            const float* const data_prev = val_prev.data(numComp_prev);

            //Displacements are in data, where:
            //data[0] corresponds to displacement in x
            //data[1] corresponds to displacement in y
            //data[2] corresponds to displacement in z            
            //Add the incremental displacement to point P
            // 
            //Note data[i] is the displacement for frame i, no the displacement
            //with respect to the previous frame
            //Thus, calculate the increment of displacement with respect to the previous frame
            //Otherwise I would need the initial position of the points for each frame
            points_cnt[P].x = points_cnt[P].x + (double)data[0] - (double)data_prev[0];
            points_cnt[P].y = points_cnt[P].y + (double)data[1] - (double)data_prev[1];
            points_cnt[P].z = points_cnt[P].z + (double)data[2] - (double)data_prev[2];
        }
    }

    return 1;
}
//This function generates the set name of CNT i, which follows the convention: CNT-i-NODES
string App_Network_From_Abaqus::Get_cnt_set_name(const int& cnt_i)const
{
    return ("CNT-"+to_string(cnt_i)+"-NODES");
}
//===========================================================================
