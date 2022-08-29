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

    //Make a copy of the GNPs in their initial position (as they were generated)
    vector<GNP> gnps0(gnps.begin(), gnps.end());

    //Set the number of cuts to 0 since there are no observation windows when reading
    //a network from file and applying displacements from abaqus
    Init->geom_sample.cut_num = 0;

    //Vector to store GNPs that lie partially outside the sample
    vector<int> gnps_outside;
    //Vector with vertices inside the sample from the GNPs that lie partially outside the sample
    vector<vector<int> > vertices_gnps_out;
    //Vector with flags that indicate which vertices are inside the sample (true)
    //and which are outside the sample (false)
    vector<vector<bool> > vertex_flags;

    //Generate vector with GNPs partially outside the sample
    if (!Get_gnps_partially_outside_sample(Init->geom_sample, gnps, gnps_outside, vertices_gnps_out, vertex_flags))
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
    //hout << "open odb=" << Init->simu_para.odb_file.c_str() << endl;

    //Access the root assembly
    odb_Assembly& root_assy = odb.rootAssembly();
    //hout << "root_assy" << endl;

    //Make sure the step indicated in the input file is in the odb file
    if (!Is_step_in_odb(odb, Init->simu_para.step_name))
    {
        hout << "Error in Nanoparticle_resistor_network_from_odb when calling Is_step_in_odb." << endl;
        return 0;
    }

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
        //Also perform again window extraction for the GNPs only
        if (i)
        {
            //Read Abaqus database and apply displacements
            ct0 = time(NULL);
            //hout<<"Apply_displacements_from_Abaqus"<<endl;
            if (!Apply_displacements_from_Abaqus(Init->simu_para.particle_type, (int)structure_cnt.size(), structure_cnt, gnps_outside, vertices_gnps_out, vertex_flags, gnps0, root_assy, allFramesInStep[i-1], allFramesInStep[i], Init->geom_sample, points_cnt, gnps))
            {
                hout << "Error when applying displacement to nanoparticles." << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Apply all displacements for current frame time: " << (int)(ct1 - ct0) << " secs." << endl;

            //Clear the GNP points
            points_gnp.clear();

            //Clear the GNP vectors that are filled with extracting an observation window
            Cutwins->boundary_gnp.clear();
            Cutwins->boundary_gnp_pts.clear();
            Cutwins->gnps_inside.clear();

            //Find the GNPs and GNP points at boundaries
            ct0 = time(NULL);
            //hout<<"Extract_observation_window (GNPs only)"<<endl;
            //For the case of reading data from an Abaqus database, window is always 0
            if (!Cutwins->Extract_observation_window(0, "GNP_cuboids", Init->geom_sample, Init->geom_sample.sample, Init->nanotube_geo, gnps, structure_cnt, vector<double>(), vector<Point_3D>(), vector<vector<int> >(), shell_gnps, structure_gnp, points_gnp)) {
                hout << "Error when extracting observation window " << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Extract observation window time: " << (int)(ct1 - ct0) << " secs." << endl;
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
        //Set flag for ignoring error in mexed contacts to 1
        //This means, error is ignored
        HoKo->ignore_eror_cnt_inside_gnp = 1;
        ct0 = time(NULL);
        if (!HoKo->Determine_clusters_and_percolation(i, Init->geom_sample.sample, Init->simu_para, Init->cutoff_dist, Init->vis_flags, Cutwins->cnts_inside, Contacts->sectioned_domain_cnts, structure_cnt, points_cnt, radii, Cutwins->boundary_cnt, Cutwins->gnps_inside, Contacts->sectioned_domain_gnps, gnps, Cutwins->boundary_gnp, structure_gnp, points_gnp)) {
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
int App_Network_From_Abaqus::Get_gnps_partially_outside_sample(const Geom_sample& geom_sample, const vector<GNP>& gnps, vector<int>& gnps_outside, vector<vector<int> >& vertices_gnps_out, vector<vector<bool> > &vertex_flags)const
{
    //Iterate over the GNPs
    for (int i = 0; i < (int)gnps.size(); i++)
    {
        //Vector to store vertices inside the sample
        vector<int> vertices_tmp;

        //Vector to store the vertex flags of current GNP
        vector<bool> vf_tmp(8, false);

        //Check if all eight vertices are inside the sample
        for (int j = 0; j < 8; j++)
        {
            //Check if vertex j of GNP i is inside the sample
            if (!gnps[i].vertices[j].is_outside_cuboid(geom_sample.sample))
            {
                //Add vertex number to temporary vector
                vertices_tmp.push_back(j);

                //Set the vertex flag as true
                vf_tmp[j] = true;
            }
        }

        //Check if not all vertices are inside the sample
        if (vertices_tmp.size() != 8)
        {
            //Check that there are at least three vertices inside the sample
            if (vertices_tmp.size() < 3)
            {
                hout << "Error in Get_gnps_partially_outside_sample." << endl;
                hout << "There is a GNP that has less than three vertices inside the sample." << endl;
                hout << "GNPs should have at least three vertices inside the sample." << endl;
                string tmp_str = (vertices_tmp.size() == 1) ? " vertex " : " vertices ";
                hout << "GNP #" << i << " has " << vertices_tmp.size() << tmp_str << "inside the sample." << endl;
                return 0;
            }

            //GNP i has some vertices outside the sample, so add GNP i to the vector
            //of GNPs partially outside the sample
            gnps_outside.push_back(i);

            //Add the vertices of GNP i that are inside the sample to the vector vertices_gnps_out
            vertices_gnps_out.push_back(vertices_tmp);

            //Add vertex flags of GNP i to the vector vertex_flags
            vertex_flags.push_back(vf_tmp);

            /*hout << "Inside vertices of GNP " << i << endl << "\t";
            for (size_t k = 0; k < vertices_tmp.size(); k++) {
                hout << vertices_tmp[k] << ' ';
            } hout << endl;*/
        }
        /*hout << "Vertices of GNP " << i << endl << "\t";
        for (size_t k = 0; k < vertices_tmp.size(); k++) {
            hout << vertices_tmp[k] << ' ';
        } hout << endl;*/
    }

    return 1;
}
//This function determines if the step indicated in the input file is in the odb file
int App_Network_From_Abaqus::Is_step_in_odb(const odb_Odb& odb, const string& step_name)const
{
    //Check if there is at least one step
    if (odb.steps().size() < 1)
    {
        hout << "Error in Is_step_in_odb: There are no steps in the odb file." << endl;
        return 0;
    }

    //Check that the step exists in the odb file
    odb_StepRepositoryIT stepIter(odb.steps());
    //hout << "from input file:" << step_name << endl;
    for (stepIter.first(); !stepIter.isDone(); stepIter.next())
    {
        if (step_name == stepIter.currentKey().CStr())
            return 1;
        //hout << stepIter.currentKey().CStr() << endl;
    }
    //hout << endl;

    hout << "Error in Is_step_in_odb: Step name indicated in input file was not found in odb file." << endl;
    hout << "Step name in input file is: " << step_name << "." << endl;
    hout << "This is a list of all steps in odb file:" << endl;
    for (stepIter.first(); !stepIter.isDone(); stepIter.next())
    {
        hout << stepIter.currentKey().CStr() << endl;
    }

    return 0;
}
//This functions adds the displacements to the CNTs, GNPs and sample
int App_Network_From_Abaqus::Apply_displacements_from_Abaqus(const string& particle_type, const int& n_cnts, const vector<vector<long int> >& structure, const vector<int>& gnps_outside, const vector<vector<int> >& vertices_gnps_out, const vector<vector<bool> >& vertex_flags, vector<GNP>& gnps0, odb_Assembly& root_assy, odb_Frame& previous_frame, odb_Frame& current_frame, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<GNP>& gnps)const
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

    //Check if there are CNTs
    if (particle_type == "CNT_wires" || particle_type == "CNT_deposit" || particle_type == "GNP_CNT_mix") 
    {
        time_t it0 = time(NULL);
        //Apply displacements to CNT points
        //hout << "Apply_displacements_to_cnts" << endl;
        if (!Apply_displacements_to_cnts(structure, root_assy, previous_fieldU, current_fieldU, n_cnts, points_cnt))
        {
            hout << "Error in Apply_displacements_from_Abaqus when calling Apply_displacements_to_cnts" << endl;
            return 0;
        }

        //If mixed particles, output time for CNTs
        if (particle_type == "GNP_CNT_mix") 
        {
            time_t it1 = time(NULL);
            hout << "Apply displacements for CNTs: " << (int)(it1 - it0) << " secs." << endl;
        }
    }

    //Check if there are GNPs
    if (particle_type == "GNP_cuboids" || particle_type == "GNP_CNT_mix") 
    {
        time_t it0 = time(NULL);
        //Apply displacements to GNP vertices
        //hout << "Apply_displacements_to_gnps" << endl;
        if (!Apply_displacements_to_gnps(gnps_outside, vertices_gnps_out, vertex_flags, gnps0, root_assy, current_fieldU, gnps))
        {
            hout << "Error in Apply_displacements_from_Abaqus when calling Apply_displacements_to_gnps" << endl;
            return 0;
        }

        //If mixed particles, output time for CNTs
        if (particle_type == "GNP_CNT_mix")
        {
            time_t it1 = time(NULL);
            hout << "Apply displacements for GNPs: " << (int)(it1 - it0) << " secs." << endl;
        }
    }

    return 1;
}
//This function gets the displacements of two of the sample nodes, applies this displacemets
//to the corresponding vertices and calculates the geometry of the deformed sample
int App_Network_From_Abaqus::Apply_displacements_to_sample(odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, Geom_sample& geom_sample)const
{
    //Name of the set for the lower left corner (it is hard coded in the python scritp too)
    string set0 = "MATRIX0";

    //Variables to store displacement vector of the current frame 
    //with respect to the previous frame 
    Point_3D disp;

    //Get the displacement for set0
    if (!Get_displacement_change_from_single_node_set(set0, root_assy, previous_fieldU, current_fieldU, disp)) {
        hout << "Error in Apply_displacements_to_sample when calling Get_displacement_change_from_single_node_set (set: "<<set0<<")" << endl;
        return 0;
    }

    //Manually adjust the displacement if some CNT elements from Abaqus lie outside the RVE
    //disp.x = disp.x * 3.0 / 3.02; disp.z = disp.z * 3.0 / 3.02;

    //Update lower left corner of sample
    //hout << "P0=" << geom_sample.sample.poi_min.str() << " disp=" << disp.str() << " P0+dips=" << (geom_sample.sample.poi_min + disp).str() << endl;
    geom_sample.sample.poi_min = geom_sample.sample.poi_min + disp;

    //Name of the set for the opposite corner (it is hard coded in the python scritp too)
    string set1 = "MATRIX1";

    //Get the displacement for set1
    if (!Get_displacement_change_from_single_node_set(set1, root_assy, previous_fieldU, current_fieldU, disp)) {
        hout << "Error in Apply_displacements_to_sample when calling Get_displacement_change_from_single_node_set (set: " << set1 << ")" << endl;
        return 0;
    }

    //Manually adjust the displacement if some CNT elements from Abaqus lie outside the RVE
    //disp.x = disp.x * 3.0 / 3.02; disp.z = disp.z * 3.0 / 3.02;
    
    //Update the maximum coordinates of the sample
    //hout << "corner=(" << geom_sample.sample.max_x << ", " << geom_sample.sample.max_y << ", " << geom_sample.sample.max_z << ") disp=" << disp.str() << " corner+disp=(" << geom_sample.sample.max_x + disp.x << ", " << geom_sample.sample.max_y + disp.y << ", " << geom_sample.sample.max_z + disp.z << ")" << endl;
    geom_sample.sample.max_x = geom_sample.sample.max_x + disp.x;
    geom_sample.sample.max_y = geom_sample.sample.max_y + disp.y;
    geom_sample.sample.max_z = geom_sample.sample.max_z + disp.z;

    //Update the dimensions of the sample along each direction
    geom_sample.sample.len_x = geom_sample.sample.max_x - geom_sample.sample.poi_min.x;
    geom_sample.sample.wid_y = geom_sample.sample.max_y - geom_sample.sample.poi_min.y;
    geom_sample.sample.hei_z = geom_sample.sample.max_z - geom_sample.sample.poi_min.z;

    //Update the sample's volume
    geom_sample.volume = geom_sample.sample.len_x * geom_sample.sample.wid_y * geom_sample.sample.hei_z;

    return 1;
}
//For the case of a set that contains a single node, this function obtaines the displacement
//at that node with respect to the previous frame
int App_Network_From_Abaqus::Get_displacement_change_from_single_node_set(const string& setname, odb_Assembly& root_assy, odb_FieldOutput& previous_fieldU, odb_FieldOutput& current_fieldU, Point_3D& disp)const
{
    //Variables to store current and previous displacements
    Point_3D current_disp, previous_disp;

    //Get the displacement for the current frame for set
    if (!Get_displacement_from_single_node_set(setname, root_assy, current_fieldU, current_disp)) {
        hout << "Error in Get_displacement_change_from_single_node_set when calling Get_displacement_from_single_node_set (current_disp)" << endl;
        return 0;
    }

    //Get the displacement for the previous frame for set
    if (!Get_displacement_from_single_node_set(setname, root_assy, previous_fieldU, previous_disp)) {
        hout << "Error in Get_displacement_change_from_single_node_set when calling Get_displacement_from_single_node_set (previous_disp)" << endl;
        return 0;
    }

    //Calculate the displacement with respect to the previous frame
    disp = current_disp - previous_disp;

    return 1;
}
//For the case of a set that contains a single node, this function obtaines the displacement
//at that node
int App_Network_From_Abaqus::Get_displacement_from_single_node_set(const string& setname, odb_Assembly& root_assy, odb_FieldOutput& fieldU, Point_3D& disp)const
{
    //Access set from root assembly
    odb_Set& sub_set = root_assy.nodeSets()[setname.c_str()];

    //Get the displacement objects of the set
    odb_FieldOutput set_disp = fieldU.getSubset(sub_set);

    //Get the sequence of values of the displacement object for matrix1
    const odb_SequenceFieldValue& vals = set_disp.values();

    //Check the size of vals is 1 (since it contains one node)
    if (vals.size() != 1)
    {
        hout << "Error in Get_displacement_from_single_node_set. The number of node set " << setname << " is not 1. Size is " << vals.size() << endl;
        return 0;
    }

    //Get the acutal values
    const odb_FieldValue val1 = vals[0];
    //Output node label
    //cout << "  Node: " << val1.nodeLabel() << endl;
    //Get the data of the displacements
    int numComp = 0; //This integer is needed to call data() in the line below
    const float* const data1 = val1.data(numComp);
    //hout << "disp1=" << data1[0] << " " << data1[1] << " " << data1[2] << endl;

    //Use the values stored in data1 as the components of the point disp
    disp.set((double)data1[0], (double)data1[1], (double)data1[2]);

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
        //hout << "Set " << set_name << endl;

        //Variables to store displacement objects
        odb_FieldOutput cnt_disp, cnt_disp_prev;

        //Sometimes abaqus does not generate the set, so catch that error if this happens
        try {
            //Access set from root assembly
            odb_Set& cnt_set = root_assy.nodeSets()[set_name.c_str()];

            //Get the displacement objects of the set
            cnt_disp = current_fieldU.getSubset(cnt_set);
            //hout << "current_fieldU" << endl;
            cnt_disp_prev = previous_fieldU.getSubset(cnt_set);
            //hout << "previous_fieldU" << endl;
        }
        catch (...) {
            hout << "Error in Apply_displacements_to_cnts." << endl;
            hout << "Error while accessing set " << set_name.c_str() << endl;
            return 0;
        }
        //hout << "nodeSets" << endl;

        //Get the values of the displacement object
        const odb_SequenceFieldValue& vals = cnt_disp.values();
        //hout << "cnt_disp" << endl;
        const odb_SequenceFieldValue& vals_prev = cnt_disp_prev.values();
        //hout << "cnt_disp_prev" << endl;
        //Output size of values
        //cout << "values.size=" << vals.size() <<" structure[CNT="<<i<<"].size()="<< structure[i].size() << endl;

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
            //Thus, I need to go in reverse order when scanning the structure vector
            //With this definition of idx, I can start at the last element of 
            //structure[i] and stop at its first element
            int idx = v_size1 - j;

            //Get current point number
            long int P = structure[i][idx];
            //hout << "P=" << P <<" idx="<< idx << endl;
            
            //Get current values
            const odb_FieldValue val = vals[j];
            //Get previous values
            const odb_FieldValue val_prev = vals_prev[j];
            //Output node label
            //hout << "  Node: " << val.nodeLabel() << endl;
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
            //hout << "P udated" << endl;
        }
        //hout << "for-end" << endl;
    }

    return 1;
}
//This function generates the set name of CNT i, which follows the convention: CNT-i-NODES
string App_Network_From_Abaqus::Get_cnt_set_name(const int& cnt_i)const
{
    return ("CNT-"+to_string(cnt_i)+"-NODES");
}
//This fucntion pplies displacements to GNPs
int App_Network_From_Abaqus::Apply_displacements_to_gnps(const vector<int>& gnps_outside, const vector<vector<int> >& vertices_gnps_out, const vector<vector<bool> >& vertex_flags, const vector<GNP>& gnps0, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU, vector<GNP>& gnps)const
{
    //Variable to store the index within vector gnps_outside 
    //of next GNP that is partially outside the sample
    int idx_gnp_out = 0;

    //Get the number of GNPs
    int n_gnps = (int)gnps.size();

    //Get the number of GNPs partially inside the sample
    int n_gnps_out = (int)gnps_outside.size();

    //Variable to store integers 0 to 7 to loop over all vertices of a GNP that
    //is completely inside the sample
    vector<int> all_vertices = All_gnp_vertices();

    //Flags for full GNPs (used when checking if points are inside the reconstructed GNP)
    vector<bool> all_true(8, true);

    //Object to reconstruct a GNP
    GNP_Reconstruction GR;

    //Iterate over the GNPs
    for (int i = 0; i < n_gnps; i++)
    {
        //Reset GNP i to its initial location
        gnps[i] = gnps0[i];

        //Check if the GNP is partially inside or completely inside
        if (i == gnps_outside[idx_gnp_out])
        {
            //GNP is partially inside the sample, so use the vector that contains
            //the vertices that are inside the sample
            //hout << endl << "Apply_displacements_to_gnp_vertices (1) GNP=" << i << endl;
            if (!Apply_displacements_to_gnp_vertices(vertices_gnps_out[idx_gnp_out], root_assy, current_fieldU, gnps[i]))
            {
                hout << "Error in Apply_displacements_to_gnps when calling Apply_displacements_to_gnp_vertices (1)" << endl;
                return 0;
            }

            //Reconstruct GNP
            //hout << "GR.Reconstruct_gnp 1" << endl;
            if (!GR.Reconstruct_gnp(vertices_gnps_out[idx_gnp_out], vertex_flags[idx_gnp_out], gnps[i]))
            {
                hout << "Error in Apply_displacements_to_gnps when calling GR.Reconstruct_gnp (partial)" << endl;
                return 0;
            }

            //Increase the index of the next available GNP in vector gnps_outside
            //Ensure that the value of the index does not exceed the number of
            //GNPs partially inside the sampe (n_gnps_out)
            idx_gnp_out = (idx_gnp_out + 1) % n_gnps_out;
        }
        else
        {
            //GNP is competely inside the sample, so use the vector that contains
            //all the GNP vertices
            //hout << endl << "Apply_displacements_to_gnp_vertices (2) GNP=" << i << endl;
            if (!Apply_displacements_to_gnp_vertices(all_vertices, root_assy, current_fieldU, gnps[i]))
            {
                hout << "Error in Apply_displacements_to_gnps when calling Apply_displacements_to_gnp_vertices (2)" << endl;
                return 0;
            }

            //Reconstruct GNP
            //Use an empty boolean vector since in the case of a full GNP the 
            //boolean vector is not needed
            //hout << "GR.Reconstruct_gnp 2" << endl;
            if (!GR.Reconstruct_gnp(all_vertices, all_true, gnps[i]))
            {
                hout << "Error in Apply_displacements_to_gnps when calling GR.Reconstruct_gnp (full)" << endl;
                return 0;
            }
        }
    }

    return 1;
}
//This function fills a vector with integers 0 to 7 to loop over all vertices of a GNP that
//lies completely inside the sample
vector<int> App_Network_From_Abaqus::All_gnp_vertices()const
{
    vector<int> idxs(8, 0);

    //Iterate over vertices 1 to 7 since 0 is already taken care of
    for (int i = 1; i < 8; i++) {
        idxs[i] = i;
    }

    return idxs;
}
//This function adds displacements to GNP vertices for the case when all vertices
//are inside the sample
int App_Network_From_Abaqus::Apply_displacements_to_gnp_vertices(const vector<int>& vertices, odb_Assembly& root_assy, odb_FieldOutput& current_fieldU, GNP& gnp_i)const
{
    //Iterate over the vertices of the GNP
    for (size_t j = 0; j < vertices.size(); j++)
    {
        //Get vertex number
        int v = vertices[j];
        //hout << "v=" << v << endl;

        //Get the set name for the node v in GNP i
        string setname = Get_gnp_set_name(gnp_i.flag, v);
        //hout << setname << endl;

        //Point to store the displacement vector with respect to the previous frame
        Point_3D disp;

        //Get the displacement for the current frame for set
        if (!Get_displacement_from_single_node_set(setname, root_assy, current_fieldU, disp)) {
            hout << "Error in Apply_displacements_to_gnp_vertices when calling Get_displacement_from_single_node_set" << endl;
            return 0;
        }

        //Update location of vertex v in GNP i
        gnp_i.vertices[v] = gnp_i.vertices[v] + disp;
    }

    return 1;
}
//This function generates the set name of GNP gnp_i, that contains the node vertex
//Set naming follows the convention: GS-i-_N-vertex
string App_Network_From_Abaqus::Get_gnp_set_name(const int& gnp_i, const int& vertex)const
{
    return ("GS-" + to_string(gnp_i) + "_N-" + to_string(vertex));
}
