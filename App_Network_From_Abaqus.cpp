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

    //Shell vectors (used to remove nanoparticles when reducing observation window size)
    vector<vector<int> > shells_cnts;

    //----------------------------------------------------------------------
    //Read CNTs and GNPs from csv file (generated by Generate_Network) and 
    //fill the CNT and GNP variables
    //Also read the sample geometry
    ct0 = time(NULL);
    if (!Generate_nanoparticle_network_from_file(Init->simu_para, Init->vis_flags, Init->geom_sample, points_cnt, radii, structure_cnt, gnps))
    {
        hout << "Error when generating a nanoparticle network from a file" << endl;
        return 0;
    }
    ct1 = time(NULL);
    hout << "Network from file time: " << (int)(ct1 - ct0) << " secs." << endl << endl;

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
            if (!Apply_displacements_from_Abaqus(Init->simu_para.particle_type, (int)structure_cnt.size(), structure_cnt, root_assy, allFramesInStep[i], Init->geom_sample, points_cnt, gnps))
            {
                hout << "Error when applying displacement to nanoparticles." << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Apply displacements forcurrent frame time: " << (int)(ct1 - ct0) << " secs." << endl;
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

        //GNP structure, each structure_gnp[i] referes to the points in GNP_i
        vector<vector<long int> > structure_gnp;

        //----------------------------------------------------------------------
        //Determine the local networks in cutoff windows
        Cutoff_Wins* Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        ct0 = time(NULL);
        //hout<<"Extract_observation_window"<<endl;
        if (!Cutwins->Extract_observation_window(i, Init->simu_para.particle_type, Init->geom_sample, window_geo, Init->nanotube_geo, gnps, structure_cnt, radii, points_cnt, shells_cnts, shell_gnps, structure_gnp, points_gnp)) {
            hout << "Error when extracting observation window " << i + 1 << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Extract observation window time: " << (int)(ct1 - ct0) << " secs." << endl;

        //----------------------------------------------------------------------
        //Determine the local networks inside the cutoff windows
        Contact_grid* Contacts = new Contact_grid;
        ct0 = time(NULL);
        if (!Contacts->Generate_contact_grid(i, Init->simu_para.particle_type, Init->geom_sample, window_geo, Cutwins->cnts_inside, points_cnt, structure_cnt, Cutwins->gnps_inside, gnps)) {
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
        delete Cutwins;
        delete HoKo;

        it1 = time(NULL);
        hout << "Frame " << i << " time: " << (int)(it1 - it0) << " secs." << endl;
    }

    //Close Abaqus database
    odb.close();

    //Finalize usage of Abaqus C++ API
    odb_finalizeAPI();

    return 1;
}
//This function reads the data from a csv file to generate a nanoparticle network
int App_Network_From_Abaqus::Generate_nanoparticle_network_from_file(const Simu_para& simu_para, const Visualization_flags& vis_flags, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure, vector<GNP>& gnps)const
{
    //Read sample geometry
    if (!Read_sample_geometry(geom_sample))
    {
        hout << "Error in Generate_nanoparticle_network_from_file when reading the sample geometry from file" << endl;
        return 0;
    }

    //Check the type of nanoparticle 
    if (simu_para.particle_type == "CNT_wires" || simu_para.particle_type == "CNT_deposit" || simu_para.particle_type == "GNP_CNT_mix") {

        //Read the structure and CNT points
        if (!Read_cnt_data_from_csv(points_cnt, radii, structure))
        {
            hout << "Error in Generate_nanoparticle_network_from_file when calling Read_cnt_data_from_csv." << endl;
            return 0;

        }
    }
    else if (simu_para.particle_type == "GNP_cuboids" || simu_para.particle_type == "GNP_CNT_mix") {

        //Read the GNP geometry
        if (!Read_gnp_data_from_csv(geom_sample.sample, gnps))
        {
            hout << "Error in Generate_nanoparticle_network_from_file when calling Read_gnp_data_from_csv." << endl;
            return 0;
        }
    }
    else if (simu_para.particle_type == "Hybrid_particles") {
        hout << "Hybrid particles not yet implemented" << endl;
        return 0;
    }
    else {
        hout << "Error in Generate_nanoparticle_network_from_file: the type of particles should be one of the following: CNT_wires, CNT_deposit, GNP_cuboids, or GNP_CNT_mix. Input value was: " << simu_para.particle_type << endl;
        return 0;
    }

    //---------------------------------------------------------------------------
    //Check if visualization files were requested for read nanoparticles
    if (vis_flags.generated_nanoparticles) {

        VTK_Export vtk_exp;

        //Export generated CNTs if any
        if (points_cnt.size()) {
            vtk_exp.Export_from_cnt_structure(points_cnt, structure, "cnts_read.vtk");
        }

        //Export generated GNPs if any
        if (gnps.size()) {
            vtk_exp.Export_gnps(gnps, "gnps_read.vtk");
        }

        //Export the sample geometry
        vtk_exp.Export_cuboid(geom_sample.sample, "sample.vtk");
    }

    return 1;
}
//This function reads the CNT data from csv files (cnt_struct.csv and cnt_coordinates.csv) 
int App_Network_From_Abaqus::Read_cnt_data_from_csv(vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure)const
{
    //Open the structure file
    ifstream struc_file;
    struc_file.open("cnt_struct.csv");
    if (!struc_file) { 
        hout << "Failed to open CNT structure file cnt_struct.csv." << endl;  
        return 0; 
    }

    //String to read lines from file
    string line;

    //Integers to store the number of CNTs
    int n_cnts, n_cnts_check;

    //Read the first line line form the file and store it in a string stream
    getline(struc_file, line);
    stringstream ss_n_cnts(line);

    //Read the number of CNTs from the string stream and ignore the comma
    ss_n_cnts >> n_cnts; ss_n_cnts.ignore();
    ss_n_cnts >> n_cnts_check;

    //Check that both numbers in the first line of the structure file are the same
    if (n_cnts != n_cnts_check) {
        hout << "Error in Read_cnt_data_from_csv. The first line of file cnt_struct.csv indicates different number of CNTs: " << n_cnts << " and " << n_cnts_check << endl;
        return 0;
    }

    //Set the size of the structure vector to store the indicated number of CNTs
    structure.assign(n_cnts, vector<long int>());

    //Set the size of the redius vector
    radii.assign(n_cnts, 0.0);

    //Variable to store the total number of CNT points in the network
    long int Np = 0;

    //Read the structure
    for (int i = 0; i < n_cnts; i++)
    {
        //Check if end-of-file has been reached
        if (struc_file.eof())
        {
            hout << "Error in Read_cnt_data_from_csv. The end-of-file of cnt_struct.csv has been reached before reading all CNT data." << endl;
            return 0;
        }

        //Variable to store the number of points in CNT i
        int np;

        //Read a line form the file and store it in a string stream
        getline(struc_file, line);
        stringstream ss(line);

        //Read the number of CNT points from the string stream
        ss >> np;
        //Ignore the comma in the line
        ss.ignore();

        //Read the radius from the string stream
        ss >> radii[i];

        //Fill structure vector for CNT i
        structure[i].assign(np, Np);
        //Increase the number of points since the current values is already used 
        //for the first point of CNT i
        Np++;

        //Iterate overt the remaining points in CNT i to fill the structure[i] vector
        //index i = 0 is ignored since that is the value with which 
        for (int j = 1; j < np; j++)
        {
            //Set the j-th point in CNT i
            structure[i][j] = Np;

            //Increase the count for the number of CNT points
            Np++;
        }
        //hout << structure[i].size() << " " << radii[i] << endl;
    }

    //Open the point coordinates file
    ifstream  coord_file;
    coord_file.open("cnt_coordinates.csv");
    if (!coord_file) {
        hout << "Failed to open point cordinates file cnt_coordinates.csv." << endl;
        return 0;
    }

    //Set the size for the points vector
    points_cnt.assign(Np, Point_3D());

    //Read the point coordinates
    for (long int i = 0; i < Np; i++)
    {
        //Check if end-of-file has been reached
        if (coord_file.eof())
        {
            hout << "Error in Read_cnt_data_from_csv. The end-of-file of cnt_coordinates.csv has been reached before reading all CNT data." << endl;
            return 0;
        }

        //Read a line form the file and store it in a string stream
        getline(coord_file, line);
        stringstream ss(line);

        //Read the point coordinates from the string stream while ignoring the commas in between
        ss >> points_cnt[i].x;
        ss.ignore();
        ss >> points_cnt[i].y;
        ss.ignore();
        ss >> points_cnt[i].z;
        
    }

    //Close files
    struc_file.close();
    coord_file.close();

    //Output a message with the number of CNTs and points read
    hout << endl << "A total of " << n_cnts << " CNTs and " << Np << " points were read." << endl;

    return 1;
}
//This function reads GNP data from a csv file (gnp_data.csv)
int App_Network_From_Abaqus::Read_gnp_data_from_csv(const cuboid& sample_geom, vector<GNP>& gnps)const
{
    //Open the file with the GNP geometric data
    ifstream gnp_file;
    gnp_file.open("gnp_data.csv");
    if (!gnp_file) {
        hout << "Failed to open file with GNP geometric data gnp_data.csv." << endl;
        return 0;
    }

    //Generate_Network object to generate the variables needed for each GNP
    Generate_Network GN;

    //String to read lines from file
    string line;

    //Read the file
    while (getline(gnp_file, line))
    {
        //Read a line from the file and it in a string stream
        stringstream ss(line);

        //GNP to store the values read from the file
        GNP new_gnp;

        //Variables to store the angles
        double theta, phi;

        //Temporary variable to ignore repetition of GNP side-length
        //For some reason using ss.ignore() three times in a roll results in an incorrect
        //reding of the csv file
        double l_tmp;

        //Read the values while ignoring the commas
        ss >> new_gnp.l; ss.ignore(); 
        ss >> l_tmp;  ss.ignore();
        ss >> new_gnp.t; ss.ignore();
        ss >> theta; ss.ignore();
        ss >> phi; ss.ignore();
        ss >> new_gnp.center.x; ss.ignore();
        ss >> new_gnp.center.y; ss.ignore();
        ss >> new_gnp.center.z; 
        //hout << new_gnp.l << " " << new_gnp.t << " " << theta << " " << phi << " " << new_gnp.center.str() << endl;

        //Calculate rotation matrix
        new_gnp.rotation = GN.Get_transformation_matrix(theta, phi);

        //Calculate (or approximate) GNP volume
        int i;
        for (i = 0; i < 8; i++)
        {
            if (!GN.Is_point_inside_cuboid(sample_geom, new_gnp.vertices[i])) {

                //GNP is partially outside the sample, so approximate volume
                double gnp_vol = 0.0;
                if (!GN.Approximate_gnp_volume_inside_sample(sample_geom, new_gnp, gnp_vol)) {
                    hout << "Error in Read_gnp_data_from_csv when calling Approximate_gnp_volume_inside_sample." << endl;
                    return 0;
                }
                new_gnp.volume = gnp_vol;

                //Terminate the for loop
                break;
            }
        }

        //Check if all vertices were inside the sample
        if (i == 7)
        {
            //All vertices were inside the sample, thus calcualte the volume of the whole GNP
            new_gnp.volume = new_gnp.l * new_gnp.l * new_gnp.t;
        }

        //Obtain coordinates of GNP vertices
        if (!GN.Obtain_gnp_vertex_coordinates(new_gnp)) {
            hout << "Error in Read_gnp_data_from_csv when calling Obtain_gnp_vertex_coordinates." << endl;
            return 0;
        }

        //Get the plane equations for the six faces
        if (!GN.Update_gnp_plane_equations(new_gnp)) {
            hout << "Error in Read_gnp_data_from_csv when calling Update_gnp_plane_equations." << endl;
            return 0;
        }

        //Add GNP to the vector of GNPs
        gnps.push_back(new_gnp);
    }

    //Close file
    gnp_file.close();

    //Output a message with the total number of GNPs read
    hout << endl << "A total of " << gnps.size() << " GNPs were read." << endl;

    return 1;
}
//This function reads the sample geomtry from a csv file (sample_geom.csv)
int App_Network_From_Abaqus::Read_sample_geometry(Geom_sample& geom_sample)const
{
    //Open the file with the sample geometry data
    ifstream sample_file;
    sample_file.open("sample_geom.csv");
    if (!sample_file) {
        hout << "Failed to open file with sample geometry data sample_geom.csv." << endl;
        return 0;
    }

    //String to read lines from file
    string line;

    //Read the first line form the file and store it in a string stream
    getline(sample_file, line);
    stringstream ss_point(line);

    //Read the lower left corner of the sample while ignoring the commas in between
    ss_point >> geom_sample.sample.poi_min.x;
    ss_point.ignore();
    ss_point >> geom_sample.sample.poi_min.y;
    ss_point.ignore();
    ss_point >> geom_sample.sample.poi_min.z;

    //Read the second line form the file and store it in a string stream
    getline(sample_file, line);
    stringstream ss_size(line);

    //Read the dimensions of the sample along each direction while ignoring the commas in between
    ss_size >> geom_sample.sample.len_x;
    ss_size.ignore();
    ss_size >> geom_sample.sample.wid_y;
    ss_size.ignore();
    ss_size >> geom_sample.sample.hei_z;

    //Calculate the sample's volume
    geom_sample.volume = geom_sample.sample.len_x * geom_sample.sample.wid_y * geom_sample.sample.hei_z;

    //Calculate the coordinates of the sample's boundaries opposite to those given by the coordinates of origin
    geom_sample.sample.max_x = geom_sample.sample.poi_min.x + geom_sample.sample.len_x;
    geom_sample.sample.max_y = geom_sample.sample.poi_min.y + geom_sample.sample.wid_y;
    geom_sample.sample.max_z = geom_sample.sample.poi_min.z + geom_sample.sample.hei_z;

    //Close file
    sample_file.close();

    return 1;
}
//This functions adds the displacements to the CNTs, GNPs and sample
int App_Network_From_Abaqus::Apply_displacements_from_Abaqus(const string& particle_type, const int& n_cnts, const vector<vector<long int> >& structure, odb_Assembly& root_assy, odb_Frame& current_frame, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<GNP>& gnps)const
{
    //Access displacement field ("U") in the current frame
    //hout << "fieldU" << endl;
    odb_FieldOutput& fieldU = current_frame.fieldOutputs()["U"];

    //Apply displacements to sample
    //hout << "Apply_displacements_to_sample" << endl;
    if (!Apply_displacements_to_sample(root_assy, fieldU, geom_sample))
    {
        hout << "Error in Apply_displacements_from_Abaqus when calling Apply_displacements_to_sample" << endl;
        return 0;
    }

    //Check the type of nanoparticle 
    if (particle_type == "CNT_wires" || particle_type == "CNT_deposit" || particle_type == "GNP_CNT_mix") {

        //Apply displacements to CNT points
        //hout << "Apply_displacements_to_cnts" << endl;
        if (!Apply_displacements_to_cnts(structure, root_assy, fieldU, n_cnts, points_cnt))
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
int App_Network_From_Abaqus::Apply_displacements_to_sample(odb_Assembly& root_assy, odb_FieldOutput& fieldU, Geom_sample& geom_sample)const
{
    //Name of the set for the lower left corner (it is hard coded in the python scritp too)
    string set0 = "MATRIX0";

    //Name of the set for the opposite corner (it is hard coded in the python scritp too)
    string set1 = "MATRIX1";

    //Access set0 from root assembly
    odb_Set& matrix0 = root_assy.nodeSets()[set0.c_str()];

    //Get the displacement object of the set
    odb_FieldOutput matrix0_disp = fieldU.getSubset(matrix0);

    //Get the sequence of values of the displacement object for matrix0
    const odb_SequenceFieldValue& vals0 = matrix0_disp.values();

    //Check the size of vals is 1 (since it contains one node)
    if (vals0.size() != 1)
    {
        hout << "Error in Apply_displacements_to_sample. The number of node set " << set0 << " is not 1. Size is " << vals0.size() << endl;
        return 0;
    }

    //Get the acutal values
    const odb_FieldValue val0 = vals0[0];
    //Output node label
    //hout << "  Node: " << val0.nodeLabel() << endl;
    //Get the data of the displacements
    int numComp = 0; //This integer is needed to call data() in the line below
    const float* const data0 = val0.data(numComp);

    //Update lower left corner of sample
    geom_sample.sample.poi_min = geom_sample.sample.poi_min + Point_3D((double)data0[0], (double)data0[1], (double)data0[2]);

    //Access set1 from root assembly
    odb_Set& matrix1 = root_assy.nodeSets()[set1.c_str()];

    //Get the displacement object of the set
    odb_FieldOutput matrix1_disp = fieldU.getSubset(matrix1);

    //Get the sequence of values of the displacement object for matrix1
    const odb_SequenceFieldValue& vals1 = matrix1_disp.values();

    //Check the size of vals is 1 (since it contains one node)
    if (vals1.size() != 1)
    {
        hout << "Error in Apply_displacements_to_sample. The number of node set " << set1 << " is not 1. Size is " << vals1.size() << endl;
        return 0;
    }

    //Get the acutal values
    const odb_FieldValue val1 = vals1[0];
    //Output node label
    //cout << "  Node: " << val1.nodeLabel() << endl;
    //Get the data of the displacements
    numComp = 0; //This integer is needed to call data() in the line below
    const float* const data1 = val1.data(numComp);

    //Update the maximum coordinates of the sample, which are the maximum coordinates of the sample
    geom_sample.sample.max_x = geom_sample.sample.max_x + (double)data1[0];
    geom_sample.sample.max_y = geom_sample.sample.max_y + (double)data1[1];
    geom_sample.sample.max_z = geom_sample.sample.max_z + (double)data1[2];

    //Update the dimensions of the sample along each direction
    geom_sample.sample.len_x = geom_sample.sample.max_x - geom_sample.sample.poi_min.x;
    geom_sample.sample.wid_y = geom_sample.sample.max_y - geom_sample.sample.poi_min.y;
    geom_sample.sample.hei_z = geom_sample.sample.max_z - geom_sample.sample.poi_min.z;

    //Update the sample's volume
    geom_sample.volume = geom_sample.sample.len_x * geom_sample.sample.wid_y * geom_sample.sample.hei_z;

    return 1;
}
//This function applies displacements to CNT points
int App_Network_From_Abaqus::Apply_displacements_to_cnts(const vector<vector<long int> >& structure, odb_Assembly& root_assy, odb_FieldOutput& fieldU, const int& n_cnts, vector<Point_3D>& points_cnt)const
{

    //Iterate over the CNTs in the sample
    for (size_t i = 0; i < structure.size(); i++)
    {
        //Get the name of the set
        //Numbering of sets starts at 1
        string set_name = Get_cnt_set_name((int)i+1);

        //Access set from root assembly
        odb_Set& cnt_set = root_assy.nodeSets()[set_name.c_str()];

        //Get the displacement object of the set
        odb_FieldOutput cnt_disp = fieldU.getSubset(cnt_set);

        //Get the values of the displacement object
        const odb_SequenceFieldValue& vals = cnt_disp.values();
        //Output size of values
        //cout << "values.size=" << vals.size() << endl;

        //Check there is the same number of points in the CNT and in the set
        if (structure[i].size() != vals.size())
        {
            hout << "Error in Apply_displacements_to_cnts. The number of points in CNT " << i << " (" << structure[i].size() << " points) is different from the number of points in set " << set_name << " ("<< vals.size() << " points)." << endl;
            return 0;
        }

        //Iterate over the values, i.e., the points in the current CNT
        //Note that nodes are actually in reverse order given the way
        //CNTs are generated (and meshed) in Abaqus
        for (int j = vals.size() - 1; j >= 0 ; j--)
        {
            //Get current point number
            long int P = structure[i][j];
            
            //Get current values
            const odb_FieldValue val = vals[j];
            //Output node label
            //cout << "  Node: " << val.nodeLabel() << endl;
            //Get the data of the displacements
            int numComp = 0; //This integer is needed to call data() in the line below
            const float* const data = val.data(numComp);

            //Displacements are in data, where:
            //data[0] corresponds to displacement in x
            //data[1] corresponds to displacement in y
            //data[2] corresponds to displacement in z            
            //Add these displacement to point P
            points_cnt[P].x = points_cnt[P].x + (double)data[0];
            points_cnt[P].y = points_cnt[P].y + (double)data[1];
            points_cnt[P].z = points_cnt[P].z + (double)data[2];
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
