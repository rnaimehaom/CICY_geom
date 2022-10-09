//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate electrical resistivity of a network read from dissplacements saved in binary files (obtained from an Abaqus output odb file)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "App_Network_From_Displacements.h"

//Generate 3D nanoparticle network, turn it into a resitor network and find its electrical conductivity
int App_Network_From_Displacements::Nanoparticle_resistor_network_from_displacements(Input* Init)const
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
        hout << "Error in Nanoparticle_resistor_network_from_displacements when calling Generate_nanoparticle_network_from_file." << endl;
        return 0;
    }
    delete Reader;
    ct1 = time(NULL);
    hout << "Network generation from file time: " << (int)(ct1 - ct0) << " secs." << endl;

    //Set the number of cuts to 0 since there are no observation windows when reading
    //a network from file and applying displacements from abaqus
    Init->geom_sample.cut_num = 0;

    //Vector with vertices inside the sample from the GNPs that lie partially outside the sample
    vector<vector<int> > vertices_gnps_in;

    //Generate vector with GNPs partially outside the sample
    if (!Get_gnps_partially_outside_sample(Init->geom_sample, gnps, vertices_gnps_in))
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

    //Iterate over the number of frames
    for (int i = 0; i < Init->simu_para.n_frames; i++)
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
            //hout<<"Apply_displacements_from_files"<<endl;
            if (!Apply_displacements_from_files(i, Init->simu_para.particle_type, structure_cnt, vertices_gnps_in, Init->geom_sample, points_cnt, gnps))
            {
                hout << "Error in Nanoparticle_resistor_network_from_displacements when calling Apply_displacements_from_files." << endl;
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
            vector<double> empty_double;
            vector<Point_3D> empty_point;
            vector<vector<int> > empty_int;
            if (!Cutwins->Extract_observation_window(0, "GNP_cuboids", Init->geom_sample, Init->geom_sample.sample, Init->nanotube_geo, gnps, structure_cnt, empty_double, empty_point, empty_int, shell_gnps, structure_gnp, points_gnp))
            {
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

    return 1;
}
//This function finds the vertices inside the sample of the GNPs that partially lie outside the sample
int App_Network_From_Displacements::Get_gnps_partially_outside_sample(const Geom_sample& geom_sample, const vector<GNP>& gnps, vector<vector<int> >& vertices_gnps_in)const
{
    //Iterate over the GNPs
    for (size_t i = 0; i < gnps.size(); i++)
    {
        //Vector to store the vertices inside the sample
        vector<int> v_inside;

        //Iterate over the number of GNP vertices
        for (int j = 0; j < 8; j++)
        {

            //Check if vertex j is inside the sample
            if (!gnps[i].vertices[j].is_outside_cuboid(geom_sample.sample))
            {
                //Vertex j is inside the sample, so add it to the vector of
                //vertices inside the sample
                v_inside.push_back(j);
                //cout << "    v=" << j << " inside" << endl;
            }
        }

        //Check if there were vertices outside the sample
        if (v_inside.size() == 8)
        {
            //All vertices are inside the sample, so add an empty vector
            vertices_gnps_in.push_back(vector<int>());
        }
        else
        {
            //Not all vertices are inside the sample,so add the vector v_inside
            vertices_gnps_in.push_back(v_inside);
        }
    }

    return 1;
}
//This function applies displacements which are read from a binary file
int App_Network_From_Displacements::Apply_displacements_from_files(const int& frame, const string& particle_type, const vector<vector<long int> >& structure_cnt, const vector<vector<int> >& vertices_gnps_in, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<GNP>& gnps)const
{
    //Apply displacements to sample
    //hout << "Apply_displacements_to_sample" << endl;
    if (frame >= 2)
    {
        //For frames 2 and onwards apply incremental displacements
        //hout << "Apply_incremental_displacements_to_sample" << endl;
        if (!Apply_incremental_displacements_to_sample(frame, geom_sample))
        {
            hout << "Error in Apply_displacements_from_files when calling Apply_incremental_displacements_to_sample." << endl;
            return 0;
        }
    }
    else 
    {
        //From frame 1 apply full displacement
        //hout << "Apply_displacements_to_sample" << endl;
        if (!Apply_displacements_to_sample(frame, geom_sample))
        {
            hout << "Error in Apply_displacements_from_files when calling Apply_displacements_to_sample." << endl;
            return 0;
        }
    }

    //Apply displacements to CNTs if any
    if (particle_type == "CNT_wires" || particle_type == "CNT_deposit" || particle_type == "GNP_CNT_mix")
    {
        time_t it0 = time(NULL);
        if (!Apply_displacements_to_cnts(frame, structure_cnt, points_cnt))
        {
            hout << "Error in Apply_displacements_from_files when calling Apply_displacements_to_cnts." << endl;
            return 0;
        }

        //If mixed particles, output time for CNTs
        if (particle_type == "GNP_CNT_mix")
        {
            time_t it1 = time(NULL);
            hout << "Apply displacements for CNTs: " << (int)(it1 - it0) << " secs." << endl;
        }
    }

    //Apply displacements to GNPs if any
    if (particle_type == "GNP_cuboids" || particle_type == "GNP_CNT_mix")
    {
        time_t it0 = time(NULL);
        //Apply displacements to GNP vertices
        //hout << "Apply_displacements_to_gnps" << endl;
        if (!Apply_displacements_to_gnps(frame, vertices_gnps_in, gnps))
        {
            hout << "Error in Apply_displacements_from_files when calling Apply_displacements_to_gnps" << endl;
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
//This functions applies displacements to the sample
//Displacements are read from a file
int App_Network_From_Displacements::Apply_displacements_to_sample(const int& frame, Geom_sample& geom_sample)const
{
    //Open the file with sample displacements for the given frame
    ifstream sample_file("Matrix_disp_F" + to_string(frame) + ".dat", ios::in | ios::binary);
    if (!sample_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with matrix displacements (Matrix_disp_F" + to_string(frame) + ".dat)." << endl;
        return 0;
    }

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Variables to store displacements 
    double dx = 0.0, dy = 0.0, dz = 0.0;

    //Read the coordinates for lower left corner (minimum coordinates)
    sample_file.read((char*)&dx, double_size);
    sample_file.read((char*)&dy, double_size);
    sample_file.read((char*)&dz, double_size);

    //Apply displacements
    geom_sample.sample.poi_min.x += dx;
    geom_sample.sample.poi_min.y += dy;
    geom_sample.sample.poi_min.z += dz;

    //Read the coordinates for upper right corner (maxium coordinates)
    sample_file.read((char*)&dx, double_size);
    sample_file.read((char*)&dy, double_size);
    sample_file.read((char*)&dz, double_size);

    //Apply displacements
    geom_sample.sample.max_x += dx;
    geom_sample.sample.max_y += dy;
    geom_sample.sample.max_z += dz;
    
    //Close file as no more data can be extracted form it
    sample_file.close();

    //Update the geom_sample object

    //Update the dimensions of the sample along each direction
    geom_sample.sample.len_x = geom_sample.sample.max_x - geom_sample.sample.poi_min.x;
    geom_sample.sample.wid_y = geom_sample.sample.max_y - geom_sample.sample.poi_min.y;
    geom_sample.sample.hei_z = geom_sample.sample.max_z - geom_sample.sample.poi_min.z;

    //Update the sample's volume
    geom_sample.volume = geom_sample.sample.len_x * geom_sample.sample.wid_y * geom_sample.sample.hei_z;

    return 1;
}
//This functions applies incremental displacements to the sample
//Displacements are read from a file
int App_Network_From_Displacements::Apply_incremental_displacements_to_sample(const int& frame, Geom_sample& geom_sample)const
{
    //Open the file with sample displacements for the current frame
    ifstream current_sample_file("Matrix_disp_F" + to_string(frame) + ".dat", ios::in | ios::binary);
    if (!current_sample_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with matrix displacements (Matrix_disp_F" + to_string(frame) + ".dat)." << endl;
        return 0;
    }

    //Open the file with sample displacements for the previous frame
    ifstream prev_sample_file("Matrix_disp_F" + to_string(frame-1) + ".dat", ios::in | ios::binary);
    if (!prev_sample_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with matrix displacements (Matrix_disp_F" + to_string(frame-1) + ".dat)." << endl;
        return 0;
    }

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Apply displacements to lower left corner (minimum coordinates)
    geom_sample.sample.poi_min.x += Get_difference(double_size, prev_sample_file, current_sample_file);
    geom_sample.sample.poi_min.y += Get_difference(double_size, prev_sample_file, current_sample_file);
    geom_sample.sample.poi_min.z += Get_difference(double_size, prev_sample_file, current_sample_file);

    //Apply displacements to upper right corner (maximum coordinates)
    geom_sample.sample.max_x += Get_difference(double_size, prev_sample_file, current_sample_file);
    geom_sample.sample.max_y += Get_difference(double_size, prev_sample_file, current_sample_file);
    geom_sample.sample.max_z += Get_difference(double_size, prev_sample_file, current_sample_file);

    //Close files as no more data can be extracted form them
    current_sample_file.close();
    prev_sample_file.close();

    //Update the geom_sample object

    //Update the dimensions of the sample along each direction
    geom_sample.sample.len_x = geom_sample.sample.max_x - geom_sample.sample.poi_min.x;
    geom_sample.sample.wid_y = geom_sample.sample.max_y - geom_sample.sample.poi_min.y;
    geom_sample.sample.hei_z = geom_sample.sample.max_z - geom_sample.sample.poi_min.z;

    //Update the sample's volume
    geom_sample.volume = geom_sample.sample.len_x * geom_sample.sample.wid_y * geom_sample.sample.hei_z;

    return 1;
}
//This function gets the diddrence between two quantities in different frames
double App_Network_From_Displacements::Get_difference(const streamsize& double_size, ifstream& previous, ifstream& current)const
{
    //Variables to store values read from file
    double d_curr = 0.0, d_prev = 0.0;

    //Read values
    previous.read((char*)&d_prev, double_size);
    current.read((char*)&d_curr, double_size);
    hout << "d_curr=" << d_curr << " d_prev=" << d_prev << endl;

    //Return the difference
    return (d_curr - d_prev);
}
//This function applies displacements to CNTs
int App_Network_From_Displacements::Apply_displacements_to_cnts(const int& frame, const vector<vector<long int> >& structure_cnt, vector<Point_3D>& points_cnt)const
{
    //Open the file with sample displacements for the current frame
    //hout << "current_file" << endl;
    ifstream current_file("CNTs_disp_F" + to_string(frame) + ".dat", ios::in | ios::binary);
    if (!current_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with CNT displacements (CNTs_disp_F" + to_string(frame) + ".dat) for current frame." << endl;
        return 0;
    }

    //Open the file with sample displacements for the previous frame
    //hout << "prev_file" << endl;
    ifstream prev_file("CNTs_disp_F" + to_string(frame - 1) + ".dat", ios::in | ios::binary);
    //Ignore for frame 1
    if (frame >= 2 && !prev_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with CNT displacements (CNTs_disp_F" + to_string(frame-1) + ".dat) for previous frame." << endl;
        return 0;
    }

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Variables to store displacemants
    double dx = 0.0, dy = 0.0, dz = 0.0;

    //Iterate over all CNTs
    //hout << "Iterate over all CNTs" << endl;
    for (size_t i = 0; i < structure_cnt.size(); i++)
    {
        //Check if end-of-file has been reached
        if (current_file.eof())
        {
            hout << "Error in Apply_displacements_to_cnts. The end-of-file of CNTs_disp_F" + to_string(frame) + ".dat has been reached before reading all CNT data (current frame)." << endl;
            return 0;
        }
        //Check if end-of-file has been reached (ignore for frame 1)
        if (frame >= 2 && prev_file.eof())
        {
            hout << "Error in Apply_displacements_to_cnts. The end-of-file of CNTs_disp_F" + to_string(frame-1) + ".dat has been reached before reading all CNT data (previous frame)." << endl;
            return 0;
        }

        //Iterate over the points of CNT i
        for (size_t j = 0; j < structure_cnt[i].size(); j++)
        {
            //Get the displacements for point j in CNT i
            if (frame >= 2)
            {
                //For frames 2 and onwards get the incremental displacement
                dx = Get_difference(double_size, prev_file, current_file);
                dy = Get_difference(double_size, prev_file, current_file);
                dz = Get_difference(double_size, prev_file, current_file);
            }
            else
            {
                //For frame 1 get the displacement for the frame
                current_file.read((char*)&dx, double_size);
                current_file.read((char*)&dy, double_size);
                current_file.read((char*)&dz, double_size);
            }

            //Get CNT point number
            long int Pj = structure_cnt[i][j];

            //Apply displacements to Pj
            points_cnt[Pj].x += dx;
            points_cnt[Pj].y += dy;
            points_cnt[Pj].z += dz;
        }
    }

    //Close files
    prev_file.close();
    current_file.close();

    return 1;
}
//This function applies displacements to GNPs
int App_Network_From_Displacements::Apply_displacements_to_gnps(const int& frame, const vector<vector<int> >& vertices_gnps_in, vector<GNP>& gnps)const
{
    //Open the file with sample displacements for the current frame
    //hout << "current_file" << endl;
    ifstream current_file("GNPs_disp_F" + to_string(frame) + ".dat", ios::in | ios::binary);
    if (!current_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with GNP displacements (GNPs_disp_F" + to_string(frame) + ".dat) for current frame." << endl;
        return 0;
    }

    //Open the file with sample displacements for the previous frame
    //hout << "prev_file" << endl;
    ifstream prev_file("GNPs_disp_F" + to_string(frame - 1) + ".dat", ios::in | ios::binary);
    //Ignore for frame 1
    if (frame >= 2 && !prev_file) {
        hout << "Error in Apply_displacements_to_sample: Failed to open file with GNP displacements (GNPs_disp_F" + to_string(frame - 1) + ".dat) for previous frame." << endl;
        return 0;
    }

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Create vector with numbers from 0 to 8
    vector<int> all_vertices = All_gnp_vertices();

    //Iterate over all GNPs
    //hout << "Iterate over all GNPs" << endl;
    for (size_t i = 0; i < gnps.size(); i++)
    {
        //Check if end-of-file has been reached
        if (current_file.eof())
        {
            hout << "Error in Apply_displacements_to_cnts. The end-of-file of GNPs_disp_F" + to_string(frame) + ".dat has been reached before reading all GNP data (current frame)." << endl;
            return 0;
        }
        //Check if end-of-file has been reached (ignore for frame 1)
        if (frame >= 2 && prev_file.eof())
        {
            hout << "Error in Apply_displacements_to_cnts. The end-of-file of GNPs_disp_F" + to_string(frame - 1) + ".dat has been reached before reading all GNP data (previous frame)." << endl;
            return 0;
        }

        //Get the vector of vertices
        vector<int> vertices = (vertices_gnps_in[i].empty()) ? all_vertices : vertices_gnps_in[i];
        //hout << "gnp_i=" << i << " vertices.size=" << vertices.size() << endl;

        //Apply displacements to vertices
        //hout << "Apply_displacements_to_vertices" << endl;
        if (!Apply_displacements_to_vertices(frame, double_size, vertices, gnps[i], prev_file,current_file))
        {
            hout << "Error in Apply_displacements_to_gnps when calling Apply_displacements_to_vertices." << endl;
            return 0;
        }
    }

    //Close files
    prev_file.close();
    current_file.close();

    return 1;
}
//This function generates a vector with the numbers from 1 to 8
vector<int> App_Network_From_Displacements::All_gnp_vertices()const
{
    vector<int> tmp(8, 0);

    //Set values 1 to 7
    tmp[1] = 1;
    tmp[2] = 2;
    tmp[3] = 3;
    tmp[4] = 4;
    tmp[5] = 5;
    tmp[6] = 6;
    tmp[7] = 7;

    return tmp;
}
//This function applies the displacements to the vertices of a single GNP
int App_Network_From_Displacements::Apply_displacements_to_vertices(const int& frame, const streamsize& double_size, const vector<int>& vertices, GNP& gnp_i, ifstream& previous, ifstream& current)const
{
    //Variables to store displacemants
    double dx = 0.0, dy = 0.0, dz = 0.0;

    //Variables to accumulate an average
    double avgx = 0.0, avgy = 0.0, avgz = 0.0;

    //Variable to iterate over vertices inside the sample
    int k = 0;

    //Get the number of vertices
    int nv = (int)vertices.size();

    //Vector to store vertices outside the sample
    vector<int> v_out;

    //Iterate over the number of vertices
    for (int j = 0; j < 8; j++)
    {
        //Get the current vertex in the vertices vector
        //hout << "j=" << j << endl;
        int v = vertices[k];
        //hout << "   v=" << v << " k=" << k << endl;
        
        //Check if vertex i is present
        if (v == j)
        {
            //Vertex j is inside the sample
            
            //Get the displacements for vertex v in gnp_i
            if (frame >= 2)
            {
                //For frames 2 and onwards get the incremental displacement
                dx = Get_difference(double_size, previous, current);
                dy = Get_difference(double_size, previous, current);
                dz = Get_difference(double_size, previous, current);
            }
            else
            {
                //For frame 1 get the displacement for the frame
                current.read((char*)&dx, double_size);
                current.read((char*)&dy, double_size);
                current.read((char*)&dz, double_size);
            }
            //hout << "dx=" << dx << " dy=" << dy << " dz=" << dz << endl;

            //Accumulate displacement
            avgx += dx;
            avgy += dy;
            avgz += dz;

            //Apply displacements
            gnp_i.vertices[v].x += dx;
            gnp_i.vertices[v].y += dy;
            gnp_i.vertices[v].z += dz;

            //Get to the next vertex in the next loop
            //If k goes above the number of vertices, reset it back to 0
            k = (k + 1) % nv;
            //hout << "   k_new=" << k << endl;
        }
        else 
        {
            //Vertex is outside the sample so add it to the corresponding vector
            v_out.push_back(j);
        }
    }

    //If there are vertices outside the sample, take the average displacement
    if (v_out.size())
    {
        //Get the number of vertices inside the sample as a double
        double nv = (double)vertices.size();
        avgx = avgx / nv;
        avgy = avgy / nv;
        avgz = avgz / nv;
        //hout << "avgx=" << avgx << " avgy=" << avgy << " avgz=" << avgz << endl;
    }

    //Iterate over the vertices outside the sample and apply average displacements
    //hout << "Apply average displacments" << endl;
    for (size_t j = 0; j < v_out.size(); j++)
    {
        //Get vertex outside the sample
        int v = v_out[j];
        //hout << "v=" << v << endl;

        //Apply average displacements
        gnp_i.vertices[v].x += avgx;
        gnp_i.vertices[v].y += avgy;
        gnp_i.vertices[v].z += avgz;
    }

    return 1;
}