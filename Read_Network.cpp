//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Read 3D nanoparticle network from file
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Read_Network.h"

//This function reads the data from a csv file to generate a nanoparticle network
int Read_Network::Generate_nanoparticle_network_from_file(const Simu_para& simu_para, const Visualization_flags& vis_flags, Geom_sample& geom_sample, vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure, vector<GNP>& gnps)const
{
    //Check a valid particle is to be read from file
    if (simu_para.particle_type == "Hybrid_particles") {
        hout << "Hybrid particles not yet implemented" << endl;
        return 0;
    }
    else if (simu_para.particle_type != "CNT_wires" && simu_para.particle_type != "CNT_deposit" && simu_para.particle_type != "GNP_CNT_mix" && simu_para.particle_type != "GNP_cuboids") {
        hout << "Error in Generate_nanoparticle_network_from_file: the type of particles should be one of the following: CNT_wires, CNT_deposit, GNP_cuboids, or GNP_CNT_mix. Input value was: " << simu_para.particle_type << endl;
        return 0;
    }

    //Read sample geometry
    if (!Read_sample_geometry(geom_sample))
    {
        hout << "Error in Generate_nanoparticle_network_from_file when reading the sample geometry from file" << endl;
        return 0;
    }

    //Check if there are CNTs
    if (simu_para.particle_type == "CNT_wires" || simu_para.particle_type == "CNT_deposit" || simu_para.particle_type == "GNP_CNT_mix") {

        //Check the type of file to read
        if (simu_para.file_type == "csv")
        {
            //Read CNT data from csv file
            if (!Read_cnt_data_from_csv(points_cnt, radii, structure))
            {
                hout << "Error in Generate_nanoparticle_network_from_file when calling Read_cnt_data_from_csv." << endl;
                return 0;
            }
        }
        else if (simu_para.file_type == "dat")
        {
            //Read CNT data from binary file (.dat)
            if (!Read_cnt_data_from_dat(points_cnt, radii, structure))
            {
                hout << "Error in Generate_nanoparticle_network_from_file when calling Read_cnt_data_from_dat." << endl;
                return 0;
            }
        }
        else {
            hout << "Error in Generate_nanoparticle_network_from_file. Invalid file type to read CNT network from. Valid options are 'csv' or 'dat' only. Input was:" << simu_para.file_type << endl;
            return 0;
        }

        /* /Approximate CNT volume fraction
        double cnt_vol = 0.0;
        for (size_t i = 0; i < structure.size(); i++)
        {
            //Length of CNT i
            double l_cnt = 0.0;
            //Point P0
            long int P0 = structure[i][0];
            for (size_t j = 1; j < structure[i].size(); j++)
            {
                //Point P1
                long int P1 = structure[i][j];
                //Accumulate distance from point j-i (P0) to point j (P1)
                l_cnt = l_cnt + points_cnt[P0].distance_to(points_cnt[P1]);
                //Update P0 for the next iteration over j
                P0 = P1;
            }
            //Accumulate CNT volume
            cnt_vol = cnt_vol + PI * radii[i] * radii[i] * l_cnt;
        }
        //Output approximated volume fraction
        hout << "vol frac=" << cnt_vol / geom_sample.volume << endl;*/
    }
    
    //Check if there are GNPs
    if (simu_para.particle_type == "GNP_cuboids" || simu_para.particle_type == "GNP_CNT_mix") {

        //Check the type of file to read
        if (simu_para.file_type == "csv")
        {
            //Read the GNP geometry form a csv file
            if (!Read_gnp_data_from_csv(geom_sample.sample, gnps))
            {
                hout << "Error in Generate_nanoparticle_network_from_file when calling Read_gnp_data_from_csv." << endl;
                return 0;
            }
        }
        else if (simu_para.file_type == "dat")
        {
            //Read the GNP geometry form a binary file (.dat)
            if (!Read_gnp_data_from_dat(geom_sample.sample, gnps))
            {
                hout << "Error in Generate_nanoparticle_network_from_file when calling Read_gnp_data_from_dat." << endl;
                return 0;
            }
        }
        else {
            hout << "Error in Generate_nanoparticle_network_from_file. Invalid file type to read GNP network from. Valid options are 'csv' or 'dat' only. Input was:" << simu_para.file_type << endl;
            return 0;
        }
    }

    //---------------------------------------------------------------------------
    //Check if visualization files were requested for read nanoparticles
    if (vis_flags.generated_nanoparticles) {

        VTK_Export VTK;

        //Export generated CNTs if any
        if (points_cnt.size()) {
            if (vis_flags.generated_nanoparticles == 2) {

                //Generate one visualization file per CNT
                for (size_t i = 0; i < structure.size(); i++) {
                    string filename = "cnt_" + to_string(i) + ".vtk";
                    VTK.Export_from_cnt_structure(points_cnt, vector<vector<long int> > (1, structure[i]), filename);
                }
            }
            else if (vis_flags.generated_nanoparticles == 1 || vis_flags.generated_nanoparticles == 3) {

                //Generate one visualization file with all the CNTs
                VTK.Export_from_cnt_structure(points_cnt, structure, "cnts_read.vtk");
            }
        }

        //Export generated GNPs if any
        if (gnps.size()) {

            if (vis_flags.generated_nanoparticles == 1)
            {
                //Generate one visualization file with all the GNPs
                VTK.Export_gnps(gnps, "gnps_read.vtk");
            }
            else if (vis_flags.generated_nanoparticles == 2)
            {
                //Generate one visualization file per GNP
                VTK.Export_gnps_single_files(gnps, "gnps_read");
            }
        }

        //Export the sample geometry
        VTK.Export_cuboid(geom_sample.sample, "sample.vtk");
    }

    return 1;
}
//This function reads the sample geomtry from a csv file (sample_geom.csv)
int Read_Network::Read_sample_geometry(Geom_sample& geom_sample)const
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
//This function reads the CNT data from csv files (cnt_struct.csv and cnt_coordinates.csv) 
int Read_Network::Read_cnt_data_from_csv(vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure)const
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

        //Iterate over the remaining points in CNT i to fill the structure[i] vector
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

    //Variable to count the number of points (again)
    long int Nq = 0;

    //Read the point coordinates
    //In order to easily update the CNT number, iterate over the structure
    for (size_t i = 0; i < structure.size(); i++) {
        for (size_t j = 0; j < structure[i].size(); j++) {

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
            ss >> points_cnt[Nq].x;
            ss.ignore();
            ss >> points_cnt[Nq].y;
            ss.ignore();
            ss >> points_cnt[Nq].z;

            //Set the CNT flag with the CNT number (i is the CNT number)
            points_cnt[Nq].flag = (int)i;

            //Increase the count of points
            Nq++;
        }
    }

    //Check that the two counts of the number of points yield the same result
    if (Np != Nq)
    {
        hout << "Error in Read_cnt_data_from_csv. The two counts of the number of points is different. Np=" << Np << " Nq=" << Nq << endl;
        return 0;
    }

    //Close files
    struc_file.close();
    coord_file.close();

    //Output a message with the number of CNTs and points read
    hout << endl << "A total of " << n_cnts << " CNTs and " << Np << " points were read." << endl;

    return 1;
}
//This function reads the CNT data from binary files (cnt_struct.dat and cnt_coordinates.dat) 
int Read_Network::Read_cnt_data_from_dat(vector<Point_3D>& points_cnt, vector<double>& radii, vector<vector<long int> >& structure)const
{
    //Open the structure file
    ifstream struc_file("cnt_struct.dat", ios::in | ios::binary);
    if (!struc_file) {
        hout << "Failed to open CNT structure file cnt_struct.dat." << endl;
        return 0;
    }
    //Get the size of a double
    streamsize double_size = sizeof(double);
    //Get the size of an int
    streamsize int_size = sizeof(int);

    //Integer to store the number of CNTs
    int n_cnts = 0;

    //Read the number of CNTs, which is the first entry in the file
    struc_file.read((char*)&n_cnts, int_size);

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
            hout << "Error in Read_cnt_data_from_dat. The end-of-file of cnt_struct.dat has been reached before reading all CNT data." << endl;
            return 0;
        }

        //Variable to store the number of points in CNT i
        int np;

        //Read the number of CNT points from the file
        struc_file.read((char*)&np, int_size);

        //Read the radius from the file
        struc_file.read((char*)&radii[i], double_size);

        //Fill structure vector for CNT i
        structure[i].assign(np, Np);
        //Increase the number of points since the current values is already used 
        //for the first point of CNT i
        Np++;

        //Iterate over the remaining points in CNT i to fill the structure[i] vector
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
    ifstream  coord_file("cnt_coordinates.dat", ios::in | ios::binary);
    if (!coord_file) {
        hout << "Failed to open point cordinates file cnt_coordinates.dat." << endl;
        return 0;
    }

    //Set the size for the points vector
    points_cnt.assign(Np, Point_3D());

    //Variable to count the number of points (again)
    long int Nq = 0;

    //Read the point coordinates
    //In order to easily update the CNT number, iterate over the structure
    for (size_t i = 0; i < structure.size(); i++) {
        for (size_t j = 0; j < structure[i].size(); j++) {

            //Check if end-of-file has been reached
            if (coord_file.eof())
            {
                hout << "Error in Read_cnt_data_from_dat. The end-of-file of cnt_coordinates.dat has been reached before reading all CNT data." << endl;
                return 0;
            }

            //Read the point coordinates from the file
            coord_file.read((char*)&points_cnt[Nq].x, double_size);
            coord_file.read((char*)&points_cnt[Nq].y, double_size);
            coord_file.read((char*)&points_cnt[Nq].z, double_size);

            //Set the CNT flag with the CNT number (i is the CNT number)
            points_cnt[Nq].flag = (int)i;

            //Increase the count of points
            Nq++;
        }
    }

    //Check that the two counts of the number of points yield the same result
    if (Np != Nq)
    {
        hout << "Error in Read_cnt_data_from_csv. The two counts of the number of points is different. Np=" << Np << " Nq=" << Nq << endl;
        return 0;
    }

    //Close files
    struc_file.close();
    coord_file.close();

    //Output a message with the number of CNTs and points read
    hout << endl << "A total of " << n_cnts << " CNTs and " << Np << " points were read." << endl;

    return 1;
}
//This function reads GNP data from a csv file (gnp_data.csv)
int Read_Network::Read_gnp_data_from_csv(const cuboid& sample_geom, vector<GNP>& gnps)const
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

    //Variable to count the GNPs
    int n_gnps = 0;

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

        //Set the GNP flag
        new_gnp.flag = n_gnps;

        //Increase the count of GNPS
        n_gnps++;

        //Add GNP to the vector of GNPs
        gnps.push_back(new_gnp);
    }

    //Close file
    gnp_file.close();

    //Output a message with the total number of GNPs read
    hout << endl << "A total of " << gnps.size() << " GNPs were read." << endl;

    return 1;
}
//This function reads GNP data from a binary file (gnp_data.dat)
int Read_Network::Read_gnp_data_from_dat(const cuboid& sample_geom, vector<GNP>& gnps)const
{
    //Open the file with the GNP geometric data
    ifstream gnp_file;
    gnp_file.open("gnp_data.dat", ios::in | ios::binary);
    if (!gnp_file) {
        hout << "Failed to open file with GNP geometric data gnp_data.dat." << endl;
        return 0;
    }

    //Generate_Network object to generate the variables needed for each GNP
    Generate_Network GN;

    //Variable to store the number of points in CNT i
    int n_gnps;

    //Read the number of GNPs
    gnp_file.read((char*)&n_gnps, sizeof(int));

    //Set the size of the GNP vector
    gnps.assign(n_gnps, GNP());

    //Get the size of a double
    streamsize double_size = sizeof(double);

    //Iterate over the number of GNPs to read the data for each GNP
    for (int i = 0; i < n_gnps; i++)
    {
        //Check if end-of-file has been reached
        if (gnp_file.eof())
        {
            hout << "Error in Read_gnp_data_from_dat. The end-of-file of gnp_data.dat has been reached before reading all GNP data." << endl;
            return 0;
        }

        //Set the flag of current GNP
        gnps[i].flag = i;

        //Read the GNP geometry
        gnp_file.read((char*)&gnps[i].l, double_size);
        gnp_file.read((char*)&gnps[i].t, double_size);

        //Variables to store the angles
        double theta, phi;

        //Read the rotation angles
        gnp_file.read((char*)&theta, double_size);
        gnp_file.read((char*)&phi, double_size);

        //Calculate rotation matrix
        gnps[i].rotation = GN.Get_transformation_matrix(theta, phi);

        //Read the coordinates of the GNP's centroid
        gnp_file.read((char*)&gnps[i].center.x, double_size);
        gnp_file.read((char*)&gnps[i].center.y, double_size);
        gnp_file.read((char*)&gnps[i].center.z, double_size);

        //Calculate (or approximate) GNP volume
        int j;
        for (j = 0; j < 8; j++)
        {
            if (!GN.Is_point_inside_cuboid(sample_geom, gnps[i].vertices[j])) {

                //GNP is partially outside the sample, so approximate volume
                double gnp_vol = 0.0;
                if (!GN.Approximate_gnp_volume_inside_sample(sample_geom, gnps[i], gnp_vol)) {
                    hout << "Error in Read_gnp_data_from_dat when calling Approximate_gnp_volume_inside_sample." << endl;
                    return 0;
                }
                gnps[i].volume = gnp_vol;

                //Terminate the for loop
                break;
            }
        }

        //Check if all vertices were inside the sample
        if (j == 7)
        {
            //All vertices were inside the sample, thus calcualte the volume of the whole GNP
            gnps[i].volume = gnps[i].l * gnps[i].l * gnps[i].t;
        }

        //Obtain coordinates of GNP vertices
        if (!GN.Obtain_gnp_vertex_coordinates(gnps[i])) {
            hout << "Error in Read_gnp_data_from_dat when calling Obtain_gnp_vertex_coordinates." << endl;
            return 0;
        }

        //Get the plane equations for the six faces
        if (!GN.Update_gnp_plane_equations(gnps[i])) {
            hout << "Error in Read_gnp_data_from_dat when calling Update_gnp_plane_equations." << endl;
            return 0;
        }
    }

    //Close file
    gnp_file.close();

    //Output a message with the total number of GNPs read
    hout << endl << "A total of " << gnps.size() << " GNPs were read." << endl;

    return 1;
}