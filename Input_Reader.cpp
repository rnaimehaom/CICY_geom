//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Read input data and save it into data structures
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_input_file(ifstream &infile)
{

	hout << "Reading input file..." << endl;
    
    //Read line by line until reaching the end-of-file
    while(!infile.eof())
    {
        istringstream istr(Get_Line(infile));
        if(infile.eof()) break;
        string str_temp;
        istr >> str_temp;
        
        //Skip over empty lines
        if(!str_temp.empty()) {
            if(str_temp=="Application") {
                if(!Read_application(app_name, infile)) return 0;
            }
            else if(str_temp=="Simulation_Parameters") {
                if(!Read_simulation_parameters(simu_para, infile)) return 0;
            }
            else if(str_temp=="Sample_Geometry") {
                if(!Read_sample_geometry(geom_sample, infile)) return 0;
            }
            else if(str_temp=="Cutoff_Distances") {
                if(!Read_cutoff_distances(cutoff_dist, infile)) return 0;
            }
            else if(str_temp=="CNT_Geometry") {
                if(!Read_nanotube_geo_parameters(nanotube_geo, infile)) return 0;
            }
            else if(str_temp=="GNP_Geometry") {
                if(!Read_gnp_geo_parameters(gnp_geo, infile)) return 0;
            }
            else if(str_temp=="Electrical_Parameters") {
                if(!Read_electrical_parameters(electric_para, infile)) return 0;
            }
            else if(str_temp=="Visualization_Flags") {
                if(!Read_visualization_flags(vis_flags, infile)) return 0;
            }
            else if(str_temp=="Output_Data_Flags") {
                if(!Read_output_data_flags(out_flags, infile)) return 0;
            }
            else
            {
                hout << "Error in Read_input_file: Keyword \"" << str_temp << "\" is not defined." << endl;
                return 0;
            }
        }
    }

	hout << "Finished reading input file!" << endl;
    
    //Check if warning messages are needed
    if(!app_name.mark) {
        Warning_message("Application");
    }
    if(!simu_para.mark) {
        Warning_message("Simulation_Parameters");
    }
    if(!geom_sample.mark) {
        Warning_message("Sample_Geometry");
    }
    if(!cutoff_dist.mark) {
        Warning_message("Cutoff_Distances");
    }
    if(!nanotube_geo.mark && simu_para.particle_type != "GNP_cuboids") {
        Warning_message("CNT_Geometry");
    }
    if(!gnp_geo.mark && simu_para.particle_type != "CNT_wires") {
        Warning_message("GNP_Geometry");
    }
    if(!electric_para.mark) {
        Warning_message("Electrical_Parameters");
    }
    if(!vis_flags.mark) {
        Warning_message("Visualization_Flags");
    }
    if(!out_flags.mark) {
        Warning_message("Output_Data_Flags");
    }
    
    /*
    if (simu_para.particle_type=="CNT_wires"||simu_para.particle_type=="GNP_cuboids") {
        hout<<"CNT volume="<<nanotube_geo.volume<<" vol frac="<<nanotube_geo.volume/geom_sample.volume<<" weight="<<nanotube_geo.weight<<" wt frac="<<nanotube_geo.weight/((geom_sample.volume-nanotube_geo.volume)*geom_sample.matrix_density + nanotube_geo.weight)<<endl;
        hout<<"GNP volume="<<gnp_geo.volume<<" vol frac="<<gnp_geo.volume/geom_sample.volume<<" weight="<<gnp_geo.weight<<" wt frac="<<gnp_geo.weight/((geom_sample.volume-gnp_geo.volume)*geom_sample.matrix_density + gnp_geo.weight)<<endl;
    }
    else {
        double matrix_weight = (geom_sample.volume-nanotube_geo.volume-gnp_geo.volume)*geom_sample.matrix_density;
        double den = matrix_weight + nanotube_geo.weight + gnp_geo.weight;
        hout<<"CNT volume="<<nanotube_geo.volume<<" vol frac="<<nanotube_geo.volume/geom_sample.volume<<" weight="<<nanotube_geo.weight<<" wt frac="<<nanotube_geo.weight/den<<endl;
        hout<<"GNP volume="<<gnp_geo.volume<<" vol frac="<<gnp_geo.volume/geom_sample.volume<<" weight="<<gnp_geo.weight<<" wt frac="<<gnp_geo.weight/den<<endl;
    }*/
    
    return 1;
}
//---------------------------------------------------------------------------
//Initialize data
int Input::Data_Initialization()
{
	//Initialize name of simulation
	app_name.keywords = "Application";
	app_name.mark = false;
	app_name.str = "3D_Electrical_Network";

	//Initialize paramters of simulation
	simu_para.keywords = "Simulation_Parameters";
	simu_para.mark = false;
	simu_para.simu_name = "Test";
	simu_para.sample_num = 1;
	simu_para.create_read_network = "Create_Network";
    simu_para.simulation_scope = 0;
    simu_para.resistances[0] = 1;
    simu_para.resistances[1] = 1;
    simu_para.resistances[2] = 1;
    simu_para.tolerance = 1e-10;
    simu_para.MAX_ATTEMPTS_CNT = 5;
    simu_para.MAX_ATTEMPTS_GNP = 1;
    
	//Initialize the geometric parameters of the sample
	geom_sample.keywords = "Sample_Geometry";
	geom_sample.mark = false;
    geom_sample.sample.poi_min.x = 0.0;
	geom_sample.sample.poi_min.y = 0.0;
	geom_sample.sample.poi_min.z = 0.0;
	geom_sample.sample.poi_min.flag = 0;
	geom_sample.sample.len_x = 1.0;
	geom_sample.sample.wid_y = 1.0;
	geom_sample.sample.hei_z = 1.0;
    geom_sample.ex_dom_cnt.poi_min.x = 0.0;
	geom_sample.ex_dom_cnt.poi_min.y = 0.0;
	geom_sample.ex_dom_cnt.poi_min.z = 0.0;
	geom_sample.ex_dom_cnt.poi_min.flag = 0;
    geom_sample.ex_dom_cnt.len_x = 1.0;
    geom_sample.ex_dom_cnt.wid_y = 1.0;
    geom_sample.ex_dom_cnt.hei_z = 1.0;
    geom_sample.ex_dom_gnp.poi_min.x = 0.0;
    geom_sample.ex_dom_gnp.poi_min.y = 0.0;
    geom_sample.ex_dom_gnp.poi_min.z = 0.0;
    geom_sample.ex_dom_gnp.poi_min.flag = 0;
    geom_sample.ex_dom_gnp.len_x = 1.0;
    geom_sample.ex_dom_gnp.wid_y = 1.0;
    geom_sample.ex_dom_gnp.hei_z = 1.0;
	geom_sample.volume = geom_sample.sample.len_x*geom_sample.sample.wid_y*geom_sample.sample.hei_z;
	geom_sample.matrix_density = 1.0;
	geom_sample.gs_minx = 1.0;
	geom_sample.gs_miny = 1.0;
	geom_sample.gs_minz = 1.0;
	geom_sample.win_max_x = 1.0;
	geom_sample.win_max_y = 1.0;
	geom_sample.win_max_z = 1.0;
	geom_sample.win_min_x = 1.0;
	geom_sample.win_min_y = 1.0;
	geom_sample.win_min_z = 1.0;
	geom_sample.win_delt_x = 1.0;
	geom_sample.win_delt_y = 1.0;
	geom_sample.win_delt_z = 1.0;

	//Initialize the geometric paramters of nanotubes
	nanotube_geo.keywords = "CNT_Geometry";
	nanotube_geo.mark = false;
	nanotube_geo.dir_distrib_type = "random";
	nanotube_geo.ini_phi = 0.0;
	nanotube_geo.ini_theta = 0.0;
	nanotube_geo.omega_b = 1.5707963267948966;
    nanotube_geo.omega_a = - nanotube_geo.omega_b;
	nanotube_geo.step_length = 0.01;
	nanotube_geo.len_distrib_type = "uniform";
	nanotube_geo.len_min = 1.0;
	nanotube_geo.len_max = 1.0;
	nanotube_geo.rad_distrib_type = "uniform";
	nanotube_geo.rad_min = 0.005;
	nanotube_geo.rad_max = 0.005;
	nanotube_geo.criterion = "vol";
	nanotube_geo.volume_fraction = 0.0;
	nanotube_geo.volume = 0.0;
	nanotube_geo.weight_fraction = 0.0;
	nanotube_geo.weight = 0.0;
    
    //Initialize the geometric paramters of GNPs
    gnp_geo.keywords = "GNP_Geometry";
    gnp_geo.mark = false;
    gnp_geo.criterion = "vol";
    gnp_geo.growth_type = "independent";
    gnp_geo.orient_distrib_type = "random";
    gnp_geo.size_distrib_type = "uniform";
    gnp_geo.thick_distrib_type = "uniform";
    gnp_geo.len_min = 1.0;
    gnp_geo.len_max = 1.0;
    gnp_geo.t_min = 0.03;
    gnp_geo.t_max = 0.03;
    gnp_geo.mass_ratio = 1.0;
    gnp_geo.volume_fraction = 0.0;
    gnp_geo.volume = 0.0;
    gnp_geo.weight_fraction = 0.0;
    gnp_geo.weight = 0.0;
    gnp_geo.density = 2.25;

	//Initialize cutoff distances
	cutoff_dist.keywords = "Cutoff_Distances";
	cutoff_dist.mark = false;
	cutoff_dist.tunneling_dist = 0.0018;
	cutoff_dist.van_der_Waals_dist = 0.00034;
    cutoff_dist.min_points = 5;

	//Initialize electrical parameters
	electric_para.keywords = "Electrical_Parameters";
	electric_para.mark = false;
	electric_para.applied_voltage = 1.0;
	electric_para.resistivity_CNT = 0.001;
    
    //Initialize visualization flags (do not export anything)
    vis_flags.keywords = "Visualization_Flags";
    vis_flags.mark = false;
    vis_flags.generated_nanoparticles = 0;
    vis_flags.clusters = 0;
    vis_flags.percolated_clusters = 0;
    vis_flags.backbone = 0;
    vis_flags.triangulations = 0;
    vis_flags.sample_domain = 0;
    vis_flags.window_domain = 0;
    
    //Initialize output data flags (do not output any data)
    out_flags.keywords = "Output_Data_Flags";
    out_flags.mark = false;
    out_flags.cnt_gnp_flag = 0;
    out_flags.gnp_data = 0;
    out_flags.cnt_data = 0;

	hout << "    Data initialization done" <<endl<<endl;

	return 1;
}
//---------------------------------------------------------------------------
//Error message for a keywork not found
void Input::Warning_message(const string &str)
{
    hout<<"Warning: Keyword \""<<str<<"\" was not found. Deafault parameters will be used."<< endl;
}
//---------------------------------------------------------------------------
//Error message for a keywork already input
void Input::Warning_message_already_input(const string &str)
{
    hout << "Warning: Keyword \"" << str << "\" has already been input." << endl;
}
//---------------------------------------------------------------------------
//Read the application type
int Input::Read_application(App_name &app_name, ifstream &infile)
{
	if(app_name.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(app_name.keywords);
        return 0;
	}
	else app_name.mark = true;

	istringstream istr(Get_Line(infile));
	istr >> app_name.str;

	return 1;
}
//---------------------------------------------------------------------------
//Read the simulation parameters
int Input::Read_simulation_parameters(Simu_para &simu_para, ifstream &infile)
{
	if(simu_para.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(simu_para.keywords);
		return 0;
	}
	else simu_para.mark = true;
    
    //Read the name of simulation
    istringstream istr0(Get_Line(infile));
	istr0 >> simu_para.simu_name;
    
    //Read the number of samples
	istringstream istr1(Get_Line(infile));
	istr1 >> simu_para.sample_num;
	if(simu_para.sample_num<1)	 {
        hout << "Error: the number of samples cannot be less than 1. Input was: "<<simu_para.sample_num<< endl; return 0; }
    
    //----------------------------------------------------------------------
    //Define the type of particles inside the sample:
    //CNT_wires for CNTs
    //GNP_cuboids for GNPs
    //Hybrid_particles for GNP-CNT hybrid particles
    //GNP_CNT_mix for mixed CNT and GNP without hybridizing
    istringstream istr_perticle_type(Get_Line(infile));
    istr_perticle_type >> simu_para.particle_type;
    if (simu_para.particle_type == "Hybrid_particles") {
        hout << "Hybrid particles are not implemented yet" << endl;
        return 0;
    }
    else if (simu_para.particle_type != "CNT_wires" && 
        simu_para.particle_type != "CNT_deposit" && 
        simu_para.particle_type != "GNP_cuboids" && 
        simu_para.particle_type != "GNP_CNT_mix") {
        hout << "Error: the type of particles shoud be one of the following: CNT_wires, CNT_deposit, GNP_cuboids or GNP_CNT_mix. Input was: "<<simu_para.particle_type<< endl;
        return 0;
    }
    
    //Read keyword for creating a new network or reading a network from a file
	istringstream istr2(Get_Line(infile));
	istr2 >> simu_para.create_read_network;
    if(simu_para.create_read_network!="Create_Network" && 
        simu_para.create_read_network != "Read_Network" && 
        simu_para.create_read_network != "Read_Network_Abq" && 
        simu_para.create_read_network!="Network_From_Seeds")
    {
        hout << "Error: Invalid keyword. Valid options are 'Create_Network', 'Read_Network', 'Read_Network_Abq', or 'Network_From_Seeds'. Input was: " << simu_para.create_read_network << endl; 
        return 0; 
    }
    
    //----------------------------------------------------------------------
    //If the network is generated from some seeds, then read them
    if (simu_para.create_read_network=="Network_From_Seeds") {
        
        //Different types of fillers need different number of seeds
        //Thus, check the filler type and read the necessary number of seeds
        if (simu_para.particle_type!="GNP_cuboids") {
            //If the particle type is not GNP_cuboids, then CNTs are generated for sure
            
            //7 seeds are required for CNTs
            int n_seeds = 7;
            simu_para.CNT_seeds.assign(n_seeds, 0);
            
            //Read the line with the seeds
            istringstream istr_network_seeds(Get_Line(infile));
            //Add seeds to the vector
            for (int i = 0; i < n_seeds; i++) {
                istr_network_seeds >> simu_para.CNT_seeds[i];
                //hout<<"seed["<<i<<"]="<<simu_para.CNT_seeds[i]<<endl;
            }
        }
        if (simu_para.particle_type!="CNT_wires" && simu_para.particle_type!="CNT_deposit") {
            //If the particle type is not CNT_wires nor CNT_deposit, then GNPs are generated for sure
            
            //6 seeds are required for GNPs
            int n_seeds = 6;
            simu_para.GNP_seeds.assign(n_seeds, 0);
            
            //Read the line with the seeds
            istringstream istr_network_seeds(Get_Line(infile));
            //Add seeds to the vector
            for (int i = 0; i < n_seeds; i++) {
                istr_network_seeds >> simu_para.GNP_seeds[i];
                //hout << simu_para.GNP_seeds[i] << endl;
            }
        }
    }
    else if (simu_para.create_read_network == "Read_Network_Abq")
    {
        //The network is read from files and displacements from an Abaqus database

        //Read the path to the odb file (the database)
        istringstream istr_odb_file(Get_Line(infile));
        istr_odb_file >> simu_para.odb_file;

        //Read the name of the step in the Abaqus simulation
        istringstream istr_step(Get_Line(infile));
        istr_step >> simu_para.step_name;

        //When displacements are read from Abaqus, the network is read from a csv file
        //since this is the file that was used to generate the network in Abaqus
        simu_para.file_type = "csv";
    }
    else if (simu_para.create_read_network == "Read_Network")
    {
        //The network is read from files

        //Read the type of file to be read
        istringstream istr_file_type(Get_Line(infile));
        istr_file_type >> simu_para.file_type;

        if (simu_para.file_type != "csv" && simu_para.file_type != "dat")
        {
            hout << "Error when reading simulation parameters: file type for reading the network is invalid. Only options are 'csv' and 'dat'. Input was: " << simu_para.file_type << endl;
            return 0;
        }
    }

    //Read the nanoparticle content
    // 
    //Read the content
    istringstream istr_content(Get_Line(infile));
    istr_content >> simu_para.criterion;
    //
    if (simu_para.criterion == "vol") {
        istr_content >> simu_para.volume_fraction;

        if (simu_para.volume_fraction > 1 || simu_para.volume_fraction < Zero) {
            hout << "Error: The volume fraction must be between 0 and 1. Input was: " << simu_para.volume_fraction << endl; return 0;
        }
    }
    else if (simu_para.criterion == "wt")
    {
        istr_content >> simu_para.weight_fraction;

        if (simu_para.weight_fraction > 1 || simu_para.weight_fraction < Zero) {
            hout << "Error: The weight fraction must be between 0 and 1. Input was:" << simu_para.weight_fraction << endl; return 0;
        }
    }
    else {
        hout << "Error: The content of nanoparticles can only be specified in volume (vol) or weight (wt) fraction. Input was: " << simu_para.criterion << endl; return 0;
    }

    //If mixed or hybrid particles were selected, then check if mass ratio or CNT density on GNPs
    //is specified
    if (simu_para.particle_type == "Hybrid_particles" || simu_para.particle_type == "GNP_CNT_mix") {
        istringstream istr_mixed(Get_Line(infile));
        istr_mixed >> simu_para.mixed;

        //Check if mass ratio is defined
        if (simu_para.mixed == "mass_ratio") {

            /*if (simu_para.criterion != "wt") {
                hout << "Error: If nanoparticle content is specified in weight fraction, then CNT/GNP mass ratio must specified, not volume faction"<< endl; return 0;
            }*/

            //Read the mass ratio
            istr_mixed >> simu_para.mass_ratio;
            if (simu_para.mass_ratio < Zero) {
                hout << "Error: The CNT/GNP mass ratio must be positive. Input was: " << simu_para.mass_ratio << endl; return 0;
            }
        }
        //Check if volume ratio is defined
        else if (simu_para.mixed == "volume_ratio") {

            if (simu_para.criterion != "vol") {
                hout << "Error: If nanoparticle content is specified in volume fraction, then CNT/GNP volume ratio must specified, not weight faction" << endl; return 0;
            }

            //Read the volume ratio
            istr_mixed >> simu_para.volume_ratio;
            if (simu_para.volume_ratio < Zero) {
                hout << "Error: The CNT/GNP volume ratio must be positive. Input was: " << simu_para.volume_ratio << endl; return 0;
            }
        }
        //Check if CNT density on GNPs is defined
        else if (simu_para.mixed == "density") {

            //This option is only valid for hybrid particles, so double check the particle type
            if (simu_para.particle_type != "Hybrid_particles") {
                hout << "Error: The CNT density on GNPs must can only specified for hybrid particles. It was specified for mixed particles." << endl; return 0;
            }
            //This option is only valid if content is specified as volume fraction
            if (simu_para.criterion != "vol") {
                hout << "Error: If CNT density on GNPs is specified, then the nanoparticle content must be specified in volume fraction." << endl; return 0;
            }

            //Read the CNT density
            istr_mixed >> simu_para.cnt_gnp_density;
            if (simu_para.cnt_gnp_density < Zero) {
                hout << "Error: The CNT density on GNPs must be positive. Input was: " << simu_para.cnt_gnp_density << endl; return 0;
            }
        }
        else {
            hout << "Error: The relative content of nanoparticles can only be specified in mass ratio (mass_ratio), volume ratio (volume_ratio), or CNT density on GNPs (density). Input was: " << simu_para.mixed << endl; return 0;
        }
    }
    
    //Flag that defines whether the non-penetrating or penetrating model is used
    istringstream istr_pm_flag(Get_Line(infile));
    istr_pm_flag >> simu_para.penetration_model_flag;
    if (simu_para.penetration_model_flag > 1 || simu_para.penetration_model_flag < 0) {
        hout << "Error: Invalid value for penetration model flag. Valid options are 0 (penetrating model) or 1 (non-penetrating model). Input was: " << simu_para.penetration_model_flag << endl; return 0;
    }
    //If CNT_depost is chose, override the flag
    if (simu_para.particle_type == "CNT_deposit") {
        simu_para.penetration_model_flag = 1;
    }
    
    //Flag to define the scope of the simulation
    // 0: Full simulations, i.e., generate a nanoparticle network and calculate its resistance
    //    When set to 0, three more flags are needed to calculate the resistance on each direction(x, y, z)
    //    One flag per direction, 1 for calculating resistance in that direction or 0 for not calculating it
    //    Example: 0 1 0 1 means calculate resistance of the network(first 0) along directions x and z 
    //    (first and second 1) and do not calculate the resistance along y(second 0)
    //    Example : 0 0 0 0 is equivalent to inputting only "1"
    // 1 : Avoid calculating the resistance of the network but calculate the fractions of 
    //     nanoparticles that belong to the backbone
    // 2 : Generate the network and only determine if there is percolation or not
    //     This is sent as a message in the output file
    // 3 : Only generate the nanoparticle network
    //     This is useful when a network needs to be generated and then exported into Abaqus
    istringstream istr3(Get_Line(infile));
    istr3 >> simu_para.simulation_scope;
    if (simu_para.simulation_scope < 0 || simu_para.simulation_scope > 3) {
        hout << "Error: Invalid value for network resistance calculation. Valid options are between 0 and 3. Input was: " << simu_para.simulation_scope << endl; 
        return 0;
    }
    //If only the network is to be generated and the network is read from file, then terminate
    //with an error as this case is not useful
    if (simu_para.simulation_scope == 3 && simu_para.create_read_network == "Read_Network") {
        hout << "Error: Cannot choose 3 (only generate network) as the flag for the simulation scope when the nanoparticle network is read from file." << endl;
        return 0;
    }
    //If only the network is to be generated and the network is read from file, then terminate
    //with an error as this case is not useful
    if ((simu_para.simulation_scope == 3 || simu_para.simulation_scope == 1)&& app_name.str == "Network_From_Abaqus") {
        hout << "Error: Cannot choose 1 or 3 as the flag for the simulation scope when the application is Network_From_Abaqus." << endl;
        return 0;
    }
    //Check if the network resistance needs to be calculated
    //When choosing 0, three more flags are needed to calculate the resistance on each direction (x, y, z)
    //One flag per direction, 1 for calculating resistance in that direction or 0 for not calculating it
    //Example: 0 1 0 1 means calculate resistance of the network (first 0) along directions x and z (first and
    //second 1) and do not calculate the resistance along y (second 0)
    //Example: 0 0 0 0 is equivalent to inputting only "1"
    if (!simu_para.simulation_scope) {
        //Read the flags for each direction
        istr3 >> simu_para.resistances[0];
        istr3 >> simu_para.resistances[1];
        istr3 >> simu_para.resistances[2];
        
        //Check the case in which the three flags are zero
        if (!(simu_para.resistances[0] + simu_para.resistances[1] + simu_para.resistances[2])) {
            
            //This is the same case as setting the simulation_scope flag to 1,
            //Thus, the simulation_scope is set to 1 since that case is simpler
            simu_para.simulation_scope = 1;
        }
    }
    
    //Tolerance for converge in the conjugate gradient algorithm
    //This value represents the reduction in the initial residual
    istringstream istr_tol(Get_Line(infile));
    istr_tol >> simu_para.tolerance;

    //Maximum number of iterations:
    //for generating a CNT point
    //for relocating a GNP
    istringstream istr_att(Get_Line(infile));
    //Check the particle type
    if (simu_para.particle_type == "CNT_wires" ||
        simu_para.particle_type == "CNT_deposit" ||
        simu_para.particle_type == "GNP_CNT_mix")
    {
        //There are CNTs, so read the maximum number of attempts for generating a CNT point
        istr_att >> simu_para.MAX_ATTEMPTS_CNT;
    }
    if (simu_para.particle_type == "GNP_cuboids" ||
        simu_para.particle_type == "GNP_CNT_mix")
    {
        //There are GNPs, so read the maximum number of attempts for relocating a GNP
        istr_att >> simu_para.MAX_ATTEMPTS_GNP;
    }
    
    
	return 1;
}
//---------------------------------------------------------------------------
//Read sample geometry
int Input::Read_sample_geometry(Geom_sample &geom_sample, ifstream &infile)
{
	if(geom_sample.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(geom_sample.keywords);
        return 0;
	}
	else geom_sample.mark = true;
    
    //----------------------------------------------------------------------
	//Read the lower-left corner of the sample, its length, width and height
	istringstream istr0(Get_Line(infile));
	istr0 >> geom_sample.sample.poi_min.x >> geom_sample.sample.poi_min.y >> geom_sample.sample.poi_min.z;
	istr0 >> geom_sample.sample.len_x >> geom_sample.sample.wid_y >> geom_sample.sample.hei_z;
	if(geom_sample.sample.len_x<=0||geom_sample.sample.wid_y<=0||geom_sample.sample.hei_z<=0)
	{
		hout << "Error: the dimensions of the sample along each direction should be positive." << endl;
		return 0;
	}
    //Calculate the sample's volume
	geom_sample.volume = geom_sample.sample.len_x*geom_sample.sample.wid_y*geom_sample.sample.hei_z;
    //Calculate the coordinates of the sample's boundaries opposite to those given by the coordinates of origin
    geom_sample.sample.max_x = geom_sample.sample.poi_min.x + geom_sample.sample.len_x;
    geom_sample.sample.max_y = geom_sample.sample.poi_min.y + geom_sample.sample.wid_y;
    geom_sample.sample.max_z = geom_sample.sample.poi_min.z + geom_sample.sample.hei_z;
    
	//----------------------------------------------------------------------
	//Read the maxmimum and minimum side lengths of the observation window and decrement in x, y and z directions
    istringstream istr1(Get_Line(infile));
	istr1 >> geom_sample.win_max_x >> geom_sample.win_max_y >> geom_sample.win_max_z;
	istringstream istr2(Get_Line(infile));
	istr2 >> geom_sample.win_delt_x >> geom_sample.win_delt_y >> geom_sample.win_delt_z;
	istringstream istr3(Get_Line(infile));
	istr3 >> geom_sample.win_min_x >> geom_sample.win_min_y >> geom_sample.win_min_z;

	if(geom_sample.win_max_x<=Zero||geom_sample.win_max_y<=Zero||geom_sample.win_max_z<=Zero||
	   geom_sample.win_max_x>geom_sample.sample.len_x||geom_sample.win_max_y>geom_sample.sample.wid_y||geom_sample.win_max_z>geom_sample.sample.hei_z)
	{
		hout << "Error: the maximum side lenght of the observation window in each direction (win_max) should be positive and smaller than the size of the sample." << endl;
		return 0;
	}
    
	if(geom_sample.win_min_x<=Zero||geom_sample.win_min_y<=Zero||geom_sample.win_min_z<=Zero||
	   geom_sample.win_min_x>geom_sample.win_max_x||geom_sample.win_min_y>geom_sample.win_max_y||geom_sample.win_min_z>geom_sample.win_max_z)
	{
		hout << "Error: the minimum side lenght of the observation window in each direction (win_min) should be positive and smaller than the size of the sample." << endl;
		return 0;
	}
	if(geom_sample.win_delt_x<=Zero||geom_sample.win_delt_y<=Zero||geom_sample.win_delt_z<=Zero)
	{
		hout << "Error: the decrement of the observation window in each direction (win_delt) should be positive." << endl;
		return 0;
	}

	//Details: +Zero for reducing numerical error when dividing by a small number
	int num[3] = {	(int)((geom_sample.win_max_x-geom_sample.win_min_x + Zero)/geom_sample.win_delt_x),
							(int)((geom_sample.win_max_y-geom_sample.win_min_y + Zero)/geom_sample.win_delt_y),
							(int)((geom_sample.win_max_z-geom_sample.win_min_z + Zero)/geom_sample.win_delt_z)	};

	if(num[0]!=num[1]||num[0]!=num[2])
	{
		hout << "Error: the number of obseravtion windows is different on each direction: " << endl;
        hout << "Windows on x = " << num[0] << endl;
        hout << "Windows on y = " << num[1] << endl;
        hout << "Windows on z = " << num[2] << endl;
		return 0;
	}
    else
    {
        //This code might be able to make the observation windows work when
        //the step size is non-integer multiple of (win_max-win_min)
        //However somewhere there is some kind of memory error
        //---AMC Aug 4, 2017
        /*/
        //Check that the last step is the smallest observation window
        int test_x = abs(geom_sample.win_max_x-num[0]*geom_sample.win_delt_x - geom_sample.win_min_x) > Zero;
        int test_y = abs(geom_sample.win_max_y-num[1]*geom_sample.win_delt_y - geom_sample.win_min_y) > Zero;
        int test_z = abs(geom_sample.win_max_z-num[2]*geom_sample.win_delt_z - geom_sample.win_min_z) > Zero;
        
        
        if (test_x || test_y || test_z) {
            //If along any direction the last step results in an observation window smaller than the minimum, then add one more
            num[0] = num[0] + 1;
            hout << "Increased number of observation windows" <<endl;
            return 0;
        }//*/
        geom_sample.cut_num = num[0];
        
    }

	//----------------------------------------------------------------------
	//Read the minimum size for background grids (looking for contact points)
	istringstream istr4(Get_Line(infile));
	istr4 >> geom_sample.gs_minx >> geom_sample.gs_miny >> geom_sample.gs_minz;
	if(geom_sample.gs_minx<=Zero||geom_sample.gs_miny<=Zero||geom_sample.gs_minz<=Zero)
	{
		hout << "Error: The length of the background grids along each direction should be positive." << endl;
        if (geom_sample.gs_minx <= Zero) 
            hout << "Length along x-axis is negative or too close to zero: " << geom_sample.gs_minx << endl;
        if (geom_sample.gs_miny <= Zero)
            hout << "Length along y-axis is negative or too close to zero: " << geom_sample.gs_miny << endl;
        if (geom_sample.gs_minz <= Zero)
            hout << "Length along z-axis is negative or too close to zero: " << geom_sample.gs_minz << endl;
		return 0;
	}
    if (geom_sample.gs_minx > geom_sample.sample.len_x || 
        geom_sample.gs_miny > geom_sample.sample.wid_y ||
        geom_sample.gs_minz > geom_sample.sample.hei_z)
    {
        hout << "Error: The size of the background grids cannot be larger than the sample." << endl;
        if (geom_sample.gs_minx > geom_sample.sample.len_x)
            hout << "Background grids larger along the x direction. Sample length is" << geom_sample.sample.len_x << " background grid length is " << geom_sample.gs_minx << endl;
        if (geom_sample.gs_miny > geom_sample.sample.wid_y)
            hout << "Background grids larger along the y direction. Sample length is" << geom_sample.sample.wid_y << " background grid length is " << geom_sample.gs_miny << endl;
        if (geom_sample.gs_minz > geom_sample.sample.hei_z)
            hout << "Background grids larger along the z direction. Sample length is" << geom_sample.sample.hei_z << " background grid length is " << geom_sample.gs_minz << endl;
    }
	if((int)(geom_sample.win_max_x/geom_sample.gs_minx)>500||
			  (int)(geom_sample.win_max_y/geom_sample.gs_miny)>500||
			  (int)(geom_sample.win_max_z/geom_sample.gs_minz)>500)
	{
		hout << "Error: The number of divisions for background grids is too large (>500). When the number of divisions for background grids in any direction is large (>500) there is a memory error. " << endl;
		return 0;	
	}
    
    //----------------------------------------------------------------------
    //Read the density of the polymer matrix
    istringstream istr5(Get_Line(infile));
    istr5 >> geom_sample.matrix_density;
    if (geom_sample.matrix_density < Zero) {
        hout << "Error: The polymer matrix density must be non-zero. Input was:"<<geom_sample.matrix_density<< endl;
        return 0;
    }

	return 1;
}
//---------------------------------------------------------------------------
//Read cutoff distances
int Input::Read_cutoff_distances(Cutoff_dist &cutoff_dist, ifstream &infile)
{
    if(cutoff_dist.mark)
    {
        //Output a message that the keyword has already been iput
        Warning_message_already_input(cutoff_dist.keywords);
        return 0;
    }
    else cutoff_dist.mark = true;
    
    //Read the the cutoff distances in the following order:
    //van der Waals distance (in microns)
    //cutoff for tunneling (in microns)
    istringstream istr(Get_Line(infile));
    istr >> cutoff_dist.van_der_Waals_dist >> cutoff_dist.tunneling_dist >> cutoff_dist.min_points;
    //hout<<"van_der_Waals_dist="<<cutoff_dist.van_der_Waals_dist<<" tunneling_dist="<<cutoff_dist.tunneling_dist<<endl;
    if (cutoff_dist.van_der_Waals_dist < Zero) {
        hout << "Error: van der Waals distance must be greater than zero. Input was: "<< cutoff_dist.van_der_Waals_dist << endl;
        return 0;
    }
    if (cutoff_dist.tunneling_dist < Zero) {
        hout << "Error: tunneling cutoff distance must be greater than zero. Input was: "<< cutoff_dist.tunneling_dist << endl;
        return 0;
    }
    if (cutoff_dist.min_points < 0) {
        hout << "Error: min_points must be greater than zero. Input was: "<< cutoff_dist.min_points << endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Read the geometric parameters of nanotubes
int Input::Read_nanotube_geo_parameters(Nanotube_Geo &nanotube_geo, ifstream &infile)
{
	if(nanotube_geo.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(nanotube_geo.keywords);
        return 0;
	}
	else nanotube_geo.mark = true;

	//----------------------------------------------------------------------
	//Read the initial growth direction type (random or specific) of CNTs
	istringstream istr0(Get_Line(infile));
	istr0 >> nanotube_geo.dir_distrib_type;
	if(nanotube_geo.dir_distrib_type!="random"&&nanotube_geo.dir_distrib_type!="specific"){ hout << "Error: the type of growth direction must be either random or specific. Input was: " <<nanotube_geo.dir_distrib_type<< endl;	return 0; }
	if(nanotube_geo.dir_distrib_type=="specific")
	{
		istr0  >> nanotube_geo.ini_theta >> nanotube_geo.ini_phi;
		if(nanotube_geo.ini_theta<Zero||nanotube_geo.ini_theta>PI)
		{
			hout << "Error: Specific growth direction was chosen, but the value of theta is not in the valid range of (0, PI). Input was: " <<nanotube_geo.ini_theta<< endl;
			return 0;
		}
        if(nanotube_geo.ini_phi<0||nanotube_geo.ini_phi>=2*PI)
        {
            hout << "Error: Specific growth direction was chosen, but the value of phi is not in the valid range of (0, 2PI). Input was: " <<nanotube_geo.ini_phi<< endl;
            return 0;
        }
	}

	//----------------------------------------------------------------------
	//Read the normal distribution range [-omega, omega] of the growth direction
    istringstream istr_omega(Get_Line(infile));
    istr_omega >> nanotube_geo.omega_a >> nanotube_geo.omega_b;
    if(nanotube_geo.omega_b > 0.5*PI){
        hout << "Error: The maximum value of angle omega is PI/2. Input was: "<<nanotube_geo.omega_a << endl;
        return 0;
    }
    if(nanotube_geo.omega_b < -0.5*PI){
        hout << "Error: The minimum value of angle omega is -PI/2. Input was: "<<nanotube_geo.omega_b << endl;
        return 0;
    }
	//----------------------------------------------------------------------
	//Read the step length (in microns) of nanotube growth
	istringstream istr2(Get_Line(infile));
	istr2 >> nanotube_geo.step_length;
	if(nanotube_geo.step_length<=Zero||
	   nanotube_geo.step_length>=0.25*geom_sample.sample.len_x||
	   nanotube_geo.step_length>=0.25*geom_sample.sample.wid_y||
       nanotube_geo.step_length>=0.25*geom_sample.sample.hei_z) {
        hout << "Error: The step length must be positive and up to 0.25 times the side length of the sample. Input was: "<<nanotube_geo.step_length << endl;	return 0; }
    
	//----------------------------------------------------------------------
    //Read the distribution type (uniform or normal) of the nanotube length (in microns) and its minimum and maximum values
    istringstream istr3(Get_Line(infile));
    istr3 >> nanotube_geo.len_distrib_type;
    if(nanotube_geo.len_distrib_type!="uniform"&&nanotube_geo.len_distrib_type!="normal"){
        hout << "Error: The distribution of the nanotube length can only be either normal or uniform. Input was: "<<nanotube_geo.len_distrib_type << endl;	return 0; }
    istr3 >> nanotube_geo.len_min >> nanotube_geo.len_max;
    if(nanotube_geo.len_min<Zero||nanotube_geo.len_max<Zero||nanotube_geo.len_max<nanotube_geo.len_min){
        hout << "Error: The nanotube length must be non-negative and the minimum lengths must be smaller than the maximum length. Input for minimum was "<< nanotube_geo.len_min<<" and for maximum was "<<nanotube_geo.len_max<<endl; return 0; }
    
    //----------------------------------------------------------------------
    //Read the minimum CNT length ('length') or number of points ('points') to keep a
    //CNT close to the boundary
    istringstream istr_min_length(Get_Line(infile));
    istr_min_length>>nanotube_geo.min_length_type;
    if (nanotube_geo.min_length_type=="length") {
        //Get the input as a double, since it is a length
        double tmp;
        istr_min_length>>tmp;
        //Convert the length to the equivalent number of points
        nanotube_geo.min_points = 1 + (int)(tmp/nanotube_geo.step_length);
    }
    else if (nanotube_geo.min_length_type=="points") {
        istr_min_length>>nanotube_geo.min_points;
    }
    else {
        hout << "Error: The minimum CNT length can only be specified as 'length' or 'points'. Input was: "<<nanotube_geo.min_length_type<<endl; return 0;
    }
    if (nanotube_geo.min_points <= 0) {
        hout << "Error: The minimum CNT length can only be specified as a non-zero length or number of points. Minimum number of points was: "<<nanotube_geo.min_length_type<<endl; return 0;
    }
    
	//----------------------------------------------------------------------
	//Read the distribution type (uniform or normal) of the nanotube radius (in microns) and its minimum and maximum values
    istringstream istr4(Get_Line(infile));
    istr4 >> nanotube_geo.rad_distrib_type;
    if(nanotube_geo.rad_distrib_type!="uniform"&&nanotube_geo.rad_distrib_type!="normal"){
        hout << "Error: The distribution of the nanotube radius can only be either normal or uniform. Input was: "<<nanotube_geo.rad_distrib_type<< endl;	return 0; }
    istr4 >> nanotube_geo.rad_min >> nanotube_geo.rad_max;
    if(nanotube_geo.rad_min<Zero||nanotube_geo.rad_max<Zero||nanotube_geo.rad_max<nanotube_geo.rad_min||
	   nanotube_geo.rad_min>3*nanotube_geo.step_length||nanotube_geo.rad_max>0.05*nanotube_geo.len_min)
	{
        hout << "Error: The radius must be non-negative, the minimum radius must be smaller than the maximum radius, the minimum radius must be smaller than 3*step_length and the maximum radius must be smaller than 0.05*len_min. Input for minimum was "<<nanotube_geo.rad_min<<" and for maximum was "<<nanotube_geo.rad_max<< endl; return 0; }
    
    //---------------------------------------------------------------------------------------
    //Read the density of the CNTs (in gr/micron3)
    istringstream istr5(Get_Line(infile));
    istr5 >> nanotube_geo.density;
    if (nanotube_geo.density <= Zero) {
        hout << "Error: The nanotube density has to be greater than zero. Input was: " << nanotube_geo.density<< endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------------------
    //If only CNTs are to be generated, then copy directly the CNT content
    if (simu_para.particle_type == "CNT_wires" || simu_para.particle_type == "CNT_deposit") {
        
        //Copy the criterion for measuring content
        nanotube_geo.criterion = simu_para.criterion;
        
        if(nanotube_geo.criterion=="vol")
        {
            nanotube_geo.volume_fraction = simu_para.volume_fraction;
            hout << "    The CNT volume fraction is "<< nanotube_geo.volume_fraction << endl;
            
            //Calcualte the total volume of (required) CNTs in the sample
            nanotube_geo.volume = nanotube_geo.volume_fraction*geom_sample.volume;
        }
        else if(nanotube_geo.criterion=="wt")
        {
            nanotube_geo.weight_fraction = simu_para.weight_fraction;
            hout << "    The CNT weight fraction is " << nanotube_geo.weight_fraction << endl;
            
            //Calculate the volume fraction of (required) CNTs in the sample
            double wt_dens = nanotube_geo.weight_fraction/nanotube_geo.density;
            nanotube_geo.volume_fraction = wt_dens/(wt_dens + (1.0 - nanotube_geo.weight_fraction)/geom_sample.matrix_density);
            //hout<<"nanotube_geo.volume_fraction="<<nanotube_geo.volume_fraction<<endl;
            
            //Calcualte the total volume of (required) CNTs in the sample
            nanotube_geo.volume = nanotube_geo.volume_fraction*geom_sample.volume;
            //hout<<"nanotube_geo.volume="<<nanotube_geo.volume<<endl;
            
            //Calculate the total weight of (required) CNTs in the sample
            nanotube_geo.weight = nanotube_geo.volume*nanotube_geo.density;
            //hout<<"nanotube_geo.weight="<<nanotube_geo.weight<<endl;
        }
    }
    //If mixed or hybrid particles are generated, CNT content is calculated after reading the
    //input value for GNP density
    
    //Calculate the extension of the domain along each direction
    //In the case of the CNT deposit, the boundary layer is drastically reduced to
    //reduce memmory usage. This because the boundary layer is also divided into
    //subregions to identify non-penetrating points
    double l_ext = (simu_para.particle_type=="CNT_deposit")? 0.2*nanotube_geo.len_max: nanotube_geo.len_max;
    double l_ext_half = 0.5*l_ext;
    
    //Get the geometry of the extended domain for CNTs
    geom_sample.ex_dom_cnt.poi_min.x = geom_sample.sample.poi_min.x - l_ext_half;
    geom_sample.ex_dom_cnt.poi_min.y = geom_sample.sample.poi_min.y - l_ext_half;
    geom_sample.ex_dom_cnt.len_x = geom_sample.sample.len_x + l_ext;
    geom_sample.ex_dom_cnt.wid_y = geom_sample.sample.wid_y + l_ext;
    //In the case of a CNT deposit, the lowest z-coordinate is tha same as the sample's, 
    //but the length is only incresaed hald the extended length (compared to the sample's)
    if (simu_para.particle_type=="CNT_deposit") {
        geom_sample.ex_dom_cnt.poi_min.z = geom_sample.sample.poi_min.z;
        //geom_sample.ex_dom_cnt.hei_z = geom_sample.sample.hei_z + l_ext_half;
    }
    else {
        geom_sample.ex_dom_cnt.poi_min.z = geom_sample.sample.poi_min.z - l_ext_half;
    }
    //In the case of the z direction, keep the maximum z-coordinate the same
    geom_sample.ex_dom_cnt.hei_z = geom_sample.sample.hei_z + l_ext;

    //Calculate maximum coordinates
    geom_sample.ex_dom_cnt.max_x = geom_sample.ex_dom_cnt.poi_min.x +  geom_sample.ex_dom_cnt.len_x;
    geom_sample.ex_dom_cnt.max_y = geom_sample.ex_dom_cnt.poi_min.y +  geom_sample.ex_dom_cnt.wid_y;
    geom_sample.ex_dom_cnt.max_z = geom_sample.ex_dom_cnt.poi_min.z +  geom_sample.ex_dom_cnt.hei_z;

    //Determine the overlapping of the overlapping sub-regions for CNTs depending on the particle type
    //This section could be read even when only GNPs are used
    if (simu_para.particle_type == "CNT_deposit" || simu_para.particle_type == "CNT_wires")
    {
        //There are only CNTs, then use the calculated overlapping for CNTs
        geom_sample.gs_overlap_cnt = 2 * nanotube_geo.rad_max + cutoff_dist.tunneling_dist;
    }
    else
    {
        //There are no CNTs, so set the overlapping to zero so it does not interfere with 
        //the overlapping for GNPs
        geom_sample.gs_overlap_cnt = 0.0;
    }
    
	return 1;
}
//---------------------------------------------------------------------------
//Readi the geometric parameters of GNPs
int Input::Read_gnp_geo_parameters(GNP_Geo &gnp_geo, ifstream &infile)
{
    if(gnp_geo.mark)
    {
        //Output a message that the keyword has already been iput
        Warning_message_already_input(gnp_geo.keywords);
        return 0;
    }
    else gnp_geo.mark = true;
    
    //----------------------------------------------------------------------
    //Depending on the network type, read:
    //the type of CNT growth on GNP surfaces
    //OR
    //the criterion for minimum GNP volume inside the sample
    istringstream istr0(Get_Line(infile));
    
    //---------------------------------------------------------------------------------------
    //Check the network being generated
    if (simu_para.particle_type == "GNP_cuboids" || simu_para.particle_type == "GNP_CNT_mix") {
        
        //----------------------------------------------------------------------
        //Save the criterion for minimum GNP volume inside the sample
        istr0 >> gnp_geo.vol_in;
        if (gnp_geo.vol_in != "all_in" && gnp_geo.vol_in != "min_in" && gnp_geo.vol_in != "no_min") {
            hout << "Error: The criterion for minimum GNP volume inside the sample can only be all_in, min_in or no_min. Input was: " << gnp_geo.vol_in << endl;
            return 0;
        }
        //If a minimum volume is specified, read it
        if (gnp_geo.vol_in == "min_in" ) {
            istr0 >> gnp_geo.min_vol_in;
            
            //Check the minimum fraction of volume is not a negative number or too small
            if (gnp_geo.min_vol_in < Zero) {
                
                //If the absolute value is greater than zero, then a negative number was input
                if (abs(gnp_geo.min_vol_in) > Zero) {
                    hout << "Error: The minimum fraction of volume of a GNP inside the sample is negative. Input was: " << gnp_geo.min_vol_in << endl;
                    return 0;
                }
                else {
                    
                    //If the absolute value is still less than zero but non-negative,
                    //then a very small number was input
                    //In such case, reset keyword to no_min
                    gnp_geo.vol_in = "no_min";
                    
                    //Send a warning for this reset value
                    hout << "Warning: The input minimum fraction of volume of a GNP inside the sample is too small. Keyword was changed from 'min_in' to 'no_min'."<< endl;
                }
            }
            //Check that the minimum fraction of volume is not greater than 1
            else if (gnp_geo.min_vol_in - 1 > Zero) {
                hout << "Error: The minimum fraction of volume of a GNP inside the sample cannot be larger than 1. Input was:" << gnp_geo.min_vol_in << endl;
                return 0;
            }
        }
    }
    else if (simu_para.particle_type == "Hybrid_particles") {
        
        //----------------------------------------------------------------------
        //Save the type of CNT growth on GNP surfaces
        istr0 >> gnp_geo.growth_type;
        if (gnp_geo.growth_type != "parallel" && gnp_geo.growth_type != "independent") {
            hout << "Error: The growth type of the CNTs on the GNP surface should be either parallel or independent. Input was: " << gnp_geo.growth_type << endl;
            return 0;
        }
    }
    
    //----------------------------------------------------------------------
    //Read the GNP orientation type (random or specific)
    istringstream istr1(Get_Line(infile));
    istr1 >> gnp_geo.orient_distrib_type;
    if(gnp_geo.orient_distrib_type!="random"&&gnp_geo.orient_distrib_type!="specific"){
        hout << "Error: The GNP orientation type must be either random or specific. Input was: " << gnp_geo.orient_distrib_type << endl;
        return 0; }
    if(gnp_geo.orient_distrib_type=="specific")
    {
        istr1  >> gnp_geo.ini_theta >> gnp_geo.ini_phi;
        if(gnp_geo.ini_theta<Zero||gnp_geo.ini_theta>PI)
        {
            hout << "Error: The range specified for theta is not in the valid range of (0, PI)." << endl;
            return 0;
        }
        if(gnp_geo.ini_phi<Zero||gnp_geo.ini_phi>=2*PI)
        {
            hout << "Error: The range specified phi is not in the valid range of (0, 2PI)." << endl;
            return 0;
        }
    }
    
    //----------------------------------------------------------------------
    //Read the distribution type (uniform or normal) of the GNP side length (in microns) and maximum and minimum values
    istringstream istr3(Get_Line(infile));
    istr3 >> gnp_geo.size_distrib_type;
    if(gnp_geo.size_distrib_type!="uniform"&&gnp_geo.size_distrib_type!="normal"){
        hout << "Error: The distribution of the length should be either normal or uniform." << endl;    return 0; }
    istr3 >> gnp_geo.len_min >> gnp_geo.len_max;
    if(gnp_geo.len_min<Zero||gnp_geo.len_max<Zero||gnp_geo.len_max<gnp_geo.len_min){
        hout << "Error: The maximum and minimum values of GNP side length must be non-negative and the minimum value must be smaller than the maximum value. Input for minimum was "<<gnp_geo.len_min<<" and for maximum was "<<gnp_geo.len_max<< endl; return 0; }
    
    //----------------------------------------------------------------------
    //Define the distribution type (uniform or normal) of the GNP thickness (in microns) and maximum and minimum values
    istringstream istr4(Get_Line(infile));
    istr4 >> gnp_geo.thick_distrib_type;
    if(gnp_geo.thick_distrib_type!="uniform"&&gnp_geo.thick_distrib_type!="normal"){
        hout << "Error: The distribution of the GNP thickness should be either normal or uniform. Input was: "<<gnp_geo.thick_distrib_type<< endl;    return 0; }
    istr4 >> gnp_geo.t_min >> gnp_geo.t_max;
    if(gnp_geo.t_min<Zero||gnp_geo.t_max<Zero||gnp_geo.t_max<gnp_geo.t_min) {
        hout << "Error: The maximum and minimum values of GNP thickness must be non-negative and the minimum value must be smaller than the maximum value. Input for minimum was "<<gnp_geo.t_min<<" and for maximum was "<<gnp_geo.t_max<< endl; return 0; }
    
    //----------------------------------------------------------------------
    //Define the GNP density (in gr/micron3)
    istringstream istr6(Get_Line(infile));
    istr6 >> gnp_geo.density;
    if (gnp_geo.density <= Zero) {
        hout << "Error: The GNP density has to be a value greater than zero. Input was:" << nanotube_geo.density<< endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------------------
    //If only GNPs are to be generated, then copy directly the total carbon content
    if (simu_para.particle_type == "GNP_cuboids") {
        
        //Copy the criterion for measuring content
        gnp_geo.criterion = simu_para.criterion;
        
        if (gnp_geo.criterion=="vol")
        {
            
            gnp_geo.volume_fraction = simu_para.volume_fraction;
            hout << "    The GNP volume fraction is "<< gnp_geo.volume_fraction << endl;
            
            //Calculate the actual volume of GNPs
            gnp_geo.volume = gnp_geo.volume_fraction*geom_sample.volume;
        }
        else if (gnp_geo.criterion=="wt")
        {

            gnp_geo.weight_fraction = simu_para.weight_fraction;
            hout << "    The GNP weight fraction is "<< gnp_geo.weight_fraction << endl;
            
            //Calculate the volume fraction of (required) GNPs in the sample
            double wt_dens = gnp_geo.weight_fraction/gnp_geo.density;
            gnp_geo.volume_fraction = wt_dens/(wt_dens + (1.0 - gnp_geo.weight_fraction)/geom_sample.matrix_density);
            
            //Calcualte the total volume of (required) GNPs in the sample
            gnp_geo.volume = gnp_geo.volume_fraction*geom_sample.volume;
            
            //Calculate the total weight of (required) GNPs in the sample
            gnp_geo.weight = gnp_geo.volume*gnp_geo.density;
        }
    }
    //If mixed or hybrid particles are generated, calculate the fraction using the mass/volume ratio
    else if (simu_para.particle_type == "GNP_CNT_mix" || simu_para.particle_type == "Hybrid_particles") {
        
        //Check what is the criterion for mixed/hybrid particles
        if (simu_para.mixed == "mass_ratio") {
            
            //Check if the total carbon content is given as volume or weight fraction
            if (simu_para.criterion=="wt") {
                
                //Set the criterion as weight fraction
                nanotube_geo.criterion = "wt";
                
                //Calculate the weight fraction of (required) CNTs in the sample
                nanotube_geo.weight_fraction = simu_para.mass_ratio*simu_para.weight_fraction/(simu_para.mass_ratio + 1.0);
                hout<<"nanotube_geo.weight_fraction="<<nanotube_geo.weight_fraction<<endl;
                
                //Calculate the weight fraction of (required) GNPs in the sample
                gnp_geo.weight_fraction = simu_para.weight_fraction/(simu_para.mass_ratio + 1.0);
                hout<<"gnp_geo.weight_fraction="<<gnp_geo.weight_fraction<<endl;
                
                //Calculate the denominator used in calculating volume fractions
                double den = (1.0 - simu_para.weight_fraction)/geom_sample.matrix_density;
                den += nanotube_geo.weight_fraction/nanotube_geo.density;
                den += gnp_geo.weight_fraction/gnp_geo.density;
                
                //Calculate the volume fraction of (required) CNTs in the sample
                nanotube_geo.volume_fraction = nanotube_geo.weight_fraction/(nanotube_geo.density*den);
                
                //Calculate the volume fraction of (required) GNPs in the sample
                gnp_geo.volume_fraction = gnp_geo.weight_fraction/(gnp_geo.density*den);
                
            }
            else if (simu_para.criterion=="vol") {
                
                //Calculate the volume fraction of (required) CNTs in the sample
                nanotube_geo.volume_fraction = simu_para.mass_ratio*gnp_geo.density*simu_para.volume_fraction/(nanotube_geo.density + simu_para.mass_ratio*gnp_geo.density);
                
                //Calculate the volume fraction of (required) GNPs in the sample
                gnp_geo.volume_fraction = simu_para.volume_fraction - nanotube_geo.volume_fraction;
                
            }
            else {
                hout<<"Error when calculating CNT and GNP content from mass ratio: Invalid criterion for measuring content. Criterion can only be \"vol\" or \"wt\". Input was:"<<simu_para.criterion<<endl;
                return 0;
            }
            
            //Calcualte the total volume of (required) CNTs in the sample
            nanotube_geo.volume = nanotube_geo.volume_fraction*geom_sample.volume;
            
            //Calculate the total weight of (required) CNTs in the sample
            nanotube_geo.weight = nanotube_geo.volume*nanotube_geo.density;
            
            //Calcualte the total volume of (required) GNPs in the sample
            gnp_geo.volume = gnp_geo.volume_fraction*geom_sample.volume;
            
            //Calculate the total weight of (required) GNPs in the sample
            gnp_geo.weight = gnp_geo.volume*gnp_geo.density;
        }
        else if (simu_para.mixed == "volume_ratio") {
            
            //Set the criterion as volume fraction
            nanotube_geo.criterion = "vol";
            
            //Calculate the volume of CNTs from volume ratio
            nanotube_geo.volume = simu_para.volume_ratio*simu_para.volume_fraction*geom_sample.volume/(simu_para.volume_ratio + 1.0);
            
            //Calculate the volume of GNPs from volume ratio
            gnp_geo.volume = simu_para.volume_fraction*geom_sample.volume/(simu_para.volume_ratio + 1.0);
        }
        //If density is selected, the amount of CNTs is determined at generation time
        else if (simu_para.mixed == "density") {
            
            //Set the criterion as volume fraction
            gnp_geo.criterion = "vol";
            
            //Calculate the volume ratio using the smallest CNT geometry and largest GNP geometry
            double vol_ratio = 2*simu_para.cnt_gnp_density*PI*nanotube_geo.rad_min*nanotube_geo.rad_min*nanotube_geo.len_min/gnp_geo.t_max;
            
            //Calculate the GNP volume using the calculated volume ratio
            gnp_geo.volume = simu_para.volume_fraction*geom_sample.volume/(vol_ratio + 1.0);
        }
    }
    
    
    
    //---------------------------------------------------------------------------------------
    //Get the geometry of the extended domain for GNPs
    double len_max_halved = gnp_geo.len_max/2.0;
    geom_sample.ex_dom_gnp.poi_min = geom_sample.sample.poi_min - Point_3D(len_max_halved,len_max_halved,len_max_halved);
    geom_sample.ex_dom_gnp.len_x = geom_sample.sample.len_x + len_max_halved;
    geom_sample.ex_dom_gnp.wid_y = geom_sample.sample.wid_y + len_max_halved;
    geom_sample.ex_dom_gnp.hei_z = geom_sample.sample.hei_z + len_max_halved;
    geom_sample.ex_dom_gnp.max_x = geom_sample.ex_dom_gnp.poi_min.x + geom_sample.ex_dom_gnp.len_x;
    geom_sample.ex_dom_gnp.max_y = geom_sample.ex_dom_gnp.poi_min.y + geom_sample.ex_dom_gnp.wid_y;
    geom_sample.ex_dom_gnp.max_z = geom_sample.ex_dom_gnp.poi_min.z + geom_sample.ex_dom_gnp.hei_z;

    //Determine the overlapping of the overlapping sub-regions for GNPs depending on the particle type
    //This section could be read even when only CNTs are used
    if (simu_para.particle_type == "GNP_cuboids" || simu_para.particle_type == "GNP_CNT_mix")
    {
        //There are only GNPs, then use the calculated overlapping for GNPs
        geom_sample.gs_overlap_gnp = geom_sample.gs_minx / (sqrt(8.0));
    }
    else 
    {
        //There are no GNPs, so set the overlapping to zero so it does not interfere with 
        //the overlapping for CNTs
        geom_sample.gs_overlap_gnp = 0.0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Read electrical properties of materials
int Input::Read_electrical_parameters(Electric_para &electric_para, ifstream &infile)
{
    if(electric_para.mark)
    {
        //Output a message that the keyword has already been iput
        Warning_message_already_input(electric_para.keywords);
        return 0;
    }
    else electric_para.mark = true;
    
    //Read applied voltage
    istringstream istr0(Get_Line(infile));
    istr0 >> electric_para.applied_voltage;
    if (electric_para.applied_voltage <= Zero) {
        hout<<"Error: Voltage cannot be negative. Input was: "<<electric_para.applied_voltage<<endl;
        return 0;
    }
    
    //Carbon nanotube resistivity
    istringstream istr1(Get_Line(infile));
    istr1 >> electric_para.resistivity_CNT;
    if (electric_para.resistivity_CNT <= Zero) {
        hout<<"Error: CNT resistivity cannot be negative. Input was: "<<electric_para.resistivity_CNT<<endl;
        return 0;
    }
    
    //GNP resistivities
    istringstream istr2(Get_Line(infile));
    istr2 >> electric_para.resistivity_GNP_t >> electric_para.resistivity_GNP_surf;
    //hout << electric_para.resistivity_GNP_t <<" ,"<<electric_para.resistivity_GNP_surf<<endl;
    if (electric_para.resistivity_GNP_t <= Zero || electric_para.resistivity_GNP_surf <= Zero) {
        hout<<"Error: GNP resistivities cannot be negative. Input was: "<<electric_para.resistivity_GNP_t<<", "<<electric_para.resistivity_GNP_surf<<endl;
        return 0;
    }
    
    //Polymer matrix resistivity
    istringstream istr3(Get_Line(infile));
    istr3 >> electric_para.resistivity_matrix;
    if (electric_para.resistivity_matrix <= Zero) {
        hout<<"Error: Polymer resistivity cannot be negative. Input was: "<<electric_para.resistivity_matrix<<endl;
        return 0;
    }
    //hout << "electric_para.resistivity_matrix = " << electric_para.resistivity_matrix << endl;
    //Electrical constants: electron charge (C), permitivity of vaccum (F/m), CNT work function (V), dielectric constant of polymer
    
    //Electrical constants to calculate junction resistance
    istringstream istr4(Get_Line(infile));
    //Check the type of junction resistance
    istr4 >> electric_para.junction_type;
    if (electric_para.junction_type == "constant") {
        
        //Read the constant value for junction resistance
        istr4 >> electric_para.junction_resistance;
        if (electric_para.junction_resistance < Zero) {
            hout<<"Error: Junction resistance cannot be negative. Input was: "<<electric_para.junction_resistance<<endl;
            return 0;
        }
    }
    else if (electric_para.junction_type == "exponential") {
        
        //Electrical constants for exponential function:
        //Planck’s constant (m2kg/s), electron charge (C), electron mass (Kg), height barrier (eV)
        istr4 >> electric_para.h_plank >> electric_para.e_charge >> electric_para.e_mass >> electric_para.lambda_barrier;
        if (electric_para.h_plank<=Zero || electric_para.e_charge<=Zero || electric_para.e_mass<=Zero || electric_para.lambda_barrier<=Zero) {
            hout<<"Error: Electric parameters cannot be negative. Input was: "<<electric_para.h_plank<<", "<<electric_para.e_charge<<", "<<electric_para.e_mass<<", "<<electric_para.lambda_barrier<<endl;
            return 0;
        }

        //Calculate a squared root term that is used twice in the equation
        double sqrt_tmp = sqrt(2.0 * electric_para.e_mass * electric_para.lambda_barrier * electric_para.e_charge);

        //Calculate C1
        electric_para.C1 = 10.0 * electric_para.h_plank * electric_para.h_plank / (electric_para.e_mass * electric_para.e_mass * sqrt_tmp);

        //Calculate C2
        electric_para.C2 = 4000.0 * PI * sqrt_tmp / electric_para.h_plank;
    }
    else {
        hout << "Error: the junction type is neither 'constant' nor 'exponential'." << endl;
        hout << "Input was: " << electric_para.junction_type << endl;
        return 0;
    }
    //hout << electric_para.h_plank <<", "<< electric_para.e_charge <<", "<< electric_para.e_mass <<", "<< electric_para.lambda_barrier << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//Read flags for visualization files
int Input::Read_visualization_flags(Visualization_flags &vis_flags, ifstream &infile)
{
	if(vis_flags.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(vis_flags.keywords);
        return 0;
	}
	else vis_flags.mark = true;
    
    //Flag to export generated nanoparticles:
    // 0: do not export visualization files
    // 1: export visualization files (if mixed or hybrid particles generated,
    //    then two files are generated, one for CNTs and one for GNPs)
    istringstream istr0(Get_Line(infile));
	istr0 >> vis_flags.generated_nanoparticles;
    if (vis_flags.generated_nanoparticles<0||vis_flags.generated_nanoparticles>3) {
        hout<<"Error: Flag to export generated CNTs can only be 0, 1, 2, or 3. Input was: "<<vis_flags.generated_nanoparticles<<endl;
        return 0;
    }
    //The flag for the sample domain has the same value as the flag for generated_nanoparticles
    vis_flags.sample_domain = vis_flags.generated_nanoparticles;
    
    //Flag to export clusters as obtained from the Hoshen-Kopelman algorithm:
    // 0: do not export visualization files
    // 1: export visualization files
    istringstream istr2(Get_Line(infile));
    istr2 >> vis_flags.clusters;
    if (vis_flags.clusters<0||vis_flags.clusters>1) {
        hout<<"Error: Flag to export clusters (from HK76 algorithm) can only be an integer, 0 or 1. Input was: "<<vis_flags.clusters<<endl;
        return 0;
    }

    //Flag to export percolated clusters:
    // 0: do not export visualization files
    // 1: export visualization files
    istringstream istr3(Get_Line(infile));
    istr3 >> vis_flags.percolated_clusters;
    if (vis_flags.percolated_clusters<0||vis_flags.percolated_clusters>1) {
        hout<<"Error: Flag to export percolated clusters can only be an integer, 0 or 1. Input was: "<<vis_flags.percolated_clusters<<endl;
        return 0;
    }
    
    //Flag to export the backbone:
    // 0: do not export visualization files
    // 1: export visualization files
    //If set to 1, up to six files are generated:
    //   Three files if there are CNTs: backbone, dead branches and isolated
    //   Three files if there are GNPs: backbone, dead (attached to the percolated cluster) and isolated
    istringstream istr4(Get_Line(infile));
    istr4 >> vis_flags.backbone;
    if (vis_flags.backbone<0||vis_flags.backbone>1) {
        hout<<"Error: Flag to export the backbone clusters can only be an integer, 0 or 1. Input was: "<<vis_flags.backbone<<endl;
        return 0;
    }
    
    //Flag to export triangulations if set to 1
    // 0: do not export triangulations files
    // 1: export triangulations
    istringstream istr5(Get_Line(infile));
    istr5 >> vis_flags.triangulations;
    if (vis_flags.triangulations<0||vis_flags.triangulations>1) {
        hout<<"Error: Flag to export triangulation can only be 0 or 1. Input was: "<<vis_flags.triangulations<<endl;
        return 0;
    }
    
    //Update the flag for exporting the observation window if needed
    if (vis_flags.clusters || vis_flags.percolated_clusters || vis_flags.backbone) {
        vis_flags.window_domain = 1;
    }
    
	return 1;
}
//---------------------------------------------------------------------------
//Read flags for output data files
int Input::Read_output_data_flags(Output_data_flags &out_flags, ifstream &infile)
{
    if(out_flags.mark)
    {
        //Output a message that the keyword has already been iput
        Warning_message_already_input(out_flags.keywords);
        return 0;
    }
    else out_flags.mark = true;
    
    //Flag to save CNT and GNP fractions and volume separately when mixed or hybrid particles are generated
    // 0: do not export CNT and GNP volumes and fractions separately
    // 1: export CNT and GNP volumes and fractions separately
    //When this flag is set to 1, in addition to total fractions and volumes, four more files are written:
    //Volumes and fractions of CNTs only (fractions respect to the CNT volume)
    //Volumes and fractions of GNPs only (fractions respect to the GNP volume)
    //If only CNTs or only GNPs are generated, then this flag is ignored
    istringstream istr_cnt_gnp_flag(Get_Line(infile));
    istr_cnt_gnp_flag >> out_flags.cnt_gnp_flag;
    //Check it is a valid flag
    if (out_flags.cnt_gnp_flag<0||out_flags.cnt_gnp_flag>1) {
        hout<<"Error: Flag to export CNT and GNP fractions and volume separately can only be 0 or 1. Input was: "<<out_flags.cnt_gnp_flag<<endl;
        return 0;
    }
    //Reset the flag to zero if only CNTs or only GNPs are generated
    if (simu_para.particle_type == "CNT_wires"||simu_para.particle_type == "CNT_deposit"|| simu_para.particle_type == "GNP_cuboids") {
        out_flags.cnt_gnp_flag = 0;
    }
    
    //Flag to output a text file with four points so that GNPs may be generated in Abaqus
    //The coice of four points depends entirely in the way the GNPs are generated in Abaqus:
    //Three points define a plane which is used to draw the squared base of the GNP
    //The fourth point is used to define the thickness of the GNP
    istringstream istr_gnp_data(Get_Line(infile));
    istr_gnp_data >> out_flags.gnp_data;
    //Check it is a valid flag
    if (out_flags.gnp_data<0||out_flags.gnp_data>4) {
        hout<<"Error: Flag to export GNP data should be an integer between 0 and 4. Input was: "<<out_flags.gnp_data<<endl;
        return 0;
    }
    //Check if precision is needed
    //Only csv files require a precision for exporting into a file
    if (out_flags.gnp_data == 1 || out_flags.gnp_data == 3) {
        
        //Flag was set to 1, so read the precision
        istr_gnp_data >> out_flags.prec_gnp;
        
        //Check the precision is a positive integer
        if (out_flags.prec_gnp <= 0) {
            hout<<"Error: Precision used to export GNP vertices must be at least 1. Input was: "<<out_flags.prec_gnp<<endl;
            return 0;
        }
    }
    
    //Flag to output a text file with the coordinates of all CNT points
    //Two files are exported, one with the coordinates and one with the number of CNTs and
    //the number of points for each CNT
    istringstream istr_cnt_data(Get_Line(infile));
    istr_cnt_data >> out_flags.cnt_data;
    //Check it is a valid flag
    if (out_flags.cnt_data<0||out_flags.cnt_data>3) {
        hout<<"Error: Flag to export CNT data should be an integer between 0 and 3. Input was: "<<out_flags.cnt_data<<endl;
        return 0;
    }
    //Check if precision is needed
    //Only csv files require a precision for exporting into a file
    if (out_flags.cnt_data == 1 || out_flags.cnt_data == 3) {
        
        //Flag was set to 1, so read the precision
        istr_cnt_data >> out_flags.prec_cnt;
        
        //Check the precision is a positive integer
        if (out_flags.prec_cnt <= 0) {
            hout<<"Error: Precision used to export CNT points must be at least 1. Input was: "<<out_flags.prec_cnt<<endl;
            return 0;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Read the input data in a whole line (to skip over the comment line starting with a '%')
string Input::Get_Line(ifstream &infile)const
{
	string s;
	//Read the input data in a whole line
	getline(infile,s);
	//to skip over the comment line starting with a '%'
	while(!infile.eof() && s.substr(0,1)=="%")
		getline(infile,s);
	return s;
}
//===========================================================================
