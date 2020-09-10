//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
//OBJECTIVE:    Read input data and save it into data structures
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Input_Reader.h"

//---------------------------------------------------------------------------
//Read data
int Input::Read_Infile(ifstream &infile)
{

	hout << "Reading input file..." << endl;
    
    //Go through each section of the input file in order
    //int i = 0;
    
    //String to save the line read from the input file
    string str;
    
    //Application type
    str = Get_Line(infile);
    if(str=="Application") {
        if(Read_application(app_name, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Application");
    }
    
    //Parameters to run the simulation
    str = Get_Line(infile);
    if(str=="Simulation_Parameters")    {
        if(Read_simulation_parameters(simu_para, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Simulation_Parameters");
    }
    
    
    str = Get_Line(infile);
    if(str=="Sample_Geometry")    {
        if(Read_sample_geometry(geom_sample, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Sample_Geometry");
    }
    
    
    str = Get_Line(infile);
    if(str=="Nanotube_Geometry")    {
        if(Read_nanotube_geo_parameters(nanotube_geo, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Nanotube_Geometry");
    }
    
    
    str = Get_Line(infile);
    if(str=="CNT_Agglomerate_Geometry")    {
        if(Read_agg_geo_parameters(agg_geo, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("CNT_Agglomerate_Geometry");
    }
    
    str = Get_Line(infile);
    if(str=="GNP_Geometry")    {
        if(Read_gnp_geo_parameters(gnp_geo, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("GNP_Geometry");
    }
    
    str = Get_Line(infile);
    if(str=="Cutoff_Distances")    {
        if(Read_cutoff_distances(cutoff_dist, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Cutoff_Distances");
    }
    
    str = Get_Line(infile);
    if(str=="Electrical_Parameters")    {
        if(Read_electrical_parameters(electric_para, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Electrical_Parameters");
    }
    
    str = Get_Line(infile);
    if(str=="Tecplot_flags")    {
        if(Read_tecplot_flags(tec360_flags, infile)==0) return 0;
        //hout<<"i="<<i<<endl;i++;
    }
    //If not found then the default values are used
    else {
        Warning_message("Tecplot_flags");
    }

	hout << "Finished reading input file!" << endl;
    
    
    //Once the input file has been read or default values have been set,
    //calculate the CNT volume and the geometry of the boudary layer
    
    //Calculate the actual volume of nanotubes
    nanotube_geo.real_volume = nanotube_geo.volume_fraction * geom_sample.volume;
    
    //Define the geometry of the boundary layer
    geom_sample.ex_origin.x = geom_sample.origin.x - nanotube_geo.len_max;
    geom_sample.ex_origin.y = geom_sample.origin.y - nanotube_geo.len_max;
    geom_sample.ex_origin.z = geom_sample.origin.z - nanotube_geo.len_max;
    geom_sample.ex_len = geom_sample.len_x + 2*nanotube_geo.len_max;
    geom_sample.ey_wid = geom_sample.wid_y+ 2*nanotube_geo.len_max;
    geom_sample.ez_hei = geom_sample.hei_z + 2*nanotube_geo.len_max;
    
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
    simu_para.avoid_resistance = 0;
    

	//Initialize the geometric parameters of the RVE
	geom_sample.keywords = "Sample_Geometry";
	geom_sample.mark = false;
    geom_sample.origin.x = 0.0;
	geom_sample.origin.y = 0.0;
	geom_sample.origin.z = 0.0;
	geom_sample.origin.flag = 0;
	geom_sample.len_x = 1.0;
	geom_sample.wid_y = 1.0;
	geom_sample.hei_z = 1.0;
	geom_sample.ex_origin.x = 0.0;
	geom_sample.ex_origin.y = 0.0;
	geom_sample.ex_origin.z = 0.0;
	geom_sample.ex_origin.flag = 0;
	geom_sample.ex_len = 1.0;
	geom_sample.ey_wid = 1.0;
	geom_sample.ez_hei = 1.0;
	geom_sample.volume = geom_sample.len_x*geom_sample.wid_y*geom_sample.hei_z;
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
	nanotube_geo.keywords = "Nanotube_Geometry";
	nanotube_geo.mark = false;
	nanotube_geo.dir_distrib_type = "random";
	nanotube_geo.ini_pha = 0.0;
	nanotube_geo.ini_sita = 0.0;
	nanotube_geo.angle_max = 1.5707963267948966;
	nanotube_geo.step_length = 0.01;
	nanotube_geo.len_distrib_type = "uniform";
	nanotube_geo.len_min = 1.0;
	nanotube_geo.len_max = 1.0;
	nanotube_geo.rad_distrib_type = "uniform";
	nanotube_geo.rad_min = 0.005;
	nanotube_geo.rad_max = 0.005;
	nanotube_geo.criterion = "vol";
	nanotube_geo.volume_fraction = 0.0;
	nanotube_geo.real_volume = 0.0;
	nanotube_geo.weight_fraction = 0.0;
	nanotube_geo.real_weight = 0.0;
	nanotube_geo.linear_density = 5.8E-5;

	//Initialize the geometric paramters of nanotube clusters
	agg_geo.keywords = "CNT_Agglomerate_Geometry";
	agg_geo.mark = false;
	agg_geo.print_key = 1;
	agg_geo.vol_fra_criterion = 0;
	agg_geo.amin = 0.0;
	agg_geo.amax = 0.0;
	agg_geo.bmin = 0.0;
	agg_geo.cmin = 0.0;
	agg_geo.growth_probability = 0.0;
	agg_geo.volf_clust = 0.0;
	agg_geo.cnt_real_volume = 0.0;

	//Initialize cutoff distances
	cutoff_dist.keywords = "Cutoff_Distances";
	cutoff_dist.mark = false;
	cutoff_dist.tunneling_dist = 0.0018;
	cutoff_dist.van_der_Waals_dist = 0.00034;

	//Initialize electrical parameters
	electric_para.keywords = "Electrical_Parameters";
	electric_para.mark = false;
	electric_para.applied_voltage = 1.0;
	electric_para.resistivity_CNT = 0.001;
    
    //Initialize tecplot flags (do not export anything)
    tec360_flags.keywords = "Tecplot_flags";
    tec360_flags.mark = false;
    tec360_flags.generated_cnts = 0;
    tec360_flags.generated_gnps = 0;
    tec360_flags.clusters = 0;
    tec360_flags.percolated_clusters = 0;
    tec360_flags.backbone = 0;
    tec360_flags.triangulations = 0;

	hout << "    Data initialization files" <<endl<<endl;

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
int Input::Read_application(struct App_name &app_name, ifstream &infile)
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
int Input::Read_simulation_parameters(struct Simu_para &simu_para, ifstream &infile)
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
    if (simu_para.particle_type != "CNT_wires" && simu_para.particle_type != "GNP_cuboids" && simu_para.particle_type != "Hybrid_particles" && simu_para.particle_type != "GNP_CNT_mix") {
        hout << "Error: the type of particles shoud be one of the following: CNT_wires, GNP_cuboids, Hybrid_particles or GNP_CNT_mix. Input was: "<<simu_para.particle_type<< endl;
        return 0;
    }
    
    //Read keyword for creating a new network or reading a network from a file
	istringstream istr2(Get_Line(infile));
	istr2 >> simu_para.create_read_network;
	if(simu_para.create_read_network!="Create_Network"&&simu_para.create_read_network!="Read_Network"&&simu_para.create_read_network!="Network_From_Seeds")
	{
        hout << "Error: Invalid keyword. Valid options are 'Create_Network', 'Read_Network', or 'Network_From_Seeds'. Input was: " << simu_para.create_read_network << endl; return 0; }
    
    //----------------------------------------------------------------------
    //If the network is generated from some seeds, then read them
    if (simu_para.create_read_network=="Network_From_Seeds") {
        
        //Different types of fillers need different number of seed
        //Thus, check the filler type and read the necessary number of seeds
        if (simu_para.particle_type=="CNT_wires") {
            
            //7 seeds are required for CNT_wires
            simu_para.network_seeds.assign(7, 0);
            
            //Read the line with the seeds
            istringstream istr_network_seeds(Get_Line(infile));
            //Add seeds to the vector
            for (int i = 0; i < 7; i++) {
                istr_network_seeds >> simu_para.network_seeds[i];
                //hout << simu_para.network_seeds[i] << endl;
            }
        }
        else {
            hout<<"Error: Seeds can only be read with the keyword \"CNT_wires\" for nanoparticle type. Input for particle type was: "<<simu_para.particle_type<<endl; return 0;
        }
    }
    
    //Flag to avoid calculating the resistance of the network
    //This can be useful when only a geometric study or visualizations are needed
    // 1: Avoid calculating the resistance of the network
    // 0: Calculate the resistance of the network
    istringstream istr3(Get_Line(infile));
    istr3 >> simu_para.avoid_resistance;
    if (simu_para.avoid_resistance < 0 || simu_para.avoid_resistance > 1) {
        hout << "Error: Invalid value for network resistance calculation. Valid options are 0 (avoid calculating network resistance) or 1 (calculate network resistance). Input was: " << simu_para.avoid_resistance << endl; return 0;
    }
    
	return 1;
}
//---------------------------------------------------------------------------
//Read sample geometry
int Input::Read_sample_geometry(struct Geom_sample &geom_sample, ifstream &infile)
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
	istr0 >> geom_sample.origin.x >> geom_sample.origin.y >> geom_sample.origin.z;
	istr0 >> geom_sample.len_x >> geom_sample.wid_y >> geom_sample.hei_z;
	if(geom_sample.len_x<=0||geom_sample.wid_y<=0||geom_sample.hei_z<=0)
	{
		hout << "Error: the dimensions of the sample along each direction should be positive." << endl;
		return 0;
	}
	geom_sample.volume = geom_sample.len_x*geom_sample.wid_y*geom_sample.hei_z;

	//----------------------------------------------------------------------
	//Read the maxmimum and minimum side lengths of the observation window and decrement in x, y and z directions
    istringstream istr1(Get_Line(infile));
	istr1 >> geom_sample.win_max_x >> geom_sample.win_max_y >> geom_sample.win_max_z;
	istringstream istr2(Get_Line(infile));
	istr2 >> geom_sample.win_delt_x >> geom_sample.win_delt_y >> geom_sample.win_delt_z;
	istringstream istr3(Get_Line(infile));
	istr3 >> geom_sample.win_min_x >> geom_sample.win_min_y >> geom_sample.win_min_z;

	if(geom_sample.win_max_x<=0.0||geom_sample.win_max_y<=0.0||geom_sample.win_max_z<=0.0||
	   geom_sample.win_max_x>geom_sample.len_x||geom_sample.win_max_y>geom_sample.wid_y||geom_sample.win_max_y>geom_sample.hei_z)
	{
		hout << "Error: the maximum side lenght of the observation window in each direction (win_max) should be positive and smaller than the size of the sample." << endl;
		return 0;
	}
    
	if(geom_sample.win_min_x<=0.0||geom_sample.win_min_y<=0.0||geom_sample.win_min_z<=0.0||
	   geom_sample.win_min_x>geom_sample.win_max_x||geom_sample.win_min_y>geom_sample.win_max_y||geom_sample.win_min_z>geom_sample.win_max_z)
	{
		hout << "Error: the minimum side lenght of the observation window in each direction (win_min) should be positive and smaller than the size of the sample." << endl;
		return 0;
	}
	if(geom_sample.win_delt_x<=0.0||geom_sample.win_delt_y<=0.0||geom_sample.win_delt_z<=0.0)
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
		hout << "Error: the number of steps is different on each direction: " << endl;
        hout << "\tSteps on x = " << num[0] << endl;
        hout << "\tSteps on y = " << num[1] << endl;
        hout << "\tSteps on z = " << num[2] << endl;
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
	if(geom_sample.gs_minx<=0||geom_sample.gs_miny<=0||geom_sample.gs_minz<=0)
	{
		hout << "Error: the number of divisions for background grids in each direction of the sample should be positive." << endl;
		return 0;
	}
	else if((int)(geom_sample.win_max_x/geom_sample.gs_minx)>500||
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

	return 1;
}
//---------------------------------------------------------------------------
//Read the geometric parameters of nanotubes
int Input::Read_nanotube_geo_parameters(struct Nanotube_Geo &nanotube_geo, ifstream &infile)
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
		istr0  >> nanotube_geo.ini_sita >> nanotube_geo.ini_pha;
		if(nanotube_geo.ini_sita<0||nanotube_geo.ini_sita>PI)
		{
			hout << "Error: Specific growth direction was chosen, but the value of theta is not in the valid range of (0, PI). Input was: " <<nanotube_geo.ini_sita<< endl;
			return 0;
		}
        if(nanotube_geo.ini_pha<0||nanotube_geo.ini_pha>=2*PI)
        {
            hout << "Error: Specific growth direction was chosen, but the value of phi is not in the valid range of (0, 2PI). Input was: " <<nanotube_geo.ini_pha<< endl;
            return 0;
        }
	}

	//----------------------------------------------------------------------
	//Read the normal distribution range [-omega, omega] of the growth direction
    istringstream istr1(Get_Line(infile));
	istr1 >> nanotube_geo.angle_max;
	if(nanotube_geo.angle_max>0.5*PI){ hout << "Error: The angle omega is not in the valid range of (-PI/2, PI/2). Input was: "<<nanotube_geo.angle_max << endl;	 return 0; }
    
	//----------------------------------------------------------------------
	//Read the step length (in microns) of nanotube growth
	istringstream istr2(Get_Line(infile));
	istr2 >> nanotube_geo.step_length;
	if(nanotube_geo.step_length<=0||
	   nanotube_geo.step_length>=0.25*geom_sample.len_x||
	   nanotube_geo.step_length>=0.25*geom_sample.wid_y||
       nanotube_geo.step_length>=0.25*geom_sample.hei_z)
	{ hout << "Error: The step length must be positive and upt to 0.25 times the side length of the sample. Input was: "<<nanotube_geo.step_length << endl;	return 0; }

	//----------------------------------------------------------------------
    //Read the distribution type (uniform or normal) of the nanotube length (in microns) and its minimum and maximum values
    istringstream istr3(Get_Line(infile));
    istr3 >> nanotube_geo.len_distrib_type;
    if(nanotube_geo.len_distrib_type!="uniform"&&nanotube_geo.len_distrib_type!="normal"){
        hout << "Error: The distribution of the nanotube length can only be either normal or uniform. Input was: "<<nanotube_geo.len_distrib_type << endl;	return 0; }
    istr3 >> nanotube_geo.len_min >> nanotube_geo.len_max;
    if(nanotube_geo.len_min<0||nanotube_geo.len_max<0||nanotube_geo.len_max<nanotube_geo.len_min){
        hout << "Error: The nanotube length must be non-negative and the minimum lengths must be smaller than the maximum length. Input for minimum was "<< nanotube_geo.len_min<<" and for maximum was "<<nanotube_geo.len_max<<endl; return 0; }

	//----------------------------------------------------------------------
	//Read the distribution type (uniform or normal) of the nanotube radius (in microns) and its minimum and maximum values
    istringstream istr4(Get_Line(infile));
    istr4 >> nanotube_geo.rad_distrib_type;
    if(nanotube_geo.rad_distrib_type!="uniform"&&nanotube_geo.rad_distrib_type!="normal"){
        hout << "Error: The distribution of the nanotube radius can only be either normal or uniform. Input was: "<<nanotube_geo.rad_distrib_type<< endl;	return 0; }
    istr4 >> nanotube_geo.rad_min >> nanotube_geo.rad_max;
    if(nanotube_geo.rad_min<0||nanotube_geo.rad_max<0||nanotube_geo.rad_max<nanotube_geo.rad_min||
	   nanotube_geo.rad_min>3*nanotube_geo.step_length||nanotube_geo.rad_max>0.05*nanotube_geo.len_min)
	{
        hout << "Error: The radius must be non-negative, the minimum radius must be smaller than the maximum radius, the minimum radius must be smaller than 3*step_length and the maximum radius must be smaller than 0.05*len_min. Input for minimum was "<<nanotube_geo.rad_min<<" and for maximum was "<<nanotube_geo.rad_max<< endl; return 0; }
    
    //---------------------------------------------------------------------------------------
    //Read the density of the CNTs (in gm/cm3)
    istringstream istr5(Get_Line(infile));
    istr5 >> nanotube_geo.density;
    if (nanotube_geo.density <= 0) {
        hout << "Error: The nanotube density has to be greater than zero. Input was: " << nanotube_geo.density<< endl;
        return 0;
    }
    
    //----------------------------------------------------------------------
	//Read how the CNT content is measured (volume or weight fraction)
	istringstream istr6(Get_Line(infile));
    istr6 >> nanotube_geo.criterion;
	if(nanotube_geo.criterion=="vol")
	{
		istr6 >> nanotube_geo.volume_fraction;
        
		if(nanotube_geo.volume_fraction>1||nanotube_geo.volume_fraction<0){
            hout << "Error: The volume fraction must be between 0 and 1. Input was: "<<nanotube_geo.volume_fraction<<endl; return 0; }
		hout << "    The CNT volume fraction is "<< nanotube_geo.volume_fraction << endl;
        
		//The total volume of the nanotube network
		nanotube_geo.real_volume = nanotube_geo.volume_fraction*geom_sample.volume;
	}
	else if(nanotube_geo.criterion=="wt")
	{
		istr6 >> nanotube_geo.weight_fraction;
        
		if(nanotube_geo.weight_fraction>1||nanotube_geo.weight_fraction<0){ hout << "Error: The weight fraction must be between 0 and 1. Input was:"<<nanotube_geo.weight_fraction<< endl; return 0; }
		hout << "    The weight fraction is " << nanotube_geo.weight_fraction << endl;
        
        //Read the linear density of a nanotube
		istr6 >> nanotube_geo.linear_density;
		if(nanotube_geo.linear_density<0){
            hout << "Error: The linear density of nanotubes should be non-negetive. Input was: "<<nanotube_geo.linear_density<< endl; return 0; }
        
		//Calculate actual weight of nanotubes
		nanotube_geo.real_weight = nanotube_geo.weight_fraction*geom_sample.volume*geom_sample.matrix_density;
	}
	else { hout << "Error: The content of nanotubes can only be specified in volume (vol) or weight (wt) fraction. Input was: "<<nanotube_geo.criterion<< endl; return 0; }
    
	return 1;
}
//---------------------------------------------------------------------------
//Read geometric parameters for agglomerates
int Input::Read_agg_geo_parameters(struct Agglomerate_Geo &agg_geo, ifstream &infile)
{
	if(agg_geo.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(agg_geo.keywords);
        return 0;
	}
	else agg_geo.mark = true;

	istringstream istr_clust_para(Get_Line(infile));
	istr_clust_para >> agg_geo.vol_fra_criterion;
    //If a negative number is input, it is understood that there are no agglomerates
    //The same happens if a 0 is input, but fur numerical stability it is better to use
    //a negative number instead of 0
    if(agg_geo.vol_fra_criterion<0){
        hout << "    No agglomerates will be generated."<< endl;
    }
    else {
        
        if(agg_geo.vol_fra_criterion>1){
            hout << "Error: The volume fraction of agglomerates must be between 0 and 1. Input was: "<<agg_geo.vol_fra_criterion<< endl; return 0; }
        
        //If the volume fraction of agglomerates is between 0 and 1, read the rest of parameters
        
        istr_clust_para >> agg_geo.amin >> agg_geo.amax;
        if(agg_geo.amin<0||agg_geo.amax<0||agg_geo.amin>agg_geo.amax){
            hout << "Error: The maximum and minimum values for the ellipsoid's (agglomerate's) longest axis must be non-negative and the minimum value must be smaller than the maximum value. Input for minimum was "<< agg_geo.amin<<" and for maximum was "<<agg_geo.amax<<endl; return 0; }
        
        istr_clust_para >> agg_geo.bmin >> agg_geo.cmin;
        if(agg_geo.bmin<0||agg_geo.cmin<0)
        {
            hout << "Error: The minimum values for the ellipsoid's (agglomerate's) middle and shortest axes must be non-negative. Input for middle axis was "<<agg_geo.bmin<<" and for shortest axis was "<<agg_geo.cmin<< endl;
            return 0;
        }
        if(agg_geo.bmin>agg_geo.amin||agg_geo.cmin>agg_geo.amin||agg_geo.cmin>agg_geo.bmin)
        {
            hout << "Error: The minimum values for the ellipsoid's (agglomerate's) middle and shortest axes must be smaller than the minimum value of longest axis,";
            hout<<"and the minimum value of shortest axis must be smaller than the minimum value of middle axis. ";
            hout<<"Input for middle axis was "<<agg_geo.bmin<<", for shortest axis was "<<agg_geo.cmin<<" and the minimum value of longest axis was "<<agg_geo.amin<< endl;
            return 0;
        }
        
        istr_clust_para >> agg_geo.growth_probability;
        if(agg_geo.growth_probability<0||agg_geo.growth_probability>1){
            hout << "Error: The growth probability of nanotubes in agglomerates must be between 0 and 1. Input was: "<<agg_geo.growth_probability<< endl; return 0; }
        
        istr_clust_para >> agg_geo.volf_clust;
        if(agg_geo.volf_clust<0||agg_geo.volf_clust>1) {
            hout << "Error: The volume fraction of nanotubes in agglomerates must be between 0 and 1. Input was: "<<agg_geo.volf_clust<< endl; return 0; }

        istr_clust_para >> agg_geo.print_key;
        if(agg_geo.print_key<0||agg_geo.print_key>2) {
            hout << "Error: The print_key for agglomerates can only be 0, 1 or 2. Input was: "<<agg_geo.print_key<< endl; return 0; }

        //clear ellipsoid vector
        agg_geo.ellips.clear();
    }

	return 1;
}
//---------------------------------------------------------------------------
//Read cutoff distances
int Input::Read_cutoff_distances(struct Cutoff_dist &cutoff_dist, ifstream &infile)
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
	istr >> cutoff_dist.van_der_Waals_dist >> cutoff_dist.tunneling_dist;
    
	return 1;
}
//---------------------------------------------------------------------------
//Read electrical properties of materials
int Input::Read_electrical_parameters(struct Electric_para &electric_para, ifstream &infile)
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
    if (electric_para.applied_voltage <= 0) {
        hout<<"Error: Voltage cannot be negative. Input was: "<<electric_para.applied_voltage<<endl;
        return 0;
    }
    
    //Carbon nanotube resistivity
    istringstream istr1(Get_Line(infile));
    istr1 >> electric_para.resistivity_CNT;
    if (electric_para.resistivity_CNT <= 0) {
        hout<<"Error: CNT resistivity cannot be negative. Input was: "<<electric_para.resistivity_CNT<<endl;
        return 0;
    }
    
    //GNP resistivities
    istringstream istr2(Get_Line(infile));
    istr2 >> electric_para.resistivity_GNP_t >> electric_para.resistivity_GNP_surf;
    //hout << electric_para.resistivity_GNP_t <<" ,"<<electric_para.resistivity_GNP_surf<<endl;
    if (electric_para.resistivity_GNP_t <= 0 || electric_para.resistivity_GNP_surf <= 0) {
        hout<<"Error: GNP resistivities cannot be negative. Input was: "<<electric_para.resistivity_GNP_t<<", "<<electric_para.resistivity_GNP_surf<<endl;
        return 0;
    }
    
    //Polymer matrix resistivity
    istringstream istr3(Get_Line(infile));
    istr3 >> electric_para.resistivity_matrix;
    if (electric_para.resistivity_matrix <= 0) {
        hout<<"Error: Polymer resistivity cannot be negative. Input was: "<<electric_para.resistivity_matrix<<endl;
        return 0;
    }
    //hout << "electric_para.resistivity_matrix = " << electric_para.resistivity_matrix << endl;
    //Electrical constants: electron charge (C), permitivity of vaccum (F/m), CNT work function (V), dielectric constant of polymer
    //istringstream istr3(Get_Line(infile));
    //istr3 >> electric_para.e_charge >> electric_para.e0_vacuum >> electric_para.CNT_work_function >> electric_para.K_polymer;
    
    //Electrical constants:Planck’s constant (m2kg/s), electron charge (C), electron mass (Kg), height barrier (eV)
    istringstream istr4(Get_Line(infile));
    istr4 >> electric_para.h_plank >> electric_para.e_charge >> electric_para.e_mass >> electric_para.lambda_barrier;
    if (electric_para.h_plank<=0||electric_para.e_charge<=0||electric_para.e_mass<=0||electric_para.lambda_barrier<=0) {
        hout<<"Error: Electric parameters cannot be negative. Input was: "<<electric_para.h_plank<<", "<<electric_para.e_charge<<", "<<electric_para.e_mass<<", "<<electric_para.lambda_barrier<<endl;
        return 0;
    }
    //hout << electric_para.h_plank <<", "<< electric_para.e_charge <<", "<< electric_para.e_mass <<", "<< electric_para.lambda_barrier << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//Read flags for Tecplot visualization files
int Input::Read_tecplot_flags(struct Tecplot_flags &tec360_flags, ifstream &infile)
{
	if(tec360_flags.mark)
	{
		//Output a message that the keyword has already been iput
        Warning_message_already_input(tec360_flags.keywords);
        return 0;
	}
	else tec360_flags.mark = true;
    
    //Flag to export generated CNTs:
    // 0: do not export Tecplot files
    // 1: wires
    // 2: singlezone meshes
    // 3: multizone meshes
    istringstream istr0(Get_Line(infile));
	istr0 >> tec360_flags.generated_cnts;
    if (tec360_flags.generated_cnts<0||tec360_flags.generated_cnts>3) {
        hout<<"Error: Flag to export generated CNTs can only be an integer in [0,3]. Input was: "<<tec360_flags.generated_cnts<<endl;
    }
    
    //Flag to export generated GNPs
    // 0: do not export Tecplot files
    // 1: export GNPs as brick elements
    istringstream istr1(Get_Line(infile));
    istr1 >> tec360_flags.generated_gnps;
    if (tec360_flags.generated_gnps<0||tec360_flags.generated_gnps>1) {
        hout<<"Error: Flag to export generated GNPs can only be 0 or 1. Input was: "<<tec360_flags.generated_gnps<<endl;
    }
    
    //Flag to export clusters as obtained from the Hoshen-Kopelman algorithm:
    // 0: do not export Tecplot files
    // 1: wires (3D lines)
    // 2: meshes
    // 3: wires (3D lines) & meshes
    istringstream istr2(Get_Line(infile));
    istr2 >> tec360_flags.clusters;
    if (tec360_flags.clusters<0||tec360_flags.clusters>3) {
        hout<<"Error: Flag to export clusters (from HK76 algorithm) can only be an integer in [0,3]. Input was: "<<tec360_flags.clusters<<endl;
    }

    //Flag to export percolated clusters:
    // 0: do not export Tecplot files
    // 1: wires (3D lines)
    // 2: meshes
    // 3: wires (3D lines) & meshes
    istringstream istr3(Get_Line(infile));
    istr3 >> tec360_flags.percolated_clusters;
    if (tec360_flags.percolated_clusters<0||tec360_flags.percolated_clusters>3) {
        hout<<"Error: Flag to export percolated clusters can only be an integer in [0,3]. Input was: "<<tec360_flags.percolated_clusters<<endl;
    }
    
    //Flag to export the backbone:
    // 0: do not export Tecplot files
    // 1: wires (3D lines)
    // 2: meshes
    // 3: wires (3D lines) & meshes
    istringstream istr4(Get_Line(infile));
    istr4 >> tec360_flags.backbone;
    if (tec360_flags.backbone<0||tec360_flags.backbone>3) {
        hout<<"Error: Flag to export the backbone clusters can only be an integer in [0,3]. Input was: "<<tec360_flags.backbone<<endl;
    }
    
    //Flag to export triangulations if set to 1
    // 0: do not export Tecplot files
    // 1: export triangulations
    istringstream istr5(Get_Line(infile));
    istr5 >> tec360_flags.triangulations;
    if (tec360_flags.triangulations<0||tec360_flags.triangulations>1) {
        hout<<"Error: Flag to export triangulation can only be 0 or 1. Input was: "<<tec360_flags.triangulations<<endl;
    }
    
	return 1;
}
//---------------------------------------------------------------------------
//Readi the geometric parameters of GNPs
int Input::Read_gnp_geo_parameters(struct GNP_Geo &gnp_geo, ifstream &infile)
{
    if(gnp_geo.mark)
    {
        //Output a message that the keyword has already been iput
        Warning_message_already_input(gnp_geo.keywords);
        return 0;
    }
    else gnp_geo.mark = true;

    //----------------------------------------------------------------------
    //Read the type of generation of CNTs on the GNP surface (parallel or independent)
    istringstream istr0(Get_Line(infile));
    istr0 >> gnp_geo.growth_type;
    if (gnp_geo.growth_type != "parallel" && gnp_geo.growth_type != "independent") {
        hout << "Error: The growth type of the CNTs on the GNP surface should be either parallel or independent. Input was: " << gnp_geo.growth_type << endl;
        return 0;
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
        istr1  >> gnp_geo.ini_sita >> gnp_geo.ini_pha;
        if(gnp_geo.ini_sita<0||gnp_geo.ini_sita>PI)
        {
            hout << "Error: The range specified for theta is not in the valid range of (0, PI)." << endl;
            return 0;
        }
        if(gnp_geo.ini_pha<0||gnp_geo.ini_pha>=2*PI)
        {
            hout << "Error: The range specified phi is not in the valid range of (0, 2PI)." << endl;
            return 0;
        }
    }
    
    //----------------------------------------------------------------------
    //Read the step length (in microns) discretization of the GNP
    istringstream istr2(Get_Line(infile));
    istr2 >> gnp_geo.discr_step_length;
    if(gnp_geo.discr_step_length<=0||
       gnp_geo.discr_step_length>=0.25*geom_sample.len_x||
       gnp_geo.discr_step_length>=0.25*geom_sample.wid_y||
       gnp_geo.discr_step_length>=0.25*geom_sample.hei_z) {
        hout << "Error: the step length must be positive and upto 0.25 times the side length of the sample. Input was: "<<gnp_geo.discr_step_length<< endl;	return 0; }
    
    //----------------------------------------------------------------------
    //Read the distribution type (uniform or normal) of the GNP side length (in microns) and maximum and minimum values
    istringstream istr3(Get_Line(infile));
    istr3 >> gnp_geo.size_distrib_type;
    if(gnp_geo.size_distrib_type!="uniform"&&gnp_geo.size_distrib_type!="normal"){
        hout << "Error: The distribution of the length should be either normal or uniform." << endl;	return 0; }
    istr3 >> gnp_geo.len_min >> gnp_geo.len_max;
    if(gnp_geo.len_min<0||gnp_geo.len_max<0||gnp_geo.len_max<gnp_geo.len_min){
        hout << "Error: The maximum and minimum values of GNP side length must be non-negative and the minimum value must be smaller than the maximum value. Input for minimum was "<<gnp_geo.len_min<<" and for maximum was "<<gnp_geo.len_max<< endl; return 0; }
    
    //----------------------------------------------------------------------
    //Define the distribution type (uniform or normal) of the GNP thickness (in microns) and maximum and minimum values
    istringstream istr4(Get_Line(infile));
    istr4 >> gnp_geo.thick_distrib_type;
    if(gnp_geo.thick_distrib_type!="uniform"&&gnp_geo.thick_distrib_type!="normal"){
        hout << "Error: The distribution of the GNP thickness should be either normal or uniform. Input was: "<<gnp_geo.thick_distrib_type<< endl;	return 0; }
    istr4 >> gnp_geo.t_min >> gnp_geo.t_max;
    if(gnp_geo.t_min<0||gnp_geo.t_max<0||gnp_geo.t_max<gnp_geo.t_min) {
        hout << "Error: The maximum and minimum values of GNP thickness must be non-negative and the minimum value must be smaller than the maximum value. Input for minimum was "<<gnp_geo.t_min<<" and for maximum was "<<gnp_geo.t_max<< endl; return 0; }
    
    //----------------------------------------------------------------------
    //Define the CNTs/GNPs mass ratio
    istringstream istr5(Get_Line(infile));
    istr5 >> gnp_geo.mass_ratio;
    if(gnp_geo.mass_ratio <= 0 && simu_para.particle_type != "Hybrid_particles") {
        hout << "Error: The CNT/GNP mass ratio must be a value greater than zero when generating GNPs only." << endl; return 0; }
    
    //----------------------------------------------------------------------
    //Define the CNT density (in gm/cm3)
    istringstream istr6(Get_Line(infile));
    istr6 >> gnp_geo.density;
    if (gnp_geo.density <= 0) {
        hout << "Error: The GNP density has to be a value greater than zero. Input was:" << nanotube_geo.density<< endl;
        return 0;
    }
    
    //----------------------------------------------------------------------
    //Define how the GNP content is measured (volume or weight fraction)
    istringstream istr7(Get_Line(infile));
    istr7 >> gnp_geo.criterion;
    if(gnp_geo.criterion=="vol")
    {
        //Check the particle type
        //If mixed CNTs+GNPs or hybrid particles are used, then the fraction of GNPs needs to be calculated from the CNT fraction
        //Then the CNT fraction needs to be adjusted
        if (simu_para.particle_type == "GNP_CNT_mix" ) {
            
            //Calculate the GNP volume fraction
            gnp_geo.volume_fraction = nanotube_geo.density*nanotube_geo.volume_fraction/(gnp_geo.mass_ratio*gnp_geo.density + nanotube_geo.density);
            
            //Adjust the CNT volume fraction
            nanotube_geo.volume_fraction = nanotube_geo.volume_fraction - gnp_geo.volume_fraction;
            //Adjust the total volume of the nanotube network
            nanotube_geo.real_volume = nanotube_geo.volume_fraction*geom_sample.volume;
            
            hout << "    Given the CNT/GNP mass ratio of " << gnp_geo.mass_ratio << ", the CNT volume fraction was adjusted to " << nanotube_geo.volume_fraction << endl;
            hout << "    The GNP volume fraction is " << gnp_geo.volume_fraction <<endl;
        }
        else if(simu_para.particle_type == "Hybrid_particles") {
            
            //If the mass ratio is positive proceed as with mixed particles
            if (gnp_geo.mass_ratio > Zero) {
                
                //Calculate the GNP volume fraction
                gnp_geo.volume_fraction = nanotube_geo.density*nanotube_geo.volume_fraction/(gnp_geo.mass_ratio*gnp_geo.density + nanotube_geo.density);
                
                //Adjust the CNT volume fraction
                //nanotube_geo.volume_fraction = nanotube_geo.volume_fraction - gnp_geo.volume_fraction;
                
                //Adjust the total volume of the nanotube network
                //nanotube_geo.real_volume = nanotube_geo.volume_fraction*geom_sample.volume;
                
                hout << "    Given the CNT/GNP mass ratio of " << gnp_geo.mass_ratio << ", the CNT volume fraction is " << (nanotube_geo.volume_fraction - gnp_geo.volume_fraction) << endl;
                hout << "    The GNP volume fraction is " << gnp_geo.volume_fraction <<endl;
            }
            //If the mass ratio is negative, then the absolute value is used as CNT density and only calculate the fraction of GNPs
            else {
                //Calculate the mass ratio from the CNT density
                //Use the smallest CNT geometry and largest GNP thickness to minimize mas ratio M
                //A minimum M will result in a maximum volume of GNP
                double M = 2*nanotube_geo.density*PI*nanotube_geo.rad_min*nanotube_geo.rad_min*nanotube_geo.len_min*abs(gnp_geo.mass_ratio);
                M = M/(gnp_geo.density*gnp_geo.t_max);
                
                hout << "Calculated mass ratio is " << M << endl;
                
                //Calculate the GNP volume fraction
                gnp_geo.volume_fraction = nanotube_geo.volume_fraction/M;
            }
            
        }
        else {
            
            istr7 >> gnp_geo.volume_fraction;
            //Check that the volume fraction is between 0 and 1
            if(gnp_geo.volume_fraction>1||gnp_geo.volume_fraction<0) {
                hout << "Error: the volume fraction must be between 0 and 1. Input was: "<<gnp_geo.volume_fraction<< endl;
                return 0;
            }
            
        }
        
        //Calculate the actual volume of GNPs
        gnp_geo.real_volume = gnp_geo.volume_fraction*geom_sample.volume;
    }
    else if(gnp_geo.criterion=="wt")
    {
        
        //Check the particle type
        //If mixed CNTs+GNPs or hybrid particles are used, then the fraction of GNPs needs to be calculated from the CNT fraction
        //Then the CNT fraction needs to be adjusted
        if (simu_para.particle_type == "GNP_CNT_mix" || simu_para.particle_type == "Hybrid_particles") {
            //Calculate the GNP weight fraction
            gnp_geo.weight_fraction = nanotube_geo.weight_fraction/(1+gnp_geo.mass_ratio);
            
            //Adjust the CNT weight fraction
            nanotube_geo.weight_fraction = nanotube_geo.weight_fraction - gnp_geo.weight_fraction;
            //Adjust the real weight of nanotubes
            nanotube_geo.real_weight = nanotube_geo.weight_fraction*geom_sample.volume*geom_sample.matrix_density;
            
            hout << "    Given the CNT/GNP mass ratio of " << gnp_geo.mass_ratio << ", the CNT weight fraction was adjusted to " << nanotube_geo.weight_fraction << endl;
            hout << "    The GNP weight fraction is " << gnp_geo.weight_fraction <<endl;
        } else {
            
            istr7 >> gnp_geo.weight_fraction;
            if(gnp_geo.weight_fraction>1||gnp_geo.weight_fraction<0){
                hout << "Error: the weight fraction must be between 0 and 1." << endl; return 0; }
            hout << "    The weight fraction is " << gnp_geo.weight_fraction << endl;
        }
        
        //Calculate the actual weight of GNPs
        gnp_geo.real_weight = gnp_geo.weight_fraction*geom_sample.volume*geom_sample.matrix_density;
    }
    else {
        hout << "Error: the criterion of generation is neither 'vol' nor 'wt'. Input was: "<< gnp_geo.criterion << endl;
        return 0;
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
