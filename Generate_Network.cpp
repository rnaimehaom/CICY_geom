//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Generate 3D nanoparticle network
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Generate_Network.h"

//NOTE ON PERIODICTY (Jul 16 2017):
//Periodicity was removed as it has not been used in any of the updates of the code since it was uploaded to github
//If periodicity needs to be added back again, then it should deal with the following:
//1) Make the extended domain equal to the composite domain. Periodicity solves the same issue as the extended domain when it comes to homegeneity of CNT distribution.
//2) Calculate the intersection with the boundary. This operation is completely useless with a non-periodic sample so I just deleted it.
//3) When adding the new_cnt vector to the global vector of CNTs, split it into segments. This operation is completely useless with a non-periodic sample and was causing some errors when using the penetrating model so I just deleted it.

//Generate a network of nanoparticles
int Generate_Network::Generate_nanoparticle_network(const Simu_para &simu_para, const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const GNP_Geo &gnp_geo, const Cutoff_dist &cutoffs, const Visualization_flags &vis_flags, const Output_data_flags &out_flags, vector<Point_3D> &points_cnt, vector<double> &radii_out, vector<vector<long int> > &structure, vector<GNP> &gnps)const
{
    //Vector for storing the CNT points
    vector<vector<Point_3D> > cnts_points;
    //Vector for radii, internal variable
    vector<double> radii_in;
    
    double gnp_vol_tot = 0, gnp_wt_tot = 0;
    if (simu_para.particle_type == "CNT_wires") {
        
        //Generate a network defined by points
        if (!Generate_cnt_network_threads_mt(simu_para, geom_sample, nanotube_geo, cutoffs, cnts_points, radii_in)) {
            hout << "Error in generating a CNT network" << endl;
            return 0;
        }
        
    }
    else if (simu_para.particle_type == "CNT_deposit") {
        
        //Sample geometry data structure used for CNT_deposit only
        Geom_sample geom_sample_deposit = geom_sample;
        
        //Set the sample geometry to be the same as the extended domain
        geom_sample_deposit.sample = geom_sample.ex_dom_cnt;
        
        //Generate a network of deposited CNTs
        if (!Generate_cnt_deposit_mt(simu_para, geom_sample, geom_sample_deposit, nanotube_geo, cutoffs, cnts_points, radii_in)) {
            hout << "Error in generating a CNT deposit" << endl;
            return 0;
        }
    }
    else if (simu_para.particle_type == "GNP_cuboids") {
        
        //Generate a GNP network
        vector<vector<int> > sectioned_domain_gnp;
        if (!Generate_gnp_network_mt(simu_para, gnp_geo, geom_sample, cutoffs, gnps, sectioned_domain_gnp, gnp_vol_tot, gnp_wt_tot)) {
            hout << "Error in generating a GNP network" << endl;
            return 0;
        }
        
    }
    else if (simu_para.particle_type == "GNP_CNT_mix") {
        
        //Generate a GNP network
        vector<vector<int> > sectioned_domain_gnp;
        if (!Generate_gnp_network_mt(simu_para, gnp_geo, geom_sample, cutoffs, gnps, sectioned_domain_gnp, gnp_vol_tot, gnp_wt_tot)) {
            hout << "Error in generating a GNP network for mixed particles" << endl;
            return 0;
        }
        
        //Generate a CNT network within a GNP network
        if (!Generate_cnt_network_threads_among_gnps_mt(simu_para, gnp_geo, nanotube_geo, geom_sample, cutoffs, gnps, sectioned_domain_gnp, gnp_vol_tot, gnp_wt_tot, cnts_points, radii_in)) {
            hout<<"Error in generating a CNT network for mixed particles"<<endl;
            return 0;
        }
        
    }
    else if (simu_para.particle_type == "Hybrid_particles") {
        hout << "Hybrid particle network is not implmented get." << endl;
        return 0;
    }
    else {
        hout << "Error: the type of particles should be one of the following: CNT_wires, CNT_deposit, GNP_cuboids, or GNP_CNT_mix. Input value was: " << simu_para.particle_type << endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------
    //Check if visualization files were requested for generated nanoparticles
    if (vis_flags.generated_nanoparticles) {
        
        VTK_Export vtk_exp;
        
        //Export generated CNTs if any
        if (cnts_points.size()) {
            if (vis_flags.generated_nanoparticles == 2) {
                
                //Generate one visualization file per CNT
                for (size_t i = 0; i < cnts_points.size(); i++) {
                    string filename = "cnt_" + to_string(i) + ".vtk";
                    vtk_exp.Export_single_cnt(cnts_points[i], filename);
                }
            }
            else {
                
                //Generate one visualization file with all the CNTs
                vtk_exp.Export_cnts_2D_vector(cnts_points, "cnts_generated.vtk");
            }
        }
        
        //Export generated GNPs if any
        if (gnps.size()) {
            vtk_exp.Export_gnps(gnps, "gnps_generated.vtk");
        }
        
        //Export the sample geometry
        vtk_exp.Export_cuboid(geom_sample.sample, "sample.vtk");
    }
    
    //---------------------------------------------------------------------------
    //If there are CNTS, transform the 2D cnts_points into 1D cpoints and 2D cstructures
    //Also, remove the CNTs in the boundary layer
    if (simu_para.particle_type != "GNP_cuboids") {
        if(Transform_points_cnts(geom_sample, nanotube_geo, cnts_points, points_cnt, radii_in, radii_out, structure)==0) {
            hout<<"Error in Transform_points for CNTs."<<endl; return 0;
        }
        if (!Recalculate_vol_fraction_cnts(geom_sample, simu_para, nanotube_geo, points_cnt, radii_out, structure)) {
            hout<<"Error in Recalculate_vol_fraction_cnts."<<endl; return 0;
        }
        hout<<endl;
    }
    
    //Export generated CNTs after removing the ones that are outside the sample
    if (vis_flags.generated_nanoparticles == 3) {
        VTK_Export vtk_exp;
        vtk_exp.Export_from_cnt_structure(points_cnt, structure, "cnts_generated_trim.vtk");
    }
    
    //Output data files if requested
    if (!Output_data_files(geom_sample, out_flags, points_cnt, radii_out, structure, gnps)) {
        hout<<"Error in Output_data_files."<<endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Generate a network defined by points and connections
//Use the Mersenne Twister for the random number generation
int Generate_Network::Generate_cnt_network_threads_mt(const Simu_para &simu_para, const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radius)const
{
    //Initial seeds, if any are in network_seeds within geom_sample.
    //However, geom_sample cannot be modified, so copy the seeds to a new vector
    unsigned int net_seeds[7];
    if (!CNT_seeds(simu_para.CNT_seeds, net_seeds)) {
        hout<<"Error in CNT_seeds"<<endl;
        return 0;
    }
    
    //Use the seeds generated above
    std::mt19937 engine_x(net_seeds[0]);
    std::mt19937 engine_y(net_seeds[1]);
    std::mt19937 engine_z(net_seeds[2]);
    std::mt19937 engine_phi(net_seeds[3]);
    std::mt19937 engine_theta(net_seeds[4]);
    std::mt19937 engine_rand(net_seeds[5]);
    std::mt19937 engine_initial_direction(net_seeds[6]);
    
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [-1, 1].
    std::uniform_real_distribution<double> dist_initial(-1.0, 1.0);

    //---------------------------------------------------------------------------
    //Time variables to keep track of generation time
    time_t ct0, ct1;
    
    //---------------------------------------------------------------------------
    //Variables for the generated CNT volume and weight
    double vol_sum = 0;
    
    //Variable to count the generated CNT seeds
    int cnt_seed_count = 0;
    
    //Variable to coun the number of rejected CNTs
    int cnt_reject_count = 0;
    //Variable to count the number of ignored CNTs
    int cnt_ignore_count = 0;
    
    //Variable to count the number of times that a point had to be relocated
    int point_overlap_count = 0;
    //Variable to count the number of points that were overlapping other points
    int point_overlap_count_unique = 0;
    
    //---------------------------------------------------------------------------
    //Vectors for handling CNT penetration
    //global_coordinates[i][0] stores the CNT number of global point i
    //global_coordinates[i][1] stores the local point number of global point i
    vector<vector<int> > global_coordinates;
    //sectioned_domain[i] contains all the points in sub-region i.
    //Sub-region i is an overlapping subregion to check for penetrations
    vector<vector<long int> > sectioned_domain;
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    //n_subregions[2] is the number of subregions along z
    int n_subregions[3];
    //Initialize the vector sub-regions
    if (!Initialize_cnt_subregions(geom_sample, n_subregions, sectioned_domain)) {
        hout << "Error in Generate_network_threads_mt when calling Initialize_cnt_subregions" <<endl;
        return 0;
    }
    
    //Get the time when generation started
    ct0 = time(NULL);
    
    //Variable used to check when 10% is completed
    double vol_completed = 0.1;
    //Variable used to store the fraction of CNT volume indicated by vol_completed
    double vol_completed_acc = vol_completed*nanotube_geo.volume;
    
    //Boolean to terminate the main while loop, initialized to false to start the loop
    bool terminate = false;
    
    //---------------------------------------------------------------------------
    while(!terminate)
    {
        //---------------------------------------------------------------------------
        //Vector for a new nanotube
        vector<Point_3D> new_cnt;
        
        //---------------------------------------------------------------------------
        //Randomly generate a CNT length and CNT radius
        double cnt_length, cnt_rad;
        //hout<<"Get_random_value_mt 1"<<endl;
        if (!Get_length_and_radius(nanotube_geo, engine_rand, dist, cnt_length, cnt_rad)) {
            hout << "Error in Generate_network_threads_mt when calling Get_length_and_radius" <<endl;
            return 0;
        }
        
        //Calculate the number of CNT growth steps
        int step_num = (int)(cnt_length/nanotube_geo.step_length) + 1;
        
        //---------------------------------------------------------------------------
        //Cross-sectional area of the current CNT. It is used to calculate the CNT volume.
        const double cnt_cross_area = PI*cnt_rad*cnt_rad;
        
        //---------------------------------------------------------------------------
        //Randomly generate an initial direction, then generate the rotation matrix that results in that rotation
        MathMatrix multiplier(3,3);
        //hout<<"Get_initial_direction_mt"<<endl;
        if (!Get_initial_direction_mt(nanotube_geo.dir_distrib_type, nanotube_geo.ini_theta, nanotube_geo.ini_phi, engine_initial_direction, dist_initial, multiplier)) {
            hout << "Error in Generate_network_threads_mt when calling Get_initial_direction_mt" <<endl;
            return 0;
        }
        
        //---------------------------------------------------------------------------
        //Randomly generate a seed (initial point) of a CNT in the extended domain
        Point_3D new_point;
        //hout<<"Get_seed_point_mt"<<endl;
        if(Get_point_in_cuboid_mt(geom_sample.ex_dom_cnt, new_point, engine_x, engine_y, engine_z, dist)==0) return 0;
        
        //Calculate a cutoff for penetration
        double rad_p_dvdw = cnt_rad + cutoffs.van_der_Waals_dist;
        //Calculate some cutoffs for self-penetration
        double cnt_cutoff = cnt_rad + rad_p_dvdw;
        double cnt_cutoff2 = cnt_cutoff*cnt_cutoff;
        
        //Map to find self-penetrating points
        map<int, vector<int> > subr_point_map;
        //hout<<"Check_penetration step_num="<<step_num<<endl;
        //Check if penetration model is being used
        if (simu_para.penetration_model_flag) {
            
            //Flag for the status of checking penetrations
            int penetration_check = 0;
            
            //Check overlapping of the intial point
            int counter = 1;
            
            //Iteratively search for a seed point that does not penetrate other CNTs
            while (counter <= MAX_ATTEMPTS && penetration_check == 0) {
                
                //Check if the seed point penetrates other CNTs and move it to a valid position
                penetration_check = Check_penetration(geom_sample, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, rad_p_dvdw, cnt_cutoff, cnt_cutoff2, subr_point_map, point_overlap_count, point_overlap_count_unique, new_point);
                
                //Check for errors
                if (penetration_check == -1) {
                    hout<<"Error when attempting to find a valid position for a seed point"<<endl;
                    return 0;
                }
                
                //Check if point could not be placed in a valid position
                else if (penetration_check == 0) {
                    
                    //Point coud not be re-located, so genereta a new point
                    if(!Get_point_in_cuboid_mt(geom_sample.ex_dom_cnt, new_point, engine_x, engine_y, engine_z, dist)) {
                        hout<<"Error when generating a new CNT point after checking for penetrations"<<endl;
                        return 0;
                    }
                }
                
                //Increase the counter
                counter++;
            }
            
            //hout << "Seed deleted" << endl;
            if (counter == MAX_ATTEMPTS && penetration_check == 0) {
                hout<<"Too many attempts to resolve overlapping of an intial CNT point ("<<counter;
                hout<<" attempts). Seed point="<<new_point.str()<<endl;
                return 0;
            }
        }
        
        //Add the CNT seed to the current CNT vector
        new_cnt.push_back(new_point);
        
        //---------------------------------------------------------------------------
        //Count the numer of CNT seeds generated
        cnt_seed_count++;
        //hout << "Seed="<<cnt_seed_count<<endl;
        int max_seed = 1E9;
        if(cnt_seed_count>max_seed)
        {
            hout << "The number of generated seeds is lager than "<<max_seed<<", but the nanotube generation still fails to acheive the requested volume fraction." << endl;
            return 0;
        }
        
        //---------------------------------------------------------------------------
        //hout<<"Growth of CNT"<<endl;
        
        //Get the location of the seed
        bool is_prev_in_sample = Is_point_inside_cuboid(geom_sample.sample, new_cnt[0]);
        
        //Variable to count the number of points of the new CNT that are inside the sample
        int points_in = (is_prev_in_sample)? 1: 0;
        
        //Variable to store the length of the current CNT that is inside the sample
        double cnt_len = 0.0;
        
        //Add the seed point to the overlapping subregions it belongs to using the map
        if (!Add_cnt_point_to_overlapping_regions_map(geom_sample, new_cnt[0], 0, is_prev_in_sample, n_subregions, subr_point_map)) {
            hout<<"Error when adding the seed point to the map of subregions"<<endl;
            return 0;
        }
        
        //Start generation of a CNT siven the generated seed
        for(int i=0; i<step_num; i++)
        {
            //Generates a new CNT point such that the new CNT segment has a random orientation
            //hout << "Get_direction_and_point " << i << endl;
            if (!Get_direction_and_point(nanotube_geo, multiplier, new_point, engine_theta, engine_phi, dist)) {
                hout<<"Error in generating a new CNT point (not a seed) at iteration i="<<i<<endl;
                return 0;
            }
            
            //Check if the new point, cnt_poi, is inside the extended domain
            if(Is_point_inside_cuboid(geom_sample.ex_dom_cnt, new_point))
            {
                //If the point is inside the extended domain, then check for penetration if non-penetrating
                //model is used, or just add it if penetrating model is used
                
                int penetration_check = 1;
                //Check if the penetrating model is used
                if (simu_para.penetration_model_flag) {
                    
                    //Penetration model is used, so check for penetrating points
                    //hout<<"penetration_check i="<<i<<endl;
                    penetration_check = Check_penetration(geom_sample, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, rad_p_dvdw, cnt_cutoff, cnt_cutoff2, subr_point_map, point_overlap_count, point_overlap_count_unique, new_point);
                    
                    //Check for error, which in this case is -1
                    if (penetration_check == -1) {
                        hout<<"Error when checking penetration for a CNT point (not a seed)"<<endl;
                        return 0;
                    }
                }
                
                //Check for penetration again but now to decide whether the point is added
                //or the CNT deleted
                //hout << "Check penetration i="<<i<<endl;
                if (!simu_para.penetration_model_flag || penetration_check == 1) {
                    
                    //---------------------------------------------------------------------------
                    //If the penetrating model is used or if the new point, cnt_poi, was placed
                    //in a valid position (i.e., without penetrating another CNT) then add the point
                    //to the generated CNTs
                    
                    //Calculate the segment length inside the sample and add it to the total CNT length
                    bool is_new_inside_sample;
                    cnt_len = cnt_len + Length_inside_sample(geom_sample.sample, new_cnt.back(), new_point, is_prev_in_sample, is_new_inside_sample);
                    
                    //If the new point, cnt_poi, is inside the sample, then increase the number
                    //of points inside the sample of the new CNT
                    if (is_new_inside_sample) {
                        points_in++;
                    }
                    
                    //For the next iteration of the for loop, cnt_poi will become previous point,
                    //so update the boolean is_prev_in_sample for the next iteration
                    is_prev_in_sample = is_new_inside_sample;
                    
                    //Add the new_point to the overlapping subregions it belongs to using the map
                    //hout << "Add_cnt_point_to_overlapping_regions_map" << endl;
                    if (!Add_cnt_point_to_overlapping_regions_map(geom_sample, new_point, (int)new_cnt.size(), is_new_inside_sample, n_subregions, subr_point_map)) {
                        hout<<"Error when adding a point to the map of subregions"<<endl;
                        return 0;
                    }
                    
                    //Add the new point to the current CNT
                    new_cnt.push_back(new_point);
                    
                    //---------------------------------------------------------------------------
                    //Check if the target volume fraction has been reached
                    if( (vol_sum + cnt_len*cnt_cross_area) >= nanotube_geo.volume) {
                        
                        //Set the terminate variable to true so that the main while-loop is terminated
                        terminate = true;
                        
                        //Break the for-loop so that the current CNT stops its growth
                        break;
                    }
                    
                } else {
                    //---------------------------------------------------------------------------
                    //If the penetrating point could not be accommodated, then delete the current CNT
                    //so that a new one is generated
                    //hout << "Penetrating point could not be accommodated" <<endl;
                    
                    //Clear the new_cnt vector so that it is not added to the rest of CNTs
                    new_cnt.clear();
                    
                    //Increase the count of rejected cnts
                    cnt_reject_count++;
                    
                    //Break the for-loop to terminate the CNT growth
                    break;
                }
                //hout << "done" << endl;
            }
            else {
                
                //If the point is outside the extended domain, break the for-loop
                //to terminate the CNT growth
                break;
            }
            //for-loop ends here
        }
        
        //---------------------------------------------------------------------------
        //Store or ignore the CNT points
        //hout<<"Store_or_ignore_new_cnt"<<endl;
        if (!Store_or_ignore_new_cnt_using_map(simu_para.penetration_model_flag, points_in, cnt_len, cnt_rad, cnt_cross_area, new_cnt, cnts_points, cnts_radius, subr_point_map, sectioned_domain, global_coordinates, vol_sum, cnt_ignore_count)) {
            hout<<"Error when storing or ignoring a new CNT"<<endl;
            return 0;
        }
        
        //Get the time to check progress
        ct1 = time(NULL);
        
        //Check progress
        if (!Check_progress("CNT", (int)(ct1-ct0), nanotube_geo.volume, vol_sum, vol_completed, vol_completed_acc)) {
            hout<<"Error when calculating the percentage of CNT volume generated"<<endl;
            return 0;
        }
        
        //while-loop ends here
    }
    
    //Output the CNT content generated
    if(nanotube_geo.criterion == "wt") {
        
        //Calculate matrix weight
        double matrix_weight = (geom_sample.volume - vol_sum)*geom_sample.matrix_density;
        
        //Calculate the CNT weight
        double cnt_weight = vol_sum*nanotube_geo.density;
        
        hout << endl << "The weight fraction of generated CNTs is: " << cnt_weight/(matrix_weight + cnt_weight);
        hout << ", the target weight fraction was " << nanotube_geo.weight_fraction << endl << endl;

    } else if(nanotube_geo.criterion == "vol") {
        hout << endl << "The volume fraction of generated CNTs was " << vol_sum/geom_sample.volume;
        hout << ", the target volume fraction was " << nanotube_geo.volume_fraction << endl << endl;
    }
    
    hout << "There were " << point_overlap_count_unique << " overlapping points and ";
    hout << point_overlap_count << " overlaps, " << cnt_reject_count << " CNTs were rejected and "<<cnt_ignore_count<<" were ignored." << endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//This function copies the seeds for the random number generators given in the input file
//If seeds were not given, seeds are generated
int Generate_Network::CNT_seeds(const vector<unsigned int> &CNT_seeds, unsigned int net_seeds[])const
{
    //---------------------------------------------------------------------------
    //Random_device is used to generate seeds for the Mersenne twister
    std::random_device rd;

    //If seeds have been specified, copy them
    if (CNT_seeds.size()) {
        for (size_t i = 0; i < CNT_seeds.size(); i++) {
            net_seeds[i] = CNT_seeds[i];
        }
    }
    //If seeds have not been specified, generate them
    else {
        
        //Generate all the new seeds and print them to the output file
        net_seeds[0] = rd();
        //hout << "seed x: "<<net_seeds[0]<<endl;
        net_seeds[1] = rd();
        //hout << "seed y: "<<net_seeds[1]<<endl;
        net_seeds[2] = rd();
        //hout << "seed z: "<<net_seeds[2]<<endl;
        net_seeds[3] = rd();
        //hout << "seed phi: "<<net_seeds[3]<<endl;
        net_seeds[4] = rd();
        //hout << "seed theta: "<<net_seeds[4]<<endl;
        net_seeds[5] = rd();
        //hout << "seed rand: "<<net_seeds[5]<<endl;
        net_seeds[6] = rd();
        //hout << "seed init dir: "<<net_seeds[6]<<endl;
    }
    
    //Output the seeds
    hout<<"CNT seeds:"<<endl;
    for (int i = 0; i < 7; i++) {
        hout<<net_seeds[i]<<' ';
    }
    hout<<endl<<endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//This functions initializes the vectors n_subregions and sectioned_domain
//
//The n_subregions vector is defined to avoid calculating the number of sub-regions for every point when the functions
//Default_region and Add_to_subregion are called. Thus, saving computational time.
//n_subregions[0] is the number of subregions along x
//n_subregions[1] is the number of subregions along y
//n_subregions[2] is the number of subregions along z
//
//The vector sectioned_domain contains the sub-regions to look for overlapping
//It is initialized with the number of sub-regions in the sample
int Generate_Network::Initialize_cnt_subregions(const Geom_sample &sample_geom, int n_subregion[], vector<vector<long int> > &sectioned_domain)const
{
    //Calculate the number of variables along each direction
    //Make sure there is at least one subregion along each direction
    //
    //Number of subregions along x
    n_subregion[0] = max(1, (int)(sample_geom.sample.len_x/sample_geom.gs_minx));
    //Number of subregions along y
    n_subregion[1] = max(1, (int)(sample_geom.sample.wid_y/sample_geom.gs_miny));
    //Number of subregions along z
    n_subregion[2] = max(1, (int)(sample_geom.sample.hei_z/sample_geom.gs_minz));
    
    //Initialize sectioned_domain
    sectioned_domain.assign(n_subregion[0]*n_subregion[1]*n_subregion[2], vector<long int>());
    
    return 1;
}
//Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const
int Generate_Network::Get_length_and_radius(const Nanotube_Geo &nanotube_geo, mt19937 &engine, uniform_real_distribution<double> &dist, double &cnt_length, double &cnt_rad)const
{
    //---------------------------------------------------------------------------
    //Randomly generate a CNT length
    //hout<<"Get_random_value_mt CNT length"<<endl;
    if(Get_random_value_mt(nanotube_geo.len_distrib_type, engine, dist, nanotube_geo.len_min, nanotube_geo.len_max, cnt_length)==0) {
        hout << "Error in Generate_network_threads_mt when calling Get_random_value_mt CNT length" <<endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------
    //Randomly generate a CNT radius
    //hout<<"Get_random_value_mt CNT radius"<<endl;
    if(!Get_random_value_mt(nanotube_geo.rad_distrib_type, engine, dist, nanotube_geo.rad_min, nanotube_geo.rad_max, cnt_rad)) {
        hout << "Error in Generate_network_threads_mt when calling Get_random_value_mt CNT radius" <<endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Check if the current CNT is penetrating another CNT, i.e. is the new point is overlapping other point
//1: a) No penetration
//   b) No need to check for penetration (point is in boundary layer or there are no other points in the same sub-region)
//   c) There was penetration but it was succesfully resolved
//0: There was penetration but could not be resolved
int Generate_Network::Check_penetration(const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<double> &radii, const vector<Point_3D> &cnt_new, const int n_subregions[], const double &rad_p_dvdw, const double &cnt_cutoff, const double &cnt_cutoff2, const map<int, vector<int> > &subr_point_map, int &point_overlap_count, int &point_overlap_count_unique, Point_3D &new_point)const
{
    //Get the sub-region the point belongs to
    //hout<<"Get_subregion 0"<<endl;
    int subregion = Get_cnt_point_subregion(geom_sample, n_subregions, new_point);
    
    //If the sub-region is -1, then the point is in the boundary layer, so there is no need to check penetration
    if (subregion == -1) {
        return 1;
    }
    //This vector will store the coordintes of the points that the input "point" is penetrating
    vector<Point_3D> affected_points;
    //This vector stores the distance at which the two points should be
    vector<double> cutoffs_p;
    //This vector stores the distance at which the two points actually are
    vector<double> distances;
    
    //Iteration 0:
    //Check if there are any penetrations in the corresponding sub-region
    //hout<<"Get_penetrating_points 0"<<endl;
    Get_penetrating_points(cnts, global_coordinates, sectioned_domain[subregion], radii, rad_p_dvdw, new_point, affected_points, cutoffs_p, distances);
    //hout<<"1 affected_points.size="<<affected_points.size()<<endl;
    
    //Check if there are any penetration within the CNT
    //hout<<"Get_penetrating_points_within_cnt 0"<<endl;
    //hout<<"point="<<point.str()<<endl;
    Get_penetrating_points_within_cnt(subregion, cnt_cutoff, cnt_cutoff2, new_point, cnt_new, subr_point_map, affected_points, cutoffs_p, distances);
    //hout<<"2 affected_points.size="<<affected_points.size()<<endl;
    
    //Update the counter of overlapping points only when an overlapping point was found the first time
    if (affected_points.size()) {
        point_overlap_count_unique++;
    }
    
    //I move the point up to max_attempts times. If there is still penetration then I delete it
    for (int attempts = 0; attempts < MAX_ATTEMPTS; attempts++) {
        
        //--------------------------------------------------------------------------------------------
        //Check if there are any penetrating points
        if (affected_points.size()) {
            //Update the counter of overlaps
            point_overlap_count++;
            
            //Find the new point
            //hout << "Point " << global_coordinates.size()-1+cnt_new.size() << " in CNT " << cnts.size() << " is overlapping." <<endl;
            //hout << "Moved a point from initial position "<<new_point.str()<<" attempts="<<attempts<<endl;
            //Check if seed point
            if (cnt_new.empty()) {
                
                //Point is a seed, so then use functions that move the point
                //hout << "Move_point" << endl;
                Move_point(cutoffs_p, distances, affected_points, new_point);
            }
            else {
                
                //If point is not a seed, use functions that rotate the CNT segment
                //hout << "Move_point_by_totating_cnt_segment" << endl;
                if (!Move_point_by_totating_cnt_segment(nanotube_geo.step_length, cnt_new, cutoffs_p, distances, affected_points, new_point)) {
                    hout<<"Error in Check_penetration when calling Move_point_by_totating_cnt_segment"<<endl;
                    return -1;
                }
            }
            //hout << "Moved a point to final position "<<new_point.str()<<endl;
            
            //Check that the new point is within the permited orientation repect to the previous segment
            //hout<<"Check_segment_orientation"<<endl;
            if (!Check_segment_orientation(new_point, cnt_new)) {
                //hout << "Deleted CNT number " << cnts.size() << " of size " << cnt_new.size();
                //hout << " (the point is not in a valid orientation)" << endl;
                //When not in a valid position it cannot be moved again so a new CNT is needed
                return 0;
            }
            
            //Need to update point sub-region as it could be relocated to a new sub-region
            //hout<<"Get_subregion 1 point="<< new_point.str()<<endl;
            subregion = Get_cnt_point_subregion(geom_sample, n_subregions, new_point);
            //Check if after moving the point it is now in the boundary layer
            if (subregion == -1) {
                //If the point is now in the boundary layer, terminate the function
                //there is no need to continue checking
                return 1;
            }
            
            //Need to clear the vectors affected_points, contact_coordinates and temporal_contacts so they are used again with the new point
            affected_points.clear();
            cutoffs_p.clear();
            distances.clear();
            
            //Check if there are any penetrations in the corresponding sub-region
            //for the point in its new position
            //hout<<"Get_penetrating_points sectioned_domain.size="<<sectioned_domain.size()<<" sectioned_domain["<<subregion<<"].size="<<sectioned_domain[subregion].size()<<endl;
            Get_penetrating_points(cnts, global_coordinates, sectioned_domain[subregion], radii, rad_p_dvdw, new_point, affected_points, cutoffs_p, distances);
            //hout<<"3 affected_points.size="<<affected_points.size()<<endl;
            
            //Check if there are any penetration within the CNT
            //hout<<"Get_penetrating_points_within_cnt"<<endl;
            //hout<<"point="<<point.str()<<endl;
            Get_penetrating_points_within_cnt(subregion, cnt_cutoff, cnt_cutoff2, new_point, cnt_new, subr_point_map, affected_points, cutoffs_p, distances);
            //hout<<"4 affected_points.size="<<affected_points.size()<<endl;
        } else {
            
            //If the size of affected_points is zero there are no penetrating points
            //Terminate the function with 1 to indicate succesful relocation of point
            return 1;
        }
        //hout<<"for end"<<endl;
    }
    
    //If the last iteration is reached and there are still affected points
    //then the point could not be accommodated, so terminate with 0, i.e.,
    //it will terminate with 0 when affected_points is not empty
    return affected_points.empty();
}
//---------------------------------------------------------------------------
//This function returns the subregion a point belongs to
int Generate_Network::Get_cnt_point_subregion(const Geom_sample &geom_sample, const int n_subregions[], const Point_3D &point)const
{
    if (Is_point_inside_cuboid(geom_sample.sample, point)) {
        
        //These variables will give me the region cordinates of the region that a point belongs to
        int a, b, c;
        
        //Calculate the subregion-coordinates
        Get_subregion_coordinates(geom_sample, n_subregions, point, a, b, c);
        
        //Return the subregion number
        return Calculate_t(a, b, c, n_subregions[0], n_subregions[1]);
    } else {
        //If the point is in the boundary layer, then there is no need to calculate its sub-region
        return -1;
    }
}
//This function finds the coordinates of a subregions where point belongs
void Generate_Network::Get_subregion_coordinates(const Geom_sample &geom_sample, const int n_subregions[], const Point_3D &point, int &a, int &b, int &c)const
{
    //Calculate the subregion-coordinates
    a = (int)((point.x-geom_sample.sample.poi_min.x)/geom_sample.gs_minx);
    //Limit the value of a as it has to go from 0 to n_subregions[0]-1
    if (a == n_subregions[0]) a--;
    b = (int)((point.y-geom_sample.sample.poi_min.y)/geom_sample.gs_miny);
    //Limit the value of b as it has to go from 0 to n_subregions[1]-1
    if (b == n_subregions[1]) b--;
    c = (int)((point.z-geom_sample.sample.poi_min.z)/geom_sample.gs_minz);
    //Limit the value of c as it has to go from 0 to n_subregions[2]-1
    if (c == n_subregions[2]) c--;
}
//---------------------------------------------------------------------------
//This function iterates over a sub-region and determines if there are any penetrating points
//If there are penetrating points, they are stored in the vector affected_points
void Generate_Network::Get_penetrating_points(const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<long int> &subregion_vec, const vector<double> &radii, const double &rad_plus_dvdw, Point_3D &point, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const
{
    //They are just intermediate variables and I only use them to make the code more readable
    long int coord2;
    int P2, CNT2;
    //cutoff_p is used for the cutoff between two points (of different CNTs), distance is the actual distance
    //between those two points
    double cutoff_p, distance;
    
    //hout << "Check0 " << subregion_vec.size() << ' ';
    //hout << "CNT1=" << point.flag <<" ("<< point.x<<", "<<point.y<< ", " <<point.z<<") r1=" << cnt_rad<<endl;
    for (long int i = 0; i < (long int)subregion_vec.size(); i++) {
        
        //hout << "Check1 " ;
        coord2 = subregion_vec[i];
        //hout << "coord2="<<coord2<<" global_coordinates.size()="<<global_coordinates.size()<<endl;
        
        CNT2 = global_coordinates[coord2][0];
        P2 = global_coordinates[coord2][1];
        
        //hout << "Check4 ";
        cutoff_p = rad_plus_dvdw + radii[CNT2];
        
        //Check if the second point is in the cube of size 2cutoff_p and centered in P1
        //This is easier and faster to check than calculating the distance from poin to point every time
        if ( (cnts[CNT2][P2].x<point.x+cutoff_p)&&(cnts[CNT2][P2].x>point.x-cutoff_p)&&(cnts[CNT2][P2].y<point.y+cutoff_p)&&(cnts[CNT2][P2].y>point.y-cutoff_p)&&(cnts[CNT2][P2].z<point.z+cutoff_p)&&(cnts[CNT2][P2].z>point.z-cutoff_p) ) {
            
            distance = point.distance_to(cnts[CNT2][P2]);
            
            //If it is inside the cube, then it is worth to take the time to calculate the distance from point ot point
            if (distance < cutoff_p) {
                affected_points.push_back(cnts[CNT2][P2]);
                cutoffs_p.push_back(cutoff_p);
                distances.push_back(distance);
                //hout << "Penetrating point"<< affected_points.size()<< " CNT2=" << CNT2 << "("<<cnts[CNT2].size()<<" points) P2=" << P2 << " r2=" << radii[CNT2] << " (" << cnts[CNT2][P2].x << ", " << cnts[CNT2][P2].y << ", " << cnts[CNT2][P2].z << ") "<< endl;
            }
        }
        //hout << "Check5 ";
    }
    //hout << "Check6 " << endl;
}
//---------------------------------------------------------------------------
//This function checks if there are points within the same CNT that may be penetrating
//NOTE: cutoff = (2*cnt_rad + d_vdw), cutoff2 = (2*cnt_rad + d_vdw)^2
void Generate_Network::Get_penetrating_points_within_cnt(const int &subregion, const double &cnt_cutoff, const double &cnt_cutoff2, const Point_3D &point, const vector<Point_3D> &cnt_new, const map<int, vector<int> > &subr_point_map, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances)const
{
    //Search for the subregion in the map and save it in the iterator variable
    map<int, vector<int> >::const_iterator it = subr_point_map.find(subregion);
    
    //Check if:
    //there is a map from subregion to a vector of point numbers
    //AND
    //there are points within the CNT in the subregion (i.e., elements in the vector)
    if (it != subr_point_map.end() && subr_point_map.at(subregion).size()) {
        
        //Calculate the index of the previous to last point in the CNT
        int prev_last_idx = (int)cnt_new.size() - 2;
        
        //There are points in the sub-region, so check if any of them penetrate
        //the newly generated point
        for (size_t i = 0; i < subr_point_map.at(subregion).size(); i++) {
            
            //Get the point number within cnt_new
            int Pi = subr_point_map.at(subregion)[i];
            
            //Ignore the last two points in the CNT as both of them are always below the
            //cutoff distance
            if (Pi < prev_last_idx) {
                
                //Calculate the squared distance to point
                double dist = point.squared_distance_to(cnt_new[Pi]);
                
                //Check if point is below the cutoff
                if (dist + Zero - cnt_cutoff2< Zero) {
                    
                    //point is too close to point Pi, so add it to the corresponding vectors
                    affected_points.push_back(cnt_new[Pi]);
                    cutoffs_p.push_back(cnt_cutoff);
                    distances.push_back(sqrt(dist));
                    //hout<<"aff_pt"<<cnt_new[Pi].str()<<" dist="<<sqrt(dist)<<" cnt_cutoff="<<cnt_cutoff<<endl;
                    //hout<<"Pi="<<Pi<<" dist="<<dist<<" cnt_cutoff2="<<cnt_cutoff2<<endl;
                    //hout<<"dist + Zero - cnt_cutoff2="<<dist + Zero - cnt_cutoff2<<endl;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------
//This function rotates a CNT segment according to the number of points it is penetrating
//The "rotation" is achieved by moving the point to a location where it does not
//penetrate any point anymore and keeps all (cutoff) distances
int Generate_Network::Move_point_by_totating_cnt_segment(const double &step, const vector<Point_3D> &cnt_new, const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &new_point)const
{
    //The number of overlapings will determine how the new point is moved
    if (affected_points.size() == 1){
        
        //Handle the case with only one penetrating point
        //hout << "ovelappings == 1"<<endl;
        if (!One_penetrating_point(step, cutoffs[0], affected_points[0], cnt_new.back(), new_point)) {
            hout<<"Error in Move_point_by_totating_cnt_segment when calling One_penetrating_point"<<endl;
            return 0;
        }
    } /*else if (affected_points.size() == 2){
        
        //Handle the case with two penetrating points
        hout << "ovelappings == 2"<<endl;
        if (!Two_penetrating_points(step, cutoffs, affected_points, cnt_new.back(), new_point)) {
            hout<<"Error in Move_point_by_totating_cnt_segment when calling Two_penetrating_points"<<endl;
            return 0;
        }
    } */
    else if (affected_points.size() >= 2) {
        
        //Handle the case with three penetrating points
        //This actually finds the two closest point and calls the function
        //that handles the case with two penetrating points
        /*hout << "ovelappings >= 3 (ovelappings="<<affected_points.size()<<")"<<endl;
        if (!Three_or_more_penetrating_points(step, cutoffs, distances, affected_points, cnt_new.back(), new_point)) {
            hout<<"Error in Move_point_by_totating_cnt_segment when calling Three_or_more_penetrating_points"<<endl;
            return 0;
        }*/
        //hout << "ovelappings >= 2 (ovelappings="<<affected_points.size()<<")"<<endl;
        if (!Two_or_more_penetrating_points(step, cutoffs, distances, affected_points, cnt_new.back(), new_point)) {
            hout<<"Error in Move_point_by_totating_cnt_segment when calling Two_or_more_penetrating_points"<<endl;
            return 0;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//This function finds the new location for a point when it penetrates only one point
int Generate_Network::One_penetrating_point(const double &d_s, const double &d_c, const Point_3D &S, const Point_3D &P0, Point_3D &N)const
{
    //Calculate unit vector from last point in CNT (P0) towars penetrating point (S)
    Point_3D u_hat = (S - P0).unit();
    
    //Calculate unit vector from segment P0S towards point N' (new location of N)
    Point_3D v_hat = ((u_hat.cross(N - P0)).cross(u_hat)).unit();
    
    //Calculate distance of segment P0S
    double d1 = P0.distance_to(S);
    
    //This variable is used twice
    double d_s2 = d_s*d_s;
    
    //Calculate dstance from S to projection of P0N' on segment P0S
    double l1 = (d_s2 + d1*d1 - d_c*d_c)*0.5/d1;
    //hout<<"d_s="<<d_s<<" d1="<<d1<<" d_c="<<d_c<<endl;
    //hout<<"d_s2="<<d_s2<<" d1*d1="<<d1*d1<<" d_c*d_c="<<d_c*d_c<<endl;
    
    //Calculate distance from segment P0S to point N'
    double h = sqrt(d_s2 - l1*l1);
    //hout<<"l1="<<l1<<" h="<<h<<endl;
    
    //Calculate the new location of N
    N = P0 + u_hat*l1 + v_hat*h;
    
    return 1;
}
//---------------------------------------------------------------------------
//This function finds the two closest points and calls the function that moves a point that
//penetrates two other points
//When a point overlaps three or more points, it becomes too difficult to find the new location
//Hence, this function that finds the two closest points
//The two closest point are chosen since those would be the more critical ones
int Generate_Network::Two_or_more_penetrating_points(const double &d_s, const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, const Point_3D &P0, Point_3D &N)const
{
    //Use the distances vector to find the two closest points
    //Index of the element with the closest point
    //Initialize with the index of the first affected point
    int i1 = 0;
    
    //Distance to the closest point
    //Initialize with the distance to the first affected point
    double d1 = distances[0];
    
    //Scan the rest of the affected points
    for (int i = 1; i < (int)distances.size(); i++) {
        
        //First check if the distances[i] is smaller than d2
        if (distances[i] - d1 < Zero) {
            
            //Update the new closest distance
            d1 = distances[i];
            
            //Update the index of the closest point
            i1 = i;
        }
    }
    
    //Now, call the function that moves a point that penetrates one point
    //hout<<"Reduced to One_penetrating_point"<<endl;
    if (!One_penetrating_point(d_s, cutoffs[i1], affected_points[i1], P0, N)) {
        hout<<"Error in Two_or_more_penetrating_points when calling One_penetrating_point"<<endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//This function finds the new location for a point when it penetrates two points
int Generate_Network::Two_penetrating_points(const double &d_s, const vector<double> &cutoffs, const vector<Point_3D> &affected_points, const Point_3D &P0, Point_3D &N)const
{
    //hout<<"affected_points.size="<<affected_points.size()<<endl;
    //hout<<"affected_points[0]="<<affected_points[0].str()<<endl;
    //hout<<"affected_points[1]="<<affected_points[1].str()<<endl;
    //hout<<"P0="<<P0.str()<<endl;
    //Calculate length of segment P0S1 and its square
    double d1 = P0.distance_to(affected_points[0]);
    double d1_2 = d1*d1;
    hout<<"d1="<<d1<<endl;
    
    //Calculate length of segment P0S2 and its square
    double d2 = P0.distance_to(affected_points[1]);
    double d2_2 = d2*d2;
    hout<<"d2="<<d2<<endl;
    hout<<"d12="<<affected_points[0].distance_to(affected_points[1])<<endl;
    hout<<"d_s="<<d_s<<" cutoffs[0]="<<cutoffs[0]<<" cutoffs[1]="<<cutoffs[1]<<endl;
    
    //Calculate squrare of d_s
    double ds_2 = d_s*d_s;
    
    //Calculate x-coordinate of N' (new location of N) in simplified geometry
    double x = (ds_2 + d1_2 - cutoffs[0]*cutoffs[0])*0.5/d1;
    hout<<"x="<<x<<endl;
    
    //Calculate unit vector in the direction from P0 to S2
    Point_3D w_hat = (affected_points[1] - P0).unit();
    hout<<"w_hat="<<w_hat.str()<<endl;
    
    //Calculate unit vector in the direction from P0 to S1
    Point_3D u_hat = (affected_points[0] - P0).unit();
    //hout<<"P0="<<P0.str()<<" N="<<N.str()<<endl;
    //hout<<"S1="<<affected_points[0].str()<<" S2="<<affected_points[1].str()<<endl;
    hout<<"u_hat="<<u_hat.str()<<endl;
    hout<<"cos(alpha)="<<u_hat.dot(w_hat)<<endl;
    
    //Calculate signed distance l1
    double l1 = d2*u_hat.dot(w_hat);
    
    //Calculate l2
    double l2 = sqrt(d2_2 - l1*l1);
    hout<<"l1="<<l1<<" l2="<<l2<<" d2_check="<<sqrt(l1*l1 + l2*l2)<<endl;
    
    //Calculate y-coordinate of N' (new location of N) in simplified geometry
    double y = (ds_2 + d2_2 - 2*l1*x - cutoffs[1]*cutoffs[1])*0.5/l2;
    hout<<"y="<<y<<" = ("<<ds_2<<" + "<<d2_2<<" - "<<2*l1*x<<" - "<<cutoffs[1]*cutoffs[1]<<")/"<<2*l2<<endl;
    hout<<"y= "<<(ds_2 + d2_2 - cutoffs[1]*cutoffs[1])/(2*l2)<<" - "<<l1/l2<<"x"<<endl;
    
    //Calculate z-coordinate of N' (new location of N) in simplified geometry
    //This z-coordinate is the height of the tetrahedron
    double h = sqrt(ds_2 - x*x - y*y);
    hout<<"ds_2 - x*x - y*y="<<ds_2<<" - "<<x*x<<" - "<<y*y<<" = "<<ds_2 - x*x - y*y<<endl;
    hout<<"h="<<h<<endl;
    
    //Unit vector from triangle P0S1S2 towards N' (new location of N)
    Point_3D v_hat = u_hat.cross(w_hat);
    //Make sure v_hat goes towards N (and thus also N')
    if (v_hat.dot(N - P0) < Zero) {
        v_hat = v_hat*(-1);
    }
    hout<<"v_hat="<<v_hat.str()<<endl;
    
    //Calculate unit vector from segment P0S1 towards S2
    Point_3D r_hat = (u_hat.cross(w_hat)).cross(u_hat);
    hout<<"r_hat="<<r_hat.str()<<endl;
    
    //Calculate new location of N (i.e., N')
    N = u_hat*x + r_hat*y + v_hat*h;
    
    return 1;
}
//---------------------------------------------------------------------------
//This function finds the two closest points and calls the function that moves a point that
//penetrates two other points
//When a point overlaps three or more points, it becomes too difficult to find the new location
//Hence, this function that finds the two closest points
//The two closest point are chosen since those would be the more critical ones
int Generate_Network::Three_or_more_penetrating_points(const double &d_s, const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, const Point_3D &P0, Point_3D &N)const
{
    //Use the distances vector to find the two closest points
    //i1 and i2 will be the indices of the closest affected_points,
    //they will be initialized with the first two
    //i1 will have the index of the closest point, while i2 will have the
    //index of the second closest point
    int i1, i2;
    
    //d1 and d2 will be the distances to the closest affected_points,
    //they will be initialized with the first two
    //d1 will have the distance to the closest point,
    //while d2 will have the distance to the second closest point
    double d1, d2;
    
    //Sort the first two distances
    if (distances[0] - distances[1] < Zero) {
        i1 = 0;
        d1 = distances[0];
        i2 = 1;
        d2 = distances[1];
    } else {
        i1 = 1;
        d1 = distances[1];
        i2 = 0;
        d2 = distances[0];
    }
    
    //Once the closest points have been initialized, then scan the rest of the affected points
    for (int i = 2; i < (int)distances.size(); i++) {
        
        //First check if the distances[i] is smaller than d2
        if (distances[i] < d2) {
            
            //In this case, distances[i] is one of the two closest points
            //Now I need to check against d1 in case distances[i] is the closest point now
            if (distances[i] < d1) {
                
                //distances[i] is the closest point
                //update the distances
                d2 = d1;
                d1 = distances[i];
                
                //update the indices
                i2 = i1;
                i1 = i;
            } else {
                
                //distances[i] is the second closest point
                //update the distances
                d2 = distances[i];
                i2 = i;
            }
        }
    }
    
    //Generate the necessary vectors for the case of two overlapping points
    vector<Point_3D> two_affected_points(2);
    two_affected_points[0] = affected_points[i1];
    two_affected_points[1] = affected_points[i2];
    vector<double> two_cutoffs(2);
    two_cutoffs[0] = cutoffs[i1];
    two_cutoffs[1] = cutoffs[i2];
    
    //Now, call the function that moves a point that penetrates two points
    hout<<"Two_penetrating_points"<<endl;
    if (!Two_penetrating_points(d_s, two_cutoffs, two_affected_points, P0, N)) {
        hout<<"Error in Three_or_more_penetrating_points when calling Two_penetrating_points"<<endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//This function moves a point according to the number of points it is overlapping
void Generate_Network::Move_point(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const
{
    //The number of overlapings will determine how the new point is moved
    //However, first I need to eliminate invalid-points, which are the points of perfect overlapping, i.e.
    //points in the exact same location
    int overlappings = (int)affected_points.size();
    if (overlappings == 1){
        //Handle the case with only one overlapping point
        //hout << "ovelappings == 1"<<endl;
        One_overlapping_point(cutoffs, distances, affected_points, point);
    } else if (overlappings == 2){
        //Handle the case with two overlapping points
        //hout << "ovelappings == 2"<<endl;
        Two_overlapping_points(cutoffs, affected_points, point);
    } else if (overlappings >= 3) {
        //Handle the case with three overlapping points
        //This actually finds the two closest point and calls the function
        //that handles the case with two overlapping points
        //hout << "ovelappings >= 3"<<endl;
        Three_or_more_overlapping_points(cutoffs, distances, affected_points, point);
    }
}
//---------------------------------------------------------------------------
//This function finds the new location for an overlapping point when it overlaps only one point
void Generate_Network::One_overlapping_point(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const
{
    //When there is overlapping with one point only, then this is the simplest and easiest case
    //Just move the point in the direction form P_old to P_new a distance cutoff from P_old
    
    //Check that point is not the same as the overlapping point
    if (distances[0] < Zero) {
        
        //The points are the same or too close to be at the same location
        
        //Create a random direction
        Point_3D rand_dir = Point_3D( (double)(rand()%200 - 99), (double)(rand()%200 - 99), (double)(rand()%200 - 99));
        rand_dir.make_unit();
        
        //Move it into a random direction
        point = affected_points[0] + rand_dir*(cutoffs[0]+Zero);
    }
    //Proceed as usual
    else {
        
        //Calculate direction unit vector
        //distances[0] already has the lenght of the vector form point to P
        Point_3D direction = (point - affected_points[0])/distances[0];
        
        //So the new point is P_old + d*n. d is the cutoff and n the unit vector
        //The Zero is to avoid machine precision errors. Without it, when comparing
        //the new point with the other points in the same region, the program was judging
        //them to be below the cutoff for the van der Waals distance. Even though they were
        //in the limit. After adding this Zero that issue was eliminated
        point = affected_points[0] + direction*(cutoffs[0]+Zero);
    }
}
//---------------------------------------------------------------------------
//This function finds the new location for an overlapping point when it overlaps two points
void Generate_Network::Two_overlapping_points(const vector<double> &cutoffs, const vector<Point_3D> &affected_points, Point_3D &new_point)const
{
    //Calculate vectors AP and AB
    //P is new_point
    //A is affected_points[0]
    //B is affected_points[1]
    Point_3D AP = new_point - affected_points[0];
    Point_3D AB = affected_points[1] - affected_points[0];
    
    //Calculate length of AB
    double d_AB = AB.length();
    
    //Calcuate distance from A to midpoint of segment AB, i.e., point Q
    double dp = AP.dot(AB)/d_AB;
    
    //Calculate location of point Q
    Point_3D Q = affected_points[0] + AB*dp/d_AB;
    
    //Calculate distance from Q to P
    //Instead of using the length of AP, I need the length of AP_new, i.e., the length of
    //the vector that goes from A to P in its new position
    //I don't know the new position yet, but I know the distance which is stored in cutoffs[0]
    //which is the one I use here
    double d = sqrt(cutoffs[0]*cutoffs[0] - dp*dp);
    
    //Calculate vector from segment AB towards point P
    Point_3D N = (AB.cross(AP)).cross(AB);
    //Make N a unitary vector
    N.make_unit();
    
    //Calculate position of new_point
    new_point = Q + N*(d + Zero);
    //The Zero is to avoid machine precision errors. Without it, when comparing
    //the new point with the other points in the same region, the program was
    //judging them to be below the cutoff for the van der Waals distance.
    //Even though they were in the limit. After adding this Zero that issue was eliminated
}
//---------------------------------------------------------------------------
//This function finds the two closest points and calls the function that moves a point that overlaps two other points
//When a point overlaps three or more points, it becomes too difficult to find the new location
//Hence, this function that finds the two closest points
//The two closest point are chosen since those would be the more critical ones
void Generate_Network::Three_or_more_overlapping_points(const vector<double> &cutoffs, const vector<double> &distances, const vector<Point_3D> &affected_points, Point_3D &point)const
{
    //Use the distances vector to find the two closest points
    //i1 and i2 will be the indices of the closest affected_points, they will be initialized with the first two
    //i1 will have the index of the closest point, while i2 will have the index of the second closest point
    int i1, i2;
    //d1 and d2 will be the distances to the closest affected_points, they will be initialized with the first two
    //d1 will have the distance to the closest point, while d2 will have the distance to the second closest point
    double d1, d2;
    //Sort the first two distances
    if (distances[0] < distances[1]) {
        i1 = 0;
        d1 = distances[0];
        i2 = 1;
        d2 = distances[1];
    } else {
        i1 = 1;
        d1 = distances[1];
        i2 = 0;
        d2 = distances[0];
    }
    
    //Once the closest points have been initialized, then scan the rest of the affected points
    for (int i = 2; i < (int)distances.size(); i++) {
        //First check if the distances[i] is smaller than d2
        if (distances[i] < d2) {
            //In this case, distances[i] is one of the two closest point
            //Now I need to check against d1 in case distances[i] is the closest point now
            if (distances[i] < d1) {
                //distances[i] is the closest point
                //update the distances
                d2 = d1;
                d1 = distances[i];
                //update the indices
                i2 = i1;
                i1 = i;
            } else {
                //distances[i] is the second closest point
                //update the distances
                d2 = distances[i];
                i2 = i;
            }
        }
    }
    
    //Generate the necessary vectors for the case of two overlapping points
    vector<Point_3D> two_affected_points(2);
    two_affected_points[0] = affected_points[i1];
    two_affected_points[1] = affected_points[i2];
    vector<double> two_cutoffs(2);
    two_cutoffs[0] = cutoffs[i1];
    two_cutoffs[1] = cutoffs[i2];
    
    //Now, call the function that moves a point that overlaps two points
    Two_overlapping_points(two_cutoffs, two_affected_points, point);
}
//---------------------------------------------------------------------------
//This function checks that "point" is within the bounds of the segment orientation.
//The criterion is just checking the point is not more than pi/2 respect with the previous
//segment. In the limiting case we have a straight triangle. So I calculate the hypotenuse.
//I also measure the distance between "point" and the second before that.
//If the distance  between points is less than the hypotenuse, then it has an
//incorrect orientation
int Generate_Network::Check_segment_orientation(const Point_3D &point, const vector<Point_3D> &cnt_new)const
{
    //If at least two points have already been generated, then check if the new point has a valid orientation
    if (cnt_new.size()>=2) {
        int last = (int)cnt_new.size()-1;
        
        //From the cosines law:
        //c^2 = a^2 + b^2 -2*a*b*cos(g)
        //where g is the angle between two consecutive segments, i.e., the angle between a and b
        
        //A valid angle is g > pi/2, so from the cosines law:
        //2*a*b*cos(g) = a^2 + b^2 - c^2 < 0 when g > pi/2
        //2*a*b*cos(g) = a^2 + b^2 - c^2 > 0 when g < pi/2
        
        //Thus an invalid angle happens when a^2 + b^2 - c^2 > 0
        //Check if this is the case
        
        double a2 = point.squared_distance_to(cnt_new.back());
        double b2 = cnt_new.back().squared_distance_to(cnt_new[last-1]);
        double c2 = point.squared_distance_to(cnt_new[last-1]);
        
        double angle_discriminant = a2 + b2 - c2;
        
        if (angle_discriminant > Zero) {
            //The point is not in a valid position
            return 0;
        }
        else {
            //The point is in a valid position
            return 1;
        }
    }
    else {
        //If point is the first or second point, its orientation does not matter
        return 1;
    }
}
//This function generates the new point of a CNT, such that the new segment has a random orientation
int Generate_Network::Get_direction_and_point(const Nanotube_Geo &nanotube_geo, MathMatrix &multiplier, Point_3D &cnt_poi, mt19937 &engine_theta, mt19937 &engine_phi, uniform_real_distribution<double> &dist)const
{
    //Randomly generate a direction in the spherical coordinates
    //To have the positive Z-axis to be a central axis
    //Then, the radial angle, theta, obeys a normal distribution (theta \in fabs[(-omega,+omega)]) and the zonal angle, phi, obeys a uniform distribution (phi \in (0,2PI))
    double cnt_theta = 0, cnt_phi = 0;
    if(!Get_direction_normal_distribution(nanotube_geo, cnt_theta, cnt_phi, engine_theta, engine_phi, dist)){
        hout<<"Error in Get_direction_and_point when calling Get_direction_normal_distribution"<<endl;
        return 0;
    }
    
    //Calculate the rotation matrix for current segment
    multiplier = multiplier*Get_transformation_matrix(cnt_theta, cnt_phi);
    
    //hout << "i="<<i<<"<"<<step_num<<endl;
    cnt_poi = cnt_poi + Get_new_point(multiplier, nanotube_geo.step_length);
    //1 means that point is not the intial point
    cnt_poi.flag = 1;
    
    return 1;
}
//---------------------------------------------------------------------------
//Calculate the effective portion (length) which falls into the given region defined by a cuboid
double Generate_Network::Length_inside_sample(const cuboid &sample, const Point_3D &prev_point, const Point_3D &new_point, const bool &is_prev_inside_sample, bool &is_new_inside_sample)const
{
    //Check if the new point is inside the sample
    is_new_inside_sample = Is_point_inside_cuboid(sample, new_point);
    
    //Check how to calculate the length inside the sample depending on the location (either
    //inside or outside) of the two points
    if (is_prev_inside_sample) {
        if (is_new_inside_sample) {
            
            //Both points are inside the sample, so return the length between the two points
            return prev_point.distance_to(new_point);
        }
        else {
            
            //The previous point is inside the sample, while the new point is outside the sample
            //Find the intersecting point at the boundary
            Point_3D boundary = Find_intersection_at_boundary(sample, new_point, prev_point);
            
            //Return the length from the boundary point towards the inside point
            return prev_point.distance_to(boundary);
        }
    }
    else {
        if (is_new_inside_sample) {
            
            //The previous point is outside the sample, while the new point is inside the sample
            //Find the intersecting point at the boundary
            Point_3D boundary = Find_intersection_at_boundary(sample, prev_point, new_point);
            
            //Return the length from the boundary point towards the inside point
            return new_point.distance_to(boundary);
            
        }
        else {

            //Both points are outside, thus there is zero length inside the sample
            return 0.0;
        }
    }
}
//This function adds a point to the overlapping subregions it belongs to using a map
int Generate_Network::Add_cnt_point_to_overlapping_regions_map(const Geom_sample &geom_sample, const Point_3D &new_point, const int &local_num, const int &is_new_inside_sample, const int n_subregions[], map<int, vector<int> > &subr_point_map)const
{
    //A point is added only if it is in the composite domain, i.e., inside the sample
    //If the point is in the boundary layer, overlapping is not important
    if (is_new_inside_sample) {
        
        //Create array for loop over overlaping regions
        int flags[2][3];
        
        //Initialize flags for overlaping regions
        int fx = 0;
        int fy = 0;
        int fz = 0;
        
        //Fill the array of flags
        Get_overlapping_flags(geom_sample, new_point, n_subregions, flags, fx, fy, fz);
        
        //In this loop I check all regions a point can belong to when it is in an overlaping zone
        for (int ii = 0; ii < 2; ii++) {
            if (!fx) ii++; //if flag is zero, do this loop only once
            for (int jj = 0; jj < 2; jj++) {
                if (!fy) jj++; //if flag is zero, do this loop only once
                for (int kk = 0; kk < 2; kk++) {
                    if (!fz) kk++; //if flag is zero, do this loop only once
                    
                    //hout <<"a="<<a<<" fx="<<fx<<" b="<<b<<" fy="<<fy<<" c="<<c<<" fz="<<fz;
                    //Calculate the region number
                    int t = Calculate_t(flags[ii][0], flags[jj][1], flags[kk][2], n_subregions[0], n_subregions[1]);
                    
                    //Add the map from subregion to local point number: t->local_num
                    subr_point_map[t].push_back(local_num);
                }
            }
        }
    }
    
    return 1;
}
//that is used to deal with self penetration
//This function checks if the newly generated CNT needs to be stored or ignore
int Generate_Network::Store_or_ignore_new_cnt_using_map(const int &penetration_model_flag, const int &points_in, const double &cnt_len, const double &cnt_rad, const double &cnt_cross_area, const vector<Point_3D> &new_cnt, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radius, const map<int, vector<int> > &subr_point_map, vector<vector<long int> > &sectioned_domain, vector<vector<int> > &global_coordinates, double &vol_sum, int &cnt_ignore_count)const
{
    //Store the CNT points
    //hout << "Store CNT ";
    //Check that is a valid CNT
    //If size is 0, the CNT was rejected
    //If size is 1, then only the CNT seed was generated in a valid position but the second point
    //could not be placed in a valid position (i.e., without interpenetrating another CNT)
    if(new_cnt.size() >= 2 && points_in) {
        
        //Update the volume and weigth of generated CNTs
        vol_sum += cnt_len*cnt_cross_area;
        
        //If the new_cnt vector has at least two points, then it can be added to the rest of the points
        cnts_points.push_back(new_cnt);
        
        //Add the radius to the vector of radii
        cnts_radius.push_back(cnt_rad);
        
        //Perform these operations when the non-overlapping model is used
        if (penetration_model_flag) {
            
            //Get the size of the global cooordinates vector before add the points in the new CNT
            long int gc_size = global_coordinates.size();
            
            //Iterate over all elements in the map
            for (map<int, vector<int> >::const_iterator it = subr_point_map.begin(); it != subr_point_map.end(); ++it) {
                
                //Get the current subregion number
                int subregion = it->first;
                
                //Iterate over the points in the vector of the current map
                for (size_t i = 0; i < it->second.size(); i++) {
                    
                    //Get the local point number of the current CNT
                    long int Pi = it->second[i];
                    
                    //Add point to its overlapping subregion
                    //The global number is the size of the sectioned domain saved in gc_size
                    //plus the local number Pi
                    sectioned_domain[subregion].push_back(gc_size+Pi);
                }
            }
            
            //Variable needed for updating global_coordinates
            vector<int> empty(2);
            
            //Get the CNT number
            int CNT = (int)cnts_points.size()-1;
            
            //Scan all points in the new CNT to add global coordinates
            for(int ii = 0; ii < (int)new_cnt.size(); ii++) {
                
                //Add global coordinate
                empty[0] = CNT;
                empty[1] = ii;
                global_coordinates.push_back(empty);
            }
        }
    }
    else if (!points_in){
        
        //The CNT can be ignored as it is completely outside the sample
        cnt_ignore_count++;
    }
    //hout << "done"<<endl;
    
    return 1;
}
int Generate_Network::Store_or_ignore_new_cnt(const Geom_sample &geom_sample, const int &penetration_model_flag, const int &points_in, const double &cnt_len, const double &cnt_rad, const double &cnt_cross_area, const vector<Point_3D> &new_cnt, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radius, const int n_subregions[], vector<vector<long int> > &sectioned_domain, vector<vector<int> > &global_coordinates, double &vol_sum, int &cnt_ignore_count)const
{
    //Store the CNT points
    //hout << "Store CNT ";
    //Check that is a valid CNT
    //If size is 0, the CNT was rejected
    //If size is 1, then only the CNT seed was generated in a valid position but the second point
    //could not be placed in a valid position (i.e., without interpenetrating another CNT)
    if(new_cnt.size() >= 2 && points_in) {
        //Update the volume and weigth of generated CNTs
        vol_sum += cnt_len*cnt_cross_area;
        
        //If the new_cnt vector has at least two points, then it can be added to the rest of the points
        cnts_points.push_back(new_cnt);
        
        //Add the radius to the vector of radii
        cnts_radius.push_back(cnt_rad);
        
        //Perform these operations when the non-overlapping model is used
        if (penetration_model_flag) {
            
            //Variable needed for updating global_coordinates
            vector<int> empty(2);
            
            //Get the CNT number
            int CNT = (int)cnts_points.size()-1;
            
            //Scan all points in the CNT that was just generated
            for(int ii = 0; ii < (int)new_cnt.size(); ii++) {
                
                //Add global coordinate
                empty[0] = CNT;
                empty[1] = ii;
                global_coordinates.push_back(empty);
                //Add point to an overlapping region in the vector sectioned_domain
                Add_cnt_point_to_overlapping_regions(geom_sample, new_cnt[ii], (long int)global_coordinates.size()-1, n_subregions, sectioned_domain);
                
            }
        }
    }
    else if (!points_in){
        
        //The CNT can be ignored as it is completely outside the sample
        cnt_ignore_count++;
    }
    //hout << "done"<<endl;
    
    return 1;
}
//---------------------------------------------------------------------------
//This function adds a point to a region so penetration can be checked
void Generate_Network::Add_cnt_point_to_overlapping_regions(const Geom_sample &geom_sample, const Point_3D &point, const long int &global_num, const int n_subregions[], vector<vector<long int> > &sectioned_domain)const
{
    //A point is added only if it is in the composite domain
    //If the point is in the boundary layer, overlapping is not important
    if (Is_point_inside_cuboid(geom_sample.sample, point)) {
        
        //Create array for loop over overlaping regions
        int flags[2][3];
        
        //Initialize flags for overlaping regions
        int fx = 0;
        int fy = 0;
        int fz = 0;
        
        //Fill the array of flags
        Get_overlapping_flags(geom_sample, point, n_subregions, flags, fx, fy, fz);
        
        //In this loop I check all regions a point can belong to when it is in an overlaping zone
        for (int ii = 0; ii < 2; ii++) {
            if (!fx) ii++; //if flag is zero, do this loop only once
            for (int jj = 0; jj < 2; jj++) {
                if (!fy) jj++; //if flag is zero, do this loop only once
                for (int kk = 0; kk < 2; kk++) {
                    if (!fz) kk++; //if flag is zero, do this loop only once
                    //hout <<"a="<<a<<" fx="<<fx<<" b="<<b<<" fy="<<fy<<" c="<<c<<" fz="<<fz;
                    int t = Calculate_t(flags[ii][0], flags[jj][1], flags[kk][2], n_subregions[0], n_subregions[1]);
                    //hout<<" t="<<t<<" sectioned_domain["<<t<<"].size()="<<sectioned_domain[t].size();
                    sectioned_domain[t].push_back(global_num);
                    //hout<<'.'<<endl;
                }
            }
        }
    }
}
//---------------------------------------------------------------------------
//This function determines the individual flags and fills the array of flags that
//help determine if a point is in more than one overlapping region
void Generate_Network::Get_overlapping_flags(const Geom_sample &geom_sample, const Point_3D &point, const int n_subregions[], int flags[][3], int &fx, int &fy, int &fz)const
{
    //Save coordinates of the point
    double x = point.x;
    double y = point.y;
    double z = point.z;
    
    //These variables will give me the region cordinates of the region that a point belongs to
    int a, b, c;
    
    //Calculate the subregion-coordinates
    Get_subregion_coordinates(geom_sample, n_subregions, point, a, b, c);
    
    //These variables are the coordinates of the lower corner of the RVE that defines its geometry
    double xmin = geom_sample.sample.poi_min.x;
    double ymin = geom_sample.sample.poi_min.y;
    double zmin = geom_sample.sample.poi_min.z;
    
    //Coordinates of non-overlaping region the point belongs to
    double x1 = (double)a*geom_sample.gs_minx +  xmin;
    double x2 = x1 + geom_sample.gs_minx;
    double y1 = (double)b*geom_sample.gs_miny +  ymin;
    double y2 = y1 + geom_sample.gs_miny;
    double z1 = (double)c*geom_sample.gs_minz +  zmin;
    double z2 = z1 + geom_sample.gs_minz;
    
    //Assign value of flag according to position of point
    //The first operand eliminates the periodicity on the boundary
    if ((a > 0) && (x >= x1) && (x <= x1+geom_sample.gs_overlap_cnt))
        fx = -1;
    else if ((a < n_subregions[0]-1) && (x >= x2-geom_sample.gs_overlap_cnt) && (x <= x2 ))
        fx = 1;
    if ((b > 0) && (y >= y1) && (y <= y1+geom_sample.gs_overlap_cnt))
        fy = -1;
    else if ((b < n_subregions[1]-1) && (y >= y2-geom_sample.gs_overlap_cnt) && (y <= y2 ))
        fy = 1;
    if ((c > 0) && (z >= z1) && (z <= z1+geom_sample.gs_overlap_cnt))
        fz = -1;
    else if ((c < n_subregions[2]-1) && (z >= z2-geom_sample.gs_overlap_cnt) && (z <= z2 ))
        fz = 1;
    
    //Fill the array of flags for overlapping regions
    //int flags[2][3] = { {a+fx, b+fy, c+fz}, {a, b, c}};
    flags[0][0] = a+fx;
    flags[0][1] = b+fy;
    flags[0][2] = c+fz;
    flags[1][0] = a;
    flags[1][1] = b;
    flags[1][2] = c;
}
//---------------------------------------------------------------------------
//Calculates the region to which a point corresponds
int Generate_Network::Calculate_t(const int &a, const int &b, const int &c, const int &sx, const int &sy)const
{
    return a + b*sx + c*sx*sy;
}
//---------------------------------------------------------------------------
//Transform the 2D cnts_points into 1D cpoints and 2D cstructures
int Generate_Network::Transform_points_cnts(const Geom_sample &geom_sample, const Nanotube_Geo &nano_geo, vector<vector<Point_3D> > &cnts_points, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures)const
{
    //Variable to count the point numbers
    long int point_count = 0;
    
    //Variable to count the CNT numbers
    int cnt_count = 0;
    
    //Variable to count the generated points
    long int n_points = 0;
    
    //Choose the way in which the output vectors are generated depending on the particle type
    for(int i=0; i<(int)cnts_points.size(); i++)
    {
        
        if (!Add_cnts_inside_sample(geom_sample, nano_geo, i, cnts_points[i], cpoints, radii_in, radii_out, cstructures, point_count, cnt_count)) {
            hout<<"Error when adding CNTs to structure."<<endl;
            return 0;
        }
        
        //Update the number of generated points
        n_points += (long int)cnts_points[i].size();
         
        //Free some memory
        cnts_points[i].clear();
    }
    
    if (cnts_points.size()) {
        //hout<<"Radii="<<radii_in.size()<<endl;
        hout << "There are " << cnts_points.size() << " CNTs with "<<n_points<<" points at least partially inside the sample before trimming."<<endl;
        hout<<"There are "<<cnt_count<<" CNTs with "<<cpoints.size() << " points inside the sample."<<endl;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
int Generate_Network::Add_cnts_inside_sample(const Geom_sample &geom_sample, const Nanotube_Geo &nano_geo, const int &CNT_old, vector<Point_3D> &cnt, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures, long int &point_count, int &cnt_count)const
{
    //Indices to define the beginning and ending of a segment of a CNT that is inside a sample
    int start = 0;
    int end = 0;
    
    //Index of the last point inside the sample
    int last_inside = 0;
    
    //Number of points in the current CNT
    int cnt_points = (int)cnt.size();
    
    //Provisionally the minimum number of points to consider a CNT is defined here
    int min_points = nano_geo.min_points;
    
    //Check where is the first point of the CNT
    bool is_first_inside_sample = Is_point_inside_cuboid(geom_sample.sample, cnt[0]);
    
    //Scan all remaining points in the current CNT
    for (int i = 1; i < cnt_points; i++) {
        
        //Check if the point is inside the sample
        if (Is_point_inside_cuboid(geom_sample.sample, cnt[i])) {
            
            //Update the last inside point
            last_inside = i;
        }
        else {
            
            //End index is the current looping index
            end = i;
            
            //Count the number of consecutive points inside the sample
            //If start is zero and it is inside the sample, then the number of points is end - start
            //because end is always outside
            //If start is not zero, then it is always outside, and since end is also outside, then
            //the number of points is end - start - 1
            int n_points = (start == 0 && is_first_inside_sample)? end - start: end - start - 1;
            
            //Check if there are enough points to consider this a whole CNT and include it in the analysis
            if (n_points >= min_points) {
                
                //Check if there are are enough points and, if so, add the current CNT segment to the data structures
                if (!Add_cnt_segment(geom_sample, is_first_inside_sample, start, end, CNT_old, cnt, cpoints, radii_in, radii_out, cstructures, point_count, cnt_count)) {
                    hout<<"Error when adding a CNT segment."<<endl;
                    return 0;
                }
            }
            //else{hout<<"Segment with "<<n_points<<" points, start="<<start<<" end="<<end<<" cnt_count="<<cnt_count<<endl;}
            
            //Reset the start index
            start = i;
        }
    }
    
    //Check if the last point of the CNT was inside the sample
    //Also check that the segment has the minimum number of CNT points required
    //for it to be consedered a CNT
    if (last_inside == cnt_points-1 && (last_inside - start + 1) >= min_points) {
        
        //Set end index as one after the last valid index
        //In this way, when adding CNTs to the structure, the last CNT point is added
        end = cnt_points;
        
        //If the last point of the CNT was inside the sample, then add a new segment
        //This was not done becuase, in the for loop, a segement is added only when it finds a point
        //outside the sample
        //Then, check if there are are enough points and, if so, add the current CNT segment to the data structures
        if (!Add_cnt_segment(geom_sample, is_first_inside_sample, start, end, CNT_old, cnt, cpoints, radii_in, radii_out, cstructures, point_count, cnt_count)) {
            hout<<"Error when adding a CNT segment."<<endl;
            return 0;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
int Generate_Network::Add_cnt_segment(const Geom_sample &geom_sample, const bool &is_first_inside_sample, const int &start, const int &end, const int &CNT_old, vector<Point_3D> &cnt, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures, long int &point_count, int &cnt_count)const
{
    
    //Temporary vector to add the point numbers to the structure vector
    vector<long int> struct_temp;
    
    //Variable to recalculate start index if needed
    int new_start = start;
    
    //Check if start is inside the boundary, this is needed to determine percolation
    //This happens when:
    //1) Always when start is non-zero. Whenever start index is not zero, the
    //   CNT segment starts outside the sample and crosses the sample boundary.
    //2) The only case when start is zero but is outside the sample is when the first
    //   CNT point is outside the sample and the second CNT point is inside the sample.
    //   For this case use the flag is_first_inside_sample
    if (start != 0 || (start == 0 && !is_first_inside_sample)) {
        
        //Note that start index is a point outside the sample
        //hout<<"start="<<start<<endl;
        //hout<<"P_start = ("<<cnt[start].x<<", "<<cnt[start].y<<", "<<cnt[start].z<<") cnt_count="<<cnt_count<<endl;
        if (!Add_boundary_point(geom_sample, cnt[start], cnt[start+1], cnt_count, cpoints, struct_temp, point_count)) {
            hout<<"Error in Add_boundary_point when adding a point at the start of the segment."<<endl;
        }
        
        //Start is outside, so the first inside point is start+1
        //Thus update start
        new_start++;
    }
    
    //Add the CNT points of the segment found to the 1D vector
    //Note that end index is actually one more than the last index of the segment
    for(int j = new_start; j < end; j++) {
        
        //Change the flag of current point to be that of its CNT number
        cnt[j].flag = cnt_count;
        
        //Add current point
        cpoints.push_back(cnt[j]);
        
        //Add the point number to the structure vector
        struct_temp.push_back(point_count);
        
        //Increase the count of points
        point_count++;
    }
    
    //Check if end index is not the last point of the CNT (since end index is actually one more
    //than the last index of the segment, then we have to check against the largest posisble index+1
    //which is the number of points in the CNT)
    //If the last index is the last point of the CNT, then the CNT segement ends inside the sample.
    //If the last index is not the last point of the CNT, then the CNT segment crosses the boundary
    //In such case, a point at the boundary is added. This is needed to determine percolation
    if (end != (int)cnt.size()) {
        
        //hout<<"end="<<end<<" points="<<cnt.size()<<endl;
        //hout<<"P_end = ("<<cnt[end].x<<", "<<cnt[end].y<<", "<<cnt[end].z<<") cnt_count="<<cnt_count<<endl;
        if (!Add_boundary_point(geom_sample, cnt[end], cnt[end-1], cnt_count, cpoints, struct_temp, point_count)) {
            hout<<"Error in Add_boundary_point when adding a point at the end of the segment."<<endl;
        }
    }
    
    //Add the corresponding radius to the output darii vector
    radii_out.push_back(radii_in[CNT_old]);
    
    //Add the temporary structure vector
    cstructures.push_back(struct_temp);
    
    //A CNT segment was added, i.e., a new CNT was added, thus increase the CNT count
    cnt_count++;
    
    return 1;
}
//---------------------------------------------------------------------------
int Generate_Network::Add_boundary_point(const Geom_sample &geom_sample, const Point_3D &p_outside, const Point_3D &p_inside, const int &cnt_count, vector<Point_3D> &cpoints, vector<long int> &struct_temp, long int &point_count)const
{
    //Find the coordinates of the point between the ouside (p_outside) and inside (p_inside) that
    //that is located at the sample boundary (one of the faces)
    Point_3D boundary = Find_intersection_at_boundary(geom_sample.sample, p_outside, p_inside);
    
    //Update the points flag
    boundary.flag = cnt_count;
    
    //Add the point to the structure vectors
    cpoints.push_back(boundary);
    
    //Add the point number to the structure vector
    struct_temp.push_back(point_count);
    
    //Increase the count of points
    point_count++;
    
    return 1;
}
//---------------------------------------------------------------------------
//Depending on the location of the outside point, the line defined by consecutive inside and outside points
//might have up to three intersections with the planes that define the boundaries of the sample.
//Of course, a point might be on the plane of a boundary, but ouside of that boundary.
//Thus, if there are multiple intersections with the planes, we need to check whether or not
//the interstion is at an actual boundary
//This function finds that intersecting point at an actual boundary
Point_3D Generate_Network::Find_intersection_at_boundary(const cuboid &sample, const Point_3D &p_outside, const Point_3D &p_inside)const
{
    //The line segment defined by p_outside and p_inside is given by:
    //P = p_outside + lambda*T
    //where T = p_inside - p_outside
    //In this way P = p_outside when lambda = 0 and P = p_inside when lambda = 1
    //Variable to store the point T = p_inside - p_outside
    Point_3D T = p_inside - p_outside;
    
    //Variable to store the point at the intersection of the line segment (between p_outside and p_inside)
    //and the boundary
    Point_3D boundary;
    
    //Lambda function to calculate the lambda coefficient, since I only use it three times here and is a
    //simple calculation, so I rather use a lambda function instead of declaring a new proper function
    //as part of the class
    auto calc_lambda = [](auto x_plane, auto x_out, auto x_T) {return (x_plane - x_out)/x_T;};
    
    //Go through each boundary and find the lambda value for those that are intersected for the
    //line segment that goes from p_outside to p_inside
    
    //Check if any of the x-boundaries is intersected
    double lambda_x = - 1.0;
    //x-left boundary
    if ( (p_outside.x - sample.poi_min.x) < Zero ) {
        
        //Calculate the lambda value
        lambda_x = calc_lambda(sample.poi_min.x, p_outside.x, T.x);
    }
    //x-right boundary
    else if ( (sample.max_x - p_outside.x) < Zero ) {
        
        //Calculate the lambda value
        lambda_x = calc_lambda(sample.max_x, p_outside.x, T.x);
    }
    //hout<<"lambda_x="<<lambda_x<<endl;
    
    //Variable to save a new value of lambda, if needed
    //If a new lambda turns out to be larger, then the old lambda needs to be updated
    //Since we need to find a value larger than lambda, which is in [0,1], then
    //new_lambda is initialized with a value smaller than lambda.
    //A negative value ensures this new_lambda will be smaller than any value lambda could get
    
    //Check if any of the y-boundaries is intersected
    double lambda_y = -1.0;
    //y-left boundary
    if ( (p_outside.y - sample.poi_min.y) < Zero ) {
        
        //Calculate the lambda value
        lambda_y = calc_lambda(sample.poi_min.y, p_outside.y, T.y);
    }
    //y-right boundary
    else if ( (sample.max_y - p_outside.y) < Zero ) {
        
        //Calculate the lambda value
        lambda_y = calc_lambda(sample.max_y, p_outside.y, T.y);
    }
    //hout<<"lambda_y="<<lambda_y<<endl;
    
    //Check if any of the z-boundaries is intersected
    double lambda_z = -1.0;
    //z-left boundary
    if ( (p_outside.z - sample.poi_min.z) < Zero ) {
        
        //Calculate the lambda value
        lambda_z = calc_lambda(sample.poi_min.z, p_outside.z, T.z);
    }
    //z-right boundary
    else if ( (sample.max_z - p_outside.z) < Zero ) {
        
        //Calculate the lambda value
        lambda_z = calc_lambda(sample.max_z, p_outside.z, T.z);
    }
    //hout<<"lambda_z="<<lambda_z<<endl;
    
    //Variable to store the coefficient lambda to parameterize the line segment between
    //p_outside (lambda = 0) and p_inside (lambda = 1)
    //The lambda needed is the largest one
    double lambda = max(lambda_x, max(lambda_y, lambda_z));
    
    //Calculate the new point using the largest lambda
    boundary = p_outside + T*lambda;
    
    //Check if a coordinate is close to zero, in which case set it to zero directly
    //to reduce numerical errors
    if (abs(boundary.x) < Zero) {
        boundary.x = 0.0;
    }
    if (abs(boundary.y) < Zero) {
        boundary.y = 0.0;
    }
    if (abs(boundary.z) < Zero) {
        boundary.z = 0.0;
    }
    
    //hout<<"P_outside = ("<<p_outside.x<<", "<<p_outside.y<<", "<<p_outside.z<<")"<<endl;
    //hout<<"P_inside = ("<<p_inside.x<<", "<<p_inside.y<<", "<<p_inside.z<<")"<<endl;
    //hout<<"P_T = ("<<T.x<<", "<<T.y<<", "<<T.z<<")"<<endl;
    //hout<<"P_boundary = ("<<boundary.x<<", "<<boundary.y<<", "<<boundary.z<<")"<<endl<<endl;
    
    return boundary;
}
//---------------------------------------------------------------------------
int Generate_Network::Recalculate_vol_fraction_cnts(const Geom_sample &geom_sample, const Simu_para &simu_para, const Nanotube_Geo &nano_geo, const vector<Point_3D> &cpoints, const vector<double> &radii, const vector<vector<long int> > &cstructures)const
{
    //Variable to store the volume of CNTs
    double cnt_vol = 0.0;
    
    //Iterate over all CNT in the structure vector
    for (int i = 0; i < (int)cstructures.size(); i++) {
        
        //Initialize the length of the CNT
        double cnt_len = 0.0;
        
        /*if (!Is_point_inside_cuboid(geom_sample.sample, cpoints[cstructures[i][0]])) {
            hout<<"Point outside the sample: cpoints["<<cstructures[i][0]<<"(i="<<i<<", j="<<0<<")]="<<cpoints[cstructures[i][0]].str()<<" CNT with "<<cstructures[i].size()<<" points"<<endl;
        }*/
        
        //Calculate the total length of current CNT_i
        for (int j = 1; j < (int)cstructures[i].size(); j++) {
            
            //Get the previous and current points
            long int prev = cstructures[i][j-1];
            long int curr = cstructures[i][j];
            
            //Add the length of two consecutive points
            cnt_len = cnt_len + cpoints[prev].distance_to(cpoints[curr]);
            
            /*if (!Is_point_inside_cuboid(geom_sample.sample, cpoints[curr])) {
                hout<<"A point is outside the sample: cpoints["<<curr<<"(i="<<i<<", j="<<j<<")]="<<cpoints[curr].str()<<endl;
            }*/
        }
        
        //Calculate the volume from the CNT length and add it to the total volume
        cnt_vol = cnt_vol + PI*radii[i]*radii[i]*cnt_len;
    }
    
    if (simu_para.criterion == "vol") {
        hout<<"CNT volume fraction after removing CNT outside the sample and small segments is: "<<cnt_vol/geom_sample.volume<<endl;
    }
    else {
        double cnt_wt = cnt_vol*nano_geo.density;
        double matrix_wt = (geom_sample.volume - cnt_vol)*geom_sample.matrix_density;
        hout<<"CNT weight fraction after removing CNT outside the sample and small segments is: "<<cnt_wt/(cnt_wt+matrix_wt)<<endl;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//CNT deposit
int Generate_Network::Generate_cnt_deposit_mt(const Simu_para &simu_para, const Geom_sample &geom_sample, const Geom_sample &geom_sample_deposit, const Nanotube_Geo &nanotube_geo, const Cutoff_dist &cutoffs, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radii)const
{
    //Initial seeds, if any are in network_seeds within geom_sample.
    //However, geom_sample cannot be modified, so copy the seeds to a new vector
    unsigned int net_seeds[7];
    if (!CNT_seeds(simu_para.CNT_seeds, net_seeds)) {
        hout<<"Error in CNT_seeds"<<endl;
        return 0;
    }
    
    //Use the seeds generated above
    std::mt19937 engine_x(net_seeds[0]);
    std::mt19937 engine_y(net_seeds[1]);
    std::mt19937 engine_z(net_seeds[2]);
    std::mt19937 engine_phi(net_seeds[3]);
    std::mt19937 engine_theta(net_seeds[4]);
    std::mt19937 engine_rand(net_seeds[5]);
    std::mt19937 engine_initial_direction(net_seeds[6]);
    
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [-1, 1].
    std::uniform_real_distribution<double> dist_initial(-1.0, 1.0);

    //Time variables to keep track of generation time
    time_t ct0, ct1;
    
    //Variables for the generated CNT volume and weight
    double vol_sum = 0;
    
    //Variable to count the generated CNT seeds
    int cnt_seed_count = 0;
    
    //Variable to coun the number of rejected CNTs
    int cnt_reject_count = 0;
    
    //Vectors for handling CNT penetration
    //global_coordinates[i][0] stores the CNT number of global point i
    //global_coordinates[i][1] stores the local point number of global point i
    vector<vector<int> > global_coordinates;
    //sectioned_domain[i] contains all the points in sub-region i.
    //Sub-region i is an overlapping subregion to check for penetrations
    vector<vector<long int> > sectioned_domain;
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    //n_subregions[2] is the number of subregions along z
    int n_subregions[3];
    //Initialize the vector sub-regions
    if (!Initialize_cnt_subregions_extended_domain(geom_sample, n_subregions, sectioned_domain)) {
        hout << "Error in Generate_network_threads_mt when calling Initialize_cnt_subregions_extended_domain" <<endl;
        return 0;
    }
    
    //Get the time when generation started
    ct0 = time(NULL);
    
    //Check when 10% is completed
    double vol_completed = 0.1;
    //Variable used to store the fraction of CNT volume indicated by vol_completed
    double vol_completed_acc = vol_completed*nanotube_geo.volume;
    
    //Boolean to terminate the main while loop, initialized to false to start the loop
    bool terminate = false;
    
    //Generate CNTs while the total volume required is not reached
    while(!terminate)
    {
        //Vector for a new nanotube
        vector<Point_3D> new_cnt;
        
        //Randomly generate a CNT length and CNT radius
        double cnt_length, cnt_rad;
        //hout<<"Get_random_value_mt 1"<<endl;
        if (!Get_length_and_radius(nanotube_geo, engine_rand, dist, cnt_length, cnt_rad)) {
            hout << "Error in Generate_network_threads_mt when calling Get_length_and_radius" <<endl;
            return 0;
        }
        
        //Calculate the number of CNT growth steps
        int step_num = (int)(cnt_length/nanotube_geo.step_length) + 1;
        
        //Cross-sectional area of the current CNT. It is used to calculate the CNT volume.
        const double cnt_cross_area = PI*cnt_rad*cnt_rad;
        
        //Randomly generate a CNT seed (initial point) at the bottom of the extended domain
        Point_3D new_point;
        if (!Get_point_in_xy_plane_mt(geom_sample.ex_dom_cnt, new_point, engine_x, engine_y, dist)) {
            hout<<"Error when generating a point in the xy plane"<<endl;
            return 0;
        }
        
        //Get the sub-region the point belongs to
        //hout<<"Get_subregion"<<endl;
        int subregion = Get_subregion_for_cnt_seed_in_deposit(geom_sample, sectioned_domain, n_subregions, new_point);
        
        //Find the upmost position where the CNT seed can be placed
        if (!Find_upmost_position_for_seed(geom_sample, cnts_points, cnts_radii, global_coordinates, sectioned_domain, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, subregion, new_point)) {
            hout<<"Error in Generate_network_threads_mt when calling Find_upmost_position_for_seed"<<endl;
            return 0;
        }
        
        //Add the CNT seed to the current CNT vector
        new_cnt.push_back(new_point);
        
        //Count the numer of CNT seeds generated
        cnt_seed_count++;
        //hout << "Seed="<<cnt_seed_count<<endl;
        int max_seed = 1E9;
        if(cnt_seed_count > max_seed)  {
            hout<<"The number of generated seeds is lager than "<<max_seed<<", but the nanotube generation still fails to acheive the requested volume fraction."<<endl;
            return 0;
        }
        
        //Randomly generate an initial direction, then generate the rotation matrix
        //that results in that rotation
        MathMatrix mult_2d(2,2);
        if (!Get_initial_direction_2d(mult_2d, engine_theta, dist)) {
            hout<<"Error when generating an initial 2D direction."<<endl;
            return 0;
        }
        
        //hout<<"Growth of CNT"<<endl;
        
        //Get the location of the seed
        bool is_prev_in_sample = Is_point_inside_cuboid(geom_sample.sample, new_cnt[0]);
        
        //Variable to store the length of the current CNT that is inside the sample
        double cnt_len = 0.0;
        
        //Start generation of a CNT given the generated seed
        for (int i = 0; i < step_num; i++) {
            
            //Generate a point in a random direction in the plane xy and find its highest position
            //NOTE: due to the nature of the deposit, the non-penetrating model is always used
            int status = Find_highest_position_for_new_point_iteratively(geom_sample_deposit, nanotube_geo, cnts_points, cnts_radii, global_coordinates, sectioned_domain, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, new_cnt, new_point, mult_2d, engine_theta, dist);
            
            //Check the value of status
            if (status == -1) {
                hout<<"Error when finding the highest position of a new point (not a seed)"<<endl;
                return 0;
            }
            else if (!status) {
                
                //A point coud not be generated, so discard current CNT
                
                //Clear the new_cnt vector so that it is not added to the rest of CNTs
                new_cnt.clear();
                
                //Increase the count of rejected cnts
                cnt_reject_count++;
                
                //Break the for-loop to terminate the CNT growth
                break;
            }
            
            //Check if new_point is inside the extended domain
            if(Is_point_inside_cuboid(geom_sample.ex_dom_cnt, new_point)) {
                
                //Calculate the segment length inside the sample and add it to the total CNT length
                bool is_new_inside_sample;
                cnt_len = cnt_len + Length_inside_sample(geom_sample.sample, new_cnt.back(), new_point, is_prev_in_sample, is_new_inside_sample);
                
                //For the next iteration of the for loop, cnt_poi will become previous point,
                //so update the boolean is_prev_in_sample for the next iteration
                is_prev_in_sample = is_new_inside_sample;
                
                //Add the new point to the current CNT
                new_cnt.push_back(new_point);
                
                //Check if the target volume fraction has been reached
                if( (vol_sum + cnt_len*cnt_cross_area) >= nanotube_geo.volume) {
                    
                    //Set the terminate variable to true so that the main while-loop is terminated
                    terminate = true;
                    
                    //Break the for-loop so that the current CNT stops its growth
                    break;
                }
                
            }
            else {
                
                //If the point is outside the extended domain, break the for-loop
                //to terminate the CNT growth
                break;
            }
            
            //for-loop ends here
        }
        
        //Store or ignore the CNT points
        //hout<<"Store_or_ignore_new_cnt"<<endl;
        int cnt_ignore_count = 0;
        //To add points to global coordinates, need to use the sample geometry of the deposit
        if (!Store_or_ignore_new_cnt(geom_sample_deposit, simu_para.penetration_model_flag, (int)new_cnt.size(), cnt_len, cnt_rad, cnt_cross_area, new_cnt, cnts_points, cnts_radii, n_subregions, sectioned_domain, global_coordinates, vol_sum, cnt_ignore_count)) {
            hout<<"Error when storing or ignoring a new CNT"<<endl;
            return 0;
        }
        
        //Get the time to check progress
        ct1 = time(NULL);
        
        //Check progress
        if (!Check_progress("CNT", (int)(ct1-ct0), nanotube_geo.volume, vol_sum, vol_completed, vol_completed_acc)) {
            hout<<"Error when calculating the percentage of CNT volume generated"<<endl;
            return 0;
        }
        
    }
    
    //Output the CNT content generated
    if(nanotube_geo.criterion == "wt") {
        
        //Calculate matrix weight
        double matrix_weight = (geom_sample.volume - vol_sum)*geom_sample.matrix_density;
        
        //Calculate the CNT weight
        double cnt_weight = vol_sum*nanotube_geo.density;
        
        hout << endl << "The weight fraction of generated CNTs is: " << cnt_weight/(matrix_weight + cnt_weight);
        hout << ", the target weight fraction was " << nanotube_geo.weight_fraction << endl << endl;

    } else if(nanotube_geo.criterion == "vol") {
        hout << endl << "The volume fraction of generated CNTs was " << vol_sum/geom_sample.volume;
        hout << ", the target volume fraction was " << nanotube_geo.volume_fraction << endl << endl;
    }
    
    hout << "There were " << cnt_reject_count << " CNTs rejected." << endl;
    
    return 1;
}
//This functions initializes the vectors n_subregions and sectioned_domain
//
//The n_subregions vector is defined to avoid calculating the number of sub-regions for every point
//n_subregions[0] is the number of subregions along x
//n_subregions[1] is the number of subregions along y
//n_subregions[2] is the number of subregions along z
//
//The vector sectioned_domain contains the sub-regions to look for overlapping
//It is initialized with the number of sub-regions in the sample
int Generate_Network::Initialize_cnt_subregions_extended_domain(const Geom_sample &sample_geom, int n_subregion[], vector<vector<long int> > &sectioned_domain)const
{
    //Calculate the number of variables along each direction
    //Make sure there is at least one subregion along each direction
    //
    //Number of subregions along x
    n_subregion[0] = max(1, (int)(sample_geom.ex_dom_cnt.len_x/sample_geom.gs_minx));
    //Number of subregions along y
    n_subregion[1] = max(1, (int)(sample_geom.ex_dom_cnt.wid_y/sample_geom.gs_miny));
    //Number of subregions along z
    n_subregion[2] = max(1, (int)(sample_geom.ex_dom_cnt.hei_z/sample_geom.gs_minz));
    
    //Initialize sectioned_domain
    sectioned_domain.assign(n_subregion[0]*n_subregion[1]*n_subregion[2], vector<long int>());
    
    return 1;
}
//Randomly generate a CNT seed (intial point of a CNT) in the xy plane of the extended domain
int Generate_Network::Get_point_in_xy_plane_mt(const cuboid &cub, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, uniform_real_distribution<double> &dist)const
{
    point.x = cub.poi_min.x + cub.len_x*dist(engine_x);
    
    point.y = cub.poi_min.y + cub.wid_y*dist(engine_y);
    
    point.z = cub.poi_min.z;
    
    //Zero flag denotes this point is the initial point of a CNT
    point.flag = 0;
    
    return 1;
}
//This function returns the subregion a point belongs to
//Here, the boundary layer is also included in the subregions
int Generate_Network::Get_subregion_for_cnt_seed_in_deposit(const Geom_sample &geom_sample, const vector<vector<long int> > &sectioned_domain, const int n_subregions[], const Point_3D &point)const
{
    //For a CNT deposit, only x and y coordinates are randomly assigned
    //Thus, first I need the a- and b- coordinates of the subregions
    
    //Calculate the region-coordinates a and b
    int a = (int)((point.x-geom_sample.ex_dom_cnt.poi_min.x)/geom_sample.gs_minx);
    //Limit the value of a as it has to go from 0 to n_subregions[0]-1
    if (a == n_subregions[0]) a--;
    int b = (int)((point.y-geom_sample.ex_dom_cnt.poi_min.y)/geom_sample.gs_miny);
    //Limit the value of b as it has to go from 0 to n_subregions[1]-1
    if (b == n_subregions[1]) b--;
    
    //Now, I need to find the c-coordinate of the subregions
    //To do this, I need to start at the highest c in the (a,b) column
    //From there I go down and find the first non-empty subregion
    //The case i=0 can be ignored in the loop, as if that coordinate is reached
    //it does not matter if there are points or not, that is the c-coordinate of the subregion
    for (int i = n_subregions[2]-1; i >= 1; i--) {
        
        //Check if subregion wihth c-coordinate equal to i is empty
        int subregion = a + (b*n_subregions[0]) + (i*n_subregions[0]*n_subregions[1]);
        if (sectioned_domain[subregion].size()) {
            
            //There are elements in subregion i
            //Thus i is the c-coordinate, so calculate the subregion using i as the c-coordinate
            return (a + (b*n_subregions[0]) + (i*n_subregions[0]*n_subregions[1]));
        }
    }
    
    //If all subregions were empty, then return the subregion calculated only with a and b
    return (a + (b*n_subregions[0]) );
}
//This function finds the upmost position for a CNT seed
int Generate_Network::Find_upmost_position_for_seed(const Geom_sample &geom_sample, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radii, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const int n_subregions[], const double &cnt_rad, const double &d_vdW, const int &subregion, Point_3D &new_point)const
{
    //Flag to determine if any 2D overlapping was found
    bool overlap_2d = false;
    
    //Check if there is a 2D overlapping
    if (!Check_2d_overlapping(sectioned_domain[subregion], cnts_points, cnts_radii, global_coordinates, cnt_rad+d_vdW, new_point, overlap_2d)) {
        hout<<"Error in Find_upmost_position_for_seed when calling Check_2d_overlapping"<<endl;
        return 0;
    }
    
    //Check if there was no overlapping and decide the next step
    if (!overlap_2d) {
        
        //Check if there is a "layer below" in the sectioned domain
        if (subregion < n_subregions[0]*n_subregions[1]) {
            
            //There is no overlapping point and there are no more subregions below
            //the current one, so set the point on the surface of the sample
            //To do this, set the z-coordinate to be the minimum z-coordinate of the sample
            //plus the radius
            new_point.z = geom_sample.sample.poi_min.z + cnt_rad;
        }
        else {
            
            //Go to the subregion in the layer below by calling this function recursively
            //The layer below is found by substracting n_subregions[0]*n_subregions[1] to
            //the current subregion
            if (!Find_upmost_position_for_seed(geom_sample, cnts_points, cnts_radii, global_coordinates, sectioned_domain, n_subregions, cnt_rad, d_vdW, subregion-n_subregions[0]*n_subregions[1], new_point)) {
                hout<<"Error in Find_upmost_position_for_seed when calling itself recursively"<<endl;
                return 0;
            }
        }
    }
    
    return 1;
}
//This function determines if there is a point in the subregion that overlaps the new point
//in the projection on the xy plane
int Generate_Network::Check_2d_overlapping(const vector<long int> &subregion, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnt_radii, const vector<vector<int> > &global_coordinates, const double &rad_p_dvdw, Point_3D &new_point, bool &overlap_2d)const
{
    //Lambda function to check the distance of two points considering only the
    //x- and y-coordinates
    auto dist_xy2 = [](const Point_3D &P1, const Point_3D &P2){
        double Dx = P1.x - P2.x, Dy = P1.y - P2.y;
        return (Dx*Dx + Dy*Dy);
    };
    
    //Variables to store the z-coordinate of the overlaping point in the xy plane
    double zcoord = 0;
    
    //Scan all points in the subregion
    for (int i = 0; i < (int)subregion.size(); i++) {
        
        //Get the global coordinate
        long int gc = subregion[i];
        
        //Get CNT and point number
        int CNTi = global_coordinates[gc][0];
        int Pj = global_coordinates[gc][1];
        
        //Calculate the squared distance in the xy plane
        double d_xy2 = dist_xy2(new_point, cnts_points[CNTi][Pj]);
        double d_xy = sqrt(d_xy2);
        
        //Calculate the minimum separation needed between the points
        double min_sep = rad_p_dvdw + cnt_radii[CNTi];
        
        //Check if new_point overlaps Pj in the projection on the xy plane
        //Zero is added to eliminate points that are at a distance min_sep
        //(or too close to that distance)
        if (d_xy - min_sep + Zero < Zero) {
            
            //Calculate the total height, i.e., the z-coordinate of new_point if placed "above" Pj
            double z_tot = cnts_points[CNTi][Pj].z + sqrt(min_sep*min_sep - d_xy2);
            
            //Check if this is the first overlapping point found,
            //i.e., the overlapping flag is still false
            if (overlap_2d) {
                
                //At least one overlapping point has been found before
                //Check if the new z-coordinate is at a higher position
                if (z_tot > zcoord) {
                    
                    //Update the z-coordinate and overlapping point
                    zcoord = z_tot;
                }
            }
            //The overlapping flag is false, so this is the first overlapping point
            else {
                
                //Set the z-coordinte and the overlapping point to be the same as that of
                //the current point (P at CNTi)
                zcoord = z_tot;
                
                //Set the overlapping flag to true
                overlap_2d = true;
            }
        }
    }
    
    //Check again if 2D overlapping was found
    if (overlap_2d) {
        
        //Update the z-coordinate of new_point
        new_point.z = zcoord;
    }
    
    return 1;
}
//This function generates a 2D rotation matrix for the initial direction
//The angle theta has a uniform distribution in [0, 2PI]
// R(theta) = |cos(theta)  -sin(theta)|
//          |sin(theta)   cos(theta)|
int Generate_Network::Get_initial_direction_2d(MathMatrix &M, mt19937 &engine_theta, uniform_real_distribution<double> &dist)const
{
    //Generate a rotaion angle theta
    //theta satisfies a uniform distribution in (0, 2PI)
    double theta = 2.0*PI*dist(engine_theta);
    
    //Fill the ratation matrix with a rotation by an anlge phi
    M.element[0][0] = cos(theta);
    M.element[1][0] = sin(theta);
    M.element[1][1] = M.element[0][0];
    M.element[0][1] = -M.element[1][0];
    
    return 1;
}
//This function generates a random direction and a point in the plane xy
int Generate_Network::Get_direction_and_point_2d(const Nanotube_Geo &nanotube_geo, MathMatrix &M, Point_3D &new_point, mt19937 &engine_theta, uniform_real_distribution<double> &dist)const
{
    //Matrix for the new direction
    MathMatrix M_new(2,2);
    
    //Randomly generate a direction
    if(!Get_direction_2d(nanotube_geo, M_new, engine_theta, dist)){
        hout<<"Error in Get_direction_and_point_2d when calling Get_direction_2d"<<endl;
        return 0;
    }
    
    //Calculate the rotation matrix for current segment
    M = M*M_new;
    
    //Get location of new point
    new_point = new_point + Get_new_point_2d(M, nanotube_geo.step_length);
    
    //1 means that point is not the intial point
    new_point.flag = 1;
    
    return 1;
}
Point_3D Generate_Network::Get_new_point_2d(const MathMatrix &M, const double &step)const
{
    //Rotate the 2D vector [step; 0] using the rotation matrix M
    //Note that the matrix M is a 2D rotation matrix, however the output is a 2D point
    //Thus, the z-coordinate of the point is set to 0
    return Point_3D(M.element[0][0]*step,M.element[1][0]*step,0);
}
//This function generates a 2D rotation matrix for the initial direction
//The angle theta has a normal distribution in [-omega, +omega]
// R(theta) = |cos(theta)  -sin(theta)|
//            |sin(theta)   cos(theta)|
int Generate_Network::Get_direction_2d(const Nanotube_Geo &nanotube_geo, MathMatrix &M, mt19937 &engine_theta, uniform_real_distribution<double> &dist)const
{
    //Initialize theta with zero
    double theta = 0.0;
    
    //Get a value for theta that follows a normal distribution in (omega_a, omega_b)
    if (!Get_random_value_mt("normal", engine_theta, dist, nanotube_geo.omega_a, nanotube_geo.omega_b, theta)) {
        hout<<"Error in Get_direction_2d when calling Get_random_value_mt"<<endl;
        return 0;
    }
    
    //Fill the rotation matrix with a rotation by an anlge theta
    M.element[0][0] = cos(theta);
    M.element[0][1] = -M.element[1][0];
    M.element[1][0] = sin(theta);
    M.element[1][1] = M.element[0][0];
    
    return 1;
}
//This function finds the highest position for a new CNT point (not a seed point)
//The highest position is found by rotating the segment from the last point in the new_cnt
//vector to new_point
int Generate_Network::Find_highest_position_for_new_point_iteratively(const Geom_sample &geom_sample_deposit, const Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radii, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const int n_subregions[], const double &cnt_rad, const double &d_vdW, const vector<Point_3D> &new_cnt, Point_3D &new_point, MathMatrix &M, mt19937 &engine_theta, uniform_real_distribution<double> &dist)const
{
    //Variable to count the number of attempts to find the highest position of new_point
    int attempts = 0;
    
    //Iterate while the maximum number of attepts is not reached
    while (attempts < MAX_ATTEMPTS) {
        
        //Matrix for the new direction
        MathMatrix M_new(2,2);
        //hout<<"attempts="<<attempts<<endl;
        
        //Randomly generate a direction. This direction is given by the rotation matrix M_new
        if(!Get_direction_2d(nanotube_geo, M_new, engine_theta, dist)){
            hout<<"Error in Find_highest_position_for_new_point_iteratively when calling Get_direction_2d"<<endl;
            return -1;
        }
        
        //Get a temporary new point
        if (!Get_temporary_new_point_2d(M, M_new, nanotube_geo.step_length, new_cnt.back(), new_point)) {
            hout<<"Error in Find_highest_position_for_new_point_iteratively when calling Get_temporary_new_point_2d"<<endl;
            return -1;
        }
        
        //Find the highest position of new_point
        if (!Find_highest_position_for_new_point(geom_sample_deposit, cnts_points, cnts_radii, global_coordinates, sectioned_domain, n_subregions, cnt_rad, d_vdW, nanotube_geo.step_length, new_cnt, new_point)) {
            hout<<"Error in Find_highest_position_for_new_point_iteratively when calling Find_upmost_position_for_new_point"<<endl;
            return -1;
        }
        
        //Check if there are at least two points in new_cnt (so that there are at leat two CNT
        //semgents in total)
        if (new_cnt.size() >= 2) {
            
            //Calculate the two vectors needed to determine if the last two CNT segments make
            //a valid angle
            Point_3D u = new_point - new_cnt.back();
            Point_3D v = new_cnt.back() - new_cnt[new_cnt.size()-2];
            
            //Check that the last two segments make an angle >= PI/2
            //This is equivalent to the vectors u and v having a positive dot product (or zero)
            if (u.dot(v) >= Zero) {
                
                //The segment makes a valid angle, so update them accumulated rotation
                //matrix and terminate the function with 1
                M = M*M_new;
                return 1;
            }
        }
        else {
            
            //If there is only one point in new_cnt, any position of new_point is valid
            //Thus terminate the function
            return 1;
        }
        
        //If this part of the while-loop is reached, then new_point is in an invalid position
        //Increase the number of attempts for the next iteration where a new point is generated
        attempts++;
    }
    
    return 0;
}
//This function generates a new_point using the previous 2D rotation matrix M and the new 2D
//rotation matrix M_new
//It is a temporary new point beacause new_point may be discarded to generate another one and,
//because of this, the accumulated rotation matix is not calculated nor used to obtain
//new_point. The accumulated rotation matix is calculated once new_point is determined to
//have a valid position
int Generate_Network::Get_temporary_new_point_2d(const MathMatrix &M, const MathMatrix &M_new, const double &step, const Point_3D &P0, Point_3D &new_point)const
{
    //Calculate the x-coordinate of new_point
    new_point.x = P0.x + step*(M.element[0][0]*M_new.element[0][0] + M.element[0][1]*M_new.element[1][0]);
    
    //Calculate the y-coordinate of new_point
    new_point.y = P0.y + step*(M.element[1][0]*M_new.element[0][0] + M.element[1][1]*M_new.element[1][0]);
    
    //Set the z-coordinate of new_point to be the same as that of P0
    new_point.z = P0.z;
    
    return 1;
}
//This function finds the upmost position for a new CNT point (not a seed point)
int Generate_Network::Find_highest_position_for_new_point(const Geom_sample &geom_sample, const vector<vector<Point_3D> > &cnts_points, const vector<double> &cnts_radii, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const int n_subregions[], const double &cnt_rad, const double &d_vdW, const double &step, const vector<Point_3D> &new_cnt, Point_3D &new_point)const
{
    //This variable is defined for code readability
    double floor = geom_sample.sample.poi_min.z;
    
    //Variable to store the upmost overlapping point
    Point_3D p_ovrlp(0.,0.0,floor-cnt_rad);
    
    //Get the new_point subregion
    //For the case of CNT deposit, the cuboids for the sample and extended domain are
    //the same
    //Thus, the function for the random distribution of CNTs can be used for a new CNT point
    //that is not a seed
    //hout<<"Get_cnt_point_subregion"<<endl;
    int subregion = Get_cnt_point_subregion(geom_sample, n_subregions, new_point);
    if (subregion == -1) {
        //new_point is outside the boundary layer, so leave it where it is
        //hout<<"new_point="<<new_point.str()<<endl;
        return 1;
    }
    
    //Radius of new CNT plus van der Waals distance
    double rad_p_dvdw = cnt_rad + d_vdW;
    
    //Calculate the minimum z-coordinate of the sphere with radius cnt_rad with center at new_point
    double zcoord = new_point.z - rad_p_dvdw;
    
    //Variables to rotate a point towards the coordinates of the torus
    double dx = new_point.x - new_cnt.back().x;
    double dy = new_point.y - new_cnt.back().y;
    double h = sqrt(dx*dx + dy*dy);
    double hdx = h/dx;
    double hdy = h/dy;
    
    //Calculate the vector of the previous segment in new_cnt
    Point_3D v_vec;
    
    if (new_cnt.size() >= 2) {
        
        //Use the last two points in the new_cnt vector to obtain P
        v_vec = new_cnt.back() - new_cnt[new_cnt.size()-2];
    }
    else {
        
        //Set vector P to have the same direction as the segment formed by the seed point
        //(which is the only point in the new_cnt vector) and new_point
        v_vec = new_point - new_cnt.back();
    }
    //hout<<endl<<"new_point="<<new_point.str()<<endl;
    
    //Scan all points in the subregion
    for (int i = 0; i < (int)sectioned_domain[subregion].size(); i++) {
        
        //Get the global coordinate
        long int gc = sectioned_domain[subregion][i];
        
        //Get CNT and point number
        int CNTi = global_coordinates[gc][0];
        int Pj = global_coordinates[gc][1];
        //hout<<"i="<<i<<" CNTi="<<CNTi<<" Pj="<<Pj<<endl;
        
        //Check if current point in the subregion:
        //is high enough to penetrate new_point
        //AND
        //is in the same direction as P_vec
        if (cnts_points[CNTi][Pj].z > zcoord - cnts_radii[CNTi] && v_vec.dot(cnts_points[CNTi][Pj] - new_cnt.back()) > Zero) {
            
            Point_3D rotated_p(0.0,0.0,floor-cnt_rad);
            
            //Check if there is penetration of point Pj with new_point
            if (cnts_points[CNTi][Pj].distance_to(new_point) < rad_p_dvdw + cnts_radii[CNTi]) {
                
                //Move new_point to avoid penetration by rotating the new segment
                //Save the result in rotated_p
                //hout<<"     d_NS="<<cnts_points[CNTi][Pj].distance_to(new_point)<<endl;
                if (!Rotate_cnt_segment(rad_p_dvdw+cnts_radii[CNTi], step, cnts_points[CNTi][Pj], new_cnt.back(), new_point, rotated_p)) {
                    hout<<"Error in Find_upmost_position_for_new_point when calling Rotate_cnt_segment"<<endl;
                    return 0;
                }
                //hout<<"     rotated_p="<<rotated_p.str()<<endl;
                //hout<<"     d_NS_new="<<cnts_points[CNTi][Pj].distance_to(rotated_p)<<endl;
            }
            //Check if there is possible penetration of point Pj with new_point if
            //new_point was rotated
            else {
                
                //Calculate distance to torus to find possible penetration
                double dist_torus = 0.0;
                //R is the radius of the donut, i.e., the step length
                //r is the radius of the circle that is sweeped along the center circle
                //of the torus, i.e., the CNT radius
                if (!Calculate_distance_to_torus(cnts_points[CNTi][Pj] - new_cnt.back(), hdx, hdy, step, cnt_rad, dist_torus)) {
                    hout<<"Error in Find_upmost_position_for_new_point when calling Calculate_distance_to_torus"<<endl;
                    return 0;
                }
                //hout<<"     dist_torus="<<dist_torus<<" dist_torus - d_vdW - cnts_radii[CNTi] ="<<dist_torus - d_vdW - cnts_radii[CNTi] <<" cnt_rad="<<cnt_rad<<" r_CNTi="<<cnts_radii[CNTi]<<endl;
                
                //Check if the distance to the torus results in penetration
                if (dist_torus - d_vdW - cnts_radii[CNTi] < Zero) {
                    
                    //If so, find the new position of the point by rotating
                    //the new segment (save it in rotated_p)
                    //hout<<"Rotate_cnt_segment (torus)"<<endl;
                    if (!Rotate_cnt_segment(rad_p_dvdw+cnts_radii[CNTi], step, cnts_points[CNTi][Pj], new_cnt.back(), new_point, rotated_p)) {
                        hout<<"Error in Find_upmost_position_for_new_point when calling Rotate_cnt_segment"<<endl;
                        return 0;
                    }
                    //hout<<"     torus_p="<<rotated_p.str()<<endl;
                }
            }
            
            //Check if new_loc needs updating, which happens when:
            //the z-coordinate of rotated_p is above floor level
            //AND
            //the z-coordinate of rotated_p is higher than the previous value of p_ovrlp
            if (rotated_p.z > floor && rotated_p.z > p_ovrlp.z) {
                
                //Update p_ovrlp
                //hout<<"     Update, rotated_p="<<rotated_p.str()<<endl;
                p_ovrlp = rotated_p;
            }
        }
    }
    
    //Check if penetration was found
    //This happens when the z-coordinate of p_ovrlp is above floor level
    //hout<<"p_ovrlp="<<p_ovrlp.str()<<endl;
    //hout<<"p_ovrlp.z="<<p_ovrlp.z<<" floor="<<floor<<endl;
    //hout<<"p_ovrlp.z - floor="<<p_ovrlp.z - floor<<endl;
    if (p_ovrlp.z - floor > Zero) {
        
        //Penetration was found, so just substitute new_point by p_ovrlp
        //hout<<"     p_ovrlp="<<p_ovrlp.str()<<endl;
        new_point = p_ovrlp;
    }
    else {
        
        //If no overlapping point was found, check if the point should be left with the same
        //z-coordinate of the previous point (i.e., at floor level) or it should be hanging
        if (new_point.z - floor - cnt_rad - Zero > Zero) {
            
            //The CNT point needs to be hanging, determine the new position
            if (!Find_hanging_position(geom_sample.sample, cnt_rad, step, new_cnt, new_point)) {
                hout<<"Error in Find_upmost_position_for_new_point when calling Find_hanging_position"<<endl;
                return 0;
            }
            //hout<<"     p_hang="<<new_point.str()<<endl;
        }
        //else, the point is already at floor level, so leave it as it is
    }
    
    return 1;
}
//This function rotates the new CNT segment 30 degrees (PI/6 radians) around an axis
int Generate_Network::Rotate_cnt_segment_around_axis(const Point_3D &u, const Point_3D &prev, Point_3D &new_point)const
{
    //Calculate axis of rotation, which is u.cross(z), where z is the unit vector along z
    Point_3D r(u.y,-u.x,0.0);
    
    //Calculate the rotation of vector u
    //cos(30) = cos(PI/6) = 0.866025403
    //sin(30) = sin(PI/6) = 0.5
    Point_3D u_rot = u*0.866025403 +r.cross(u)*0.5;
    
    //Calculate new_point from the rotated segment (u_rot) and the previous point in the CNT
    new_point = prev + u_rot;
    
    return 1;
}
//Function that calculates the distance from a point P0 to the torus that is formed by
//all the positions that new_point might take if rotated along the plane where the new
//segment lies and is perpendicular to the xy plane
int Generate_Network::Calculate_distance_to_torus(const Point_3D &P0, const double &hdx, const double &hdy, const double &R, const double &r, double &dist)const
{
    //Map P0 to the coordinate system of the torus
    Point_3D P0_T(hdx*P0.x + hdy*P0.z, hdy*P0.x - hdx*P0.z, P0.y);
    
    //Calculate some quantitites
    //x0^2 + y0^2
    double l_xy2 = P0_T.x*P0_T.x + P0_T.y*P0_T.y;
    //R^2
    double R2 = R*R;
    
    //Check in which region P0 is located
    if (l_xy2 - R2 >= Zero || abs(P0_T.z) - r < Zero) {
        
        //P0 is outside the circle or radius R
        //OR
        //P0 is inside the circle of radius R but its z-coordinate is in [-r, r]
        //Thus P0 is in region 2
        
        //Calculate coordinates of auxiliary circle
        double tmp = R/sqrt(l_xy2);
        double xc = P0_T.x*tmp;
        double yc = P0_T.y*tmp;
        
        //Calculate the distance form P0_T to the center of the auxiliary circle
        double Dx = P0_T.x - xc, Dy = P0_T.y - yc;
        dist = Dx*Dx + Dy*Dy + P0_T.z*P0_T.z;
        dist = sqrt(dist) - r;
    }
    else {
        
        //P0 is in region 1
        dist = abs(P0_T.z) - r;
    }
    
    return 1;
}
//If new_point is penetrating another point (p_point), this function moves new_point to a new
//position by "rotating" the CNT segment that goes from prev_point to new_point
//The rotation is actually done by calculating the new position of new_point constraining the
//plane that contains the CNT segment (from prev_point to new_point) that is perpendicular
//to the plane xy
int Generate_Network::Rotate_cnt_segment(const double &d_new_p, const double &step, const Point_3D &p_point, const Point_3D &prev_point, const Point_3D &new_point, Point_3D &rotated_p)const
{
    //Get the coordinates of new_point taking prev_point as the origin
    Point_3D N = new_point - prev_point;
    
    //Get the penetrating point in the coordinates where prev_point is the origin
    Point_3D S = p_point - prev_point;
    
    //Get the components of the unit vector on the xy plane
    double len_xy = sqrt(N.x*N.x + N.y*N.y);
    double ux = N.x/len_xy;
    double uy = N.y/len_xy;
    
    //Step length squared is used more than once
    double step2 = step*step;
    
    //Calculate "k" constants
    double k1 = -2.0*ux*S.x - 2.0*uy*S.y;
    double k2 = -2*S.z;
    double k3 = step2 + S.length2() - d_new_p*d_new_p;
    
    //Some intermediate quantities
    //Squared k's
    double k1_2 = k1*k1;
    double k2_2 = k2*k2;
    double k3_2 = k3*k3;
    //Length of u
    double len_u_2 = ux*ux + uy*uy;
    //denominator
    double den = k1_2 + k2_2*len_u_2;
    //Common term
    double common = -k1*k3/den;
    //Numerator of "conjugate" term
    double num = k2*sqrt(step2*k1_2 - len_u_2*(k3_2 - step2*k2_2));
    //"Conjugate" term
    double con = num/den;
    
    //Calculate the two possible values for the length along the unit vector on the xy plane
    double a1 = common + con;
    double a2 = common - con;
    
    //Find the largest positive a
    double a;
    if (a1 > Zero && a2 > Zero) {
        
        //Set a to have the value of the largest between a1 and a2
        a = max(a1, a2);
    }
    else {
        
        //One of a1 and a2 is negative, set a to have the positive value
        a = (a1 < Zero)? a2: a1;
    }
    
    //Using a, set the x- and y-coordinates of rotated_p
    rotated_p.x = a*ux;
    rotated_p.y = a*uy;
    
    //Now that a is defined, we can calculate the z-coordinate of rotated_p
    rotated_p.z = sqrt(step2 - a*a*len_u_2);
    
    //Demap rotated_p, so add the previous point
    rotated_p = rotated_p + prev_point;
    //hout<<"     a1="<<a1<<" a2="<<a2<<endl;
    //hout<<"     p_a1="<<(Point_3D(a1*ux,a1*uy,sqrt(step2 - a1*a1*len_u_2))+ prev_point).str()<<endl;
    //hout<<"     p_a2="<<(Point_3D(a2*ux,a2*uy,sqrt(step2 - a2*a2*len_u_2))+ prev_point).str()<<endl;
    
    return 1;
}
//This function finds the position of a new point when the point needs to be hanging
//The new segment hangs half the angle that the previous segment makes with the xy plane
int Generate_Network::Find_hanging_position(const cuboid &sample, const double &cnt_rad, const double &step, const vector<Point_3D> &new_cnt, Point_3D &new_point)const
{
    //Calculate vector u, it's magnitude and u_hat (unit vector along u)
    //u is the vector of the new segment
    Point_3D u = new_point - new_cnt.back();
    //Point_3D u_hat = u.unit();
    //double u_len = u.length();
    
    //Check if there is a previous segment, i.e., if there are at least two points in new_cnt
    if (new_cnt.size() >= 2) {
        
        //Calculate rotation vector
        //Point_3D r = u_hat.cross(Point_3D(0.0,0.0,1.0));
        //hout<<"new_point="<<new_point.str()<<endl;
        //hout<<"d0="<<new_point.distance_to(new_cnt.back());
        
        //Calculate vector v (the vector of the previous segment) and its length
        Point_3D v = new_cnt.back() - new_cnt[new_cnt.size() - 2];
        double v_len = v.length();
        
        //Obtain half the angle of v with the xy plane using the arcsine function
        //Using the arcsine helps recover the sign of the angle
        //According to cplusplus.com, asin returns "Principal arc sine of x,
        //in the interval [-pi/2,+pi/2] radians."
        //PI/10 is around 0.31415926536
        double theta_half = asin(v.z/v_len)*0.5 - 0.31415926536;
        
        //Calculate sine and cosine of theta_half
        double sinT = sin(theta_half);
        double cosT = cos(theta_half);
        
        //Calculate factor that comes from the length of rotation vector
        double rh_len = sqrt(u.x*u.x + u.y*u.y);
        
        //Calculate new vector u
        Point_3D u_new = u*cosT;
        u_new.z = u_new.z + rh_len*sinT;
        
        //Update the location of new_point
        new_point = u_new + new_cnt.back();
    }
    else {
        
        //There is no previous segment, so check if the point is already at floor level
        //Calculate the distance below the previous point where new_point
        //should be located
        //1/sqrt(2)=0.707106781, for a 45 degree (PI/4) angle
        //0.5 for a 30 degree (PI/6) angle
        double dz = 0.5*step;
        
        //Check if moving new_point a distance dz along the direction -z would result
        //in new_point being above floor level (floor = sample.poi_min.z)
        if (new_point.z - dz - cnt_rad - sample.poi_min.z > Zero) {
            
            //Move new_point a distance dz along the direction -z
            
            //First calculate the temporary position along the current segment
            new_point = new_cnt.back() + (new_point - new_cnt.back()).unit()*dz;
            
            //Update the z coordinate
            new_point.z = new_point.z - dz;
        }
        //else: There is no space for the CNT to hang a 30 degree (PI/6) angle,
        //so leave it where it is
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Generate a random value through a probability distribution function
int Generate_Network::Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const
{
    //Check if limits are correctly defined
    if(min > max) {
        hout<<"Error in Get_random_value_mt. The minimum value is larger than the maximum: min = "<<min<<", max = "<<max<<endl;
        return 0;
    }
    
    //Check if the interval has 0 length
    if (abs(max - min) < Zero) {
        //In this case, value is either of the limits
        //To be consistent with the formulation below, value is set equal to min
        value = min;
        return 1;
    }
    
    //Check if uniform distribution
    if(dist_type=="uniform")
    {
        value = (max-min)*dist(engine) + min;
    }
    //Check if normal distribution
    else if(dist_type=="normal")
    {
        double sum=0;
        for(int i=0; i<12; i++)
        {
            sum = sum + dist(engine);
        }
        value = (max-min)*sum/12.0 + min;
    }
    //Wrong input for distribution type
    else {
        hout<<"Error in Get_random_value_mt. The distribution type is not valid. Input was: "<<dist_type<<endl;
        return 0;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//Randomly generate a seed (intial point of a CNT or GNP center) in the given cuboid
int Generate_Network::Get_point_in_cuboid_mt(const cuboid &cub, Point_3D &point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const
{
    
    point.x = cub.poi_min.x + cub.len_x*dist(engine_x);
    
    point.y = cub.poi_min.y + cub.wid_y*dist(engine_y);
    
    point.z = cub.poi_min.z + cub.hei_z*dist(engine_z);
    
    //Zero flag denotes this point is the initial point of a CNT
    point.flag = 0;
    
    return 1;
}
//---------------------------------------------------------------------------
//Randomly generate an initial direction, then generate the rotation matrix that results in that rotation
int Generate_Network::Get_initial_direction_mt(const string &dir_distrib_type, const double &ini_theta, const double &ini_phi, mt19937 &engine_inital_direction, uniform_real_distribution<double> &dist_initial, MathMatrix &rotation)const
{
    if(dir_distrib_type=="random")
    {
        //Choose three random numbers between -1 and 1,
        //they are the components of the vector that will define the initial direction
        double a = dist_initial(engine_inital_direction);
        double b = dist_initial(engine_inital_direction);
        double c = dist_initial(engine_inital_direction);
        
        //Check that a, b and c are not all zero
        if (abs(a) < Zero && abs(b) < Zero && abs(c) < Zero) {
            //Then transform this to a simple case where a = b = c
            a = 1.0;
            b = 1.0;
            c = 1.0;
        }
        
        //Calculate the length of the vector v = (a,b,c)
        double v_length = sqrt(a*a + b*b + c*c);
        
        //This quantity is used three times:
        double quantity = sqrt(a*a + b*b);
        
        //Calculate the trigonometric functions of the angles theta and phi
        double cos_phi = a/quantity;
        double sin_phi = b/quantity;
        double cos_theta = c/v_length;
        double sin_theta = quantity/v_length;
        
        //Fill the elements of the rotation matrix
        rotation.element[0][0] = cos_phi*cos_theta;
        rotation.element[0][1] = -sin_phi;
        rotation.element[0][2] = cos_phi*sin_theta;
        
        rotation.element[1][0] = sin_phi*cos_theta;
        rotation.element[1][1] = cos_phi;
        rotation.element[1][2] = sin_phi*sin_theta;
        
        rotation.element[2][0] = -sin_theta;
        rotation.element[2][2] = cos_theta;
    }
    else if(dir_distrib_type=="specific")
    {
        //initialize variables with the initial direction
        double cnt_theta = ini_theta;
        double cnt_phi = ini_phi;
        
        //Use the probability of a random number to be even to use the opposite direction half the time
        if( engine_inital_direction()%2==0 )
        {
            //Invert the direction
            cnt_theta = PI - ini_theta;
            cnt_phi = PI + ini_phi;
        }
        //Get the rotation matrix
        rotation = Get_transformation_matrix(cnt_theta, cnt_phi);
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//This function generates the angles theta and phi
//Theta defines the curvature of the CNT, phi defines the growth direction
int Generate_Network::Get_direction_normal_distribution(const Nanotube_Geo &nanotube_geo, double &cnt_theta, double &cnt_phi, mt19937 &engine_theta, mt19937 &engine_phi, uniform_real_distribution<double> &dist)const
{
    //Get a value for theta that follows a normal distribution in (omega_a, omega_b)
    if (!Get_random_value_mt("normal", engine_theta, dist, nanotube_geo.omega_a, nanotube_geo.omega_b, cnt_theta)) {
        hout<<"Error in Get_direction_normal_distribution when calling Get_random_value_mt"<<endl;
        return 0;
    }
    //angle theta is always positive
    cnt_theta = abs(cnt_theta);
    
    //phi satisfies a uniform distribution in (0, 2PI)
    cnt_phi = 2.0*PI*dist(engine_phi);
    
    return 1;
}
//---------------------------------------------------------------------------
//Get rotation matrix from two angles
MathMatrix Generate_Network::Get_transformation_matrix(const double &theta, const double &phi)const
{
    //M = M_phi*M_theta
    //          |cos(phi) -sin(phi) 0|
    // M_phi  = |sin(phi)  cos(phi) 0|
    //          |   0         0     1|
    //
    //           | cos(theta)  0  sin(theta)|
    // M_theta = |      0      1      0     |
    //           |-sin(theta)  0  cos(theta)|
    //Calculate the matrix elements directly, instead of multiplying two rotation matrices
    MathMatrix M(3,3);
    M.element[0][0] = cos(phi)*cos(theta);
    M.element[0][1] = -sin(phi);
    M.element[0][2] = cos(phi)*sin(theta);
    
    M.element[1][0] = sin(phi)*cos(theta);
    M.element[1][1] = cos(phi);
    M.element[1][2] = sin(phi)*sin(theta);
    
    M.element[2][0] = -sin(theta);
    M.element[2][2] = cos(theta);
    
    return M;
}
//---------------------------------------------------------------------------
//Calculate the coordinates of the new CNT point (transformation of coordinates)
Point_3D Generate_Network::Get_new_point(MathMatrix &Matrix, const double &Rad)const
{
    //Point = Matrix*v
    //v = [0; 0; Rad]
    //Calculate the new point directly
    Point_3D Point(Matrix.element[0][2]*Rad, Matrix.element[1][2]*Rad, Matrix.element[2][2]*Rad);
    
    return Point;
}
//This function checks if a point is inside a cuboid
int Generate_Network::Is_point_inside_cuboid(const cuboid &cub, const Point_3D &point)const
{
    if(point.x<cub.poi_min.x||point.x>cub.max_x||
       point.y<cub.poi_min.y||point.y>cub.max_y||
       point.z<cub.poi_min.z||point.z>cub.max_z) {
              //Point is outside cuboid, return false (0)
              return 0;
          }
    
    return 1;
}
//This function is used to display the percentage of generated volume in increments
//that are multiple of 10%
int Generate_Network::Check_progress(const string &particle, const int &elapsed_time, const double &target_vol, const double &generated_vol, double &vol_completed, double &vol_completed_acc)const
{
    if (generated_vol > vol_completed_acc) {
        
        //Find the actual increment
        while (generated_vol > vol_completed_acc) {
            
            //Increase the volume used to check progress
            vol_completed = vol_completed + 0.1;
            vol_completed_acc = vol_completed*target_vol;
        }
        
        //Output elapsed time
        hout << "Completed " << (vol_completed - 0.1)*100 << " % of " << particle << " target volume. Elapsed time: " << elapsed_time << " secs." <<endl;
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//Generate a GNP network
int Generate_Network::Generate_gnp_network_mt(const Simu_para &simu_para, const GNP_Geo &gnp_geo, const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, vector<GNP> &gnps, vector<vector<int> > &sectioned_domain, double &gnp_vol_tot, double &gnp_wt_tot)const
{

    //geom_sample cannot be modified, so copy the seeds to an array if they were specified in
    //the input file otherwise generate them
    unsigned int net_seeds[6];
    if (!GNP_seeds(simu_para.GNP_seeds, net_seeds)) {
        hout<<"Error in GNP_seeds"<<endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------
    //Generate differnet engines for different variables
    std::mt19937 engine_x(net_seeds[0]);
    std::mt19937 engine_y(net_seeds[1]);
    std::mt19937 engine_z(net_seeds[2]);
    std::mt19937 engine_l(net_seeds[3]);
    std::mt19937 engine_t(net_seeds[4]);
    std::mt19937 engine_orientation(net_seeds[5]);
    
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [-1, 1].
    std::uniform_real_distribution<double> dist_orientation(-1.0, 1.0);
    
    //---------------------------------------------------------------------------
    //Total volume of generated GNPs, initialized at zero
    gnp_vol_tot = 0;
    
    //Variable to count the number of GNPs that were deleted due to penetration
    int gnp_reject_count = 0;
    //Variable to count the number of GNPs ignored because they were not even partially inside the sample
    int gnp_ignored_count = 0;
    
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    //n_subregions[2] is the number of subregions along z
    int n_subregion[] = {0,0,0};
    //Calculate the number of subregins along each direction and intialize the sectioned_domain vector
    //hout<<"Initialize_gnp_subregions"<<endl;
    if (!Initialize_gnp_subregions(geom_sample, n_subregion, sectioned_domain)) {
        hout<<"Error in Initialize_gnp_subregions"<<endl;
        return 1;
    }
    
    //Time variables to keep track of generation time
    time_t ct0, ct1;
    //Get the time when generation started
    ct0 = time(NULL);
    
    //Variable used to check when 10% is completed
    double vol_completed = 0.1;
    //Variable used to store the fraction of CNT volume indicated by vol_completed
    double vol_completed_acc = vol_completed*gnp_geo.volume;
    
    //---------------------------------------------------------------------------
    while( gnp_vol_tot < gnp_geo.volume )
    {
        //---------------------------------------------------------------------------
        //Generate a GNP
        GNP gnp;
        //Location of the GNP center
        Point_3D gnp_poi;
        
        //---------------------------------------------------------------------------
        //Randomly generate a direction as the orientation of the GNP
        //hout<<"Get_initial_direction_mt"<<endl;
        if (!Get_initial_direction_mt(gnp_geo.orient_distrib_type, gnp_geo.ini_theta, gnp_geo.ini_phi, engine_orientation, dist_orientation, gnp.rotation)) {
            hout << "Error in Generate_gnp_network_mt when calling Get_initial_direction_mt" << endl;
            return 0;
        }
        
        //---------------------------------------------------------------------------
        //Randomly generate a point inside the extended domain, this will be the displacement
        //applied to the GNP, i.e, its random location
        //hout<<"Get_seed_point_mt"<<endl;
        if(Get_point_in_cuboid_mt(geom_sample.ex_dom_gnp, gnp.center, engine_x, engine_y, engine_z, dist)==0) return 0;
        //hout<<"gnp.center="<<gnp.center.str()<<endl;
        
        //Generate the GNP and all its remaining geometric parameters (vertices, planes, normal vectors)
        //hout<<"Generate_gnp"<<endl;
        if (!Generate_gnp(gnp_geo, gnp, engine_l, engine_t, dist)) {
            hout << "Error in Generate_gnp_network_mt when calling Generate_gnp" << endl;
            return 0;
        }
        
        //Flag to determine if new gnp was rejected
        bool rejected = false;
        
        //Variable to store all subregions a GNP belongs to
        set<int> subregions_gnp;
        
        //Check if we are allowing penetrations
        if (simu_para.penetration_model_flag) {
            //Penetrations are not allowed
            
            //Check if the GNP penetrates another GNP
            //hout<<"Deal_with_gnp_interpenetrations"<<endl;
            if (!Deal_with_gnp_interpenetrations(geom_sample, cutoffs, gnps, n_subregion, gnp, subregions_gnp, sectioned_domain, rejected)) {
                hout<<"Error in Deal_with_gnp_interpenetrations"<<endl;
                return 0;
            }
        }
        
        //Check if the new GNP was rejected
        //If the penetrating model is used, then rejected is always false and the volume of gnp_new
        //is always added
        if (rejected) {
            
            //The GNP was rejected so increase the count of rejected GNPs
            //hout<<"GNP WAS REJECTED"<<endl;
            gnp_reject_count++;
        }
        else {
            //The GNP was not rejected, so update generated volume/weight and add it to the
            //vector of GNPs
            
            //Flag to determine if a GNP can be added or not into the vector of GNPs
            //This will depend on the minimum required volume of a GNP inside the sample
            bool is_enough_inside = false;
            
            //Check what is the criterion for the minimum GNP volume inside the sample
            if (gnp_geo.vol_in == "all_in") {
                
                //Check that all vertices of the GNP are inside the sample
                //hout<<"Check_all_gnp_vertices_are_inside_sample"<<endl;
                if (!Check_all_gnp_vertices_are_inside_sample(geom_sample, gnp, is_enough_inside)) {
                    hout<<"Error in Check_all_gnp_vertices_are_inside_sample."<<endl;
                    return 0;
                }
                //hout<<"is_enough_inside="<<is_enough_inside<<endl;
            }
            else {
                
                //Calculate (or approximate) the generated volume
                //hout<<"Calculate_generated_gnp_vol_and_update_total_vol"<<endl;
                if (!Calculate_generated_gnp_vol(gnp_geo, geom_sample, gnp, is_enough_inside)) {
                    hout<<"Error in Calculate_generated_gnp_vol."<<endl;
                    return 0;
                }
                //hout<<"gnp.volume="<<gnp.volume<<endl;
            }
            
            //Only add the GNP if it is at least partially inside
            if (is_enough_inside) {
                
                //Current GNP will not be rejected nor ignored, so add it to the
                //corresponding subregions
                if (!Add_valid_gnp_to_subregions((int)gnps.size(), subregions_gnp, sectioned_domain)) {
                    hout<<"Error in Add_valid_gnp_to_subregions after displacement."<<endl;
                    return 0;
                }
                
                //Update the GNP flag
                gnp.flag = (int)gnps.size();
                
                //hout << "Total volume = " << gnp_vol_tot << ". Approximated volume = " << gnp.volume <<endl;
                //Update the total volume
                gnp_vol_tot = gnp_vol_tot + gnp.volume;
                
                //Add the current particle to the vector of particles
                gnps.push_back(gnp);
            }
            else {
                
                //The GNP is completely outside the sample or not of it enough inside the sample
                //so increase the count of ignored GNPs
                gnp_ignored_count++;
            }
        }//End of if for rejected flag
        
        //Get the time to check progress
        ct1 = time(NULL);
        
        //Check progress
        if (!Check_progress("GNP", (int)(ct1-ct0), gnp_geo.volume, gnp_vol_tot, vol_completed, vol_completed_acc)) {
            hout<<"Error when calculating the percentage of GNP volume generated"<<endl;
            return 0;
        }
        
        //while-loop ends here
    }
    
    if (simu_para.particle_type == "GNP_cuboids" || simu_para.particle_type == "GNP_CNT_mix") {
        
        //Print the amount of generated GNPs when GNPs only or mixed fillers are generated
        if(gnp_geo.criterion == "wt") {
            
            //Update the weight from the volume
            gnp_wt_tot = gnp_vol_tot*gnp_geo.density;
            
            //Calculate matrix weight
            double matrix_weight = (geom_sample.volume - gnp_vol_tot)*geom_sample.matrix_density;
            
            hout << endl << "The weight fraction of generated GNPs inside the sample is: " << gnp_wt_tot/(matrix_weight + gnp_wt_tot) << ", the target weight fraction was " << gnp_geo.weight_fraction << endl;
        }
        else if(gnp_geo.criterion == "vol") {
            
            hout << endl << "The volume fraction of generated GNPs inside the sample is: " << gnp_vol_tot/geom_sample.volume;
            hout << ", the target volume fraction was " << gnp_geo.volume_fraction << endl;
        }
        
        hout << "There are " << (int)gnps.size()<< " GNPs inside the sample domain. ";
        //If there were rejected GNPs, also output that information
        if (gnp_reject_count || gnp_ignored_count) {
            hout<<"In total "<<(int)gnps.size()+gnp_ignored_count+gnp_reject_count<<" GNPs were generated, from which "<<gnp_reject_count<<" were rejected and "<<gnp_ignored_count<<" were ignored.";
        }
        hout<< endl << endl ;
    }
    
    return 1;
}
//This function copies the seeds for the random number generators given in the input file
//If seeds were not given, seeds are generated
int Generate_Network::GNP_seeds(const vector<unsigned int> &GNP_seeds, unsigned int net_seeds[])const
{
    //---------------------------------------------------------------------------
    //Random_device is used to generate seeds for the Mersenne twister
    std::random_device rd;

    //If seeds have been specified, copy them
    if (GNP_seeds.size()) {
        for (size_t i = 0; i < GNP_seeds.size(); i++) {
            net_seeds[i] = GNP_seeds[i];
        }
    }
    //If seeds have not been specified, generate them
    else {
        //Generate all the new seeds
        net_seeds[0] = rd();
        //hout << "seed x: "<<net_seeds[0]<<endl;
        net_seeds[1] = rd();
        //hout << "seed y: "<<net_seeds[1]<<endl;
        net_seeds[2] = rd();
        //hout << "seed z: "<<net_seeds[2]<<endl;
        net_seeds[3] = rd();
        //hout << "seed length: "<<net_seeds[3]<<endl;
        net_seeds[4] = rd();
        //hout << "seed thickness: "<<net_seeds[4]<<endl;
        net_seeds[5] = rd();
        //hout << "seed orientation: "<<net_seeds[5]<<endl;
    }
    
    //Output the seeds
    hout<<"GNP seeds:"<<endl;
    for (int i = 0; i < 6; i++) {
        hout<<net_seeds[i]<<' ';
    }
    
    hout<<endl;
    return 1;
}
int Generate_Network::Initialize_gnp_subregions(const Geom_sample &sample_geom, int n_subregion[], vector<vector<int> > &sectioned_domain)const
{   
    //Calculate the number of variables along each direction
    //Make sure there is at least one subregion along each direction
    //
    //Number of subregions along x
    n_subregion[0] = max(1, (int)(sample_geom.sample.len_x/sample_geom.gs_minx));
    //Number of subregions along y
    n_subregion[1] = max(1, (int)(sample_geom.sample.wid_y/sample_geom.gs_miny));
    //Number of subregions along z
    n_subregion[2] = max(1, (int)(sample_geom.sample.hei_z/sample_geom.gs_minz));
    
    //Initialize sectioned_domain
    sectioned_domain.assign(n_subregion[0]*n_subregion[1]*n_subregion[2], vector<int>());
    
    return 1;
}
//This fuction generates a GNP with a square base and generates other usefuld geometric parameters:
//its volume, the equations of the planes for the faces, thu unit normal vectors to the planes,
//and the vertices
int Generate_Network::Generate_gnp(const GNP_Geo &gnp_geo, GNP &gnp, mt19937 &engine_l, mt19937 &engine_t, uniform_real_distribution<double> &dist)const
{
    //Define side length for the squared GNP surface
    gnp.l = gnp_geo.len_min + (gnp_geo.len_max - gnp_geo.len_min)*dist(engine_l);
    
    //Thickness
    gnp.t = gnp_geo.t_min + (gnp_geo.t_max - gnp_geo.t_min)*dist(engine_t);
    
    //Calculate volume
    gnp.volume = gnp.l*gnp.l*gnp.t;
    
    //Calcualte the coordinates of the GNP vertices
    if (!Obtain_gnp_vertex_coordinates(gnp)) {
        hout << "Error in Generate_gnp when calling Obtain_gnp_vertex_coordinates." << endl;
        return 0;
    }
    
    //Get the plane equations for the six faces
    if (!Update_gnp_plane_equations(gnp)) {
        hout<<"Error in Generate_gnp when calling Update_gnp_plane_equations."<<endl;
        return 0;
    }
    
    return 1;
}
//This function calculates the vertex coordinates of a GNP
//Here, it is assumed that the length, thickness, rotation matrix, and center are known
int Generate_Network::Obtain_gnp_vertex_coordinates(GNP& gnp)const
{
    //The GNP in local coordinates is assumed to have its center at the origin
    //Thus its vertices are at +/- l/2 in x and y, and +/- t in z in local coordinates
    //After applying the displacement and rotation (that are already in the gnp variable),
    //vertices are mapped to the global coordinates
    //
    //Top surface
    gnp.vertices[0] = (Point_3D(gnp.l / 2, gnp.l / 2, gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[1] = (Point_3D(-gnp.l / 2, gnp.l / 2, gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[2] = (Point_3D(-gnp.l / 2, -gnp.l / 2, gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[3] = (Point_3D(gnp.l / 2, -gnp.l / 2, gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    //Bottom surface
    gnp.vertices[4] = (Point_3D(gnp.l / 2, gnp.l / 2, -gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[5] = (Point_3D(-gnp.l / 2, gnp.l / 2, -gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[6] = (Point_3D(-gnp.l / 2, -gnp.l / 2, -gnp.t / 2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[7] = (Point_3D(gnp.l / 2, -gnp.l / 2, -gnp.t / 2)).rotation(gnp.rotation, gnp.center);

    return 1;
}
//This function updates the plane equations for a fully generated GNP
int Generate_Network::Update_gnp_plane_equations(GNP &gnp)const
{
    //Calculate the plane normals (as unit vectors) and their equations
    //Normal is calcualted as (P2-P1)x(P3-P1)
    //
    //Top face
    gnp.faces[0] = Plane_3D(gnp.vertices[0], gnp.vertices[1], gnp.vertices[3]);
    //Bottom face
    gnp.faces[1] = Plane_3D(gnp.vertices[4], gnp.vertices[7], gnp.vertices[5]);
    //Front face
    gnp.faces[2] = Plane_3D(gnp.vertices[0], gnp.vertices[3], gnp.vertices[4]);
    //Right face
    gnp.faces[3] = Plane_3D(gnp.vertices[0], gnp.vertices[4], gnp.vertices[1]);
    //Back face
    gnp.faces[4] = Plane_3D(gnp.vertices[1], gnp.vertices[5], gnp.vertices[2]);
    //Left face
    gnp.faces[5] = Plane_3D(gnp.vertices[2], gnp.vertices[6], gnp.vertices[3]);
    
    return 1;
}
//This function checks whether the newly generated GNP penetrates another GNP
//If there is interpenetration, then the GNP is moved
//If the attempts to relocate the GNP exceed the number of maximum attempts, the GNP is rejected
int Generate_Network::Deal_with_gnp_interpenetrations(const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, const vector<GNP> &gnps, const int n_subregions[], GNP &gnp_new, set<int> &subregions_gnp, vector<vector<int> > &sectioned_domain, bool &rejected)const
{
    //Variable to count the number of attempts
    int attempts = 0;
    
    //Flag to determine if gnp_new was moved
    bool displaced;
    
    //hout<<endl<<"GNP_new ="<<gnps.size()<<" center="<<gnp_new.center.str()<<endl;
    //Keep moving gnp_new as long as the number of attempts does not exceed the maximum allowed
    while (attempts <= MAX_ATTEMPTS) {
        
        //Empty the subregions set
        subregions_gnp.clear();
        
        //Find the subregions gnp_new occupies
        //hout<<"attempts="<<attempts<<endl;
        if (!Get_gnp_subregions(geom_sample, gnp_new, n_subregions, subregions_gnp)) {
            hout<<"Error in Get_gnp_subregions"<<endl;
            return 0;
        }
        
        //Scan the subregions the GNP occupies to determine if it might be close enough to another GNP
        //Get all GNPs it might penetrate
        set<int> gnp_set;
        //hout<<"Get_gnps_in_subregions"<<endl;
        if (!Get_gnps_in_subregions(sectioned_domain, subregions_gnp, gnp_set)) {
            hout<<"Error in Get_gnps_in_subregions"<<endl;
            return 0;
        }
        //hout<<"gnp_set.size="<<gnp_set.size()<<endl;
        
        //Check if close enough GNPs were found
        if (!gnp_set.empty()) {
            
            //Determine if gnp_new needs to be moved and, if so, move it
            //hout<<"Move_gnps_if_needed"<<endl;
            if (!Move_gnps_if_needed(cutoffs, gnps, gnp_set, gnp_new, displaced)) {
                hout<<"Error in Move_gnps_if_needed."<<endl;
                return 0;
            }
            
            //Check if the GNP_new was not moved
            if (!displaced) {
                
                //If gnp_new was not moved, then it is in a valid position
                //Thus, set the rejected flag to false
                rejected = false;
                
                return 1;
            }
            //If the GNP was moved, then one more iteration is needed to check that the new
            //position is a valid position
            //else{hout<<"displaced gnp_new.center="<<gnp_new.center.str()<<endl;}
        }
        else {
            
            //If there are no GNPs that could be too close or interpenetrating gnp_new,
            //then gnp_new is in a valid position
            //Thus, set the rejected flag to false
            rejected = false;
            
            //terminate the function
            return 1;
        }
        
        //Increase the number of attempts
        attempts++;
    }
    
    //If this part of the code is reached, then gnp_new was moved and we need to check if it is
    //in a valid positon. However, since the maximum number of attempts was reached then gnp_new
    //is to be rejected. So set the rejected flag to true
    rejected = true;
    
    return 1;
}
//This function finds all the subregions the GNP_new occupies and those subregions where there might be
//close enough GNPs below the van der Waals distance
int Generate_Network::Get_gnp_subregions(const Geom_sample &geom_sample, const GNP &gnp_new, const int n_subregions[], set<int> &subregions)const
{
    //Number of points to discretize the GNP along the x direction (of the GNP local coordinates)
    int n_points_x = max(2, 1 + (int)(gnp_new.l/geom_sample.gs_minx));
    
    //Number of points to discretize the GNP along the y direction (of the GNP local coordinates)
    int n_points_y = max(2, 1 + (int)(gnp_new.l/geom_sample.gs_miny));
    
    //Number of points to discretize the GNP along the z direction (of the GNP local coordinates)
    //Make sure there are at least two points in the discretization along z
    int n_points_z = max(2, 1 + (int)(gnp_new.t/geom_sample.gs_minz));
    
    //Iterate over all points in the discretization
    //Index k moves the point in the discretization along z
    for (int k = 0; k < n_points_z; k++) {
        
        //Calculate the start and end points along the y-direction
        
        //Their locations are proportional to index k and the number of points
        //in the discretization along z
        double lambda_z = (double)k/(n_points_z-1);
        
        //These points move along edges 3-7 and 0-4 of the GNP
        Point_3D start1_y = gnp_new.vertices[7] + (gnp_new.vertices[3] - gnp_new.vertices[7])*lambda_z;
        Point_3D end1_y = gnp_new.vertices[4] + (gnp_new.vertices[0] - gnp_new.vertices[4])*lambda_z;
        Point_3D diff1_y = end1_y - start1_y;
        
        //These points move along edges 2-6 and 1-5 of the GNP
        Point_3D start2_y = gnp_new.vertices[6] + (gnp_new.vertices[2] - gnp_new.vertices[6])*lambda_z;
        Point_3D end2_y = gnp_new.vertices[5] + (gnp_new.vertices[1] - gnp_new.vertices[5])*lambda_z;
        Point_3D diff2_y = end2_y - start2_y;
        
        //Index j moves the point in the discretization along y
        for (int j = 0; j < n_points_y; j++) {
            
            //Calculate the start and end points along the x-direction
            
            //Their locations are proportional to index j and the number of points
            //in the discretization along y
            double lambda_y = (double)j/(n_points_y-1);
            
            //These points move along the y-direction from start_y to end_y
            Point_3D start_x = start1_y + diff1_y*lambda_y;
            Point_3D end_x = start2_y + diff2_y*lambda_y;
            Point_3D diff_x = end_x - start_x;

            //Index i moves the point in the discretization along x
            for (int i = 0; i < n_points_x; i++) {
                
                //Calculate the position of new_point (the point in the discretization)
                
                //The location of new_point is proportional to index i and the number of points
                //in the discretization along x
                double lambda_x = (double)i/(n_points_x-1);
                
                //Location of new_point
                Point_3D new_point = start_x + diff_x*lambda_x;
                

                //Check if new_point is inside the sample
                if (Is_point_inside_cuboid(geom_sample.sample, new_point)) {
                    
                    //If the point is inside the sample, then add the subregion the GNP point occupies
                    if (!Add_gnp_subregions_to_set_for_gnp_point(geom_sample, new_point, n_subregions, subregions)) {
                        hout<<"Error in Add_gnp_subregions_to_set_for_gnp_point"<<endl;
                        return 0;
                    }
                }
            }
            
        }
    }
    
    return 1;
}
int Generate_Network::Add_gnp_subregions_to_set_for_gnp_point(const Geom_sample &geom_sample, const Point_3D &new_point, const int n_subregions[], set<int> &subregions)const
{
    //Minimum and maximum coordinates of the overlapping subregions
    //Initialize them with the maximum and minimum values, respectively, so they can be updated
    int min_a = n_subregions[0], max_a = 0;
    int min_b = n_subregions[1], max_b = 0;
    int min_c = n_subregions[2], max_c = 0;
    
    //Calculate the region-coordinates a, b, c
    int a = (int)((new_point.x-geom_sample.sample.poi_min.x)/geom_sample.gs_minx);
    //Limit the value of a as it has to go from 0 to n_subregions[0]-1
    if (a == n_subregions[0]) a--;
    int b = (int)((new_point.y-geom_sample.sample.poi_min.y)/geom_sample.gs_miny);
    //Limit the value of b as it has to go from 0 to n_subregions[1]-1
    if (b == n_subregions[1]) b--;
    int c = (int)((new_point.z-geom_sample.sample.poi_min.z)/geom_sample.gs_minz);
    //Limit the value of c as it has to go from 0 to n_subregions[2]-1
    if (c == n_subregions[2]) c--;
    //hout<<"a="<<a<<" b="<<b<<" c="<<c<<endl;
    
    //Update the maximum and minimum values of coordinates
    if (a < min_a) { min_a = a; }
    if (a > max_a) { max_a = a; }
    if (b < min_b) { min_b = b; }
    if (b > max_b) { max_b = b; }
    if (c < min_c) { min_c = c; }
    if (c > max_c) { max_c = c; }
    
    //Calculate the coordinates of non-overlaping region the point belongs to
    double x1 = (double)a*geom_sample.gs_minx +  geom_sample.sample.poi_min.x;
    double x2 = x1 + geom_sample.gs_minx;
    double y1 = (double)b*geom_sample.gs_miny +  geom_sample.sample.poi_min.y;
    double y2 = y1 + geom_sample.gs_miny;
    double z1 = (double)c*geom_sample.gs_minz +  geom_sample.sample.poi_min.z;
    double z2 = z1 + geom_sample.gs_minz;
    
    //Update subregion coordinates in case the point is in the overlapping part of the subregion
    //The first operand eliminates the periodicity on the boundary
    if ((a > 0) && (new_point.x >= x1) && (new_point.x <= x1+geom_sample.gs_overlap_gnp)) {
        //The point is in the "left" side of the overlapping region (along the x axis),
        //then the minimum a-coordinate needs to be updated again
        min_a = min_a - 1;
    }
    else if ((a < n_subregions[0]-1) && (new_point.x >= x2-geom_sample.gs_overlap_gnp) && (new_point.x <= x2 )) {
        //The point is in the "right" side of the overlapping region (along the x axis),
        //then the maximum a-coordinate needs to be updated again
        max_a = max_a + 1;
    }
    if ((b > 0) && (new_point.y >= y1) && (new_point.y <= y1+geom_sample.gs_overlap_gnp)) {
        //The point is in the "left" side of the overlapping region (along the y axis),
        //then the minimum a-coordinate needs to be updated again
        min_b = min_b - 1;
    }
    else if ((b < n_subregions[1]-1) && (new_point.y >= y2-geom_sample.gs_overlap_gnp) && (new_point.y <= y2 )){
        //The point is in the "right" side of the overlapping region (along the y axis),
        //then the maximum a-coordinate needs to be updated again
        max_b = max_b + 1;
    }
    if ((c > 0) && (new_point.z >= z1) && (new_point.z <= z1+geom_sample.gs_overlap_gnp)) {
        //The point is in the "left" side of the overlapping region (along the z axis),
        //then the minimum a-coordinate needs to be updated again
        min_c = min_c - 1;
    }
    else if ((c < n_subregions[2]-1) && (new_point.z >= z2-geom_sample.gs_overlap_gnp) && (new_point.z <= z2 )) {
        //The point is in the "right" side of the overlapping region (along the z axis),
        //then the maximum a-coordinate needs to be updated again
        max_c = max_c + 1;
    }
    
    //Generate all combinations of the coordinates to find the subregions the GNP occupies
    //as given the coordinates of new_point
    /*hout<<"min_a="<<min_a<<" max_a="<<max_a<<endl;
    hout<<"min_b="<<min_b<<" max_b="<<max_b<<endl;
    hout<<"min_c="<<min_c<<" max_c="<<max_c<<endl;
    hout<<"subregions:"<<endl;*/
    for (int ii = min_a; ii <= max_a; ii++) {
        for (int jj = min_b; jj <= max_b; jj++) {
            for (int kk = min_c; kk <= max_c; kk++) {
                
                //Calculate the subregion index
                int idx = ii + jj*n_subregions[0] + kk*n_subregions[0]*n_subregions[1];
                //hout<<idx<<' ';
                
                //Add the index to the vector of subregions
                subregions.insert(idx);
            }
        }
    }
    //hout<<endl;
    
    return 1;
}
//From a given list of subregions, this function fids all GNPs in those subregions
//A set is the putput value to avoid repetition of GNP numbers
int Generate_Network::Get_gnps_in_subregions(const vector<vector<int> > &sectioned_domain, const set<int> &subregions, set<int> &gnp_set)const
{
    //hout<<"sectioned_domain.size="<<sectioned_domain.size()<<endl;
    for (set<int>::iterator i = subregions.begin(); i != subregions.end(); i++) {
        
        //Current sub-region
        int s = *i;
        
        //hout<<"   sectioned_domain["<<s<<"].size="<<sectioned_domain[s].size()<<endl;
        for (size_t j = 0; j < sectioned_domain[s].size(); j++) {
            
            //Current GNP
            int GNP = sectioned_domain[s][j];
            
            //Add current GNP to set
            gnp_set.insert(GNP);
        }
    }
    
    return 1;
}
//This function moves a GNP if needed
int Generate_Network::Move_gnps_if_needed(const Cutoff_dist &cutoffs, const vector<GNP> &gnps, set<int> &gnp_set, GNP &gnp_new, bool &displaced) const
{
    //Create a new Collision_detection object to determine whether two GNPs interpenetrate
    //eachother or not, and their distances if they do not interpenetrate each other
    Collision_detection GJK_EPA;
    
    //Initialize the displaced flag to false
    displaced = false;
    
    //Variable to store the two largest displacements
    vector<double> disps;
    //Variable to store the displacement vectors of the two largest displacements
    vector<Point_3D> disps_vec;
    
    //int el = 0;
    for (set<int>::iterator i = gnp_set.begin(); i!=gnp_set.end(); i++) {
        
        //Get current GNP number
        int GNP_i = *i;
        //hout<<"GNP_i="<<GNP_i<<" gnps.size()="<<gnps.size()<<endl;
        /*if (gnps.size() == 33 && GNP_i == 9) {
            hout<<"===================="<<endl;
            hout<<"===================="<<endl;
            hout<<"GNP_i:"<<endl;
            for (int i = 0; i < 8; i++) {
                hout<<gnps[GNP_i].vertices[i].x<<' '<<gnps[GNP_i].vertices[i].y<<' '<<gnps[GNP_i].vertices[i].z<<endl;
            }
            hout<<"GNP_new:"<<endl;
            for (int i = 0; i < 8; i++) {
                hout<<gnp_new.vertices[i].x<<' '<<gnp_new.vertices[i].y<<' '<<gnp_new.vertices[i].z<<endl;
            }
        }//*/
        
        //Vector to stor the simplex that encloses the origin in case of interpenetration
        vector<Point_3D> simplex;
        
        //Flags for penetration (p_flag) and touching (t_flag)
        bool p_flag = false;
        bool t_flag = false;
        
        //Check if GNP and GNP_new penetrate each other
        //hout<<"GJK_EPA.GJK"<<endl;
        if (!GJK_EPA.GJK(gnps[GNP_i], gnp_new, simplex, p_flag, t_flag)) {
            hout<<"Error in Move_gnps_if_needed when calling GJK"<<endl;
            return 0;
        }
        //hout<<"GJK_EPA.GJK end"<<endl;

        //Variables to store the penetration depth (PD) and direction vector (N) along which
        //gnp_new needs to move
        double PD;
        Point_3D N;
        
        if (p_flag) {
            //hout<<"Penetration with GNP "<<GNP_i<<endl;
            
            //There is penetration, so then use EPA to find the penetration depth PD and direction vector N
            if (!GJK_EPA.EPA(gnps[GNP_i].vertices, gnp_new.vertices, simplex, N, PD)) {
                hout<<"Error in Move_gnps_if_needed when calling EPA"<<endl;
                return 0;
            }
            //hout<<"PD="<<PD<<" normal="<<N.str()<<endl;
            
            //Add to the vector of displacements
            if (!Add_to_vector_of_displacements(PD + cutoffs.van_der_Waals_dist, N, disps, disps_vec)) {
                hout<<"Error in Move_gnps_if_needed when calling Add_to_vector_of_displacements (1)"<<endl;
                return 0;
            }
        }
        else {
            //There is no penetration, so check wether they are touching or not
            
            if (t_flag) {
                //hout<<"Touch with GNP "<<GNP_i<<endl;
                
                //Find the simpleces that share the same plane/line/vertex and the direction in which
                //gnp_new should be moved
                if (!Find_direction_of_touching_gnps(GJK_EPA, gnps[GNP_i], gnp_new, N)) {
                    hout<<"Error in Move_gnps_if_needed when calling Find_direction_of_touching_gnps"<<endl;
                    return 0;
                }
                
                //Since the GNPs are touching, then the distance to be moved is just
                //the van der Waals distance
                //Add to the vector of displacements
                if (!Add_to_vector_of_displacements(cutoffs.van_der_Waals_dist, N, disps, disps_vec)) {
                    hout<<"Error in Move_gnps_if_needed when calling Add_to_vector_of_displacements (2)"<<endl;
                    return 0;
                }
            }
            else {
                //hout<<"No penetration with GNP "<<GNP_i<<endl;
                
                //Find the distance between the two GNPs and make sure they are separated at least the
                //van der Waals distance
                if (!GJK_EPA.Distance_and_direction_from_simplex_to_origin(simplex, N, PD)) {
                    hout<<"Error in Move_gnps_if_needed when calling Distance_and_direction_estimation"<<endl;
                    return 0;
                }
                
                //Check that the separation is at least the van der Waals distance
                double new_dist = cutoffs.van_der_Waals_dist - PD;
                if (new_dist > Zero) {
                    
                    //If the separation is less than the van der Waals distance,
                    //then gnp_new needs to be moved
                    //Add to the vector of displacements
                    if (!Add_to_vector_of_displacements(new_dist, N, disps, disps_vec)) {
                        hout<<"Error in Move_gnps_if_needed when calling Add_to_vector_of_displacements (3)"<<endl;
                        return 0;
                    }
                }
            }
        }
    }
    
    //Check if a displacement is needed
    if (disps.size()) {
        
        //Set the displace flag to true
        displaced = true;
        
        //hout<<"GNP moved disp_tot="<<disp_tot.str()<<endl;
        //A displacement is needed, then move gnp_new
        if (!Move_gnp_two_displacements(disps_vec, gnp_new)) {
            hout<<"Error Move_gnps_if_needed in when calling Move_gnp_two_displacements."<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function adds a new displacement to the vector of displacements and makes sure only the
//two largest displacements are kept
//Element with index 0 is the largest displacement
//It also keeps the displaments vector that correspond to the largest displacements
int Generate_Network::Add_to_vector_of_displacements(const double &disp, const Point_3D &N, vector<double> &disps, vector<Point_3D> &disps_vec)const
{
    //Check the size of the vector of diplacements
    if (disps.size() <= 1) {
        
        //Since the vectors have 1 or 0 elements, first just add the displacement
        disps.push_back(disp);
        disps_vec.push_back(N);
        
        //If the vector has only one element, terminate the function as there is nothing
        //more to check
        if (disps.size() == 1)
            return 1;
        
        //If not empty check if the elements need swapping
        if (disps[0] - disps[1] < Zero) {
            
            //Element in index 1 is larger, so swap the elements
            double tmp = disps[0];
            disps[0] = disps[1];
            disps[1] = tmp;
            
            //Swap the displacement vectors
            Point_3D tmp_p = disps_vec[0];
            disps_vec[0] = disps_vec[1];
            disps_vec[1] = tmp_p;
        }
    }
    else {
        
        //There are already two elements in the vectors, so check if the new displacement
        //is larger than any of the ones in the vector disps
        
        //Check if the new displacement is larger than the largest displacement
        if (disp > disps[0]) {
            
            //Copy element [0] into element [1]
            disps[1] = disps[0];
            disps_vec[1] = disps_vec[0];
            
            //Save new displament and vector into element [0]
            disps[0] = disp;
            disps_vec[0] = N;
        }
        //New displacement is not larger than the largest displacement
        //So check if the new displacement is larger than the second largest displacement
        else if (disp > disps[1]) {
            
            //Substitute element [1] by the new displacement and vector
            disps[1] = disp;
            disps_vec[1] = N;
        }
        //New displacement is not larger than the second largest displacement
        //So it is not added to the vectors
    }
    
    return 1;
}
//This function determines the direction in wich a GNP should be moved in case it is touching
//another GNP
int Generate_Network::Find_direction_of_touching_gnps(Collision_detection &GJK_EPA, const GNP &gnpA, const GNP &gnpB, Point_3D &N)const
{
    //Find a point V in the Minkowski sum
    Point_3D V = gnpA.center - gnpB.center;
    
    //Squared length of the vector from V to origin
    double len_V2 = V.length2();
    
    //Lambda function to obtain the distance from V to the vector OQ
    //The '&' lets me use variables in scope
    auto dist2 = [&](const Point_3D &Q_) {
        //Vector OQ = Q - 0 = Q
        //Vector VO = 0 - V = -V
        //Due to the dot product being squared in the return value,
        //it does not matter to use -V or V in this dot product
        double dot_ = V.dot(Q_);
        return (len_V2 - dot_*dot_/Q_.length2());
    };
    
    //Get the support vector in the direction from V to origin
    Point_3D Q = GJK_EPA.Support_AB(V*(-1), gnpA.vertices, gnpB.vertices);
    
    //Calculate the distance from V to the line segment OQ, this will be the reference distance
    double dist_ref;
    //Check if Q is the origin
    if (Q.length2() < Zero) {
        
        //Q is the origin, so the reference distance is just the distance from V to the origin
        dist_ref = len_V2;
    }
    else {
        
        //Q is not the origin, so use the lambda function to calculate the distance from V
        //to the line segment OQ
        dist_ref = dist2(Q);
    }
    
    
    //Sets to save the simplices in A and B that are touching
    set<int> simplexA, simplexB;
    
    //Iterate over all vertices of the Minkowski sum
    //i iterates over the vertices of A
    for (int i = 0; i < 8; i++) {
        //j iterates over the vertices of B
        for (int j = 0; j < 8; j++) {
            
            //Calculate the point in the Minkowski sum
            Q = gnpA.vertices[i] - gnpB.vertices[j];
            
            //Check if Q is the origin
            if (Q.length2() < Zero) {
                
                //Q is the origin, so i and j should be added to the simplices
                simplexA.insert(i);
                simplexB.insert(j);
            }
            
            //Q is not the origin, so check if the distance from V to line segment OQ
            //is the same (or almost the same) as the reference distance
            else if (abs(dist_ref - dist2(Q)) < Zero) {
                
                //The vertices in gnpA and gnpB are in the edge or face of the Minkowski sum
                //that crosses the origin, so add them to their corresponding simplices
                simplexA.insert(i);
                simplexB.insert(j);
            }
        }
    }
    
    //Check if there is a face in any of the simplices found
    if (simplexA.size() == 4 || simplexB.size() == 4) {
        
        //Find one face and the normal to that face
        if (simplexA.size() == 4) {
            
            //Get the first three vertices of simplex A
            set<int>::iterator it = simplexA.begin();
            int A0 = *it; ++it;
            int A1 = *it; ++it;
            int A2 = *it;
            
            //gnpA has a face as part of the touch, get the normal
            N = (gnpA.vertices[A1] - gnpA.vertices[A0]).cross(gnpA.vertices[A2] - gnpA.vertices[A0]);
        }
        else {
            
            //Get the first three vertices of simplex B
            set<int>::iterator it = simplexB.begin();
            int B0 = *it; ++it;
            int B1 = *it; ++it;
            int B2 = *it;
            
            //Only gnpB has a plane as part of the touch
            N = (gnpB.vertices[B1] - gnpB.vertices[B0]).cross(gnpB.vertices[B2] - gnpB.vertices[B0]);
        }
        
        //Make the normal a unit vector
        N.make_unit();
        
        //Make sure N goes in the direction from gnpA to gnpB
        if (N.dot(gnpB.center - gnpA.center) < Zero) {
            //Revert the direction
            N = N*(-1);
        }
    }
    else {
        
        //No face is touching, in such case gnpB can be moved in the direction from centerA to centerB
        N = (gnpB.center - gnpA.center).unit();
    }
    
    return 1;
}
int Generate_Network::Move_gnp_two_displacements(const vector<Point_3D> &disps_vec, GNP &gnp)const
{
    //Move the GNP according the the displacement in index [0]
    if (!Move_gnp(disps_vec[0], gnp)) {
        hout<<"Error in Move_gnp_two_displacements when calling Move_gnp."<<endl;
        return 0;
    }
    
    //Check if there is a second displacement
    if (disps_vec.size() == 2) {
        
        //Calculate the orthogonal component of disps_vec[1] with respect to disps_vec[0]
        
        //Calculate dot product
        double dot_p = disps_vec[0].dot(disps_vec[1]);
        
        //Calculate the squared length of disps_vec[0]
        double d0_2 = disps_vec[0].length2();
        
        //Calculate projection of disps_vec[1] onto disps_vec[0]
        Point_3D proj = disps_vec[0]*(dot_p/d0_2);
        
        //Calculate the orthoganal component we are looking for
        Point_3D orth = disps_vec[1] - proj;
        
        //Move GNP according the the displacement given by orth
        if (!Move_gnp(orth, gnp)) {
            hout<<"Error in Move_gnp_two_displacements when calling Move_gnp (orth)."<<endl;
            return 0;
        }
    }
    
    return 1;
}
//This function moves a GNP on a given direction and updates the plane equations of its faces
int Generate_Network::Move_gnp(const Point_3D &displacement, GNP &gnp)const
{
    //Move the GNP center
    gnp.center = gnp.center + displacement;
    
    //Move the vertices
    for (int i = 0; i < 8; i++) {
        gnp.vertices[i] = gnp.vertices[i] + displacement;
    }
    
    //Update the planes of the gnp faces
    if (!Update_gnp_plane_equations(gnp)) {
        hout<<"Error in Move_gnp when calling Update_gnp_plane_equations."<<endl;
        return 0;
    }
    
    return 1;
}
//This function adds a GNP to all the subregions it occupies (as given by the subregions vector)
int Generate_Network::Add_valid_gnp_to_subregions(const int &gnp_new_idx, const set<int> &subregions, vector<vector<int> > &sectioned_domain)const
{
    //Also, since gnp_new was successfully generated, add it to all corresponding sub regions
    for (set<int>::iterator i = subregions.begin(); i != subregions.end(); i++) {
        
        //Current subregion
        int s = *i;
        
        //Add GNP number to current subregion
        sectioned_domain[s].push_back(gnp_new_idx);
    }
    
    return 1;
}
//This function checks if GNP is completely inside the sample
int Generate_Network::Check_all_gnp_vertices_are_inside_sample(const Geom_sample &sample_geom, GNP &gnp, bool &is_enough_inside)const
{
    //Iterate over the eight vertices of the GNP and check if all of them are inside the sample
    for (int i = 0; i < 8; i++) {
        
        //Check if current vertex is outside the sample
        if (!Is_point_inside_cuboid(sample_geom.sample, gnp.vertices[i])) {
            
            //A vertex was found to be outside the sample
            //Break the loop as it is not necessary to continue cheking the rest of the vertices
            //Flag is_enough_inside is already false, so there is no need to update it
            //Thus, just terminate the function
            return 1;
        }
    }
    
    //The eight GNP vertices were inside the sample, so set the flag is_enough_inside to true
    is_enough_inside = true;
    
    return 1;
}
//This function calculates the generated volume of a GNP and adds it to the global variables
int Generate_Network::Calculate_generated_gnp_vol(const GNP_Geo gnp_geom, const Geom_sample &sample_geom, GNP &gnp, bool &is_enough_inside)const
{
    //---------------------------------------------------------------------------
    //Add the volume and weight corresponding to the GNP
    double gnp_vol = gnp.volume;
    
    //Get the location of the first vertex
    int loc_0 = Is_point_inside_cuboid(sample_geom.sample, gnp.vertices[0]);
    
    //Iterate over the eight vertices of the GNP and check if all of them have the same
    //location as the first vertex
    int i = 0;
    for (i = 1; i < 8; i++) {
        
        //Check if current vertex has the same location as the first vertex
        if (loc_0 != Is_point_inside_cuboid(sample_geom.sample, gnp.vertices[i])) {
            
            //Vertex i has a different location than vertex 0
            //This means that the GNP is partially inside the sample
            //Thus, approximate the GNP volume inside the sample
            if (!Approximate_gnp_volume_inside_sample(sample_geom.sample, gnp, gnp_vol)) {
                hout << "Error in Calculate_generated_gnp_vol when calling Approximate_gnp_volume_inside_sample." << endl;
                return 0;
            }
            
            //Check the criterion for the minimum GNP volume inside the sample
            if (gnp_geom.vol_in == "no_min" || gnp_vol/gnp.volume - gnp_geom.min_vol_in >= Zero) {
                
                //If no minimum specified there is at least some volume inside the sample
                //If specified, the minimum required GNP volume is acutually inside the sample
                
                //Thus, update the flag that indicates that the minimun required GNP volume
                //is actually inside the sample
                is_enough_inside = true;
                
                //Update the GNP's volume
                gnp.volume = gnp_vol;
            }
            
            //Terminate the function as it is not necessary to continue
            //cheking the rest of the vertices
            return 1;
        }
    }
    
    //If all vertices had the same location as the first one, this means that the GNP is
    //either completely inside or completely outside the sample
    //Check if the GNP is completely inside the sample
    if (loc_0) {
        
        //Update the flag that indicates that the minimun required GNP volume is acutually
        //inside the sample
        //No need to update GNP volume as the variable already has the volume of the GNP
        is_enough_inside = true;
    }
    
    return 1;
}
//This function approximates the volume of a GNP inside the sample depending of the number of points in
//the discretization of its middle plane that are inside the sample
int Generate_Network::Approximate_gnp_volume_inside_sample(const cuboid &sample_geom, const GNP &gnp, double &gnp_vol)const
{
    //---------------------------------------------------------------------------
    //Number of points per side to approximate volume
    int n_points = 20;
    
    //Variable to count the number of points inside the sample
    int points_in = 0;
    
    //Midpoint between vertices 3 and 7
    Point_3D corner = (gnp.vertices[3] + gnp.vertices[7])/2.0;
    //Midpoint between vertices 0 and 4
    Point_3D corner_right = (gnp.vertices[0] + gnp.vertices[4])/2.0;
    //Midpoint between vertices 2 and 6
    Point_3D corner_up = (gnp.vertices[2] + gnp.vertices[6])/2.0;
    //Midpoint between vertices 1 and 5
    Point_3D corner_up_right = (gnp.vertices[1] + gnp.vertices[5])/2.0;
    
    //Direction vector from corner to corner_right
    Point_3D dir_right = (corner_right - corner);
    //Direction vector from corner_up to corner_up_right
    Point_3D dir_right_up = (corner_up_right - corner_up);
    
    //Iterate over the total number of points
    //Index i moves new_point to the right
    for (int i = 0; i < n_points; i++) {
        
        //Get the starting point in the bottom edge
        double lambda_right = (double)i/(n_points-1);
        Point_3D start = corner + dir_right*lambda_right;
        
        //Get the ending point in the top edge
        Point_3D end = corner_up + dir_right_up*lambda_right;
        
        //Get the direction up for the new point
        Point_3D dir_up_new = end - start;
        
        //Index j moves new_point up
        for (int j = 0; j < n_points; j++) {
            
            //Get the value of the point that will move over the surface
            double lambda_up = (double)j/(n_points-1);
            Point_3D new_point = start + dir_up_new*lambda_up;
            
            //Check if new point is inside the sample
            if (Is_point_inside_cuboid(sample_geom, new_point)) {
                //The point is inside the sample, so increase the number of points inside the sample
                points_in++;
            }
        }
    }
    
    //Check if all points where outside the sample
    if (points_in != 0) {
        
        //The gnp is partially outside the sample, so approximate the volume inside
        gnp_vol = gnp_vol*((double)points_in/(n_points*n_points));
    }
    
    
    
    return 1;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function generates CNTs for a mixed nanoparticle network
int Generate_Network::Generate_cnt_network_threads_among_gnps_mt(const Simu_para &simu_para, const GNP_Geo &gnp_geo, const Nanotube_Geo &nanotube_geo, const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const double &gnp_vol_tot, const double &gnp_wt_tot, vector<vector<Point_3D> > &cnts_points, vector<double> &cnts_radius)const
{
    //Initial seeds, if any are in network_seeds within geom_sample.
    //However, geom_sample cannot be modified, so copy the seeds to a new vector
    unsigned int net_seeds[7];
    if (!CNT_seeds(simu_para.CNT_seeds, net_seeds)) {
        hout<<"Error in CNT_seeds"<<endl;
        return 0;
    }
    
    //Use the seeds generated above
    std::mt19937 engine_x(net_seeds[0]);
    std::mt19937 engine_y(net_seeds[1]);
    std::mt19937 engine_z(net_seeds[2]);
    std::mt19937 engine_phi(net_seeds[3]);
    std::mt19937 engine_theta(net_seeds[4]);
    std::mt19937 engine_rand(net_seeds[5]);
    std::mt19937 engine_initial_direction(net_seeds[6]);
    
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [0, 1].
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    // "Filter" MT's output to generate double values, uniformly distributed on the closed interval [-1, 1].
    std::uniform_real_distribution<double> dist_initial(-1.0, 1.0);

    //---------------------------------------------------------------------------
    //Time variables to keep track of generation time
    time_t ct0, ct1;
    
    //---------------------------------------------------------------------------
    //Variables for the generated CNT volume and weight
    double vol_sum = 0;
    
    //Variable to count the generated CNT seeds
    int cnt_seed_count = 0;
    
    //Variable to coun the number of rejected CNTs
    int cnt_reject_count = 0;
    //Variable to count the number of ignored CNTs
    int cnt_ignore_count = 0;
    
    //Variable to count the number of times that a point had to be relocated
    int point_overlap_count = 0;
    //Variable to count the number of points that were overlapping other points
    int point_overlap_count_unique = 0;
    
    //---------------------------------------------------------------------------
    //Vectors for handling CNT penetration
    //global_coordinates[i][0] stores the CNT number of global point i
    //global_coordinates[i][1] stores the local point number of global point i
    vector<vector<int> > global_coordinates;
    //sectioned_domain[i] contains all the points in sub-region i.
    //Sub-region i is an overlapping subregion to check for penetrations
    vector<vector<long int> > sectioned_domain;
    //n_subregions[0] is the number of subregions along x
    //n_subregions[1] is the number of subregions along y
    //n_subregions[2] is the number of subregions along z
    int n_subregions[3];
    //Initialize the vector sub-regions
    //hout<<"Initialize_cnt_subregions"<<endl;
    Initialize_cnt_subregions(geom_sample, n_subregions, sectioned_domain);
    
    //Get the time when generation started
    ct0 = time(NULL);
    //Check when 10% is completed
    double vol_completed = 0.1;
    //Variable used to store the fraction of CNT volume indicated by vol_completed
    double vol_completed_acc = vol_completed*nanotube_geo.volume;
    
    //Boolean to terminate the main while loop, initialized to false to start the loop
    bool terminate = false;
    
    //---------------------------------------------------------------------------
    while(!terminate)
    {
        //hout<<"CNTs generated:"<<cnts_points.size()<<endl;
        //---------------------------------------------------------------------------
        //Vector for a new nanotube
        vector<Point_3D> new_cnt;
        
        //---------------------------------------------------------------------------
        //Randomly generate a CNT length and CNT radius
        double cnt_length, cnt_rad;
        //hout<<"Get_random_value_mt 1"<<endl;
        if (!Get_length_and_radius(nanotube_geo, engine_rand, dist, cnt_length, cnt_rad)) {
            hout << "Error in Generate_network_threads_mt when calling Get_length_and_radius" <<endl;
            return 0;
        }
        
        //Calculate the number of CNT growth steps
        int step_num = (int)(cnt_length/nanotube_geo.step_length) + 1;
        
        //---------------------------------------------------------------------------
        //Cross-sectional area of the current CNT. It is used to calculate the CNT volume.
        const double cnt_cross_area = PI*cnt_rad*cnt_rad;
        
        //---------------------------------------------------------------------------
        //Randomly generate an initial direction, then generate the rotation matrix that results in that rotation
        MathMatrix multiplier(3,3);
        //hout<<"Get_initial_direction_mt"<<endl;
        if (!Get_initial_direction_mt(nanotube_geo.dir_distrib_type, nanotube_geo.ini_theta, nanotube_geo.ini_phi, engine_initial_direction, dist_initial, multiplier)) {
            hout << "Error in Generate_network_threads_mt when calling Get_initial_direction_mt" <<endl;
            return 0;
        }
        
        //Calculate a cutoff for penetration
        double rad_p_dvdw = cnt_rad + cutoffs.van_der_Waals_dist;
        //Calculate some cutoffs for self-penetration
        double cnt_cutoff = cnt_rad + rad_p_dvdw;
        double cnt_cutoff2 = cnt_cutoff*cnt_cutoff;
        
        //Map to find self-penetrating points
        map<int, vector<int> > subr_point_map;
        
        //---------------------------------------------------------------------------
        //Generate an initial point of the CNT in the extended domain
        //hout<<"Generate_initial_point_of_cnt"<<endl;
        Point_3D new_point;
        if (Generate_initial_point_of_cnt(geom_sample, simu_para, cnts_points, cnts_radius, new_cnt, rad_p_dvdw, cnt_cutoff, cnt_cutoff2, nanotube_geo.step_length, subr_point_map, global_coordinates, sectioned_domain, gnps, sectioned_domain_gnp, n_subregions, new_point, engine_x, engine_y, engine_z, dist) != 1) {
            hout << "Error in Generate_network_threads_mt when calling Generate_initial_point_of_cnt" <<endl;
            return 0;
        }
        
        //Add the CNT seed to the current CNT vector
        new_cnt.push_back(new_point);
        
        //---------------------------------------------------------------------------
        //Increase the count of initial CNT points generated
        cnt_seed_count++;
        //hout << "Seed="<<cnt_seed_count<<endl;
        int max_seed = 1E9;
        if(cnt_seed_count>max_seed) {
            hout << "The number of generated seeds is lager than "<<max_seed<<", but the nanotube generation still fails to acheive the requested volume fraction." << endl;
            return 0;
        }
        
        //Get the location of the initial point
        bool is_prev_in_sample = Is_point_inside_cuboid(geom_sample.sample, new_cnt[0]);
        
        //Set the flag of the point depending on its position
        //1: inside the sample
        //0: outside the sample
        new_cnt[0].flag = (is_prev_in_sample)? 1: 0;
        
        //Variable to count the number of points of the new CNT that are inside the sample
        //It is initialize with the same value as the flag
        //int points_in = (is_prev_in_sample)? 1: 0;
        int points_in = new_cnt[0].flag;
        
        //Variable to store the length of the current CNT that is inside the sample
        double cnt_len = 0.0;
        
        //Add the seed point to the overlapping subregions it belongs to using the map
        if (!Add_cnt_point_to_overlapping_regions_map(geom_sample, new_cnt[0], 0, is_prev_in_sample, n_subregions, subr_point_map)) {
            hout<<"Error when adding the seed point to the map of subregions"<<endl;
            return 0;
        }
        
        //Start generation of a CNT siven the generated seed
        for(int i = 0; i < step_num; i++)
        {
            //Generates a new CNT point such that the new CNT segment has a random orientation
            if (!Get_direction_and_point(nanotube_geo, multiplier, new_point, engine_theta, engine_phi, dist)) {
                hout<<"Error in generating a new CNT point (not a seed) at iteration i="<<i<<endl;
                return 0;
            }
            
            //Check if the new point, cnt_poi, is inside the extended domain
            if(Is_point_inside_cuboid(geom_sample.ex_dom_cnt, new_point))
            {
                //If the point is inside the extended domain, then check for penetration if non-penetrating
                //model is used, or just add it if penetrating model is used
                
                //Check for penetration
                //hout << "Check penetration "<<endl;
                int mixed_interpenetration = 1;
                //if (penetration_model_flag) {
                if (simu_para.penetration_model_flag) {
                    
                    //Get the status on mixed penetration
                    //hout<<"Check_mixed_interpenetration"<<endl;
                    mixed_interpenetration = Check_mixed_interpenetration(geom_sample, cnts_points, cnts_radius, new_cnt, rad_p_dvdw, cnt_cutoff, cnt_cutoff2, nanotube_geo.step_length, subr_point_map, global_coordinates, sectioned_domain, gnps, sectioned_domain_gnp, n_subregions, new_point);
                    if (mixed_interpenetration == -1) {
                        hout << "Error in Generate_network_threads_mt when calling Check_mixed_interpenetration (CNT initial point)" <<endl;
                        return 0;
                    }
                    
                }
                
                if (!simu_para.penetration_model_flag || mixed_interpenetration) {
                //if (!penetration_model_flag || mixed_interpenetration) {
                    
                    //---------------------------------------------------------------------------
                    //If the penetrating model is used or if new_point was placed
                    //in a valid position (i.e., without penetrating another CNT or GNP)
                    //then add the point to the generated CNTs
                    
                    //Calculate the segment length inside the sample and add it to the total CNT length
                    bool is_new_inside_sample;
                    //hout<<"Length_inside_sample"<<endl;
                    cnt_len = cnt_len + Length_inside_sample(geom_sample.sample, new_cnt.back(), new_point, is_prev_in_sample, is_new_inside_sample);
                    
                    //If new_point is inside the sample, then increase the number
                    //of points inside the sample of the new CNT
                    if (is_new_inside_sample) {
                        points_in++;
                        
                        //Set the flag of the point equal to 1 to indicate the point is inside the sample
                        new_point.flag = 1;
                    }
                    else {
                        //Set the flag of the point equal to 0 to indicate the point is in the boundary layer
                        new_point.flag = 0;
                    }
                    
                    //For the next iteration of the for loop, cnt_poi will become previous point,
                    //so update the boolean is_prev_in_sample for the next iteration
                    is_prev_in_sample = is_new_inside_sample;
                    
                    //Add the new_point to the overlapping subregions it belongs to using the map
                    if (!Add_cnt_point_to_overlapping_regions_map(geom_sample, new_point, (int)new_cnt.size(), is_new_inside_sample, n_subregions, subr_point_map)) {
                        hout<<"Error when adding a point to the map of subregions"<<endl;
                        return 0;
                    }
                    
                    //Add the new point to the current CNT
                    new_cnt.push_back(new_point);
                    
                    //---------------------------------------------------------------------------
                    //Check if the target volume fraction has been reached
                    if( (vol_sum + cnt_len*cnt_cross_area) >= nanotube_geo.volume) {
                        
                        //Set the terminate variable to true so that the main while-loop is terminated
                        terminate = true;
                        
                        //Break the for-loop so that the current CNT stops its growth
                        break;
                    }
                    
                } else {
                    //---------------------------------------------------------------------------
                    //If the penetrating point could not be accommodated, then delete the current CNT
                    //so that a new one is generated
                    //hout << "Penetrating point could not be accommodated" <<endl;
                    
                    //Clear the new_cnt vector so that it is not added to the rest of CNTs
                    new_cnt.clear();
                    
                    //Increase the count of rejected cnts
                    cnt_reject_count++;
                    //hout <<"cnt_reject_count="<<cnt_reject_count<< endl;
                    
                    //Break the for-loop to terminate the CNT growth
                    break;
                }
                //hout << "done" << endl;
            }
            else {
                
                //If the point is outside the extended domain, break the for-loop
                //to terminate the CNT growth
                break;
            }
            //for-loop ends here
        }
        
        //---------------------------------------------------------------------------
        //Store or ignore the CNT points
        //hout<<"Store_or_ignore_new_cnt"<<endl;
        /*if (!Store_or_ignore_new_cnt(geom_sample, simu_para.penetration_model_flag, points_in, cnt_len, cnt_rad, cnt_cross_area, new_cnt, cnts_points, cnts_radius, n_subregions, sectioned_domain, global_coordinates, vol_sum, cnt_ignore_count)) {
            hout<<"Error when storing or ignoring a new CNT"<<endl;
            return 0;
        }*/
        if (!Store_or_ignore_new_cnt_using_map(simu_para.penetration_model_flag, points_in, cnt_len, cnt_rad, cnt_cross_area, new_cnt, cnts_points, cnts_radius, subr_point_map, sectioned_domain, global_coordinates, vol_sum, cnt_ignore_count)) {
            hout<<"Error when storing or ignoring a new CNT"<<endl;
            return 0;
        }
        
        //Get the time to check progress
        ct1 = time(NULL);
        
        //Check progress
        if (!Check_progress("CNT", (int)(ct1-ct0), nanotube_geo.volume, vol_sum, vol_completed, vol_completed_acc)) {
            hout<<"Error when calculating the percentage of CNT volume generated"<<endl;
            return 0;
        }
        
        //while-loop ends here
    }
    
    //Output the CNT content generated
    if(nanotube_geo.criterion == "wt") {
        
        //Calculate matrix weight
        double matrix_weight = (geom_sample.volume - vol_sum)*geom_sample.matrix_density;
        
        //Calculate the CNT weight
        double cnt_weight = vol_sum*nanotube_geo.density;
        
        hout << endl << "The weight fraction of generated CNTs is: " << cnt_weight/(matrix_weight + cnt_weight);
        hout << ", the target weight fraction was " << nanotube_geo.weight_fraction << endl << endl;

    } else if(nanotube_geo.criterion == "vol") {
        hout << endl << "The volume fraction of generated CNTs was " << vol_sum/geom_sample.volume;
        hout << ", the target volume fraction was " << nanotube_geo.volume_fraction << endl << endl;
    }
    
    hout << "There were " << point_overlap_count_unique << " overlapping points and ";
    hout << point_overlap_count << " overlaps, " << cnt_reject_count << " CNTs were rejected and "<<cnt_ignore_count<<" were ignored." << endl;
    
    
    return 1;
}
//This function generates the initial point of a CNT
//If non-penetrating model is used, then it also makes sure that the initial point of a CNT
//is in a valid position (i.e., without interpenetrations)
int Generate_Network::Generate_initial_point_of_cnt(const Geom_sample &geom_sample, const Simu_para &simu_para, const vector<vector<Point_3D> > &cnts, const vector<double> &radii, vector<Point_3D> &new_cnt, const double &rad_plus_dvdw, const double &cnt_cutoff, const double &cnt_cutoff2, const double &step_length, const map<int, vector<int> > &subr_point_map, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const int n_subregions[], Point_3D &new_point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const
{
    //Variable to count the attempts
    int attempts = 0;
    
    while (attempts <= MAX_ATTEMPTS) {
        
        //Genereta a new seed point
        //hout<<"attempts seed="<<attempts<<endl;
        if(!Get_point_in_cuboid_mt(geom_sample.ex_dom_cnt, new_point, engine_x, engine_y, engine_z, dist)){
            hout << "Error in Generate_network_threads_mt when calling Get_point_in_cuboid_mt (CNT initial point)" <<endl;
            return -1;
        }
        
        if (simu_para.penetration_model_flag) {
            
            //If non-penetrating model is used, then make sure the intial point of the CNT
            //is in a valid position
            //hout<<"Check_mixed_interpenetration"<<endl;
            int interpenetration = Check_mixed_interpenetration(geom_sample, cnts, radii, new_cnt, rad_plus_dvdw, cnt_cutoff, cnt_cutoff2, step_length, subr_point_map, global_coordinates, sectioned_domain_cnt, gnps, sectioned_domain_gnp, n_subregions, new_point);
            if (interpenetration == -1) {
                hout << "Error in Generate_network_threads_mt when calling Check_mixed_interpenetration (CNT initial point)" <<endl;
                return -1;
            }
            else if (interpenetration) {
                
                //There is no interpenetration, so loop can be terminated
                return 1;
            }
        }
        else {
            
            //There is no interpenetration, so loop can be terminated
            return 1;
        }
        
        //Increase the number of attempts
        attempts++;
    }
    
    //There is still penetration so terminate with 0 and output error message
    hout << "Error: Too many attempts to resolve overlapping of an intial CNT point . ";
    hout << new_point.x << ' ' << new_point.y << ' ' << new_point.z << endl;
    return 0;
}
//This function checks if a CNT point is penetrating a CNT or GNP
int Generate_Network::Check_mixed_interpenetration(const Geom_sample &geom_sample, const vector<vector<Point_3D> > &cnts, const vector<double> &radii, vector<Point_3D> &new_cnt, const double &rad_plus_dvdw, const double &cnt_cutoff, const double &cnt_cutoff2, const double &step_length, const map<int, vector<int> > &subr_point_map, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const int n_subregions[], Point_3D &new_point)const
{
    //Vector that stores the coordintes of the points that the input point "new_point" is penetrating
    vector<Point_3D> affected_points;
    //Vector that stores the distance at which the two points should be
    vector<double> cutoffs_p;
    //Vector that stores the distance at which the two points actually are
    vector<double> distances;
    
    //Variable to count the number of attempts
    int attempts = 0;
    
    //Loop while the maximum numbe of attempts has not been reached
    while (attempts <= MAX_ATTEMPTS) {
        
        //Get the subregion of CNT p, GNP &affected_gnpoint
        //hout<<"attempts="<<attempts<<endl;
        int subregion = Get_cnt_point_subregion(geom_sample, n_subregions, new_point);
        //hout<<"new_point="<<new_point.str()<<" subregion="<<subregion<<endl;
        if (subregion == -1) {
            //If the sub-region is -1, then the point is in te boundary layer,
            //so there is no need to check penetration
            return 1;
        }
        
        //Check if point is too close or penetrating a GNP
        //hout<<"Get_gnp_penetrating_points"<<endl;
        GNP affected_gnp;
        if (!Get_gnp_penetrating_points(gnps, sectioned_domain_gnp[subregion], rad_plus_dvdw, new_point, affected_points, cutoffs_p, distances, affected_gnp)) {
            hout<<"Error in Check_mixed_interpenetration when calling Get_gnp_penetrating_points"<<endl;
            return -1;
        }
        //hout<<"sectioned_domain_gnp[subregion].size="<<sectioned_domain_gnp[subregion].size()<<endl;
        
        //Check for a special case
        if (!affected_points.empty() && affected_points[0].flag == -1) {
            //This is a special case when the point is inside a GNP
            //Since the point is inside a GNP, there is no need to check other points
            //because there cannot be penetration with other CNTs, this would have
            //been avoided due to using the non-penetrating model
            
            //Go to the special case when the CNT point is inside the GNP
            //hout<<"Deal_with_point_inside_gnp"<<endl;
            if (!Deal_with_point_inside_gnp(geom_sample, affected_gnp, rad_plus_dvdw, new_cnt, new_point)) {
                //new_point could not be accomodated so reject it
                return 0;
            }
            
            //Check if the segment has a valid orientation, only when:
            //new_cnt has at least two points
            //AND
            //new_cnt.back() is inside the sample
            //hout<<"Is_point_inside_cuboid"<<endl;
            if (new_cnt.size() >= 2 && Is_point_inside_cuboid(geom_sample.sample, new_cnt.back())) {
                
                //hout<<"Check_segment_orientation 1"<<endl;
                if (!Check_segment_orientation(new_point, new_cnt)) {
                    //When not in a valid position it cannot be moved again so a new CNT is needed
                    //hout<<"Point not in a valid position (special case)"<<endl;
                    return 0;
                }
            }
            
        }
        else {
            
            //Check if point is too close or penetrating another CNT
            //hout<<"Get_penetrating_points"<<endl;
            Get_penetrating_points(cnts, global_coordinates, sectioned_domain_cnt[subregion], radii, rad_plus_dvdw, new_point, affected_points, cutoffs_p, distances);
            
            //Check if there are any penetration within the CNT
            //hout<<"Get_penetrating_points_within_cnt"<<endl;
            Get_penetrating_points_within_cnt(subregion, cnt_cutoff, cnt_cutoff2, new_point, new_cnt, subr_point_map, affected_points, cutoffs_p, distances);
            //hout<<"affected_points.empty()="<<affected_points.empty()<<endl;
            
            //Check if there are any penetrating points
            if (!affected_points.empty()) {
                
                //To determine how to move the point, check if it is a seed point
                if (new_cnt.empty()) {
                    
                    //Point is a seed, so then use functions that move the point
                    //hout<<"Move_point"<<endl;
                    Move_point(cutoffs_p, distances, affected_points, new_point);
                }
                else {
                    
                    //Point is not a seed so use functions that rotate the CNT segment
                    //hout<<"Move_point_by_totating_cnt_segment cnt.back="<<((new_cnt.empty())?" ":new_cnt.back().str())<<endl;
                    //hout<<"new_point="<<new_point.str()<<endl;
                    if (!Move_point_by_totating_cnt_segment(step_length, new_cnt, cutoffs_p, distances, affected_points, new_point)) {
                        hout<<"Error in Check_mixed_interpenetration when calling Move_point_by_totating_cnt_segment"<<endl;
                        return -1;
                    }
                    //hout<<"new_point="<<new_point.str()<<endl;
                }
                
                //Check if the segment has a valid orientation
                //hout<<"Check_segment_orientation 2"<<endl;
                if (!Check_segment_orientation(new_point, new_cnt)) {
                    //When not in a valid position it cannot be moved again so a new CNT is needed
                    //hout<<"Point not in a valid position"<<endl;
                    return 0;
                }
            }
            else {
                
                //The point does not need to be moved so terminate with true
                return 1;
            }
        }
        
        //Clear the vectors affected_points, contact_coordinates and temporal_contacts
        //so they are used in the next iteration
        affected_points.clear();
        cutoffs_p.clear();
        distances.clear();
        
        //Increase the number of attempts
        attempts++;
    }
    
    //The point could not be accomodated so terminate with false
    //hout<<"Point could not be accomodated"<<endl;
    return 0;
}
//This function finsd the GNPs that are too close to a CNT point penetrating the GNP
//The points in the GNP that are in contact with the CNT are also found
int Generate_Network::Get_gnp_penetrating_points(const vector<GNP> &gnps, const vector<int> &subregion_gnp, const double &cutoff, const Point_3D &new_point, vector<Point_3D> &affected_points, vector<double> &cutoffs_p, vector<double> &distances, GNP &affected_gnp)const
{
    
    //Scan all GNPs in the subregion
    for (int i = 0; i < subregion_gnp.size(); i++) {
        
        //Get the GNP number of the current subregion
        int GNP_i = subregion_gnp[i];
        
        //Find the distance from new_point to GNP  and the closest point in the GNP
        double distance = 0;
        //hout<<"Get_gnp_point_closest_to_point i="<<i<<endl;
        Point_3D P = Get_gnp_point_closest_to_point(gnps[GNP_i], new_point, distance);
        
        //Check if the point is inside the GNP
        if (P.flag == -1) {
            
            //The point is inside of the GNP, so set this as the only point in affected points
            affected_points.clear();
            affected_points.push_back(P);
            
            //Get the affected GNP
            affected_gnp = gnps[GNP_i];
            
            //Terminate the function as this case is handled on its own
            return 1;
        }
        
        //Check if distance is below the cutoff
        if (distance < cutoff) {
            
            //Add the current point, distance and cutoff to the output vectors
            affected_points.push_back(P);
            cutoffs_p.push_back(cutoff);
            distances.push_back(distance);
        }
    }
    
    return 1;
}
//This function finds the point on a GNP closest to a given point
Point_3D Generate_Network::Get_gnp_point_closest_to_point(const GNP &gnp, const Point_3D &P, double &distance) const
{
    
    //These vectors are used more than once
    Point_3D V0P = P - gnp.vertices[0];
    Point_3D V1P = P - gnp.vertices[1];
    
    //Lambda function to calculate projection on a plane
    auto calculate_pojection = [](const Point_3D &P, const Point_3D &N, const double &dist) {
        //Calcualte point in plane
        Point_3D Q = P - N*dist;
        //Set flag to 0
        Q.flag = 0;
        return Q;
    };
    
    //Check if P is above plane 0
    if (gnp.faces[0].N.dot(V0P) > Zero) {
        
        //Call the function to find the point P and distance to closest simplex
        //above the plane of face 0
        //hout<<"P is above face 0"<<endl;
        return Get_point_closest_to_large_gnp_face(gnp, 0, 1, 2, 3, 0, P, distance);
    }
    else {
        
        //Check if P is above plane 1
        if (gnp.faces[1].N.dot(P - gnp.vertices[4]) > Zero) {
            
            //Call the function to find the point P and distance to closest simplex
            //above the plane of face 1
            //hout<<"P is above face 1"<<endl;
            return Get_point_closest_to_large_gnp_face(gnp, 4, 5, 6, 7, 1, P, distance);
        }
        else {
            
            //Check if point P is above face 2
            if (gnp.faces[2].N.dot(V0P) > Zero) {

                //hout<<"P is above face 2"<<endl;
                //Face V0V3V7V4 (2), edge V0V4, or edge V3V7 is closest
                if (gnp.faces[3].N.dot(V0P) > Zero) {
                    
                    //Edge V0V4 is closest
                    //hout<<"Edge V0V4 is closest"<<endl;
                    return Distance_from_point_to_edge(P, gnp.vertices[0], gnp.vertices[4], distance);
                }
                else {
                    if (gnp.faces[5].N.dot(P - gnp.vertices[3]) > Zero) {
                        
                        //Edge V3V7 is closest
                        //hout<<"Edge V3V7 is closest"<<endl;
                        return Distance_from_point_to_edge(P, gnp.vertices[3], gnp.vertices[7], distance);
                    }
                    else {
                        
                        //Face V0V3V7V4 (2) is closest
                        //hout<<"Face V0V3V7V4 (2) is closest"<<endl;
                        distance = gnp.faces[2].distance_to(P);

                        //Calcualte point in plane
                        return calculate_pojection(P, gnp.faces[2].N, distance);
                    }
                }
            }
            else {
                
                //Check if point P is above face 4
                if (gnp.faces[4].N.dot(V1P) > Zero) {
                    
                    //Face V1V2V6V5 (4), edge V1V5, or edge V2V6 is closest
                    //hout<<"point P is above face 4"<<endl;
                    if (gnp.faces[3].N.dot(V1P) > Zero) {
                        
                        //Edge V1V5 is closest
                        //hout<<"Edge V1V5 is closest"<<endl;
                        return Distance_from_point_to_edge(P, gnp.vertices[1], gnp.vertices[5], distance);
                    }
                    else {
                        if (gnp.faces[5].N.dot(P - gnp.vertices[2]) > Zero) {
                            
                            //Edge V2V6 is closest
                            //hout<<"Edge V2V6 is closest"<<endl;
                            return Distance_from_point_to_edge(P, gnp.vertices[2], gnp.vertices[6], distance);
                        }
                        else {
                            
                            //Face V1V2V6V5 (4) is closest
                            //hout<<"Face V1V2V6V5 (4) is closest"<<endl;
                            distance = gnp.faces[4].distance_to(P);
                            
                            //Calcualte point in plane
                            return calculate_pojection(P, gnp.faces[4].N, distance);
                        }
                    }
                }
                else {
                    
                    //Face V0V1V5V4 or face V2V3V7V6 is closest or P is inside
                    if (gnp.faces[3].N.dot(V0P) > Zero) {
                        
                        //Face V0V1V5V4 (3) is closest
                        //hout<<"Face V0V1V5V4 (3) is closest"<<endl;
                        distance = gnp.faces[3].distance_to(P);
                        
                        //Calcualte point in plane
                        return calculate_pojection(P, gnp.faces[3].N, distance);
                    }
                    else {
                        if (gnp.faces[5].N.dot(P - gnp.vertices[2]) > Zero) {
                            
                            //Face V2V3V7V6 (5) is closest
                            //hout<<"Face V2V3V7V6 (5) is closest"<<endl;
                            distance = gnp.faces[5].distance_to(P);
                            
                            //Calcualte point in plane
                            return calculate_pojection(P, gnp.faces[5].N, distance);
                        }
                        else {
                            
                            //Point is inside GNP
                            //hout<<"Point is inside GNP"<<endl;
                            //Create a point with flag equal to -1
                            Point_3D out;
                            out.flag = -1;
                            
                            //Set the distance equal to zero so that the point is identified as
                            //being below the cutoff
                            distance = 0.0;
                            
                            return out;
                        }
                    }
                }
            }
        }
    }
}
Point_3D Generate_Network::Get_point_closest_to_large_gnp_face(const GNP &gnp, const int &V0, const int &V1, const int &V2, const int &V3, const int &F, const Point_3D &P, double &distance)const
{
    //These vectors are used more than once
    Point_3D V0P = P - gnp.vertices[V0];
    Point_3D V1P = P - gnp.vertices[V1];
    
    //P is above plane 0
    //Now check if it is outside on the side of face 2
    if (gnp.faces[2].N.dot(V0P) > Zero) {
        
        //P might be closest to V0V3, V0 or V3
        //hout<<"outside on the side of face 2"<<endl;
        if (gnp.faces[3].N.dot(V0P) > Zero) {
            
            //P is closest to V0
            //Calcualte the distance between points
            //hout<<"P is closest to V0"<<endl;
            distance = P.distance_to(gnp.vertices[V0]);
            
            //Create a point with flag zero
            Point_3D Q = gnp.vertices[V0];
            Q.flag = 0;
            
            return Q;
        }
        else {
            if (gnp.faces[5].N.dot(P - gnp.vertices[V3]) > Zero) {
                
                //P is closest to V3
                //Calcualte the distance between points
                //hout<<"P is closest to V3"<<endl;
                distance = P.distance_to(gnp.vertices[V3]);
                
                //Create a point with flag zero
                Point_3D Q = gnp.vertices[V3];
                Q.flag = 0;
                
                return Q;
            }
            else {
                
                //P is closest to V0V3
                //hout<<"P is closest to V0V3"<<endl;
                return Distance_from_point_to_edge(P, gnp.vertices[V0], gnp.vertices[V3], distance);
            }
        }
    }
    else {

        //P is above plane of face 0
        //Now check if it is outside on the side of face 4
        if (gnp.faces[4].N.dot(V1P) > Zero) {
            
            //P might be closest to V1V2, V1 or V2
            //hout<<"P is outside on the side of face 4"<<endl;
            if (gnp.faces[3].N.dot(V1P) > Zero) {
                
                //P is closest to V1
                //Calcualte the distance between points
                //hout<<"P is closest to V1"<<endl;
                distance = P.distance_to(gnp.vertices[V1]);
                
                //Create a point with flag zero
                Point_3D Q = gnp.vertices[V1];
                Q.flag = 0;
                
                return Q;
            }
            else {
                if (gnp.faces[5].N.dot(P - gnp.vertices[V2]) > Zero) {
                    
                    //P is closest to V2
                    //Calcualte the distance between points
                    //hout<<"P is closest to V2"<<endl;
                    distance = P.distance_to(gnp.vertices[V2]);
                    
                    //Create a point with flag zero
                    Point_3D Q = gnp.vertices[V2];
                    Q.flag = 0;
                    
                    return Q;
                }
                else {
                    
                    //P is closest to V1V2
                    //hout<<"P is closest to V1V2"<<endl;
                    return Distance_from_point_to_edge(P, gnp.vertices[V1], gnp.vertices[V2], distance);
                    
                }
            }
        }
        else {
            
            //Face 0 (V0V1V2V3), V0V1 or V2V3 is closest
            if (gnp.faces[3].N.dot(V0P) > Zero) {
                
                //V0V1 is closest
                //hout<<"P is closest to V0V1"<<endl;
                return Distance_from_point_to_edge(P, gnp.vertices[V0], gnp.vertices[V1], distance);
            }
            else {
                if (gnp.faces[5].N.dot(P - gnp.vertices[V3]) > Zero) {
                    
                    //V2V3 is closest
                    //hout<<"P is closest to V2V3"<<endl;
                    return Distance_from_point_to_edge(P, gnp.vertices[V2], gnp.vertices[V3], distance);
                }
                else {
                    
                    //Face V0V1V2V3 (0) is closest
                    //hout<<"P is closest to V0V1V2V3 (0)"<<endl;
                    //Calculate distance from plane to face
                    distance = gnp.faces[F].distance_to(P);
                    
                    //Calcualte point in plane
                    Point_3D P_proj = P - gnp.faces[F].N*distance;
                    
                    //Set the flag to 0
                    P_proj.flag = 0;
                    
                    return P_proj;
                }
            }
        }
    }
    
}
//This function calculates the distance from an edge to a point P, it also calculates the point on that
//edge closest to P
Point_3D Generate_Network::Distance_from_point_to_edge(const Point_3D &P, const Point_3D &A, const Point_3D &B, double &distance)const
{
    //Calcualte vector AB
    Point_3D AB = B - A;
    
    //Calculate vector AP
    Point_3D AP = P - A;
    
    //Calculate the proyeted distance on edge V1V2 (squared quantity)
    double dot_ = AP.dot(AB);
    double dp2 = dot_*dot_/AB.length2();
    
    //Calculate the distance from AB to the point P
    distance = sqrt(AP.length2() - dp2);
    
    //Calculate the position of the point closest to P in edge AB
    Point_3D Q = A + AB.unit()*sqrt(dp2);
    
    //Set the flag to zero
    Q.flag = 0;
    
    return Q;
}
//This function finds the new position of a CNT point when it is inside a GNP
int Generate_Network::Deal_with_point_inside_gnp(const Geom_sample &geom_sample, const GNP &gnp, const double &cutoff, vector<Point_3D> &new_cnt, Point_3D &new_point)const
{
    //Check if this is the first point of the CNT
    if (new_cnt.empty()) {
        
        //Find the closest face and relocate it
        //hout<<"Find_closest_face_and_relocate"<<endl;
        new_point = Find_closest_face_and_relocate(gnp, new_point, cutoff);
        
        //Nothing more to do, just terminate the function
        return 1;
    }
    
    //Check if the last point in new_cnt:
    //is inside the sample (this ensures that the CNT point is outside the GNP)
    //OR
    //outside the GNP (this is evaluated when the last point in new_cnt is in the boundary layer)
    if (new_cnt.back().flag || !gnp.Is_point_inside_gnp(new_cnt.back())) {
        
        //Find new position of new_point when new_cnt.back() is outside the GNP but
        //new_point is inside the GNP
        new_point = Find_new_position_one_point_inside_gnp(gnp, new_cnt.back(), new_point, cutoff);
    }
    else {
        
        //Here, the point new_cnt.back() is in the boundary layer and inside a GNP
        
        //Attempt to relocate the points new_cnt.back() and new_point
        if (!Find_new_position_two_points_inside_gnp(geom_sample, gnp, cutoff, new_cnt, new_point)) {
            //The two points could not be succesfully relocated so terminate with 0
            //to reject the CNT
            return 0;
        }
        
    }
    
    //If this part of the code was reached, then new_point was successfully relocated
    return 1;
}
//Given a CNT point inside a GNP, this function moves it outside the GNP towards the closest face
Point_3D Generate_Network::Find_closest_face_and_relocate(const GNP &gnp, const Point_3D &new_point, const double &cutoff)const
{
    //Array to save the faces closest to P
    int closest[] = {0,0,0};
    
    //Vector from centroid of GNP to point new_point
    Point_3D CP = new_point - gnp.center;
    
    //Check if face 0 or 1 is closest
    //hout<<"gnp.faces[0].N.dot(CP)="<<gnp.faces[0].N.dot(CP)<<endl;
    if (gnp.faces[0].N.dot(CP) < Zero) {
        
        //Closest to face 1
        closest[0] = 1;
    }
    
    //Check if face 2 or 4 is closest
    //hout<<"gnp.faces[2].N.dot(CP)="<<gnp.faces[2].N.dot(CP)<<endl;
    if (gnp.faces[2].N.dot(CP) > Zero) {
        
        //Closest to face 2
        closest[1] = 2;
    }
    else {
        
        //Closest to face 4
        closest[1] = 4;
    }
    
    //Check if face 3 or 5 is closest
    //hout<<"gnp.faces[3].N.dot(CP)="<<gnp.faces[3].N.dot(CP)<<endl;
    if (gnp.faces[3].N.dot(CP) > Zero) {
        
        //Closest to face 3
        closest[2] = 3;
    }
    else {
        
        //Closest to face 5
        closest[2] = 5;
    }
    
    //Initialize the output distance with the distance to the first face
    double dist = gnp.faces[closest[0]].distance_to(new_point);
    
    //Index of minimum
    int idx = closest[0];
    
    //Iterate over the other two planes
    for (int i = 1; i < 3; i++) {
        
        //Current plane
        int pl = closest[i];
        
        //Distance from P to plane
        double new_dist = gnp.faces[pl].distance_to(new_point);
        
        if (new_dist < dist) {
            dist = new_dist;
            idx = i;
        }
    }
    
    /*for (int i = 0; i < 6; i++) {
        hout<<"P to plane "<<i<<"="<<gnp.faces[i].distance_to(P)<<endl;
    }
    hout<<"idx="<<idx<<endl;*/
    
    //Calculate the new location of new_point
    Point_3D P = new_point + gnp.faces[idx].N*(cutoff+dist);
    
    return P;
}
//This function finds the new position of a CNT point when it is inside the GNP
//Here P1 is the point outside the GNP (new_cnt.back()) and P2 is the point inside the GNP
Point_3D Generate_Network::Find_new_position_one_point_inside_gnp(const GNP &gnp, const Point_3D &P1, const Point_3D &P2, const double &cutoff)const
{
    
    //Array of vertices to calculate dot products
    int V[] = {0,4,0,0,1,2};
    
    //Vector to store the faces with positive dot product
    vector<int> close_faces;
    
    //Find the faces close to P1
    for (int i = 0; i < 6; i++) {
        
        //Check the sign of the dot product
        if (gnp.faces[i].N.dot(P1 - gnp.vertices[V[i]]) > 0) {
            
            //P1 is above the plane of face i, so add it to the vector of close faces
            close_faces.push_back(i);
        }
    }
    
    //Calculate the vector P1P2
    Point_3D P1P2 = P2 - P1;
    
    //Lambda (anonymous) function to calculate the lambda value
    auto calc_lambda = [](const Point_3D &N, const double &d, const Point_3D &P1, const Point_3D &P1P2){
        return ((-N.dot(P1) - d)/(N.dot(P1P2)));
    };
    
    //Variable to store the index of the intersected face
    int idx = 0;
    
    //Variable to store the lambda value
    double lambda = 0.0;
    
    //Find the lambda value of all remaining faces and choose the largest
    //The face with the largest lambda is the intersected face
    for (int i = 0; i < (int)close_faces.size(); i++) {
        
        //Calculate new lambda
        double new_lambda = calc_lambda(gnp.faces[close_faces[i]].N, gnp.faces[close_faces[i]].coef[3], P1, P1P2);
        
        //Check if an update is needed
        //NOTE: the lambda we need is in [0,1]
        if (new_lambda - 1.0 <= Zero && new_lambda >= Zero && new_lambda > lambda) {
            lambda = new_lambda;
            idx = close_faces[i];
        }
    }
    
    //Calculate the distance fom the point inside the GNP to the intersected face
    double dist = gnp.faces[idx].distance_to(P2);
    
    //Move the point P2 towrads the closest face a distance dist+cutoff, so the it
    //takes its new position
    Point_3D P_out = P2 + gnp.faces[idx].N*(dist+cutoff);
    
    return P_out;
}
//This function finds the position of a CNT point when it is inside the GNP and the previous point
//in the CNT is in the boundary layer
//Here P1 is the point in the boundary layer (new_cnt.back()) and P2 is the new point inside the GNP
int Generate_Network::Find_new_position_two_points_inside_gnp(const Geom_sample &geom_sample, const GNP &gnp, const double &cutoff, vector<Point_3D> &new_cnt, Point_3D &new_point)const
{
    //Get the GNP faces that intersect the boundary
    vector<int> faces = Find_gnp_faces_intersecting_boundary(geom_sample, gnp);
    
    //Find the face closest to the CNT point new_cnt.back()
    //Initialize variables with the first face in the vector
    int idx = faces[0];
    double dist = gnp.faces[idx].distance_to(new_cnt.back());
    
    //Scan the rest of the faces
    for (int i = 1; i < (int)faces.size(); i++) {
        
        //Get the distance from point to plane
        double new_dist = gnp.faces[faces[i]].distance_to(new_cnt.back());
        
        //Check if an update is needed
        if (new_dist < dist) {
            dist = new_dist;
            idx = faces[i];
        }
    }
    
    //Move the point in the direction towards the closest face
    new_cnt.back() = new_cnt.back() + gnp.faces[idx].N*(dist+cutoff);
    
    //Check if the point new_cnt.back() is still outside
    if (Is_point_inside_cuboid(geom_sample.sample, new_cnt.back())) {
        //The point is now inside the sample, this will make a repositioning more difficult,
        //so terminate with 0 so that the CNT is rejected
        return 0;
    }
    
    //The point new_cnt.back() is still outside, so move the new point
    dist = gnp.faces[idx].distance_to(new_point);
    new_point = new_point + gnp.faces[idx].N*(dist+cutoff);
    
    return 1;
}
//This function finds the faces of a GNP that are intersecting a sample boundary
vector<int> Generate_Network::Find_gnp_faces_intersecting_boundary(const Geom_sample &geom_sample, const GNP &gnp)const
{
    //Vector to store the intersecting faces
    vector<int> faces;
    
    //Array to store the location of each GNP vertex
    int locations[8];
    
    //Fill the locations array
    for (int i = 0; i < 8; i++) {
        locations[i] = Is_point_inside_cuboid(geom_sample.sample, gnp.vertices[i]);
    }
    
    //Create a lambda function
    auto add_to_vector = [](const int locations[], const int &F, const int &V1, const int &V2, const int &V3, const int &V4, vector<int> &faces){

        int sum = locations[V1] + locations[V2] + locations[V3] + locations[V4];
        if (sum != 0 && sum != 4) {
            faces.push_back(F);
        }
    };
    
    //Check sums face by face and add face if needed
    //Face 0
    add_to_vector(locations,0,0,1,2,3,faces);
    //Face 1
    add_to_vector(locations,1,4,5,6,7,faces);
    //Face 2
    add_to_vector(locations,2,3,7,4,0,faces);
    //Face 3
    add_to_vector(locations,3,0,1,5,4,faces);
    //Face 4
    add_to_vector(locations,4,1,2,6,5,faces);
    //Face 5
    add_to_vector(locations,5,2,3,7,6,faces);
    
    return faces;
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//This function generates output data files if requested
int Generate_Network::Output_data_files(const Geom_sample &geom_sample, const Output_data_flags &out_flags, vector<Point_3D> &points_cnt, vector<double> &radii_out, vector<vector<long int> > &structure, vector<GNP> &gnps)const
{
    //Create printer object in case it might be used
    Printer Pr;
    
    //Output the sample geometry data if any data is requested
    if (out_flags.gnp_data || out_flags.cnt_data) {
        
        //Open file
        ofstream otec("sample_geom.csv");
        
        //Output corner of sample
        otec<<geom_sample.sample.poi_min.str(out_flags.prec_cnt)<<endl;
        
        //Output length of sample along each side
        otec<<geom_sample.sample.len_x<<", "<<geom_sample.sample.wid_y<<", "<<geom_sample.sample.hei_z<<endl;
        
        //Close file
        otec.close();
    }
    
    //Check if output data files were requested for generated GNPs
    if (out_flags.gnp_data == 1) {
        
        //Print the GNP needed to generate them in Abaqus
        Pr.Print_gnp_data(gnps, out_flags.prec_gnp, "gnp_data.csv");
    }
    else if (out_flags.gnp_data == 2) {
        
        //Print the four vertices of a GNP needed to generate them in Abaqus
        Pr.Print_4_vertices_gnps(gnps, out_flags.prec_gnp, "gnp_vertices.csv");
    }
    
    //Check if output data files were requested for generated CNTs
    if (out_flags.cnt_data) {
        
        //Print the CNT points needed to generated them in Abaqus
        Pr.Print_cnt_points_and_structure(geom_sample.sample, structure, points_cnt, radii_out, out_flags.prec_cnt, "cnt_coordinates.csv", "cnt_struct.csv");
    }
    
    return 1;
}
