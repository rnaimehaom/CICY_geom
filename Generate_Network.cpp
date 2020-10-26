//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
int Generate_Network::Generate_nanoparticle_network(const Simu_para &simu_para, const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const GNP_Geo &gnp_geo, const Cutoff_dist &cutoffs, const Tecplot_flags &tec360_flags, vector<Point_3D> &cpoints, vector<double> &cnts_radius_out, vector<vector<long int> > &cstructures, vector<GNP> &gnps)const
{
    //Vector of storing the CNT points
    vector<vector<Point_3D> > cnts_points;
    //Vector for radii, internal variable
    vector<double> cnts_radius_in;
    
    double gnp_vol_tot = 0, gnp_wt_tot = 0;
    if (simu_para.particle_type == "CNT_wires") {
        
        //Generate a network defined by points
        if (!Generate_cnt_network_threads_mt(simu_para, geom_sample, nanotube_geo, cutoffs, cnts_points, cnts_radius_in)) {
            hout << "Error in generating a CNT network" << endl;
            return 0;
        }
        
    } else if (simu_para.particle_type == "GNP_cuboids") {
        
        //Generate a GNP network
        vector<vector<int> > sectioned_domain_gnp;
        if (!Generate_gnp_network_mt(simu_para, gnp_geo, geom_sample, cutoffs, gnps, sectioned_domain_gnp, gnp_vol_tot, gnp_wt_tot)) {
            hout << "Error in generating a GNP network" << endl;
            return 0;
        }
        
    } else if (simu_para.particle_type == "GNP_CNT_mix") {
        
        //Generate a GNP network
        vector<vector<int> > sectioned_domain_gnp;
        if (!Generate_gnp_network_mt(simu_para, gnp_geo, geom_sample, cutoffs, gnps, sectioned_domain_gnp, gnp_vol_tot, gnp_wt_tot)) {
            hout << "Error in generating a GNP network for mixed particles" << endl;
            return 0;
        }
        
        //Generate a CNT network within a GNP network
        if (!Generate_cnt_network_threads_among_gnps_mt(simu_para, gnp_geo, nanotube_geo, geom_sample, cutoffs, gnps, sectioned_domain_gnp, gnp_vol_tot, gnp_wt_tot, cnts_points, cnts_radius_in)) {
            hout<<"Error in generating a CNT network for mixed particles"<<endl;
            return 0;
        }
        
    } else if (simu_para.particle_type == "Hybrid_particles") {
    } else {
        hout << "Error: the type of particles shoud be one of the following: CNT_wires, GNP_cuboids, Hybrid_particles or GNP_CNT_mix. Input value was: " << simu_para.particle_type << endl;
        return 0;
    }
    
    //---------------------------------------------------------------------------
    //Check if Tecplot visualization files were requested for CNTs
    if (tec360_flags.generated_cnts) {
        
        //A new Tecplot_Export object
        Tecplot_Export *Tecexpt = new Tecplot_Export;
        
        if (tec360_flags.generated_cnts == 1) {
            
            //Export the CNT network as threads in Tecplot
            if(Tecexpt->Export_network_threads(geom_sample.sample, cnts_points)==0) return 0;
        }
        else {
            
            //Export the CNT network as 3D "cylinders" in Tecplot
            if(Tecexpt->Export_cnt_network_meshes(tec360_flags.generated_cnts, geom_sample.sample, cnts_points, cnts_radius_in)==0) return 0;
        }
        
        delete Tecexpt;
    }
    
    //---------------------------------------------------------------------------
    //If there are CNTS, transform the 2D cnts_points into 1D cpoints and 2D cstructures
    //Also, remove the CNTs in the boundary layer
    if (simu_para.particle_type != "GNP_cuboids") {
        if(Transform_points_cnts(geom_sample, nanotube_geo, cnts_points, cpoints, cnts_radius_in, cnts_radius_out, cstructures)==0) {
            hout<<"Error in Transform_points for CNTs."<<endl; return 0;
        }
        if (!Recalculate_vol_fraction_cnts(geom_sample, simu_para, nanotube_geo, cpoints, cnts_radius_out, cstructures)) {
            hout<<"Error in Recalculate_vol_fraction_cnts."<<endl; return 0;
        }
        hout<<endl;
    }
    
    //---------------------------------------------------------------------------
    //Check if Tecplot visualization files were requested for CNTs
    if (tec360_flags.generated_cnts) {
        
        //A new Tecplot_Export object
        Tecplot_Export *Tecexpt = new Tecplot_Export;
        
        if (tec360_flags.generated_cnts == 1) {
            
            string filename = "CNT_wires_sample.dat";
            ofstream otec(filename.c_str());
            otec << "TITLE = \"Wires\"" << endl;
            otec << "VARIABLES = X, Y, Z" << endl;
            
            //---------------------------------------------------------------------------
            //Export 3D nanotube threads
            if(Tecexpt->Export_nano_threads(otec, cpoints, cstructures)==0) return 0;
            otec.close();
        }
        else {
            
            //Export the CNT network as 3D "cylinders" in Tecplot
            //if(Tecexpt->Export_cnt_network_meshes(tec360_flags.generated_cnts, geom_sample.sample, cnts_points, cnts_radius_in)==0) return 0;
        }
        
        delete Tecexpt;
    }
    
    //---------------------------------------------------------------------------
    //Check if Tecplot visualization files were requested for GNPs
    if (tec360_flags.generated_gnps) {
        
        //A new Tecplot_Export object
        Tecplot_Export *Tecexpt = new Tecplot_Export;
        
        ofstream otec("Sample.dat");
        otec << "TITLE =\"Sample\"" << endl;
        otec << "VARIABLES = X, Y, Z" << endl;
        if(Tecexpt->Export_cuboid(otec, geom_sample.sample)==0) return 0;
        otec.close();
        
        ofstream otec2("GNP_cuboids.dat");
        otec2 << "TITLE =\"GNPs\"" << endl;
        otec2 << "VARIABLES = X, Y, Z" << endl;
        string zone_name = "as_generated";
        if(Tecexpt->Export_randomly_oriented_gnps(otec2, gnps)==0) return 0;
        otec2.close();
        
        ofstream otec3("Ext_domain.dat");
        otec3 << "TITLE =\"Extended Domain\"" << endl;
        otec3 << "VARIABLES = X, Y, Z" << endl;
        if(Tecexpt->Export_cuboid(otec3, geom_sample.ex_dom_cnt)==0) return 0;
        otec3.close();
        
        delete Tecexpt;
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
    //Check when 10% is completed
    double vol_completed = 0.1;
    
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
        Point_3D cnt_poi;
        //hout<<"Get_seed_point_mt"<<endl;
        if(Get_point_in_cuboid_mt(geom_sample.ex_dom_cnt, cnt_poi, engine_x, engine_y, engine_z, dist)==0) return 0;
        
        //Check overlapping of the intial point
        int counter = 1;
        //hout<<"Check_penetration (while)"<<endl;
        while (simu_para.penetration_model_flag && !Check_penetration(geom_sample, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, point_overlap_count, point_overlap_count_unique, cnt_poi)) {
            if(Get_point_in_cuboid_mt(geom_sample.ex_dom_cnt, cnt_poi, engine_x, engine_y, engine_z, dist)==0) return 0;
            cnt_seed_count++;                    //record the number of seed generations
            //hout << "Seed deleted" << endl;
            if (counter == MAX_ATTEMPTS) {
                hout << "Too many attempts to resolve overlapping of an intial CNT point (" << counter << " attempts). ";
                hout << cnt_poi.x << ' ' << cnt_poi.y << ' ' << cnt_poi.z << endl;
                return 0;
            }
            counter ++;
        }
        
        //Add the CNT seed to the current CNT vector
        new_cnt.push_back(cnt_poi);
        
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
        
        //Start generation of a CNT siven the generated seed
        for(int i=0; i<step_num; i++)
        {
            //Randomly generate a direction in the spherical coordinates
            //To have the positive Z-axis to be a central axis
            //Then, the radial angle, theta, obeys a normal distribution (theta \in fabs[(-omega,+omega)]) and the zonal angle, phi, obeys a uniform distribution (phi \in (0,2PI))
            double cnt_theta = 0, cnt_phi = 0;
            if(Get_normal_direction_mt(nanotube_geo.angle_max, cnt_theta, cnt_phi, engine_theta, engine_phi, dist)==0) return 0;
            
            //Calculate the rotation matrix for current segment
            multiplier = multiplier*Get_transformation_matrix(cnt_theta, cnt_phi);
            
            //hout << "i="<<i<<"<"<<step_num<<endl;
            cnt_poi = cnt_poi + Get_new_point(multiplier, nanotube_geo.step_length);
            //1 means that point is not the intial point
            cnt_poi.flag = 1;
            
            //Check if the new point, cnt_poi, is inside the extended domain
            if(Is_point_inside_cuboid(geom_sample.ex_dom_cnt, cnt_poi))
            {
                //If the point is inside the extended domain, then check for penetration if non-penetrating
                //model is used, or just add it if penetrating model is used
                
                //Check for penetration
                //hout << "Check penetration ";
                if (!simu_para.penetration_model_flag || Check_penetration(geom_sample, nanotube_geo, cnts_points, global_coordinates, sectioned_domain, cnts_radius, new_cnt, n_subregions, cnt_rad, cutoffs.van_der_Waals_dist, point_overlap_count, point_overlap_count_unique, cnt_poi)) {
                    
                    //---------------------------------------------------------------------------
                    //If the penetrating model is used or if the new point, cnt_poi, was placed
                    //in a valid position (i.e., without penetrating another CNT) then add the point
                    //to the generated CNTs
                    
                    //Calculate the segment length inside the sample and add it to the total CNT length
                    bool is_new_inside_sample;
                    cnt_len = cnt_len + Length_inside_sample(geom_sample, new_cnt.back(), cnt_poi, is_prev_in_sample, is_new_inside_sample);
                    
                    //If the new point, cnt_poi, is inside the sample, then increase the number
                    //of points inside the sample of the new CNT
                    if (is_new_inside_sample) {
                        points_in++;
                    }
                    
                    //For the next iteration of the for loop, cnt_poi will become previous point,
                    //so update the boolean is_prev_in_sample for the next iteration
                    is_prev_in_sample = is_new_inside_sample;
                    
                    //Add the new point to the current CNT
                    new_cnt.push_back(cnt_poi);
                    
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
        //Store the CNT points
        //hout << "Store CNT ";
        //Check that is a valid CNT
        //If size is 0, the CNT was rejected
        //If size is 1, then only the CNT seed was generated in a valid position but the second point
        //could not be placed in a valid position (i.e., without interpenetrating another CNT)
        if(new_cnt.size() >= 2 && points_in)
        {
            //Update the volume and weigth of generated CNTs
            vol_sum += cnt_len*cnt_cross_area;
            
            //If the new_cnt vector has at least two points, then it can be added to the rest of the points
            cnts_points.push_back(new_cnt);
            
            //Add the radius to the vector of radii
            cnts_radius.push_back(cnt_rad);
            
            //Perform these operations when the non-overlapping model is used
            if (simu_para.penetration_model_flag) {
                
                //Variable needed for updating global_coordinates
                vector<int> empty;
                
                int CNT = (int)cnts_points.size()-1;
                
                //Scan all points in the CNT that was just generated
                for(int ii=0; ii<(int)new_cnt.size(); ii++)
                {
                    
                    //Add global coordinate
                    global_coordinates.push_back(empty);
                    global_coordinates.back().push_back(CNT);
                    global_coordinates.back().push_back(ii);
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
        
        //Check progress
        if (vol_sum > vol_completed*nanotube_geo.volume) {
            //Get the time
            ct1 = time(NULL);
            
            //Output elapsed time
            hout << "Completed " << vol_completed*100 << " % of target volume. Elapsed time: " << (int)(ct1-ct0) << " secs." <<endl;
            
            //When the next 10% is completed send another message
            vol_completed = vol_completed + 0.1;
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
    n_subregion[0] = max(1, (int)(sample_geom.len_x/sample_geom.gs_minx));
    //Number of subregions along y
    n_subregion[1] = max(1, (int)(sample_geom.wid_y/sample_geom.gs_miny));
    //Number of subregions along z
    n_subregion[2] = max(1, (int)(sample_geom.hei_z/sample_geom.gs_minz));
    
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
int Generate_Network::Check_penetration(const Geom_sample &geom_sample, const Nanotube_Geo &nanotube_geo, const vector<vector<Point_3D> > &cnts, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain, const vector<double> &radii, const vector<Point_3D> &cnt_new, const int n_subregions[], const double &cnt_rad, const double &d_vdw, int &point_overlap_count, int &point_overlap_count_unique, Point_3D &point)const
{
    //Get the sub-region the point belongs to
    //hout<<"Get_subregion 0"<<endl;
    int subregion = Get_cnt_point_subregion(geom_sample, n_subregions, point);
    
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
    //Need to keep the count of the number of attempts outside the for-loop
    int attempts;
    //I move the point up to max_attempts times. If there is still penetration then I delete it
    for (attempts = 0; attempts <= MAX_ATTEMPTS; attempts++) {
        
        //Check if there are any penetrations in the corresponding sub-region
        //hout<<"Get_penetrating_points"<<endl;
        Get_penetrating_points(cnts, global_coordinates, sectioned_domain[subregion], radii, cnt_rad+d_vdw, point, affected_points, cutoffs_p, distances);
        
        //--------------------------------------------------------------------------------------------
        //Check if there are any penetrating points
        if (affected_points.size()) {
            //Update the counter of overlaps
            point_overlap_count++;
            
            //Update the counter of overlapping points only when an overlapping point was found the first time
            //i.e. attempts = 0
            if (!attempts) {
                point_overlap_count_unique++;
            }
            
            //If this is the last iteration and there are still affected points then the point could not be accommodated
            if (attempts == MAX_ATTEMPTS) {
                //hout << "Deleted CNT number " << cnts.size() << " of size " << cnt_new.size();
                //hout << " (reached maximum number of attempts for relocation)" << endl;//*/
                return 0;
            }
            
            //Find the new point
            //hout << "Point " << global_coordinates.size()-1+cnt_new.size() << " in CNT " << cnts.size() << " is overlapping." <<endl;
            //hout << "Moved a point from initial position (" << point.x << ", " << point.y << ", " << point.z << ")." << endl;
            Move_point(cnt_new, point, cutoffs_p, distances, affected_points);
            //hout << "Moved a point to final position (" << point.x << ", " << point.y << ", " << point.z << ")." << endl;
            
            //Check that the new point is within the permited orientation repect to the previous segment
            //hout<<"Check_segment_orientation"<<endl;
            if (!Check_segment_orientation(point, cnt_new)) {
                //hout << "Deleted CNT number " << cnts.size() << " of size " << cnt_new.size();
                //hout << " (the point is not in a valid orientation)" << endl;//*/
                //When not in a valid position it cannot be moved again so a new CNT is needed
                return 0;
            }
            
            //Need to update point sub-region as it could be relocated to a new sub-region
            //hout<<"Get_subregion 1"<<endl;
            subregion = Get_cnt_point_subregion(geom_sample, n_subregions, point);
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
            
        } else {
            //if the size of affected_points is zero, then terminate the function
            return 1;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
//This function returns the subregion a point belongs to
int Generate_Network::Get_cnt_point_subregion(const Geom_sample &geom_sample, const int n_subregions[], const Point_3D &point)const
{
    if (Is_point_inside_cuboid(geom_sample.sample, point)) {
        //These variables will give me the region cordinates of the region that a point belongs to
        int a, b, c;
        //Calculate the region-coordinates
        a = (int)((point.x-geom_sample.origin.x)/geom_sample.gs_minx);
        //Limit the value of a as it has to go from 0 to n_subregions[0]-1
        if (a == n_subregions[0]) a--;
        b = (int)((point.y-geom_sample.origin.y)/geom_sample.gs_miny);
        //Limit the value of b as it has to go from 0 to n_subregions[1]-1
        if (b == n_subregions[1]) b--;
        c = (int)((point.z-geom_sample.origin.z)/geom_sample.gs_minz);
        //Limit the value of c as it has to go from 0 to n_subregions[2]-1
        if (c == n_subregions[2]) c--;
        
        return (a + (b*n_subregions[0]) + (c*n_subregions[0]*n_subregions[1]));
    } else {
        //If the point is in the boundary layer, then there is no need to calculate its sub-region
        return -1;
    }
}
//---------------------------------------------------------------------------
//This functions iterates over a sub-region and determines if there are any penetrating points
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
//This function moves a point according to the number of points it is overlapping
void Generate_Network::Move_point(const vector<Point_3D> &cnt_new, Point_3D &point, vector<double> &cutoffs, vector<double> &distances, vector<Point_3D> &affected_points)const
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
    
    //Get the penetrating point. Its coordinates are in the first (and only) element of affected_points
    Point_3D P = affected_points[0];
    //Calculate direction unit vector. distances[0] already has the lenght of the vector form point to P
    Point_3D direction = (point - P)/distances[0];
    //So the new point is P_old + d*n. d is the cutoff and n the unit vector
    point = P + direction*(cutoffs[0]+Zero); //The Zero is to avoid machine precision errors. Without it, when comparing
    //the new point with the other points in the same region, the program was judging them to be below the cutoff
    //for the van der Waals distance. Even though they were in the limit. After adding this Zero that issue
    //was eliminated
}
//---------------------------------------------------------------------------
//This function finds the new location for an overlapping point when it overlaps two points
void Generate_Network::Two_overlapping_points(const vector<double> &cutoffs, const vector<Point_3D> &affected_points, Point_3D &point)const
{
    //Point variables
    Point_3D P, Q, R, P1, P2;
    //distance variables
    double a, b, c, d;
    
    //Get the penetrating points.
    P1 = affected_points[0];
    //hout << "(" << affected_points[0].x << ", " << affected_points[0].y << ", " << affected_points[0].z << ")." << endl;
    P2 = affected_points[1];
    //hout << "(" << affected_points[1].x << ", " << affected_points[1].y << ", " << affected_points[1].z << ")." << endl;
    //Calculate P vector
    P = P2 - P1;
    //hout<<"P="<<P.str()<<endl;
    //Calculate Q vector
    Q = point - P1;
    //hout<<"Q="<<Q.str()<<endl;
    //Calculate normal vector PxQ
    R = (P.cross(Q)).cross(P);
    //hout<<"R="<<R.str()<<endl;
    //Sides of the triangle
    a = cutoffs[0];
    b = cutoffs[1];
    c = P1.distance_to(P2);
    //Distance from P1 to M
    d = (b*b - a*a - c*c)/(-2*c);
    //hout<<"a="<<a<<" b="<<b<<" c="<<c<<" d="<<d<<endl;
    //Make P a unit vector
    P = P/sqrt(P.dot(P));
    //Make R a unit vector
    R = R/sqrt(R.dot(R));
    //Calculate new position
    point = P1 + P*(d + Zero) + R*(sqrt(a*a - d*d)+Zero);//The Zero is to avoid machine precision errors. Without it, when comparing
    //the new point with the other points in the same region, the program was judging them to be below the cutoff
    //for the van der Waals distance. Even though they were in the limit. After adding this Zero that issue
    //was eliminated
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
    vector<Point_3D> two_affected_points;
    two_affected_points.push_back(affected_points[i1]);
    two_affected_points.push_back(affected_points[i2]);
    vector<double> two_cutoffs;
    two_cutoffs.push_back(cutoffs[i1]);
    two_cutoffs.push_back(cutoffs[i2]);
    
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
//---------------------------------------------------------------------------
//Calculate the effective portion (length) which falls into the given region defined by a cuboid
double Generate_Network::Length_inside_sample(const Geom_sample &geom_sample, const Point_3D &prev_point, const Point_3D &new_point, const bool &is_prev_inside_sample, bool &is_new_inside_sample)const
{
    //Check if the new point is inside the sample
    is_new_inside_sample = Is_point_inside_cuboid(geom_sample.sample, new_point);
    
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
            Point_3D boundary = Find_intersection_at_boundary(geom_sample, new_point, prev_point);
            
            //Return the length from the boundary point towards the inside point
            return prev_point.distance_to(boundary);
        }
    }
    else {
        if (is_new_inside_sample) {
            
            //The previous point is outside the sample, while the new point is inside the sample
            //Find the intersecting point at the boundary
            Point_3D boundary = Find_intersection_at_boundary(geom_sample, prev_point, new_point);
            
            //Return the length from the boundary point towards the inside point
            return new_point.distance_to(boundary);
            
        }
        else {

            //Both points are outside, thus there is zero length inside the sample
            return 0.0;
        }
    }
}
//---------------------------------------------------------------------------
//This function adds a point to a region so penetration can be checked
void Generate_Network::Add_cnt_point_to_overlapping_regions(const Geom_sample &geom_sample, Point_3D point, long int global_num, const int n_subregions[], vector<vector<long int> > &sectioned_domain)const
{
    //A point is added only if it is in the composite domain
    //If the point is in the boundary layer, overlapping is not important
    if (Is_point_inside_cuboid(geom_sample.sample, point)) {
        //Save coordinates of the point
        double x = point.x;
        double y = point.y;
        double z = point.z;
        
        //These variables will give me the region cordinates of the region that a point belongs to
        int a, b, c;
        //Calculate the region-coordinates
        a = (int)((x-geom_sample.origin.x)/geom_sample.gs_minx);
        //Limit the value of a as it has to go from 0 to n_subregions[0]-1
        if (a == n_subregions[0]) a--;
        b = (int)((y-geom_sample.origin.y)/geom_sample.gs_miny);
        //Limit the value of b as it has to go from 0 to n_subregions[1]-1
        if (b == n_subregions[1]) b--;
        c = (int)((z-geom_sample.origin.z)/geom_sample.gs_minz);
        //Limit the value of c as it has to go from 0 to n_subregions[2]-1
        if (c == n_subregions[2]) c--;
        
        //These variables are the coordinates of the lower corner of the RVE that defines its geometry
        double xmin = geom_sample.origin.x;
        double ymin = geom_sample.origin.y;
        double zmin = geom_sample.origin.z;
        
        //Coordinates of non-overlaping region the point belongs to
        double x1 = (double)a*geom_sample.gs_minx +  xmin;
        double x2 = x1 + geom_sample.gs_minx;
        double y1 = (double)b*geom_sample.gs_miny +  ymin;
        double y2 = y1 + geom_sample.gs_miny;
        double z1 = (double)c*geom_sample.gs_minz +  zmin;
        double z2 = z1 + geom_sample.gs_minz;
        
        //Initialize flags for overlaping regions
        int fx = 0;
        int fy = 0;
        int fz = 0;
        
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
        
        //Create array for loop over overlaping regions
        int temp[2][3] = { {a+fx, b+fy, c+fz}, {a, b, c}};
        int t;
        
        //In this loop I check all regions a point can belong to when it is in an overlaping zone
        for (int ii = 0; ii < 2; ii++) {
            if (!fx) ii++; //if flag is zero, do this loop only once
            for (int jj = 0; jj < 2; jj++) {
                if (!fy) jj++; //if flag is zero, do this loop only once
                for (int kk = 0; kk < 2; kk++) {
                    if (!fz) kk++; //if flag is zero, do this loop only once
                    //hout <<"a="<<a<<" fx="<<fx<<" b="<<b<<" fy="<<fy<<" c="<<c<<" fz="<<fz;
                    t = Calculate_t(temp[ii][0],temp[jj][1],temp[kk][2],n_subregions[0],n_subregions[1]);
                    //hout<<" t="<<t<<" sectioned_domain["<<t<<"].size()="<<sectioned_domain[t].size();
                    sectioned_domain[t].push_back(global_num);
                    //hout<<'.'<<endl;
                }
            }
        }


    }
}
//---------------------------------------------------------------------------
//Calculates the region to which a point corresponds
int Generate_Network::Calculate_t(int a, int b, int c, int sx, int sy)const
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
int Generate_Network::Add_cnts_inside_sample(const struct Geom_sample &geom_sample, const struct Nanotube_Geo &nano_geo, const int &CNT_old, vector<Point_3D> &cnt, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures, long int &point_count, int &cnt_count)const
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
    
    /*
    if (!is_first_inside_sample) {
        vector<vector<Point_3D> > tmp(1,cnt);
        Tecplot_Export tec;
        string filename = "CNT_" + to_string(CNT_old) + ".dat";
        tec.Export_network_threads(geom_sample.sample, tmp, filename);
    }*/
    
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
                if (!Add_cnt_segment(geom_sample, is_first_inside_sample, start, end, min_points, CNT_old, cnt, cpoints, radii_in, radii_out, cstructures, point_count, cnt_count)) {
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
    if (last_inside == cnt_points-1) {
        
        //Set end index as one after the last valid index
        //In this way, when adding CNTs to the structure, the last CNT point is added
        end = cnt_points;
        
        //If the last point of the CNT was inside the sample, then add a new segment
        //This was not done becuase, in the for loop, a segement is added only when it finds a point
        //outside the sample
        //Then, check if there are are enough points and, if so, add the current CNT segment to the data structures
        if (!Add_cnt_segment(geom_sample, is_first_inside_sample, start, end, min_points, CNT_old, cnt, cpoints, radii_in, radii_out, cstructures, point_count, cnt_count)) {
            hout<<"Error when adding a CNT segment."<<endl;
            return 0;
        }
    }
    
    return 1;
}
//---------------------------------------------------------------------------
int Generate_Network::Add_cnt_segment(const Geom_sample &geom_sample, const bool &is_first_inside_sample, const int &start, const int &end, const int &min_points, const int &CNT_old, vector<Point_3D> &cnt, vector<Point_3D> &cpoints, vector<double> &radii_in, vector<double> &radii_out, vector<vector<long int> > &cstructures, long int &point_count, int &cnt_count)const
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
    //Find_intersection_at_boundary(const struct Geom_sample &geom_sample, const Point_3D &p_outside, const Point_3D &p_inside)
    Point_3D boundary = Find_intersection_at_boundary(geom_sample, p_outside, p_inside);
    
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
Point_3D Generate_Network::Find_intersection_at_boundary(const Geom_sample &geom_sample, const Point_3D &p_outside, const Point_3D &p_inside)const
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
    if ( (p_outside.x - geom_sample.origin.x) < Zero ) {
        
        //Calculate the lambda value
        lambda_x = calc_lambda(geom_sample.origin.x, p_outside.x, T.x);
    }
    //x-right boundary
    else if ( (geom_sample.x_max - p_outside.x) < Zero ) {
        
        //Calculate the lambda value
        lambda_x = calc_lambda(geom_sample.x_max, p_outside.x, T.x);
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
    if ( (p_outside.y - geom_sample.origin.y) < Zero ) {
        
        //Calculate the lambda value
        lambda_y = calc_lambda(geom_sample.origin.y, p_outside.y, T.y);
    }
    //y-right boundary
    else if ( (geom_sample.y_max - p_outside.y) < Zero ) {
        
        //Calculate the lambda value
        lambda_y = calc_lambda(geom_sample.y_max, p_outside.y, T.y);
    }
    //hout<<"lambda_y="<<lambda_y<<endl;
    
    //Check if any of the z-boundaries is intersected
    double lambda_z = -1.0;
    //z-left boundary
    if ( (p_outside.z - geom_sample.origin.z) < Zero ) {
        
        //Calculate the lambda value
        lambda_z = calc_lambda(geom_sample.origin.z, p_outside.z, T.z);
    }
    //z-right boundary
    else if ( (geom_sample.z_max - p_outside.z) < Zero ) {
        
        //Calculate the lambda value
        lambda_z = calc_lambda(geom_sample.z_max, p_outside.z, T.z);
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
//---------------------------------------------------------------------------
//Generate a random value through a probability distribution function
int Generate_Network::Get_random_value_mt(const string &dist_type, mt19937 &engine, uniform_real_distribution<double> &dist, const double &min, const double &max, double &value)const
{
    //Check if limits are correctly defined
    if(min>max) { hout << "Error, the minimum value is larger than the maximum value (Get_random_value)!" << endl; return 0; }
    
    //Check if the interval has 0 length
    if (max == min) {
        //In this case, value is either of the limits
        //To be consistent with the formulation below, value is set equal to min
        value = min;
        return 1;
    }
    
    if(dist_type=="uniform")    //uniform distribution
    {
        value = (max-min)*dist(engine) + min;
    }
    else if(dist_type=="normal")    //normal distribution
    {
        double sum=0;
        for(int i=0; i<12; i++)
        {
            sum = sum + dist(engine);
        }
        value = (max-min)*sum/12.0 + min;
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
int Generate_Network::Get_normal_direction_mt(const double &omega, double &cnt_theta, double &cnt_phi, mt19937 &engine_theta, mt19937 &engine_phi, uniform_real_distribution<double> &dist)const
{
    
    //theta centers around 0 and obeys a normal distribution in (-omega, +omega)
    double sum=0;
    for(int i=0; i<12; i++)
    {
        sum = sum + dist(engine_theta);
    }
    cnt_theta = fabs(omega*(sum/6 - 1));
    
    //phi satisfies a uniform distribution in (0, 2PI)
    cnt_phi = 2.0*PI*dist(engine_phi);//*/
    
    return 1;
}
//---------------------------------------------------------------------------
//Transform angles into matrix
MathMatrix Generate_Network::Get_transformation_matrix(const double &theta, const double &phi)const
{
    //M = M_phi*M_theta
    //          |cos(phi) -sin(phi) 0|
    // M_phi  = |sin(phi)  cos(phi) 0|
    //          |   0         0     1|
    //
    //          | cos(theta)  0  sin(theta)|
    // M_theta = |     0      1      0    |
    //          |-sin(theta)  0  cos(theta)|
    //Calculate the matrix elements directly, instead of multiplying two matrices
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
        //Tecplot_Export tec;
        //string str =  "GNP_"+to_string(gnps.size())+".dat";
        //tec.Export_singlegnp(gnp, str);
        
        //Flag to determine if new gnp was rejected
        bool rejected = false;
        
        //Check if we are allowing penetrations
        if (simu_para.penetration_model_flag) {
            //Penetrations are not allowed
            
            //Check if the GNP penetrates another GNP
            //hout<<"Deal_with_gnp_interpenetrations"<<endl;
            if (!Deal_with_gnp_interpenetrations(geom_sample, cutoffs, gnps, n_subregion, gnp, sectioned_domain, rejected)) {
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
            bool is_all_outside = false;
            
            //Calculate (or approximate) the generated volume
            //hout<<"Calculate_generated_gnp_vol_and_update_total_vol"<<endl;
            if (!Calculate_generated_gnp_vol_and_update_total_vol(gnp_geo, geom_sample, gnp, gnp_vol_tot, is_all_outside)) {
                hout<<"Error in Calculate_generated_gnp_volume."<<endl;
                return 0;
            }
            
            //Only add the GNP if it is at least partially inside
            if (!is_all_outside) {
                
                //Update the GNP flag
                gnp.flag = (int)gnps.size();
                
                //Add the current particle to the vector of particles
                gnps.push_back(gnp);
            }
            else {
                
                //The GNP is completely outside so increase the count of ignored GNPs
                gnp_ignored_count++;
                
            }
            
        }
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
    n_subregion[0] = max(1, (int)(sample_geom.len_x/sample_geom.gs_minx));
    //Number of subregions along y
    n_subregion[1] = max(1, (int)(sample_geom.wid_y/sample_geom.gs_miny));
    //Number of subregions along z
    n_subregion[2] = max(1, (int)(sample_geom.hei_z/sample_geom.gs_minz));
    
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
    
    //The GNP in local coordinates it is assumed to have its center at the origin, thus its vertices
    //are at +/- l/2 in x and y, and +/- t in z in local coordinates
    //After applying the displacement and rotation (that are already in the gnp variable),
    //vertices are mapped to the global coordinates
    //
    //Top surface
    gnp.vertices[0] = (Point_3D( gnp.l/2, gnp.l/2,  gnp.t/2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[1] = (Point_3D(-gnp.l/2, gnp.l/2,  gnp.t/2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[2] = (Point_3D(-gnp.l/2,-gnp.l/2,  gnp.t/2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[3] = (Point_3D( gnp.l/2,-gnp.l/2,  gnp.t/2)).rotation(gnp.rotation, gnp.center);
    //Bottom surface
    gnp.vertices[4] = (Point_3D( gnp.l/2, gnp.l/2, -gnp.t/2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[5] = (Point_3D(-gnp.l/2, gnp.l/2, -gnp.t/2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[6] = (Point_3D(-gnp.l/2,-gnp.l/2, -gnp.t/2)).rotation(gnp.rotation, gnp.center);
    gnp.vertices[7] = (Point_3D( gnp.l/2,-gnp.l/2, -gnp.t/2)).rotation(gnp.rotation, gnp.center);
    
    //Get the plane equations for the six faces
    if (!Update_gnp_plane_equations(gnp)) {
        hout<<"Error in Generate_gnp when calling Update_gnp_plane_equations."<<endl;
        return 0;
    }
    
    return 1;
}
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
int Generate_Network::Deal_with_gnp_interpenetrations(const Geom_sample &geom_sample, const Cutoff_dist &cutoffs, const vector<GNP> &gnps, const int n_subregions[], GNP &gnp_new, vector<vector<int> > &sectioned_domain, bool &rejected)const
{
    //Variable to count the number of attempts
    int attempts = 0;
    
    //Flag to determine if gnp_new was moved
    bool displaced;
    
    //Varibale to store all GNP subregions
    set<int> subregions;
    
    //hout<<endl<<"GNP_new ="<<gnps.size()<<" center="<<gnp_new.center.str()<<endl;
    //Keep moving gnp_new as long as the number of attempts does not exceed the maximum allowed
    while (attempts <= MAX_ATTEMPTS) {
        
        //Empty the subregions set
        subregions.clear();
        
        //Find the subregions gnp_new occupies
        //hout<<"attempts="<<attempts<<endl;
        if (!Get_gnp_subregions(geom_sample, gnp_new, n_subregions, subregions)) {
            hout<<"Error in Get_gnp_subregions"<<endl;
            return 0;
        }
        
        //Scan the subregions the GNP occupies to determine if it might be close enough to another GNP
        //Get all GNPs it might penetrate
        set<int> gnp_set;
        //hout<<"Get_gnps_in_subregions"<<endl;
        if (!Get_gnps_in_subregions(sectioned_domain, subregions, gnp_set)) {
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
                
                if (!Add_valid_gnp_to_subregions((int)gnps.size(), subregions, sectioned_domain)) {
                    hout<<"Error in Add_valid_gnp_to_subregions after displacement."<<endl;
                    return 0;
                }
                
                return 1;
            }
            //If the GNP was moved, then one more iteration is needed to check that the new
            //position is a valid position
            //else{hout<<"displaced gnp_new.center="<<gnp_new.center.str()<<endl;}
        }
        else {
            
            //If no GNPs that could be to close or interpenetrating gnp_new, then gnp_new is in a valid
            //position
            //Thus, set the rejected flag to false
            rejected = false;
            
            //Also, since gnp_new was successfully generated, add it to all corresponding sub regions
            if (!Add_valid_gnp_to_subregions((int)gnps.size(), subregions, sectioned_domain)) {
                hout<<"Error in Add_valid_gnp_to_subregions (no displacement needed)."<<endl;
                return 0;
            }
            
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
    int a = (int)((new_point.x-geom_sample.origin.x)/geom_sample.gs_minx);
    //Limit the value of a as it has to go from 0 to n_subregions[0]-1
    if (a == n_subregions[0]) a--;
    int b = (int)((new_point.y-geom_sample.origin.y)/geom_sample.gs_miny);
    //Limit the value of b as it has to go from 0 to n_subregions[1]-1
    if (b == n_subregions[1]) b--;
    int c = (int)((new_point.z-geom_sample.origin.z)/geom_sample.gs_minz);
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
    double x1 = (double)a*geom_sample.gs_minx +  geom_sample.origin.x;
    double x2 = x1 + geom_sample.gs_minx;
    double y1 = (double)b*geom_sample.gs_miny +  geom_sample.origin.y;
    double y2 = y1 + geom_sample.gs_miny;
    double z1 = (double)c*geom_sample.gs_minz +  geom_sample.origin.z;
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
    
    //Variable to store the total displacement (if the GNP_new needs to be moved)
    Point_3D disp_tot = Point_3D(0.0,0.0,0.0);
    
    //Initialize the displaced flag to false
    displaced = false;
    
    //int el = 0;
    for (set<int>::iterator i = gnp_set.begin(); i!=gnp_set.end(); i++) {
        
        //Get current GNP number
        int GNP_i = *i;
        //hout<<"GNP_i="<<GNP_i<<" el="<<el<<endl;el++;
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
            Tecplot_Export tec;
            string str =  "GNPA.dat";
            tec.Export_singlegnp(gnps[GNP_i], str);
            str = "GNPB.dat";
            tec.Export_singlegnp(gnp_new, str);
        }//*/
        
        //Vector to stor the simplex that encloses the origin in case of interpenetration
        vector<Point_3D> simplex;
        
        //Flags for penetration (p_flag) and touching (t_flag)
        bool p_flag = false;
        bool t_flag = false;
        
        //Check if GNP and GNP_new penetrate each other
        if (!GJK_EPA.GJK(gnps[GNP_i], gnp_new, simplex, p_flag, t_flag)) {
            hout<<"Error in Move_gnps_if_needed when calling GJK"<<endl;
            return 0;
        }

        //Variables to store the penetration depth (PD) and direction vector (N) along which
        //gnp_new needs to move
        double PD;
        Point_3D N;
        
        if (p_flag) {
            //hout<<"Penetration with GNP "<<GNP_i<<endl;
            
            //There is penetration, so then use EPA to fint the penetration depth PD and direction vector N
            if (!GJK_EPA.EPA(gnps[GNP_i].vertices, gnp_new.vertices, simplex, N, PD)) {
                hout<<"Error in Move_gnps_if_needed when calling EPA"<<endl;
                return 0;
            }
            //hout<<"PD="<<PD<<" normal="<<N.str()<<endl;
            
            //Update the total dispacement
            disp_tot = disp_tot + N*(PD + cutoffs.van_der_Waals_dist);
            
            //Change the displaced flag if needed
            if (!displaced) {
                displaced =  true;
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
                disp_tot = disp_tot + N*cutoffs.van_der_Waals_dist;
                
                //Change the displaced flag if needed
                if (!displaced) {
                    displaced =  true;
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
                    
                    //If the separation is less than the van der Waals distance, then gnp_new needs
                    //to be moved, so update the total displacement
                    disp_tot = disp_tot + N*new_dist;
                    
                    //Change the displaced flag if needed
                    if (!displaced) {
                        displaced =  true;
                    }
                }
            }
        }
    }
    
    //Check if a displacement is needed
    if (displaced) {
        
        //hout<<"GNP moved disp_tot="<<disp_tot.str()<<endl;
        //A displacement is needed, then move gnp_new
        if (!Move_gnp(disp_tot, gnp_new)) {
            hout<<"Error in Move_gnp."<<endl;
            return 0;
        }
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
        double dot_ = V.dot(Q_);
        return (len_V2 + dot_*dot_/Q_.length2());
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
    
    
    //Vectors to save the simplices in A and B that are touching
    vector<int> simplexA, simplexB;
    
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
                simplexA.push_back(i);
                simplexB.push_back(j);
            }
            
            //Q is not the origin, so check if the distance from V to line segment OQ
            //is the same (or almost the same) as the reference distance
            else if (abs(dist_ref - dist2(Q)) < Zero) {
                
                //The vertices in gnpA and gnpB are in the edge or face of the Minkowski sum
                //that crosses the origin, so add them to their corresponding simplices
                simplexA.push_back(i);
                simplexB.push_back(j);
            }
        }
    }
    
    //Check if there is a face in any of the simplices found
    if (simplexA.size() == 4 || simplexB.size() == 4) {
        
        //Find one face and the normal to that face
        if (simplexA.size() == 4) {
            
            //gnpA has a face as part of the touch, get the normal
            N = (gnpA.vertices[simplexA[1]] - gnpA.vertices[simplexA[0]]).cross(gnpA.vertices[simplexA[2]] - gnpA.vertices[simplexA[0]]);
        }
        else {
            
            //Only gnpB has a plane as part of the touch
            N = (gnpB.vertices[simplexB[1]] - gnpB.vertices[simplexB[0]]).cross(gnpB.vertices[simplexB[2]] - gnpB.vertices[simplexB[0]]);
        }
        
        //Make the normal a unit vector
        N.make_unit();
        
        //Make sure N goes in the direction from gnpA to gnpB
        if (N.dot(gnpB.center - gnpA.center) < Zero) {
            //Rever the direction
            N = N*(-1);
        }
    }
    else {
        
        //No face is touching, in such case gnpB can be moved in the direction from centerA to centerB
        N = (gnpB.center - gnpA.center).unit();
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
//This function calculates the generated volume of a GNP and adds it to the global variables
int Generate_Network::Calculate_generated_gnp_vol_and_update_total_vol(const GNP_Geo gnp_geom, const Geom_sample &sample_geom, const GNP &gnp, double &gnp_vol_tot, bool &is_all_outside)const
{
    //---------------------------------------------------------------------------
    //Add the volume and weight corresponding to the GNP
    double gnp_vol = gnp.volume;
    
    //Iterate over the eight vertices of the GNP and check if all of them are inside the sample
    for (int i = 0; i < 8; i++) {
        
        //Check if current vertex is outside the sample
        if (!Is_point_inside_cuboid(sample_geom.sample, gnp.vertices[i])) {
            
            //Vertex i is outside the sample, then approximate the GNP volume
            if (!Approximate_gnp_volume_inside_sample(sample_geom, gnp, gnp_vol, is_all_outside)) {
                hout << "Error in Calculate_generated_volume when calling Approximate_gnp_volume_inside_sample." << endl;
                return 0;
            }
            
            //Break the loop as it is not necessary to continue cheking the rest of the vertices
            break;
        }
    }
    
    //hout << "Total volume = " << gnp_vol_tot << ". Approximated volume = " << gnp_vol <<endl;
    //Update the total volume
    gnp_vol_tot = gnp_vol_tot + gnp_vol;
    
    return 1;
}
//This function approximates the volume of a GNP inside the sample depending of the number of points in
//the discretization of its middle plane that are inside the sample
int Generate_Network::Approximate_gnp_volume_inside_sample(const Geom_sample &sample_geom, const GNP &gnp, double &gnp_vol, bool &is_all_outside)const
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
            if (Is_point_inside_cuboid(sample_geom.sample, new_point)) {
                //The point is inside the sample, so increase the number of points inside the sample
                points_in++;
            }
        }
    }
    
    //Check if all points where outside the sample
    if (points_in == 0) {
        
        //All points are outside the sample, so set the is_all_outside flag to true
        is_all_outside = true;
        
        //Directly set the GNP volume to zero
        gnp_vol = 0.0;
    }
    else {
        
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
        
        //---------------------------------------------------------------------------
        //Generate an initial point of the CNT in the extended domain
        //hout<<"Generate_initial_point_of_cnt"<<endl;
        Point_3D new_point;
        if (Generate_initial_point_of_cnt(geom_sample, simu_para, cnts_points, cnts_radius, new_cnt, cnt_rad+cutoffs.van_der_Waals_dist, global_coordinates, sectioned_domain, gnps, sectioned_domain_gnp, n_subregions, new_point, engine_x, engine_y, engine_z, dist) != 1) {
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
        
        //Start generation of a CNT siven the generated seed
        for(int i = 0; i < step_num; i++)
        {
            //Randomly generate a direction in the spherical coordinates
            //To have the positive Z-axis to be a central axis
            //Then, the radial angle, theta, obeys a normal distribution (theta \in fabs[(-omega,+omega)]) and the zonal angle, phi, obeys a uniform distribution (phi \in (0,2PI))
            double cnt_theta = 0, cnt_phi = 0;
            //hout<<"Get_normal_direction_mt"<<endl;
            if(Get_normal_direction_mt(nanotube_geo.angle_max, cnt_theta, cnt_phi, engine_theta, engine_phi, dist)==0) return 0;
            
            //Calculate the rotation matrix for current segment
            multiplier = multiplier*Get_transformation_matrix(cnt_theta, cnt_phi);
            
            //hout << "i="<<i<<"<"<<step_num<<endl;
            new_point = new_point + Get_new_point(multiplier, nanotube_geo.step_length);
            //1 means that point is not the intial point
            new_point.flag = 1;
            
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
                    mixed_interpenetration = Check_mixed_interpenetration(geom_sample, cnts_points, cnts_radius, new_cnt, cnt_rad+cutoffs.van_der_Waals_dist, global_coordinates, sectioned_domain, gnps, sectioned_domain_gnp, n_subregions, new_point);
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
                    cnt_len = cnt_len + Length_inside_sample(geom_sample, new_cnt.back(), new_point, is_prev_in_sample, is_new_inside_sample);
                    
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

                    //if(cnts_points.size() == 300){hout<<"P_(i="<<new_cnt.size()<<")="<<new_point.str()<<endl;}
                    
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
        //Store the CNT points
        //hout << "Store CNT ";
        //Check that is a valid CNT
        //If size is 0, the CNT was rejected
        //If size is 1, then only the CNT seed was generated in a valid position but the second point
        //could not be placed in a valid position (i.e., without interpenetrating another CNT)
        if(new_cnt.size() >= 2 && points_in)
        {
            //Update the volume and weigth of generated CNTs
            vol_sum += cnt_len*cnt_cross_area;
            
            //If the new_cnt vector has at least two points, then it can be added to the rest of the points
            cnts_points.push_back(new_cnt);
            
            //Add the radius to the vector of radii
            cnts_radius.push_back(cnt_rad);
            
            //Perform these operations when the non-overlapping model is used
            if (simu_para.penetration_model_flag) {
                
                //Variable needed for updating global_coordinates
                vector<int> empty;
                
                int CNT = (int)cnts_points.size()-1;
                
                //Scan all points in the CNT that was just generated
                for(int ii=0; ii<(int)new_cnt.size(); ii++)
                {
                    
                    //Add global coordinate
                    global_coordinates.push_back(empty);
                    global_coordinates.back().push_back(CNT);
                    global_coordinates.back().push_back(ii);
                    //Add point to an overlapping region in the vector sectioned_domain
                    //hout<<"Add_cnt_point_to_overlapping_regions"<<endl;
                    Add_cnt_point_to_overlapping_regions(geom_sample, new_cnt[ii], (long int)global_coordinates.size()-1, n_subregions, sectioned_domain);
                    
                }
            }
        }
        else if (!points_in){
            
            //The CNT can be ignored as it is completely outside the sample
            cnt_ignore_count++;
        }
        //hout << "done"<<endl;
        
        //Check progress
        if (vol_sum > vol_completed*nanotube_geo.volume) {
            //Get the time
            ct1 = time(NULL);
            
            //Output elapsed time
            hout << "Completed " << vol_completed*100 << " % of target volume. Elapsed time: " << (int)(ct1-ct0) << " secs." <<endl;
            
            //When the next 10% is completed send another message
            vol_completed = vol_completed + 0.1;
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
int Generate_Network::Generate_initial_point_of_cnt(const Geom_sample &geom_sample, const Simu_para &simu_para, const vector<vector<Point_3D> > &cnts, const vector<double> &radii, vector<Point_3D> &new_cnt, const double &rad_plus_dvdw, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const int n_subregions[], Point_3D &new_point, mt19937 &engine_x, mt19937 &engine_y, mt19937 &engine_z, uniform_real_distribution<double> &dist)const
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
            int interpenetration = Check_mixed_interpenetration(geom_sample, cnts, radii, new_cnt, rad_plus_dvdw, global_coordinates, sectioned_domain_cnt, gnps, sectioned_domain_gnp, n_subregions, new_point);
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
int Generate_Network::Check_mixed_interpenetration(const Geom_sample &geom_sample, const vector<vector<Point_3D> > &cnts, const vector<double> &radii, vector<Point_3D> &new_cnt, const double &rad_plus_dvdw, const vector<vector<int> > &global_coordinates, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, const int n_subregions[], Point_3D &new_point)const
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
            
            //Move the point for this special case
            if (!Deal_with_point_inside_gnp(geom_sample, affected_gnp, rad_plus_dvdw, new_cnt, new_point)) {
                //new_point could not be accomodated so reject it
                return 0;
            }
            
            //Check if the segment has a valid orientation, only when new_cnt.back() is inside the sample
            if (Is_point_inside_cuboid(geom_sample.sample, new_cnt.back())) {
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
            
            //Check if there are any penetrating points
            if (!affected_points.empty()) {
                
                //Move the point depending on the case
                //hout<<"Move_point"<<endl;
                Move_point(new_cnt, new_point, cutoffs_p, distances, affected_points);
                
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
        new_point = Find_closest_face_and_relocate(gnp, new_point, cutoff);
        
        //Nothing more to do, just terminate the function
        return 1;
    }
    
    //Check if the last point in new_cnt:
    //is inside the sample (this ensures that the CNT point is outside the GNP)
    //OR
    //outside the GNP (this is evaluated when the last point in new_cnt is in the boundary layer)
    if (new_cnt.back().flag || !Is_point_inside_gnp(gnp, new_cnt.back())) {
        
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
//This function checks if a point is inside a GNP
bool Generate_Network::Is_point_inside_gnp(const GNP &gnp, const Point_3D &P)const
{
    //This vector is used multiple times
    Point_3D V4P = P - gnp.vertices[4];
    
    //Check if P is between faces 0 and 1
    if (gnp.faces[0].N.dot(P - gnp.vertices[0]) < Zero && gnp.faces[1].N.dot(V4P) < Zero) {
        
        //Check if P is between faces 2 and 4
        if (gnp.faces[2].N.dot(V4P) < Zero && gnp.faces[4].N.dot(P - gnp.vertices[1]) < Zero) {
            
            //Check if P is between faces 3 and 5
            if (gnp.faces[3].N.dot(V4P) < Zero && gnp.faces[5].N.dot(P - gnp.vertices[2]) < Zero) {
                
                //Point P is inside the GNP so return true
                return true;
            }
        }
    }
    
    //If this part of the code is reached, then the point P is outside the GNP
    return false;
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