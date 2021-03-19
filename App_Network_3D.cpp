//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Create 3D nanoparticle network. Find backbone, dead branches and isolated clusters. Calculate electrical conductivity of sample (and observation windows)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "App_Network_3D.h"

//Generate 3D nanoparticle network, turn it into a resitor network and find its electrical conductivity
int App_Network_3D::Generate_nanoparticle_resistor_network(Input *Init)const
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
    
    //Shell vectors (used to remove nanoparticles when reduding observation window size)
    vector<vector<int> > shells_cnts;
    
    //----------------------------------------------------------------------
    //Network Generation with overlapping
    hout << "Generating nanoparticle network......" << endl;
    ct0 = time(NULL);
    Generate_Network *Generator = new Generate_Network;
    if (!Generator->Generate_nanoparticle_network(Init->simu_para, Init->geom_sample, Init->nanotube_geo, Init->gnp_geo, Init->cutoff_dist, Init->vis_flags, Init->out_flags, points_cnt, radii, structure_cnt, gnps)) {
        hout<<"Error in Generate_nanoparticle_network."<<endl;
        return 0;
    }
    delete Generator;
    ct1 = time(NULL);
    hout << "Nanotube network generation time: " << (int)(ct1-ct0) <<" secs." << endl;
    //Printer *P = new Printer;
    //P->Print_1d_vec(gnps_point, "gnps_point_00.txt");
    //delete P;
    
    //----------------------------------------------------------------------
    //Vector for GNP shells
    vector<Shell> shell_gnps(gnps.size());
    ct0 = time(NULL);
    Shells *SH = new Shells;
    if (!SH->Generate_shells(Init->geom_sample, points_cnt, gnps, shells_cnts, shell_gnps)) {
        hout << "Error when generating shells" << endl;
        return 0;
    }
    delete SH;
    ct1 = time(NULL);
    hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
    //Variable to store the geometry of the observation window
    cuboid window_geo;
    
    for(int i=0; i<=Init->geom_sample.cut_num; i++)
    {
        hout << "============================================================================"<<endl;
        hout << "============================================================================"<<endl;
        hout << "Iteration " << i+1 << " of " << Init->geom_sample.cut_num+1 << endl;
        time_t it0, it1;
        it0 = time(NULL);
        
        //Update observation window geometry
        //hout<<"Update observation window geometry"<<endl;
        if (!Update_obseravtion_window_geometry(i, Init->vis_flags.window_domain, Init->geom_sample, Init->simu_para.particle_type, window_geo)) {
            hout<<"Error when updating the geometry for observation window "<<i<<endl;
            return 0;
        }
        
        //----------------------------------------------------------------------
        //Determine the local networks in cutoff windows
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        ct0 = time(NULL);
        //hout<<"Extract_observation_window"<<endl;
        if(!Cutwins->Extract_observation_window(i, Init->simu_para.particle_type, Init->geom_sample, window_geo, Init->nanotube_geo, gnps, structure_cnt, radii, points_cnt, shells_cnts, shell_gnps, structure_gnp, points_gnp)) {
            hout << "Error when extracting observation window "<< i+1 << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Extract observation window time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //----------------------------------------------------------------------
        //Determine the local networks inside the cutoff windows
        Contact_grid *Contacts = new Contact_grid;
        ct0 = time(NULL);
        if (!Contacts->Generate_contact_grid(i, Init->simu_para.particle_type, Init->geom_sample, window_geo, Cutwins->cnts_inside, points_cnt, structure_cnt, Cutwins->gnps_inside, gnps)) {
            hout << "Error when generating contact grid" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Generate contact grid time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //----------------------------------------------------------------------
        //Hoshen-Kopelman algorithm
        Hoshen_Kopelman *HoKo = new Hoshen_Kopelman;
        ct0 = time(NULL);
        if (!HoKo->Determine_clusters_and_percolation(i, Init->simu_para, Init->cutoff_dist, Init->vis_flags, Cutwins->cnts_inside, Contacts->sectioned_domain_cnts, structure_cnt, points_cnt, radii, Cutwins->boundary_cnt, Cutwins->gnps_inside, Contacts->sectioned_domain_gnps, gnps, Cutwins->boundary_gnp, structure_gnp, points_gnp)) {
            hout << "Error when finding clusters and determining percolation" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Find clusters and determine percolation: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Contacts are not needed anymore, so delete the object
        delete Contacts;
        
        //Loop over the different clusters so that the direct electrifying algorithm is aplied on each cluster
        Electrical_analysis *EA = new Electrical_analysis;
        if (!EA->Perform_analysis_on_clusters(i, window_geo, Init->simu_para, Init->electric_para, Init->cutoff_dist, Init->vis_flags, Init->out_flags, HoKo, Cutwins, structure_cnt, points_cnt, radii, points_gnp, structure_gnp, gnps)) {
            hout << "Error when performing electrical analysis" << endl;
            return 0;
        }
        
        //Delete objects to free memory
        delete EA;
        delete Cutwins;
        delete HoKo;
        
        it1 = time(NULL);
        hout << "Iteration "<<i+1<<" time: "<<(int)(it1-it0)<<" secs."<<endl;
    }
    
    return 1;
}
//Update the geometry of the observation window
int App_Network_3D::Update_obseravtion_window_geometry(const int &window, const int &window_domain, const Geom_sample &sample_geo, const string particle_type, cuboid &window_geo)const
{
    //Check the particle type
    //If a deposit of CNTs is generated, then the observation window will not change along z
    if (particle_type != "CNTdeposit") {
        
        //Update the dimension along z of the current observation window
        window_geo.hei_z = sample_geo.win_max_z - ((double)window)*sample_geo.win_delt_z;
        
        //Upadte the z-coordinate of the lower corner of the observation window
        window_geo.poi_min.z = sample_geo.sample.poi_min.z + (sample_geo.sample.hei_z - window_geo.hei_z)/2;
    }
    
    //Update dimensions of the current observation window
    window_geo.len_x = sample_geo.win_max_x - ((double)window)*sample_geo.win_delt_x;
    window_geo.wid_y = sample_geo.win_max_y - ((double)window)*sample_geo.win_delt_y;
    //hout<<"window_geo.len_x="<<window_geo.len_x<<" window_geo.wid_y="<<window_geo.wid_y<<" window_geo.hei_z="<<window_geo.hei_z<<endl;
    
    //These variables are the coordinates of the lower corner of the observation window
    window_geo.poi_min.x = sample_geo.sample.poi_min.x + (sample_geo.sample.len_x - window_geo.len_x)/2;
    window_geo.poi_min.y = sample_geo.sample.poi_min.y + (sample_geo.sample.wid_y - window_geo.wid_y)/2;
    //hout<<"window_geo.poi_min.x="<<window_geo.poi_min.x<<" window_geo.poi_min.y="<<window_geo.poi_min.y<<" window_geo.poi_min.z="<<window_geo.poi_min.z<<endl;
    
    //Boundaries with maximum values of coordinates
    window_geo.max_x = window_geo.poi_min.x + window_geo.len_x;
    window_geo.max_y = window_geo.poi_min.y + window_geo.wid_y;
    window_geo.max_z = window_geo.poi_min.z + window_geo.hei_z;
    //hout<<"window_geo.max_x="<<window_geo.max_x<<" window_geo.max_y="<<window_geo.max_y<<" window_geo.max_z="<<window_geo.max_z<<endl;
    
    //Export the window geometry if needed
    if (window_domain) {
        string str = "window_" + to_string(window) + ".vtk";
        VTK_Export VTK_E;
        VTK_E.Export_cuboid(window_geo, str);
    }
    
    return 1;
}
//===========================================================================
