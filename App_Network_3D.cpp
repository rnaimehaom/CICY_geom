//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GS networks
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
    
    //Deprecated:
    vector<GCH> hybrid_particles;
    vector<vector<long int> > gnps_structure;
    
    //----------------------------------------------------------------------
    //Network Generation with overlapping
    hout << "Generating nanoparticle network......" << endl;
    ct0 = time(NULL);
    Generate_Network *Generator = new Generate_Network;
    if (!Generator->Generate_nanoparticle_network(Init->simu_para, Init->geom_sample, Init->nanotube_geo, Init->gnp_geo, Init->cutoff_dist, Init->vis_flags, points_cnt, radii, structure_cnt, gnps)) {
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
    Shells *Shell = new Shells;
    if (!Shell->Generate_shells(Init->geom_sample, points_cnt, gnps, shells_cnts, shell_gnps)) {
        hout << "Error when generating shells" << endl;
        return 0;
    }
    delete Shell;
    ct1 = time(NULL);
    hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
    //Variable to store the geometry of the observation window
    struct Geom_sample window_geom;
    cuboid window_geo;
    
    for(int i=0; i<=Init->geom_sample.cut_num; i++)
    {
        hout << "============================================================================"<<endl;
        hout << "============================================================================"<<endl;
        hout << "Iteration " << i+1 << " of " << Init->geom_sample.cut_num+1 << endl;
        time_t it0, it1;
        it0 = time(NULL);
        
        //Update observation window geometry
        if (!Update_obseravtion_window_geometry(i, Init->geom_sample, window_geo)) {
            hout<<"Error when updating the geometry for observation window "<<i<<endl;
            return 0;
        }
        
        //----------------------------------------------------------------------
        //These vectors are used to export tecplot files
        vector<long int> empty;
        vector<vector<long int> > all_dead_indices(7,empty);
        vector<vector<long int> > all_percolated_indices(7,empty);
        vector<int> empty_int;
        vector<vector<int> > all_percolated_gnp(7,empty_int);
        vector<vector<int> > all_dead_gnps(7,empty_int);
        
        //----------------------------------------------------------------------
        //Determine the local networks in cutoff windows
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        //From this function I get the internal variables cnts_inside and boundary_cnt
        ct0 = time(NULL);
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
        if (!HoKo->Determine_clusters_and_percolation(Init->simu_para, Init->cutoff_dist, Cutwins->cnts_inside, Contacts->sectioned_domain_cnts, structure_cnt, points_cnt, radii, Cutwins->boundary_cnt, Cutwins->gnps_inside, Contacts->sectioned_domain_gnps, gnps, Cutwins->boundary_gnp, structure_gnp, points_gnp)) {
            hout << "Error when finding clusters and determining percolation" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Find clusters and determine percolation: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Contacts are not needed anymore, so delete the object
        delete Contacts;
        
        //----------------------------------------------------------------------
        //These vectors are used to store the fractions of the different families in the current observation window
        //families_lengths has 8 elements because of the 7 percolated families and the non-percoalted CNTs, the same is true for fractions
        vector<double> families_lengths(8,0);
        vector<double> fractions(8,0);
        vector<double> branches_lengths(7,0);
        
        //Loop over the different clusters so that the direct electrifying algorithm is aplied on each cluster
        Electrical_analysis *Electric_A = new Electrical_analysis;/*
        if (HoKo->clusters_cnt.size() || HoKo->clusters_gch.size()) {
            
            //Perform the electrical analysis to obtain the backbone and calculate the electrical resistance
            ct0 = time(NULL);
            if (!Electric_A->Perform_analysis_on_clusters(Init->simu_para.avoid_resistance, HoKo->family, HoKo, Cutwins, structure_cnt, points_cnt, radii, gnps_structure, points_gnp, window_geom, Init->electric_para, Init->cutoff_dist, hybrid_particles, all_dead_indices, all_percolated_indices, all_dead_gnps, all_percolated_gnp)) {
                hout << "Error when performing electrical analysis" << endl;
                return 0;
            }
            ct1 = time(NULL);
            hout << "Perform electrical analysis time: "<<(int)(ct1-ct0)<<" secs."<<endl;
            
        } else {
            
            //Calculate the resistances and resistivities of the polymer matrix
            vector<double> matrix_resistances;
            if (!Electric_A->Calculate_matrix_resistances(Init->electric_para.resistivity_matrix, window_geom, matrix_resistances)) {
                hout << "Error when calculating matrix resistances for a sample without percolated clusters" << endl;
                return 0;
            }
            
            //Calculate the resistances and resistivities along each direction
            vector<double> resistivities;
            //paralel_resistors is initialized with three empty vectors
            vector<vector<double> > paralel_resistors(3,resistivities);
            if (!Electric_A->Calculate_resistances_and_resistivities(window_geom, matrix_resistances, paralel_resistors, Electric_A->resistors, resistivities)) {
                hout << "Error when calculating matrix resistivities for a sample without percolated clusters" << endl;
                return 0;
            }
            
            //Append resistors to a file
            Printer *P = new Printer;
            P->Append_1d_vec(Electric_A->resistors, "resistors.txt");
            P->Append_1d_vec(resistivities, "resistivities.txt");
            delete P;
            
            hout << "There are no percolated clusters" << endl;
        }*/
        delete Electric_A;
        
        //Calculate the fractions of CNTs that belong to each family and save them to a file
        Clusters_fractions *Fracs = new Clusters_fractions;
        ct0 = time(NULL);
        if (!Fracs->Calculate_fractions(Init->geom_sample, Cutwins->cnts_inside, Cutwins->gnps_inside, structure_cnt, points_cnt, radii, HoKo->isolated, hybrid_particles, HoKo->isolated_gch, all_dead_indices, all_percolated_indices, all_dead_gnps, all_percolated_gnp)) {
            hout << "Error when calculating clusters fractions" << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Calculate fractions time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Delete objects to free memory
        delete Fracs;
        delete Cutwins;

        //----------------------------------------------------------------------
        //Check if visualization files were requested for the backbone
        
        //Delete objects to free memory
        delete HoKo;
        
        //----------------------------------------------------------------------
        //Check if visualization files were requested for the triangulations
        
        it1 = time(NULL);
        hout << "Iteration "<<i+1<<" time: "<<(int)(it1-it0)<<" secs."<<endl;

        
    }
    
    
    return 1;
}
//Update the geometry of the observation window
int App_Network_3D::Update_obseravtion_window_geometry(const int &window, const Geom_sample &sample_geo, cuboid &window_geo)const
{
    //Dimensions of the current observation window
    window_geo.len_x = sample_geo.win_max_x - ((double)window)*sample_geo.win_delt_x;
    window_geo.wid_y = sample_geo.win_max_y - ((double)window)*sample_geo.win_delt_y;
    window_geo.hei_z = sample_geo.win_max_z - ((double)window)*sample_geo.win_delt_z;
    
    //These variables are the coordinates of the lower corner of the observation window
    window_geo.poi_min.x = sample_geo.sample.poi_min.x + (sample_geo.sample.len_x - window_geo.len_x)/2;
    window_geo.poi_min.y = sample_geo.sample.poi_min.y + (sample_geo.sample.wid_y - window_geo.wid_y)/2;
    window_geo.poi_min.z = sample_geo.sample.poi_min.z + (sample_geo.sample.hei_z - window_geo.hei_z)/2;
    
    //Boundaries with maximum values of coordinates
    window_geo.max_x = window_geo.poi_min.x + window_geo.len_x;
    window_geo.max_y = window_geo.poi_min.y + window_geo.wid_y;
    window_geo.max_z = window_geo.poi_min.z + window_geo.hei_z;
    
    return 1;
}
//===========================================================================
