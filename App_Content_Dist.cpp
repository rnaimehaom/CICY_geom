//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate the volume or weight fraction for each observation window
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "App_Content_Dist.h"

//Generate 3D nanoparticle network, turn it into a resitor network and find its electrical conductivity
int App_Content_Dist::Calculate_content_on_each_window(Input *Init)const
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
    //Check if network is generated or read from file
    if (Init->simu_para.create_read_network == "Create_Network" ||
        Init->simu_para.create_read_network == "Network_From_Seeds") {

        //----------------------------------------------------------------------
        //Network Generation with overlapping
        hout << "Generating nanoparticle network......" << endl;
        ct0 = time(NULL);
        Generate_Network* Generator = new Generate_Network;
        if (!Generator->Generate_nanoparticle_network(Init->simu_para, Init->geom_sample, Init->nanotube_geo, Init->gnp_geo, Init->cutoff_dist, Init->vis_flags, Init->out_flags, points_cnt, radii, structure_cnt, gnps)) {
            hout << "Error in Generate_nanoparticle_network." << endl;
            return 0;
        }
        delete Generator;
        ct1 = time(NULL);
        hout << "Network generation time: " << (int)(ct1 - ct0) << " secs." << endl;
        //Printer *P = new Printer;
        //P->Print_1d_vec(gnps_point, "gnps_point_00.txt");
        //delete P;
    }
    else if (Init->simu_para.create_read_network == "Read_Network") {

        //----------------------------------------------------------------------
        //Network Generation from file
        hout << "Generating nanoparticle network from file......" << endl;
        ct0 = time(NULL);
        Read_Network* Reader = new Read_Network;
        if (!Reader->Generate_nanoparticle_network_from_file(Init->simu_para, Init->vis_flags, Init->geom_sample, points_cnt, radii, structure_cnt, gnps)) {
            hout << "Error in Generate_nanoparticle_network_from_file." << endl;
            return 0;
        }
        delete Reader;
        ct1 = time(NULL);
        hout << "Nanoparticle network generation from file time: " << (int)(ct1 - ct0) << " secs." << endl;
    }
    else {
        hout << "Error in Calculate_content_on_each_window. Invalid keyword for generating nanoparticle network using application " << Init->app_name.str << "." << endl;
        hout << "Only valid options are 'Create_Network', 'Read_Network' or 'Network_From_Seeds'. Input was:" << Init->simu_para.create_read_network << endl;
        return 0;
    }
    
    //----------------------------------------------------------------------
    //Vector for GNP shells
    vector<Shell> shell_gnps(gnps.size());
    ct0 = time(NULL);
    Shells *SH = new Shells;
    if (!SH->Generate_shells(Init->simu_para.particle_type, Init->geom_sample, points_cnt, gnps, shells_cnts, shell_gnps)) {
        hout << "Error when generating shells" << endl;
        return 0;
    }
    delete SH;
    ct1 = time(NULL);
    hout << "Generate shells and structure time: "<<(int)(ct1-ct0)<<" secs."<<endl;//*/
    
    for(int i=0; i<=Init->geom_sample.cut_num; i++)
    {
        hout << "============================================================================"<<endl;
        hout << "============================================================================"<<endl;
        hout << "Iteration " << i+1 << " of " << Init->geom_sample.cut_num+1 << endl;
        time_t it0, it1;
        it0 = time(NULL);
        
        //Update observation window geometry
        //hout<<"Update observation window geometry"<<endl;
        cuboid window_geo = Init->geom_sample.update_observation_window_geometry(i, Init->simu_para.particle_type);
        //hout<<"window_geo = "<<window_geo.str()<<endl;
        
        //----------------------------------------------------------------------
        //Extract the observation window
        Cutoff_Wins *Cutwins = new Cutoff_Wins;
        ct0 = time(NULL);
        //hout<<"Extract_observation_window"<<endl;
        if(!Cutwins->Extract_observation_window(i, Init->simu_para.particle_type, Init->geom_sample, window_geo, Init->nanotube_geo, gnps, structure_cnt, radii, points_cnt, shells_cnts, shell_gnps, structure_gnp, points_gnp)) {
            hout << "Error when extracting observation window "<< i+1 << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Extract observation window time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //----------------------------------------------------------------------
        //Calculate the CNTand/or GNP content in the observation window
        ct0 = time(NULL);
        if (!Calculate_content_in_window(window_geo, Init->simu_para.particle_type, points_cnt, radii, structure_cnt, Cutwins->cnts_inside, gnps,  Cutwins->gnps_inside)) {
            hout << "Error when calculating content in observation window "<< i+1 << endl;
            return 0;
        }
        ct1 = time(NULL);
        hout << "Calculate content in window time: "<<(int)(ct1-ct0)<<" secs."<<endl;
        
        //Delete object
        delete Cutwins;
        
        it1 = time(NULL);
        hout << "Iteration "<<i+1<<" time: "<<(int)(it1-it0)<<" secs."<<endl;
    }
    
    return 1;
}
//This function calculate the volume fraction of CNTs and/or GNPs on a given observation window
int App_Content_Dist::Calculate_content_in_window(const cuboid &window_geo, const string &particle_type, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &structure, const vector<int> &cnts_inside, const vector<GNP> &gnps, const vector<int> &gnps_inside)const
{
    //Calculate are of the face at the lower xy-plane
    double window_face_area = window_geo.len_x*window_geo.wid_y;
    
    //Calculate the volume of the window
    double window_vol = window_face_area*window_geo.hei_z;
    
    //Printer object to output to file
    Printer Pr;
    
    //Check if there are CNTs
    if (particle_type != "GNP_cuboids") {
        
        //Initialize the CNT volume at zero
        double cnt_vol = 0.0;
        
        //Calculate the volume of the CNTs inside the window
        for (size_t i = 0; i < cnts_inside.size(); i++) {
            
            //Get the CNT number
            int CNTi = cnts_inside[i];
            
            //Initialize variable to store the length of CNTo
            double cnti_len = 0.0;
            
            //Iterate over the points in the CNT
            for (size_t j = 0; j < structure[CNTi].size() - 1; j++) {
                
                //Get points j and k=j+1
                long int Pj = structure[CNTi][j];
                long int Pk = structure[CNTi][j+1];
                
                //Add the distance from Pj to Pk to the length of CNTi
                cnti_len = cnti_len + points_cnt[Pj].distance_to(points_cnt[Pk]);
            }
            
            //Calculate the volume of CNTi and add it to the total CNT volume
            cnt_vol = cnt_vol + cnti_len*PI*radii[CNTi]*radii[CNTi];
        }
        
        //Output to file depeneding if there is a deposit or not
        if (particle_type != "CNT_deposit") {
            
            //Output the volume fraction to the file for CNT volume fractions
            Pr.Append(cnt_vol/window_vol, "cnt_vol_fracs.txt");
        }
        else {
            
            //Output the volume per unit area to the file for CNT volume per unit area
            Pr.Append(cnt_vol/window_face_area, "cnt_vol_per_area.txt");
        }
    }
    
    //Check if there are GNPs
    if (particle_type != "CNT_wires" && particle_type != "CNT_deposit") {
        
        //Initialize the CNT volume at zero
        double gnp_vol = 0.0;
        
        //Calculate the volume of the GNPs inside the window
        for (size_t i = 0; i < gnps_inside.size(); i++) {
            
            //Get the GNP number
            int GNPi = gnps_inside[i];
            
            //Add the volume of GNPi to the total volume
            gnp_vol = gnp_vol + gnps[GNPi].volume;
        }
        
        //Output the volume fraction to the file for CNT volume fractions
        Pr.Append(gnp_vol/window_vol, "gnp_vol_fracs.txt");
    }
    
    return 1;
}
