//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Calculate electrical resistivity of a network read from a file (csv or binary)
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "App_Network_From_File.h"

//Generate 3D nanoparticle network, turn it into a resitor network and find its electrical conductivity
int App_Network_From_Abaqus::Generate_nanoparticle_resistor_network_from_file(Input* Init)const
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
    //Network Generation from file
    hout << "Generating nanoparticle network from file......" << endl;
    

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

    for (int i = 0; i <= Init->geom_sample.cut_num; i++)
    {
        hout << "============================================================================" << endl;
        hout << "============================================================================" << endl;
        hout << "Iteration " << i + 1 << " of " << Init->geom_sample.cut_num + 1 << endl;
        time_t it0, it1;
        it0 = time(NULL);

        //Update observation window geometry
        //hout<<"Update observation window geometry"<<endl;
        cuboid window_geo = Init->geom_sample.update_observation_window_geometry(i, Init->simu_para.particle_type);
        //hout<<"window_geo = "<<window_geo.str()<<endl;

        //Export the window geometry if needed
        if (Init->vis_flags.window_domain) {
            string str = "window_" + to_string(i) + ".vtk";
            VTK_Export VTK_E;
            VTK_E.Export_cuboid(window_geo, str);
        }

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
        hout << "Iteration " << i + 1 << " time: " << (int)(it1 - it0) << " secs." << endl;
    }

    return 1;
}