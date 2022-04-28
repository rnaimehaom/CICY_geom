//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Generate a grid to facilitate findic contact points
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Contact_grid.h"

/*
 
 This function groups the points inside the observation window into sub-regions.
 
 The task of finding overlapping points can be computationally expensive. To reduce computational cost, the sample was divided into smaller cubic sub-regions. The cubic sample of size $[a \times a \times a ]$ was divided into $m$ segments in each direction resulting in $m^3$ cubic sub-regions. Then, overlapping points are searched for only on each smaller sub-region rather than in the whole sample. It may happen that two overlapping points belong to different sub-regions. In order to take into account these overlapping points, each sub-region is ``extended". If a point lies exactly at a boundary, then an overlapping point can be at a maximum distance equal to $r_{max}+r_{max}+d_W = 2r_{max}+d_W$ from the boundary. Here, $r_{max}$ is the maximum radii of the CNTs. Then, each boundary plane of each cubic sub-region is translated externally to a parallel plane at a distance $2r_{max}+d_W$.
 
 */

int Contact_grid::Generate_contact_grid(const int &window, const string &particle_type, const Geom_sample &sample_geom, const cuboid &window_geom, const vector<int> &cnts_inside, vector<Point_3D> &points_cnt, const vector<vector<long int> > &structure, const vector<int> &gnps_inside, const vector<GNP> &gnps)
{
    
    //Number of subregions on each direction
    int sx = 1 + (int)(window_geom.len_x/sample_geom.gs_minx);
    int sy = 1 + (int)(window_geom.wid_y/sample_geom.gs_miny);
    int sz = 1 + (int)(window_geom.hei_z/sample_geom.gs_minz);
    int n_regions[] = {sx,sy,sz};
    //hout << "sx = " << sx << '\t' << "lx = " << window_geom.len_x << '\t' << "dx = " << sample_geom.gs_minx << "\n";
    //hout << "sy = " << sy << '\t' << "ly = " << window_geom.wid_y << '\t' << "dy = " << sample_geom.gs_miny << "\n";
    //hout << "sz = " << sz << '\t' << "lz = " << window_geom.hei_z << '\t' << "dz = " << sample_geom.gs_minz << "\n";
    
    //Variable to store the sizes of the subregions on each direction
    double l_regions[] = {sample_geom.gs_minx, sample_geom.gs_miny, sample_geom.gs_minz};
    
    //Calculate the largest overlapping between the one for CNTs and the one for GNPs
    double max_overlapping = max(sample_geom.gs_overlap_cnt, sample_geom.gs_overlap_gnp);
    //max_overlapping = max_overlapping * 1.05;
    //hout << "max_overlapping=" << max_overlapping << endl;
    
    //Check if the number of regions needs to be adjusted
    //hout << "Adjust_regions_if_needed" << endl;
    if (!Adjust_regions_if_needed(max_overlapping, sample_geom, window_geom, n_regions, l_regions)) {
        hout << "Error in Generate_contact_grid when calling Adjust_regions_if_needed" << endl;
        return 0;
    }

    //Calculate the total number of subregions
    int tot_regions = n_regions[0] * n_regions[1] * n_regions[2];
    if (tot_regions <= 0) {
        hout << "Error in Fill_sectioned_domain_gnps. Invalid number of subregions: " << tot_regions << endl;
        return 0;
    }

    //hout<<"There are "<< tot_regions <<" overlapping sub-regions."<<endl;
    //hout<<"\tn_regions[0]="<<n_regions[0]<<" n_regions[1]="<<n_regions[1]<<" n_regions[2]="<<n_regions[2]<<endl;
    //hout<<"\tl_regions[0]="<<l_regions[0]<<" l_regions[1]="<<l_regions[1]<<" l_regions[2]="<<l_regions[2]<<endl;
    
    //Fill the vector for the sectioned domain for CNTs only when there are CNTs in the structure
    //i.e., when the structure is not made of only GNPs
    if (particle_type != "GNP_cuboids") {
        
        //hout<<"Fill_sectioned_domain_cnts"<<endl;
        if (!Fill_sectioned_domain_cnts(window_geom, cnts_inside, structure, points_cnt, sample_geom.gs_overlap_cnt, n_regions, l_regions, tot_regions)) {
            hout << "Error in Generate_contact_grid when calling Generate_sectioned_domain_cnts" << endl;
            return 0;
        }
    }
    
    //Fill the vector for the sectioned domain for discrete GNPs only when there are GNPs in the structure
    //Equivalently, this is done when the structure is not made of only CNTs
    if (particle_type != "CNT_wires" && particle_type != "CNT_deposit") {
        
        //Fill the sectioned domain corresponding to GNP numbers
        //hout<<"Fill_sectioned_domain_gnps"<<endl;
        if (!Fill_sectioned_domain_gnps(window_geom, gnps, gnps_inside, sample_geom.gs_overlap_gnp, n_regions, l_regions, tot_regions)) {
            hout << "Error in Generate_contact_grid when calling Generate_sectioned_domain_cnts" << endl;
            return 0;
        }
    }
    
    return 1;
}
//This function checks if the subregions are too small, and if so, their size is adjusted
int Contact_grid::Adjust_regions_if_needed(const double &overlapping, const Geom_sample &sample_geom, const cuboid &window_geom, int n_regions[], double l_regions[])
{
    //Check that the regions are not too small for the maximum cutoff distance 2r_max+tunnel
    //If they are, then change the number of sections to the maximum possible
    if (sample_geom.gs_minx < 2*overlapping) {
        l_regions[0] = 2.0 * (overlapping + Zero);
        n_regions[0] = 1 + (int)(window_geom.len_x/l_regions[0]);
        hout << "Modified the number of sections along x. " << "sx = " << n_regions[0] << '\t' << "dx = " << l_regions[0] << endl;
    }
    if (sample_geom.gs_miny < 2*overlapping) {
        l_regions[1] = 2.0 * (overlapping + Zero);
        n_regions[1] = 1 + (int)(window_geom.wid_y/ l_regions[1]);
        hout << "Modified the number of sections along y. " << "sy = " << n_regions[1] << '\t' << "dy = " << l_regions[1] << endl;
    }
    if (sample_geom.gs_minz < 2*overlapping) {
        l_regions[2] = 2.0 * (overlapping + Zero);
        n_regions[2] = 1 + (int)(window_geom.hei_z/ l_regions[2]);
        hout << "Modified the number of sections along z. " << "sz = " << n_regions[2] << '\t' << "dz = " << l_regions[2] << endl;
    }
    return 1;
}
//This function fills the sectioned domain vector for CNTs
int Contact_grid::Fill_sectioned_domain_cnts(const cuboid &window_geom, const vector<int> &cnts_inside, const vector<vector<long int> > &structure, const vector<Point_3D> &points_cnt, const double &overlapping, const int n_regions[], const double l_regions[], const int& tot_regions)
{
    //There will be n_regions[0]*n_regions[1]*n_regions[2] different regions
    sectioned_domain_cnts.clear();
    sectioned_domain_cnts.assign(tot_regions, vector<long int>());
    //hout<<"window_geom="<<window_geom.str()<<endl;
    //hout<<"sectioned_domain_cnts.size="<<sectioned_domain_cnts.size()<<endl;
    
    //First loop over the particles inside the box, then loop over the points inside each particle
    for (int i = 0; i < (int)cnts_inside.size(); i++) {
        //hout<<"cnts_inside["<<i<<"]="<<cnts_inside[i]<<' ';
        
        //Current CNT number
        int CNT = cnts_inside[i];
        //hout<<"structure["<<CNT<<"].size()="<<structure[CNT].size()<<endl;
        
        //Scan each point in the particle
        for (int j = 0; j < (int)structure[CNT].size(); j++) {
            
            //Current point number
            long int P = structure[CNT][j];
            //hout<<"P=structure["<<CNT<<"]["<<j<<"]="<<structure[CNT][j]<< endl;
            
            //Calculate the region-coordinates
            int a, b, c;
            //hout << "points_cnt[P="<<P<<"]="<<points_cnt[P].str()<< endl;
            if (!Calculate_region_coordinates(window_geom, points_cnt[P], l_regions, a, b, c)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Calculate_region_coordinates" << endl;
                return 0;
            }
            //hout<<"a="<<a<<" b="<<b<<" c="<<c<<endl;
            
            //Initialize flags for overlaping regions
            int f_regions[] = {0,0,0};
            
            //Assign value of flag according to position of point
            //hout<<"Calculate_overlapping_flags"<<endl;
            if (!Calculate_overlapping_flags(window_geom, points_cnt[P], overlapping, a, b, c, n_regions, l_regions, f_regions)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Calculate_postion_flags" << endl;
                return 0;
            }
            
            //Assign the point P to the correspoding region or regions
            //hout<<"Assign_point_to_regions_cnts"<<endl;
            if (!Assign_point_to_regions_cnts(a, b, c, f_regions, n_regions, tot_regions, P)) {
                hout << "Error in Generate_sectioned_domain_cnts when calling Assign_point_to_region" << endl;
                hout << "P=" << points_cnt[P].str() << endl;
                hout << "CNT=" << CNT << " P#=" << P << endl;
                hout << "window_geom=" << window_geom.str() << endl;
                return 0;
            }
            //hout<<"for-j"<<endl;
        }
        //hout<<"for-i"<<endl;
    }
    
    return 1;
}
//This function calculates the subregion coordinates for a given point
int Contact_grid::Calculate_region_coordinates(const cuboid &window_geom, const Point_3D &point, const double l_regions[], int &a, int &b, int &c)
{
    
    //Calculate the region-coordinates
    a = (int)((point.x - window_geom.poi_min.x)/l_regions[0]);
    b = (int)((point.y - window_geom.poi_min.y)/l_regions[1]);
    c = (int)((point.z - window_geom.poi_min.z)/l_regions[2]);
    
    return 1;
}
//
int Contact_grid::Calculate_overlapping_flags(const cuboid &window_geom, const Point_3D &point, const double &overlapping, const int &a, const int &b, const int &c, const int n_regions[], const double l_regions[], int f_regions[])
{
    //Coordinates of non-overlaping region the point belongs to
    double x1 = (double)a*l_regions[0] + window_geom.poi_min.x;
    double x2 = x1 + l_regions[0];
    double y1 = (double)b*l_regions[1] +  window_geom.poi_min.y;
    double y2 = y1 + l_regions[1];
    double z1 = (double)c*l_regions[2] +  window_geom.poi_min.z;
    double z2 = z1 + l_regions[2];
    
    //Assign value of flag according to position of point
    //The first operand eliminates the periodicity on the boundary
    if ((a > 0) && (point.x >= x1) && (point.x <= x1+overlapping))
        f_regions[0] = -1;
    else if ((a < n_regions[0]-1) && (point.x >= x2-overlapping) && (point.x <= x2 ))
        f_regions[0] = 1;
    if ((b > 0) && (point.y >= y1) && (point.y <= y1+overlapping))
        f_regions[1] = -1;
    else if ((b < n_regions[1]-1) && (point.y >= y2-overlapping) && (point.y <= y2 ))
        f_regions[1] = 1;
    if ((c > 0) && (point.z >= z1) && (point.z <= z1+overlapping))
        f_regions[2] = -1;
    else if ((c < n_regions[2]-1) && (point.z >= z2-overlapping) && (point.z <= z2 ))
        f_regions[2] = 1;
    
    return 1;
}
//
int Contact_grid::Assign_point_to_regions_cnts(const int& a, const int& b, const int& c, const int f_regions[], const int n_regions[], const int& tot_regions, const long int& P)
{
    //Calculate default rregion
    int t = Calculate_t(a, b, c, n_regions[0], n_regions[1]);

    //Assign point to default region
    //hout<<"t="<<t<<endl;
    if (t >= tot_regions || t < 0)
    {
        hout << "Error in Assign_point_to_regions_cnts. Point belongs to a subregion outside the range of subregions." << endl;
        hout << "Subregion number (t) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
        hout << "a=" << a << " b=" << b << " c=" << c << endl;
        return 0;
    }
    sectioned_domain_cnts[t].push_back(P);

    //Check if P needs to be added to more regions due to their overlapping
    //Check the flag for the subregion's a-coordinate
    if (f_regions[0] != 0) {

        //Calculate the region
        t = Calculate_t(a + f_regions[0], b, c, n_regions[0], n_regions[1]);
        //hout<<"\tt0="<<t<<endl;
        if (t >= tot_regions || t < 0)
        {
            hout << "Error in Assign_point_to_regions_cnts. Point belongs to a subregion outside the range of subregions." << endl;
            hout << "Subregion number (t0) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
            hout << "a=" << a << " b=" << b << " c=" << c << endl;
            hout << "f_regions[0]=" << f_regions[0] << " f_regions[1]=" << f_regions[1] << " f_regions[2]=" << f_regions[2] << endl;
            return 0;
        }

        //Assign point to region
        sectioned_domain_cnts[t].push_back(P);

    }
    //Check the flag for the subregion's b-coordinate
    if (f_regions[1] != 0) {

        //Calculate the region
        //hout<<"\tt1="<<t<<endl;
        t = Calculate_t(a, b + f_regions[1], c, n_regions[0], n_regions[1]);
        if (t >= tot_regions || t < 0)
        {
            hout << "Error in Assign_point_to_regions_cnts. Point belongs to a subregion outside the range of subregions." << endl;
            hout << "Subregion number (t1) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
            hout << "a=" << a << " b=" << b << " c=" << c << endl;
            hout << "f_regions[0]=" << f_regions[0] << " f_regions[1]=" << f_regions[1] << " f_regions[2]=" << f_regions[2] << endl;
            return 0;
        }

        //Assign point to region
        sectioned_domain_cnts[t].push_back(P);

    }
    //Check the flag for the subregion's c-coordinate
    if (f_regions[2] != 0) {

        //Calculate the region
        //hout<<"\tt2="<<t<<endl;
        t = Calculate_t(a, b, c + f_regions[2], n_regions[0], n_regions[1]);
        if (t >= tot_regions || t < 0)
        {
            hout << "Error in Assign_point_to_regions_cnts. Point belongs to a subregion outside the range of subregions." << endl;
            hout << "Subregion number (t2) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
            hout << "a=" << a << " b=" << b << " c=" << c << endl;
            hout << "f_regions[0]=" << f_regions[0] << " f_regions[1]=" << f_regions[1] << " f_regions[2]=" << f_regions[2] << endl;
            return 0;
        }

        //Assign point to region
        sectioned_domain_cnts[t].push_back(P);

    }

    return 1;
}
//Calculates the region to which a point corresponds
int Contact_grid::Calculate_t(const int &a, const int &b, const int &c, const int &sx, const int &sy)
{
    return a + b*sx + c*sx*sy;
}
//
int Contact_grid::Fill_sectioned_domain_gnps(const cuboid &window_geom, const vector<GNP> &gnps, const vector<int> &gnps_inside, const double &overlapping, const int n_regions[], const double l_regions[], const int& tot_regions)
{
    //Initialize the sectioned domain with the number of regions
    sectioned_domain_gnps.clear();
    sectioned_domain_gnps.assign(tot_regions, vector<int>());
    //hout << "sectioned_domain" << endl;
    
    //First loop over the GNPs inside the sample, then loop over the points inside each GNP
    for (int i = 0; i < (int)gnps_inside.size(); i++) {
        
        //Get the GNP number
        //hout<<"gnps_inside[i="<<i<<"]="<<gnps_inside[i]<<endl;
        int GNPi = gnps_inside[i];
        
        //Discretize the GNP and add it to all the corresponding subregions
        //hout << "Fill_sectioned_domain_single_gnp" << endl;
        if (!Fill_sectioned_domain_single_gnp(window_geom, gnps[GNPi], overlapping, n_regions, l_regions, tot_regions)) {
            hout<<"Error in Fill_sectioned_domain_gnps when calling Fill_sectioned_domain_single_gnp"<<endl;
            return 0;
        }
        
    }
    
    return 1;
}
//This function finds all the subregions the GNP_new occupies and those subregions where there might be
//close enough GNPs below the van der Waals distance
int Contact_grid::Fill_sectioned_domain_single_gnp(const cuboid &window_geom, const GNP &gnp, const double &overlapping, const int n_regions[], const double l_regions[], const int& tot_regions)
{
    //Number of points to discretize the GNP along the x direction (of the GNP local coordinates)
    int n_points_x = max(20, 1 + (int)(gnp.l/l_regions[0]));
    
    //Number of points to discretize the GNP along the y direction (of the GNP local coordinates)
    int n_points_y = max(20, 1 + (int)(gnp.l/l_regions[1]));
    
    //Number of points to discretize the GNP along the z direction (of the GNP local coordinates)
    //Make sure there are at least two points in the discretization along z
    int n_points_z = max(2, 1 + (int)(gnp.t/l_regions[2]));

    //Calculate
    
    //Iterate over all points in the discretization
    //Index k moves the point in the discretization along z
    for (int k = 0; k < n_points_z; k++) {
        
        //Calculate the start and end points along the y-direction
        //hout<<"k="<<k<<endl;
        
        //Their locations are proportional to index k and the number of points
        //in the discretization along z
        double lambda_z = (double)k/(n_points_z-1);
        
        //These points move along edges 3-7 and 0-4 of the GNP
        Point_3D start1_y = gnp.vertices[7] + (gnp.vertices[3] - gnp.vertices[7])*lambda_z;
        Point_3D end1_y = gnp.vertices[4] + (gnp.vertices[0] - gnp.vertices[4])*lambda_z;
        Point_3D diff1_y = end1_y - start1_y;
        
        //These points move along edges 2-6 and 1-5 of the GNP
        Point_3D start2_y = gnp.vertices[6] + (gnp.vertices[2] - gnp.vertices[6])*lambda_z;
        Point_3D end2_y = gnp.vertices[5] + (gnp.vertices[1] - gnp.vertices[5])*lambda_z;
        Point_3D diff2_y = end2_y - start2_y;
        
        //Index j moves the point in the discretization along y
        for (int j = 0; j < n_points_y; j++) {
            
            //Calculate the start and end points along the x-direction
            //hout<<"**j="<<j<<endl;
            
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
                //hout<<"****i="<<i<<endl;
                
                //The location of new_point is proportional to index i and the number of points
                //in the discretization along x
                double lambda_x = (double)i/(n_points_x-1);
                
                //Location of new_point
                Point_3D point = start_x + diff_x*lambda_x;

                //Check if new_point is inside the sample
                if( !(point.x<=window_geom.poi_min.x||point.x>=window_geom.max_x||
                      point.y<=window_geom.poi_min.y||point.y>=window_geom.max_y||
                      point.z<=window_geom.poi_min.z||point.z>=window_geom.max_z) ) {
                    
                    //If the point is inside the window, then calcualte the region coordinates of the point
                    int a, b, c;
                    //hout<<"Calculate_region_coordinates"<<endl;
                    if (!Calculate_region_coordinates(window_geom, point, l_regions, a, b, c)) {
                        hout<<"Error in Fill_sectioned_domain_single_gnp when calling Calculate_region_coordinates"<<endl;
                        return 0;
                    }
                    //hout<<"a="<<a<<" b="<<b<<" c="<<c<<endl;
                    
                    //Initialize flags for overlaping regions
                    int f_regions[] = {0,0,0};
                    
                    //hout<<"Calculate_overlapping_flags"<<endl;
                    if (!Calculate_overlapping_flags(window_geom, point, overlapping, a, b, c, n_regions, l_regions, f_regions)) {
                        hout<<"Error in Fill_sectioned_domain_single_gnp when calling Calculate_overlapping_flags"<<endl;
                        return 0;
                    }
                    
                    //Add to the subregion(s) the GNP point occupies
                    //hout<<"Assign_point_to_regions_gnps point="<<point.str()<<endl;
                    if (!Assign_point_to_regions_gnps(a, b, c, f_regions, n_regions, tot_regions, gnp.flag)) {
                        hout<<"Error in Fill_sectioned_domain_single_gnp when calling Assign_point_to_regions_gnps"<<endl;
                        return 0;
                    }
                }
            }
        }
    }
    
    return 1;
}
//
int Contact_grid::Assign_point_to_regions_gnps(const int &a, const int &b, const int &c, const int f_regions[], const int n_regions[], const int& tot_regions, const int &gnp_i)
{
    //Calculate default region
    int t = Calculate_t(a, b, c, n_regions[0], n_regions[1]);
    if (t >= tot_regions || t < 0)
    {
        hout << "Error in Assign_point_to_regions_gnps. Point belongs to a subregion outside the range of subregions." << endl;
        hout << "Subregion number (t) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
        hout << "GNPi=" << gnp_i << endl;
        hout << "a=" << a << " b=" << b << " c=" << c << endl;
        return 0;
    }
    
    //Assign GNP number of point to default region if not already added
    if (sectioned_domain_gnps[t].empty() || sectioned_domain_gnps[t].back() != gnp_i) {
        //hout << "GNP " << gnp_i << " in region " << t << endl;
        sectioned_domain_gnps[t].push_back(gnp_i);
    }
    
    //Check if P needs to be added to more regions due to their overlapping
    //Check the flag for the subregion's a-coordinate
    if (f_regions[0] != 0) {
        
        //Calculate the region
        t = Calculate_t(a+f_regions[0], b, c, n_regions[0], n_regions[1]);
        if (t >= tot_regions || t < 0)
        {
            hout << "Error in Assign_point_to_regions_gnps. Point belongs to a subregion outside the range of subregions." << endl;
            hout << "Subregion number (t0) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
            hout << "GNPi=" << gnp_i << endl;
            hout << "a=" << a << " b=" << b << " c=" << c << endl;
            hout << "f_regions[0]=" << f_regions[0] << " f_regions[1]=" << f_regions[1] << " f_regions[2]=" << f_regions[2] << endl;
            return 0;
        }
        
        //Assign GNP number of point to region if not already added
        if (sectioned_domain_gnps[t].empty() || sectioned_domain_gnps[t].back() != gnp_i) {
            //hout << "GNP " << gnp_i << " also in region " << t << endl;
            sectioned_domain_gnps[t].push_back(gnp_i);
        }
        
    }
    //Check the flag for the subregion's b-coordinate
    if (f_regions[1] != 0) {
        
        //Calculate the region
        t = Calculate_t(a, b+f_regions[1], c, n_regions[0], n_regions[1]);
        if (t >= tot_regions || t < 0)
        {
            hout << "Error in Assign_point_to_regions_gnps. Point belongs to a subregion outside the range of subregions." << endl;
            hout << "Subregion number (t1) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
            hout << "GNPi=" << gnp_i << endl;
            hout << "a=" << a << " b=" << b << " c=" << c << endl;
            hout << "f_regions[0]=" << f_regions[0] << " f_regions[1]=" << f_regions[1] << " f_regions[2]=" << f_regions[2] << endl;
            return 0;
        }
        
        //Assign GNP number of point to region if not already added
        if (sectioned_domain_gnps[t].empty() || sectioned_domain_gnps[t].back() != gnp_i) {
            //hout << "GNP " << gnp_i << " also in region " << t << endl;
            sectioned_domain_gnps[t].push_back(gnp_i);
        }
        
    }
    //Check the flag for the subregion's c-coordinate
    if (f_regions[2] != 0) {
        
        //Calculate the region
        t = Calculate_t(a, b, c+f_regions[2], n_regions[0], n_regions[1]);
        if (t >= tot_regions || t < 0)
        {
            hout << "Error in Assign_point_to_regions_gnps. Point belongs to a subregion outside the range of subregions." << endl;
            hout << "Subregion number (t2) is " << t << ". Maximum number of subregions is " << tot_regions << endl;
            hout << "GNPi=" << gnp_i << endl;
            hout << "a=" << a << " b=" << b << " c=" << c << endl;
            hout << "f_regions[0]=" << f_regions[0] << " f_regions[1]=" << f_regions[1] << " f_regions[2]=" << f_regions[2] << endl;
            return 0;
        }
        
        //Assign GNP number of point to region if not already added
        if (sectioned_domain_gnps[t].empty() || sectioned_domain_gnps[t].back() != gnp_i) {
            //hout << "GNP " << gnp_i << " also in region " << t << endl;
            sectioned_domain_gnps[t].push_back(gnp_i);
        }
        
    }
    
    return 1;
}
