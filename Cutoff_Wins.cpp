//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Cut out an observation window
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#include "Cutoff_Wins.h"

//This function removes the points that are outside the observation window.
//The vector cnts_inside is created so only the CNTs inside the obseration window are considered in other functions
int Cutoff_Wins::Extract_observation_window(const int &window, const string &particle_type, const Geom_sample &sample_geo, const cuboid &window_geo, const Nanotube_Geo &cnts_geo, vector<GNP> &gnps, vector<vector<long int> > &structure_cnt, vector<double> &radii, vector<Point_3D> &points_cnt, vector<vector<int> > &shells_cnt, vector<Shell> &shells_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp)
{
    
    //Vector to save initial seeds
    vector<long int> seeds;
    //Save the initial points of the CNTs that are attached to the GNP
    /* /Check if particle type is the hybrid
    //hout << "Save seeds" << endl;
    if (particle_type == "Hybrid_particles") {
        if (!Save_seeds(hybrid_particles, structure, seeds)) {
            hout << "Error in Extract_observation_window when calling Save_seeds" << endl;
            return 0;
        }
    }*/
    
    //Trim the CNTs if there are CNTs in the structure
    //Check if the generated structure has CNTs, this happens when the particle type is not GNPs
    if (particle_type != "GNP_cuboids") {

        //hout<<"Trim_boundary_cnts"<<endl;
        if (!Trim_boundary_cnts(window, particle_type, sample_geo, window_geo, cnts_geo, points_cnt, structure_cnt, shells_cnt, radii)) {
            hout << "Error in Extract_observation_window when calling Trim_boundary_cnts" << endl;
            return 0;
        }
        
        //hout<<"Fill_cnts_inside"<<endl;
        //Fill the vector cnts_inside
        if (!Fill_cnts_inside(structure_cnt)) {
            hout << "Error in Extract_observation_window when calling Fill_cnts_inside" << endl;
            return 0;
        }
    }
    
    //Flag to determine if hybrid particles are used
    //This flag determines if GNP contacts with the boundary are ignored
    int hybrid_flag = 0;
    
    //Update the CNTs attached to the GNP
    //Check if particle type is hybrid
    if (particle_type == "Hybrid_particles") {
        
        //hout << "Compare_seeds" << endl;
        //Compare the initial points of the CNTs attached to the GNPs
        /* /If they are different, that means that the CNT is not attached to the GNP anymore
        if (!Compare_seeds(hybrid_particles, structure, seeds)) {
            hout << "Error in Extract_observation_window when calling Compare_seeds" << endl;
            return 0;
        }*/
        
        //Set the flag to 1 as hybrid particles are being used
        hybrid_flag = 1;
    }
    
    //Remove GNPs that are outside the observation window
    //Check if the structure has GNPs, this happens when the particle type is not CNT
    if (particle_type != "CNT_wires" && particle_type != "CNT_deposit") {
        
        //Assign the correct size to the GNP structure
        structure_gnp.assign(gnps.size(), vector<long int>());
    
        //Fill the vector gnps_inside
        //hout<<"Fill_gnps_inside"<<endl;
        if (!Fill_gnps_inside(window, window_geo, shells_gnp, gnps, structure_gnp, points_gnp)) {
            hout << "Error in Extract_observation_window when calling Fill_gnps_inside" << endl;
            return 0;
        }
    }
    
    return 1;
}
int Cutoff_Wins::Trim_boundary_cnts(const int &window, const string &particle_type, const Geom_sample &sample_geo, const cuboid &window_geo, const Nanotube_Geo &cnts, vector<Point_3D> &points_cnt, vector<vector<long int> > &structure_cnt, vector<vector<int> > &shells_cnt, vector<double> &radii)
{
    //String to save the location of a point (inside the window, outside the window, or at a boundary)
    string point_location;

    //Cutoff for the step length
    double step_cutoff = cnts.step_length * 0.5;
    
    //Initialize the vectors of CNTs and CNT points at boundaries
    boundary_cnt.assign(6, vector<int>());
    boundary_cnt_pts.assign(6, vector<long int>());
    
    //Variable to reduce computations when adding CNTs to a shell
    //vars_shells[0] = midpoints
    //vars_shells[1] = boundary_layer
    //vars_shells[2] = core
    //vars_shells[3] = half_step
    double vars_shells[4][3] = {
        {sample_geo.sample.poi_min.x+sample_geo.sample.len_x/2.0,
        sample_geo.sample.poi_min.y+sample_geo.sample.wid_y/2.0,
        sample_geo.sample.poi_min.z+sample_geo.sample.hei_z/2.0},
        {sample_geo.sample.poi_min.x+(sample_geo.sample.len_x-sample_geo.win_max_x)/2.0,
        sample_geo.sample.poi_min.y+(sample_geo.sample.wid_y-sample_geo.win_max_y)/2.0,
        sample_geo.sample.poi_min.z+(sample_geo.sample.hei_z-sample_geo.win_max_z)/2.0},
        {sample_geo.sample.poi_min.x+(sample_geo.sample.len_x-sample_geo.win_min_x)/2.0,
        sample_geo.sample.poi_min.y+(sample_geo.sample.wid_y-sample_geo.win_min_y)/2.0,
        sample_geo.sample.poi_min.z+(sample_geo.sample.hei_z-sample_geo.win_min_z)/2.0},
        {sample_geo.win_delt_x/2.0, sample_geo.win_delt_y/2.0, sample_geo.win_delt_z/2.0}
    };
    
    //Scan all CNTs in the current shell
    for (long int i = 0; i < (long int)shells_cnt[window].size(); i++) {
        
        //Get the current CNT
        int CNT = shells_cnt[window][i];
        //hout << "CNT=" << CNT << " shell=" << window << " i=" << i << endl;
        
        //Generate the percolation layer cuboid
        cuboid layer_geom;
        /*if (!Get_percolation_layer_cuboid(radii[CNT], window_geo, layer_geom)) {
            hout<<"Error in Trim_boundary_cnts when calling Get_percolation_layer_cuboid"<<endl;
            return 0;
        }*/
        
        //Indices to define the beginning and ending of a segment of a CNT that is inside a sample
        int start = 0;
        int end = 0;
        
        //Variables to store the index of the first and last point in the first segement
        int first_idx = 0;
        int last_idx = 0;
        
        //Variable to store the number of CNT segments of the current CNT
        int segments = 0;
        
        //Index of the last point inside the sample
        int last_inside = 0;
        
        //Number of points in the current CNT
        int cnt_points = (int)structure_cnt[CNT].size();
        //hout << "cnt_points=" << cnt_points << endl;
        
        //Check if there are points in the CNT
        //Because there a minimum number of points may be specified to determine
        //if a CNT is inored, a CNT may be in a shell but have size 0
        if (cnt_points)
        {

            //Initialize the point number with the first point in the CNT
            long int P1 = structure_cnt[CNT][0];

            //Check where is the first point of the CNT
            //flag = 1 because this point is an endpoint
            //Calculate vector V
            //Point_3D V = (points_cnt[P1] - points_cnt[structure_cnt[CNT][1]]).unit();
            //string is_first_inside_sample = Where_is_with_layer(points_cnt[P1], window_geo, layer_geom, 1, radii[CNT], V);
            string is_first_inside_sample = Where_is(points_cnt[P1], window_geo);
            //hout << "P0=" << points_cnt[P1].str() << endl;
            //hout<<"P0="<<P1<<" idx=0 CNT="<<CNT<<" loc="<<is_first_inside_sample<<endl;

            //Scan all remaning points in the current CNT
            for (int j = 1; j < cnt_points; j++) {

                //Get the current point number
                P1 = structure_cnt[CNT][j];

                //Check if this is the last point
                /*if (j == cnt_points - 1) {
                    //It is the last point
                    //Recalculate vector V
                    long int P2 = structure_cnt[CNT][j-1];
                    V = (points_cnt[P1] - points_cnt[P2]).unit();
                    //flag = 1 because this point is an endpoint
                    point_location = Where_is_with_layer(points_cnt[P1], window_geo, layer_geom, 1, radii[CNT], V);
                }
                else {
                    point_location = Where_is_with_layer(points_cnt[P1], window_geo, layer_geom);
                }*/
                point_location = Where_is(points_cnt[P1], window_geo);
                //hout << "P1=" << points_cnt[P1].str() << endl;
                //hout << "P1=" << P1 << " idx=" << j << " CNT=" << CNT << " loc=" << point_location << endl;

                //Check if the point is inside the window
                if (point_location != "outside") {

                    //Update the last inside point
                    last_inside = j;

                }
                //The point is outside the window or at a boundary, then a segment might be added
                else {

                    //End index is the current looping index
                    end = j;

                    //Count the number of consecutive points inside the sample
                    //If start is zero and it is inside the sample or at the boundary, then the number of
                    //points in the CNT segment is end - start because end is outside (always in this for-loop)
                    //If start is zero and it is outside the sample, then the number of points in
                    //the CNT segment is end - start - 1 because both end and start are outside
                    //If start is not zero, then it is always outside, and since end is also outside, then
                    //the number of points in the CNT segment is end - start - 1
                    //Thus, the number of points in the CNT segment is end - start only when start
                    //is zero and not outside the sample.
                    //Otherwise the number of points is end - start - 1
                    int n_points = (start == 0 && is_first_inside_sample != "outside") ? end - start : end - start - 1;

                    //Check if there are enough points to consider this a whole CNT and include it in the analysis
                    if (n_points >= cnts.min_points) {

                        //Add the current segment to the structure
                        //hout << "Add_cnt_segment_to_structure loop" << endl;
                        if (!Add_cnt_segment_to_structure(particle_type, step_cutoff, window_geo, layer_geom, vars_shells, start, end, cnts.min_points, CNT, point_location, points_cnt, structure_cnt, shells_cnt, radii, segments, first_idx, last_idx)) {
                            hout << "Error when adding a CNT segment (Add_cnt_segment_to_structure 1)." << endl;
                            return 0;
                        }
                    }

                    //Reset the start index
                    start = j;
                }
            }

            //Check if the last index of the vector structure[CNT] was inside the sample
            //hout<<"last_inside="<<last_inside<<" cnt_points="<<cnt_points<<endl;
            //Also check that the segment has the minimum number of CNT points required
            //for it to be consedered a CNT
            if (last_inside == cnt_points - 1 && (last_inside - start + 1) >= cnts.min_points) {

                //Set end index as the last valid index
                end = cnt_points - 1;

                //If the last point of the CNT was inside the sample, then add a new segment
                //This was not done becuase, in the for loop, a segement is added only when it finds a point
                //outside the sample
                //Add the current segment to the structure
                //hout << "Add_cnt_segment_to_structure" << endl;
                if (!Add_cnt_segment_to_structure(particle_type, step_cutoff, window_geo, layer_geom, vars_shells, start, end, cnts.min_points, CNT, point_location, points_cnt, structure_cnt, shells_cnt, radii, segments, first_idx, last_idx)) {
                    hout << "Error when adding a CNT segment (Add_cnt_segment_to_structure 2)." << endl;
                    return 0;
                }
            }

            //Check if there is at least one segment
            //hout<<"Check if there is at least one segment"<<endl;
            if (segments > 0) {

                //Move the points to the front of the CNT if the start index of the first segment is not zero
                if (first_idx != 0) {
                    for (int k = first_idx; k <= last_idx; k++) {
                        structure_cnt[CNT][k - first_idx] = structure_cnt[CNT][k];
                    }

                    //Update the last idx
                    last_idx = last_idx - first_idx;
                }

                //If there were multiple segments, then remove the points from the current CNT
                //that are after the last index of the first segment and now belog to other CNTs
                for (int k = last_idx + 1; k < cnt_points; k++) {
                    structure_cnt[CNT].pop_back();
                }
            }
            //If there are no segments clear the points in the structure so that the CNT is not included
            //in the vector of CNTs inside
            else if (!segments) {
                structure_cnt[CNT].clear();
            }
            
        }
        
    }
    
    return 1;
}
//This function calculates the geometry for the percolation layer cuboid
int Cutoff_Wins::Get_percolation_layer_cuboid(const double &cnt_rad, const cuboid &window_geo, cuboid &layer_geom)
{
    //For the corner of the window, increase each coordinate a distance equal to the CNT radius
    layer_geom.poi_min = window_geo.poi_min + Point_3D(cnt_rad, cnt_rad, cnt_rad);
    
    //Reduce the maximum coordinates of the window by a distance equal to the CNT radius
    layer_geom.max_x = window_geo.max_x - cnt_rad;
    layer_geom.max_y = window_geo.max_y - cnt_rad;
    layer_geom.max_z = window_geo.max_z - cnt_rad;
    
    //For completeness reduce the length of each side of the window cuboid a distance equal to
    //twice the CNT radius (i.e., the CNT diameter)
    double cnt_diam = 2*cnt_rad;
    layer_geom.len_x = window_geo.len_x - cnt_diam;
    layer_geom.wid_y = window_geo.wid_y - cnt_diam;
    layer_geom.hei_z = window_geo.hei_z - cnt_diam;
    
    return 1;
}
//===========================================================================
int Cutoff_Wins::Add_cnt_segment_to_structure(const string &particle_type, const double& step_cutoff, const cuboid &window_geo, const cuboid &layer_geom, const double var_shells[][3], const int &start, const int &end, const int &min_points, const int &CNT, const string &end_point_loc, vector<Point_3D> &points_cnt, vector<vector<long int> > &structure_cnt, vector<vector<int> > &shells_cnt, vector<double> &radii, int &segments, int &first_idx, int &last_idx)
{
    //Variables used in case the start and end indices need to change
    int new_start = start;
    int new_end = end;
    
    //Variables for the (possibly) outside and inside points for the start of the segment
    long int p_out_start = structure_cnt[CNT][new_start];
    long int p_ins_start = structure_cnt[CNT][new_start+1];
    //hout<<"Start=("<<points_cnt[p_out_start].x<<", "<<points_cnt[p_out_start].y<<", "<<points_cnt[p_out_start].z<<") CNT="<<CNT<<endl;
    
    //Double check where is the start point of the current segment, if the first point is:
    //inside: do not add it to the boundary vectors
    //boundary: add it to the boundary vectors
    //outside: calculate the boundary point, subtitute the outside point by the boundary point,
    //          then add it to the boundary vectors
    string start_point_loc;
   /* if (p_out_start == structure_cnt[CNT][0]) {
        //p_out_start is the first point of the CNT
        //Recalculate vector V
        long int P2 = structure_cnt[CNT][1];
        Point_3D V = (points_cnt[p_out_start] - points_cnt[P2]).unit();
        //flag = 1 because this point is an endpoint
        start_point_loc = Where_is_with_layer(points_cnt[p_out_start], window_geo, layer_geom, 1, radii[CNT], V);
    }
    else {
        start_point_loc = Where_is_with_layer(points_cnt[p_out_start], window_geo, layer_geom);
    }*/
    start_point_loc = Where_is(points_cnt[p_out_start], window_geo);
    if (start_point_loc != "inside") {
        
        if (start_point_loc == "outside") {
            
            //Subtitute the point
            //The point at the boundary will be stored in points_cnt[p_out_start]
            if (!Substitute_boundary_point(window_geo, points_cnt[p_ins_start], points_cnt[p_out_start])) {
                hout<<"Error when substituting boundary point (start)"<<endl;
                return 0;
            }
        }
        
        //Choose the CNT to be added
        int CNT_num = (segments)? (int)structure_cnt.size(): CNT;

        //Check if the length of the segment at the boundary is below the step cutoff
        if (points_cnt[p_out_start].distance_to(points_cnt[p_out_start + 1]) - step_cutoff < Zero) {
            //hout << "d_start=" << points_cnt[p_out_start].distance_to(points_cnt[p_out_start + 1]) << endl;
            //The segment at the boundary is below the step cutoff,
            //so substitute this point by the following point
            //recall the boundary point is stored at points_cnt[p_out_start]
            points_cnt[p_out_start + 1] = points_cnt[p_out_start];
            //Update the start of the CNT
            new_start++;

            //Update the index p_out_start that contains the boundary point
            p_out_start++;
            //hout << "d_start_new=" << points_cnt[p_out_start].distance_to(points_cnt[p_out_start + 1]) << endl<<endl;
        }
        
        //Add to the boundary vectors, this happens when the first point is either outside or at a boundary
        if (!Add_cnt_point_to_boundary_vectors(window_geo, points_cnt[p_out_start], p_out_start, CNT_num)) {
            hout<<"Error in Add_cnt_segment_to_structure when caling Add_cnt_point_to_boundary_vectors (1)"<<endl;
            return 0;
        }
    }
    
    //Variables for the (possibly) outside and inside points for the end of the segment
    long int p_out_end = structure_cnt[CNT][new_end];
    long int p_ins_end = structure_cnt[CNT][new_end - 1];
    //hout<<"End=("<<points_cnt[p_out_end].x<<", "<<points_cnt[p_out_end].y<<", "<<points_cnt[p_out_end].z<<") CNT="<<CNT<<endl;
    
    //Double check where is the end point of the current segment, if the end point is:
    //inside: do not add it to the boundary vectors
    //boundary: add it to the boundary vectors
    //outside: calculate the boundary point, subtitute the outside point by the boundary point,
    //          then add it to the boundary vectors
    if (end_point_loc != "inside") {
        
        if (end_point_loc == "outside") {
            
            //Subtitute the point
            if (!Substitute_boundary_point(window_geo, points_cnt[p_ins_end], points_cnt[p_out_end])) {
                hout<<"Error when substituting boundary point (end)"<<endl;
                return 0;
            }
        }
        
        //Choose the CNT to be added
        int CNT_num = (segments)? (int)structure_cnt.size(): CNT;

        //Check if the length of the segment at the boundary is below the step cutoff
        if (points_cnt[p_out_end - 1].distance_to(points_cnt[p_out_end]) - step_cutoff < Zero) {
            //hout << "d_end=" << points_cnt[p_out_end - 1].distance_to(points_cnt[p_out_end]) << endl;
            //The segment at the boundary is below the step cutoff,
            //so substitute this point by the previous point
            //recall the boundary point is stored at points_cnt[p_out_end]
            points_cnt[p_out_end - 1] = points_cnt[p_out_end];

            //Update the end of the CNT
            new_end--;

            //Update the index p_out_end that contains the boundary point
            p_out_end--; 
            //hout << "d_end_new=" << points_cnt[p_out_end - 1].distance_to(points_cnt[p_out_end]) << endl<<endl;
        }
        
        //Add to the boundary vectors, this happens when the last point is either outside or at a boundary
        if (!Add_cnt_point_to_boundary_vectors(window_geo, points_cnt[p_out_end], p_out_end, CNT_num)) {
            hout<<"Error in Add_cnt_segment_to_structure when caling Add_cnt_point_to_boundary_vectors (2)"<<endl;
            return 0;
        }
    }
    
    //At this point all points in the segment are inside the observation window (or its boundary)
    //Check if this is the first segment
    if (segments == 0) {
        
        //This is the first segment, so save the two indices of the first segment
        first_idx = new_start;
        last_idx = new_end;
        //hout<<"segment 1, new_start="<<new_start<<" end="<<end<<endl;
    }
    else {
        
        //A new CNT is added, so get the new CNT number
        int new_CNT = (int)structure_cnt.size();
        //hout<<"New segment added CNT="<<CNT<<" new CNT="<<new_CNT<<endl;
        //hout<<"New segment added, new_start="<<new_start<<" end="<<end<<endl;
        
        //If this is not the first segment, then a new CNT needs to be added to the structure
        
        //Temporary vector to add a new CNT to the structure
        vector<long int> struct_temp;
        
        //This bg variable is used to add the new CNT into the corresponding shell or shells
        Shells shells = Shells();
        
        //Add the CNT points of the new segment to the 1D vector
        //After dealing with the boundary points, ALL CNTs are inside the shell, and thus
        //ALL points are added to the structure
        for(int j = new_start; j <= new_end; j++) {
            
            //Get the current point number
            long int P = structure_cnt[CNT][j];
            
            //Change the flag of current point to be that of the new CNT number
            points_cnt[P].flag = new_CNT;
            
            //Add the point number to the structure vector
            struct_temp.push_back(P);
            
            //Update the shell of the new CNT points
            shells.Add_to_cnt_shells(particle_type, var_shells[0], var_shells[1], var_shells[2], var_shells[3], points_cnt[P], (int)shells_cnt.size(), shells_cnt);
        }
        
        //Update the radii vector
        //The new CNT is just a segment of the old one, so they should have the same radius
        radii.push_back(radii[CNT]);
        
        //Add the new CNT to the structure
        structure_cnt.push_back(struct_temp);
    }
    
    //Update the number of segments
    //hout<<"segments++"<<endl;
    segments++;
    
    return 1;
}
int Cutoff_Wins::Substitute_boundary_point(const cuboid &window_geo, const Point_3D &p_inside, Point_3D &p_outside)
{
    //The line segment defined by p_outside and p_inside is given by:
    //P = p_outside + lambda*T
    //where T = p_inside - p_outside
    //In this way P = p_outside when lambda = 0 and P = p_inside when lambda = 1
    
    //Variable to store the point T = p_inside - p_outside
    Point_3D T = p_inside - p_outside;
    
    //Lambda function to calculate the lambda coefficient, since I only use it multiple times here and is a
    //simple calculation I rather use a lambda function instead of declaring a new proper function
    auto calc_lambda = [](auto x_plane, auto x_out, auto x_T) {return (x_plane - x_out)/x_T;};
    
    //Go through each boundary and, if the segment defined by p_inside and p_outside
    //crosses the boundary, calculate the lambda for that boundary
    
    //Check if any of the x-boundaries is intersected
    double lambda_x = -1.0;
    //x-left boundary
    if ( (p_outside.x - window_geo.poi_min.x) < Zero ) {
        
        //Calculate the lambda value
        lambda_x = calc_lambda(window_geo.poi_min.x, p_outside.x, T.x);
    }
    //x-right boundary
    else if ( (window_geo.max_x - p_outside.x) < Zero ) {
        
        //Calculate the lambda value
        lambda_x = calc_lambda(window_geo.max_x, p_outside.x, T.x);
    }
    //hout<<"lambda_x="<<lambda_x<<endl;
    
    //Check if any of the y-boundaries is intersected
    double lambda_y = -1.0;
    //y-left boundary
    if ( (p_outside.y - window_geo.poi_min.y) < Zero ) {
        
        //Calculate the lambda value
        lambda_y = calc_lambda(window_geo.poi_min.y, p_outside.y, T.y);
    }
    //y-right boundary
    else if ( (window_geo.max_y - p_outside.y) < Zero ) {
        
        //Calculate the lambda value
        lambda_y = calc_lambda(window_geo.max_y, p_outside.y, T.y);
    }
    //hout<<"lambda_y="<<lambda_y<<endl;
    
    //Check if any of the z-boundaries is intersected
    double lambda_z = -1.0;
    //z-left boundary
    if ( (p_outside.z - window_geo.poi_min.z) < Zero ) {
        
        //Calculate the lambda value
        lambda_z = calc_lambda(window_geo.poi_min.z, p_outside.z, T.z);
    }
    //z-right boundary
    else if ( (window_geo.max_z - p_outside.z) < Zero ) {
        
        //Calculate the lambda value
        lambda_z = calc_lambda(window_geo.max_z, p_outside.z, T.z);
    }
    //hout<<"lambda_z="<<lambda_z<<endl;
    
    //Get the largest lambda to calculate the intersection at the boundary
    double lambda = max(lambda_x, max(lambda_y, lambda_z));
    
    //hout<<"P_outside = ("<<p_outside.x<<", "<<p_outside.y<<", "<<p_outside.z<<")"<<endl;
    //hout<<"P_inside = ("<<p_inside.x<<", "<<p_inside.y<<", "<<p_inside.z<<")"<<endl;
    //hout<<"P_T = ("<<T.x<<", "<<T.y<<", "<<T.z<<")"<<endl;
    
    //Calculate the point at the boundary
    p_outside = p_outside + T*lambda;
    //hout<<"P_intersection = ("<<p_outside.x<<", "<<p_outside.y<<", "<<p_outside.z<<")"<<endl<<endl;
    
    //Copy the flag, i.e., the CNT number
    p_outside.flag = p_inside.flag;
    
    return 1;
}
//This function checks in which of these three location a point is placed (considering that
//there is a percolation layer):
//outside the observation window
//inside the observation window
//at the boundary of the observation window
string Cutoff_Wins::Where_is(const Point_3D &point, const cuboid &window_geo)
{
    //If the point is close enough to the boudary of the observation window, then consider this point
    //to be at the boundary
    if (point.is_at_cuboid_boundary(window_geo))
        return "boundary";
    //If the point is outside the observation window then any these conditions needs to be true
    else if (point.is_outside_cuboid(window_geo))
        return "outside";
    //If the point is not outside the observation window nor at the boundary, then it's inside
    else
        return "inside";
}
string Cutoff_Wins::Where_is_with_layer(const Point_3D &point, const cuboid &window_geo, const cuboid &layer_geom, const int &flag, const double &cnt_rad, const Point_3D &V)
{
    //Check if the point is outside the sample
    if (point.is_outside_cuboid(window_geo)) {
        
        //In this case the point is outside
        return "outside";
    }
    else if (point.is_outside_cuboid(layer_geom)) {
        
        //In this case the point is inside the sample but outside the percolation layer,
        //thus the point is within the percolation layer
        //Check if the flag is set to 1, in which case point is either the first or
        //last point of the CNT
        if (flag) {
            
            //point is the first or last point of the CNT, so check if the point orientation
            //meets the criteria to be considered inside or at boundary
            if (Check_criteria_for_cnt_at_boundary(point.x, window_geo.poi_min.x, window_geo.max_x, layer_geom.poi_min.x, layer_geom.max_x, V.x, cnt_rad) ||
                Check_criteria_for_cnt_at_boundary(point.y, window_geo.poi_min.y, window_geo.max_y, layer_geom.poi_min.y, layer_geom.max_y, V.y, cnt_rad) ||
                Check_criteria_for_cnt_at_boundary(point.z, window_geo.poi_min.z, window_geo.max_z, layer_geom.poi_min.z, layer_geom.max_z, V.z, cnt_rad)) {
                return "boundary";
            }
        }
        else {
            
            //point is not one of the endpoints of the CNT, so since point is inside the
            //boundary layer and is not an endpoint then point is considered to be at a boundary
            //hout<<"middle point considered to be at boundary"<<endl;
            return "boundary";
        }
    }
    
    //If none of the options above, then the point is definitely inside
    return "inside";
}
//This function checks the criteria in two opposite boundaries for an endpoint
//that is in the percolation layer
int Cutoff_Wins::Check_criteria_for_cnt_at_boundary(const double &p_coord, const double &min_boundary, const double &max_boundary, const double &min_layer, const double &max_layer, const double &V_cood, const double &cnt_rad)
{
    //Left boundary
    if (min_boundary - p_coord < Zero && p_coord - min_layer < Zero) {
        
        //The normal to the boundary is N = x or y or z
        //Check:
        //The dot product N.dot(V) is negative, which reduces to V_cood = Vx or Vy or Vz
        //being negative
        //AND
        //cos(alpha) = sqrt(1 - (N.dot(V))^2 ) is less than the threshold d/cnt_rad
        //where d is the distance to the boundary
        //also note that sqrt(1 - (N.dot(V))^2 ) reduces to sqrt(1 - (V_cood)^2 )
        if (V_cood < Zero && sqrt(1 - V_cood*V_cood) - (p_coord - min_boundary)/cnt_rad < Zero) {
            
            //point is in a position such that the CNT touches the boundary
            return 1;
        }
    }
    //Right boundary
    else if (max_layer < p_coord && p_coord < max_boundary) {
        
        //The normal to the boundary is N = -x or -y or -z
        //Check:
        //The dot product N.dot(V) is negative, which reduces to V_cood = Vx or Vy or Vz
        //being positive
        //AND
        //cos(alpha) = sqrt(1 - (N.dot(V))^2 ) is less than the threshold d/cnt_rad
        //where d is the distance to the boundary
        //also note that sqrt(1 - (N.dot(V))^2 ) reduces to sqrt(1 - (V_cood)^2 )
        if (V_cood > Zero && sqrt(1 - V_cood*V_cood) - (max_boundary - p_coord)/cnt_rad < Zero) {
            
            //point is in a position such that the CNT touches the boundary
            return 1;
        }
        
    }
    
    //-----------------------------------------------------------------------
    //For reference this is the code I wrote originally for the x-boundary
    //-----------------------------------------------------------------------
    //Find a boundary the point might be at
    /* /x-boundary
    if (window_geo.poi_min.x - point.x < Zero && point.x - layer_geo.poi_min.x < Zero) {
        
        //The normal to the boundary is N = x
        //Check:
        //The dot product N.dot(V) is negative, which reduces to Vx being negative
        //AND
        //cos(alpha) = sqrt(1 - (N.dot(V))^2 ) is less than the threshold d/cnt_rad
        //where d is the distance to the boundary
        //also note that sqrt(1 - (N.dot(V))^2 ) reduces to sqrt(1 - (Vx)^2 )
        if (V.x < Zero && sqrt(1 - V.x*V.x) - (point.x - window_geo.poi_min.x)/cnt_rad < Zero) {
            
            //point is in a position such that the CNT touches the boundary
            return "boundary";
        }
    }
    //Opposite x-boundary
    else if (layer_geo.max_x < point.x && point.x < window_geo.max_x) {
        
        //The normal to the boundary is N = -x
        //Check:
        //The dot product N.dot(V) is negative, which reduces to Vx being positive
        //AND
        //cos(alpha) = sqrt(1 - (N.dot(V))^2 ) is less than the threshold d/cnt_rad
        //where d is the distance to the boundary
        //also note that sqrt(1 - (N.dot(V))^2 ) reduces to sqrt(1 - (Vx)^2 )
        if (V.x > Zero && sqrt(1 - V.x*V.x) - (window_geo.max_x - point.x)/cnt_rad < Zero) {
            
            //point is in a position such that the CNT touches the boundary
            return "boundary";
        }
        
    }*/
    
    return 0;
}
string Cutoff_Wins::Where_is_with_boundary(const Point_3D &point, const cuboid &window_geo, int &boundary)
{
    double x = point.x;
    double y = point.y;
    double z = point.z;
    
    //Lambda function to determine if a coordinate is boundaded
    auto is_bounded = [](const double &coord, const double &bound_min, const double &bound_max){
        return (bound_min - coord <= Zero && coord - bound_max <= Zero);
    };
    
    //Check if coordinates are bounded by the sample
    //In case a coordinate is not bounded, then return outside directly
    bool is_x_bounded = is_bounded(x, window_geo.poi_min.x, window_geo.max_x);
    if (!is_x_bounded) {
        //hout<<"unbounded x"<<point.str()<<endl;
        return "outside";
    }
    bool is_y_bounded = is_bounded(y, window_geo.poi_min.y, window_geo.max_y);
    if (!is_y_bounded) {
        //hout<<"unbounded y"<<point.str()<<endl;
        return "outside";
    }
    bool is_z_bounded = is_bounded(z, window_geo.poi_min.z, window_geo.max_z);
    if (!is_z_bounded) {
        //hout<<"unbounded z"<<point.str()<<endl;
        return "outside";
    }
    
    //At this point, all coordinates are bounded
    //Thus, the point is either inside or at a boundary
    //Then, check if the point is close enough to the boudary of the observation window so that
    //it is considered to be at the boundary
    if (abs(x - window_geo.poi_min.x) < Zero) {
        boundary = 4;
        return "boundary";
    }
    else if(abs(x - window_geo.max_x) < Zero) {
        boundary = 2;
        return "boundary";
    }
    else if (abs(y - window_geo.poi_min.y) < Zero){
        boundary = 5;
        return "boundary";
    }
    else if(abs(y - window_geo.max_y) < Zero){
        boundary = 3;
        return "boundary";
    }
    else if (abs(z - window_geo.poi_min.z) < Zero){
        boundary = 1;
        return "boundary";
    }
    else if(abs(z - window_geo.max_z) < Zero) {
        boundary = 0;
        return "boundary";
    }
    //If the point is bounded and not at a boundary of the observation window, then it's inside
    else
        return "inside";
}
//Add a point to the corrsponding boundary vector
int Cutoff_Wins::Add_cnt_point_to_boundary_vectors(const cuboid &layer_geom, const Point_3D &P, const long int &P_num, const int &CNT_num)
{
    //Get the boundary of the point
    //int boundary = In_which_boundary(P, layer_geom);
    int boundary = -1;
    string location = Where_is_with_boundary(P, layer_geom, boundary);
    
    //Check for error
    if (boundary == -1) {
        //hout<<"Error: point should be at boundary but it is not."<<endl;
        hout<<"Error: point should be at boundary but it is not. Point location is: "<<location<<endl;
        hout<<"\tP="<<P.str()<<", P_num="<<P_num<<endl;
        return 0;
    }
    
    //Add the point and CNT number to the corresponding boundary
    //Only add the CNT if it has not been already added
    if (boundary_cnt[boundary].empty() || boundary_cnt[boundary].back() != CNT_num) {
        boundary_cnt[boundary].push_back(CNT_num);
    }
    boundary_cnt_pts[boundary].push_back(P_num);
    
    return 1;
}
//Find in which boundary a point is located, considering there is a percolation layer
int Cutoff_Wins::In_which_boundary(const Point_3D &P, const cuboid &layer_geom)
{
    //To determine the boundary number, it suffices to determine in which side of the
    //parcolation layer the point is located
    if (P.x - layer_geom.poi_min.x <= Zero) {
        //boundary = 4;
        return 4;
    }
    else if(layer_geom.max_x - P.x <= Zero) {
        //boundary = 2;
        return 2;
    }
    else if (P.y - layer_geom.poi_min.y <= Zero){
        //boundary = 5;
        return 5;
    }
    else if(layer_geom.max_y - P.y <= Zero){
        //boundary = 3;
        return 3;
    }
    else if (P.z - layer_geom.poi_min.z <= Zero){
        //boundary = 1;
        return 1;
    }
    else if(layer_geom.max_z - P.z <= Zero) {
        //boundary = 0;
        return 0;
    }
    
    return -1;
}
//Fill the vector cnts_inside
int Cutoff_Wins::Fill_cnts_inside(const vector<vector<long int> > &structure)
{
    //Flag to determine invalid CNTs
    int flag = 0;
    //Scan all CNTs in the structure
    for (int i = 0; i < (int)structure.size(); i++) {
        
        //A CNT needs at least two points
        //However, CNTs that only have two points might cause issues if they have a junction
        //This junction will be at a boundary, and if it happens to be a boundary with
        //prescribed boundary conditions, the junction might cause numerical issues
        //For instance, it can result in having a voltage at a node that is larger than
        //the presecribed one
        //Thus, all CNTs with 2 points are ignored and only those with 3 points or more are
        //considered in the simulation
        if (structure[i].size() > 2) {
            cnts_inside.push_back(i);
        }
        else if (structure[i].size() == 1) {
            hout<<"Error in Extract_observation_window. CNT "<<i<<" has only one point. A CNT must have at least 2 points."<<endl;
            flag = 1;
        }
    }
    
    //If the flag was set, then there were CNTs with one point
    //The function is not terminated at the first CNTs with one point found so that
    //all these CNTs can be displayed in the output file
    return !flag;
}
//Function that fills the vector gnps_inside
int Cutoff_Wins::Fill_gnps_inside(const int &window, const cuboid &window_geo, const vector<Shell> &shells_gnp, vector<GNP> &gnps, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp)
{
    //Initialize boundary vectors
    boundary_gnp.assign(6, vector<int>());
    boundary_gnp_pts.assign(6, vector<long int>());
    
    //Itertate over all GNPs
    for (int i = 0; i < (int)gnps.size(); i++) {
        
        //Check if GNP i is inside the window
        //hout<<"GNPi="<<i<<" GNP.flag="<<gnps[i].flag<<endl;
        if (window < shells_gnp[i].shell_min) {
            
            //The GNP is completely inside the observation window, so add it the the vector of gnps_inside
            gnps_inside.push_back(i);
        }
        else if (window >= shells_gnp[i].shell_min && window <= shells_gnp[i].shell_max) {
            
            //The GNP might be partially outside, so check if boundary points need to be added
            //hout<<"Find_gnp_boundary_points"<<endl;
            if (!Find_gnp_boundary_points(window_geo, gnps[i], structure_gnp, points_gnp)) {
                hout<<"Error in Fill_gnps_inside when calling Find_gnp_boundary_points"<<endl;
                return 0;
            }
            //hout<<endl<<endl;
        }
        //If window > shells_gnp[i].shell_max, then the GNP is outside the observation window
        //so it is ignored and not included in the vector of gnps_inside
    }
    
    return 1;
}
//This function determines if a GNP is partially inside an observation window or inside of it
//If it is partially inside, the boundary points are calculated
int Cutoff_Wins::Find_gnp_boundary_points(const cuboid &window_geo, GNP &gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp)
{
    //Vector to accumulate all boundary points vertices
    vector<vector<Point_3D> > points_acc(6);
    
    //Vector to store vertices inside the window
    vector<int> inside_v;
    
    //Vector to keep track of vertex locations
    vector<string> locations(8, "");
    
    //Find the locations of all vertices of the GNP
    for (int i = 0; i < 8; i++) {
        
        //Variable for the boundary location
        int boundary_l;
        
        //Get the location of vertex i
        //hout<<"Get the location of vertex i="<<i<<" P="<<gnp.vertices[i].str()<<endl;
        locations[i] = Where_is_with_boundary(gnp.vertices[i], window_geo, boundary_l);
        
        //Increase a counter if needed
        if (locations[i] == "inside") {
            inside_v.push_back(i);
        }
        else if (locations[i] == "boundary") {
            
            //Accumulate a vertex on a boundary on its corresponding boundary
            //hout<<"boundary_l="<<boundary_l<<endl;
            points_acc[boundary_l].push_back(gnp.vertices[i]);
        }
    }
    
    //Check if there are any vertices inside the window
    //hout<<"Check if there are any vertices inside the window"<<endl;
    if (inside_v.size() > 0) {
        
        //The GNP is at least partially inside
        //So first add it to the vector of GNPs inside
        gnps_inside.push_back(gnp.flag);
        
        //If the GNP is completely inside, then there no boundary points to calculate
        //So just terminate the function
        if (inside_v.size() == 8) {
            return 1;
        }
        
        //Accumulate all points at window boundaries
        //hout<<"Accumulate_boundary_points_due_to_intersections"<<endl;
        if (!Accumulate_boundary_points_due_to_intersections(window_geo, gnp, locations, points_acc)) {
            hout<<"Error in Find_gnp_boundary_points when calling Accumulate_boundary_points_due_to_intersections"<<endl;
            return 0;
        }
        
        //Calculate the average point of the points accumulated at each boundary and add it
        //to the vector points_gnp
        //hout<<"points_acc.size()="<<points_acc.size()<<endl;
        for (int i = 0; i < (int)points_acc.size(); i++) {
            
            //Check if any point was accumulated at boundary i
            //hout<<"points_acc[i="<<i<<"].size()="<<points_acc[i].size()<<endl;
            if (points_acc[i].size()) {
                
                //Point to store the average
                Point_3D P_avg(0,0,0);
                P_avg.flag = gnp.flag;
                
                //Add all points in the boundary
                for (int j = 0; j < (int)points_acc[i].size(); j++) {
                    //hout<<"points_acc[i][j]="<<points_acc[i][j].str()<<" is inside GNP?="<<gnp.Is_point_inside_gnp(points_acc[i][j])<<endl;
                    P_avg = P_avg + points_acc[i][j];
                }
                
                //Add GNP number to corresponding boundary
                boundary_gnp[i].push_back(gnp.flag);
                
                //Get GNP point number
                long int P_gnp_num = (long int)points_gnp.size();
                
                //Add the GNP point number for boundary points
                boundary_gnp_pts[i].push_back(P_gnp_num);
                
                //Add the GNP point number to the GNP structure
                structure_gnp[gnp.flag].push_back(P_gnp_num);
                
                //Add the average to the boundary vector
                points_gnp.push_back(P_avg/((double)points_acc[i].size()));
                //hout<<"P_avg="<<points_gnp.back().str()<<" is inside GNP?="<<gnp.Is_point_inside_gnp(points_gnp.back())<<endl;
            }
        }
        
        //Update the GNP volume
        Generate_Network GN;
        //Calculate the total GNP volume
        double gnp_vol = gnp.l*gnp.l*gnp.t;
        //hout<<"GN.Approximate_gnp_volume_inside_sample"<<endl;
        if (!GN.Approximate_gnp_volume_inside_sample(window_geo, gnp, gnp_vol)) {
            //hout<<"Error in Find_gnp_boundary_points when calling GN.Approximate_gnp_volume_inside_sample"<<endl;
            return 0;
        }
        gnp.volume = gnp_vol;
        
    }
    //If there were no vertices inside the window, even if there are vertices at the boundary,
    //the GNP can be considered to be outside the window
    //In such case, it is not added to the vector gnps_inside
    
    return 1;
}
//This function accumulates all boundary points due to intersections of GNP edges with window faces,
//and intersections due to window edges with GNP faces
//If a window vertex is inside the GNP, then it is added to the corresponding boundaries
int Cutoff_Wins::Accumulate_boundary_points_due_to_intersections(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, vector<vector<Point_3D> > &points_acc)
{
    //Find intersections of GNP edges with window boundaries
    //hout<<"Find_intersections_of_gnp_edges_with_window_boundaries"<<endl;
    if (!Find_intersections_of_gnp_edges_with_window_boundaries(window_geo, gnp, locations, points_acc)) {
        hout<<"Error in Accumulate_boundary_points_due_to_intersections when calling Find_intersections_of_GNP_edges_with_window_boundaries"<<endl;
        return 0;
    }
    
    //Find intersections of window edges with GNP surfaces
    //hout<<"Find_intersections_of_window_edges_with_gnp_faces"<<endl;
    if (!Find_intersections_of_window_edges_with_gnp_faces(window_geo, gnp, points_acc)) {
        hout<<"Error in Accumulate_boundary_points_due_to_intersections when calling Find_intersections_of_window_edges_with_gnp_faces"<<endl;
        return 0;
    }
    
    //Check if a sample vertex is inside the GNP and, if so, accumulate it
    //into the corresponding boundaries
    //hout<<"Find_window_vertex_inside_gnp"<<endl;
    if (!Find_window_vertex_inside_gnp(window_geo, gnp, points_acc)) {
        hout<<"Error in Accumulate_boundary_points_due_to_intersections when calling Find_window_vertex_inside_gnp"<<endl;
        return 0;
    }
    
    return 1;
}
//This function finds the intersections of GNP edges with faces of the observation window
int Cutoff_Wins::Find_intersections_of_gnp_edges_with_window_boundaries(const cuboid &window_geo, const GNP &gnp, const vector<string> &locations, vector<vector<Point_3D> > &points_acc)
{
    //Array of GNP edges
    Edge gnp_ed[] = {
        Edge(0,1),Edge(0,3),Edge(0,4),
        Edge(2,1),Edge(2,3),Edge(2,6),
        Edge(5,1),Edge(5,4),Edge(5,6),
        Edge(7,3),Edge(7,4),Edge(7,6)
    };
    
    
    //Iterate over all edges
    for (int i = 0; i < 12; i++) {
        
        //Get the vertices of the edge
        int v1 = gnp_ed[i].v1;
        int v2 = gnp_ed[i].v2;
        
        //Check if the GNP vertices go inside-outside
        if ( (locations[v1] == "inside" && locations[v2] == "outside") ||
             (locations[v2] == "inside" && locations[v1] == "outside") ) {
            
            //A boundary point is needed, which is initialized at a GNP vertex
            Point_3D P_out = (locations[v1] == "inside")? gnp.vertices[v2] : gnp.vertices[v1];
            Point_3D P_in = (locations[v1] == "inside")? gnp.vertices[v1] : gnp.vertices[v2];
            
            //Substitute P (which has the coordinates of the outside vertex)
            //by the intersection at the boundary
            //hout<<"P="<<P_out.str();
            if (!Substitute_boundary_point(window_geo, P_in, P_out)) {
                hout<<"Error in Find_boundary_points_partially_inside_case when calling Substitute_boundary_point"<<endl;
                return 0;
            }
            //hout<<" P_substituted="<<P_out.str();
            
            //Check in which boundary P is located
            int P_boundary;
            string P_loc = Where_is_with_boundary(P_out, window_geo, P_boundary);
            //hout<<" P_loc="<<P_loc<<endl;
            
            //Double check that P is actually at the boundary
            if (P_loc != "boundary") {
                hout<<"Error in Find_boundary_points_partially_inside_case calculating point at boundary."<<endl;
                hout<<"Point P is not at bundary but it should be since one vertex is inside the window and the second vertex is outside."<<endl;
                hout<<"Location of P is: "<<P_loc<<", P="<<P_out.str()<<endl;
                hout<<"V1="<<gnp.vertices[v1].str()<<". V2="<<gnp.vertices[v2].str()<<endl;
                return 0;
            }
            
            //Accumulate P into the corresponding boundary
            points_acc[P_boundary].push_back(P_out);
        }
        
        //Check if the vertices go outside-outside
        else if (locations[v1] == "outside" && locations[v2] == "outside") {
            
            //Go to the case where two outside vertices intersect a sample boundary
            vector<Point_3D> Pts;
            if (!Find_two_intersections_of_gnp_edges_with_window(window_geo, gnp, gnp.vertices[v1], gnp.vertices[v2], Pts)) {
                hout<<"Error in Find_boundary_points_partially_inside_case when calling Find_two_intersections_of_gnp_edges_with_window"<<endl;
                return 0;
            }
            
            //If intersections were found, add them to the accumulator
            if (Pts.size()) {
                for (int k = 0; k < (int)Pts.size(); k++) {
                    //Add to the corresponding boundary as indicated by the point flag
                    //hout<<"Pts[k="<<k<<"]="<<Pts[k].str()<<" flag="<<Pts[k].flag<<endl;
                    points_acc[Pts[k].flag].push_back(Pts[k]);
                }
            }
        }
    }
    
    return 1;
}
//This function checks if any of the corners of the window is inside a GNP
int Cutoff_Wins::Find_two_intersections_of_gnp_edges_with_window(const cuboid &window_geo, const GNP &gnp, const Point_3D &V1, const Point_3D &V2, vector<Point_3D> &Pts)
{
    //Lambda function to check if a window boundary is between two vertex coordinates
    auto find_calc_intersection = [&](const double &V1_coord, const double &V2_coord, const double &window_coord){
        
        //Check if the boundary coordinate is between the vertex coordinates
        if ((V1_coord - window_coord < Zero && V2_coord - window_coord > Zero) ||
            (V2_coord - window_coord < Zero && V1_coord - window_coord > Zero)) {
            
            //Calculate the value of lambda at the intersection with the boundary
            double lambda = (window_coord - V1_coord)/(V2_coord - V1_coord);
            
            //Calculate the point at window_coord
            Point_3D P = V1 + (V2 - V1)*lambda;
            
            //Variable to store the boundary
            int P_loc;
            
            //Check if P is actually at a face, not just the plane
            //hout<<"Check if P is actually at a face, not just the plane"<<endl;
            if (Where_is_with_boundary(P, window_geo, P_loc) == "boundary") {
                
                //Add P to the vector of points since it is actually at a face
                //Also set the point flag equal to the boundary number
                P.flag = P_loc;
                Pts.push_back(P);
            }
        }
    };
    
    //Go boundary by boundary, and check if there is an intersection
    //x-boundaries
    find_calc_intersection(V1.x, V2.x, window_geo.poi_min.x);
    find_calc_intersection(V1.x, V2.x, window_geo.max_x);
    //y-boundaries
    find_calc_intersection(V1.y, V2.y, window_geo.poi_min.y);
    find_calc_intersection(V1.y, V2.y, window_geo.max_y);
    //z-boundaries
    find_calc_intersection(V1.z, V2.z, window_geo.poi_min.z);
    find_calc_intersection(V1.z, V2.z, window_geo.max_z);
    
    return 1;
}
//This function finds the intersections of window edges with faces of a GNP
int Cutoff_Wins::Find_intersections_of_window_edges_with_gnp_faces(const cuboid &window_geo, const GNP &gnp, vector<vector<Point_3D> > &points_acc)
{
    //Array of window vertices
    Point_3D window_vertices[] = {
        Point_3D(window_geo.max_x, window_geo.max_y, window_geo.max_z),     //0
        window_geo.poi_min + Point_3D(0,window_geo.wid_y,window_geo.hei_z), //1
        window_geo.poi_min + Point_3D(0,0,window_geo.hei_z),                //2
        window_geo.poi_min + Point_3D(window_geo.len_x,0,window_geo.hei_z), //3
        window_geo.poi_min + Point_3D(window_geo.len_x,window_geo.wid_y,0), //4
        window_geo.poi_min + Point_3D(0,window_geo.wid_y,0),                //5
        window_geo.poi_min,                                                 //6
        window_geo.poi_min + Point_3D(window_geo.len_x,0,0),                //7
    };
    
    //Array of window edges
    Edge wes[] = {
        Edge(0,1),Edge(0,3),Edge(0,4),
        Edge(2,1),Edge(2,3),Edge(2,6),
        Edge(5,1),Edge(5,4),Edge(5,6),
        Edge(7,3),Edge(7,4),Edge(7,6)
    };
    
    //Array of boundaries to be added
    int boundaries[][2] = {
        {0,3}, {0,2}, {2,3},
        {0,4}, {0,5}, {4,5},
        {3,4}, {1,3}, {1,4},
        {2,5}, {1,2}, {1,5}
    };
    
    //Go through all edges of the window
    for (int i = 0; i < 12; i++) {
        
        //Check if current edge intersects a GNP face
        //hout<<"edge="<<wes[i].str()<<endl;
        Point_3D P;
        if (Does_edge_intersect_gnp(wes[i], window_vertices, gnp, P)) {
            
            //Current Edge intersects a GNP face
            //Intersection is stored in P
            //Add P to the correspoding boundaries
            int b1 = boundaries[i][0];
            int b2 = boundaries[i][1];
            
            //P is added to two boundaries because and edge is the intersection
            //of two boundaries
            points_acc[b1].push_back(P);
            points_acc[b2].push_back(P);
            //hout<<"P="<<P.str()<<" b1="<<b1<<" b2="<<b2<<" edge="<<wes[i].str()<<endl;
            //hout<<"window_vertices["<<wes[i].v1<<"]="<<window_vertices[wes[i].v1].str()<<" window_vertices["<<wes[i].v2<<"]="<<window_vertices[wes[i].v2].str()<<endl;
        }
    }
    
    return 1;
}
//This function determines if an edge defined by two points intersects a GNP face
int Cutoff_Wins::Does_edge_intersect_gnp(const Edge &edg, const Point_3D window_vertices[], const GNP &gnp, Point_3D &P)
{
    //Adjacency matrix of GNP faces
    int adj_face[][4] = {
        {2,3,4,5},
        {2,3,4,5},
        {0,1,3,5},
        {0,1,2,4},
        {0,1,3,5},
        {0,1,2,4},
    };
    
    //Vertices on each plane to calculate dot products
    int V[] = {0,4,0,0,1,2};
    
    //Vertex numbers of the window edge
    int V1 = edg.v1;
    int V2 = edg.v2;
    
    //Lambda (anonymous) function to calculate the lambda value of a line defined by two points
    auto calc_lambda = [](const Point_3D &N, const double &d, const Point_3D &P1, const Point_3D &P1P2){
        return ((-N.dot(P1) - d)/(N.dot(P1P2)));
    };
    
    //Go through all GNP faces
    for (int i = 0; i < 6; i++) {
        
        //Check wich side of GNP face i are the vertices of the edge
        int v1_loc = gnp.faces[i].N.dot(window_vertices[V1] - gnp.vertices[V[i]]) > Zero;
        int v2_loc = gnp.faces[i].N.dot(window_vertices[V2] - gnp.vertices[V[i]]) > Zero;
        //hout<<"face i="<<i<<" v1_loc="<<v1_loc<<" v2_loc="<<v2_loc<<endl;
        
        //If signs are different, then edge V1V2 intersects the plane of GNP face i
        if (v2_loc != v1_loc) {
            
            //Signs are different, so edge V1V2 intersects the plane of GNP face i
            //Now get the intersection
            
            //First calculate the lambda of the line equation
            //Calculate the vector P1P2
            Point_3D V1V2 = window_vertices[V2] - window_vertices[V1];
            double lambda = calc_lambda(gnp.faces[i].N, gnp.faces[i].coef[3], window_vertices[V1], V1V2);
            
            //Calculate the point at the plane of GNP face i
            P = window_vertices[V1] + V1V2*lambda;
            //hout<<"P.N+d="<<P.dot(gnp.faces[i].N)+gnp.faces[i].coef[3]<<endl;
            
            bool ignore_flag = false;
            
            //Check if P is bounded by the adjacent faces of GNP face i
            for (int j = 0; j < 4; j++) {
                
                //Get adjacent face j
                int face_j = adj_face[i][j];
                
                //Check if P is "above" face_j
                //hout<<"gnp.faces[face_j="<<face_j<<"].N.dot(P - V[face_j])="<<gnp.faces[face_j].N.dot(P - gnp.vertices[V[face_j]])<<endl;
                if (gnp.faces[face_j].N.dot(P - gnp.vertices[V[face_j]]) > Zero) {
                    
                    //Set the ignore flag to true
                    ignore_flag = true;
                    
                    //break the for-loop with index j
                    break;
                }
            }
            
            //Check if the intersection was valid or it should be ignored
            if (!ignore_flag) {
                
                //The intersection is valid and should not be ignored, then return 1
                return 1;
            }
        }
    }
    
    //If this part of the code is reached, then no intersections or invalid intersections were found
    return 0;
}
//This function finds if a window vertex is inside a GNP and accumulates it into the corresponding boundaries
int Cutoff_Wins::Find_window_vertex_inside_gnp(const cuboid &window_geo, const GNP &gnp, vector<vector<Point_3D> > &points_acc)
{
    //Array of window vertices
    Point_3D window_vertices[] = {
        Point_3D(window_geo.max_x, window_geo.max_y, window_geo.max_z),     //0
        window_geo.poi_min + Point_3D(0,window_geo.wid_y,window_geo.hei_z), //1
        window_geo.poi_min + Point_3D(0,0,window_geo.hei_z),                //2
        window_geo.poi_min + Point_3D(window_geo.len_x,0,window_geo.hei_z), //3
        window_geo.poi_min + Point_3D(window_geo.len_x,window_geo.wid_y,0), //4
        window_geo.poi_min + Point_3D(0,window_geo.wid_y,0),                //5
        window_geo.poi_min,                                                 //6
        window_geo.poi_min + Point_3D(window_geo.len_x,0,0),                //7
    };
    
    //Array of boundaries that need to be added if a window vertex is inside the GNP
    int boundaries[][3] = {
        {0,2,3},
        {0,3,4},
        {0,4,5},
        {0,2,5},
        {1,2,3},
        {1,3,4},
        {1,4,5},
        {1,2,5},
    };
    
    //Iterate over all window vertices
    for (int i = 0; i < 8; i++) {
        
        //Check if vertex i is inside the GNP
        if (gnp.Is_point_inside_gnp(window_vertices[i])) {
            
            //Add the window vertex to three boundaries
            int b = boundaries[i][0];
            points_acc[b].push_back(window_vertices[i]);
            b = boundaries[i][1];
            points_acc[b].push_back(window_vertices[i]);
            b = boundaries[i][2];
            points_acc[b].push_back(window_vertices[i]);
            //hout<<"window_v="<<window_vertices[i].str()<<" b1="<<boundaries[i][0]<<" b2="<<boundaries[i][1]<<" b3="<<boundaries[i][2]<<endl;
            
            //More than one vertex cannot be inside the GNP, so terminate the function
            return 1;
        }
    }
    
    return 1;
}
