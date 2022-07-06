//===========================================================================
//SOFTWARE:     3D geometric model for CNT and GNP networks
//OBJECTIVE:    Implementation of Hoshen-Kopelman Algorithm
//AUTHOR:       Angel Mora
//E-MAIL:       angel.mora@cicy.mx
//===========================================================================

#ifndef HOSHENKOPELMAN_H
#define HOSHENKOPELMAN_H

#include "time.h"
#include <algorithm>
#include <map>
#include <set>
#include "Input_Reader.h"
#include "Printer.h"
#include "Collision_detection.h"
#include "Generate_Network.h"

//-------------------------------------------------------
class Hoshen_Kopelman
{
public:
    
    //Variables
    //Cluster vectors to be used by other classes
    vector<vector<int> > clusters_cnt;
    vector<vector<int> > isolated_cnt;
    //Cluster vectors for GNPs
    vector<vector<int> > clusters_gnp;
    vector<vector<int> > isolated_gnp;
    //CNT elements
    vector<set<long int> > elements_cnt;
    //Junctions
    vector<Junction> junctions_cnt;
    vector<Junction> junctions_gnp;
    vector<Junction> junctions_mixed;
    //Junctions grouped as clusters
    vector<vector<int> > cluster_cnt_junctions;
    vector<vector<int> > cluster_gnp_junctions;
    vector<vector<int> > cluster_mix_junctions;
    //Vector with the family for each cluster
    //Families are represented by an integer as follows:
    //0 for X
    //1 for Y
    //2 for Z
    //3 for XY
    //4 for XZ
    //5 for YZ
    //6 for XYZ
    vector<int> family;
    
    //Temporary map for determining if a pair of veritces is an edge
    map<int, vector<int> > edge_map;

    //Flag for ignoring errors in mixed contacts where a CNT point is inside a GNP
    int ignore_eror_cnt_inside_gnp;
    
    //Constructor
    Hoshen_Kopelman()
    {
        //Initialize flag for ignoring errors in mixed contacts
        //where a CNT point is inside a GNP
        //0 means do not ignore error
        ignore_eror_cnt_inside_gnp = 0;
    };
    
    //
    int Determine_clusters_and_percolation(const int &iter, const cuboid &sample, const Simu_para &simu_para, const Cutoff_dist &cutoffs, const Visualization_flags &vis_flags, const vector<int> &cnts_inside, const vector<vector<long int> > &sectioned_domain_cnt, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<int> > &boundary_cnt, const vector<int> &gnps_inside, const vector<vector<int> > &sectioned_domain_gnp, vector<GNP> &gnps, const vector<vector<int> > &boundary_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp);
    int Make_cnt_clusters(const string& particle_type, const vector<Point_3D> &points_cnt, const vector<double> &radii, const Cutoff_dist &cutoffs, const vector<vector<long int> > &sectioned_domain_cnt, const vector<vector<long int> > &structure_cnt, vector<int> &labels_cnt, int &n_labels_cnt);
    int Label_cnts_in_window(const vector<Point_3D> &points_cnt, const vector<double> &radii, const Cutoff_dist& cutoffs, const vector<vector<long int> > &sectioned_domain_cnt, map<long int, map<int, long int> > &point_contacts, map<long int, map<int, double> > &point_contacts_dist, vector<map<int, set<long int> > > &contact_elements, vector<int> &labels_cnt, vector<int> &labels_labels_cnt);
    int Cleanup_labels(vector<int> &labels_labels, vector<int> &labels, int &n_labels);
    int Renumber_gnp_labels(const int &cnt_labels, vector<int> &labels_gnp);
    int Compress_cnt_cnt_contact_segments(const Cutoff_dist &cutoffs, const map<long int, map<int, long int> > &point_contacts, const map<long int, map<int, double> > &point_contacts_dist, const vector<map<int, set<long int> > > &contact_elements);
    int Complete_cnt_elements(const vector<vector<long int> > &structure_cnt);
    int Make_gnp_clusters(const cuboid& sample, const vector<int> &gnps_inside, const vector<vector<int> > &sectioned_domain_gnp, const vector<GNP> &gnps, const Cutoff_dist& cutoffs, const int &n_labels_cnt, int &n_total_labels, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp, vector<int> &labels_gnp);
    int Label_gnps_in_window(const cuboid& sample, const vector<int> &gnps_inside, const vector<vector<int> > &sectioned_domain_gnp, const vector<GNP> &gnps, const Cutoff_dist& cutoffs, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp, vector<int> &labels_gnp, vector<int> &labels_labels_gnp);
    int Get_distance_between_GNPs(const Cutoff_dist& cutoffs, const vector<GNP>& gnps, const int& GNPa, const int& GNPb, bool& p_flag, Point_3D& disp_tot, Point_3D& N, double& dist, GNP& gnp_B);
    int Find_junction_points_for_gnps(const cuboid &sample, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, vector<Point_3D> &points_gnp, Point_3D& PointA, Point_3D& PointB);
    int Find_closest_simplices_of_gnps_in_contact(const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, vector<int> &simplexA, vector<int> &simplexB, int &v_sumA, int &v_sumB);
    int Find_junction_points_in_gnps(const vector<int> &simplexA, const vector<int> &simplexB, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, const int &face_sumA, const int &face_sumB, Point_3D &PointA, Point_3D &PointB);
    int Find_point_b_for_vertex_in_simplex_a(const vector<int> &simplexB, const GNP &GNP_B, const Point_3D &N, const double &distance, const Point_3D &PointA, Point_3D &PointB);
    int Find_point_b_for_edge_in_simplex_a(const vector<int> &simplexA, const vector<int> &simplexB, const int &face_sum, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, Point_3D &PointA, Point_3D &PointB);
    int Find_point_for_edge_edge_contact(const GNP &GNP_A, const GNP &GNP_B, const Point_3D edgeA[], const vector<int> &simplexB, const Point_3D &disp, Point_3D &PointA, Point_3D &PointB);
    int Find_point_for_edge_face_contact(const int &face_sum, const GNP &GNP_A, const GNP &GNP_B, const Point_3D edgeA[], const Point_3D &disp, Point_3D &PointA, Point_3D &PointB);
    int Get_edges_of_face(const int &face_sum, vector<Edge> &edges, vector<int> &normals);
    Point_3D Get_intersection_with_face_edge(const int &countA, const int arrA[], const Point_3D &A0, const Point_3D &A1, const vector<Edge> &edges, const GNP GNP_B);
    double Lambda_of_two_lines(const Point_3D &P_in, const Point_3D &P_out, const Point_3D &edge1, const Point_3D &edge2);
    int Find_point_b_for_face_in_simplex_a(const vector<int> &simplexA, const vector<int> &simplexB, const int &face_sumA, const int &face_sumB, const GNP &GNP_A, const GNP &GNP_B, const Point_3D &N, const double &distance, Point_3D &PointA, Point_3D &PointB);
    int Get_vertices_inside_face(const vector<int> &simplexA, const vector<Edge> &edgesB, const vector<int> &normalsB, const GNP &GNP_A, const GNP &GNP_B, vector<int> &verticesA_inside);
    int Average_point_with_two_interserctions(const vector<Edge> &edgesA, const vector<int> &normalsA, const GNP &gnpA, const vector<Edge> &edgesB, const vector<int> &normalsB, const GNP &gnpB, const Point_3D &disp, const int &vB, Point_3D &PB);
    int Add_gnp_junction(const vector<GNP>& gnps, const int& GNPa, const int& GNPb, const bool& p_flag, const Point_3D& PointA, const Point_3D& PointB, const double& dist, const double& tol_gnp, vector<Point_3D>& points_gnp, vector<vector<long int> >& structure_gnp);
    int Add_gnp_point_to_vectors_if_not_repeated(const Point_3D& P, const double& tol_gnp, const bool& p_flag, long int& Pj, vector<Point_3D>& points_gnp, vector<long int>& structure_gnp);
    int Make_mixed_clusters(const int &n_labels, const Cutoff_dist &cutoffs, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, vector<int> &labels_cnt, vector<int> &labels_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp, int &n_clusters);
    int Find_adjacent_labels(const int& n_labels, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const Cutoff_dist &cutoffs, const vector<Point_3D> &points_cnt, const vector<double> &radii, const vector<vector<long int> > &sectioned_domain_cnt, const vector<GNP> &gnps, const vector<vector<int> > &sectioned_domain_gnp, vector<Point_3D> &contact_points_gnp, map<long int, map<int, long int> > &point_contacts, map<long int, map<int, double> > &point_contacts_dist, vector<map<int, set<long int> > > &gnp_cnt_point_contacts, vector<int>& labels_mixed, vector<int>& labels_labels_mixed);
    int Initialize_mixed_labels(const vector<int>& labels_cnt, const vector<int>& labels_gnp, vector<int>& labels_mixed, vector<int>& labels_labels_mixed);
    int Compress_mixed_contacts(const Cutoff_dist &cutoffs, map<long int, map<int, long int> > &point_contacts, map<long int, map<int, double> > &point_contacts_dist, vector<map<int, set<long int> > > &gnp_cnt_point_contacts, const vector<Point_3D> &contact_points_gnp, vector<vector<long int> > &structure_gnp, vector<Point_3D> &points_gnp);
    int Add_mixed_junction(const Cutoff_dist& cutoffs, const long int& Pi_junc, const long int& Pj_junc, const int& CNTj, const int& GNPi, const double& d_junc_min, const vector<Point_3D>& contact_points_gnp, vector<vector<long int> >& structure_gnp, vector<Point_3D>& points_gnp);
    int Update_cnt_and_gnp_labels(const vector<int>& labels_mixed, vector<int>& labels_cnt, vector<int>& labels_gnp);
    int Make_particle_clusters(const int &n_clusters, const vector<int> &particles_inside, vector<int> &labels, vector<vector<int> > &isolated, vector<vector<int> > &clusters_particles);
    int Find_percolated_clusters(const int &n_clusters, const vector<vector<int> > &boundary_cnt, const vector<vector<int> > &boundary_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const vector<vector<int> > &clusters_cnt_tmp, const vector<vector<int> > &clusters_gnp_tmp, map<int,int> &percoalted_labels);
    int Find_clusters_connected_to_boundaries(const vector<vector<int> > &boundary_cnt, const vector<vector<int> > &boundary_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const int boundary_pairs[][2], vector<set<int> > &percolated_dirs);
    int Add_clusters_in_boundary(const vector<int> &boundary_particles, const vector<int> &labels, set<int> &boundary1);
    int Add_percolated_direction(const int &d, const vector<int> &boundary_particles, const vector<int> &labels, const set<int> &boundary1, vector<set<int> > &percolated_dirs);
    int Determine_family_of_percolated_clusters(const int &n_clusters, const vector<vector<int> > &clusters_cnt_tmp, const vector<vector<int> > &clusters_gnp_tmp, const vector<set<int> > &percolated_dirs, map<int,int> &percolated_labels);
    int Family_map(const int &perc_sum, int &family);
    int Group_junctions(const vector<Point_3D> &points_cnt, const vector<Point_3D> &points_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const map<int,int> &percolated_labels);
    int Group_junctions_same_particle(const vector<Point_3D> &points, const vector<int> &labels, const vector<Junction> &junctions, const map<int,int> &percolated_labels, vector<vector<int> > &cluster_junction);
    int Group_junctions_mix_particle(const vector<Point_3D> &points_cnt, const vector<Point_3D> &points_gnp, const vector<int> &labels_cnt, const vector<int> &labels_gnp, const map<int,int> &percolated_labels);
    int Remove_gnp_points_in_non_relevant_boundaries(vector<vector<long int> >& structure_gnp, vector<GNP>& gnps);
    int Get_vector_with_relevant_boundaries(const int& family, vector<int> &relevant_boundaries);
    int Check_relevant_boundary_and_remove_if_needed(const int& n_rb, const int& gnp_j, const vector<GNP>& gnps, const vector<int>& relevant_boundaries, vector<pair<long int, int> >& relevant_points, vector<vector<long int> >& structure_gnp);
    bool Is_in_relevant_boundary(const int& n_rb, const int& bk, const vector<int>& relevant_boundaries);
    //-------------------------------------------------------
    //HK76
    int HK76(const int &CNT1, const int &CNT2, int &new_label, vector<int> &labels, vector<int> &labels_labels);
    int Find_root(const int &L, vector<int> &labels_labels);
    int Merge_labels(const int &root1, const int &root2, vector<int> &labels_labels);
    //-------------------------------------------------------
    //Visualization files
    int Export_clusters(const int &percolation, const int &iter, const vector<vector<int> > &clusters_cnt, const vector<vector<long int> > &structure_cnt, const vector<Point_3D> &points_cnt, const vector<vector<int> > &clusters_gnp, const vector<GNP> &gnps);
    int Combine_into_one_cluster(const vector<vector<int> > &clusters, vector<int> &cluster);
    
    
private:
    
};
//-------------------------------------------------------
#endif
//===========================================================================
