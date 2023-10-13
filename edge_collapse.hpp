#ifndef EDGE_COLLAPSE_H_
#define EDGE_COLLAPSE_H_

#include "header.hpp"
#include "raster.hpp"

#include <chrono>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <limits>

#include <CGAL/Point_set_3.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_svd.h>
#endif

#ifdef CGAL_USE_GMP
#include <CGAL/Gmpzf.h>
typedef CGAL::Gmpzf ET;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float ET;
#endif

typedef CGAL::Point_set_3<Exact_predicates_kernel::Point_3> Point_set;
typedef CGAL::Eigen_svd::Vector                             Eigen_vector;
typedef CGAL::Eigen_svd::Matrix                             Eigen_matrix;
typedef CGAL::Quadratic_program<ET>                         Program;
typedef CGAL::Quadratic_program_solution<ET>                Solution;

typedef CGAL::Search_traits_3<Exact_predicates_kernel>      Traits_base;
typedef CGAL::Search_traits_adapter<Point_set::Index, Point_set::Point_map, Traits_base> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>      Neighbor_search;
typedef Neighbor_search::Tree                               Point_tree;
typedef CGAL::AABB_tree<CGAL::AABB_traits<K, CGAL::AABB_triangle_primitive<K, std::vector<K::Triangle_3>::iterator>>> Triangle_tree;

namespace SMS = CGAL::Surface_mesh_simplification;

#include "edge_collapse.hpp"

void add_label(Surface_mesh &mesh, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face);

struct LindstromTurk_param {
	float volume_preservation;
	float boundary_preservation;
	float volume_optimisation;
	float boundary_optimization;
	float triangle_shape_optimization;
	float label_preservation;

	LindstromTurk_param(float volume_preservation,
						float boundary_preservation,
						float volume_optimisation,
						float boundary_optimization,
						float triangle_shape_optimization,
						float label_preservation);
};

class Custom_placement {
	const LindstromTurk_param &params;
	std::map<Surface_mesh::Halfedge_index, K::FT> &costs;
	const Point_set &point_cloud;
	std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face;

	public:
		Custom_placement (const LindstromTurk_param &params, std::map<Surface_mesh::Halfedge_index, K::FT> &costs, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face);

		boost::optional<SMS::Edge_profile<Surface_mesh>::Point> operator()(const SMS::Edge_profile<Surface_mesh>& profile) const;
};

class Custom_cost {
	const K::FT alpha, beta, gamma;
	const Point_set &point_cloud;
	std::map<Surface_mesh::Halfedge_index, K::FT> &costs;
	std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face;

	public:
		Custom_cost (const K::FT alpha, const K::FT beta, const K::FT gamma, std::map<Surface_mesh::Halfedge_index, K::FT> &costs, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face);

		boost::optional<SMS::Edge_profile<Surface_mesh>::FT> operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const;
};

class Cost_stop_predicate {
	public:

		Cost_stop_predicate(const float cost);

		bool operator()(const SMS::Edge_profile<Surface_mesh>::FT & current_cost, const SMS::Edge_profile<Surface_mesh> &, const SMS::Edge_profile<Surface_mesh>::edges_size_type, const SMS::Edge_profile<Surface_mesh>::edges_size_type) const;

	private:
		const float cost;
};

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh> {
	private:
		int i_collecte = 0;
		Surface_mesh &mesh;
		const Surface_mesh_info &mesh_info;
		std::chrono::time_point<std::chrono::system_clock> start_collecte;
		std::chrono::time_point<std::chrono::system_clock> start_collapse;
		bool output[30] = {false};

		const Point_set &point_cloud;
		std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face;
		std::set<Point_set::Index> points_to_be_change;

	public:
		My_visitor(Surface_mesh &mesh, const Surface_mesh_info &mesh_info, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face);

		void OnStarted (Surface_mesh&);

		void OnCollected(const SMS::Edge_profile<Surface_mesh>&, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >&);

		void OnSelected (const SMS::Edge_profile<Surface_mesh>&, boost::optional< SMS::Edge_profile<Surface_mesh>::FT > cost, const SMS::Edge_profile<Surface_mesh>::edges_size_type initial_edge_count, const SMS::Edge_profile<Surface_mesh>::edges_size_type current_edge_count);

		void OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&);

		void OnCollapsed (const SMS::Edge_profile<Surface_mesh>&, const Surface_mesh::Vertex_index vd);
};

#endif  /* !EDGE_COLLAPSE_H_ */
