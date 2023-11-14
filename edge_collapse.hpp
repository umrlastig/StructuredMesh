#ifndef EDGE_COLLAPSE_H_
#define EDGE_COLLAPSE_H_

#include "header.hpp"

#include <chrono>
#include <vector>
#include <set>

#include <CGAL/Point_set_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

typedef CGAL::Point_set_3<Exact_predicates_kernel::Point_3> Point_set;

namespace SMS = CGAL::Surface_mesh_simplification;

K::FT get_mean_point_per_area(Surface_mesh &mesh, const Point_set &point_cloud);

struct LindstromTurk_param {
	float volume_preservation;
	float boundary_preservation;
	float volume_optimisation;
	float boundary_optimization;
	float triangle_shape_optimization;
	float label_preservation;
	float semantic_border_optimization;

	LindstromTurk_param(float volume_preservation,
						float boundary_preservation,
						float volume_optimisation,
						float boundary_optimization,
						float triangle_shape_optimization,
						float label_preservation,
						float semantic_border_optimization);
};

class Custom_placement {
	const LindstromTurk_param &params;
	Surface_mesh::Property_map<Surface_mesh::Halfedge_index, K::FT> placement_costs;
	const Point_set &point_cloud;

	public:
		Custom_placement (const LindstromTurk_param &params, Surface_mesh &mesh, const Point_set &point_cloud);

		boost::optional<SMS::Edge_profile<Surface_mesh>::Point> operator()(const SMS::Edge_profile<Surface_mesh>& profile) const;
};

class Custom_cost {
	const K::FT alpha, beta, gamma, delta;
	const K::FT min_point_per_area;
	const Point_set &point_cloud;
	Surface_mesh::Property_map<Surface_mesh::Halfedge_index, K::FT> placement_costs;

	public:
		Custom_cost (const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT delta, K::FT min_point_per_area, Surface_mesh &mesh, const Point_set &point_cloud);

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
		std::chrono::time_point<std::chrono::system_clock> start_collecte;
		std::chrono::time_point<std::chrono::system_clock> start_collapse;
		bool output[30] = {false};
		CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;
		std::set<Point_set::Index> points_to_be_change;
		
		const K::FT alpha, beta;
		const K::FT min_point_per_area;
		Surface_mesh &mesh;
		const Surface_mesh_info &mesh_info;
		Point_set &point_cloud;

		Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
		Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
		Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
		Point_set::Property_map<unsigned char> point_cloud_label;

	public:
		My_visitor(const K::FT alpha, const K::FT beta, const K::FT min_point_per_area, Surface_mesh &mesh, const Surface_mesh_info &mesh_info, Point_set &point_cloud);

		void OnStarted (Surface_mesh&);

		void OnFinished (Surface_mesh&);

		void OnCollected(const SMS::Edge_profile<Surface_mesh>&, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >&);

		void OnSelected (const SMS::Edge_profile<Surface_mesh>&, boost::optional< SMS::Edge_profile<Surface_mesh>::FT > cost, const SMS::Edge_profile<Surface_mesh>::edges_size_type initial_edge_count, const SMS::Edge_profile<Surface_mesh>::edges_size_type current_edge_count);

		void OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&);

		void OnCollapsed (const SMS::Edge_profile<Surface_mesh>&, const Surface_mesh::Vertex_index vd);
};

#endif  /* !EDGE_COLLAPSE_H_ */
