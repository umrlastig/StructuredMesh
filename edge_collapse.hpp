#ifndef EDGE_COLLAPSE_H_
#define EDGE_COLLAPSE_H_

#include "header.hpp"
#include "timer.hpp"

#include <chrono>
#include <vector>
#include <set>

#include <CGAL/Point_set_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

typedef Exact_predicates_kernel Point_set_kernel;

typedef CGAL::Point_set_3<Point_set_kernel::Point_3> Point_set;

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

struct CollapseDataElement{
	Surface_mesh::Halfedge_index halfedge;
	unsigned char label;
	K::FT cost;
	std::list<Point_set::Index> points;
};

struct CollapseData {
	K::FT cost;
	std::list<CollapseDataElement> elements;
};

class Custom_placement {
	const LindstromTurk_param &params;
	Surface_mesh::Property_map<Surface_mesh::Edge_index, CollapseData> collapse_datas;
	const Point_set &point_cloud;

	public:
		Custom_placement (const LindstromTurk_param &params, Surface_mesh &mesh, const Point_set &point_cloud);

		boost::optional<SMS::Edge_profile<Surface_mesh>::Point> operator()(const SMS::Edge_profile<Surface_mesh>& profile) const;
};

class Custom_cost {
	const LindstromTurk_param &params;
	const K::FT alpha, beta, gamma, delta;
	const K::FT min_point_per_area;
	const Point_set &point_cloud;
	Surface_mesh::Property_map<Surface_mesh::Edge_index, CollapseData> collapse_datas;

	public:
		Custom_cost (const LindstromTurk_param &params, const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT delta, K::FT min_point_per_area, Surface_mesh &mesh, const Point_set &point_cloud);

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
		TimerUtils::Timer collected_timer;
		TimerUtils::Timer collapsing_timer;
		bool output[30] = {false};
		CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;
		
		const LindstromTurk_param &params;
		const K::FT alpha, beta, gamma;
		const K::FT min_point_per_area;
		Surface_mesh &mesh;
		const Surface_mesh_info &mesh_info;
		Point_set &point_cloud;

		Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
		Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
		Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
		Surface_mesh::Property_map<Surface_mesh::Edge_index, CollapseData> collapse_datas;
		Point_set::Property_map<unsigned char> point_cloud_label;

	public:
		My_visitor(const LindstromTurk_param &params, const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT min_point_per_area, Surface_mesh &mesh, const Surface_mesh_info &mesh_info, Point_set &point_cloud);

		void OnStarted (Surface_mesh&);

		void OnFinished (Surface_mesh&);

		void OnCollected(const SMS::Edge_profile<Surface_mesh>&, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >&);

		void OnSelected (const SMS::Edge_profile<Surface_mesh>&, boost::optional< SMS::Edge_profile<Surface_mesh>::FT > cost, const SMS::Edge_profile<Surface_mesh>::edges_size_type initial_edge_count, const SMS::Edge_profile<Surface_mesh>::edges_size_type current_edge_count);

		void OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&);

		void OnCollapsed (const SMS::Edge_profile<Surface_mesh>&, const Surface_mesh::Vertex_index vd);
};

#endif  /* !EDGE_COLLAPSE_H_ */
