#include "bridge.hpp"

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include "ceres/ceres.h"

#include <cmath>

namespace PMP = CGAL::Polygon_mesh_processing;

std::pair<Surface_mesh::Face_index, Point_2> point_on_path_border(const Surface_mesh &mesh, Surface_mesh::Face_index face, std::list<K::Segment_2> segments);


void add_road(
		std::map<CGAL::Straight_skeleton_2<K>::Halfedge_handle, K::FT> &roads,
		K::FT distance,
		CGAL::Straight_skeleton_2<K>::Vertex_handle vertex) {

	if (distance < 50) {
		auto he = vertex->halfedge_around_vertex_begin();
		do {
			if ((*he)->is_inner_bisector() && (*he)->opposite()->is_inner_bisector()) {
				auto length = sqrt(CGAL::squared_distance(vertex->point(), (*he)->opposite()->vertex()->point()));
				if (roads.count(*he) == 0 && roads.count((*he)->opposite()) == 0) {
					roads[*he] = distance;
					add_road(roads, distance + length, (*he)->opposite()->vertex());
				} else if (roads.count(*he) == 0) {
					if (roads[*he] > distance) {
						roads[*he] = distance;
						add_road(roads, distance + length, (*he)->opposite()->vertex());
					}
				} else {
					if (roads[(*he)->opposite()] > distance) {
						roads[(*he)->opposite()] = distance;
						add_road(roads, distance + length, (*he)->opposite()->vertex());
					}
				}
			}
		} while (++he != vertex->halfedge_around_vertex_begin());
	}

}

std::pair<K::FT, K::FT> road_width (std::pair<skeletonPoint,skeletonPoint> link) {
	K::Vector_2 vector(link.first.point, link.second.point);
	vector /= sqrt(vector.squared_length());

	std::map<CGAL::Straight_skeleton_2<K>::Halfedge_handle, K::FT> roads1;
	std::map<CGAL::Straight_skeleton_2<K>::Halfedge_handle, K::FT> roads2;

	if (link.first.vertex != nullptr) {
		add_road(roads1, 0, link.first.vertex);
	} else {
		roads1[link.first.halfedge] = 0;
		add_road(roads1, sqrt(CGAL::squared_distance(link.first.halfedge->vertex()->point(), link.first.point)), link.first.halfedge->vertex());
		add_road(roads1, sqrt(CGAL::squared_distance(link.first.halfedge->opposite()->vertex()->point(), link.first.point)), link.first.halfedge->opposite()->vertex());
	}

	if (link.second.vertex != nullptr) {
		add_road(roads2, 0, link.second.vertex);
	} else {
		roads2[link.second.halfedge] = 0;
		add_road(roads2, sqrt(CGAL::squared_distance(link.second.halfedge->vertex()->point(), link.second.point)), link.second.halfedge->vertex());
		add_road(roads2, sqrt(CGAL::squared_distance(link.second.halfedge->opposite()->vertex()->point(), link.second.point)), link.second.halfedge->opposite()->vertex());
	}

	K::FT road_width1 = 0;
	if (roads1.size() > 0) {
		K::FT sum1 = 0;
		for (auto it = roads1.begin(); it != roads1.end(); ++it) {
			auto width = it->first->vertex()->time() + it->first->opposite()->vertex()->time();
			auto vec = K::Vector_2(it->first->vertex()->point(), it->first->opposite()->vertex()->point());
			auto length = sqrt(vec.squared_length());
			auto cos_angle = abs(CGAL::scalar_product(vec, vector) / sqrt(vec.squared_length()));
			auto coef = (cos_angle/2+0.6)*length*(51*49/50 * (1/(it->second+1)+50*50/49/51-1));
			road_width1 += coef;
			sum1 += coef/width;
		}
		road_width1 /= sum1;
	} else if (link.first.vertex != nullptr) {
		road_width1 = 2*link.first.vertex->time();
	} else {
		road_width1 = link.first.halfedge->vertex()->time() + link.first.halfedge->opposite()->vertex()->time();
	}

	K::FT road_width2 = 0;
	if (roads2.size() > 0) {
		K::FT sum2 = 0;
		for (auto it = roads2.begin(); it != roads2.end(); ++it) {
			auto width = it->first->vertex()->time() + it->first->opposite()->vertex()->time();
			auto vec = K::Vector_2(it->first->vertex()->point(), it->first->opposite()->vertex()->point());
			auto length = sqrt(vec.squared_length());
			auto cos_angle = abs(CGAL::scalar_product(vec, vector) / sqrt(vec.squared_length()));
			auto coef = (cos_angle/2+0.6)*length*(51*49/50 * (1/(it->second+1)+50*50/49/51-1));
			road_width2 += coef;
			sum2 += coef/width;
		}
		road_width2 /= sum2;
	} else if (link.second.vertex != nullptr) {
		road_width2 = 2*link.second.vertex->time();
	} else {
		road_width2 = link.second.halfedge->vertex()->time() + link.second.halfedge->opposite()->vertex()->time();
	}

	return std::pair<K::FT, K::FT>(road_width1, road_width2);

}


std::set<pathLink> link_paths(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> &medial_axes, const Raster &raster) {

	K::FT minimal_path_width = raster.coord_distance_to_grid_distance(2); // in m

	// Get label property
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::set<pathLink> result;

	for (int selected_label:  {3, 8, 9}) {
		// List path with selected label
		std::list<int> same_label_paths;
		for (int i = 0; i < paths.size(); i++) {
			if (label[paths[i].front()] == selected_label && medial_axes.count(i) == 1) {
				same_label_paths.push_back(i);
			}
		}

		for (int path1: same_label_paths) {
			for (int path2: same_label_paths) {
				if (path1 < path2) {

					// Compute pair-distance between vertices and edges.
					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Vertex_handle>, K::FT> distance_v1v2;
					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Halfedge_handle>, std::pair<K::FT, K::Point_2>> distance_v1h2;
					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Halfedge_handle>, std::pair<K::FT, K::Point_2>> distance_v2h1;

					// For vertices pairs
					for (auto v1: medial_axes.at(path1)->vertex_handles()) {
						if (v1->is_skeleton()) {
							for (auto v2: medial_axes.at(path2)->vertex_handles()) {
								if (v2->is_skeleton()) {
									distance_v1v2[std::make_pair(v1, v2)] = CGAL::squared_distance(v1->point(), v2->point());
								}
							}
						}
					}

					// For vertex on path1 and edge on path2
					for (auto v1: medial_axes.at(path1)->vertex_handles()) {
						if (v1->is_skeleton()) {
							for (auto edge2: medial_axes.at(path2)->halfedge_handles()) {
								if (edge2->vertex()->id() < edge2->opposite()->vertex()->id() && edge2->is_inner_bisector() && edge2->opposite()->is_inner_bisector()) {
									auto segment = K::Segment_2(edge2->opposite()->vertex()->point(), edge2->vertex()->point());
									auto proj = segment.supporting_line().projection(v1->point());
									if (segment.collinear_has_on(proj)) {
										distance_v1h2[std::make_pair(v1, edge2)] = std::make_pair(CGAL::squared_distance(v1->point(), proj), proj);
									}
								}
							}
						}
					}

					// For vertex on path2 and edge on path1
					for (auto v2: medial_axes.at(path2)->vertex_handles()) {
						if (v2->is_skeleton()) {
							for (auto edge1: medial_axes.at(path1)->halfedge_handles()) {
								if (edge1->vertex()->id() < edge1->opposite()->vertex()->id() && edge1->is_inner_bisector() && edge1->opposite()->is_inner_bisector()) {
									auto segment = K::Segment_2(edge1->opposite()->vertex()->point(), edge1->vertex()->point());
									auto proj = segment.supporting_line().projection(v2->point());
									if (segment.collinear_has_on(proj)) {
										distance_v2h1[std::make_pair(v2, edge1)] = std::make_pair(CGAL::squared_distance(v2->point(), proj), proj);
									}
								}
							}
						}
					}

					// For vertices pairs
					for (auto it = distance_v1v2.begin(); it != distance_v1v2.end(); ++it) {
						auto v1 = it->first.first;
						auto v2 = it->first.second;
						auto d = it->second;

						// Exit link
						Exact_predicates_kernel::Segment_2 segment(Exact_predicates_kernel::Point_2(v1->point().x(), v1->point().y()), Exact_predicates_kernel::Point_2(v2->point().x(), v2->point().y()));
						int intersect = 0;
						for (auto edge = path_polygon.at(path1).outer_boundary().edges_begin(); intersect <= 1 && edge != path_polygon.at(path1).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect++;
							}
						}
						for (auto hole: path_polygon.at(path1).holes()) {
							for (auto edge = hole.edges_begin(); intersect <= 1 && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect++;
								}
							}
						}
						if (intersect != 1) continue;

						auto he = v1->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (distance_v1v2[std::make_pair(v, v2)] < d) {
									goto exit1;
								}
								if ((*he)->vertex()->id() < (*he)->opposite()->vertex()->id()) {
									if (distance_v2h1.count(std::make_pair(v2, *he)) > 0 && distance_v2h1[std::make_pair(v2, *he)].first < d) {
										goto exit1;
									}
								} else {
									if (distance_v2h1.count(std::make_pair(v2, (*he)->opposite())) > 0 && distance_v2h1[std::make_pair(v2, (*he)->opposite())].first < d) {
										goto exit1;
									}
								}
							}
						} while (++he != v1->halfedge_around_vertex_begin());

						he = v2->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (distance_v1v2[std::make_pair(v1, v)] < d) {
									goto exit1;
								}
								if ((*he)->vertex()->id() < (*he)->opposite()->vertex()->id()) {
									if (distance_v1h2.count(std::make_pair(v1, *he)) > 0 && distance_v1h2[std::make_pair(v1, *he)].first < d) {
										goto exit1;
									}
								} else {
									if (distance_v1h2.count(std::make_pair(v1, (*he)->opposite())) > 0 && distance_v1h2[std::make_pair(v1, (*he)->opposite())].first < d) {
										goto exit1;
									}
								}
							}
						} while (++he != v2->halfedge_around_vertex_begin());

						result.insert(std::make_pair(skeletonPoint(path1, v1), skeletonPoint(path2, v2)));
						exit1:
							continue;
					}

					// For vertex on path1 and edge on path2
					for (auto it = distance_v1h2.begin(); it != distance_v1h2.end(); ++it) {
						auto v1 = it->first.first;
						auto e2 = it->first.second;
						auto d = it->second.first;
						auto p2 = it->second.second;

						// Exit link
						Exact_predicates_kernel::Segment_2 segment(Exact_predicates_kernel::Point_2(v1->point().x(), v1->point().y()), Exact_predicates_kernel::Point_2(p2.x(), p2.y()));
						int intersect = 0;
						for (auto edge = path_polygon.at(path1).outer_boundary().edges_begin(); intersect <= 1 && edge != path_polygon.at(path1).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect++;
							}
						}
						for (auto hole: path_polygon.at(path1).holes()) {
							for (auto edge = hole.edges_begin(); intersect <= 1 && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect++;
								}
							}
						}
						if (intersect != 1) continue;

						auto he = v1->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (CGAL::squared_distance(v->point(), p2) < d) {
									goto exit2;
								}
								auto segment = K::Segment_2((*he)->opposite()->vertex()->point(), (*he)->vertex()->point());
								auto proj = segment.supporting_line().projection(p2);
								if (segment.collinear_has_on(proj)) {
									goto exit2;
								}
							}
						} while (++he != v1->halfedge_around_vertex_begin());

						result.insert(std::make_pair(skeletonPoint(path1, v1), skeletonPoint(path2, e2, p2)));
						exit2:
							continue;
					}

					// For vertex on path2 and edge on path1
					for (auto it = distance_v2h1.begin(); it != distance_v2h1.end(); ++it) {
						auto v2 = it->first.first;
						auto e1 = it->first.second;
						auto d = it->second.first;
						auto p1 = it->second.second;

						// Exit link
						Exact_predicates_kernel::Segment_2 segment(Exact_predicates_kernel::Point_2(v2->point().x(), v2->point().y()), Exact_predicates_kernel::Point_2(p1.x(), p1.y()));
						int intersect = 0;
						for (auto edge = path_polygon.at(path2).outer_boundary().edges_begin(); intersect <= 1 && edge != path_polygon.at(path2).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect++;
							}
						}
						for (auto hole: path_polygon.at(path2).holes()) {
							for (auto edge = hole.edges_begin(); intersect <= 1 && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect++;
								}
							}
						}
						if (intersect != 1) continue;

						auto he = v2->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (CGAL::squared_distance(v->point(), p1) < d) {
									goto exit3;
								}
								auto segment = K::Segment_2((*he)->opposite()->vertex()->point(), (*he)->vertex()->point());
								auto proj = segment.supporting_line().projection(p1);
								if (segment.collinear_has_on(proj)) {
									goto exit3;
								}
							}
						} while (++he != v2->halfedge_around_vertex_begin());

						result.insert(std::make_pair(skeletonPoint(path1, e1, p1), skeletonPoint(path2, v2)));
						exit3:
							continue;
					}

				} else if (path1 == path2) {

					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Vertex_handle>, K::FT> distance_vv;
					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Halfedge_handle>, std::pair<K::FT, K::Point_2>> distance_vh;

					// Compute pair-distance between vertices and edges.
					for (auto v1: medial_axes.at(path1)->vertex_handles()) {
						if (v1->is_skeleton()) {
							// For vertices pairs
							for (auto v2: medial_axes.at(path1)->vertex_handles()) {
								if (v2->is_skeleton()) {
									distance_vv[std::make_pair(v1, v2)] = CGAL::squared_distance(v1->point(), v2->point());
								}
							}
							// For vertex and edge
							for (auto edge2: medial_axes.at(path1)->halfedge_handles()) {
								if (edge2->vertex()->id() < edge2->opposite()->vertex()->id() && edge2->is_inner_bisector() && edge2->opposite()->is_inner_bisector()) {
									auto segment = K::Segment_2(edge2->opposite()->vertex()->point(), edge2->vertex()->point());
									auto proj = segment.supporting_line().projection(v1->point());
									if (segment.collinear_has_on(proj)) {
										distance_vh[std::make_pair(v1, edge2)] = std::make_pair(CGAL::squared_distance(v1->point(), proj), proj);
									}
								}
							}
						}
					}

					// For vertices pairs
					for (auto it = distance_vv.begin(); it != distance_vv.end(); ++it) {
						auto v1 = it->first.first;
						auto v2 = it->first.second;
						auto d = it->second;

						if (v1->id() == v2->id()) continue;

						// Exit link
						Exact_predicates_kernel::Segment_2 segment(Exact_predicates_kernel::Point_2(v1->point().x(), v1->point().y()), Exact_predicates_kernel::Point_2(v2->point().x(), v2->point().y()));
						int intersect = 0;
						for (auto edge = path_polygon.at(path1).outer_boundary().edges_begin(); intersect <= 2 && edge != path_polygon.at(path1).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect++;
							}
						}
						for (auto hole: path_polygon.at(path1).holes()) {
							for (auto edge = hole.edges_begin(); intersect <= 2 && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect++;
								}
							}
						}
						if (intersect != 2) continue;

						auto he = v1->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (distance_vv[std::make_pair(v, v2)] < d) {
									goto exit4;
								}
								if ((*he)->vertex()->id() < (*he)->opposite()->vertex()->id()) {
									if (distance_vh.count(std::make_pair(v2, *he)) > 0 && distance_vh[std::make_pair(v2, *he)].first < d) {
										goto exit4;
									}
								} else {
									if (distance_vh.count(std::make_pair(v2, (*he)->opposite())) > 0 && distance_vh[std::make_pair(v2, (*he)->opposite())].first < d) {
										goto exit4;
									}
								}
							}
						} while (++he != v1->halfedge_around_vertex_begin());

						he = v2->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (distance_vv[std::make_pair(v1, v)] < d) {
									goto exit4;
								}
								if ((*he)->vertex()->id() < (*he)->opposite()->vertex()->id()) {
									if (distance_vh.count(std::make_pair(v1, *he)) > 0 && distance_vh[std::make_pair(v1, *he)].first < d) {
										goto exit4;
									}
								} else {
									if (distance_vh.count(std::make_pair(v1, (*he)->opposite())) > 0 && distance_vh[std::make_pair(v1, (*he)->opposite())].first < d) {
										goto exit4;
									}
								}
							}
						} while (++he != v2->halfedge_around_vertex_begin());

						result.insert(std::make_pair(skeletonPoint(path1, v1), skeletonPoint(path1, v2)));
						exit4:
							continue;
					}

					// For vertex on path1 and edge on path2
					for (auto it = distance_vh.begin(); it != distance_vh.end(); ++it) {
						auto v1 = it->first.first;
						auto e2 = it->first.second;
						auto d = it->second.first;
						auto p2 = it->second.second;

						if (v1->id() == e2->vertex()->id() || v1->id() == e2->opposite()->vertex()->id()) continue;

						// Exit link
						Exact_predicates_kernel::Segment_2 segment(Exact_predicates_kernel::Point_2(v1->point().x(), v1->point().y()), Exact_predicates_kernel::Point_2(p2.x(), p2.y()));
						int intersect = 0;
						for (auto edge = path_polygon.at(path1).outer_boundary().edges_begin(); intersect <= 2 && edge != path_polygon.at(path1).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect++;
							}
						}
						for (auto hole: path_polygon.at(path1).holes()) {
							for (auto edge = hole.edges_begin(); intersect <= 2 && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect++;
								}
							}
						}
						if (intersect != 2) continue;

						auto he = v1->halfedge_around_vertex_begin();
						do {
							auto v = (*he)->opposite()->vertex();
							if (v->is_skeleton()) {
								if (CGAL::squared_distance(v->point(), p2) < d) {
									goto exit5;
								}
								auto segment = K::Segment_2((*he)->opposite()->vertex()->point(), (*he)->vertex()->point());
								auto proj = segment.supporting_line().projection(p2);
								if (segment.collinear_has_on(proj)) {
									goto exit5;
								}
							}
						} while (++he != v1->halfedge_around_vertex_begin());

						result.insert(std::make_pair(skeletonPoint(path1, v1), skeletonPoint(path1, e2, p2)));
						exit5:
							continue;
					}
				}
			}
		}
	}

	// Remove bridges between too small path
	for (auto it = result.begin(); it != result.end(); ) {
		auto width = road_width(*it);
		if (width.first < minimal_path_width || width.second < minimal_path_width) {
			it = result.erase(it);
		} else {
			++it;
		}
	}

	Surface_mesh links;

	for(auto link: result) {
		auto z1 = raster.dsm[int(link.first.point.y())][int(link.first.point.x())];
		auto z2 = raster.dsm[int(link.second.point.y())][int(link.second.point.x())];
		auto v1 = links.add_vertex(Point_3((float) link.first.point.x(), (float) link.first.point.y(), z1));
		auto v2 = links.add_vertex(Point_3((float) link.second.point.x(), (float) link.second.point.y(), z2));
		links.add_edge(v1,v2);
	}

	bool created;
	Surface_mesh::Property_map<Surface_mesh::Edge_index, int> edge_prop;
	boost::tie(edge_prop, created) = links.add_property_map<Surface_mesh::Edge_index, int>("prop",0);
	assert(created);

	double min_x, min_y;
	raster.grid_to_coord(0, 0, min_x, min_y);

	for(auto vertex : links.vertices()) {
		auto point = links.point(vertex);
		double x, y;
		raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
		links.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
	}

	std::ofstream mesh_ofile ("links.ply");
	CGAL::IO::write_PLY (mesh_ofile, links);
	mesh_ofile.close();

	return result;

}

pathBridge::pathBridge(pathLink link): link(link), cost(0) {
	N = int(sqrt(CGAL::squared_distance(link.first.point, link.second.point)));
	xl = new double[N+1];
	xr = new double[N+1];
	z_segment = new double[N+1];
}

pathBridge::pathBridge(const pathBridge& other) : link(other.link), label(other.label), N(other.N), cost(other.cost), crossing_surface(other.crossing_surface) {
	// Copy all elements from 'other' to 'data'
	xl = new double[N+1];
	xr = new double[N+1];
	z_segment = new double[N+1];
	std::copy_n(other.xl, N + 1, xl);
	std::copy_n(other.xr, N + 1, xr);
	std::copy_n(other.z_segment, N + 1, z_segment);
}

pathBridge::~pathBridge() {
	delete [] xl;
	delete [] xr;
	delete [] z_segment;
}

// Surface solving

// regularity of the surface
struct surface_regularity {
	double coef;

	surface_regularity (double coef) : coef(coef) {}

	template <typename T>
	bool operator()(const T* const z0, const T* const z1, T* residual) const {
		residual[0] = (z0[0] - z1[0])*coef;
		return true;
	}
};

// attachment to DSM data 
class SurfaceCost : public ceres::SizedCostFunction<1, 1, 1, 1> {
	private:
		double coef;
		double cost;

		Point_2 start;
		K::Vector_2 ortho_vect;
		unsigned char label;
		double tunnel_height;
		const Surface_mesh &mesh;
		const AABB_tree &tree;
		Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_labels;

		double cost_at_point(const double j, const double* const z, double *grad) const {
			double local_cost = 0;

			auto p2d = start + j * ortho_vect;
			auto p3d = K::Point_3(p2d.x(), p2d.y(), z[0]);

			const K::Ray_3 ray_top(p3d, K::Direction_3(0, 0, 1));
			const K::Ray_3 ray_bottom(p3d, K::Direction_3(0, 0, -1));

			auto location_top = PMP::locate_with_AABB_tree(ray_top, tree, mesh);
			auto location_bottom = PMP::locate_with_AABB_tree(ray_bottom, tree, mesh);

			//if no bottom point, we are outside the bounding boxe
			if(location_bottom.first != mesh.null_face()) {

				if (location_top.first == mesh.null_face()) {
					// point is above the surface
					auto l = mesh_labels[location_bottom.first];
					auto point_bottom = PMP::construct_point(location_bottom, mesh);
					local_cost = abs(z[0] - ((double) point_bottom.z()));
					if (grad != nullptr) *grad = (z[0] > point_bottom.z()) ? 1 : -1;
					if (l != 0 && l != label) {
						if ((l != 8 && l != 9) || (label != 8 && label != 9)) {
							local_cost += cost;
						}
					}
				} else {
					auto p1_top = mesh.point(CGAL::source(CGAL::halfedge(location_top.first, mesh), mesh));
					auto p2_top = mesh.point(CGAL::target(CGAL::halfedge(location_top.first, mesh), mesh));
					auto p3_top = mesh.point(CGAL::target(CGAL::next(CGAL::halfedge(location_top.first, mesh), mesh), mesh));

					if (K::Orientation_3()(p1_top, p2_top, p3_top, p3d) == CGAL::POSITIVE) {
						// the point is above the surface
						auto l = mesh_labels[location_bottom.first];
						auto point_bottom = PMP::construct_point(location_bottom, mesh);
						local_cost = abs(z[0] - ((double) point_bottom.z()));
						if (grad != nullptr) *grad = (z[0] > point_bottom.z()) ? 1 : -1;
						if (l != 0 && l != label) {
							if ((l != 8 && l != 9) || (label != 8 && label != 9)) {
								local_cost += cost;
							}
						}

						auto point_top = PMP::construct_point(location_top, mesh);
						if (point_top.z() - z[0] < tunnel_height) {
							local_cost += (tunnel_height - (point_top.z() - z[0])) / 2;
							if (grad != nullptr) *grad += 1/2;
						}
					} else {
						// the point is under the surface
						auto point_top = PMP::construct_point(location_top, mesh);
						if (point_top.z() - z[0] < tunnel_height/2) {
							auto l = mesh_labels[location_top.first];
							local_cost = ((double) point_top.z()) - z[0];
							if (grad != nullptr) *grad = -1;
							if (l != 0 && l != label) {
								if ((l != 8 && l != 9) || (label != 8 && label != 9)) {
									local_cost += cost;
								}
							}
						} else if (point_top.z() - z[0] < tunnel_height) {
							local_cost = z[0] + tunnel_height - point_top.z();
							if (grad != nullptr) *grad = 1;
						} else {
							local_cost = 0;
							if (grad != nullptr) *grad = 0;
						}
					}
				}
			}

			return local_cost;
		}

	public:
		SurfaceCost (double coef,
					double cost,
					Point_2 start,
					K::Vector_2 ortho_vect,
					unsigned char label,
					double tunnel_height,
					const Surface_mesh &mesh,
					const AABB_tree &tree) :
		coef(coef),
		cost(cost),
		start(start),
		ortho_vect(ortho_vect),
		label(label),
		tunnel_height(tunnel_height),
		mesh(mesh),
		tree(tree) {
			// Get label property
			bool has_label;
			boost::tie(mesh_labels, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);
		}

		bool Evaluate(double const* const* parameters, double* residual, double** jacobians) const {
			const double* const xl = parameters[0];
			const double* const xr = parameters[1];
			const double* const z = parameters[2];

			if (xr[0] + xl[0] < 0) { // the road has a negative width
				residual[0] = (-xl[0] - xr[0])*coef*10;

				if (jacobians != nullptr) {
					jacobians[0][0] = - coef * 10; // dC/dxl
					jacobians[1][0] = - coef * 10; // dC/dxr
					jacobians[2][0] = 0;
				}

				return true;
			}

			double step = 0.3;

			int num_step = ceil((xl[0] + xr[0]) / step);
			if (num_step < 1) num_step = 1;
			step = (xl[0] + xr[0]) / num_step;

			if (jacobians != nullptr) {
				double grad;

				double left_cost = cost_at_point(-xl[0], z, &grad);
				residual[0] = left_cost * step / 2;
				jacobians[0][0] = coef * left_cost;
				jacobians[2][0] = grad * step / 2;

				for (int i = 1; i < num_step; i++) {
					residual[0] += cost_at_point(i * (xl[0] + xr[0]) / num_step - xl[0], z, &grad) * step;
					jacobians[2][0] += grad * step;
				}

				double right_cost = cost_at_point(xr[0], z, &grad);
				residual[0] += right_cost * step / 2;
				jacobians[1][0] = coef * right_cost ;
				jacobians[2][0] += grad * step / 2;

				jacobians[2][0] *= coef;
			} else {

				residual[0] = cost_at_point(-xl[0], z, nullptr);
				residual[0] *= step / 2;

				for (int i = 1; i < num_step; i++) {
					residual[0] += cost_at_point(i * (xl[0] + xr[0]) / num_step - xl[0], z, nullptr) * step;
				}

				double right_cost = cost_at_point(xr[0], z, nullptr);
				residual[0] += right_cost * step / 2;
			}

			residual[0] *= coef;

			return true;

		}
};

// border
struct surface_border {
	double coef;
	double border_z;

	surface_border(double coef, double border_z) : coef(coef), border_z(border_z) {}

	template <typename T>
	bool operator()(const T* const z, T* residual) const {
		residual[0] = (z[0] - border_z)*coef;
		return true;
	}
};

//regularity of the contour
/*struct contour_regularity {
	double coef;

	contour_regularity (double coef) : coef(coef) {}

	template <typename T>
	bool operator()(const T* const x0, const T* const x1, const T* const x2, T* residual) const {
		residual[0] = (x0[0] - 2.0 * x1[0] + x2[0])*coef;
		return true;
	}
};*/

struct contour_regularity {
	double coef;

	contour_regularity (double coef) : coef(coef) {}

	template <typename T>
	bool operator()(const T* const x0, const T* const x1, T* residual) const {
		residual[0] = (x0[0] - x1[0])*coef;
		return true;
	}
};

//width of the reconstructed surface
struct surface_width {
	double coef;
	double width;

	surface_width(double coef, double width) : coef(coef), width(width) {}

	template <typename T>
	bool operator()(const T* const xl, const T* const xr, T* residual) const {
		residual[0] = (xl[0] + xr[0] - width)*coef;
		return true;
	}
};

//centering of the surface on the link vertices
struct surface_centering {
	double coef;

	surface_centering(double coef) : coef(coef) {}

	template <typename T>
	bool operator()(const T* const xl, const T* const xr, T* residual) const {
		residual[0] = (xl[0] - xr[0])*coef;
		return true;
	}
};

//constraint border inside path
struct border_constraint {
	double coef;
	double max_value;

	border_constraint (double coef, double max_value) : coef(coef), max_value(max_value) {}

	template <typename T>
	bool operator()(const T* const x, T* residual) const {
		if (x[0] > max_value) {
			residual[0] = (x[0] - max_value) * coef;
		} else {
			residual[0] = ((T) 0.0);
		}
		return true;
	}
};

pathBridge bridge (pathLink link, const Surface_mesh &mesh, const AABB_tree &tree, const Raster &raster) {

	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
	Filtered_graph filtered_sm1(mesh, link.first.path, path);
	Filtered_graph filtered_sm2(mesh, link.second.path, path);

	typedef CGAL::dynamic_vertex_property_t<Point_2>               Point_2_property;
	auto projection_pmap = CGAL::get(Point_2_property(), mesh);

	for(auto v : CGAL::vertices(mesh)) {
		const Point_3& p = mesh.point(v);
		put(projection_pmap, v, Point_2(p.x(), p.y()));
	}

	auto location1 = PMP::locate(link.first.point, filtered_sm1, CGAL::parameters::vertex_point_map(projection_pmap));
	auto location2 = PMP::locate(link.second.point, filtered_sm2, CGAL::parameters::vertex_point_map(projection_pmap));
	auto point1 = PMP::construct_point(location1, mesh);
	auto point2 = PMP::construct_point(location2, mesh);

	K::Vector_2 link_vector(link.first.point, link.second.point);
	float length = sqrt(link_vector.squared_length());
	auto l = link_vector / length;
	auto n = l.perpendicular(CGAL::COUNTERCLOCKWISE);

	auto width = road_width(link);

	Point_2 left1 = link.first.point - n*width.first/2;
	Point_2 right1 = link.first.point + n*width.first/2;
	Point_2 left2 = link.second.point - n*width.second/2;
	Point_2 right2 = link.second.point + n*width.second/2;

	auto left1_border = point_on_path_border(mesh, location1.first, {K::Segment_2(link.first.point, left1)});
	auto right1_border = point_on_path_border(mesh, location1.first, {K::Segment_2(link.first.point, right1)});
	auto left2_border = point_on_path_border(mesh, location2.first, {K::Segment_2(link.second.point, left2)});
	auto right2_border = point_on_path_border(mesh, location2.first, {K::Segment_2(link.second.point, right2)});

	float dl0 = sqrt(CGAL::squared_distance(link.first.point, left1_border.second));
	float dr0 = sqrt(CGAL::squared_distance(link.first.point, right1_border.second));
	float dlN = sqrt(CGAL::squared_distance(link.second.point, left2_border.second));
	float drN = sqrt(CGAL::squared_distance(link.second.point, right2_border.second));

	int N = int(sqrt(link_vector.squared_length ()));
	pathBridge bridge(link);
	for (int i = 0; i <= bridge.N; i++) {
		bridge.xl[i] = dl0 + ((float) i)/bridge.N*(dlN-dl0);
		bridge.xr[i] = dr0 + ((float) i)/bridge.N*(drN-dr0);
		bridge.z_segment[i] = point1.z() + (point2.z() - point1.z())*((float) i)/bridge.N;
	}
	bridge.label = label[location1.first];

	float tunnel_height = 3; // in meter

	ceres::Problem problem;

	float alpha = 10; // regularity of the surface
	float beta = 1; // attachment to DSM data
	float gamma = 1; // regularity of the contour
	float delta = 2; // width of the reconstructed surface
	float epsilon = 1; // centering of the surface on the link vertices
	float zeta = 10; // border elevation
	float eta = 100; // constraint border inside path
	float theta = 25; // cost for label error

	// regularity of the surface
	for (int i = 0; i < bridge.N; i++) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<surface_regularity, 1, 1, 1>(new surface_regularity(alpha)),
			nullptr,
			bridge.z_segment + i, //z_segment[i] 
			bridge.z_segment + i + 1); // z_segment[i+1]
	}

	// attachment to DSM data
	for (int i = 0; i <= bridge.N; i++) {
		problem.AddResidualBlock(
			new SurfaceCost(beta, theta, link.first.point + ((float) i)/bridge.N*link_vector, n, bridge.label, tunnel_height, mesh, tree),
			nullptr,
			bridge.xl + i, //x^l_{i}
			bridge.xr + i, //x^r_{i}
			bridge.z_segment + i); //z_segment[i]
	}

	// border
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<surface_border, 1, 1>(new surface_border(zeta, point1.z())),
		nullptr,
		bridge.z_segment); //z_segment[0]
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<surface_border, 1, 1>(new surface_border(zeta, point2.z())),
		nullptr,
		bridge.z_segment + bridge.N); //z_segment[bridge.N]

	//regularity of the contour
	/*for (int j = 1; j < bridge.N; j++) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<contour_regularity, 1, 1, 1, 1>(new contour_regularity(gamma)),
			nullptr,
			bridge.xl + j - 1, //x^l_{j-1}
			bridge.xl + j, //x^l_{j}
			bridge.xl + j + 1); //x^l_{j+1}
	}
	for (int j = 1; j < bridge.N; j++) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<contour_regularity, 1, 1, 1, 1>(new contour_regularity(gamma)),
			nullptr,
			bridge.xr + j - 1, //x^r_{j-1}
			bridge.xr + j, //x^r_{j}
			bridge.xr + j + 1); //x^r_{j+1}
	}*/

	for (int j = 0; j < bridge.N; j++) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<contour_regularity, 1, 1, 1>(new contour_regularity(gamma)),
			nullptr,
			bridge.xl + j, //x^l_{j}
			bridge.xl + j + 1); //x^l_{j+1}
	}
	for (int j = 0; j < bridge.N; j++) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<contour_regularity, 1, 1, 1>(new contour_regularity(gamma)),
			nullptr,
			bridge.xr + j, //x^r_{j}
			bridge.xr + j + 1); //x^r_{j+1}
	}

	//width of the reconstructed surface
	for (int j = 0; j <= bridge.N; j++) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<surface_width, 1, 1, 1>(new surface_width(delta, width.first + j*(width.second - width.first)/bridge.N)),
			nullptr,
			bridge.xl + j, //x^l_{j}
			bridge.xr + j); //x^r_{j}
	}

	//centering of the surface on the link vertices
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<surface_centering, 1, 1, 1>(new surface_centering(epsilon)),
		nullptr,
		bridge.xl, //x^l_{0}
		bridge.xr); //x^r_{0}
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<surface_centering, 1, 1, 1>(new surface_centering(epsilon)),
		nullptr,
		bridge.xl + bridge.N, //x^l_{N}
		bridge.xr + bridge.N); //x^r_{N}

	//constraint border inside path
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<border_constraint, 1, 1>(new border_constraint(eta, dl0)),
		nullptr,
		bridge.xl); //x^l_{0}
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<border_constraint, 1, 1>(new border_constraint(eta, dlN)),
		nullptr,
		bridge.xl + bridge.N); //x^l_{N}
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<border_constraint, 1, 1>(new border_constraint(eta, dr0)),
		nullptr,
		bridge.xr); //x^r_{0}
	problem.AddResidualBlock(
		new ceres::AutoDiffCostFunction<border_constraint, 1, 1>(new border_constraint(eta, drN)),
		nullptr,
		bridge.xr + bridge.N); //x^r_{N}

	// solving
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.use_nonmonotonic_steps = true;
	options.logging_type = ceres::SILENT;
	options.minimizer_progress_to_stdout = true;
	ceres::Solver::Summary summary;
	Solve(options, &problem, &summary);

	//Border constraint
	if (bridge.xl[0] > dl0) bridge.xl[0] = dl0;
	if (bridge.xl[bridge.N] > dlN) bridge.xl[bridge.N] = dlN;
	if (bridge.xr[0] > dr0) bridge.xr[0] = dr0;
	if (bridge.xr[bridge.N] > drN) bridge.xr[bridge.N] = drN;

	//Width constraint
	for (int j = 0; j <= bridge.N; j++) {
		if (bridge.xl[j] + bridge.xr[j] < 0) { // x^l_j − x^r_j
			bridge.xl[j] += (-bridge.xl[j] - bridge.xr[j]) / 2;
			bridge.xr[j] += (-bridge.xl[j] - bridge.xr[j]) / 2;
		}
	}

	// Compute cost

	//regularity of the contour
	/*for (int j = 1; j < bridge.N; j++) {
		bridge.cost += pow(gamma * (bridge.xl[j-1]+bridge.xl[j+1]-2*bridge.xl[j]),2);
		bridge.cost += pow(gamma * (bridge.xr[j-1]+bridge.xr[j+1]-2*bridge.xr[j]),2);
	}*/
	for (int j = 0; j < bridge.N; j++) {
		// x^l_j − x^l_{j+1}
		bridge.cost += pow(gamma * (bridge.xl[j]-bridge.xl[j+1]),2);
		bridge.cost += pow(gamma * (bridge.xr[j]-bridge.xr[j+1]),2);
	}

	//width of the reconstructed surface
	for (int j = 0; j <= bridge.N; j++) {
		// x^l_j − x^l_{j+1}
		bridge.cost += pow(delta * ((bridge.xl[j] + bridge.xr[j]) - (width.first + j*(width.second - width.first)/bridge.N)),2);
	}

	//centering of the surface on the link vertices
	bridge.cost += pow(epsilon * (bridge.xl[0] - bridge.xr[0]),2);
	bridge.cost += pow(epsilon * (bridge.xl[bridge.N] - bridge.xr[bridge.N]),2);

	// regularity of the surface
	for (int i = 0; i < bridge.N; i++) {
		bridge.cost += pow(alpha * (bridge.z_segment[i] - bridge.z_segment[i+1]),2);
	}

	Surface_mesh bridge_surface_cost;
	std::set<Surface_mesh::Face_index> bridge_crossing;

	// attachment to DSM data
	for (int i = 0; i <= bridge.N; i++) {
		float tmp_cost = 0;

		double step = 0.3;
		int num_step = ceil((bridge.xl[i] + bridge.xr[i]) / step);
		if (num_step < 1) num_step = 1;
		step = (bridge.xl[i] + bridge.xr[i]) / num_step;

		for (int k = 0; k <= num_step; k++) {

			double j = (bridge.xl[i] + bridge.xr[i]) * k / num_step - bridge.xl[i];
			double point_cost = 0;

			auto p2d = link.first.point + ((float) i)/bridge.N*link_vector + j * n;
			auto p3d = K::Point_3(p2d.x(), p2d.y(), bridge.z_segment[i]);

			const K::Ray_3 ray_top(p3d, K::Direction_3(0, 0, 1));
			const K::Ray_3 ray_bottom(p3d, K::Direction_3(0, 0, -1));

			auto location_top = PMP::locate_with_AABB_tree(ray_top, tree, mesh);
			auto location_bottom = PMP::locate_with_AABB_tree(ray_bottom, tree, mesh);

			//if no bottom point, we are outside the bounding boxe
			if(location_bottom.first != mesh.null_face()) {

				if (location_top.first == mesh.null_face()) {
					// point is above the surface
					auto l = label[location_bottom.first];
					auto point_bottom = PMP::construct_point(location_bottom, mesh);
					point_cost = abs(bridge.z_segment[i] - point_bottom.z());
					if (l != 0 && l != bridge.label) {
						if ((l != 8 && l != 9) || (bridge.label != 8 && bridge.label != 9)) {
							point_cost += theta;
						} else {
							bridge_crossing.insert(location_bottom.first);
						}
					}
				} else {
					auto p1_top = mesh.point(CGAL::source(CGAL::halfedge(location_top.first, mesh), mesh));
					auto p2_top = mesh.point(CGAL::target(CGAL::halfedge(location_top.first, mesh), mesh));
					auto p3_top = mesh.point(CGAL::target(CGAL::next(CGAL::halfedge(location_top.first, mesh), mesh), mesh));

					if (K::Orientation_3()(p1_top, p2_top, p3_top, p3d) == CGAL::POSITIVE) {
						// the point is above the surface
						auto l = label[location_bottom.first];
						auto point_bottom = PMP::construct_point(location_bottom, mesh);
						point_cost = abs(bridge.z_segment[i] - point_bottom.z());
						if (l != 0 && l != bridge.label) {
							if ((l != 8 && l != 9) || (bridge.label != 8 && bridge.label != 9)) {
								point_cost += theta;
							} else {
								bridge_crossing.insert(location_bottom.first);
							}
						}

						auto point_top = PMP::construct_point(location_top, mesh);
						if (point_top.z() - bridge.z_segment[i] < tunnel_height) {
							point_cost += abs(tunnel_height - (point_top.z() - bridge.z_segment[i])) / 2;
						}
					} else {
						// the point is under the surface
						auto point_top = PMP::construct_point(location_top, mesh);
						if (point_top.z() - bridge.z_segment[i] < tunnel_height/2) {
							auto l = label[location_top.first];
							point_cost = abs(bridge.z_segment[i] - ((double) point_top.z()));
							if (l != 0 && l != bridge.label) {
								if ((l != 8 && l != 9) || (bridge.label != 8 && bridge.label != 9)) {
									point_cost += theta;
								} else {
									bridge_crossing.insert(location_top.first);
								}
							}
						} else if (point_top.z() - bridge.z_segment[i] < tunnel_height) {
							point_cost = abs(tunnel_height - (point_top.z() - bridge.z_segment[i]));
						}
					}
				}
			}

			if (k == 0 || k == num_step) {
				point_cost /= 2;
			}

			tmp_cost += point_cost;

			auto v1 = bridge_surface_cost.add_vertex(p3d);
			auto v2 = bridge_surface_cost.add_vertex(K::Point_3(p3d.x(), p3d.y(), p3d.z() + 50*(point_cost)));
			bridge_surface_cost.add_edge(v1, v2);
		}

		tmp_cost *= step;

		bridge.cost += pow(beta * tmp_cost, 2);
	}

	{ // surface_cost
		bool created;
		Surface_mesh::Property_map<Surface_mesh::Edge_index, int> edge_prop;
		boost::tie(edge_prop, created) = bridge_surface_cost.add_property_map<Surface_mesh::Edge_index, int>("prop",0);
		assert(created);

		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : bridge_surface_cost.vertices()) {
			auto point = bridge_surface_cost.point(vertex);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			bridge_surface_cost.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream bridge_surface_cost_name;
		bridge_surface_cost_name << "bridge_surface_cost_" << ((int) bridge.label) << "_" << link.first.path << "_" << link.second.path << "_" << link.first.point << "_" << link.second.point << ".ply";
		std::ofstream mesh_ofile (bridge_surface_cost_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, bridge_surface_cost);
		mesh_ofile.close();
	}

	// border
	bridge.cost += pow(zeta * (bridge.z_segment[0] - point1.z()),2);
	bridge.cost += pow(zeta * (bridge.z_segment[bridge.N] - point2.z()),2);

	//Compute crossing surface
	bridge.crossing_surface.clear();

	if (bridge_crossing.size() > 0) {

		CGAL::Face_filtered_graph<Surface_mesh> filtered_graph(mesh, bridge_crossing);

		CGAL::copy_face_graph(filtered_graph, bridge.crossing_surface);
	}

	if (bridge.cost > 50) {
		return bridge;
	}

	{ // Surface

		Surface_mesh bridge_mesh;
		Surface_mesh::Vertex_index Xl[bridge.N+1];
		Surface_mesh::Vertex_index Xr[bridge.N+1];

		// Add points
		for (int i = 0; i <= bridge.N; i++) {
			auto p1 = link.first.point + ((float) i)/bridge.N*link_vector - bridge.xl[i] * n;
			auto p2 = link.first.point + ((float) i)/bridge.N*link_vector + bridge.xr[i] * n;
			Xl[i] = bridge_mesh.add_vertex(Point_3(p1.x(), p1.y(), bridge.z_segment[i]));
			Xr[i] = bridge_mesh.add_vertex(Point_3(p2.x(), p2.y(), bridge.z_segment[i]));
		}

		// Add faces
		for (int i = 0; i < bridge.N; i++) {
			bridge_mesh.add_face(Xl[i], Xr[i], Xr[i+1]);
			bridge_mesh.add_face(Xl[i], Xr[i+1], Xl[i+1]);
		}

		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : bridge_mesh.vertices()) {
			auto point = bridge_mesh.point(vertex);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			bridge_mesh.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream bridge_mesh_name;
		bridge_mesh_name << "bridge_mesh_" << ((int) bridge.label) << "_" << link.first.path << "_" << link.second.path << "_" << link.first.point << "_" << link.second.point << "(" << bridge.cost << ").ply";
		std::ofstream mesh_ofile (bridge_mesh_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, bridge_mesh);
		mesh_ofile.close();

	}


	{ // Skeleton
		Surface_mesh skeleton;

		auto v1 = skeleton.add_vertex(point1);
		auto v2 = skeleton.add_vertex(point2);
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left1_border.second.x(), left1_border.second.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(left2_border.second.x(), left2_border.second.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(right1_border.second.x(), right1_border.second.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(right2_border.second.x(), right2_border.second.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left1.x(), left1.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(right1.x(), right1.y(), point1.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left2.x(), left2.y(), point2.z()));
		v2 = skeleton.add_vertex(Point_3(right2.x(), right2.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		/*auto init = skeleton.add_vertex(Point_3(Xr[0].x(), Xr[0].y(), point1.z()));
		auto temp1 = CGAL::SM_Vertex_index(init);
		for (int i = 1; i <= bridge.N; i++) {
			auto temp2 = skeleton.add_vertex(Point_3(Xr[i].x(), Xr[i].y(), point1.z()));
			skeleton.add_edge(temp1, temp2);
			temp1 = temp2;
		}
		for (int i = bridge.N; i >= 0; i--) {
			auto temp2 = skeleton.add_vertex(Point_3(Xl[i].x(), Xl[i].y(), point1.z()));
			skeleton.add_edge(temp1, temp2);
			temp1 = temp2;
		}
		skeleton.add_edge(temp1, init);*/

		Surface_mesh::Property_map<Surface_mesh::Edge_index, int> edge_prop;
		bool created;
		boost::tie(edge_prop, created) = skeleton.add_property_map<Surface_mesh::Edge_index, int>("prop",0);
		assert(created);

		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : skeleton.vertices()) {
			auto point = skeleton.point(vertex);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			skeleton.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream skeleton_name;
		skeleton_name << "bridge_" << ((int) bridge.label) << "_" << link.first.path << "_" << link.second.path << "_" << link.first.point << "_" << link.second.point << ".ply";
		std::ofstream mesh_ofile (skeleton_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, skeleton);
		mesh_ofile.close();
	}

	return bridge;

}

void close_surface_mesh(Surface_mesh &mesh) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, bool> true_face;
	bool created;
	boost::tie(true_face, created) = mesh.add_property_map<Surface_mesh::Face_index, bool>("true_face", true);
	assert(created);

	// Minimal elevation
	float min_h = 50;
	for (auto point: mesh.points()) {
		if (point.z() < min_h) {
			min_h = point.z();
		}
	}
	min_h -= 50;

	float min_x;
	float min_y;

	std::vector<Surface_mesh::Halfedge_index> borders_edge;
	CGAL::Polygon_mesh_processing::extract_boundary_cycles (mesh, std::back_inserter(borders_edge));
	assert(borders_edge.size() == 1);
	auto border_edge = borders_edge[0];

	std::vector<Point_3> points_bottom;
	std::vector<Point_3> points_top;
	std::vector<Surface_mesh::Vertex_index> vertices_bottom;

	auto halfedge = border_edge;
	do {
		auto vt = CGAL::target(halfedge, mesh);
		auto pt = mesh.point(vt);
		auto vb = mesh.add_vertex(Point_3(pt.x(), pt.y(), min_h));
		vertices_bottom.push_back(vb);
		points_top.push_back(pt);
		points_bottom.push_back(Point_3(pt.x(), pt.y(), min_h));
		halfedge = CGAL::next(halfedge, mesh);
	} while (halfedge != border_edge);

	std::vector<CGAL::Triple<int, int, int>> patch;
	patch.reserve(points_bottom.size() -2);
	CGAL::Polygon_mesh_processing::triangulate_hole_polyline (points_bottom, points_top, std::back_inserter(patch));

	for(auto face: patch) {
		auto f = mesh.add_face(vertices_bottom[face.first], vertices_bottom[face.second], vertices_bottom[face.third]);
		true_face[f] = false;
	}

	auto v0 = CGAL::target(halfedge, mesh);
	auto v1 = CGAL::target(CGAL::next(halfedge, mesh), mesh);

	auto f = mesh.add_face(v0, vertices_bottom.at(1), vertices_bottom.at(0));
	true_face[f] = false;
	f = mesh.add_face(v0, v1, vertices_bottom.at(1));
	true_face[f] = false;

	bool exist;
	boost::tie(border_edge, exist) = CGAL::halfedge(v0, vertices_bottom.at(0), mesh);
	assert(exist);

	std::vector<Surface_mesh::Face_index> patch_facets;
	CGAL::Polygon_mesh_processing::triangulate_hole(mesh, border_edge, std::back_inserter(patch_facets));
	for (auto face: patch_facets) {
		true_face[face] = false;
	}

	assert(!CGAL::Polygon_mesh_processing::does_self_intersect(mesh));
	assert(CGAL::is_closed(mesh));
	CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(mesh);
}

AABB_tree index_surface_mesh(Surface_mesh &mesh) {
	AABB_tree tree;
	PMP::build_AABB_tree(mesh, tree);
	return tree;
}

class Cost_stop_predicate {
	public:

		Cost_stop_predicate(const float cost) : cost(cost) {}

		bool operator()(const CGAL::Surface_mesh_simplification::Edge_profile<Surface_mesh>::FT & current_cost, const CGAL::Surface_mesh_simplification::Edge_profile<Surface_mesh> &, const CGAL::Surface_mesh_simplification::Edge_profile<Surface_mesh>::edges_size_type, const CGAL::Surface_mesh_simplification::Edge_profile<Surface_mesh>::edges_size_type) const {
			return current_cost > cost;
		}

	private:
		const float cost;
};

struct CorefinementVisitor : public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<Surface_mesh> {
	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	Surface_mesh::Property_map<Surface_mesh::Face_index, bool> true_face;

	int current_path;
	unsigned char current_label;
	bool current_true_face;
	const Surface_mesh *main_mesh;

	CorefinementVisitor (const Surface_mesh *mesh) : main_mesh(mesh) {
		if (main_mesh != nullptr) {
			bool has_path;
			boost::tie(path, has_path) = mesh->property_map<Surface_mesh::Face_index, int>("path");
			assert(has_path);

			bool has_label;
			boost::tie(label, has_label) = mesh->property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);

			bool has_true_face;
			boost::tie(true_face, has_true_face) = mesh->property_map<Surface_mesh::Face_index, bool>("true_face");
			assert(has_true_face);
		}
	}

	void before_subface_creations (Surface_mesh::Face_index f_split, const Surface_mesh &tm) {
		if (&tm == main_mesh) {
			current_path = path[f_split];
			current_label = label[f_split];
			current_true_face = true_face[f_split];
		} else {
			Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> bridge_label;
			bool has_label;
			boost::tie(bridge_label, has_label) = tm.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);
			current_label = bridge_label[f_split];
		}
	}

	void after_subface_created (Surface_mesh::Face_index f_new, const Surface_mesh &tm) {
		if (&tm == main_mesh) {
			path[f_new] = current_path;
			label[f_new] = current_label;
			true_face[f_new] = current_true_face;
		} else {
			Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> bridge_label;
			bool has_label;
			boost::tie(bridge_label, has_label) = tm.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);
			bridge_label[f_new] = current_label;
		}
	}

	void after_face_copy (Surface_mesh::Face_index f_src, const Surface_mesh &tm_src, Surface_mesh::Face_index f_tgt, const Surface_mesh &tm_tgt) {
		if (main_mesh != nullptr) {
			assert(&tm_tgt == main_mesh);
			assert(&tm_src != main_mesh);

			Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> bridge_label;
			bool has_label;
			boost::tie(bridge_label, has_label) = tm_src.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);
			label[f_tgt] = bridge_label[f_src];
		} else {
			Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> bridge_label1, bridge_label2;
			bool has_label;
			boost::tie(bridge_label1, has_label) = tm_src.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);
			boost::tie(bridge_label2, has_label) = tm_tgt.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);
			bridge_label2[f_tgt] = bridge_label1[f_src];
		}
	}

};

Surface_mesh compute_remove_mesh(const pathBridge &bridge, const Raster &raster) {
	float tunnel_height = 3; // in meter

	K::Vector_2 link_vector(bridge.link.first.point, bridge.link.second.point);
	float length = sqrt(link_vector.squared_length());
	auto l = link_vector / length;
	auto n = l.perpendicular(CGAL::COUNTERCLOCKWISE);

	Surface_mesh bridge_mesh;
	Surface_mesh::Vertex_index Xlb[bridge.N+1]; // X left bottom
	Surface_mesh::Vertex_index Xrb[bridge.N+1]; // X right bottom
	Surface_mesh::Vertex_index Xlt[bridge.N+1]; // X left top
	Surface_mesh::Vertex_index Xrt[bridge.N+1]; // X right top

	// Add points
	for (int i = 0; i <= bridge.N; i++) {
		auto p1 = bridge.link.first.point + ((float) i)/bridge.N*link_vector - bridge.xl[i] * n;
		auto p2 = bridge.link.first.point + ((float) i)/bridge.N*link_vector + bridge.xr[i] * n;
		Xlb[i] = bridge_mesh.add_vertex(Point_3(p1.x(), p1.y(), bridge.z_segment[i]));
		Xrb[i] = bridge_mesh.add_vertex(Point_3(p2.x(), p2.y(), bridge.z_segment[i]));
		Xlt[i] = bridge_mesh.add_vertex(Point_3(p1.x(), p1.y(), bridge.z_segment[i] + tunnel_height));
		Xrt[i] = bridge_mesh.add_vertex(Point_3(p2.x(), p2.y(), bridge.z_segment[i] + tunnel_height));
	}

	// Label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool created_label;
	boost::tie(label, created_label) = bridge_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label", 0);
	assert(created_label);

	// Add faces
	bridge_mesh.add_face(Xrt[0], Xlt[0], Xlb[0]);
	bridge_mesh.add_face(Xrt[0], Xlb[0], Xrb[0]);

	for (int i = 0; i < bridge.N; i++) {
		bridge_mesh.add_face(Xlt[i], Xrt[i], Xlt[i+1]);
		bridge_mesh.add_face(Xrt[i], Xrt[i+1], Xlt[i+1]);

		auto f1 = bridge_mesh.add_face(Xlb[i], Xrb[i+1], Xrb[i]);
		auto f2 = bridge_mesh.add_face(Xlb[i], Xlb[i+1], Xrb[i+1]);
		label[f1] = bridge.label;
		label[f2] = bridge.label;

		bridge_mesh.add_face(Xlt[i], Xlt[i+1], Xlb[i]);
		bridge_mesh.add_face(Xlt[i+1], Xlb[i+1], Xlb[i]);

		bridge_mesh.add_face(Xrt[i], Xrb[i+1], Xrt[i+1]);
		bridge_mesh.add_face(Xrt[i], Xrb[i], Xrb[i+1]);
	}

	bridge_mesh.add_face(Xlt[bridge.N], Xrb[bridge.N], Xlb[bridge.N]);
	bridge_mesh.add_face(Xlt[bridge.N], Xrt[bridge.N], Xrb[bridge.N]);

	CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(bridge_mesh);

	{ // output
		Surface_mesh output_mesh(bridge_mesh);
		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : output_mesh.vertices()) {
			auto point = output_mesh.point(vertex);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			output_mesh.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream bridge_mesh_name;
		bridge_mesh_name << "bridge_tube_" << ((int) bridge.label) << "_" << bridge.link.first.path << "_" << bridge.link.second.path << "_" << bridge.link.first.point << "_" << bridge.link.second.point << "(" << bridge.cost << ").ply";
		std::ofstream mesh_ofile (bridge_mesh_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, output_mesh);
		mesh_ofile.close();
	}

	CGAL::Cartesian_converter<Surface_mesh::Point::R, CGAL::Exact_predicates_exact_constructions_kernel> to_exact;

	Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3> exact_points;
	bool created_point;
	boost::tie(exact_points, created_point) = bridge_mesh.add_property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
	assert(created_point);

	for (auto vertex : bridge_mesh.vertices()) {
		exact_points[vertex] = to_exact(bridge_mesh.point(vertex));
	}

	return bridge_mesh;
}

Surface_mesh compute_support_mesh(const pathBridge &bridge, const Raster &raster) {
	float tunnel_height = 3; // in meter

	K::Vector_2 link_vector(bridge.link.first.point, bridge.link.second.point);
	float length = sqrt(link_vector.squared_length());
	auto l = link_vector / length;
	auto n = l.perpendicular(CGAL::COUNTERCLOCKWISE);

	Surface_mesh bridge_mesh;
	Surface_mesh::Vertex_index Xlb[bridge.N+1]; // X left bottom
	Surface_mesh::Vertex_index Xrb[bridge.N+1]; // X right bottom
	Surface_mesh::Vertex_index Xlt[bridge.N+1]; // X left top
	Surface_mesh::Vertex_index Xrt[bridge.N+1]; // X right top

	// Add points
	for (int i = 0; i <= bridge.N; i++) {
		auto p1 = bridge.link.first.point + ((float) i)/bridge.N*link_vector - bridge.xl[i] * n;
		auto p2 = bridge.link.first.point + ((float) i)/bridge.N*link_vector + bridge.xr[i] * n;
		Xlt[i] = bridge_mesh.add_vertex(Point_3(p1.x(), p1.y(), bridge.z_segment[i]));
		Xrt[i] = bridge_mesh.add_vertex(Point_3(p2.x(), p2.y(), bridge.z_segment[i]));
		Xlb[i] = bridge_mesh.add_vertex(Point_3(p1.x(), p1.y(), bridge.z_segment[i] - tunnel_height/6));
		Xrb[i] = bridge_mesh.add_vertex(Point_3(p2.x(), p2.y(), bridge.z_segment[i] - tunnel_height/6));
	}

	// Label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool created_label;
	boost::tie(label, created_label) = bridge_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label", 0);
	assert(created_label);

	// Add faces
	bridge_mesh.add_face(Xrt[0], Xlt[0], Xlb[0]);
	bridge_mesh.add_face(Xrt[0], Xlb[0], Xrb[0]);

	for (int i = 0; i < bridge.N; i++) {
		auto f1 = bridge_mesh.add_face(Xlt[i], Xrt[i], Xlt[i+1]);
		auto f2 = bridge_mesh.add_face(Xrt[i], Xrt[i+1], Xlt[i+1]);
		label[f1] = bridge.label;
		label[f2] = bridge.label;

		bridge_mesh.add_face(Xlb[i], Xrb[i+1], Xrb[i]);
		bridge_mesh.add_face(Xlb[i], Xlb[i+1], Xrb[i+1]);

		bridge_mesh.add_face(Xlt[i], Xlt[i+1], Xlb[i]);
		bridge_mesh.add_face(Xlt[i+1], Xlb[i+1], Xlb[i]);

		bridge_mesh.add_face(Xrt[i], Xrb[i+1], Xrt[i+1]);
		bridge_mesh.add_face(Xrt[i], Xrb[i], Xrb[i+1]);
	}

	bridge_mesh.add_face(Xlt[bridge.N], Xrb[bridge.N], Xlb[bridge.N]);
	bridge_mesh.add_face(Xlt[bridge.N], Xrt[bridge.N], Xrb[bridge.N]);

	CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(bridge_mesh);

	{ // output
		Surface_mesh output_mesh(bridge_mesh);
		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : output_mesh.vertices()) {
			auto point = output_mesh.point(vertex);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			output_mesh.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream bridge_mesh_name;
		bridge_mesh_name << "bridge_support_" << ((int) bridge.label) << "_" << bridge.link.first.path << "_" << bridge.link.second.path << "_" << bridge.link.first.point << "_" << bridge.link.second.point << "(" << bridge.cost << ").ply";
		std::ofstream mesh_ofile (bridge_mesh_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, output_mesh);
		mesh_ofile.close();
	}

	CGAL::Cartesian_converter<Surface_mesh::Point::R, CGAL::Exact_predicates_exact_constructions_kernel> to_exact;

	Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3> exact_points;
	bool created_point;
	boost::tie(exact_points, created_point) = bridge_mesh.add_property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
	assert(created_point);

	for (auto vertex : bridge_mesh.vertices()) {
		exact_points[vertex] = to_exact(bridge_mesh.point(vertex));
	}

	return bridge_mesh;
}

struct Extrudor {
	bool top;
	float tunnel_height = 3; // in meter
	Surface_mesh &omesh;
	Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3> exact_point;

	Extrudor(bool top, Surface_mesh &omesh) : top(top), omesh(omesh) {
		bool has_exact_points;
		boost::tie(exact_point, has_exact_points) = omesh.property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
		assert(has_exact_points);
	}

	void operator()(const Surface_mesh::Vertex_index &vin, const Surface_mesh::Vertex_index &vout) const {
		if (top) {
			exact_point[vout] += CGAL::Exact_predicates_exact_constructions_kernel::Vector_3(0, 0, tunnel_height);
		} else { // bottom
			exact_point[vout] += CGAL::Exact_predicates_exact_constructions_kernel::Vector_3(0, 0, -tunnel_height);
		}
	}
};

Surface_mesh compute_crossing_mesh(Surface_mesh mesh, const pathBridge &bridge, std::list<Point_2> *polygon, const Raster &raster) {

	// Get label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	// Compute bridge border
	K::Vector_2 link_vector(bridge.link.first.point, bridge.link.second.point);
	float length = sqrt(link_vector.squared_length());
	auto l = link_vector / length;
	auto n = l.perpendicular(CGAL::COUNTERCLOCKWISE);

	// Add points
	for (int i = 0; i <= bridge.N; i++) {
		polygon->push_back(bridge.link.first.point + ((float) i)/bridge.N*link_vector + bridge.xr[i] * n);
	}
	for (int i = bridge.N; i >= 0; i--) {
		polygon->push_back(bridge.link.first.point + ((float) i)/bridge.N*link_vector - bridge.xl[i] * n);
	}

	// Compute crossing border

	// Exact points
	CGAL::Cartesian_converter<Surface_mesh::Point::R, CGAL::Exact_predicates_exact_constructions_kernel> to_exact;
	Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3> exact_points;
	bool created_point;

	// crossing_mesh
	Surface_mesh crossing_mesh;

	boost::tie(exact_points, created_point) = crossing_mesh.add_property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
	assert(created_point);

	Extrudor extrudor1 (false, crossing_mesh);
	Extrudor extrudor2 (true, crossing_mesh);
	CGAL::Polygon_mesh_processing::extrude_mesh (bridge.crossing_surface, crossing_mesh, extrudor1, extrudor2, CGAL::parameters::all_default(), CGAL::parameters::vertex_point_map(exact_points));
	assert(!CGAL::Polygon_mesh_processing::does_self_intersect(crossing_mesh, CGAL::parameters::vertex_point_map(exact_points)));

	{ // output
		CGAL::Cartesian_converter<CGAL::Exact_predicates_exact_constructions_kernel, Surface_mesh::Point::R> to_float;
		Surface_mesh output_mesh(crossing_mesh);
		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : output_mesh.vertices()) {
			auto point = to_float(exact_points[vertex]);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			output_mesh.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream bridge_mesh_name;
		bridge_mesh_name << "crossing_mesh_" << ((int) bridge.label) << "_" << bridge.link.first.path << "_" << bridge.link.second.path << "_" << bridge.link.first.point << "_" << bridge.link.second.point << ".ply";
		std::ofstream mesh_ofile (bridge_mesh_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, output_mesh);
		mesh_ofile.close();
	}

	return crossing_mesh;
}

void add_bridge_to_mesh(Surface_mesh &mesh, const std::vector<pathBridge> &bridges, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const Raster &raster) {

	// Label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	CGAL::Cartesian_converter<Surface_mesh::Point::R, CGAL::Exact_predicates_exact_constructions_kernel> to_exact;
	CGAL::Cartesian_converter<CGAL::Exact_predicates_exact_constructions_kernel, Surface_mesh::Point::R> to_kernel;

	// Add exact points to mesh
	Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3> exact_points;
	bool created;
	boost::tie(exact_points, created) = mesh.add_property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
	assert(created);

	for (auto vertex : mesh.vertices()) {
		exact_points[vertex] = to_exact(mesh.point(vertex));
	}

	Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3> vpm1, vpm2, vpm3;
	bool has_exact_points;

	for (auto bridge: bridges) {

		std::cout << "Bridge " << bridge.link.first.path << " (" << bridge.link.first.point << ") -> " << bridge.link.second.path << " (" << bridge.link.second.point << ")\n";

		if (bridge.crossing_surface.num_faces() > 0) {

			std::list<Point_2> polygon;
			Surface_mesh mesh_crossing = compute_crossing_mesh(mesh, bridge, &polygon, raster);

			auto r_b = compute_remove_mesh(bridge, raster);
			auto s_b = compute_support_mesh(bridge, raster);

			boost::tie(vpm1, has_exact_points) = r_b.property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
			assert(has_exact_points);
			boost::tie(vpm2, has_exact_points) = s_b.property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
			assert(has_exact_points);
			boost::tie(vpm3, has_exact_points) = mesh_crossing.property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
			assert(has_exact_points);

			// Crossing relabling

			CGAL::Polygon_mesh_processing::corefine(r_b, mesh_crossing, CGAL::parameters::vertex_point_map(vpm1).visitor(CorefinementVisitor(nullptr)), CGAL::parameters::vertex_point_map(vpm3).do_not_modify(true));

			CGAL::Polygon_mesh_processing::corefine(s_b, mesh_crossing, CGAL::parameters::vertex_point_map(vpm2).visitor(CorefinementVisitor(nullptr)), CGAL::parameters::vertex_point_map(vpm3).do_not_modify(true));

			CGAL::Side_of_triangle_mesh<Surface_mesh, CGAL::Exact_predicates_exact_constructions_kernel, Surface_mesh::Property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>> inside(mesh_crossing, vpm3);

			bool has_label;

			Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label_r_b;
			boost::tie(label_r_b, has_label) = r_b.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);

			for (auto face: r_b.faces()) {
				if (label_r_b[face] == bridge.label) {
					auto p1 = vpm1[CGAL::source(CGAL::halfedge(face, r_b), r_b)];
					auto p2 = vpm1[CGAL::target(CGAL::halfedge(face, r_b), r_b)];
					auto p3 = vpm1[CGAL::target(CGAL::next(CGAL::halfedge(face, r_b), r_b), r_b)];

					if (inside(CGAL::centroid(p1, p2, p3)) == CGAL::ON_BOUNDED_SIDE) {
						label_r_b[face] = 11;
					}
				}
			}

			Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label_s_b;
			boost::tie(label_s_b, has_label) = s_b.property_map<Surface_mesh::Face_index, unsigned char>("label");
			assert(has_label);

			for (auto face: s_b.faces()) {
				if (label_s_b[face] == bridge.label) {
					auto p1 = vpm2[CGAL::source(CGAL::halfedge(face, s_b), s_b)];
					auto p2 = vpm2[CGAL::target(CGAL::halfedge(face, s_b), s_b)];
					auto p3 = vpm2[CGAL::target(CGAL::next(CGAL::halfedge(face, s_b), s_b), s_b)];

					if (inside(CGAL::centroid(p1, p2, p3)) == CGAL::ON_BOUNDED_SIDE) {
						label_s_b[face] = 11;
					}
				}
			}

			std::cout << "Remove bridge: " << CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh, r_b, mesh, CGAL::parameters::vertex_point_map(exact_points).visitor(CorefinementVisitor(&mesh)), CGAL::parameters::vertex_point_map(vpm1), CGAL::parameters::vertex_point_map(exact_points)) << "\n";

			std::cout << "Support bridge: " << CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, s_b, mesh, CGAL::parameters::vertex_point_map(exact_points).visitor(CorefinementVisitor(&mesh)), CGAL::parameters::vertex_point_map(vpm2), CGAL::parameters::vertex_point_map(exact_points)) << "\n";

			assert(!CGAL::Polygon_mesh_processing::does_self_intersect(mesh, CGAL::parameters::vertex_point_map(exact_points)));
			assert(std::cout << CGAL::Polygon_mesh_processing::does_bound_a_volume(mesh, CGAL::parameters::vertex_point_map(exact_points)));

		} else {

			auto r_b = compute_remove_mesh(bridge, raster);
			boost::tie(vpm1, has_exact_points) = r_b.property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
			assert(has_exact_points);
			std::cout << "Remove bridge: " << CGAL::Polygon_mesh_processing::corefine_and_compute_difference(mesh, r_b, mesh, CGAL::parameters::vertex_point_map(exact_points).visitor(CorefinementVisitor(&mesh)), CGAL::parameters::vertex_point_map(vpm1), CGAL::parameters::vertex_point_map(exact_points)) << "\n";

			auto s_b = compute_support_mesh(bridge, raster);
			boost::tie(vpm2, has_exact_points) = s_b.property_map<Surface_mesh::Vertex_index, CGAL::Exact_predicates_exact_constructions_kernel::Point_3>("v:exact_point");
			assert(has_exact_points);
			std::cout << "Support bridge: " << CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh, s_b, mesh, CGAL::parameters::vertex_point_map(exact_points).visitor(CorefinementVisitor(&mesh)), CGAL::parameters::vertex_point_map(vpm2), CGAL::parameters::vertex_point_map(exact_points)) << "\n";

			assert(!CGAL::Polygon_mesh_processing::does_self_intersect(mesh, CGAL::parameters::vertex_point_map(exact_points)));
			assert(std::cout << CGAL::Polygon_mesh_processing::does_bound_a_volume(mesh, CGAL::parameters::vertex_point_map(exact_points)));

		}

	}

	for (auto vertex : mesh.vertices()) {
		mesh.point(vertex) = to_kernel(exact_points[vertex]);
	}

}
