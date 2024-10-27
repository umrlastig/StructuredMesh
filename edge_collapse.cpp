#include "edge_collapse.hpp"

#include <list>
#include <limits>
#include <algorithm>
#include <filesystem>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/intersections.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

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

typedef CGAL::Eigen_svd::Vector                             Eigen_vector;
typedef CGAL::Eigen_svd::Matrix                             Eigen_matrix;
typedef CGAL::Quadratic_program<ET>                         Program;
typedef CGAL::Quadratic_program_solution<ET>                Solution;

typedef CGAL::Search_traits_3<Point_set_kernel>             Traits_base;
typedef CGAL::Search_traits_adapter<Point_set::Index, Point_set::Point_map, Traits_base> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>      Neighbor_search;
typedef Neighbor_search::Tree                               Point_tree;

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>        AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                AABB_tree;

void compute_stat(Surface_mesh &mesh, const Ablation_study &ablation, TimerUtils::Timer &timer, K::FT cost) {
	AABB_tree mesh_tree;
	PMP::build_AABB_tree(mesh, mesh_tree);
	CGAL::Cartesian_converter<Point_set_kernel, K> type_converter;

	K::FT min_distance = std::numeric_limits<K::FT>::max();
	K::FT max_distance = std::numeric_limits<K::FT>::min();
	K::FT mean_distance = 0;

	for (const auto &p: ablation.ground_truth_surface_mesh.points()) {
		auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
		if (!std::isnan(location.second[0])) {
			K::FT d = CGAL::sqrt(CGAL::squared_distance(p, PMP::construct_point(location, mesh)));
			if (d < min_distance) min_distance = d;
			if (d > max_distance) max_distance = d;
			mean_distance += d / ablation.ground_truth_surface_mesh.num_vertices();
		}
	}

	AABB_tree ground_truth_surface_mesh_tree;
	PMP::build_AABB_tree(ablation.ground_truth_surface_mesh, ground_truth_surface_mesh_tree);

	K::FT min_distance_r = std::numeric_limits<K::FT>::max();
	K::FT max_distance_r = std::numeric_limits<K::FT>::min();
	K::FT mean_distance_r = 0;

	std::vector<K::Point_3> samples;
	PMP::sample_triangle_mesh(mesh, std::back_inserter(samples), CGAL::parameters::random_seed(0));
	for (const auto &p: samples) {
		auto location = PMP::locate_with_AABB_tree(p, ground_truth_surface_mesh_tree, ablation.ground_truth_surface_mesh);
		if (!std::isnan(location.second[0])) {
			K::FT d = CGAL::sqrt(CGAL::squared_distance(p, PMP::construct_point(location, ablation.ground_truth_surface_mesh)));
			if (d < min_distance_r) min_distance_r = d;
			if (d > max_distance_r) max_distance_r = d;
			mean_distance_r += d / samples.size();
		}
	}

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool created_mesh_label, created_point_in_face;
	boost::tie(mesh_label, created_mesh_label) = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("f:label", LABEL_OTHER);

	Point_set::Property_map<unsigned char> ground_truth_point_cloud_label;
	bool has_ground_truth_point_cloud_label;
	boost::tie(ground_truth_point_cloud_label, has_ground_truth_point_cloud_label) = ablation.ground_truth_point_cloud.property_map<unsigned char>("p:label");

	if (created_mesh_label) {
		// Create point_in_face
		boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points", std::list<Point_set::Index>());

		if(created_point_in_face) {
			for (const auto &ph: ablation.ground_truth_point_cloud) {
				auto p = type_converter(ablation.ground_truth_point_cloud.point(ph));
				auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
				point_in_face[location.first].push_back(ph);
			}
		}

		add_label(mesh, ablation.ground_truth_point_cloud, 0);
	}

	// Check label error
	std::size_t num_wrong_points = 0;
	std::size_t total_num_points = 0;

	if (has_ground_truth_point_cloud_label) {
		for (const auto &ph: ablation.ground_truth_point_cloud) {
			auto p = type_converter(ablation.ground_truth_point_cloud.point(ph));
			auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
			if (mesh_label[location.first] != ground_truth_point_cloud_label[ph]) num_wrong_points += 1;
		}
		total_num_points = ablation.ground_truth_point_cloud.size();
	}

	// Compute semantic border length
	K::FT total_semanctic_contour_length = 0;

	for (const auto &edge: mesh.edges()) {
		if (!mesh.is_border(edge)) {
			auto h = mesh.halfedge(edge);
			if (mesh_label[mesh.face(h)] != mesh_label[mesh.face(mesh.opposite(h))]) {
				total_semanctic_contour_length += CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(h)), mesh.point(mesh.target(h))));
			}
		}
	}

	if (created_mesh_label) {
		mesh.remove_property_map<Surface_mesh::Face_index, unsigned char>(mesh_label);
		if (created_point_in_face) {
			mesh.remove_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>(point_in_face);
		}
	}

	std::ifstream input_file("commande_line.txt");

	if (!std::filesystem::exists("results.csv")) {
		std::ofstream results;
		results.open ("results.csv", std::ios::app);
		results << "mesh, point cloud, l1, l2, l3, l4, l5, l6, l7, c1, c2, c3, c4, ns / cs, subsample, min_point_factor, subdivide, direct_search, border_point, num_mesh_vertices, num_mesh_edges, min_distance, max_distance, mean_distance, min_distance_r, max_distance_r, mean_distance_r, num_wrong_points, total_num_points, total_semanctic_contour_length, time, final_cost\n";
	}

	std::ofstream results;
	results.open ("results.csv", std::ios::app);

	std::string line;
	std::getline(input_file, line);
	while (std::getline(input_file, line) && !line.empty()) {
		size_t equals_pos = line.find('=');
		if (equals_pos != std::string::npos) {
			std::string value = line.substr(equals_pos + 1);
			results << value << ", ";
		}
	}

	results << mesh.number_of_vertices() << ", " << mesh.number_of_edges() << ", " << min_distance << ", " << max_distance << ", " << mean_distance << ", " << min_distance_r << ", " << max_distance_r << ", " << mean_distance_r << ", " << num_wrong_points << ", " << total_num_points << ", " << total_semanctic_contour_length << ", " << timer.getElapsedTime() << ", " << cost << "\n";
}

K::FT get_mean_point_per_area(Surface_mesh &mesh, const Point_set &point_cloud) {
	K::FT total_area = 0;
	for (const auto &face: mesh.faces()) {
		auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
		total_area +=  CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
	}
	if (total_area > 0) return point_cloud.size() / total_area;
	return 0;
}

void add_label(Surface_mesh &mesh, const Point_set &point_cloud, const K::FT min_point_per_area) {
	//Update mesh label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	bool has_mesh_label;
	boost::tie(mesh_label, has_mesh_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("f:label");
	assert(has_mesh_label);

	Point_set::Property_map<unsigned char> point_cloud_label;
	bool has_point_cloud_label;
	boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_point_cloud_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = mesh.property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	std::set<Surface_mesh::Face_index> faces_with_no_label;
	for (const auto &face: mesh.faces()) {
		K::FT min_point = 0;
		if (min_point_per_area > 0) {
			auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
			K::FT face_area = CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
			min_point = min_point_per_area * face_area;
		}
		if (point_in_face[face].size() > 0) {

			if (point_in_face[face].size() > min_point) {

				int face_label[LABELS.size()] = {0};

				for (const auto &point: point_in_face[face]) {
					face_label[point_cloud_label[point]]++;
				}

				auto argmax = std::max_element(face_label, face_label+LABELS.size());
				mesh_label[face] = argmax - face_label;

			} else { // We don't have information about this face.
				mesh_label[face] = LABEL_OTHER;
			}
		} else {
			if (min_point <= 1) { // It's probable that there is no point
				faces_with_no_label.insert(face);
			} else { // We don't have information about this face.
				mesh_label[face] = LABEL_OTHER;
			}
		}
	}
	while(faces_with_no_label.size() > 0) {
		std::list<Surface_mesh::Face_index> face_to_be_removed;
		for (const auto &face: faces_with_no_label) {
			K::FT face_label[LABELS.size()] = {0};
			bool no_neighbor = true;
			for (const auto &he: mesh.halfedges_around_face(mesh.halfedge(face))) {
				if (!mesh.is_border(Surface_mesh::Edge_index(he))) {
					no_neighbor = false;
					if (faces_with_no_label.count(mesh.face(mesh.opposite(he))) == 0) {
						face_label[mesh_label[mesh.face(mesh.opposite(he))]] += CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(he)), mesh.point(mesh.target(he))));
					}
				}
			}
			if (no_neighbor) {
				mesh_label[face] = LABEL_OTHER;
				face_to_be_removed.push_back(face);
			} else {
				auto argmax = std::max_element(face_label, face_label+LABELS.size());
				if (*argmax > 0) {
					mesh_label[face] = argmax - face_label;
					face_to_be_removed.push_back(face);
				}
			}
		}
		for (const auto &face: face_to_be_removed) {
			faces_with_no_label.erase(face);
		}
		if (face_to_be_removed.size() == 0) {
			for (const auto &face_id: faces_with_no_label) {
				mesh_label[face_id] = LABEL_OTHER;
			}
			break;
		} else {
			face_to_be_removed.clear();
		}
	}
}

std::pair<K::Vector_3, K::FT> compute_SVM(std::vector<Point_set::Index> points_for_svm, std::vector<int> y, const Point_set &point_cloud, float * quality = nullptr) {

// TimerUtils::Timer timer;
// timer.start();

	assert(points_for_svm.size() == y.size());

	CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;

	Point_set_kernel::FT c = 10000;

	Program qp (CGAL::EQUAL, true, 0, true, c);

	for (std::size_t i = points_for_svm.size() - 1; i < points_for_svm.size(); i--) {
		auto vr = Point_set_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) * y[i];
		qp.set_d(i, i, CGAL::scalar_product(vr,vr));
		for (std::size_t j = 0; j < i; j++) {
			qp.set_d(i, j, CGAL::scalar_product(vr, Point_set_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[j])) * y[j]));
		}
		qp.set_c(i, -1);
		qp.set_a(i, 0,  y[i]);
	}
	Solution s = CGAL::solve_quadratic_program(qp, ET());
	//std::cerr << "s: " << s << "\n";
	if (!s.solves_quadratic_program(qp)) std::cerr << "ALERT !!!!!!!!!!!!!!!!!!!\n";
	assert (s.solves_quadratic_program(qp));

// timer.pause();
// auto inter = timer.getElapsedTime();
// timer.resume();

	Point_set_kernel::Vector_3 w(CGAL::NULL_VECTOR);
	auto value = s.variable_values_begin();
	for (std::size_t i = 0; i < points_for_svm.size(); i++) {
		w += y[i] * CGAL::to_double(*(value++)) * Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]));
	}

	Point_set_kernel::FT b = 0;
	int count = 0;
	Point_set_kernel::FT min_positive = std::numeric_limits<Point_set_kernel::FT>::max();
	Point_set_kernel::FT max_negative = std::numeric_limits<Point_set_kernel::FT>::lowest();
	for (std::size_t i = 0; i < points_for_svm.size(); i++) {
		Point_set_kernel::FT v = CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
		if (y[i] > 0 && v < min_positive) min_positive = v;
		if (y[i] < 0 && v > max_negative) max_negative = v;
	}
	if (max_negative < min_positive) {
		b += - (max_negative + min_positive) / 2;
		count = 1;
	}
	if (count == 0) {
		value = s.variable_values_begin();
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			double v = CGAL::to_double(*(value++));
			if (v > 0 && v < c) {
				b += y[i] - CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
				count++;
			}
		}
	}
	if (count == 0) {
		value = s.variable_values_begin();
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			if (CGAL::to_double(*(value++)) > 0) {
				b += y[i] - CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
				count++;
			}
		}
	}
	if (count == 0) {
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			b += y[i] - CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
			count++;
		}
	}
	assert(count > 0);
	b /= count;

	if (quality != nullptr) {
		*quality = 0;
		value = s.variable_values_begin();
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			if (CGAL::to_double(*(value++)) > 0) {
				Point_set_kernel::FT v = CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
				if (y[i] > 0 && v < -b) *quality += 1;
				if (y[i] < 0 && v > -b) *quality += 1;
			}
		}
		*quality = 1 - *quality/points_for_svm.size();
	}

	// if (*quality < 0.8) {
	// 	std::cerr << "q: " << *quality << "\n";
	// 	// timer.pause();
	// 	// std::cerr << "t1: " << inter << "\n";
	// 	// std::cerr << "t: " << timer.getElapsedTime() << "\n";
	// 	std::cerr << "w: " << w << "\n";
	// 	std::cerr << "b: " << b << "\n";

	// 	for (std::size_t i = 0; i < points_for_svm.size(); i++) {
	// 		std::cerr << "y" << i << ":\t" << y[i] << "\t";
	// 		std::cerr << "p: " << Point_set_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) << "   \t";
	// 		std::cerr << "s: " << CGAL::to_double(*(s.variable_values_begin() + i)) << "\t";
	// 		std::cerr << "w.p: " << CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) << "  \t";
	// 		std::cerr << "y(w.p + b): " << y[i]*(CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) + b) << "\n";
	// 		if (y[i]*(CGAL::scalar_product(w, Point_set_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) + b) < 0) {
	// 			std::cerr << "Error\n";
	// 		}
	// 	}
	// }

	return std::pair<K::Vector_3, K::FT>(type_converter(w), -type_converter(b));
}

std::vector<Surface_mesh::Face_index> subdivide_face(Surface_mesh& mesh, Surface_mesh::Face_index face, int max_label, const Point_set &point_cloud) {
	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = mesh.property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	std::vector<int> y;
	std::vector<Point_set::Index> points_for_svm;
	y.reserve(point_in_face[face].size());
	points_for_svm.reserve(point_in_face[face].size());
	for (const auto &point: point_in_face[face]) {
		if (label[point] == max_label) {
			points_for_svm.push_back(point);
			y.push_back(1);
		} else {
			points_for_svm.push_back(point);
			y.push_back(-1);
		}
	}

	auto plan = compute_SVM(points_for_svm, y, point_cloud);
	K::Plane_3 cut(CGAL::ORIGIN + plan.second/plan.first.squared_length()*plan.first, plan.first);

	auto h0 = mesh.halfedge(face);
	auto h1 = mesh.next(h0);
	auto h2 = mesh.next(h1);
	auto v0 = mesh.target(h0);
	auto v1 = mesh.target(h1);
	auto v2 = mesh.target(h2);

	auto vp0 = mesh.target(mesh.next(mesh.opposite(h0)));
	auto vp1 = mesh.target(mesh.next(mesh.opposite(h1)));
	auto vp2 = mesh.target(mesh.next(mesh.opposite(h2)));
	if (mesh.is_border(mesh.edge(h0))) vp0 = mesh.null_vertex();
	if (mesh.is_border(mesh.edge(h1))) vp1 = mesh.null_vertex();
	if (mesh.is_border(mesh.edge(h2))) vp2 = mesh.null_vertex();

	auto p0 = mesh.point(v0);
	auto p1 = mesh.point(v1);
	auto p2 = mesh.point(v2);

	Surface_mesh::Vertex_index vm0, vm1, vm2;
	if (CGAL::approximate_angle(cut.orthogonal_vector(), K::Plane_3(p0, p1, p2).orthogonal_vector()) > 10) {
		// the SVM plane and the face are not too parallel
		const auto result0 = CGAL::intersection(cut, K::Segment_3(p2, p0));
		if (result0) {
			if (const Point_3* p = boost::get<Point_3>(&*result0)) {
				if (*p != p0 && *p != p2) {
					vm0 = mesh.add_vertex(*p);
				} else {
					vm0 = mesh.add_vertex(CGAL::midpoint(p2, p0));
				}
			} else {
				vm0 = mesh.add_vertex(CGAL::midpoint(p2, p0));
			}
		} else {
			vm0 = mesh.add_vertex(CGAL::midpoint(p2, p0));
		}
		const auto result1 = CGAL::intersection(cut, K::Segment_3(p0, p1));
		if (result1) {
			if (const Point_3* p = boost::get<Point_3>(&*result1)) {
				if (*p != p0 && *p != p1) {
					vm1 = mesh.add_vertex(*p);
				} else {
					vm1 = mesh.add_vertex(CGAL::midpoint(p0, p1));
				}
			} else {
				vm1 = mesh.add_vertex(CGAL::midpoint(p0, p1));
			}
		} else {
			vm1 = mesh.add_vertex(CGAL::midpoint(p0, p1));
		}
		const auto result2 = CGAL::intersection(cut, K::Segment_3(p1, p2));
		if (result2) {
			if (const Point_3* p = boost::get<Point_3>(&*result2)) {
				if (*p != p1 && *p != p2) {
					vm2 = mesh.add_vertex(*p);
				} else {
					vm2 = mesh.add_vertex(CGAL::midpoint(p1, p2));
				}
			} else {
				vm2 = mesh.add_vertex(CGAL::midpoint(p1, p2));
			}
		} else {
			vm2 = mesh.add_vertex(CGAL::midpoint(p1, p2));
		}
	} else {
		vm0 = mesh.add_vertex(CGAL::midpoint(p2, p0));
		vm1 = mesh.add_vertex(CGAL::midpoint(p0, p1));
		vm2 = mesh.add_vertex(CGAL::midpoint(p1, p2));
	}
	
	K::Point_3 pp0, pp1, pp2;
	if (vp0 != mesh.null_vertex()) pp0 = mesh.point(vp0);
	if (vp1 != mesh.null_vertex()) pp1 = mesh.point(vp1);
	if (vp2 != mesh.null_vertex()) pp2 = mesh.point(vp2);

	if (vp0 != mesh.null_vertex()) CGAL::Euler::remove_face(mesh.opposite(h0), mesh);
	if (vp1 != mesh.null_vertex()) CGAL::Euler::remove_face(mesh.opposite(h1), mesh);
	if (vp2 != mesh.null_vertex()) CGAL::Euler::remove_face(mesh.opposite(h2), mesh);
	CGAL::Euler::remove_face(h0, mesh);

	if (mesh.is_removed(v0)) v0 = mesh.add_vertex(p0);
	if (mesh.is_removed(v1)) v1 = mesh.add_vertex(p1);
	if (mesh.is_removed(v2)) v2 = mesh.add_vertex(p2);
	if (vp0 != mesh.null_vertex() && mesh.is_removed(vp0)) vp0 = mesh.add_vertex(pp0);
	if (vp1 != mesh.null_vertex() && mesh.is_removed(vp1)) vp1 = mesh.add_vertex(pp1);
	if (vp2 != mesh.null_vertex() && mesh.is_removed(vp2)) vp2 = mesh.add_vertex(pp2);

	std::vector<Surface_mesh::Face_index> new_faces;
	new_faces.push_back(mesh.add_face(vm0, vm1, vm2));
	new_faces.push_back(mesh.add_face(v0, vm1, vm0));
	new_faces.push_back(mesh.add_face(v1, vm2, vm1));
	new_faces.push_back(mesh.add_face(v2, vm0, vm2));
	if (vp0 != mesh.null_vertex()) new_faces.push_back(mesh.add_face(v0, vm0, vp0));
	if (vp0 != mesh.null_vertex()) new_faces.push_back(mesh.add_face(v2, vp0, vm0));
	if (vp1 != mesh.null_vertex()) new_faces.push_back(mesh.add_face(v1, vm1, vp1));
	if (vp1 != mesh.null_vertex()) new_faces.push_back(mesh.add_face(v0, vp1, vm1));
	if (vp2 != mesh.null_vertex()) new_faces.push_back(mesh.add_face(v2, vm2, vp2));
	if (vp2 != mesh.null_vertex()) new_faces.push_back(mesh.add_face(v1, vp2, vm2));

	return new_faces;
}

std::list<std::pair <K::Vector_3, K::FT>> volume_preservation_and_optimisation (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<std::pair <K::Vector_3, K::FT>> result;

	SMS::Edge_profile<Surface_mesh>::Triangle_vector triangles = profile.triangles();
	for (const auto &triangle: triangles) {
		Point_3 p0 = get(profile.vertex_point_map(),triangle.v0);
		Point_3 p1 = get(profile.vertex_point_map(),triangle.v1);
		Point_3 p2 = get(profile.vertex_point_map(),triangle.v2);

		if (!CGAL::collinear (p0, p1, p2)) {
/*if (profile.v0().idx() == 51431 && profile.v1().idx() == 175066) {
	std::cerr << "triangle: " << p0 << ", " << p1 << ", " << p2 << "\n";
	std::cerr << "normal: " << CGAL::normal(p0, p1, p2) << "\n";
	std::cerr << "determinant: " << CGAL::determinant(K::Vector_3(Point_3(CGAL::ORIGIN), p0), K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p2)) << "\n";
}*/

			K::Vector_3 n = CGAL::normal(p0, p1, p2) / 2;
			K::FT det = CGAL::scalar_product (n, K::Vector_3(Point_3(CGAL::ORIGIN), p0));

			result.push_back(std::make_pair(n, det));
		}
	}

	return result;
}

std::list<std::pair <K::Vector_3, K::Vector_3>> boundary_preservation_and_optimisation (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<std::pair <K::Vector_3, K::Vector_3>> result;

	for (const auto &edge: profile.border_edges()) {
		Point_3 p0 = get(profile.vertex_point_map(), profile.surface_mesh().source(edge));
		Point_3 p1 = get(profile.vertex_point_map(), profile.surface_mesh().target(edge));

		//std::cerr << "edge: " << p0 << ", " << p1 << "\n";

		auto e1 = K::Vector_3(p0, p1);
		auto e2 = CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0));

		/*std::cerr << "e1: " << K::Vector_3(p0, p1) << "\n";
		std::cerr << "e2: " << CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0)) << "\n";*/

		result.push_back(std::make_pair(e1, e2));
	}

	return result;
}

std::list<K::Vector_3> triangle_shape_optimization (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<K::Vector_3> result;

	for (const auto &v: profile.link()) {
		result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), v)));
	}

	return result;

}

class Label_simple_optimization {
	const std::map<Point_set::Index, K::Point_3> &cloud_point_in_plane;
	const std::set<Point_set::Index> &points_in_faces;
	const std::map<Surface_mesh::Vertex_index, K::Point_3> &point_in_plane;
	const std::vector<Surface_mesh::Halfedge_index> &faces_border_halfedge;
	const Surface_mesh &mesh;
	const Point_set::Property_map<unsigned char> label;
	std::map<K::Point_3, K::FT> known_energy;
	std::set<std::pair<int,int>> grid_search_visit;

	public:
		Label_simple_optimization(const std::map<Point_set::Index, K::Point_3> &cloud_point_in_plane, const std::set<Point_set::Index> &points_in_faces, const std::map<Surface_mesh::Vertex_index, K::Point_3> &point_in_plane, const std::vector<Surface_mesh::Halfedge_index> &faces_border_halfedge, const Surface_mesh &mesh, const Point_set::Property_map<unsigned char> label) : cloud_point_in_plane(cloud_point_in_plane), points_in_faces(points_in_faces), point_in_plane(point_in_plane), faces_border_halfedge(faces_border_halfedge), mesh(mesh), label(label) {}

		/*K::FT energy(K::Point_3 middle) {
			if (auto search = known_energy.find(middle); search != known_energy.end()) {
				// std::cerr << "es: " << middle << " " << search->second << "\n";
				return search->second;
			} else {

				for (std::size_t j = 0; j < faces_border_halfedge.size(); j++) {
					if (CGAL::scalar_product(CGAL::orthogonal_vector (middle, mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j]))), CGAL::orthogonal_vector (mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), mesh.point(mesh.target(mesh.next(faces_border_halfedge[j]))))) < 0) {

						// std::cerr << "t_a: " << K::Triangle_3(mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), middle) << "\n";
						// std::cerr << "t_b: " << K::Triangle_3(mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), mesh.point(mesh.target(mesh.next(faces_border_halfedge[j])))) << "\n";
						
						known_energy[middle] = std::numeric_limits<K::FT>::max();
						return std::numeric_limits<K::FT>::max();
					}
				}

				CGAL::Eigen_matrix< K::FT, Eigen::Dynamic, Eigen::Dynamic > Dpt(points_in_faces.size(), faces_border_halfedge.size());
				for (std::size_t j = 0; j < faces_border_halfedge.size(); j++) {
					K::Triangle_3 triangle (middle, point_in_plane.at(mesh.source(faces_border_halfedge[j])), point_in_plane.at(mesh.target(faces_border_halfedge[j])));
					std::size_t i = 0;
					for (const auto &v: points_in_faces) {
						auto d = CGAL::sqrt(CGAL::squared_distance(triangle, cloud_point_in_plane.at(v)));
						Dpt.set(i++,j, d);
					}
				}

				// std::cerr << "Dpt\n" << Dpt << "\n";
				
				// K::FT tho = 100;
				// auto FDpt = (Dpt.array() * (- tho)).exp();

				K::FT tho = 0.0005;
				auto FDpt = (Dpt.array() < tho).select(- Dpt.array() / tho + K::FT(1), K::FT(0));

				// std::cerr << "FDpt\n" << FDpt << "\n";

				auto sums = FDpt.rowwise().sum();
				auto Apt = FDpt.colwise() / (sums + (sums == 0).cast<K::FT>());

				// std::cerr << "Apt\n" << Apt << "\n";

				std::map<unsigned char, std::size_t> classes;
				for (const auto &v: points_in_faces) {
					if(classes.count(label[v]) == 0) {
						classes[label[v]] = classes.size();
					}
				}
				Eigen::SparseMatrix<K::FT> Ipc(points_in_faces.size(), classes.size());
				Ipc.reserve(points_in_faces.size());
				std::size_t i = 0;
				for (const auto &v: points_in_faces) {
					Ipc.insert(i++, classes[label[v]]) = 1;
				}

				// std::cerr << "Ipc\n" << Ipc << "\n";

				auto sums2 = Apt.colwise().sum();
				auto Pct = (Ipc.transpose() * Apt.matrix()).array().rowwise() / (sums2 + (sums2 == 0).cast<K::FT>());

				// std::cerr << "Pct\n" << Pct << "\n";

				auto E = (Pct.array() > 0).select(- Pct * Pct.log(), Pct);

				auto e = E.sum();

				known_energy[middle] = e;

				// std::cerr << "e: " << middle << " " << e << "\n";

				return e;
			}
		}*/

		K::FT energy(K::Point_3 middle) {
			if (auto search = known_energy.find(middle); search != known_energy.end()) {
				// std::cerr << "es: " << middle << " " << search->second << "\n";
				return search->second;
			} else {

				// std::cerr << "SIZE: " << faces_border_halfedge.size() << "\n";

				for (std::size_t j = 0; j < faces_border_halfedge.size(); j++) {
					
					// K::Triangle_3 t_a (mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), middle);
					// K::Triangle_3 t_b (mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), mesh.point(mesh.target(mesh.next(faces_border_halfedge[j]))));
					// std::cerr << "t_a: " << t_a << " : " << CGAL::orthogonal_vector (mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), middle) << "\n";
					// std::cerr << "t_b: " << t_b << " : " << CGAL::orthogonal_vector (mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), mesh.point(mesh.target(mesh.next(faces_border_halfedge[j])))) << "\n";

					if (CGAL::scalar_product(CGAL::orthogonal_vector (middle, mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j]))), CGAL::orthogonal_vector (mesh.point(mesh.source(faces_border_halfedge[j])), mesh.point(mesh.target(faces_border_halfedge[j])), mesh.point(mesh.target(mesh.next(faces_border_halfedge[j]))))) < 0) {

						known_energy[middle] = std::numeric_limits<K::FT>::max();
						return std::numeric_limits<K::FT>::max();
					}
				}

				CGAL::Eigen_matrix< K::FT, Eigen::Dynamic, Eigen::Dynamic > Dpt(points_in_faces.size(), faces_border_halfedge.size());
				for (std::size_t j = 0; j < faces_border_halfedge.size(); j++) {
					K::Triangle_3 triangle (middle, point_in_plane.at(mesh.source(faces_border_halfedge[j])), point_in_plane.at(mesh.target(faces_border_halfedge[j])));
					std::size_t i = 0;
					for (const auto &v: points_in_faces) {
						auto d = CGAL::sqrt(CGAL::squared_distance(triangle, cloud_point_in_plane.at(v)));
						Dpt.set(i++,j, d);
					}
				}

				// std::cerr << "Dpt\n" << Dpt << "\n";
				
				K::FT tho = 0.01;
				K::FT alpha = 0.2;
				auto FDpt = (Dpt.array() < tho).select((Dpt.array() < tho/10).select(Dpt.array() * (alpha - 1) * K::FT(10) / tho + K::FT(1), alpha * K::FT(10/9) * (K::FT(1) - Dpt.array() / tho)), K::FT(0));

				// std::cerr << "FDpt\n" << FDpt << "\n";

				auto sums = FDpt.rowwise().sum();
				auto Apt = FDpt.colwise() / (sums + (sums == 0).cast<K::FT>());

				// std::cerr << "Apt\n" << Apt << "\n";
				
				std::map<unsigned char, std::size_t> classes;
				for (const auto &v: points_in_faces) {
					if(classes.count(label[v]) == 0) {
						classes[label[v]] = classes.size();
					}
				}
				Eigen::SparseMatrix<K::FT> Ipc(points_in_faces.size(), classes.size());
				Ipc.reserve(points_in_faces.size());
				std::size_t i = 0;
				for (const auto &v: points_in_faces) {
					Ipc.insert(i++, classes[label[v]]) = 1;
				}

				// std::cerr << "Ipc\n" << Ipc << "\n";

				auto Pct = Ipc.transpose() * Apt.matrix();

				// std::cerr << "Pct\n" << Pct << "\n";

				auto Mct = Pct.array().colwise().maxCoeff();

				auto e = Pct.sum() - Mct.sum();

				known_energy[middle] = e;

				// std::cerr << "e: " << middle << " " << e << "\n";

				return e;
			}
		}

		std::pair<K::Point_3, K::FT> linear_search(K::Point_3 middle, K::Vector_3 vec1) {
			K::FT pas = 1;
			K::FT pos = 0;
			K::FT e = energy(middle);
			// if (abs(middle.x() - 3.81392) < 0.0001) std::cerr << "p: " << middle << ": " << e << "\n";
			while (e == std::numeric_limits<K::FT>::max() && (pas > 0.01 || pas < -0.01)) {
				// std::cerr << "l\t" <<  pas << "\t" << pos << "\n";
				if (pas > 0) {
					if (pos < 1) {
						pos += pas;
					} else {
						pas = - pas;
						pos = pas;
					}
				} else {
					if (pos > -1) {
						pos += pas;
					} else {
						pas = - pas / 2;
						pos = pas;
					}
				}
				e = energy(middle + pos * vec1);
				// if (abs(middle.x() - 3.81392) < 0.0001) std::cerr << "p: " << middle + pos * vec1 << ": " << e << "\n";
			}
			return std::make_pair(middle + pos * vec1, CGAL::abs(pas));
		}

		void grid_search(std::vector<K::Point_3> &results, K::FT *r_energy, K::Vector_3 &vec1, K::Vector_3 &vec2, int pos1, int pos2) {
			// std::cerr << "g\t" << pos1 << "\t" << pos2 << "\n";
			// for (const auto &r: grid_search_visit) {
			// 	std::cerr << "\t\t" << r.first << "\t" << r.second << "\n";
			// }
			if (grid_search_visit.count(std::make_pair(pos1, pos2)) == 0) {
				grid_search_visit.insert(std::make_pair(pos1, pos2));
				auto e = energy(results[0] + pos1*vec1 + pos2*vec2);
				if (e < *r_energy) {
					results[1] = results[0] + pos1*vec1 + pos2*vec2;
					*r_energy = e;
				}
				// std::cerr << "e: " << e << "\n";
				if (e < std::numeric_limits<K::FT>::max() && abs(pos1) < 10 && abs(pos2) < 10) {
					grid_search(results, r_energy, vec1, vec2, pos1 + 1, pos2);
					grid_search(results, r_energy, vec1, vec2, pos1 - 1, pos2);
					grid_search(results, r_energy, vec1, vec2, pos1, pos2 + 1);
					grid_search(results, r_energy, vec1, vec2, pos1, pos2 - 1);
				}
			}
		}
};


std::pair<K::Point_3, std::pair<K::Vector_3, K::Vector_3>> best_position(const SMS::Edge_profile<Surface_mesh>& profile, const Point_set &point_cloud, const std::set<Point_set::Index> &points_in_faces) {

	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	K::Vector_3 ortho_plane (CGAL::NULL_VECTOR);
	if (profile.left_face_exists()) {
		ortho_plane += CGAL::orthogonal_vector(profile.surface_mesh().point(profile.v0()), profile.surface_mesh().point(profile.v1()), profile.surface_mesh().point(profile.vL()));
	}
	if (profile.right_face_exists()) {
		ortho_plane += CGAL::orthogonal_vector(profile.surface_mesh().point(profile.v1()), profile.surface_mesh().point(profile.v0()), profile.surface_mesh().point(profile.vR()));
	}
	K::Plane_3 plane(profile.surface_mesh().point(profile.v0()), ortho_plane);

	CGAL::Cartesian_converter<Point_set_kernel, K> type_converter;
	std::map<Point_set::Index, K::Point_3> cloud_point_in_plane;
	for (const auto &v: points_in_faces) {
		cloud_point_in_plane[v] = plane.projection(type_converter(point_cloud.point(v)));
	}

	std::map<Surface_mesh::Vertex_index, K::Point_3> point_in_plane;
	for (const auto &v: profile.link()) {
		point_in_plane[v] = plane.projection(profile.surface_mesh().point(v));
	}

	std::vector<Surface_mesh::Halfedge_index> faces_border_halfedge;
	for (std::size_t face_id = 0; face_id < profile.link().size(); face_id++) {
		auto he = profile.surface_mesh().halfedge(profile.link()[face_id], profile.link()[(face_id + 1) % profile.link().size()]);
		if (he != Surface_mesh::null_halfedge() && !profile.surface_mesh().is_border(he)) {
			faces_border_halfedge.push_back(he);
		}
	}

	Label_simple_optimization optim(cloud_point_in_plane, points_in_faces, point_in_plane, faces_border_halfedge, profile.surface_mesh(), label);

	std::vector<K::Point_3> results;
	results.reserve(2);
	K::Vector_3 vec1 (profile.surface_mesh().point(profile.v0()), profile.surface_mesh().point(profile.v1()));
	vec1 /= 2;
	K::FT pas = CGAL::sqrt(vec1.squared_length());
	K::Vector_3 vec2 = CGAL::cross_product(vec1, ortho_plane);
	vec2 *= pas / CGAL::sqrt(vec2.squared_length());

	// Linear seach
	auto r = optim.linear_search(CGAL::midpoint(profile.surface_mesh().point(profile.v0()), profile.surface_mesh().point(profile.v1())), vec1);

	// if (profile.v0().idx() == 7415 && profile.v1().idx() == 6860) std::cerr << "linear search: " << r.first << " (" << r.second << ")\n";

	if (r.second == 0.01) {
		return std::pair<K::Point_3, std::pair<K::Vector_3, K::Vector_3>>(K::Point_3(CGAL::ORIGIN), std::pair<K::Vector_3, K::Vector_3>(CGAL::NULL_VECTOR, CGAL::NULL_VECTOR));
	}
	results.push_back(r.first);

	// Grid search
	vec1 *= r.second;
	vec2 *= r.second;
	pas *= r.second;

	K::FT r_energy = std::numeric_limits<K::FT>::max();
	// std::cerr << "-----------------------------\n";
	optim.grid_search(results, &r_energy, vec1, vec2, 0, 0);
	/*if (optim.grid_search_visit.size() < 10) {
		optim.grid_search_visit.clear();
		vec1 /= 2;
		vec2 /= 2;
		pas /= 2;
		optim.grid_search(results, &r_energy, vec1, vec2, 0, 0);
	}*/

	if (optim.energy(r.first) <= optim.energy(results[1])) {
		results[1] = r.first;
	}

	// Cross search	
	K::FT pos1 = 0;
	K::FT pos2 = 0;
	K::FT var = 1;

	while(var > 1e-4 && results.size() < 20) {
		K::FT values[5];
		values[0] = optim.energy(results[1] + pos1 * vec1 + pos2 * vec2); // center
		values[1] = optim.energy(results[1] + (pos1 - pas) * vec1 + pos2 * vec2); // left
		values[2] = optim.energy(results[1] + (pos1 + pas) * vec1 + pos2 * vec2); // right
		values[3] = optim.energy(results[1] + pos1 * vec1 + (pos2 - pas) * vec2); // bottom
		values[4] = optim.energy(results[1] + pos1 * vec1 + (pos2 + pas) * vec2); // top

		auto argmin = std::min_element(values, values + 5);
		int p = argmin - values;

		if (p == 0) {
			pas /= 2;
			auto argmin = std::min_element(values + 1, values + 5);
			var = *argmin - values[0];
		} else if (p == 1) {
			pos1 -= pas;
			var = values[0] - values[1];
		} else if (p == 2) {
			pos1 += pas;
			var = values[0] - values[2];
		} else if (p == 3) {
			pos2 -= pas;
			var = values[0] - values[3];
		} else {
			pos2 += pas;
			var = values[0] - values[4];
		}

		results.push_back(results[1] + pos1 * vec1 + pos2 * vec2);
	}

	return std::pair<K::Point_3, std::pair<K::Vector_3, K::Vector_3>>(results.back(), std::pair<K::Vector_3, K::Vector_3>(vec1, vec2));


	// std::list<Surface_mesh::Face_index> selected_face;
	// for (const auto &face: profile.triangles()) {
	// 	auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
	// 	selected_face.push_back(fh);
	// }
	// CGAL::Face_filtered_graph<Surface_mesh> filtered_sm(profile.surface_mesh(), selected_face);

	// Surface_mesh filtered_mesh;
	// CGAL::copy_face_graph(filtered_sm, filtered_mesh);
	// AABB_tree filtered_sm_tree;
	// PMP::build_AABB_tree(filtered_mesh, filtered_sm_tree);

	// auto location = PMP::locate_with_AABB_tree(K::Ray_3(results.back(), ortho_plane), filtered_sm_tree, filtered_mesh);
	// if (location.first != Surface_mesh::Face_index()) {
	// 	auto p1 = PMP::construct_point(location, profile.surface_mesh());
	// 	auto location2 = PMP::locate_with_AABB_tree(K::Ray_3(results.back(), -ortho_plane), filtered_sm_tree, filtered_mesh);
	// 	if (location2.first != Surface_mesh::Face_index()) {
	// 		auto p2 = PMP::construct_point(location2, profile.surface_mesh());
	// 		if (CGAL::squared_distance(p1, results.back()) < CGAL::squared_distance(p2, results.back())) {
	// 			return p1;
	// 		} else {
	// 			return p2;
	// 		}
	// 	} else {
	// 		return p1;
	// 	}
	// } else {
	// 	auto location = PMP::locate_with_AABB_tree(K::Ray_3(results.back(), -ortho_plane), filtered_sm_tree, filtered_mesh);
	// 	if (location.first != Surface_mesh::Face_index()) {
	// 		return PMP::construct_point(location, profile.surface_mesh());
	// 	} else {
	// 		return results.back();
	// 	}
	// }

	/*if (profile.v0().idx() == 7415 && profile.v1().idx() == 6860) std::cerr << "final energy: " << optim.energy(results.back()) << "\n";

	if (profile.v0().idx() == 7415 && profile.v1().idx() == 6860) { // if (optim.energy(results[0]) > 1e15) {
		std::cerr << "v0: " << profile.v0() << "\n";
		std::cerr << "v1: " << profile.v1() << "\n";
		std::cerr << "p0: " << profile.p0() << "\n";
		std::cerr << "p1: " << profile.p1() << "\n";
		std::cerr << "points_in_faces.size: " << points_in_faces.size() << "\n";
		std::cerr << "p:" << results.back() << "\n";
		std::cerr << "---\n";
		auto mesh = profile.surface_mesh();
		for (std::size_t j = 0; j < faces_border_halfedge.size(); j++) {
			std::cerr <<  point_in_plane.at(mesh.source(faces_border_halfedge[j])) << " - " << point_in_plane.at(mesh.target(faces_border_halfedge[j])) << "\n";
		}
		std::cerr << "---\n";
		for (std::size_t j = 0; j < faces_border_halfedge.size(); j++) {
			auto p = mesh.point(mesh.source(faces_border_halfedge[j]));
			std::cerr <<  "\t[" << p.x() << ", " << p.y() << ", " << p.z() << "],\n";
		}
		std::cerr << "---\n";
		for (const auto &v: points_in_faces) {
			auto p = cloud_point_in_plane[v];
			std::cerr <<  "\t[" << p.x() << ", " << p.y() << ", " << p.z() << "],\n";
		}
		std::cerr << "---\n";
		for (const auto &v: points_in_faces) std::cerr << "int: " << int(label[v]) << "\n";
		std::cerr << "------------------------------------------------------------------------------------------\n";
	}
	// 	for (const auto &r: results) std::cerr << "r: " << r << ": " << optim.energy(r) << "\n";
	// 	std::cerr << "-----\n";
	// } else {
	// 	std::cerr << "---\n";
	// 	for (const auto &r: results) std::cerr << "r: " << r << ": " << optim.energy(r) << "\n";
	// 	std::cerr << "-----\n";
	// }

	if (profile.v0().idx() == 7415 && profile.v1().idx() == 6860) { // Output point cloud and division
		if (true || points_in_faces.size() > 10) {

			auto p = results.back();

			Surface_mesh output_mesh;

			bool created;
			Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> red;
			Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> green;
			Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> blue;
			boost::tie(red, created) = output_mesh.add_property_map<Surface_mesh::Vertex_index, unsigned char>("red",0);
			assert(created);
			boost::tie(green, created) = output_mesh.add_property_map<Surface_mesh::Vertex_index, unsigned char>("green",0);
			assert(created);
			boost::tie(blue, created) = output_mesh.add_property_map<Surface_mesh::Vertex_index, unsigned char>("blue",0);
			assert(created);

			CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;

			for (const auto &point: points_in_faces) {
				auto v = output_mesh.add_vertex(type_converter(point_cloud.point(point)));
				red[v] = LABELS.at(label[point]).red;
				green[v] = LABELS.at(label[point]).green;
				blue[v] = LABELS.at(label[point]).blue;
			}

			auto v = output_mesh.add_vertex(p);

			for (std::size_t face_id = 0; face_id < profile.link().size(); face_id++) {
				auto he = profile.surface_mesh().halfedge(profile.link()[face_id], profile.link()[(face_id + 1) % profile.link().size()]);
				if (he != Surface_mesh::null_halfedge() && !profile.surface_mesh().is_border(he)) {
					auto v1 = output_mesh.add_vertex(profile.surface_mesh().point(profile.surface_mesh().source(he)));
					auto v2 = output_mesh.add_vertex(profile.surface_mesh().point(profile.surface_mesh().target(he)));
					
					output_mesh.add_face(v1,v2,v);
				}
			}

			std::stringstream output_mesh_name;
			output_mesh_name << "optim_" << optim.energy(p) << "_" << profile.v0() << "_" << profile.v1() << "_" << profile.p0() << "_" << profile.p1() << "_" << points_in_faces.size() << "_p:" << p << ".ply";
			std::ofstream mesh_ofile (output_mesh_name.str().c_str(), std::ios_base::binary);
			CGAL::IO::set_binary_mode (mesh_ofile);
			CGAL::IO::write_PLY (mesh_ofile, output_mesh);
			mesh_ofile.close();

		}
	}*/

}

std::list<std::pair <K::Vector_3, K::FT>> label_preservation (const SMS::Edge_profile<Surface_mesh>& profile, const Point_set &point_cloud, const Ablation_study &ablation) {

TimerUtils::Timer timer;
timer.start();

	std::list<std::pair <K::Vector_3, K::FT>> result;

	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	Point_set::Property_map<bool> isborder;
	if (ablation.border_point) {
		bool has_isborder;
		boost::tie(isborder, has_isborder) = point_cloud.property_map<bool>("p:isborder");
		assert(has_isborder);
	}

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	std::set<Point_set::Index> points_in_faces;
	int count_collapse_label[LABELS.size()] = {0};
	for (const auto &face: profile.triangles()) {
		auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
		if (point_in_face[fh].size() > 0) {
			int count_face_label[LABELS.size()] = {0};
			for (const auto &ph: point_in_face[fh]) {
				if (!ablation.border_point || isborder[ph]) {
					points_in_faces.insert(ph);
				}
				count_face_label[label[ph]]++;
			}
			auto argmax = std::max_element(count_face_label, count_face_label+LABELS.size()); // face label
			count_collapse_label[argmax - count_face_label]++;
		}
	}

	int count_diff_label = 0;
	for (std::size_t i = 0; i < LABELS.size(); i++) {
		if (count_collapse_label[i] > 0) count_diff_label++;
	}

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();
	CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;

	if (count_diff_label < 2) {
		return result;
	} else if (count_diff_label == 2) {

		unsigned char label1, label2;
		bool first = true;
		for (std::size_t i = 0; i < LABELS.size(); i++) {
			if (count_collapse_label[i] > 0) {
				if (first) {
					label1 = i;
					first = false;
				} else {
					label2 = i;
				}
			}
		}

		std::vector<int> y;
		std::vector<Point_set::Index> points_for_svm;
		y.reserve(points_in_faces.size());
		points_for_svm.reserve(points_in_faces.size());
		bool has_label1 = false, has_label2 = false; 
		for (const auto &point: points_in_faces) {
			if (label[point] == label1) {
				points_for_svm.push_back(point);
				y.push_back(1);
				has_label1 = true;
			} else if (label[point] == label2) {
				points_for_svm.push_back(point);
				y.push_back(-1);
				has_label2 = true;
			}
		}
		if (!has_label1 || !has_label2) {
			// this is strange, why is it need ?
			//add no border point to points_in_faces and restart SVM
			points_in_faces.clear();
			for (const auto &face: profile.triangles()) {
				auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
				if (point_in_face[fh].size() > 0) {
					points_in_faces.insert(point_in_face[fh].begin(), point_in_face[fh].end());
				}
			}
			y.clear();
			points_for_svm.clear();
			y.reserve(points_in_faces.size());
			points_for_svm.reserve(points_in_faces.size());
			for (const auto &point: points_in_faces) {
				if (label[point] == label1) {
					points_for_svm.push_back(point);
					y.push_back(1);
				} else if (label[point] == label2) {
					points_for_svm.push_back(point);
					y.push_back(-1);
				}
			}
		}

		float quality;
		auto plan = compute_SVM(points_for_svm, y, point_cloud, &quality);

		if(quality > 0.75) result.push_back(std::pair<K::Vector_3, K::FT>(plan.first*squared_length, plan.second*squared_length));
		// else std::cerr << "quality: " << quality << "\n";

		// if (profile.v0().idx() == 7415 && profile.v1().idx() == 6860) { // Output point cloud and plane
		// 	if (true || points_in_faces.size() > 10) {

		// 		Surface_mesh output_mesh;

		// 		bool created;
		// 		Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> red;
		// 		Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> green;
		// 		Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> blue;
		// 		boost::tie(red, created) = output_mesh.add_property_map<Surface_mesh::Vertex_index, unsigned char>("red",0);
		// 		assert(created);
		// 		boost::tie(green, created) = output_mesh.add_property_map<Surface_mesh::Vertex_index, unsigned char>("green",0);
		// 		assert(created);
		// 		boost::tie(blue, created) = output_mesh.add_property_map<Surface_mesh::Vertex_index, unsigned char>("blue",0);
		// 		assert(created);

		// 		CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;

		// 		for (const auto &point: points_in_faces) {
		// 			auto v = output_mesh.add_vertex(type_converter(point_cloud.point(point)));
		// 			red[v] = LABELS.at(label[point]).red;
		// 			green[v] = LABELS.at(label[point]).green;
		// 			blue[v] = LABELS.at(label[point]).blue;
		// 		}

		// 		auto n = plan.first;
		// 		auto d = plan.second;

		// 		K::Plane_3 plan2 (n.x(), n.y(), n.z(), -d);
		// 		auto orig = plan2.projection(type_converter(point_cloud.point(*points_in_faces.begin())));

		// 		auto v1 = output_mesh.add_vertex(orig - 2*plan2.base1() - 2*plan2.base2());
		// 		auto v2 = output_mesh.add_vertex(orig + 4*plan2.base1());
		// 		auto v3 = output_mesh.add_vertex(orig + 4*plan2.base2());
		// 		output_mesh.add_face(v1,v2,v3);

		// 		std::stringstream output_mesh_name;
		// 		output_mesh_name << "svm_(" << profile.p0() << ")-(" << profile.p1() << "), w: " << n << ", b: " << d << ".ply";
		// 		std::ofstream mesh_ofile (output_mesh_name.str().c_str(), std::ios_base::binary);
		// 		CGAL::IO::set_binary_mode (mesh_ofile);
		// 		CGAL::IO::write_PLY (mesh_ofile, output_mesh);
		// 		mesh_ofile.close();

		// 	}
		// }

	} else {

		for (std::size_t i_label = 0; i_label < LABELS.size(); i_label++) {
			if (count_collapse_label[i_label] > 0) {

				std::vector<int> y;
				std::vector<Point_set::Index> points_for_svm;
				y.reserve(points_in_faces.size());
				points_for_svm.reserve(points_in_faces.size());
				bool has_i_label = false, has_other_label = false;
				for (const auto &point: points_in_faces) {
					if (label[point] == i_label) {
						points_for_svm.push_back(point);
						y.push_back(1);
						has_i_label = true;
					} else {
						points_for_svm.push_back(point);
						y.push_back(-1);
						has_other_label = true;
					}
				}

				if (has_i_label && has_other_label) {
					float quality;
					auto plan = compute_SVM(points_for_svm, y, point_cloud, &quality);
					if(quality > 0.75) result.push_back(std::pair<K::Vector_3, K::FT>(plan.first*squared_length, plan.second*squared_length));
					// else std::cerr << "quality: " << quality << "\n";
				}
			}
		}
	}

	if (ablation.direct_search && result.size() == 0) {
		//SVM is poor quality
		points_in_faces.clear();
		for (const auto &face: profile.triangles()) {
			auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
			if (point_in_face[fh].size() > 0) {
				points_in_faces.insert(point_in_face[fh].begin(), point_in_face[fh].end());
			}
		}

		auto p = best_position(profile, point_cloud, points_in_faces);
		if (p.first != CGAL::ORIGIN) {
			p.second.first /= CGAL::sqrt(p.second.first.squared_length());
			p.second.second /= CGAL::sqrt(p.second.second.squared_length());
			auto b1 = CGAL::scalar_product(K::Vector_3(CGAL::ORIGIN, p.first), p.second.first);
			auto b2 = CGAL::scalar_product(K::Vector_3(CGAL::ORIGIN, p.first), p.second.second);
			result.push_back(std::pair<K::Vector_3, K::FT>(p.second.first*squared_length, b1*squared_length));
			result.push_back(std::pair<K::Vector_3, K::FT>(p.second.second*squared_length, b2*squared_length));
			// result.push_back(std::pair<K::Vector_3, K::FT>(K::Vector_3(1,0,0)*squared_length, p.first.x()*squared_length));
			// result.push_back(std::pair<K::Vector_3, K::FT>(K::Vector_3(0,1,0)*squared_length, p.first.y()*squared_length));
			// result.push_back(std::pair<K::Vector_3, K::FT>(K::Vector_3(0,0,1)*squared_length, p.first.z()*squared_length));
		} /*else {
			std::cerr << "--------------------------------------------------\n";
			std::cerr << "v0: " << profile.surface_mesh().point(profile.v0()) << "\n";
			std::cerr << "v1: " << profile.surface_mesh().point(profile.v1()) << "\n";
			std::cerr << "------- point in face\n";
			for (const auto &ph: points_in_faces) {
				std::cerr << "point_cloud.point: " << point_cloud.point(ph) << "\n";
			}
			std::cerr << "------- link\n";
			for (std::size_t face_id = 0; face_id < profile.link().size(); face_id++) {
				auto he = profile.surface_mesh().halfedge(profile.link()[face_id], profile.link()[(face_id + 1) % profile.link().size()]);
				if (he != Surface_mesh::null_halfedge() && !profile.surface_mesh().is_border(he)) {
					std::cerr << "he: " << profile.surface_mesh().point(profile.surface_mesh().source(he)) << " - " << profile.surface_mesh().point(profile.surface_mesh().target(he)) << "\n";
				}
			}
		}*/
	}


timer.pause();
// std::cerr << "timer.getElapsedTime: " << timer.getElapsedTime() << "\n";
	return result;

}

std::list<K::Vector_3> semantic_border_optimization (const SMS::Edge_profile<Surface_mesh>& profile, const Point_set &point_cloud) {

	std::list<K::Vector_3> result;

	Point_set::Property_map<unsigned char> point_cloud_label;
	bool has_point_cloud_label;
	boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_point_cloud_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	bool has_mesh_label;
	boost::tie(mesh_label, has_mesh_label) = profile.surface_mesh().property_map<Surface_mesh::Face_index, unsigned char>("f:label");
	assert(has_mesh_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	for (const auto &h: profile.surface_mesh().halfedges_around_target(profile.v1_v0())) {
		if (h != profile.v1_v0() && h != profile.vL_v0() && !profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
			if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))] && point_in_face[profile.surface_mesh().face(h)].size() > 0 && point_in_face[profile.surface_mesh().face(profile.surface_mesh().opposite(h))].size() > 0) {
				result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), profile.surface_mesh().source(h))));
			}
		}
	}

	for (const auto &h: profile.surface_mesh().halfedges_around_target(profile.v0_v1())) {
		if (h != profile.v0_v1() && h != profile.vR_v1() && !profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
			if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))] && point_in_face[profile.surface_mesh().face(h)].size() > 0 && point_in_face[profile.surface_mesh().face(profile.surface_mesh().opposite(h))].size() > 0) {
				result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), profile.surface_mesh().source(h))));
			}
		}
	}

	/*if (result.size() > 0) {
		std::cerr << "result.size: " << result.size() << "\n";
		for (const auto &face: profile.triangles()) {
			auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
			std::cerr << "fh: " << fh << ": " << int(mesh_label[fh]) << "\n";
		}
	}*/

	return result;

}

LindstromTurk_param::LindstromTurk_param(
	float volume_preservation,
	float boundary_preservation,
	float volume_optimisation,
	float boundary_optimization,
	float triangle_shape_optimization,
	float label_preservation,
	float semantic_border_optimization) :
	volume_preservation(volume_preservation),
	boundary_preservation(boundary_preservation),
	volume_optimisation(volume_optimisation),
	boundary_optimization(boundary_optimization),
	triangle_shape_optimization(triangle_shape_optimization),
	label_preservation(label_preservation),
	semantic_border_optimization(semantic_border_optimization) {}

Ablation_study::Ablation_study(
	bool subdivide,
	bool direct_search,
	bool border_point,
	bool step_mesh) :
	subdivide(subdivide), direct_search(direct_search), border_point(border_point), step_mesh(step_mesh) {}

Custom_placement::Custom_placement (const LindstromTurk_param &params, Surface_mesh &mesh, const Point_set &point_cloud, const Ablation_study ablation) : params(params), point_cloud(point_cloud), ablation(ablation) {
	bool created_collapse_datas;
	boost::tie(collapse_datas, created_collapse_datas) = mesh.add_property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");
}


boost::optional<SMS::Edge_profile<Surface_mesh>::Point> Custom_placement::operator()(const SMS::Edge_profile<Surface_mesh>& profile) const {
	typedef boost::optional<SMS::Edge_profile<Surface_mesh>::Point> result_type;

	std::list<std::pair <K::Vector_3, K::FT>> r1;
	if (params.volume_preservation > 0 || params.volume_optimisation > 0) r1 = volume_preservation_and_optimisation(profile);
	std::list<std::pair <K::Vector_3, K::Vector_3>> r2;
	if (params.boundary_preservation > 0 || params.boundary_optimization > 0) r2 = boundary_preservation_and_optimisation(profile);
	std::list<K::Vector_3> r3;
	if (params.triangle_shape_optimization > 0) r3 = triangle_shape_optimization(profile);
	std::list<std::pair <K::Vector_3, K::FT>> r4;
	if (params.label_preservation > 0) r4 = label_preservation(profile, point_cloud, ablation);
	std::list<K::Vector_3> r5;
	if (params.semantic_border_optimization > 0 && r4.size() < 3) r5 = semantic_border_optimization(profile, point_cloud);

	Eigen_vector B(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0 && r2.size() > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0) + ((params.semantic_border_optimization > 0) ? 3*r5.size() : 0));
	Eigen_matrix A(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0 && r2.size() > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0) + ((params.semantic_border_optimization > 0) ? 3*r5.size() : 0), 3);

	std::size_t i = 0;

	// Volume preservation
	if (params.volume_preservation > 0) {
		K::Vector_3 n (CGAL::NULL_VECTOR);
		K::FT det (0);
		for (const auto &vpo: r1) {
			n += vpo.first;
			det += vpo.second;
		}
		A.set(i, 0, n.x()/3*params.volume_preservation);
		A.set(i, 1, n.y()/3*params.volume_preservation);
		A.set(i, 2, n.z()/3*params.volume_preservation);
		B.set(i++, det/3*params.volume_preservation);
	}

	// Volume optimisation
	if (params.volume_optimisation > 0) {
		for (const auto &vpo: r1) {
			A.set(i, 0, vpo.first.x()/3*params.volume_optimisation);
			A.set(i, 1, vpo.first.y()/3*params.volume_optimisation);
			A.set(i, 2, vpo.first.z()/3*params.volume_optimisation);
			B.set(i++, vpo.second/3*params.volume_optimisation);
		}
	}

	// Boundary preservation
	if (params.boundary_preservation > 0 && r2.size() > 0) {
		K::Vector_3 e1 (CGAL::NULL_VECTOR);
		K::Vector_3 e2 (CGAL::NULL_VECTOR);
		for (const auto &vpo: r2) {
			e1 += vpo.first;
			e2 += vpo.second;
		}
		A.set(i, 0, 0);
		A.set(i, 1, -e1.z()/2*params.boundary_preservation);
		A.set(i, 2, e1.y()/2*params.boundary_preservation);
		B.set(i++, e2.x()/2*params.boundary_preservation);
		A.set(i, 0, e1.z()/2*params.boundary_preservation);
		A.set(i, 1, 0);
		A.set(i, 2, -e1.x()/2*params.boundary_preservation);
		B.set(i++, e2.y()/2*params.boundary_preservation);
		A.set(i, 0, -e1.y()/2*params.boundary_preservation);
		A.set(i, 1, e1.x()/2*params.boundary_preservation);
		A.set(i, 2, 0);
		B.set(i++, e2.z()/2*params.boundary_preservation);
	}

	// Boundary optimisation
	if (params.boundary_optimization > 0) {
		for (const auto &vpo: r2) {
			A.set(i, 0, 0);
			A.set(i, 1, -vpo.first.z()/2*params.boundary_optimization);
			A.set(i, 2, vpo.first.y()/2*params.boundary_optimization);
			B.set(i++, vpo.second.x()/2*params.boundary_optimization);
			A.set(i, 0, vpo.first.z()/2*params.boundary_optimization);
			A.set(i, 1, 0);
			A.set(i, 2, -vpo.first.x()/2*params.boundary_optimization);
			B.set(i++, vpo.second.y()/2*params.boundary_optimization);
			A.set(i, 0, -vpo.first.y()/2*params.boundary_optimization);
			A.set(i, 1, vpo.first.x()/2*params.boundary_optimization);
			A.set(i, 2, 0);
			B.set(i++, vpo.second.z()/2*params.boundary_optimization);
		}
	}

	// Triange shape optimisation
	if (params.triangle_shape_optimization > 0) {
		for (const auto &vpo: r3) {
			A.set(i, 0, params.triangle_shape_optimization);
			A.set(i, 1, 0);
			A.set(i, 2, 0);
			B.set(i++, vpo.x()*params.triangle_shape_optimization);
			A.set(i, 0, 0);
			A.set(i, 1, params.triangle_shape_optimization);
			A.set(i, 2, 0);
			B.set(i++, vpo.y()*params.triangle_shape_optimization);
			A.set(i, 0, 0);
			A.set(i, 1, 0);
			A.set(i, 2, params.triangle_shape_optimization);
			B.set(i++, vpo.z()*params.triangle_shape_optimization);
		}
	}

	// Label preservation
	if (params.label_preservation > 0) {
		for (const auto &vpo: r4) {
			A.set(i, 0, vpo.first.x()*params.label_preservation);
			A.set(i, 1, vpo.first.y()*params.label_preservation);
			A.set(i, 2, vpo.first.z()*params.label_preservation);
			B.set(i++, vpo.second*params.label_preservation);
		}
	}

	// Semantic border optimization
	if (params.semantic_border_optimization > 0 && r4.size() < 3) {
		for (const auto &vpo: r5) {
			A.set(i, 0, params.semantic_border_optimization);
			A.set(i, 1, 0);
			A.set(i, 2, 0);
			B.set(i++, vpo.x()*params.semantic_border_optimization);
			A.set(i, 0, 0);
			A.set(i, 1, params.semantic_border_optimization);
			A.set(i, 2, 0);
			B.set(i++, vpo.y()*params.semantic_border_optimization);
			A.set(i, 0, 0);
			A.set(i, 1, 0);
			A.set(i, 2, params.semantic_border_optimization);
			B.set(i++, vpo.z()*params.semantic_border_optimization);
		}
	}

	// Solve AX=B
	auto C = B;
	CGAL::Eigen_svd::solve(A, B);

	auto R = A*B - C;

// if (profile.v0().idx() == 6984 && profile.v1().idx() == 7620) { // save cost detail
// 	std::ofstream outfile;

// 	outfile.open("placement.output", std::ios::app);
// 	if (profile.v0().idx() < profile.v1().idx()) {
// 		outfile << "edge " << profile.v0().idx() << " -> " << profile.v1().idx() << ":\n";
// 	} else {
// 		outfile << "edge " << profile.v1().idx() << " -> " << profile.v0().idx() << ":\n";
// 	}
// 	outfile << "1 - " << r1.size() << " - 3 - " << (3*r2.size()) << " - " << (3*r3.size()) << " - " << r4.size() << " - " << (3*r5.size()) << "\n";
// 	outfile << "A:\n" << A << "\n";
// 	outfile << "B:\n" << C << "\n";
// 	outfile << "R:\n" << R << "\n";
// }

	// Save cost
	Point_3 placement(B.vector()[0], B.vector()[1], B.vector()[2]);
	collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].placement_cost = (R.transpose()*R)(0,0);

	return result_type(placement);
}

Custom_cost::Custom_cost (const LindstromTurk_param &params, const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT delta, const K::FT min_point_per_area, Surface_mesh &mesh, const Point_set &point_cloud, char *next_mesh) : params(params), alpha(alpha), beta(beta), gamma(gamma), delta(delta), min_point_per_area(min_point_per_area), point_cloud(point_cloud), next_mesh(next_mesh) {
	bool created_collapse_datas;
	boost::tie(collapse_datas, created_collapse_datas) = mesh.add_property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");
}

boost::optional<SMS::Edge_profile<Surface_mesh>::FT> Custom_cost::operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
	typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;
	CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;

	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
	bool has_face_costs;
	boost::tie(face_costs, has_face_costs) = profile.surface_mesh().property_map<Surface_mesh::Face_index, K::FT>("f:cost");
	assert(has_face_costs);

	collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.clear();

	if (placement) {

		K::FT squared_distance = 0;
		int count_semantic_error = 0;
		K::FT semantic_border_length = 0;
		K::FT old_cost = 0;
		if (alpha > 0 || beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {

			Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
			bool has_point_in_face;
			boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
			assert(has_point_in_face);

			std::set<Point_set::Index> points_to_be_change;
			for (const auto &face: profile.triangles()) {
				auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
				points_to_be_change.insert(point_in_face[fh].begin(), point_in_face[fh].end());
				old_cost += face_costs[fh];
// std::cerr << "\tface " << fh << " cost\t" << face_costs[fh] << "\n";
			}

			std::vector<K::Triangle_3> new_faces;
			std::vector<Surface_mesh::Halfedge_index> new_faces_border_halfedge;
			Point_3 C = *placement;
			for (std::size_t face_id = 0; face_id < profile.link().size(); face_id++) {
				auto he = profile.surface_mesh().halfedge(profile.link()[face_id], profile.link()[(face_id + 1) % profile.link().size()]);
				if (he != Surface_mesh::null_halfedge() && !profile.surface_mesh().is_border(he)) {
					Point_3 A = get(profile.vertex_point_map(), profile.link()[face_id]);
					Point_3 B = get(profile.vertex_point_map(), profile.link()[(face_id + 1) % profile.link().size()]);
					new_faces.push_back(K::Triangle_3(A, B, C));
					new_faces_border_halfedge.push_back(he);
				}
			}

			if (new_faces.size() == 0) return result_type();

			std::vector<K::FT> new_face_cost (new_faces.size(), 0);

			// geometric error
			std::vector<std::list<Point_set::Index>> points_in_new_face (new_faces.size());
//std::cerr << "points_in_new_face.size: " << points_in_new_face.size() << "\n";
			for (const auto &ph: points_to_be_change) {
				auto point = type_converter(point_cloud.point(ph));
				K::FT min_d = std::numeric_limits<K::FT>::max();
				std::size_t closest_face = 0;
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					auto d = CGAL::squared_distance(point, new_faces[face_id]);
					if (d < min_d) {
						min_d = d;
						closest_face = face_id;
					}
				}
				points_in_new_face[closest_face].push_back(ph);
				if (alpha > 0) {
					squared_distance += min_d;
					new_face_cost[closest_face] += alpha * min_d;
				}
			}

			// semantic error
			if (beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
				Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
				bool has_mesh_label;
				boost::tie(mesh_label, has_mesh_label) = profile.surface_mesh().property_map<Surface_mesh::Face_index, unsigned char>("f:label");
				assert(has_mesh_label);

				Point_set::Property_map<unsigned char> point_cloud_label;
				bool has_label;
				boost::tie(point_cloud_label, has_label) = point_cloud.property_map<unsigned char>("p:label");
				assert(has_label);

				std::vector<unsigned char> new_face_label(new_faces.size(), LABEL_OTHER);

				// label new face and semantic error
				std::set<std::size_t> faces_with_no_label;
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					K::FT min_point = 0;
					if (min_point_per_area > 0) {
						min_point = min_point_per_area * CGAL::sqrt(new_faces[face_id].squared_area());
					}
					if (points_in_new_face[face_id].size() > 0) {

						if (points_in_new_face[face_id].size() > min_point) {

							int face_label[LABELS.size()] = {0};

							for (const auto &point: points_in_new_face[face_id]) {
								face_label[point_cloud_label[point]]++;
							}

							auto argmax = std::max_element(face_label, face_label+LABELS.size());
							count_semantic_error += points_in_new_face[face_id].size() - *argmax;
							new_face_cost[face_id] += beta * (points_in_new_face[face_id].size() - *argmax);
// if (profile.v0_v1().idx() == 1710) std::cerr << "face_id: " << face_id << " (" << profile.surface_mesh().face(new_faces_border_halfedge[face_id]) << ") (beta1 x nb_error): " << beta * (points_in_new_face[face_id].size() - *argmax) << "\n";
							new_face_label[face_id] = argmax - face_label;

						} else { // We don't have information about this face.
							new_face_label[face_id] = LABEL_OTHER;
							count_semantic_error += points_in_new_face[face_id].size();
							new_face_cost[face_id] += beta * points_in_new_face[face_id].size();
// if (profile.v0_v1().idx() == 1710) std::cerr << "face_id: " << face_id << " (" << profile.surface_mesh().face(new_faces_border_halfedge[face_id]) << ") (beta2 x nb_error): " << beta * points_in_new_face[face_id].size() << "\n";
						}
					} else {
						if (min_point <= 1) { // It's probable that there is no point
							faces_with_no_label.insert(face_id);
						} else { // We don't have information about this face.
							new_face_label[face_id] = LABEL_OTHER;
						}
					}
				}
				while(faces_with_no_label.size() > 0) {
					std::list<std::size_t> face_to_be_removed;

					for (const auto &face_id: faces_with_no_label) {
						K::FT face_label[LABELS.size()] = {0};
						bool no_neighbor = true;
						// Outside face
						auto border_he = profile.surface_mesh().opposite(new_faces_border_halfedge[face_id]);
						if (!profile.surface_mesh().is_border(border_he)) {
							no_neighbor = false;
							face_label[mesh_label[profile.surface_mesh().face(border_he)]] += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(0), new_faces[face_id].vertex(1)));
						}
						if (profile.surface_mesh().target(new_faces_border_halfedge[(face_id - 1 + new_faces_border_halfedge.size()) % new_faces_border_halfedge.size()]) == profile.surface_mesh().source(new_faces_border_halfedge[face_id])) {
							no_neighbor = false;
							if(faces_with_no_label.count((face_id - 1 + new_faces.size()) % new_faces.size()) == 0) {
								face_label[new_face_label[(face_id - 1 + new_faces.size()) % new_faces.size()]] += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(0), C));
							}
						}
						if (profile.surface_mesh().source(new_faces_border_halfedge[(face_id + 1) % new_faces_border_halfedge.size()]) == profile.surface_mesh().target(new_faces_border_halfedge[face_id])) {
							no_neighbor = false;
							if(faces_with_no_label.count((face_id + 1) % new_faces.size()) == 0) {
								face_label[new_face_label[(face_id + 1) % new_faces.size()]] += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(1), C));
							}
						}
						if (no_neighbor) {
							new_face_label[face_id] = LABEL_OTHER;
							face_to_be_removed.push_back(face_id);
						} else {
							auto argmax = std::max_element(face_label, face_label+LABELS.size());
							if (*argmax > 0) {
								new_face_label[face_id] = argmax - face_label;
								face_to_be_removed.push_back(face_id);
							}
						}
					}
					for (const auto &face_id: face_to_be_removed) {
						faces_with_no_label.erase(face_id);
					}
					if (face_to_be_removed.size() == 0) {
						for (const auto &face_id: faces_with_no_label) {
							new_face_label[face_id] = LABEL_OTHER;
						}
						break;
					} else {
						face_to_be_removed.clear();
					}
				}

				// semantic border length error
				if (gamma > 0) {
					for (const auto &h: profile.surface_mesh().halfedges_around_target(profile.v1_v0())) {
						if (!profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
							if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))]) {
								semantic_border_length -= CGAL::sqrt(K::Vector_3(get(profile.vertex_point_map(), profile.surface_mesh().source(h)), get(profile.vertex_point_map(), profile.surface_mesh().target(h))).squared_length());
							}
						}
					}

					for (const auto &h: profile.surface_mesh().halfedges_around_target(profile.v0_v1())) {
						if (h != profile.v0_v1() && !profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
							if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))]) {
								semantic_border_length -= CGAL::sqrt(K::Vector_3(get(profile.vertex_point_map(), profile.surface_mesh().source(h)), get(profile.vertex_point_map(), profile.surface_mesh().target(h))).squared_length());
							}
						}
					}

					for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
						if (profile.surface_mesh().source(new_faces_border_halfedge[(face_id + 1) % new_faces_border_halfedge.size()]) == profile.surface_mesh().target(new_faces_border_halfedge[face_id]) && new_face_label[face_id] != new_face_label[(face_id + 1) % new_faces_border_halfedge.size()]) {
							semantic_border_length += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(1), C));
						}
					}

					for(std::size_t face_id = 0; face_id < new_faces_border_halfedge.size(); face_id++) {
						if (!profile.surface_mesh().is_border(profile.surface_mesh().opposite(new_faces_border_halfedge[face_id]))) {
							bool old_diff = (mesh_label[profile.surface_mesh().face(new_faces_border_halfedge[face_id])] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(new_faces_border_halfedge[face_id]))]);
							bool new_diff = (new_face_label[face_id] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(new_faces_border_halfedge[face_id]))]);

							if (old_diff != new_diff) {
								K::FT length = CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(0), new_faces[face_id].vertex(1)));
								if (new_diff) {
									semantic_border_length += length;
								} else {
									semantic_border_length -= length;
								}
							}
						}
					}
				}

				collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.reserve(new_faces.size());
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					CollapseDataElement r;
					r.halfedge = new_faces_border_halfedge[face_id];
					r.label = new_face_label[face_id];
					r.cost = new_face_cost[face_id];
					r.points.reserve(points_in_new_face[face_id].size());
					for (const auto &ph: points_in_new_face[face_id]) r.points.push_back(ph);
					collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.push_back(r);
				}

				if (next_mesh != nullptr) {
					// Output collapse surface
					Surface_mesh output_mesh;

					bool created;
					Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> red;
					Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> green;
					Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> blue;
					boost::tie(red, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("red",0);
					assert(created);
					boost::tie(green, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("green",0);
					assert(created);
					boost::tie(blue, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("blue",0);
					assert(created);

					for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
						auto &t = new_faces[face_id];
						auto v0 = output_mesh.add_vertex(t.vertex(0));
						auto v1 = output_mesh.add_vertex(t.vertex(1));
						auto v2 = output_mesh.add_vertex(t.vertex(2));
						auto f = output_mesh.add_face(v0, v1, v2);

						red[f] = LABELS.at(new_face_label[face_id]).red;
						green[f] = LABELS.at(new_face_label[face_id]).green;
						blue[f] = LABELS.at(new_face_label[face_id]).blue;
					}

					K::FT export_cost = (- old_cost + alpha * squared_distance + beta * count_semantic_error + gamma * semantic_border_length + delta * collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].placement_cost);

					std::stringstream output_mesh_name;
					if (profile.v0().idx() < profile.v1().idx()) {
						output_mesh_name << next_mesh << "/next_mesh_" << export_cost << "_" << profile.v0().idx() << "_" << profile.v1().idx() << "_" << profile.surface_mesh().point(profile.v0()) << "_" << profile.surface_mesh().point(profile.v1()) << ".ply";
					} else {
						output_mesh_name << next_mesh << "/next_mesh_" << export_cost << "_" << profile.v1().idx() << "_" << profile.v0().idx() << "_" << profile.surface_mesh().point(profile.v1()) << "_" << profile.surface_mesh().point(profile.v0()) << ".ply";
					}
					std::stringstream output_mesh_coment;
					output_mesh_coment << "cost: " << export_cost;
					output_mesh_coment << "\nold_cost: " << old_cost;
					output_mesh_coment << "\nalpha * squared_distance: " << (alpha * squared_distance);
					output_mesh_coment << "\nbeta * count_semantic_error: " << (beta * count_semantic_error);
					output_mesh_coment << "\ngamma * semantic_border_length: " << (gamma * semantic_border_length);
					output_mesh_coment << "\ndelta * collapse_datas[].placement_cost: " << (delta * collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].placement_cost);
					std::ofstream mesh_ofile (output_mesh_name.str().c_str());
					CGAL::IO::write_PLY (mesh_ofile, output_mesh, output_mesh_coment.str());
					mesh_ofile.close();
				}
			} else {
				collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.reserve(new_faces.size());
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					CollapseDataElement r;
					r.halfedge = new_faces_border_halfedge[face_id];
					r.cost = new_face_cost[face_id];
					r.points.reserve(points_in_new_face[face_id].size());
					for (const auto &ph: points_in_new_face[face_id]) r.points.push_back(ph);
					collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.push_back(r);
				}
			}
		}

		// std::cerr << "\tcost_explain detail " << profile.v0_v1() << "\t" << old_cost << "\t" << squared_distance << "\t" << count_semantic_error << "\n";
		// std::cerr << "\tcost_explain total " << profile.v0_v1() << "\t" << (- old_cost + alpha * squared_distance + beta * count_semantic_error) << "\t" << (gamma * semantic_border_length) << "\t" << (delta * collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].placement_cost) << "\n";

		K::FT final_cost = - old_cost + alpha * squared_distance + beta * count_semantic_error + gamma * semantic_border_length + delta * collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].placement_cost;

		collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].cost = final_cost;

		return result_type(final_cost);
	}

	return result_type();
}

Cost_stop_predicate::Cost_stop_predicate(const float cost) : cost(cost) {}

bool Cost_stop_predicate::operator()(const SMS::Edge_profile<Surface_mesh>::FT & current_cost, const SMS::Edge_profile<Surface_mesh> &, const SMS::Edge_profile<Surface_mesh>::edges_size_type, const SMS::Edge_profile<Surface_mesh>::edges_size_type) const {
	return current_cost > cost;
}

My_visitor::My_visitor(const LindstromTurk_param &params, const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT min_point_per_area, Surface_mesh &mesh, const Surface_mesh_info &mesh_info, Point_set &point_cloud, const Ablation_study ablation) : params(params), alpha(alpha), beta(beta), gamma(gamma), min_point_per_area(min_point_per_area), mesh(mesh), mesh_info(mesh_info), point_cloud(point_cloud), ablation(ablation) {
	bool created_collapse_datas;
	boost::tie(collapse_datas, created_collapse_datas) = mesh.add_property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");
}

void My_visitor::OnStarted (Surface_mesh&) {

	std::cout << "Starting edge_collapse" << std::endl;
	total_timer.start();

	if (beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
		// Add label to face
		bool created_label;
		boost::tie(mesh_label, created_label) = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("f:label", LABEL_OTHER);
		assert(created_label);
	}

	if (beta > 0 || gamma > 0 || params.label_preservation > 0 || params.semantic_border_optimization > 0) {
		// Get point_cloud_label
		bool has_point_cloud_label;
		boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
		assert(has_point_cloud_label);
	}

	// Create point_in_face
	bool created_point_in_face;
	boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points", std::list<Point_set::Index>());

	// Create face_costs
	bool created_face_costs;
	boost::tie(face_costs, created_face_costs) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:cost", 0);

	if(created_point_in_face) {
		AABB_tree mesh_tree;
		PMP::build_AABB_tree(mesh, mesh_tree);
		for (const auto &ph: point_cloud) {
			auto p = type_converter(point_cloud.point(ph));
			auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
			point_in_face[location.first].push_back(ph);
			if (created_face_costs && alpha > 0) face_costs[location.first] += alpha * CGAL::squared_distance(p, PMP::construct_point(location, mesh));
		}
	} else if (created_face_costs && alpha > 0) {
		for (const auto &face: mesh.faces()) {
			auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
			K::Triangle_3 triangle_face(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++));
			for (const auto &ph: point_in_face[face]) {
				auto p = type_converter(point_cloud.point(ph));
				face_costs[face] += alpha * CGAL::squared_distance(p, triangle_face);
			}
		}
	}
	if (created_face_costs && beta > 0) {
		for (const auto &face: mesh.faces()) {
			K::FT min_point = 0;
			if (min_point_per_area > 0) {
				auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
				K::FT face_area = CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
				min_point = min_point_per_area * face_area;
			}
			if (point_in_face[face].size() > 0) {
				if (point_in_face[face].size() > min_point) {
					int face_label[LABELS.size()] = {0};
					for (const auto &point: point_in_face[face]) {
						face_label[point_cloud_label[point]]++;
					}
					auto argmax = std::max_element(face_label, face_label+LABELS.size());
					face_costs[face] += beta * (point_in_face[face].size() - *argmax);
				} else { // We don't have information about this face.
					face_costs[face] += beta * point_in_face[face].size();
				}
			}
		}
	}

	// Set isBorder
	if(params.label_preservation > 0 && ablation.border_point) {
		bool created_point_isborder;
		Point_set::Property_map<bool> isborder;
		boost::tie (isborder, created_point_isborder) = point_cloud.add_property_map<bool>("p:isborder", false);
		if (created_point_isborder) {
			int N = 10;
			Point_tree point_tree(point_cloud.begin(), point_cloud.end(), Point_tree::Splitter(), TreeTraits(point_cloud.point_map()));
			Neighbor_search::Distance tr_dist(point_cloud.point_map());
			for (const auto &point: point_cloud) {
				Neighbor_search search(point_tree, point_cloud.point(point), N, 0, true, tr_dist, false);
				unsigned char l = point_cloud_label[search.begin()->first];
				for (auto it = search.begin() + 1; it != search.end() && !isborder[point]; ++it) {
					if (point_cloud_label[it->first] != l) isborder[point] = true;
				}
			}
			std::cout << "Border points computed" << std::endl;
		} else {
			std::cout << "Border points found" << std::endl;
		}
	}

	{// Save point cloud
		if (ablation.border_point) {
			total_timer.pause();
			Point_set output_point_cloud(point_cloud);

			// Color
			bool created;
			Point_set::Property_map<unsigned char> red;
			Point_set::Property_map<unsigned char> green;
			Point_set::Property_map<unsigned char> blue;
			boost::tie(red, created) = output_point_cloud.add_property_map<unsigned char>("red",0);
			assert(created);
			boost::tie(green, created) = output_point_cloud.add_property_map<unsigned char>("green",0);
			assert(created);
			boost::tie(blue, created) = output_point_cloud.add_property_map<unsigned char>("blue",0);
			assert(created);

			if (params.label_preservation > 0) {
				Point_set::Property_map<bool> output_isborder;
				bool has_output_isborder;
				boost::tie(output_isborder, has_output_isborder) = output_point_cloud.property_map<bool>("p:isborder");
				assert(has_output_isborder);

				for (const auto &point: output_point_cloud) {
					if (output_isborder[point]) {
						red[point] = 255;
						green[point] = 255;
						blue[point] = 255;
					} else {
						red[point] = 155;
						green[point] = 155;
						blue[point] = 155;
					}
				}

				CGAL::IO::write_point_set("pc_with_border.ply", output_point_cloud);
			}

			if (beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
				Point_set::Property_map<unsigned char> output_label;
				bool has_output_label;
				boost::tie(output_label, has_output_label) = output_point_cloud.property_map<unsigned char>("p:label");
				assert(has_output_label);

				for (const auto &point: output_point_cloud) {
					red[point] = LABELS.at(output_label[point]).red;
					green[point] = LABELS.at(output_label[point]).green;
					blue[point] = LABELS.at(output_label[point]).blue;
				}
			}

			CGAL::IO::write_point_set("pc_with_color.ply", output_point_cloud);
			total_timer.resume();
		}
	}


	if (beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
		add_label(mesh, point_cloud, min_point_per_area);

		// Subdivide wrong face
		if (ablation.subdivide) {
			std::set<Surface_mesh::Face_index> face_to_divide;
			for (const auto &face: mesh.faces()) {
				face_to_divide.insert(face);
			}

			while (face_to_divide.size() > 0) {
				auto face = *face_to_divide.begin();

				if (point_in_face[face].size() > 1) {

					int face_label[LABELS.size()] = {0};

					for (const auto &point: point_in_face[face]) {
						face_label[point_cloud_label[point]]++;
					}

					auto argmax = std::max_element(face_label, face_label+LABELS.size());

					if (*argmax < point_in_face[face].size()) {
						auto h0 = mesh.halfedge(face);
						auto h1 = mesh.next(h0);
						auto h2 = mesh.next(h1);

						std::set<Point_set::Index> points_to_be_change;
						points_to_be_change.insert(point_in_face[face].begin(), point_in_face[face].end());
						if (!mesh.is_border(mesh.edge(h0))) points_to_be_change.insert(point_in_face[mesh.face(mesh.opposite(h0))].begin(), point_in_face[mesh.face(mesh.opposite(h0))].end());
						if (!mesh.is_border(mesh.edge(h1))) points_to_be_change.insert(point_in_face[mesh.face(mesh.opposite(h1))].begin(), point_in_face[mesh.face(mesh.opposite(h1))].end());
						if (!mesh.is_border(mesh.edge(h2))) points_to_be_change.insert(point_in_face[mesh.face(mesh.opposite(h2))].begin(), point_in_face[mesh.face(mesh.opposite(h2))].end());

						face_to_divide.erase(face);
						face_to_divide.erase(mesh.face(mesh.opposite(h0)));
						face_to_divide.erase(mesh.face(mesh.opposite(h1)));
						face_to_divide.erase(mesh.face(mesh.opposite(h2)));

						auto new_faces = subdivide_face(mesh, face, argmax - face_label, point_cloud);

						std::vector<K::Triangle_3> new_faces_triangle;
						for (const auto &new_face: new_faces) {
							auto r = mesh.vertices_around_face(mesh.halfedge(new_face)).begin();
							new_faces_triangle.push_back(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)));
						}

						// geometric error
						for (const auto &ph: points_to_be_change) {
							auto point = type_converter(point_cloud.point(ph));
							K::FT min_d = std::numeric_limits<K::FT>::max();
							std::size_t closest_face = 0;
							for(std::size_t face_id = 0; face_id < new_faces_triangle.size(); face_id++) {
								auto d = CGAL::squared_distance(point, new_faces_triangle[face_id]);
								if (d < min_d) {
									min_d = d;
									closest_face = face_id;
								}
							}
							point_in_face[new_faces[closest_face]].push_back(ph);
							if (alpha > 0) {
								face_costs[new_faces[closest_face]] += alpha * min_d;
							}
						}

						for (const auto &new_face: new_faces) {
							if (point_in_face[new_face].size() > 0) face_to_divide.insert(new_face);
						}

						// semantic error
						if (beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
							// label new face and semantic error
							std::set<Surface_mesh::Face_index> faces_with_no_label;
							for (const auto &face: new_faces) {
								K::FT min_point = 0;
								if (min_point_per_area > 0) {
									auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
									K::FT face_area = CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
									min_point = min_point_per_area * face_area;
								}
								if (point_in_face[face].size() > 0) {

									if (point_in_face[face].size() > min_point) {

										int face_label[LABELS.size()] = {0};

										for (const auto &point: point_in_face[face]) {
											face_label[point_cloud_label[point]]++;
										}

										auto argmax = std::max_element(face_label, face_label+LABELS.size());
										face_costs[face] += beta * (point_in_face[face].size() - *argmax);
										mesh_label[face] = argmax - face_label;

									} else { // We don't have information about this face.
										mesh_label[face] = LABEL_OTHER;
										face_costs[face] += beta * point_in_face[face].size();
									}
								} else {
									if (min_point <= 1) { // It's probable that there is no point
										faces_with_no_label.insert(face);
									} else { // We don't have information about this face.
										mesh_label[face] = LABEL_OTHER;
									}
								}
							}
							while(faces_with_no_label.size() > 0) {
								std::list<Surface_mesh::Face_index> face_to_be_removed;
								for (const auto &face: faces_with_no_label) {
									K::FT face_label[LABELS.size()] = {0};
									bool no_neighbor = true;
									for (const auto &he: mesh.halfedges_around_face(mesh.halfedge(face))) {
										if (!mesh.is_border(Surface_mesh::Edge_index(he))) {
											no_neighbor = false;
											if (faces_with_no_label.count(mesh.face(mesh.opposite(he))) == 0) {
												face_label[mesh_label[mesh.face(mesh.opposite(he))]] += CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(he)), mesh.point(mesh.target(he))));
											}
										}
									}
									if (no_neighbor) {
										mesh_label[face] = LABEL_OTHER;
										face_to_be_removed.push_back(face);
									} else {
										auto argmax = std::max_element(face_label, face_label+LABELS.size());
										if (*argmax > 0) {
											mesh_label[face] = argmax - face_label;
											face_to_be_removed.push_back(face);
										}
									}
								}
								for (const auto &face: face_to_be_removed) {
									faces_with_no_label.erase(face);
								}
								if (face_to_be_removed.size() == 0) {
									for (const auto &face_id: faces_with_no_label) {
										mesh_label[face_id] = LABEL_OTHER;
									}
									break;
								} else {
									face_to_be_removed.clear();
								}
							}
						}
					} else {
						face_to_divide.erase(face);
					}
				} else {
					face_to_divide.erase(face);
				}
			}
		}
	}

	total_timer.pause();
	mesh_info.save_mesh(mesh, "initial-mesh.ply");

	total_timer.resume();

/*{
	K::FT face_cost = 0;
	for (const auto &face: mesh.faces()) {
		face_cost += face_costs[face];
		// std::cerr << "\tface " << face << " cost\t" << face_costs[face] << "\n";
	}
	K::FT edge_cost = 0;
	for (const auto &e: mesh.edges()) {
		if (!mesh.is_border(e)) {
			auto h = mesh.halfedge(e);
			auto f1 = mesh.face(h);
			auto f2 = mesh.face(mesh.opposite(h));
			if (f1 != mesh.null_face() && f2 != mesh.null_face() && mesh_label[f1] != mesh_label[f2]) {
				edge_cost += 0.01 * CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(h)), mesh.point(mesh.target(h))));
			}
		}
	}
	std::cerr << "Initial cost\t" << (face_cost + edge_cost) << "\n";
	std::cerr << "Initial cost_explain\t\t" << face_cost << "\t" << edge_cost << "\n";
}*/

	collected_timer.start();
}

void My_visitor::OnFinished (Surface_mesh &mesh) {
	std::cout << "\rMesh simplified                                               " << std::endl;

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = mesh.property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
	bool has_face_costs;
	boost::tie(face_costs, has_face_costs) = mesh.property_map<Surface_mesh::Face_index, K::FT>("f:cost");
	assert(has_face_costs);

	Surface_mesh::Property_map<Surface_mesh::Edge_index, CollapseData> collapse_datas;
	bool has_collapse_datas;
	boost::tie(collapse_datas, has_collapse_datas) = mesh.property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");
	assert(has_collapse_datas);

	// { //output edge cost

	// 	Point_set output_point_cloud;

	// 	bool created;
	// 	Point_set::Property_map<unsigned char> red;
	// 	Point_set::Property_map<unsigned char> green;
	// 	Point_set::Property_map<unsigned char> blue;
	// 	Point_set::Property_map<float> quality;
	// 	boost::tie(red, created) = output_point_cloud.add_property_map<unsigned char>("red",0);
	// 	assert(created);
	// 	boost::tie(green, created) = output_point_cloud.add_property_map<unsigned char>("green",0);
	// 	assert(created);
	// 	boost::tie(blue, created) = output_point_cloud.add_property_map<unsigned char>("blue",0);
	// 	assert(created);
	// 	boost::tie(quality, created) = output_point_cloud.add_property_map<float>("quality",0);
	// 	assert(created);

	// 	CGAL::Cartesian_converter<K, Point_set_kernel> type_converter;

	// 	for (const auto &edge: mesh.edges()) {
	// 		K:FT cost = collapse_datas[edge].cost;

	// 		auto point = output_point_cloud.insert(type_converter(CGAL::midpoint(mesh.point(mesh.source(mesh.halfedge(edge))), mesh.point(mesh.target(mesh.halfedge(edge))))));

	// 		if (cost < 0) {
	// 			red[*point] = 255;
	// 			green[*point] = 200 - int(std::min(- cost / 10 * 200, (float) 200));
	// 			blue[*point] = 200 - int(std::min(- cost / 10 * 200, (float) 200));
	// 		} else {
	// 			red[*point] = 255 - int(std::min(cost / 50 * 255, (float) 255));
	// 			green[*point] = 255 - int(std::min(cost / 50 * 255, (float) 255));
	// 			blue[*point] = 255;
	// 		}

	// 		quality[*point] = cost;
	// 	}

	// 	std::ofstream mesh_ofile ("halfedge_collapsing_cost.ply");
	// 	CGAL::IO::write_PLY (mesh_ofile, output_point_cloud);
	// 	mesh_ofile.close();

	// }

	/*std::cerr << "Total cost\t" << total_cost << "\n";

	{
		K::FT cost = 0;
		for (const auto &face: mesh.faces()) {
			cost += face_costs[face];

			K::FT sum_squared_distance = 0;
			K::FT sum_false_point = 0;

			auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
			K::Triangle_3 triangle_face(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++));
			for (const auto &point: point_in_face[face]) {
				sum_squared_distance += CGAL::squared_distance(triangle_face, point_cloud.point(point);
				if (point_cloud_label[point] != mesh_label[face]) sum_false_point += 1;
			}

			if (abs((alpha * sum_squared_distance + beta * sum_false_point) - face_costs[face]) != 0) {

				std::cerr << "face_costs[" << face << "] = " << face_costs[face] << "\n";
				std::cerr << "alpha * sum_squared_distance + beta * sum_false_point = " << (alpha * sum_squared_distance + beta * sum_false_point) << "\n";
			}

		}
		for (const auto &e: mesh.edges()) {
			if (!mesh.is_border(e)) {
				auto h = mesh.halfedge(e);
				auto f1 = mesh.face(h);
				auto f2 = mesh.face(mesh.opposite(h));
				if (f1 != mesh.null_face() && f2 != mesh.null_face() && mesh_label[f1] != mesh_label[f2]) {
					cost += 0.01 * CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(h)), mesh.point(mesh.target(h))));
				}
			}
		}
		std::cerr << "Final cost\t" << cost << "\n";
	}

	{
		int num_point = 0;
		for (const auto &face: mesh.faces()) {
			num_point += point_in_face[face].size();
		}
		std::cerr << "Final num_point_in_face\t" << num_point << "\n";
		std::cerr << "Point cloud num_point\t" << point_cloud.size() << "\n";
	}*/

	mesh.remove_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>(point_in_face);
	mesh.remove_property_map<Surface_mesh::Face_index, K::FT>(face_costs);
	mesh.remove_property_map<Surface_mesh::Edge_index, CollapseData>(collapse_datas);
}

void My_visitor::OnCollected(const SMS::Edge_profile<Surface_mesh>&, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >&) {
	collapsing_timer.start();
	i_collecte++;
	if (i_collecte%1000 == 0) {
		auto time = collected_timer.getElapsedTime();
		std::cout << "\rCollecte: " << i_collecte << "/" << mesh.number_of_edges() << " (" << ((int) (((float) i_collecte)/mesh.number_of_edges()*100)) << "%)" << " still " << (((float) mesh.number_of_edges() - i_collecte) * time / i_collecte) << "s" << " (" << (((float) i_collecte) / time) << " op/s)" << std::flush;
	}
}

void My_visitor::OnSelected (const SMS::Edge_profile<Surface_mesh>&, boost::optional< SMS::Edge_profile<Surface_mesh>::FT > cost, const SMS::Edge_profile<Surface_mesh>::edges_size_type initial_edge_count, const SMS::Edge_profile<Surface_mesh>::edges_size_type current_edge_count) {
	if (current_edge_count%100 == 0) {
		auto time = collapsing_timer.getElapsedTime();
		std::cout << "\rCollapse: " << (initial_edge_count-current_edge_count) << "/" << initial_edge_count << " (" << ((int) (((float) (initial_edge_count-current_edge_count))/initial_edge_count*100)) << "%)" << " still " << (((float) current_edge_count) * time / (initial_edge_count-current_edge_count)) << "s" << " (" << (((float) (initial_edge_count-current_edge_count)) / time) << " op/s)";
		if (cost) {
			std::cout << " - cost: " << *cost << "     " << std::flush;
		}
	}

	if (cost && ablation.step_mesh) {
// std::cerr << "Try collapsing " << profile.v0_v1() << "\t\t" << *cost << "\n";
// c_cost = *cost;

		if (current_edge_count < last_size) {
			total_timer.pause();
			collapsing_timer.pause();
			compute_stat(mesh, ablation, total_timer, *cost);
			last_size = ((int) ((current_edge_count) / 1000)) * 1000;
			total_timer.resume();
			collapsing_timer.resume();
		}

		int cost_id = -1;
		if(!output[++cost_id] && *cost > 100) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c100.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 10) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c10.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 1) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c1.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 0.3) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.03.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 0.2) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.02.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 0.1) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.01.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 0.001) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.001.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && *cost > 0) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		}
		cost_id = 8;
		if(!output[++cost_id] && current_edge_count <= 1000000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-1000000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 250000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-250000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 100000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-100000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 50000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-50000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 40000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-40000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 39000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-39000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 10000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-10000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 5000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-5000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 2500) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-2500.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 2000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-2000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 1000) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-1000.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 500) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-500.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 400) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-400.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 300) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-300.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 100) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-100.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 50) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-50.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		} else if(!output[++cost_id] && current_edge_count <= 10) {
			total_timer.pause();
			collapsing_timer.pause();
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-10.ply");
			// compute_stat(mesh, ablation, total_timer, *cost);
			total_timer.resume();
			collapsing_timer.resume();
		}
		// std::cerr << "Used_cost: " << *cost << "\n";
	}

}

void My_visitor::OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&) {
	// Called when an edge is about to be collapsed and replaced by a vertex whose position is *placement
	if (alpha > 0 || beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
		for (const auto &face: profile.triangles()) {
			auto fh = mesh.face(mesh.halfedge(face.v0, face.v1));
			/*std::cerr << "fh: " << fh << "\n";
			Point_3 p0 = get(profile.vertex_point_map(),face.v0);
			Point_3 p1 = get(profile.vertex_point_map(),face.v1);
			Point_3 p2 = get(profile.vertex_point_map(),face.v2);
			std::cerr << "p0: " << p0 << " " << p1 << " " << p2 << "\n";
			for (const auto &v: point_in_face[fh]) {
				std::cerr << "v: " << v << " : " << point_cloud.point(v) << "\n";
			}*/

			point_in_face[fh].clear();
			face_costs[fh] = 0;
		}
	}
}

void My_visitor::OnCollapsed (const SMS::Edge_profile<Surface_mesh>& prof, const Surface_mesh::Vertex_index vd) {
	// Called when an edge has been collapsed and replaced by the vertex vd

	// change point_in_face and face_costs
	if (alpha > 0 || beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) {
		for (const auto &element: collapse_datas[Surface_mesh::Edge_index(prof.v0_v1())].elements) {
			auto face = mesh.face(element.halfedge);
			for (const auto &ph: element.points) point_in_face[face].push_back(ph);
			face_costs[face] = element.cost;
			if (beta > 0 || gamma > 0 || params.semantic_border_optimization > 0) mesh_label[face] = element.label;
		}
	}

/*{
	std::cerr << "Collapsed cost\t" << c_cost << "\n";

	K::FT face_cost = 0;
	for (const auto &face: mesh.faces()) {
		face_cost += face_costs[face];
		// std::cerr << "\tface " << face << " cost\t" << face_costs[face] << "\n";
	}
	K::FT edge_cost = 0;
	for (const auto &e: mesh.edges()) {
		if (!mesh.is_border(e)) {
			auto h = mesh.halfedge(e);
			auto f1 = mesh.face(h);
			auto f2 = mesh.face(mesh.opposite(h));
			if (f1 != mesh.null_face() && f2 != mesh.null_face() && mesh_label[f1] != mesh_label[f2]) {
				edge_cost += 0.01 * CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(h)), mesh.point(mesh.target(h))));
			}
		}
	}
	std::cerr << "Intermediate cost\t\t" << (face_cost + edge_cost) << "\n";
	std::cerr << "Intermediate cost_explain\t\t" << face_cost << "\t" << edge_cost << "\n";

}*/

// total_cost += c_cost;
}

Point_set compute_point_cloud (Surface_mesh& mesh) {
	std::list<K::Point_3> out;
	PMP::sample_triangle_mesh(mesh, std::back_inserter(out), CGAL::parameters::use_grid_sampling(true).do_sample_edges(false).do_sample_vertices(false));

	CGAL::Cartesian_converter<K, Point_set_kernel> type_converter;

	Point_set point_cloud;
	for (const auto& p: out) {
		point_cloud.insert(type_converter(p));
	}

	return point_cloud;
}

void associate_mesh_point_cloud (Surface_mesh& mesh, Point_set& point_cloud) {
	Point_set::Property_map<unsigned char> point_cloud_label;
	bool created_point_label;
	boost::tie (point_cloud_label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", LABEL_OTHER);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	bool has_mesh_label;
	boost::tie(mesh_label, has_mesh_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("f:label");
	assert(has_mesh_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool created_point_in_face;
	boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points", std::list<Point_set::Index>());
	assert(created_point_in_face);

	AABB_tree mesh_tree;
	PMP::build_AABB_tree(mesh, mesh_tree);

	CGAL::Cartesian_converter<Point_set_kernel, K> type_converter;

	for (auto &ph: point_cloud) {
		auto p = type_converter(point_cloud.point(ph));
		auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
		if (created_point_label) point_cloud_label[ph] = mesh_label[location.first];
		point_in_face[location.first].push_back(ph);
	}
}
