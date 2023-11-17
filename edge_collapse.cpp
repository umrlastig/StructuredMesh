#include "edge_collapse.hpp"

#include <list>
#include <limits>

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

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

typedef CGAL::Search_traits_3<Exact_predicates_kernel>      Traits_base;
typedef CGAL::Search_traits_adapter<Point_set::Index, Point_set::Point_map, Traits_base> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>      Neighbor_search;
typedef Neighbor_search::Tree                               Point_tree;
typedef CGAL::AABB_tree<CGAL::AABB_traits<K, CGAL::AABB_triangle_primitive<K, std::vector<K::Triangle_3>::iterator>>> Triangle_tree;

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>        AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                AABB_tree;

K::FT get_mean_point_per_area(Surface_mesh &mesh, const Point_set &point_cloud) {
	K::FT total_area = 0;
	for(auto face: mesh.faces()) {
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
	boost::tie(mesh_label, has_mesh_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
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
	for(auto face: mesh.faces()) {
		K::FT min_point = 0;
		if (min_point_per_area > 0) {
			auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
			K::FT face_area = CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
			min_point = min_point_per_area * face_area;
		}	
		if (point_in_face[face].size() > 0) {

			if (point_in_face[face].size() > min_point) {

				int face_label[LABELS.size()] = {0};

				for (auto point: point_in_face[face]) {
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
		for(auto face: faces_with_no_label) {
			K::FT face_label[LABELS.size()] = {0};
			bool no_neighbor = true;
			for (auto he: mesh.halfedges_around_face(mesh.halfedge(face))) {
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
		for (auto face: face_to_be_removed) {
			faces_with_no_label.erase(face);
		}
		face_to_be_removed.clear();
	}
}

std::list<std::pair <K::Vector_3, K::FT>> volume_preservation_and_optimisation (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<std::pair <K::Vector_3, K::FT>> result;

	SMS::Edge_profile<Surface_mesh>::Triangle_vector triangles = profile.triangles();
	for (auto triangle: triangles) {
		Point_3 p0 = get(profile.vertex_point_map(),triangle.v0);
		Point_3 p1 = get(profile.vertex_point_map(),triangle.v1);
		Point_3 p2 = get(profile.vertex_point_map(),triangle.v2);

		if (!CGAL::collinear (p0, p1, p2)) {
			/*std::cerr << "triangle: " << p0 << ", " << p1 << ", " << p2 << "\n";
			std::cerr << "normal: " << CGAL::normal(p0, p1, p2) << "\n";
			std::cerr << "determinant: " << CGAL::determinant(K::Vector_3(Point_3(CGAL::ORIGIN), p0), K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p2)) << "\n";*/
			
			K::Vector_3 n = CGAL::normal(p0, p1, p2) / 2;
			K::FT det = CGAL::scalar_product (n, K::Vector_3(Point_3(CGAL::ORIGIN), p0));

			result.push_back(std::make_pair(n, det));
		}
	}

	return result;
}

std::list<std::pair <K::Vector_3, K::Vector_3>> boundary_preservation_and_optimisation (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<std::pair <K::Vector_3, K::Vector_3>> result;

	for (auto edge: profile.border_edges()) {
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

	for (auto v: profile.link()) {
		result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), v)));
	}

	return result;

}

std::list<std::pair <K::Vector_3, K::FT>> label_preservation (const SMS::Edge_profile<Surface_mesh>& profile, const Point_set &point_cloud) {

	std::list<std::pair <K::Vector_3, K::FT>> result;

	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	Point_set::Property_map<bool> isborder;
	bool has_isborder;
	boost::tie(isborder, has_isborder) = point_cloud.property_map<bool>("p:isborder");
	assert(has_isborder);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	std::set<Point_set::Index> points_in_faces;
	int count_collapse_label[LABELS.size()] = {0};
	for(auto face: profile.triangles()) {
		auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
		if (point_in_face[fh].size() > 0) {
			int count_face_label[LABELS.size()] = {0};
			for (auto ph: point_in_face[fh]) {
				if (isborder[ph]) {
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
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

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
		for (auto point: points_in_faces) {
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
			//add no border point to points_in_faces and restart SVM
			points_in_faces.clear();
			for(auto face: profile.triangles()) {
				auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
				if (point_in_face[fh].size() > 0) {
					points_in_faces.insert(point_in_face[fh].begin(), point_in_face[fh].end());
				}
			}
			y.clear();
			points_for_svm.clear();
			y.reserve(points_in_faces.size());
			points_for_svm.reserve(points_in_faces.size());
			for (auto point: points_in_faces) {
				if (label[point] == label1) {
					points_for_svm.push_back(point);
					y.push_back(1);
				} else if (label[point] == label2) {
					points_for_svm.push_back(point);
					y.push_back(-1);
				}
			}
		}

		Exact_predicates_kernel::FT c = 1;

		Program qp (CGAL::EQUAL, true, 0, true, c);
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			auto vr = Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) * y[i];
			qp.set_d(i, i, CGAL::scalar_product(vr,vr));
			for (std::size_t j = 0; j < i; j++) {
				qp.set_d(i, j, CGAL::scalar_product(vr, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[j])) * y[j]));
			}
			qp.set_c(i, -1);
			qp.set_a(i, 0,  y[i]);
		}
		Solution s = CGAL::solve_quadratic_program(qp, ET());
		//std::cerr << s << "\n";
		if (!s.solves_quadratic_program(qp)) std::cerr << "ALERT !!!!!!!!!!!!!!!!!!!\n";
		assert (s.solves_quadratic_program(qp));

		Exact_predicates_kernel::Vector_3 w(CGAL::NULL_VECTOR);
		auto value = s.variable_values_begin();
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			w += y[i] * CGAL::to_double(*(value++)) * Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]));
		}

		Exact_predicates_kernel::FT b = 0;
		int count = 0;
		Exact_predicates_kernel::FT min_positive = std::numeric_limits<Exact_predicates_kernel::FT>::max();
		Exact_predicates_kernel::FT max_negative = std::numeric_limits<Exact_predicates_kernel::FT>::lowest();
		for (std::size_t i = 0; i < points_for_svm.size(); i++) {
			Exact_predicates_kernel::FT v = CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
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
					b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
					count++;
				}
			}
		}
		if (count == 0) {
			value = s.variable_values_begin();
			for (std::size_t i = 0; i < points_for_svm.size(); i++) {
				if (CGAL::to_double(*(value++)) > 0) {
					b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
					count++;
				}
			}
		}
		if (count == 0) {
			for (std::size_t i = 0; i < points_for_svm.size(); i++) {
				b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
				count++;
			}
		}
		assert(count > 0);
		b /= count;

		result.push_back(std::pair<K::Vector_3, K::FT>(type_converter(w)*squared_length, -type_converter(b)*squared_length));

	} else {

		for (std::size_t i_label = 0; i_label < LABELS.size(); i_label++) {
			if (count_collapse_label[i_label] > 0) {
				
				std::vector<int> y;
				std::vector<Point_set::Index> points_for_svm;
				y.reserve(points_in_faces.size());
				points_for_svm.reserve(points_in_faces.size());
				bool has_i_label = false, has_other_label = false;
				for (auto point: points_in_faces) {
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

					Exact_predicates_kernel::FT c = 1;

					Program qp (CGAL::EQUAL, true, 0, true, c);
					for (std::size_t i = 0; i < points_for_svm.size(); i++) {
						auto vr = Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) * y[i];
						qp.set_d(i, i, CGAL::scalar_product(vr,vr));
						for (std::size_t j = 0; j < i; j++) {
							qp.set_d(i, j, CGAL::scalar_product(vr, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[j])) * y[j]));
						}
						qp.set_c(i, -1);
						qp.set_a(i, 0,  y[i]);
					}
					Solution s = CGAL::solve_quadratic_program(qp, ET());
					//std::cerr << s << "\n";
					if (!s.solves_quadratic_program(qp)) std::cerr << "ALERT !!!!!!!!!!!!!!!!!!!\n";
					assert (s.solves_quadratic_program(qp));

					Exact_predicates_kernel::Vector_3 w(CGAL::NULL_VECTOR);
					auto value = s.variable_values_begin();
					for (std::size_t i = 0; i < points_for_svm.size(); i++) {
						w += y[i] * CGAL::to_double(*(value++)) * Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]));
					}

					Exact_predicates_kernel::FT b = 0;
					int count = 0;
					Exact_predicates_kernel::FT min_positive = std::numeric_limits<Exact_predicates_kernel::FT>::max();
					Exact_predicates_kernel::FT max_negative = std::numeric_limits<Exact_predicates_kernel::FT>::lowest();
					for (std::size_t i = 0; i < points_for_svm.size(); i++) {
						Exact_predicates_kernel::FT v = CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
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
								b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
								count++;
							}
						}
					}
					if (count == 0) {
						value = s.variable_values_begin();
						for (std::size_t i = 0; i < points_for_svm.size(); i++) {
							if (CGAL::to_double(*(value++)) > 0) {
								b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
								count++;
							}
						}
					}
					if (count == 0) {
						for (std::size_t i = 0; i < points_for_svm.size(); i++) {
							b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
							count++;
						}
					}
					assert(count > 0);
					b /= count;

					result.push_back(std::pair<K::Vector_3, K::FT>(type_converter(w)*squared_length, -type_converter(b)*squared_length));
				}
			}
		}
	}

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
	boost::tie(mesh_label, has_mesh_label) = profile.surface_mesh().property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_mesh_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	for (auto h: profile.surface_mesh().halfedges_around_target(profile.v1_v0())) {
		if (h != profile.v1_v0() && h != profile.vL_v0() && !profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
			if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))] && point_in_face[profile.surface_mesh().face(h)].size() > 0 && point_in_face[profile.surface_mesh().face(profile.surface_mesh().opposite(h))].size() > 0) {
				result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), profile.surface_mesh().source(h))));
			}
		}
	}

	for (auto h: profile.surface_mesh().halfedges_around_target(profile.v0_v1())) {
		if (h != profile.v0_v1() && h != profile.vR_v1() && !profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
			if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))] && point_in_face[profile.surface_mesh().face(h)].size() > 0 && point_in_face[profile.surface_mesh().face(profile.surface_mesh().opposite(h))].size() > 0) {
				result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), profile.surface_mesh().source(h))));
			}
		}
	}

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

Custom_placement::Custom_placement (const LindstromTurk_param &params, Surface_mesh &mesh, const Point_set &point_cloud) : params(params), point_cloud(point_cloud) {
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
	if (params.label_preservation > 0) r4 = label_preservation(profile, point_cloud);
	std::list<K::Vector_3> r5;
	if (params.semantic_border_optimization > 0) r5 = semantic_border_optimization(profile, point_cloud);

	Eigen_vector B(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0 && r2.size() > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0) + ((params.semantic_border_optimization > 0) ? 3*r5.size() : 0));
	Eigen_matrix A(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0 && r2.size() > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0) + ((params.semantic_border_optimization > 0) ? 3*r5.size() : 0), 3);

	std::size_t i = 0;

	// Volume preservation
	if (params.volume_preservation > 0) {
		K::Vector_3 n (CGAL::NULL_VECTOR);
		K::FT det (0);
		for (auto vpo: r1) {
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
		for (auto vpo: r1) {
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
		for (auto vpo: r2) {
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
		for (auto vpo: r2) {
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
		for (auto vpo: r3) {
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
		for (auto vpo: r4) {
			A.set(i, 0, vpo.first.x()*params.label_preservation);
			A.set(i, 1, vpo.first.y()*params.label_preservation);
			A.set(i, 2, vpo.first.z()*params.label_preservation);
			B.set(i++, vpo.second*params.label_preservation);
		}
	}

	// Semantic border optimization
	if (params.semantic_border_optimization > 0) {
		for (auto vpo: r5) {
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
	
	// Save cost
	Point_3 placement(B.vector()[0], B.vector()[1], B.vector()[2]);
	collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].cost = (R.transpose()*R)(0,0);

	return result_type(placement);
}

Custom_cost::Custom_cost (const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT delta, const K::FT min_point_per_area, Surface_mesh &mesh, const Point_set &point_cloud) : alpha(alpha), beta(beta), gamma(gamma), delta(delta), min_point_per_area(min_point_per_area), point_cloud(point_cloud) {
	bool created_collapse_datas;
	boost::tie(collapse_datas, created_collapse_datas) = mesh.add_property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");
}

boost::optional<SMS::Edge_profile<Surface_mesh>::FT> Custom_cost::operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
	typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

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
		if (alpha > 0 || beta > 0 || gamma > 0) {

			Point_set::Property_map<unsigned char> point_cloud_label;
			bool has_label;
			boost::tie(point_cloud_label, has_label) = point_cloud.property_map<unsigned char>("p:label");
			assert(has_label);

			Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
			bool has_point_in_face;
			boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points");
			assert(has_point_in_face);

			std::set<Point_set::Index> points_to_be_change;
			for(auto face: profile.triangles()) {
				auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
				points_to_be_change.insert(point_in_face[fh].begin(), point_in_face[fh].end());
				old_cost += face_costs[fh];
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

			std::vector<K::FT> new_face_cost (new_faces.size(), 0);

			// geometric error
			std::vector<std::list<Point_set::Index>> points_in_new_face (new_faces.size());
			Triangle_tree tree(new_faces.begin(), new_faces.end());
			for(auto ph: points_to_be_change) {
				auto point = type_converter(point_cloud.point(ph));
				auto point_and_id = tree.closest_point_and_primitive(point);

				points_in_new_face[point_and_id.second - new_faces.begin()].push_back(ph);
				if (alpha > 0) {
					K::FT distance = CGAL::squared_distance(point, point_and_id.first);
					squared_distance += distance;
					new_face_cost[point_and_id.second - new_faces.begin()] += alpha * distance;
				}
			}

			// semantic error
			if (beta > 0 || gamma > 0) {
				std::vector<unsigned char> new_face_label(new_faces.size(), LABEL_OTHER);
				Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
				bool has_mesh_label;
				boost::tie(mesh_label, has_mesh_label) = profile.surface_mesh().property_map<Surface_mesh::Face_index, unsigned char>("label");
				assert(has_mesh_label);

				std::set<std::size_t> faces_with_no_label;
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					K::FT min_point = 0;
					if (min_point_per_area > 0) {
						min_point = min_point_per_area * CGAL::sqrt(new_faces[face_id].squared_area());
					}	
					if (points_in_new_face[face_id].size() > 0) {

						if (points_in_new_face[face_id].size() > min_point) {

							int face_label[LABELS.size()] = {0};

							for (auto point: points_in_new_face[face_id]) {
								face_label[point_cloud_label[point]]++;
							}

							auto argmax = std::max_element(face_label, face_label+LABELS.size());
							count_semantic_error += points_in_new_face[face_id].size() - *argmax;
							new_face_cost[face_id] += beta * (points_in_new_face[face_id].size() - *argmax);
							new_face_label[face_id] = argmax - face_label;

						} else { // We don't have information about this face.
							new_face_label[face_id] = LABEL_OTHER;
							count_semantic_error += points_in_new_face[face_id].size();
							new_face_cost[face_id] += beta * points_in_new_face[face_id].size();
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

					for(auto face_id: faces_with_no_label) {
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
					for (auto face_id: face_to_be_removed) {
						faces_with_no_label.erase(face_id);
					}
					face_to_be_removed.clear();
				}

				// semantic border length error
				if (gamma > 0) {
					for (auto h: profile.surface_mesh().halfedges_around_target(profile.v1_v0())) {
						if (!profile.surface_mesh().is_border(Surface_mesh::Edge_index(h))) {
							if (mesh_label[profile.surface_mesh().face(h)] != mesh_label[profile.surface_mesh().face(profile.surface_mesh().opposite(h))]) {
								semantic_border_length -= CGAL::sqrt(K::Vector_3(get(profile.vertex_point_map(), profile.surface_mesh().source(h)), get(profile.vertex_point_map(), profile.surface_mesh().target(h))).squared_length());
							}
						}
					}

					for (auto h: profile.surface_mesh().halfedges_around_target(profile.v0_v1())) {
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

				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					CollapseDataElement r;
					r.halfedge = new_faces_border_halfedge[face_id];
					r.label = new_face_label[face_id];
					r.cost = new_face_cost[face_id];
					for(auto ph: points_in_new_face[face_id]) r.points.push_back(ph);
					collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.push_back(r);
				}
			} else {
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					CollapseDataElement r;
					r.halfedge = new_faces_border_halfedge[face_id];
					r.cost = new_face_cost[face_id];
					for(auto ph: points_in_new_face[face_id]) r.points.push_back(ph);
					collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].elements.push_back(r);
				}
			}
		}

		return result_type(- old_cost + alpha * squared_distance + beta * count_semantic_error + gamma * semantic_border_length + delta * collapse_datas[Surface_mesh::Edge_index(profile.v0_v1())].cost);
	}

	return result_type();
}

Cost_stop_predicate::Cost_stop_predicate(const float cost) : cost(cost) {}

bool Cost_stop_predicate::operator()(const SMS::Edge_profile<Surface_mesh>::FT & current_cost, const SMS::Edge_profile<Surface_mesh> &, const SMS::Edge_profile<Surface_mesh>::edges_size_type, const SMS::Edge_profile<Surface_mesh>::edges_size_type) const {
	return current_cost > cost;
}

My_visitor::My_visitor(const K::FT alpha, const K::FT beta, const K::FT min_point_per_area, Surface_mesh &mesh, const Surface_mesh_info &mesh_info, Point_set &point_cloud) : alpha(alpha), beta(beta), min_point_per_area(min_point_per_area), mesh(mesh), mesh_info(mesh_info), point_cloud(point_cloud) {
	bool created_collapse_datas;
	boost::tie(collapse_datas, created_collapse_datas) = mesh.add_property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");
}

void My_visitor::OnStarted (Surface_mesh&) {
	
	// Add label to face
	bool created_label;
	boost::tie(mesh_label, created_label) = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label", LABEL_OTHER);
	assert(created_label);

	// Get point_cloud_label
	bool has_point_cloud_label;
	boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_point_cloud_label);

	// Create point_in_face
	bool created_point_in_face;
	boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points", std::list<Point_set::Index>());

	// Create face_costs
	bool created_face_costs;
	boost::tie(face_costs, created_face_costs) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:cost", 0);

	if(created_point_in_face || created_face_costs) {
		AABB_tree mesh_tree;
		PMP::build_AABB_tree(mesh, mesh_tree);
		for(auto ph: point_cloud) {
			auto p = type_converter(point_cloud.point(ph));
			auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
			point_in_face[location.first].push_back(ph);
			if (created_face_costs) face_costs[location.first] += alpha * CGAL::squared_distance(p, PMP::construct_point(location, mesh));
		}
		if (created_face_costs) {
			for(auto face: mesh.faces()) {
				K::FT min_point = 0;
				if (min_point_per_area > 0) {
					auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
					K::FT face_area = CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
					min_point = min_point_per_area * face_area;
				}	
				if (point_in_face[face].size() > 0) {
					if (point_in_face[face].size() > min_point) {
						int face_label[LABELS.size()] = {0};
						for (auto point: point_in_face[face]) {
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
	}
	
	bool created_point_isborder;
	Point_set::Property_map<bool> isborder;
	boost::tie (isborder, created_point_isborder) = point_cloud.add_property_map<bool>("p:isborder", false);
	if(created_point_isborder) {
		// Set isBorder
		int N = 10;
		Point_tree point_tree(point_cloud.begin(), point_cloud.end(), Point_tree::Splitter(), TreeTraits(point_cloud.point_map()));
		Neighbor_search::Distance tr_dist(point_cloud.point_map());
		for (auto point: point_cloud) {
			Neighbor_search search(point_tree, point_cloud.point(point), N, 0, true, tr_dist, false);
			unsigned char l = point_cloud_label[search.begin()->first];
			for(auto it = search.begin() + 1; it != search.end() && !isborder[point]; ++it) {
				if (point_cloud_label[it->first] != l) isborder[point] = true;
			}
		}

		std::cout << "Border points found" << std::endl;
	}

	{// Save point cloud
		Point_set output_point_cloud(point_cloud);
		Point_set::Property_map<unsigned char> output_label;
		bool has_output_label;
		boost::tie(output_label, has_output_label) = output_point_cloud.property_map<unsigned char>("p:label");
		assert(has_output_label);

		Point_set::Property_map<bool> output_isborder;
		bool has_output_isborder;
		boost::tie(output_isborder, has_output_isborder) = output_point_cloud.property_map<bool>("p:isborder");
		assert(has_output_isborder);
	
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

		for (auto point: output_point_cloud) {
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

		for (auto point: output_point_cloud) {
			red[point] = LABELS[output_label[point]].red;
			green[point] = LABELS[output_label[point]].green;
			blue[point] = LABELS[output_label[point]].blue;
		}

		CGAL::IO::write_point_set("pc_with_color.ply", output_point_cloud);
	}

	Surface_mesh::Property_map<Surface_mesh::Edge_index, CollapseData> collapse_datas;
	bool created_collapse_datas;
	boost::tie(collapse_datas, created_collapse_datas) = mesh.add_property_map<Surface_mesh::Edge_index, CollapseData>("e:c_datas");

	add_label(mesh, point_cloud, min_point_per_area);
	mesh_info.save_mesh(mesh, "initial-mesh.ply");

	start_collecte = std::chrono::system_clock::now();
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

	mesh.remove_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>(point_in_face);
	mesh.remove_property_map<Surface_mesh::Face_index, K::FT>(face_costs);
	mesh.remove_property_map<Surface_mesh::Edge_index, CollapseData>(collapse_datas);
}

void My_visitor::OnCollected(const SMS::Edge_profile<Surface_mesh>&, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >&) {
	start_collapse = std::chrono::system_clock::now();
	i_collecte++;
	if (i_collecte%1000 == 0) {
		std::chrono::duration<double> diff = start_collapse - start_collecte;
		std::cout << "\rCollecte: " << i_collecte << "/" << mesh.number_of_edges() << " (" << ((int) (((float) i_collecte)/mesh.number_of_edges()*100)) << "%)" << " still " << (((float) mesh.number_of_edges() - i_collecte) * diff.count() / i_collecte) << "s" << " (" << (((float) i_collecte) / diff.count()) << " op/s)" << std::flush;
	}
}

void My_visitor::OnSelected (const SMS::Edge_profile<Surface_mesh>&, boost::optional< SMS::Edge_profile<Surface_mesh>::FT > cost, const SMS::Edge_profile<Surface_mesh>::edges_size_type initial_edge_count, const SMS::Edge_profile<Surface_mesh>::edges_size_type current_edge_count) {
	if (current_edge_count%100 == 0) {
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = end - start_collapse;
		std::cout << "\rCollapse: " << (initial_edge_count-current_edge_count) << "/" << initial_edge_count << " (" << ((int) (((float) (initial_edge_count-current_edge_count))/initial_edge_count*100)) << "%)" << " still " << (((float) current_edge_count) * diff.count() / (initial_edge_count-current_edge_count)) << "s" << " (" << (((float) (initial_edge_count-current_edge_count)) / diff.count()) << " op/s)";
		if (cost) {
			std::cout << " - cost: " << *cost << "     " << std::flush;
		}
	}

	if (cost) {
		int cost_id = -1;
		if(!output[++cost_id] && *cost > 100) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c100.ply");
		} else if(!output[++cost_id] && *cost > 10) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c10.ply");
		} else if(!output[++cost_id] && *cost > 1) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c1.ply");
		} else if(!output[++cost_id] && *cost > 0.3) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.03.ply");
		} else if(!output[++cost_id] && *cost > 0.2) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.02.ply");
		} else if(!output[++cost_id] && *cost > 0.1) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.01.ply");
		} else if(!output[++cost_id] && *cost > 0.001) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.001.ply");
		} else if(!output[++cost_id] && *cost > 0) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-c0.ply");
		}
		cost_id = 8;
		if(!output[++cost_id] && current_edge_count <= 1000000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-1000000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 250000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-250000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 100000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-100000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 50000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-50000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 40000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-40000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 39000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-39000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 10000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-10000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 5000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-5000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 2500) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-2500.ply");
		} else if(!output[++cost_id] && current_edge_count <= 2000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-2000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 1000) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-1000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 500) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-500.ply");
		} else if(!output[++cost_id] && current_edge_count <= 400) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-400.ply");
		} else if(!output[++cost_id] && current_edge_count <= 300) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-300.ply");
		} else if(!output[++cost_id] && current_edge_count <= 100) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-100.ply");
		} else if(!output[++cost_id] && current_edge_count <= 50) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-50.ply");
		} else if(!output[++cost_id] && current_edge_count <= 10) {
			output[cost_id] = true;
			mesh_info.save_mesh(mesh,"mesh-10.ply");
		}
	}

}

void My_visitor::OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&) {
	// Called when an edge is about to be collapsed and replaced by a vertex whose position is *placement
	for(auto face: profile.triangles()) {
		auto fh = mesh.face(mesh.halfedge(face.v0, face.v1));
		point_in_face[fh].clear();
		face_costs[fh] = 0;
	}
}

void My_visitor::OnCollapsed (const SMS::Edge_profile<Surface_mesh>& prof, const Surface_mesh::Vertex_index vd) {
	// Called when an edge has been collapsed and replaced by the vertex vd

	// change point_in_face and face_costs
	for (auto element: collapse_datas[Surface_mesh::Edge_index(prof.v0_v1())].elements) {
		auto face = mesh.face(element.halfedge);
		for (auto ph: element.points) point_in_face[face].push_back(ph);
		face_costs[face] = element.cost;
		if (beta > 0 || gamma > 0) mesh_label[face] = element.label;
	}
}
