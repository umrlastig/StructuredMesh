#include "header.hpp"
#include "raster.hpp"

#include <chrono>
#include <list>
#include <map>
#include <vector>
#include <set>
#include <limits>

#include <CGAL/Point_set_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/QP_models.h>
#include <CGAL/QP_functions.h>

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

namespace SMS = CGAL::Surface_mesh_simplification;

void save_mesh(const Surface_mesh &mesh, const Raster &raster, const char *filename);

std::list<std::pair <K::Vector_3, K::FT>> volume_preservation_and_optimisation (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<std::pair <K::Vector_3, K::FT>> result;

	SMS::Edge_profile<Surface_mesh>::Triangle_vector triangles = profile.triangles();
	for (auto triangle: triangles) {
		if (!((triangle.v0 == profile.v0() && (triangle.v1 == profile.v1() || triangle.v2 == profile.v1()))
			|| (triangle.v1 == profile.v0() && (triangle.v0 == profile.v1() || triangle.v2 == profile.v1()))
			|| (triangle.v2 == profile.v0() && (triangle.v0 == profile.v1() || triangle.v1 == profile.v1()))
		)) {
			Point_3 p0 = get(profile.vertex_point_map(),triangle.v0);
			Point_3 p1 = get(profile.vertex_point_map(),triangle.v1);
			Point_3 p2 = get(profile.vertex_point_map(),triangle.v2);

			if (!CGAL::collinear (p0, p1, p2)) {
				/*std::cerr << "triangle: " << p0 << ", " << p1 << ", " << p2 << "\n";
				std::cerr << "normal: " << CGAL::normal(p0, p1, p2) << "\n";
				std::cerr << "determinant: " << CGAL::determinant(K::Vector_3(Point_3(CGAL::ORIGIN), p0), K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p2)) << "\n";*/
				
				K::Vector_3 n = CGAL::normal(p0, p1, p2) / 6;
				K::FT det = CGAL::scalar_product (n, K::Vector_3(Point_3(CGAL::ORIGIN), p0));

				result.push_back(std::make_pair(n, det));
			}
		} else {
			//std::cerr << "no way" << "\n";
		}
	}

	return result;
}

std::list<std::pair <K::Vector_3, K::Vector_3>> boundary_preservation_and_optimisation (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<std::pair <K::Vector_3, K::Vector_3>> result;

	auto length = CGAL::sqrt(K::Vector_3(profile.p0(), profile.p1()).squared_length());

	for (auto edge: profile.border_edges()) {
		Point_3 p0 = get(profile.vertex_point_map(), profile.surface_mesh().source(edge));
		Point_3 p1 = get(profile.vertex_point_map(), profile.surface_mesh().target(edge));

		//std::cerr << "edge: " << p0 << ", " << p1 << "\n";

		auto e1 = K::Vector_3(p0, p1) * length / 2;
		auto e2 = CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0)) * length / 2;

		/*std::cerr << "e1: " << K::Vector_3(p0, p1) << "\n";
		std::cerr << "e2: " << CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0)) << "\n";*/

		result.push_back(std::make_pair(e1, e2));
	}

	return result;
}

std::pair<std::list<K::Vector_3>, K::FT> triangle_shape_optimization (const SMS::Edge_profile<Surface_mesh>& profile) {

	std::list<K::Vector_3> result;

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();

	for (auto v: profile.link()) {
		result.push_back(K::Vector_3(K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), v))*squared_length);
	}

	return std::make_pair(result, squared_length);

}

std::list<std::pair <K::Vector_3, K::FT>> label_preservation (const SMS::Edge_profile<Surface_mesh>& profile, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) {

	std::list<std::pair <K::Vector_3, K::FT>> result;

	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	std::set<Point_set::Index> points_in_faces;
	int count_collapse_label[LABELS.size()] = {0};
	for(auto face: profile.triangles()) {
		auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
		points_in_faces.insert(point_in_face[fh].begin(), point_in_face[fh].end());
		int count_face_label[LABELS.size()] = {0};
		for (auto ph: point_in_face[fh]) {
			count_face_label[int(label[ph])]++;
		}
		auto argmax = std::max_element(count_face_label, count_face_label+LABELS.size());
		count_collapse_label[argmax - count_face_label]++;
	}

	int count_diff_label = 0;
	for (int i = 0; i < LABELS.size(); i++) {
		if (count_collapse_label[i] > 0) count_diff_label++;
	}

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

	if (count_diff_label < 2) {
		return result;
	} else if (count_diff_label == 2) {
		
		unsigned char label1, label2;
		bool first = true;
		for (int i = 0; i < LABELS.size(); i++) {
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
		for (auto point: points_in_faces) {
			if (label[point] == label1) {
				points_for_svm.push_back(point);
				y.push_back(1);
			} else if (label[point] == label2) {
				points_for_svm.push_back(point);
				y.push_back(-1);
			}
		}

		Exact_predicates_kernel::FT c = 5;

		Program qp (CGAL::EQUAL, true, 0, true, c);
		for (int i = 0; i < points_for_svm.size(); i++) {
			auto vr = Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) * y[i];
			qp.set_d(i, i, CGAL::scalar_product(vr,vr));
			for (int j = 0; j < i; j++) {
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
		for (int i = 0; i < points_for_svm.size(); i++) {
			w += y[i] * CGAL::to_double(*(value++)) * Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]));
		}

		Exact_predicates_kernel::FT b = 0;
		int count = 0;
		Exact_predicates_kernel::FT min_positive = std::numeric_limits<Exact_predicates_kernel::FT>::max();
		Exact_predicates_kernel::FT max_negative = std::numeric_limits<Exact_predicates_kernel::FT>::lowest();
		for (int i = 0; i < points_for_svm.size(); i++) {
			Exact_predicates_kernel::FT v = CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
			if (y[i] > 0 && v < min_positive) min_positive = v;
			if (y[i] < 0 && v > max_negative) max_negative = v;
		}
		if (max_negative < min_positive) {
			b += - (max_negative + min_positive) / 2;
			count = 1;
		} else {
			std::cerr << "to_be_test " << max_negative << " " << min_positive << "\n";
		}
		if (count == 0) {
			value = s.variable_values_begin();
			for (int i = 0; i < points_for_svm.size(); i++) {
				double v = CGAL::to_double(*(value++));
				if (v > 0 && v < c) {
					b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
					count++;
				}
			}
		}
		if (count == 0) {
			value = s.variable_values_begin();
			for (int i = 0; i < points_for_svm.size(); i++) {
				if (CGAL::to_double(*(value++)) > 0) {
					b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
					count++;
				}
			}
		}
		if (count == 0) {
			for (int i = 0; i < points_for_svm.size(); i++) {
				b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
				count++;
			}
		}
		b /= count;
		
		Exact_predicates_kernel::FT w_length = CGAL::sqrt(w.squared_length());
		w /= w_length;
		b /= w_length;

		std::cerr << "w: " << w << "\n";
		std::cerr << "b: " << b << "\n";

		for (int i = 0; i < points_for_svm.size(); i++) {
			std::cerr << "y" << i << ":\t" << y[i] << "\t";
			std::cerr << Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) << "   \t";
			std::cerr << CGAL::to_double(*(s.variable_values_begin() + i)) << "\t";
			std::cerr << CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) << "  \t";
			std::cerr << y[i]*(CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) +b) << "\n";
			if (y[i]*(CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) +b) < 0) {
				std::cerr << "Error\n";
			}
		}

		result.push_back(std::pair<K::Vector_3, K::FT>(type_converter(w)*squared_length, type_converter(b)*squared_length));
	} else {

		for (int i = 0; i < LABELS.size(); i++) {
			if (count_collapse_label[i] > 0) {
				
				std::vector<int> y;
				std::vector<Point_set::Index> points_for_svm;
				y.reserve(points_in_faces.size());
				points_for_svm.reserve(points_in_faces.size());
				for (auto point: points_in_faces) {
					if (label[point] == i) {
						points_for_svm.push_back(point);
						y.push_back(1);
					} else {
						points_for_svm.push_back(point);
						y.push_back(-1);
					}
				}

				Exact_predicates_kernel::FT c = 5;

				Program qp (CGAL::EQUAL, true, 0, true, c);
				for (int i = 0; i < points_for_svm.size(); i++) {
					auto vr = Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) * y[i];
					qp.set_d(i, i, CGAL::scalar_product(vr,vr));
					for (int j = 0; j < i; j++) {
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
				for (int i = 0; i < points_for_svm.size(); i++) {
					w += y[i] * CGAL::to_double(*(value++)) * Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]));
				}

				Exact_predicates_kernel::FT b = 0;
				int count = 0;
				Exact_predicates_kernel::FT min_positive = std::numeric_limits<Exact_predicates_kernel::FT>::max();
				Exact_predicates_kernel::FT max_negative = std::numeric_limits<Exact_predicates_kernel::FT>::lowest();
				for (int i = 0; i < points_for_svm.size(); i++) {
					Exact_predicates_kernel::FT v = CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
					if (y[i] > 0 && v < min_positive) min_positive = v;
					if (y[i] < 0 && v > max_negative) max_negative = v;
				}
				if (max_negative < min_positive) {
					b += - (max_negative + min_positive) / 2;
					count = 1;
				} else {
					std::cerr << "to_be_test " << max_negative << " " << min_positive << "\n";
				}
				if (count == 0) {
					value = s.variable_values_begin();
					for (int i = 0; i < points_for_svm.size(); i++) {
						double v = CGAL::to_double(*(value++));
						if (v > 0 && v < c) {
							b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
							count++;
						}
					}
				}
				if (count == 0) {
					value = s.variable_values_begin();
					for (int i = 0; i < points_for_svm.size(); i++) {
						if (CGAL::to_double(*(value++)) > 0) {
							b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
							count++;
						}
					}
				}
				if (count == 0) {
					for (int i = 0; i < points_for_svm.size(); i++) {
						b += y[i] - CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i])));
						count++;
					}
				}
				b /= count;
				
				Exact_predicates_kernel::FT w_length = CGAL::sqrt(w.squared_length());
				w /= w_length;
				b /= w_length;

				std::cerr << "w: " << w << "\n";
				std::cerr << "b: " << b << "\n";

				for (int i = 0; i < points_for_svm.size(); i++) {
					std::cerr << "y" << i << ":\t" << y[i] << "\t";
					std::cerr << Exact_predicates_kernel::Vector_3(CGAL::ORIGIN, point_cloud.point(points_for_svm[i])) << "   \t";
					std::cerr << CGAL::to_double(*(s.variable_values_begin() + i)) << "\t";
					std::cerr << CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) << "  \t";
					std::cerr << y[i]*(CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) +b) << "\n";
					if (y[i]*(CGAL::scalar_product(w, Exact_predicates_kernel::Vector_3(CGAL::ORIGIN,point_cloud.point(points_for_svm[i]))) +b) < 0) {
						std::cerr << "Error\n";
					}
				}

				result.push_back(std::pair<K::Vector_3, K::FT>(type_converter(w)*squared_length, type_converter(b)*squared_length));
			}
		}
	}

	return result;

}

std::pair<K::Vector_3, K::FT> middle_point (const SMS::Edge_profile<Surface_mesh>& profile) {

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();

	return std::make_pair((K::Vector_3(CGAL::ORIGIN, profile.p0()) + K::Vector_3(CGAL::ORIGIN, profile.p1())) * squared_length / 2, squared_length);

}

struct LindstromTurk_param {
	float volume_preservation;
	float boundary_preservation;
	float volume_optimisation;
	float boundary_optimization;
	float triangle_shape_preservation;
	float triangle_shape_optimization;
	float label_preservation;

	LindstromTurk_param(float volume_preservation,
						float boundary_preservation,
						float volume_optimisation,
						float boundary_optimization,
						float triangle_shape_preservation,
						float triangle_shape_optimization,
						float label_preservation) :
						volume_preservation(volume_preservation),
						boundary_preservation(boundary_preservation),
						volume_optimisation(volume_optimisation),
						boundary_optimization(boundary_optimization),
						triangle_shape_preservation(triangle_shape_preservation),
						triangle_shape_optimization(triangle_shape_optimization),
						label_preservation(label_preservation) {}
};

class Custom_placement {
	const LindstromTurk_param &params;
	std::map<Surface_mesh::Halfedge_index, K::FT> &costs;
	const Point_set &point_cloud;
	std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face;

	public:
		Custom_placement (const LindstromTurk_param &params, std::map<Surface_mesh::Halfedge_index, K::FT> &costs, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) : params(params), costs(costs), point_cloud(point_cloud), point_in_face(point_in_face) {}

		boost::optional<SMS::Edge_profile<Surface_mesh>::Point> operator()(const SMS::Edge_profile<Surface_mesh>& profile) const {
			typedef boost::optional<SMS::Edge_profile<Surface_mesh>::Point> result_type;

			auto r1 = volume_preservation_and_optimisation(profile);
			auto r2 = boundary_preservation_and_optimisation(profile);
			auto r3 = triangle_shape_optimization(profile);
			auto r4 = label_preservation(profile, point_cloud, point_in_face);
			auto r5 = middle_point(profile);

			Eigen_vector B(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_preservation > 0) ? 3 : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.first.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0));
			Eigen_matrix A(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_preservation > 0) ? 3 : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.first.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0), 3);

			std::size_t i = 0;

			// Volume preservation
			if (params.volume_preservation > 0) {
				K::Vector_3 n (CGAL::NULL_VECTOR);
				K::FT det (0);
				for (auto vpo: r1) {
					n += vpo.first;
					det += vpo.second;
				}
				A.set(i, 0, n.x()*params.volume_preservation);
				A.set(i, 1, n.y()*params.volume_preservation);
				A.set(i, 2, n.z()*params.volume_preservation);
				B.set(i++, det*params.volume_preservation);
			}

			// Volume optimisation
			if (params.volume_optimisation > 0) {
				for (auto vpo: r1) {
					A.set(i, 0, vpo.first.x()*params.volume_optimisation);
					A.set(i, 1, vpo.first.y()*params.volume_optimisation);
					A.set(i, 2, vpo.first.z()*params.volume_optimisation);
					B.set(i++, vpo.second*params.volume_optimisation);
				}
			}

			// Boundary preservation
			if (params.boundary_preservation > 0) {
				K::Vector_3 e1 (CGAL::NULL_VECTOR);
				K::Vector_3 e2 (CGAL::NULL_VECTOR);
				for (auto vpo: r2) {
					e1 += vpo.first;
					e2 += vpo.second;
				}
				A.set(i, 0, 0);
				A.set(i, 1, -e1.z()*params.boundary_preservation);
				A.set(i, 2, e1.y()*params.boundary_preservation);
				B.set(i++, e2.x()*params.boundary_preservation);
				A.set(i, 0, e1.z()*params.boundary_preservation);
				A.set(i, 1, 0);
				A.set(i, 2, -e1.x()*params.boundary_preservation);
				B.set(i++, e2.y()*params.boundary_preservation);
				A.set(i, 0, -e1.y()*params.boundary_preservation);
				A.set(i, 1, e1.x()*params.boundary_preservation);
				A.set(i, 2, 0);
				B.set(i++, e2.z()*params.boundary_preservation);
			}

			// Boundary optimisation
			if (params.boundary_optimization > 0) {
				for (auto vpo: r2) {
					A.set(i, 0, 0);
					A.set(i, 1, -vpo.first.z()*params.boundary_optimization);
					A.set(i, 2, vpo.first.y()*params.boundary_optimization);
					B.set(i++, vpo.second.x()*params.boundary_optimization);
					A.set(i, 0, vpo.first.z()*params.boundary_optimization);
					A.set(i, 1, 0);
					A.set(i, 2, -vpo.first.x()*params.boundary_optimization);
					B.set(i++, vpo.second.y()*params.boundary_optimization);
					A.set(i, 0, -vpo.first.y()*params.boundary_optimization);
					A.set(i, 1, vpo.first.x()*params.boundary_optimization);
					A.set(i, 2, 0);
					B.set(i++, vpo.second.z()*params.boundary_optimization);
				}
			}

			// Triange shape preservation
			if (params.triangle_shape_preservation > 0) {
				K::Vector_3 e (CGAL::NULL_VECTOR);
				for (auto vpo: r3.first) {
					e += vpo;
				}
				e /= r3.first.size();
				A.set(i, 0, r3.second*params.triangle_shape_preservation);
				A.set(i, 1, 0);
				A.set(i, 2, 0);
				B.set(i++, e.x()*params.triangle_shape_preservation);
				A.set(i, 0, 0);
				A.set(i, 1, r3.second*params.triangle_shape_preservation);
				A.set(i, 2, 0);
				B.set(i++, e.y()*params.triangle_shape_preservation);
				A.set(i, 0, 0);
				A.set(i, 1, 0);
				A.set(i, 2, r3.second*params.triangle_shape_preservation);
				B.set(i++, e.z()*params.triangle_shape_preservation);
			}

			// Triange shape optimisation
			if (params.triangle_shape_optimization > 0) {
				for (auto vpo: r3.first) {
					A.set(i, 0, r3.second*params.triangle_shape_optimization);
					A.set(i, 1, 0);
					A.set(i, 2, 0);
					B.set(i++, vpo.x()*params.triangle_shape_optimization);
					A.set(i, 0, 0);
					A.set(i, 1, r3.second*params.triangle_shape_optimization);
					A.set(i, 2, 0);
					B.set(i++, vpo.y()*params.triangle_shape_optimization);
					A.set(i, 0, 0);
					A.set(i, 1, 0);
					A.set(i, 2, r3.second*params.triangle_shape_optimization);
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

			// Solve AX=B
			auto C = B;
			CGAL::Eigen_svd::solve(A, B);

			auto R = A*B - C;
			
			// Save cost
			Point_3 placement(B.vector()[0], B.vector()[1], B.vector()[2]);
			costs[profile.v0_v1()] = (R.transpose()*R)(0,0);

			return result_type(placement);
		}
};

class Custom_cost {
	std::map<Surface_mesh::Halfedge_index, K::FT> &costs;

	public:
		Custom_cost (std::map<Surface_mesh::Halfedge_index, K::FT> &costs) : costs(costs) {}

		boost::optional<SMS::Edge_profile<Surface_mesh>::FT> operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
			typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;

			if (placement) {
				return result_type(costs[profile.v0_v1()]);
			}

			return result_type();
		}
};

class Cost_stop_predicate {
	public:

		Cost_stop_predicate(const float cost) : cost(cost) {}

		bool operator()(const SMS::Edge_profile<Surface_mesh>::FT & current_cost, const SMS::Edge_profile<Surface_mesh> &, const SMS::Edge_profile<Surface_mesh>::edges_size_type, const SMS::Edge_profile<Surface_mesh>::edges_size_type) const {
			return current_cost > cost;
		}

	private:
		const float cost;
};

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh> {
	private:
		int i_collecte = 0;
		Surface_mesh &mesh;
		const Raster &raster;
		std::chrono::time_point<std::chrono::system_clock> start_collecte;
		std::chrono::time_point<std::chrono::system_clock> start_collapse;
		bool output[25] = {false};

		const Point_set &point_cloud;
		std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face;
		std::set<Point_set::Index> points_to_be_change;

	public:
		My_visitor(Surface_mesh &mesh, const Raster &raster, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) : mesh(mesh), raster(raster), point_cloud(point_cloud), point_in_face(point_in_face) {}

		void OnStarted (Surface_mesh &mesh) {
			start_collecte = std::chrono::system_clock::now();
		}

		void OnCollected(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >& cost) {
			start_collapse = std::chrono::system_clock::now();
			i_collecte++;
			if (i_collecte%1000 == 0) {
				std::chrono::duration<double> diff = start_collapse - start_collecte;
				std::cout << "\rCollecte: " << i_collecte << "/" << mesh.number_of_edges() << " (" << ((int) (((float) i_collecte)/mesh.number_of_edges()*100)) << "%)" << " still " << (((float) mesh.number_of_edges() - i_collecte) * diff.count() / i_collecte) << "s" << " (" << (((float) i_collecte) / diff.count()) << " op/s)" << std::flush;
			}
		}

		void OnSelected (const SMS::Edge_profile<Surface_mesh> &profile, boost::optional< SMS::Edge_profile<Surface_mesh>::FT > cost, const SMS::Edge_profile<Surface_mesh>::edges_size_type initial_edge_count, const SMS::Edge_profile<Surface_mesh>::edges_size_type current_edge_count) {
			if (current_edge_count%100 == 0) {
				auto end = std::chrono::system_clock::now();
				std::chrono::duration<double> diff = end - start_collapse;
				std::cout << "\rCollapse: " << (initial_edge_count-current_edge_count) << "/" << initial_edge_count << " (" << ((int) (((float) (initial_edge_count-current_edge_count))/initial_edge_count*100)) << "%)" << " still " << (((float) current_edge_count) * diff.count() / (initial_edge_count-current_edge_count)) << "s" << " (" << (((float) (initial_edge_count-current_edge_count)) / diff.count()) << " op/s)";
				if (cost) {
					std::cout << " - cost: " << *cost << "     " << std::flush;
				}
			}

			if (cost) {
				if(*cost > 100 && !output[0]) {
					output[0] = true;
					save_mesh(mesh,raster,"mesh-100.ply");
				} else if(*cost > 10 && !output[0]) {
					output[0] = true;
					save_mesh(mesh,raster,"mesh-10.ply");
				} else if(*cost > 0 && !output[1]) {
					output[1] = true;
					save_mesh(mesh,raster,"mesh-0.ply");
				}
				if(!output[2] && current_edge_count <= 1000000) {
					output[2] = true;
					save_mesh(mesh,raster,"mesh-1000000.ply");
				} else if(!output[3] && current_edge_count <= 250000) {
					output[3] = true;
					save_mesh(mesh,raster,"mesh-250000.ply");
				} else if(!output[4] && current_edge_count <= 100000) {
					output[4] = true;
					save_mesh(mesh,raster,"mesh-100000.ply");
				} else if(!output[11] && current_edge_count <= 50000) {
					output[11] = true;
					save_mesh(mesh,raster,"mesh-50000.ply");
				} else if(!output[5] && current_edge_count <= 10000) {
					output[5] = true;
					save_mesh(mesh,raster,"mesh-10000.ply");
				} else if(!output[10] && current_edge_count <= 5000) {
					output[10] = true;
					save_mesh(mesh,raster,"mesh-5000.ply");
				} else if(!output[6] && current_edge_count <= 1000) {
					output[6] = true;
					save_mesh(mesh,raster,"mesh-1000.ply");
				} else if(!output[7] && current_edge_count <= 100) {
					output[7] = true;
					save_mesh(mesh,raster,"mesh-100.ply");
				} else if(!output[8] && current_edge_count <= 50) {
					output[8] = true;
					save_mesh(mesh,raster,"mesh-50.ply");
				} else if(!output[9] && current_edge_count <= 10) {
					output[9] = true;
					save_mesh(mesh,raster,"mesh-10.ply");
				}
			}

		};

		void Collapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) {
			// Called when an edge is about to be collapsed and replaced by a vertex whose position is *placement
			for(auto face: profile.triangles()) {
				auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
				points_to_be_change.insert(point_in_face[fh].begin(), point_in_face[fh].end());
				point_in_face[fh].clear();
			}
		};

		void OnCollapsed (const SMS::Edge_profile<Surface_mesh> &profile, const Surface_mesh::Vertex_index vd) {
			// Called when an edge has been collapsed and replaced by the vertex vd
			CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

			for(auto ph: points_to_be_change) {
				K::FT min_d = std::numeric_limits<K::FT>::max();
				Surface_mesh::Face_index nearest_face;
				for(auto face: profile.surface_mesh().faces_around_target(profile.surface_mesh().halfedge(vd))) {
					auto r = profile.surface_mesh().vertices_around_face(profile.surface_mesh().halfedge(face)).begin();
					auto d = CGAL::squared_distance(K::Triangle_3(profile.surface_mesh().point(*r++), profile.surface_mesh().point(*r++), profile.surface_mesh().point(*r)), type_converter(point_cloud.point(ph)));
					if (d < min_d) {
						min_d = d;
						nearest_face = face;
					}
				}
				if (nearest_face.is_valid()) {
					point_in_face[nearest_face].push_back(ph);
				}
			}
		};
};

std::tuple<Surface_mesh, Surface_mesh> compute_meshes(const Raster &raster) {

	std::cout << "Terrain mesh" << std::endl;
	Surface_mesh terrain_mesh;
	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> terrain_vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			terrain_vertex_index[L][P] = terrain_mesh.add_vertex(Point_3(0.5 + P, 0.5 + L, raster.dtm[L][P]));
		}
	}
	std::cout << "Point added" << std::endl;
	// Add faces
	for (int L = 0; L < raster.ySize-1; L++) {
		for (int P = 0; P < raster.xSize-1; P++) {
			if (pow(raster.dtm[L][P]-raster.dtm[L+1][P+1], 2) < pow(raster.dtm[L+1][P]-raster.dtm[L][P+1], 2)) {
				terrain_mesh.add_face(terrain_vertex_index[L][P], terrain_vertex_index[L+1][P+1], terrain_vertex_index[L+1][P]);
				terrain_mesh.add_face(terrain_vertex_index[L][P], terrain_vertex_index[L][P+1], terrain_vertex_index[L+1][P+1]);
			} else {
				terrain_mesh.add_face(terrain_vertex_index[L][P], terrain_vertex_index[L][P+1], terrain_vertex_index[L+1][P]);
				terrain_mesh.add_face(terrain_vertex_index[L][P+1], terrain_vertex_index[L+1][P+1], terrain_vertex_index[L+1][P]);
			}
		}
	}
	std::cout << "Faces added" << std::endl;

	save_mesh(terrain_mesh, raster, "initial-terrain-mesh.ply");

	SMS::edge_collapse(terrain_mesh, Cost_stop_predicate(10));
	std::cout << "Terrain mesh simplified" << std::endl;

	save_mesh(terrain_mesh, raster, "terrain-mesh.ply");

	std::cout << "Surface mesh" << std::endl;
	Surface_mesh mesh;
	std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	
	Point_set point_cloud;
	bool created_label;
	Point_set::Property_map<unsigned char> label;
	boost::tie (label, created_label) = point_cloud.add_property_map<unsigned char>("p:label", 0);
	assert(created_label);

	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	std::vector<std::vector<Point_set::Index>> point_index(raster.ySize, std::vector<Point_set::Index>(raster.xSize, Point_set::Index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			vertex_index[L][P] = mesh.add_vertex(Point_3(0.5 + P, 0.5 + L, raster.dsm[L][P]));
			point_index[L][P] = *(point_cloud.insert(Point_set::Point_3(0.5 + P, 0.5 + L, raster.dsm[L][P])));
			label[point_index[L][P]] = raster.land_cover[L][P];
		}
	}
	std::cout << "Point added" << std::endl;
	// Add faces
	for (int L = 0; L < raster.ySize-1; L++) {
		for (int P = 0; P < raster.xSize-1; P++) {
			if (pow(raster.dsm[L][P]-raster.dsm[L+1][P+1], 2) < pow(raster.dsm[L+1][P]-raster.dsm[L][P+1], 2)) {
				auto f1 = mesh.add_face(vertex_index[L][P], vertex_index[L+1][P+1], vertex_index[L+1][P]);
				point_in_face[f1].push_back(point_index[L][P]);
				point_in_face[f1].push_back(point_index[L+1][P+1]);
				point_in_face[f1].push_back(point_index[L+1][P]);
				auto f2 = mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P+1]);
				point_in_face[f2].push_back(point_index[L][P]);
				point_in_face[f2].push_back(point_index[L][P+1]);
				point_in_face[f2].push_back(point_index[L+1][P+1]);
			} else {
				auto f1 = mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P]);
				point_in_face[f1].push_back(point_index[L][P]);
				point_in_face[f1].push_back(point_index[L][P+1]);
				point_in_face[f1].push_back(point_index[L+1][P]);
				auto f2 = mesh.add_face(vertex_index[L][P+1], vertex_index[L+1][P+1], vertex_index[L+1][P]);
				point_in_face[f2].push_back(point_index[L][P+1]);
				point_in_face[f2].push_back(point_index[L+1][P+1]);
				point_in_face[f2].push_back(point_index[L+1][P]);
			}
		}
	}
	std::cout << "Faces added" << std::endl;

	save_mesh(mesh, raster, "initial-mesh.ply");

	//Cost_stop_predicate stop(10);
	SMS::Count_stop_predicate<Surface_mesh> stop(50);
	const LindstromTurk_param params (1,1,1,1,0,0,0.01);
	std::map<Surface_mesh::Halfedge_index, K::FT> costs;
	Custom_placement pf(params, costs, point_cloud, point_in_face);
	Custom_cost cf(costs);
	SMS::Bounded_normal_change_filter<> filter;
	int r = SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(My_visitor(mesh, raster, point_cloud, point_in_face)));
	std::cout << "\rMesh simplified                                               " << std::endl;

	save_mesh(mesh, raster, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, mesh);
}
