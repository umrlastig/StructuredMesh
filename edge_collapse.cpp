#include "edge_collapse.hpp"

void add_label(Surface_mesh &mesh, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face, int min_surface) {
	//Update mesh label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	bool has_mesh_label;
	boost::tie(mesh_label, has_mesh_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_mesh_label);

	Point_set::Property_map<unsigned char> point_cloud_label;
	bool has_point_cloud_label;
	boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_point_cloud_label);

	// Compute face_area and mean_point_per_area to filter face with less point than 1 / min_surface of mean point
	std::map<Surface_mesh::Face_index, K::FT> face_area;
	K::FT mean_point_per_area = 0;
	if (min_surface > 0) {
		K::FT point_per_area = 0;
		for(auto face: mesh.faces()) {
			CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
			boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
			auto p0 = mesh.point(*(vbegin++));
			auto p1 = mesh.point(*(vbegin++));
			auto p2 = mesh.point(*(vbegin++));
			K::FT area = CGAL::sqrt(K::Triangle_3(p0, p1, p2).squared_area());
			point_per_area += point_in_face[face].size() / area;
			face_area[face] = area;
		}
		mean_point_per_area = point_per_area / mesh.number_of_faces();
	}

	for(auto face: mesh.faces()) {
		int face_label[LABELS.size()] = {0};

		for (auto point: point_in_face[face]) {
			face_label[point_cloud_label[point]]++;
		}

		int min_point = 0;
		if (min_surface > 0) {
			min_point = int(mean_point_per_area * face_area[face] / min_surface);
		}

		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		if (*argmax > min_point) {
			mesh_label[face] = argmax - face_label;
		} else {
			mesh_label[face] = LABEL_OTHER;
		}
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

std::list<std::pair <K::Vector_3, K::FT>> label_preservation (const SMS::Edge_profile<Surface_mesh>& profile, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) {

	std::list<std::pair <K::Vector_3, K::FT>> result;

	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	Point_set::Property_map<bool> isborder;
	bool has_isborder;
	boost::tie(isborder, has_isborder) = point_cloud.property_map<bool>("p:isborder");
	assert(has_isborder);

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

LindstromTurk_param::LindstromTurk_param(
	float volume_preservation,
	float boundary_preservation,
	float volume_optimisation,
	float boundary_optimization,
	float triangle_shape_optimization,
	float label_preservation) :
	volume_preservation(volume_preservation),
	boundary_preservation(boundary_preservation),
	volume_optimisation(volume_optimisation),
	boundary_optimization(boundary_optimization),
	triangle_shape_optimization(triangle_shape_optimization),
	label_preservation(label_preservation) {}

Custom_placement::Custom_placement (const LindstromTurk_param &params, std::map<Surface_mesh::Halfedge_index, K::FT> &costs, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) : params(params), costs(costs), point_cloud(point_cloud), point_in_face(point_in_face) {}


boost::optional<SMS::Edge_profile<Surface_mesh>::Point> Custom_placement::operator()(const SMS::Edge_profile<Surface_mesh>& profile) const {
	typedef boost::optional<SMS::Edge_profile<Surface_mesh>::Point> result_type;

	auto r1 = volume_preservation_and_optimisation(profile);
	auto r2 = boundary_preservation_and_optimisation(profile);
	auto r3 = triangle_shape_optimization(profile);
	auto r4 = label_preservation(profile, point_cloud, point_in_face);

	Eigen_vector B(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0 && r2.size() > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0));
	Eigen_matrix A(((params.volume_preservation > 0) ? 1 : 0) + ((params.volume_optimisation > 0) ? r1.size() : 0) + ((params.boundary_preservation > 0 && r2.size() > 0) ? 3 : 0) + ((params.boundary_optimization > 0) ? 3*r2.size() : 0) + ((params.triangle_shape_optimization > 0) ? 3*r3.size() : 0) + ((params.label_preservation > 0) ? r4.size() : 0), 3);

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

	// Solve AX=B
	auto C = B;
	CGAL::Eigen_svd::solve(A, B);

	auto R = A*B - C;
	
	// Save cost
	Point_3 placement(B.vector()[0], B.vector()[1], B.vector()[2]);
	costs[profile.v0_v1()] = (R.transpose()*R)(0,0);

	return result_type(placement);
}

Custom_cost::Custom_cost (const K::FT alpha, const K::FT beta, const K::FT gamma, std::map<Surface_mesh::Halfedge_index, K::FT> &costs, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) : alpha(alpha), beta(beta), gamma(gamma), point_cloud(point_cloud), costs(costs), point_in_face(point_in_face) {}

boost::optional<SMS::Edge_profile<Surface_mesh>::FT> Custom_cost::operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
	typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

	Point_set::Property_map<unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_label);

	if (placement) {
		
		std::set<Point_set::Index> points_to_be_change;
		for(auto face: profile.triangles()) {
			auto fh = profile.surface_mesh().face(profile.surface_mesh().halfedge(face.v0, face.v1));
			points_to_be_change.insert(point_in_face[fh].begin(), point_in_face[fh].end());
		}

		std::vector<K::Triangle_3> new_faces;
		Point_3 C = *placement;
		for (auto he : CGAL::halfedges_around_source(profile.v0(), profile.surface_mesh())) {
			if (he != profile.v0_v1() && he != profile.v0_vR()) {
				Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
				Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
				new_faces.push_back(K::Triangle_3(A, B, C));
			}
		}
		for (auto he : CGAL::halfedges_around_source(profile.v1(), profile.surface_mesh())) {
			if (he != profile.v1_v0() && he != profile.v1_vL()) {
				Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
				Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
				new_faces.push_back(K::Triangle_3(A, B, C));
			}
		}

		// geometric error
		K::FT squared_distance = 0;
		std::vector<std::list<Point_set::Index>> points_in_new_face (new_faces.size());
		Triangle_tree tree(new_faces.begin(), new_faces.end());
		for(auto ph: points_to_be_change) {
			auto point = type_converter(point_cloud.point(ph));
			auto point_and_id = tree.closest_point_and_primitive(point);

			points_in_new_face[point_and_id.second - new_faces.begin()].push_back(ph);
			squared_distance += CGAL::squared_distance(point, point_and_id.first);
		}

		// semantic error
		int count_semantic_error = 0;
		for(std::size_t i = 0; i < new_faces.size(); i++) {
			int face_label[LABELS.size()] = {0};
			int sum_face_label = 0;
			for (auto ph: points_in_new_face.at(i)) {
				sum_face_label++;
				face_label[label[ph]]++;
			}
			auto argmax = std::max_element(face_label, face_label+LABELS.size());
			count_semantic_error += sum_face_label - *argmax;
		}

		return result_type(alpha * squared_distance + beta * count_semantic_error + gamma * costs[profile.v0_v1()]);
	}

	return result_type();
}

Cost_stop_predicate::Cost_stop_predicate(const float cost) : cost(cost) {}

bool Cost_stop_predicate::operator()(const SMS::Edge_profile<Surface_mesh>::FT & current_cost, const SMS::Edge_profile<Surface_mesh> &, const SMS::Edge_profile<Surface_mesh>::edges_size_type, const SMS::Edge_profile<Surface_mesh>::edges_size_type) const {
	return current_cost > cost;
}

My_visitor::My_visitor(Surface_mesh &mesh, const Surface_mesh_info &mesh_info, const Point_set &point_cloud, std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> &point_in_face) : mesh(mesh), mesh_info(mesh_info), point_cloud(point_cloud), point_in_face(point_in_face) {}

void My_visitor::OnStarted (Surface_mesh&) {
	start_collecte = std::chrono::system_clock::now();
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
		if(*cost > 100 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c100.ply");
		} else if(*cost > 10 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c10.ply");
		} else if(*cost > 0.3 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c0.03.ply");
		} else if(*cost > 0.2 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c0.02.ply");
		} else if(*cost > 0.1 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c0.01.ply");
		} else if(*cost > 0.001 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c0.001.ply");
		} else if(*cost > 0 && !output[++cost_id]) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-c0.ply");
		}
		cost_id = 8;
		if(!output[++cost_id] && current_edge_count <= 1000000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-1000000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 250000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-250000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 100000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-100000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 50000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-50000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 10000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-10000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 5000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-5000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 2500) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-2500.ply");
		} else if(!output[++cost_id] && current_edge_count <= 2000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-2000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 1000) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-1000.ply");
		} else if(!output[++cost_id] && current_edge_count <= 500) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-500.ply");
		} else if(!output[++cost_id] && current_edge_count <= 100) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-100.ply");
		} else if(!output[++cost_id] && current_edge_count <= 50) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-50.ply");
		} else if(!output[++cost_id] && current_edge_count <= 10) {
			output[cost_id] = true;
			add_label(mesh, point_cloud, point_in_face);
			mesh_info.save_mesh(mesh,"mesh-10.ply");
		}
	}

}

void My_visitor::OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&) {
	// Called when an edge is about to be collapsed and replaced by a vertex whose position is *placement
	for(auto face: profile.triangles()) {
		auto fh = mesh.face(mesh.halfedge(face.v0, face.v1));
		points_to_be_change.insert(point_in_face[fh].begin(), point_in_face[fh].end());
		point_in_face[fh].clear();
	}
}

void My_visitor::OnCollapsed (const SMS::Edge_profile<Surface_mesh>&, const Surface_mesh::Vertex_index vd) {
	// Called when an edge has been collapsed and replaced by the vertex vd
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

	for(auto ph: points_to_be_change) {
		K::FT min_d = std::numeric_limits<K::FT>::max();
		std::map<Surface_mesh::Face_index, K::FT> nearest_face;
		K::FT threshold = 0.001;
		for(auto face: mesh.faces_around_target(mesh.halfedge(vd))) {
			if (face != mesh.null_face()) {
				auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
				auto d = CGAL::squared_distance(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)), type_converter(point_cloud.point(ph)));
				if (d < min_d) {
					if (d < min_d - threshold) {
						nearest_face.clear();
					}
					nearest_face[face] = d;
					min_d = d;
				} else if (d < min_d + threshold) {
					nearest_face[face] = d;
				}
			}
		}
		for (auto face_d: nearest_face) {
			if (face_d.second < min_d + threshold) {
				assert(face_d.first.is_valid());
				point_in_face[face_d.first].push_back(ph);
			}
		}
	}

	points_to_be_change.clear();
}

std::tuple<Surface_mesh, Surface_mesh> compute_meshes(const Raster &raster, const Surface_mesh_info &mesh_info) {

	std::cout << "Terrain mesh" << std::endl;
	Surface_mesh terrain_mesh;
	
	double x, y;
	
	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> terrain_vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			raster.grid_to_coord(P, L, x, y);
			terrain_vertex_index[L][P] = terrain_mesh.add_vertex(Point_3(x - mesh_info.x_0, y - mesh_info.y_0, raster.dtm[L][P]));
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

	// Return mesh if coords are in reverse order
	double x_0, y_0, x_1, y_1;	
	raster.grid_to_coord(0, 0, x_0, y_0);
	raster.grid_to_coord(1, 1, x_1, y_1);
	if ((x_1-x_0)*(y_1-y_0) < 0) {
		CGAL::Polygon_mesh_processing::reverse_face_orientations(terrain_mesh); 	
	}

	mesh_info.save_mesh(terrain_mesh, "initial-terrain-mesh.ply");

	SMS::edge_collapse(terrain_mesh, Cost_stop_predicate(10));
	std::cout << "Terrain mesh simplified" << std::endl;

	mesh_info.save_mesh(terrain_mesh, "terrain-mesh.ply");

	std::cout << "Surface mesh" << std::endl;
	Surface_mesh mesh;
	std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;

	Point_set point_cloud;
	bool created_point_label;
	Point_set::Property_map<unsigned char> label;
	boost::tie (label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", 0);
	assert(created_point_label);
	bool created_point_isborder;
	Point_set::Property_map<bool> isborder;
	boost::tie (isborder, created_point_isborder) = point_cloud.add_property_map<bool>("p:isborder", false);
	assert(created_point_isborder);

	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	std::vector<std::vector<Point_set::Index>> point_index(raster.ySize, std::vector<Point_set::Index>(raster.xSize, Point_set::Index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			raster.grid_to_coord(P, L, x, y);
			vertex_index[L][P] = mesh.add_vertex(Point_3(x - mesh_info.x_0, y - mesh_info.y_0, raster.dsm[L][P]));
			point_index[L][P] = *(point_cloud.insert(Point_set::Point_3(x - mesh_info.x_0, y - mesh_info.y_0, raster.dsm[L][P])));
			label[point_index[L][P]] = raster.land_cover[L][P];
		}
	}
	std::cout << "Point added" << std::endl;

	// Set isBorder
	int N = 10;
	Point_tree tree(point_cloud.begin(), point_cloud.end(), Point_tree::Splitter(), TreeTraits(point_cloud.point_map()));
	Neighbor_search::Distance tr_dist(point_cloud.point_map());
	for (auto point: point_cloud) {
		Neighbor_search search(tree, point_cloud.point(point), N, 0, true, tr_dist, false);
		unsigned char l = label[search.begin()->first];
		for(auto it = search.begin() + 1; it != search.end() && !isborder[point]; ++it) {
			if (label[it->first] != l) isborder[point] = true;
		}
	}
	std::cout << "Border points found" << std::endl;

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

	// Return mesh if coords are in reverse order
	if ((x_1-x_0)*(y_1-y_0) < 0) {
		CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh); 	
	}

	auto created_label = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label",0);
	assert(created_label.second);

	add_label(mesh, point_cloud, point_in_face);
	mesh_info.save_mesh(mesh, "initial-mesh.ply");

	Cost_stop_predicate stop(5);
	//SMS::Count_stop_predicate<Surface_mesh> stop(50);
	const LindstromTurk_param params (1,1,1,1,0.1,1);
	std::map<Surface_mesh::Halfedge_index, K::FT> costs;
	Custom_placement pf(params, costs, point_cloud, point_in_face);
	Custom_cost cf(1, 1, 0.01, costs, point_cloud, point_in_face);
	SMS::Bounded_normal_change_filter<> filter;
	SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(My_visitor(mesh, mesh_info, point_cloud, point_in_face)));
	std::cout << "\rMesh simplified                                               " << std::endl;

	add_label(mesh, point_cloud, point_in_face, 2);
	mesh_info.save_mesh(mesh, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, mesh);
}
