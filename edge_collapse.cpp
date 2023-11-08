#include "edge_collapse.hpp"

K::FT add_label(Surface_mesh &mesh, const Point_set &point_cloud, int min_surface, K::FT mean_point_per_area) {
	//Update mesh label
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	bool has_mesh_label;
	boost::tie(mesh_label, has_mesh_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_mesh_label);

	Point_set::Property_map<unsigned char> point_cloud_label;
	bool has_point_cloud_label;
	boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_point_cloud_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = mesh.property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	// Compute face_area and mean_point_per_area to filter face with less point than 1 / min_surface of mean point
	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_area;

	if (min_surface > 0) {
		bool created_face_area;
		boost::tie(face_area, created_face_area) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:area", 0);
		if (created_face_area) {
			for(auto face: mesh.faces()) {
				CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
				boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
				auto p0 = mesh.point(*(vbegin++));
				auto p1 = mesh.point(*(vbegin++));
				auto p2 = mesh.point(*(vbegin++));
				face_area[face] = CGAL::sqrt(K::Triangle_3(p0, p1, p2).squared_area());
			}
		}

		if (mean_point_per_area == 0) {
			K::FT point_per_area = 0;
			for(auto face: mesh.faces()) {
				if (face_area[face] > 0) point_per_area += point_in_face[face].size() / face_area[face];
			}
			if (mesh.number_of_faces() > 0) mean_point_per_area = point_per_area / mesh.number_of_faces();
		}
	}

	std::set<Surface_mesh::Face_index> faces_with_no_label;
	for(auto face: mesh.faces()) {
		K::FT min_point = 0;
		if (min_surface > 0) {
			min_point = mean_point_per_area * face_area[face] / min_surface;
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
	for(auto face: faces_with_no_label) {
		K::FT face_label[LABELS.size()] = {0};
		for (auto he: mesh.halfedges_around_face(mesh.halfedge(face))) {
			if (!mesh.is_border(Surface_mesh::Edge_index(he))) {
				if (faces_with_no_label.count(mesh.face(mesh.opposite(he))) == 0) {
					face_label[mesh_label[mesh.face(mesh.opposite(he))]] += CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(he)), mesh.point(mesh.target(he))));
				}
			}
		}
		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		if (*argmax > 0) {
			mesh_label[face] = argmax - face_label;
		} else { // No neighbor with known label
			mesh_label[face] = LABEL_OTHER;
		}
	}

	return mean_point_per_area;
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

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points");
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

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points");
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
	bool has_placement_costs;
	boost::tie(placement_costs, has_placement_costs) = mesh.property_map<Surface_mesh::Halfedge_index, K::FT>("h:p_cost");
	assert(has_placement_costs);
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
	placement_costs[profile.v0_v1()] = (R.transpose()*R)(0,0);

	return result_type(placement);
}

Custom_cost::Custom_cost (const K::FT alpha, const K::FT beta, const K::FT gamma, const K::FT delta, int min_surface, K::FT mean_point_per_area, Surface_mesh &mesh, const Point_set &point_cloud) : alpha(alpha), beta(beta), gamma(gamma), delta(delta), min_surface(min_surface), mean_point_per_area(mean_point_per_area), point_cloud(point_cloud) {
	bool has_face_costs;
	boost::tie(face_costs, has_face_costs) = mesh.property_map<Surface_mesh::Face_index, K::FT>("f:cost");
	assert(has_face_costs);
	
	bool has_placement_costs;
	boost::tie(placement_costs, has_placement_costs) = mesh.property_map<Surface_mesh::Halfedge_index, K::FT>("h:p_cost");
	assert(has_placement_costs);
}

boost::optional<SMS::Edge_profile<Surface_mesh>::FT> Custom_cost::operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
	typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

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

			Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
			bool has_point_in_face;
			boost::tie(point_in_face, has_point_in_face) = profile.surface_mesh().property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points");
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
				if (he != Surface_mesh::null_halfedge()) {
					Point_3 A = get(profile.vertex_point_map(), profile.link()[face_id]);
					Point_3 B = get(profile.vertex_point_map(), profile.link()[(face_id + 1) % profile.link().size()]);
					new_faces.push_back(K::Triangle_3(A, B, C));
					new_faces_border_halfedge.push_back(he);
				}
			}

			// geometric error
			std::vector<std::list<Point_set::Index>> points_in_new_face (new_faces.size());
			Triangle_tree tree(new_faces.begin(), new_faces.end());
			for(auto ph: points_to_be_change) {
				auto point = type_converter(point_cloud.point(ph));
				auto point_and_id = tree.closest_point_and_primitive(point);

				points_in_new_face[point_and_id.second - new_faces.begin()].push_back(ph);
				if (alpha > 0) squared_distance += CGAL::squared_distance(point, point_and_id.first);
			}

			// semantic error
			if (beta > 0 || gamma > 0) {
				std::vector<int> new_face_label(new_faces.size(), LABEL_OTHER);
				Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
				bool has_mesh_label;
				boost::tie(mesh_label, has_mesh_label) = profile.surface_mesh().property_map<Surface_mesh::Face_index, unsigned char>("label");
				assert(has_mesh_label);

				std::set<std::size_t> faces_with_no_label;
				for(std::size_t face_id = 0; face_id < new_faces.size(); face_id++) {
					K::FT min_point = 0;
					if (min_surface > 0) {
						min_point = mean_point_per_area * CGAL::sqrt(new_faces[face_id].squared_area()) / min_surface;
					}	
					if (points_in_new_face[face_id].size() > 0) {

						if (points_in_new_face[face_id].size() > min_point) {

							int face_label[LABELS.size()] = {0};

							for (auto point: points_in_new_face[face_id]) {
								face_label[point_cloud_label[point]]++;
							}

							auto argmax = std::max_element(face_label, face_label+LABELS.size());
							count_semantic_error += points_in_new_face[face_id].size() - *argmax;
							new_face_label[face_id] = argmax - face_label;

						} else { // We don't have information about this face.
							new_face_label[face_id] = LABEL_OTHER;
						}
					} else {
						if (min_point <= 1) { // It's probable that there is no point
							faces_with_no_label.insert(face_id);
						} else { // We don't have information about this face.
							new_face_label[face_id] = LABEL_OTHER;
						}
					}
				}
				for(auto face_id: faces_with_no_label) {
					K::FT face_label[LABELS.size()] = {0};
					// Outside face
					auto border_he = profile.surface_mesh().opposite(new_faces_border_halfedge[face_id]);
					if (!profile.surface_mesh().is_border(border_he)) {
						face_label[mesh_label[profile.surface_mesh().face(border_he)]] += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(0), new_faces[face_id].vertex(1)));
					}
					if (profile.surface_mesh().target(new_faces_border_halfedge[(face_id - 1 + new_faces_border_halfedge.size()) % new_faces_border_halfedge.size()]) == profile.surface_mesh().source(new_faces_border_halfedge[face_id]) && faces_with_no_label.count((face_id - 1 + new_faces.size()) % new_faces.size()) == 0) {
						face_label[new_face_label[(face_id - 1 + new_faces.size()) % new_faces.size()]] += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(0), C));
					}
					if (profile.surface_mesh().source(new_faces_border_halfedge[(face_id + 1) % new_faces_border_halfedge.size()]) == profile.surface_mesh().target(new_faces_border_halfedge[face_id]) && faces_with_no_label.count((face_id + 1) % new_faces.size()) == 0) {
						face_label[new_face_label[(face_id + 1) % new_faces.size()]] += CGAL::sqrt(CGAL::squared_distance(new_faces[face_id].vertex(1), C));
					}
					auto argmax = std::max_element(face_label, face_label+LABELS.size());
					if (*argmax > 0) {
						new_face_label[face_id] = argmax - face_label;
					} else { // No neighbor with known label
						new_face_label[face_id] = LABEL_OTHER;
					}
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
			}
		}

		return result_type(- old_cost + alpha * squared_distance + beta * count_semantic_error + gamma * semantic_border_length + delta * placement_costs[profile.v0_v1()]);
	}

	return result_type();
}

Cost_stop_predicate::Cost_stop_predicate(const float cost) : cost(cost) {}

bool Cost_stop_predicate::operator()(const SMS::Edge_profile<Surface_mesh>::FT & current_cost, const SMS::Edge_profile<Surface_mesh> &, const SMS::Edge_profile<Surface_mesh>::edges_size_type, const SMS::Edge_profile<Surface_mesh>::edges_size_type) const {
	return current_cost > cost;
}

My_visitor::My_visitor(const K::FT alpha, const K::FT beta, int min_surface, K::FT mean_point_per_area, Surface_mesh &mesh, const Surface_mesh_info &mesh_info, const Point_set &point_cloud) : alpha(alpha), beta(beta), min_surface(min_surface), mean_point_per_area(mean_point_per_area), mesh(mesh), mesh_info(mesh_info), point_cloud(point_cloud) {
	bool has_point_in_face;
	boost::tie(point_in_face, has_point_in_face) = mesh.property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points");
	assert(has_point_in_face);

	bool has_face_costs;
	boost::tie(face_costs, has_face_costs) = mesh.property_map<Surface_mesh::Face_index, K::FT>("f:cost");
	assert(has_face_costs);
}

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

}

void My_visitor::OnCollapsing (const SMS::Edge_profile<Surface_mesh> &profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>&) {
	// Called when an edge is about to be collapsed and replaced by a vertex whose position is *placement
	for(auto face: profile.triangles()) {
		auto fh = mesh.face(mesh.halfedge(face.v0, face.v1));
		points_to_be_change.insert(point_in_face[fh].begin(), point_in_face[fh].end());
		point_in_face[fh].clear();
		face_costs[fh] = 0;
	}
}

void My_visitor::OnCollapsed (const SMS::Edge_profile<Surface_mesh>&, const Surface_mesh::Vertex_index vd) {
	// Called when an edge has been collapsed and replaced by the vertex vd
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> mesh_label;
	bool has_mesh_label;
	boost::tie(mesh_label, has_mesh_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_mesh_label);

	Point_set::Property_map<unsigned char> point_cloud_label;
	bool has_point_cloud_label;
	boost::tie(point_cloud_label, has_point_cloud_label) = point_cloud.property_map<unsigned char>("p:label");
	assert(has_point_cloud_label);

	// change point_in_face and face_costs
	for(auto ph: points_to_be_change) {
		K::FT min_d = std::numeric_limits<K::FT>::max();
		Surface_mesh::Face_index nearest_face;
		for(auto face: mesh.faces_around_target(mesh.halfedge(vd))) {
			if (face != mesh.null_face()) {
				auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
				auto d = CGAL::squared_distance(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)), type_converter(point_cloud.point(ph)));
				if (d < min_d) {
					nearest_face = face;
					min_d = d;
				}
			}
		}
		assert(nearest_face != mesh.null_face());
		point_in_face[nearest_face].push_back(ph);
		face_costs[nearest_face] += alpha * min_d;
		if (prof.v0_v1().idx() == 144548) std::cerr << ph << " " << type_converter(point_cloud.point(ph)) << " (alpha x p_distance): " << alpha * min_d << "\n";
	}

	std::set<Surface_mesh::Face_index> faces_with_no_label;
	for(auto face: mesh.faces_around_target(mesh.halfedge(vd))) {
		if (face != mesh.null_face()) {
			K::FT min_point = 0;
			if (min_surface > 0 && mean_point_per_area != 0) {
				CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
				boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
				auto p0 = mesh.point(*(vbegin++));
				auto p1 = mesh.point(*(vbegin++));
				auto p2 = mesh.point(*(vbegin++));
				K::FT area = CGAL::sqrt(K::Triangle_3(p0, p1, p2).squared_area());
				min_point = mean_point_per_area * area / min_surface;
			}	
			if (point_in_face[face].size() > 0) {

				if (point_in_face[face].size() > min_point) {

					int face_label[LABELS.size()] = {0};

					for (auto point: point_in_face[face]) {
						face_label[point_cloud_label[point]]++;
					}

					auto argmax = std::max_element(face_label, face_label+LABELS.size());
					mesh_label[face] = argmax - face_label;
					face_costs[face] += beta * (point_in_face[face].size() - *argmax);

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
	}
	for(auto face: faces_with_no_label) {
		K::FT face_label[LABELS.size()] = {0};
		for (auto he: mesh.halfedges_around_face(mesh.halfedge(face))) {
			if (!mesh.is_border(Surface_mesh::Edge_index(he))) {
				if (faces_with_no_label.count(mesh.face(mesh.opposite(he))) == 0) {
					face_label[mesh_label[mesh.face(mesh.opposite(he))]] += CGAL::sqrt(CGAL::squared_distance(mesh.point(mesh.source(he)), mesh.point(mesh.target(he))));
				}
			}
		}
		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		if (*argmax > 0) {
			mesh_label[face] = argmax - face_label;
		} else { // No neighbor with known label
			mesh_label[face] = LABEL_OTHER;
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

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	bool created_point_in_face;
	boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points", std::vector<Point_set::Index>());
	assert(created_point_in_face);

	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
	bool created_face_costs;
	boost::tie(face_costs, created_face_costs) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:cost", 0);
	assert(created_face_costs);

	Point_set point_cloud;
	bool created_point_label;
	Point_set::Property_map<unsigned char> point_cloud_label;
	boost::tie (point_cloud_label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", LABEL_OTHER);
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
			point_cloud_label[point_index[L][P]] = raster.land_cover[L][P];
		}
	}
	std::cout << "Point added" << std::endl;

	// Set isBorder
	int N = 10;
	Point_tree tree(point_cloud.begin(), point_cloud.end(), Point_tree::Splitter(), TreeTraits(point_cloud.point_map()));
	Neighbor_search::Distance tr_dist(point_cloud.point_map());
	for (auto point: point_cloud) {
		Neighbor_search search(tree, point_cloud.point(point), N, 0, true, tr_dist, false);
		unsigned char l = point_cloud_label[search.begin()->first];
		for(auto it = search.begin() + 1; it != search.end() && !isborder[point]; ++it) {
			if (point_cloud_label[it->first] != l) isborder[point] = true;
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
				face_costs[f1] = 0;
				auto f2 = mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P+1]);
				point_in_face[f2].push_back(point_index[L][P]);
				point_in_face[f2].push_back(point_index[L][P+1]);
				point_in_face[f2].push_back(point_index[L+1][P+1]);
				face_costs[f2] = 0;
			} else {
				auto f1 = mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P]);
				point_in_face[f1].push_back(point_index[L][P]);
				point_in_face[f1].push_back(point_index[L][P+1]);
				point_in_face[f1].push_back(point_index[L+1][P]);
				face_costs[f1] = 0;
				auto f2 = mesh.add_face(vertex_index[L][P+1], vertex_index[L+1][P+1], vertex_index[L+1][P]);
				point_in_face[f2].push_back(point_index[L][P+1]);
				point_in_face[f2].push_back(point_index[L+1][P+1]);
				point_in_face[f2].push_back(point_index[L+1][P]);
				face_costs[f2] = 0;
			}
		}
	}
	std::cout << "Faces added" << std::endl;

	float alpha = 1, beta = 1;
	for(auto face: mesh.faces()) {
		if (point_in_face[face].size() > 0) {
			int face_label[LABELS.size()] = {0};
			int sum_face_label = 0;
			for (auto ph: point_in_face[face]) {
				sum_face_label++;
				face_label[point_cloud_label[ph]]++;
			}
			auto argmax = std::max_element(face_label, face_label+LABELS.size());
			face_costs[face] = beta * (sum_face_label - *argmax);
		}
	}
	std::cout << "Faces cost" << std::endl;

	// Return mesh if coords are in reverse order
	if ((x_1-x_0)*(y_1-y_0) < 0) {
		CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh); 	
	}

	auto created_label = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label", LABEL_OTHER);
	assert(created_label.second);

	int min_surface = 2;
	K::FT mean_point_per_area = add_label(mesh, point_cloud, min_surface);
	mesh_info.save_mesh(mesh, "initial-mesh.ply");

	Surface_mesh::Property_map<Surface_mesh::Halfedge_index, K::FT> placement_costs;
	bool created_placement_costs;
	boost::tie(placement_costs, created_placement_costs) = mesh.add_property_map<Surface_mesh::Halfedge_index, K::FT>("h:p_cost", 0);
	assert(created_placement_costs);

	Cost_stop_predicate stop(5);
	//SMS::Count_stop_predicate<Surface_mesh> stop(50);
	const LindstromTurk_param params (1,1,1,1,0.1,1,1);
	Custom_placement pf(params, mesh, point_cloud);
	Custom_cost cf(alpha, beta, 1, 1, min_surface, mean_point_per_area, mesh, point_cloud);
	My_visitor mv (alpha, beta, min_surface, mean_point_per_area, mesh, mesh_info, point_cloud);
	SMS::Bounded_normal_change_filter<> filter;
	SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(mv));
	std::cout << "\rMesh simplified                                               " << std::endl;

	add_label(mesh, point_cloud, min_surface, mean_point_per_area);
	mesh_info.save_mesh(mesh, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, mesh);
}
