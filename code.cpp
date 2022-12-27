#include "header.hpp"
#include "raster.hpp"

std::list<Polygon> get_LOD0_from_shapefile(char *path) {
	GDALAllRegister();

	GDALDataset *dataset;

	std::list<Polygon> polygons;

	dataset = (GDALDataset*) GDALOpenEx(path, GDAL_OF_VECTOR, NULL, NULL, NULL );
	if( dataset == NULL ) {
		std::cerr << "Unable to open " << path << "." << std::endl;
		return polygons;
	}

	for( OGRLayer *layer: dataset->GetLayers() ) {
		for( const auto& feature: *layer ) {
			if ((*feature)["num_class"].GetInteger() == 4) {
				// Get only the buildings
				OGRGeometry *geometry = feature->GetGeometryRef();

				if (wkbFlatten(geometry->getGeometryType()) == wkbPolygon) {

					OGRPolygon *shp_polygon = geometry->toPolygon();
					Polygon cgal_polygon;
					std::istringstream wkt(shp_polygon->exportToWkt());
					CGAL::IO::read_polygon_WKT (wkt, cgal_polygon);
					polygons.push_back(cgal_polygon);

				} else if (wkbFlatten(geometry->getGeometryType()) == wkbMultiPolygon) {

					OGRMultiPolygon *shp_multi_polygon = geometry->toMultiPolygon();
					for (OGRPolygon *shp_polygon: *shp_multi_polygon) {
						Polygon cgal_polygon;
						std::istringstream wkt(shp_polygon->exportToWkt());
						CGAL::IO::read_polygon_WKT (wkt, cgal_polygon);
						polygons.push_back(cgal_polygon);
					}

				}
			}
		}
	}

	return polygons;
}

float single_face_cost(const Raster &raster, const Point_3 &p0, const Point_3 &p1, const Point_3 &p2) {
	float nz = ((-p0.x() + p1.x()) * (-p0.y() + p2.y()) - (-p0.x() + p2.x()) * (-p0.y() + p1.y()));
	if (nz == 0) {
		// flat triangle
		return 0;
	}

	std::list<std::pair<int,int>> pixels = raster.triangle_to_pixel(p0, p1, p2);

	// Entropy
	int face_label[LABELS.size()] = {0};
	int sum_face_label = 0;
	for (auto pixel : pixels) {
		if (raster.land_cover[pixel.second][pixel.first] > -1) {
			sum_face_label++;
			face_label[raster.land_cover[pixel.second][pixel.first]]++;
		}
	}
	float entropy = 0;
	if (sum_face_label > 0) {
		for (int i = 0; i < LABELS.size(); i++) {
			if (face_label[i] > 0) {
				entropy += ((float) face_label[i])*log((float) face_label[i]);
			}
		}
		entropy = log((float) sum_face_label) - entropy/((float) sum_face_label);
	}

	// Least squares
	float least_squares = 0;
	float least = 0;
	if (pixels.size() != 0) {
		K::Plane_3 plane(p0, p1, p2);
		for (auto pixel : pixels) {
			float px = 0.5 + pixel.first;
			float py = 0.5 + pixel.second;
			float pz = - (plane.a() * px + plane.b() * py + plane.d()) / plane.c();
			least_squares += pow(raster.dsm[pixel.second][pixel.first] - pz,2);
			least += abs(raster.dsm[pixel.second][pixel.first] - pz);
		}
		least_squares /= pixels.size();
		least /= pixels.size();
	}

	// Verticality
	float verticality = 0;
	if (pixels.size() == 0) {
		float surface = pow(K::Triangle_3(p0, p1, p2).squared_area (), 0.5);
		verticality = abs(nz)/(2*surface);
	}

	// Eccentricity
	Point_3 S = CGAL::centroid(p0,p1,p2);
	float M = (pow(p2.x()-S.x(),2) + pow(p2.y()-S.y(),2) + pow(p2.z()-S.z(),2) + (pow(p1.x()-p0.x(),2) + pow(p1.y()-p0.y(),2) + pow(p1.z()-p0.z(),2))/3)/4;
	float N = sqrtf(pow((p2.y()-S.y())*(p1.z()-p0.z())-(p2.z()-S.z())*(p1.y()-p0.y()),2) + pow((p2.z()-S.z())*(p1.x()-p0.x())-(p2.x()-S.x())*(p1.z()-p0.z()),2) + pow((p2.x()-S.x())*(p1.y()-p0.y())-(p2.y()-S.y())*(p1.x()-p0.x()),2))/(4 * sqrtf(3));
	float eccentricity = M*M - 4*N*N;
	eccentricity = sqrtf(1-(M-eccentricity)/(M+eccentricity));

	return 0 * least_squares + 1 * least + 1 * entropy + 0 * verticality + 0 * eccentricity;
}

float face_cost(const Raster &raster, const Point_3 &p0, const Point_3 &p1, const Point_3 &p2) {
	float surface = pow(K::Triangle_3(p0, p1, p2).squared_area (), 0.5);
	return surface * (1 + single_face_cost(raster, p0, p1, p2));
}

Point_3 best_point(const Raster &raster, K::FT x, K::FT y, const SMS::Edge_profile<Surface_mesh>& profile) {

	float z;
	float d = 0;
	float D = 0;
	std::list<std::pair<float,float>> values;
	float t = 0;
	int count = 0;

	//Foreach face
	for (auto he : CGAL::halfedges_around_source(profile.v0(), profile.surface_mesh())) {
		if (he != profile.v0_v1() && he != profile.v0_vR()) {
			Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
			Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
			float i2 = ((-A.x() + B.x()) * (-A.y() + y) - (-A.x() + x) * (-A.y() + B.y()));
			if (i2 != 0) {
				for (auto pixel : raster.triangle_to_pixel(A, B, Point_3(x, y, 0))) { //for each pixel
					float px = 0.5 + pixel.first;
					float py = 0.5 + pixel.second;
					float pz = raster.dsm[pixel.second][pixel.first];

					float i1 = ((A.x() - B.x()) * (-A.y() + py) + (-A.x() + px) * (-A.y() + B.y()));

					d += i1 * (((A.x() - px) * (A.y() * B.z() - A.z() * B.y() + A.z() * y - B.z() * y) + (A.y() - py) * (-A.x() * B.z() + A.z() * B.x() - A.z() * x + B.z() * x) + (A.z() - pz) * i2) / (i2 * i2));
					D += (i1 * i1) / (i2 * i2);
					values.push_back(std::pair<float,float>(((A.x() - px) * (A.y() * B.z() - A.z() * B.y() + A.z() * y - B.z() * y) + (A.y() - py) * (-A.x() * B.z() + A.z() * B.x() - A.z() * x + B.z() * x) + (A.z() - pz) * i2) / i1, abs(i1/i2)));
					t += abs(i1/i2);
					count ++;
				}
			}
		}
	}
	for (auto he : CGAL::halfedges_around_source(profile.v1(), profile.surface_mesh())) {
		if (he != profile.v1_v0() && he != profile.v1_vL()) {
			Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
			Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
			float i2 = ((-A.x() + B.x()) * (-A.y() + y) - (-A.x() + x) * (-A.y() + B.y()));
			if (i2 != 0) {
				for (auto pixel : raster.triangle_to_pixel(A, B, Point_3(x, y, 0))) { //for each pixel
					float px = 0.5 + pixel.first;
					float py = 0.5 + pixel.second;
					float pz = raster.dsm[pixel.second][pixel.first];

					float i1 = ((A.x() - B.x()) * (-A.y() + py) + (-A.x() + px) * (-A.y() + B.y()));

					d += i1 * (((A.x() - px) * (A.y() * B.z() - A.z() * B.y() + A.z() * y - B.z() * y) + (A.y() - py) * (-A.x() * B.z() + A.z() * B.x() - A.z() * x + B.z() * x) + (A.z() - pz) * i2) / (i2 * i2));
					D += (i1 * i1) / (i2 * i2);
					values.push_back(std::pair<float,float>(((A.x() - px) * (A.y() * B.z() - A.z() * B.y() + A.z() * y - B.z() * y) + (A.y() - py) * (-A.x() * B.z() + A.z() * B.x() - A.z() * x + B.z() * x) + (A.z() - pz) * i2) / i1, abs(i1/i2)));
					t += abs(i1/i2);
					count ++;
				}
			}
		}
	}

	if (count > 0) {
		z = d / D;

		values.sort([](std::pair<float,float> a, std::pair<float,float> b) {
			return a.first > b.first;
		});

		auto p = values.begin();
		float s = 0;
		while (s+p->second < t/2) {
			s += p->second;
			p++;
		}
		z = p->first;

	} else {
		if (profile.p1().x() != profile.p0().x()) {
			z = (x - profile.p0().x())/(profile.p1().x() - profile.p0().x())*(profile.p1().z() - profile.p0().z());
		} else {
			z = (y - profile.p0().y())/(profile.p1().y() - profile.p0().y())*(profile.p1().z() - profile.p0().z());
		}
	}

	return Point_3(x,y,z);
}

class Custom_placement {
	const Raster &raster;

	public:
		Custom_placement (const Raster &raster) : raster(raster) {}

		boost::optional<SMS::Edge_profile<Surface_mesh>::Point> operator()(const SMS::Edge_profile<Surface_mesh>& profile) const {
			typedef boost::optional<SMS::Edge_profile<Surface_mesh>::Point> result_type;

			K::Vector_3 vector(profile.p0(), profile.p1());
			Point_3 p[5] = {profile.p0(), profile.p0() + 0.25 * vector, profile.p0() + 0.5 * vector, profile.p0() + 0.75 * vector, profile.p1()};
			float cost[5] = {0};

			/*for (int j = 0; j < 5; j++) {
				p[j] = best_point(raster, p[j].x(), p[j].y(), profile);
			}*/

			for (auto he : CGAL::halfedges_around_source(profile.v0(), profile.surface_mesh())) {
				if (he != profile.v0_v1() && he != profile.v0_vR()) {
					Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
					Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
					for (int j = 0; j < 5; j++) {
						cost[j] += face_cost(raster, A, B, p[j]);
					}
				}
			}
			for (auto he : CGAL::halfedges_around_source(profile.v1(), profile.surface_mesh())) {
				if (he != profile.v1_v0() && he != profile.v1_vL()) {
					Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
					Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
					for (int j = 0; j < 5; j++) {
						cost[j] += face_cost(raster, A, B, p[j]);
					}
				}
			}

			for (int i = 0; i < 2; i++) {
				int min_cost = std::min_element(cost, cost + 5) - cost;

				if (min_cost == 0 || min_cost == 1) {
					p[4] = p[2];
					cost[4] = cost[2];
					p[2] = p[1];
					cost[2] = cost[1];
				} else if (min_cost == 2) {
					p[0] = p[1];
					cost[0] = cost[1];
					p[4] = p[3];
					cost[4] = cost[3];
				} else {
					p[0] = p[2];
					cost[0] = cost[2];
					p[2] = p[3];
					cost[2] = cost[3];
				}

				vector = K::Vector_3(p[0], p[4]);
				p[1] = p[0] + 0.25 * vector;
				p[3] = p[0] + 0.75 * vector;
				//p[1] = best_point(raster, p[1].x(), p[1].y(), profile);
				//p[3] = best_point(raster, p[1].x(), p[1].y(), profile);

				cost[1] = 0;
				cost[3] = 0;

				for (auto he : CGAL::halfedges_around_source(profile.v0(), profile.surface_mesh())) {
					if (he != profile.v0_v1() && he != profile.v0_vR()) {
						Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
						Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
						cost[1] += face_cost(raster, A, B, p[1]);
						cost[3] += face_cost(raster, A, B, p[3]);
					}
				}
				for (auto he : CGAL::halfedges_around_source(profile.v1(), profile.surface_mesh())) {
					if (he != profile.v1_v0() && he != profile.v1_vL()) {
						Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
						Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
						cost[1] += face_cost(raster, A, B, p[1]);
						cost[3] += face_cost(raster, A, B, p[3]);
					}
				}
			}

			int min_cost = std::min_element(cost, cost + 5) - cost;

			//std::cout << "Placement: (" << profile.p0() << ") - (" << profile.p1() << ") -> (" << p[min_cost] << ")\n";
			Point_3 placement = p[min_cost];

			for (auto he : CGAL::halfedges_around_source(profile.v0(), profile.surface_mesh())) {
				if (he != profile.v0_v1() && he != profile.v0_vR()) {
					Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
					Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
					if (CGAL::orientation(Point_2(A.x(),A.y()), Point_2(B.x(),B.y()), Point_2(profile.p0().x(),profile.p0().y())) != CGAL::orientation(Point_2(A.x(),A.y()), Point_2(B.x(),B.y()), Point_2(placement.x(),placement.y()))) {
						return result_type();
					}
				}
			}
			for (auto he : CGAL::halfedges_around_source(profile.v1(), profile.surface_mesh())) {
				if (he != profile.v1_v0() && he != profile.v1_vL()) {
					Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
					Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
					if (CGAL::orientation(Point_2(A.x(),A.y()), Point_2(B.x(),B.y()), Point_2(profile.p1().x(),profile.p1().y())) != CGAL::orientation(Point_2(A.x(),A.y()), Point_2(B.x(),B.y()), Point_2(placement.x(),placement.y()))) {
						return result_type();
					}
				}
			}

			return result_type(placement);
		}
};

class Custom_cost {
	const Raster &raster;

	public:
		Custom_cost (const Raster &raster) : raster(raster) {}

		boost::optional<SMS::Edge_profile<Surface_mesh>::FT> operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
			typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;

			if (placement) {

				float old_cost = 0;
				float new_cost = 0;

				SMS::Edge_profile<Surface_mesh>::Triangle_vector triangles = profile.triangles();
				for (auto triange: triangles) {
					Point_3 A = get(profile.vertex_point_map(),triange.v0);
					Point_3 B = get(profile.vertex_point_map(),triange.v1);
					Point_3 C = get(profile.vertex_point_map(),triange.v2);
					old_cost += face_cost(raster, A, B, C);
				}

				Point_3 C = *placement;
				for (auto he : CGAL::halfedges_around_source(profile.v0(), profile.surface_mesh())) {
					if (he != profile.v0_v1() && he != profile.v0_vR()) {
						Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
						Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
						new_cost += face_cost(raster, A, B, C);
					}
				}
				for (auto he : CGAL::halfedges_around_source(profile.v1(), profile.surface_mesh())) {
					if (he != profile.v1_v0() && he != profile.v1_vL()) {
						Point_3 A = get(profile.vertex_point_map(),CGAL::target(he, profile.surface_mesh()));
						Point_3 B = get(profile.vertex_point_map(),CGAL::target(CGAL::next(he, profile.surface_mesh()), profile.surface_mesh()));
						new_cost += face_cost(raster, A, B, C);
					}
				}

				//std::cout << "Cost: " << (new_cost - old_cost) << "\n";

				return result_type(new_cost - old_cost);
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

void save_mesh(const Surface_mesh &mesh, const Raster &raster, const char *filename) {
	typedef CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3> OutputMesh;
	OutputMesh output_mesh;

	std::unordered_map<boost::graph_traits<Surface_mesh>::face_descriptor, boost::graph_traits<OutputMesh>::face_descriptor> f2f;
	CGAL::copy_face_graph (mesh, output_mesh, CGAL::parameters::face_to_face_output_iterator(std::inserter(f2f, f2f.end())));

	// Label
	OutputMesh::Property_map<OutputMesh::Face_index, unsigned char> output_label;
	bool created;
	boost::tie(output_label, created) = output_mesh.add_property_map<OutputMesh::Face_index, unsigned char>("label",0);
	assert(created);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	if (has_label) {
		for (auto face : mesh.faces()) {
			output_label[f2f[face]] = label[face];
		}
	}

	// Color
	OutputMesh::Property_map<OutputMesh::Face_index, unsigned char> red;
	OutputMesh::Property_map<OutputMesh::Face_index, unsigned char> green;
	OutputMesh::Property_map<OutputMesh::Face_index, unsigned char> blue;
	boost::tie(red, created) = output_mesh.add_property_map<OutputMesh::Face_index, unsigned char>("red",0);
	assert(created);
	boost::tie(green, created) = output_mesh.add_property_map<OutputMesh::Face_index, unsigned char>("green",0);
	assert(created);
	boost::tie(blue, created) = output_mesh.add_property_map<OutputMesh::Face_index, unsigned char>("blue",0);
	assert(created);

	// Entropy
	OutputMesh::Property_map<OutputMesh::Face_index, float> quality;
	boost::tie(quality, created) = output_mesh.add_property_map<OutputMesh::Face_index, float>("quality",0);
	assert(created);

	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	if (has_path) {
		for (auto face : mesh.faces()) {
			quality[f2f[face]] = path[face];
		}
	}

	for (auto face : output_mesh.faces()) {
		if (!has_label) {
			int face_label[LABELS.size()] = {0};

			CGAL::Vertex_around_face_iterator<OutputMesh> vbegin, vend;
			boost::tie(vbegin, vend) = vertices_around_face(output_mesh.halfedge(face), output_mesh);
			for (auto pixel : raster.triangle_to_pixel(output_mesh.point(*(vbegin++)), output_mesh.point(*(vbegin++)), output_mesh.point(*(vbegin++)))) {
				if (raster.land_cover[pixel.second][pixel.first] > -1) {
					face_label[raster.land_cover[pixel.second][pixel.first]]++;
				}
			}

			auto argmax = std::max_element(face_label, face_label+LABELS.size());
			output_label[face] = argmax - face_label;
		}

		red[face] = LABELS[output_label[face]].red;
		green[face] = LABELS[output_label[face]].green;
		blue[face] = LABELS[output_label[face]].blue;

		if (!has_path) {
			CGAL::Vertex_around_face_iterator<OutputMesh> vbegin, vend;
			boost::tie(vbegin, vend) = vertices_around_face(output_mesh.halfedge(face), output_mesh);
			auto pa = output_mesh.point(*(vbegin++));
			auto pb = output_mesh.point(*(vbegin++));
			auto pc = output_mesh.point(*(vbegin++));
			quality[face] = single_face_cost(raster, Point_3(pa.x(), pa.y(), pa.z()), Point_3(pb.x(), pb.y(), pb.z()), Point_3(pc.x(), pc.y(), pc.z()));
		}
	}

	double min_x, min_y;
	raster.grid_to_coord(0, 0, min_x, min_y);

	char *temp;
	const char *options_json[] = { "MULTILINE=NO", NULL };
	raster.get_crs().exportToPROJJSON(&temp, options_json);
	std::string crs_as_json(temp);
	CPLFree(temp);

	crs_as_json = regex_replace(crs_as_json, std::regex(",\"id\":\\{\"authority\":\"EPSG\",\"code\":[^\\}]+\\}"), "");
	crs_as_json = regex_replace(crs_as_json, std::regex("Lambert-93"), "Custom");

	std::smatch matches;
	regex_search(crs_as_json, matches, std::regex("\"name\":\"Northing at false origin\",\"value\":([^,]+),"));
	crs_as_json = regex_replace(crs_as_json, std::regex("\"name\":\"Northing at false origin\",\"value\":([^,]+),"), "\"name\":\"Northing at false origin\",\"value\":" + std::to_string(std::stoi(matches[1]) - min_y) + ",");

	regex_search(crs_as_json, matches, std::regex("\"name\":\"Easting at false origin\",\"value\":([^,]+),"));
	crs_as_json = regex_replace(crs_as_json, std::regex("\"name\":\"Easting at false origin\",\"value\":([^,]+),"), "\"name\":\"Easting at false origin\",\"value\":" + std::to_string(std::stoi(matches[1]) - min_x) + ",");

	OGRSpatialReference output_crs;
	output_crs.SetFromUserInput(crs_as_json.c_str());
	output_crs.AutoIdentifyEPSG();

	/*const char *options_wkt[] = { "MULTILINE=NO", "FORMAT=WKT2", NULL };
	output_crs.exportToWkt(&temp, options_wkt);*/ output_crs.exportToProj4(&temp); // WKT format is too long for MeshLab
	std::string crs_as_wkt(temp);
	CPLFree(temp);

	for(auto vertex : output_mesh.vertices()) {
		auto point = output_mesh.point(vertex);
		double x, y;
		raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
		output_mesh.point(vertex) = CGAL::Simple_cartesian<double>::Point_3(x-min_x, y-min_y, (double) point.z());
	}

	std::ofstream mesh_ofile (filename, std::ios_base::binary);
	CGAL::IO::set_binary_mode (mesh_ofile);
	CGAL::IO::write_PLY (mesh_ofile, output_mesh, "crs " + crs_as_wkt);
	mesh_ofile.close();
}


struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh> {
	private:
		int i_collecte = 0;
		const Surface_mesh *mesh;
		const Raster *raster;
		std::chrono::time_point<std::chrono::system_clock> start_collecte;
		std::chrono::time_point<std::chrono::system_clock> start_collapse;
		bool output[5] = {false};

	public:
		My_visitor(const Surface_mesh *mesh, const Raster *raster) : mesh(mesh), raster(raster) {}

		void OnStarted (Surface_mesh &mesh) {
			start_collecte = std::chrono::system_clock::now();
		}

		void OnCollected(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional< SMS::Edge_profile<Surface_mesh>::FT >& cost) {
			start_collapse = std::chrono::system_clock::now();
			i_collecte++;
			if (i_collecte%1000 == 0) {
				std::chrono::duration<double> diff = start_collapse - start_collecte;
				std::cout << "\rCollecte: " << i_collecte << "/" << mesh->number_of_edges() << " (" << ((int) (((float) i_collecte)/mesh->number_of_edges()*100)) << "%)" << " still " << (((float) mesh->number_of_edges() - i_collecte) * diff.count() / i_collecte) << "s" << " (" << (((float) i_collecte) / diff.count()) << " op/s)" << std::flush;
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
				if(*cost > 1e-4 && !output[0]) {
					output[0] = true;
					save_mesh(*mesh,*raster,"mesh-1e-4.ply");
				} else if(*cost > 0 && !output[1]) {
					output[1] = true;
					save_mesh(*mesh,*raster,"mesh-0.ply");
				}
				if(current_edge_count <= 100000 && !output[2]) {
					output[2] = true;
					save_mesh(*mesh,*raster,"mesh-100000.ply");
				} else if(current_edge_count <= 10000 && !output[3]) {
					output[3] = true;
					save_mesh(*mesh,*raster,"mesh-10000.ply");
				} else if(current_edge_count <= 5000 && !output[4]) {
					output[4] = true;
					save_mesh(*mesh,*raster,"mesh-5000.ply");
				} else if(current_edge_count <= 1000000 && !output[5]) {
					output[5] = true;
					save_mesh(*mesh,*raster,"mesh-1000000.ply");
				}
			}

		}

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
	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			vertex_index[L][P] = mesh.add_vertex(Point_3(0.5 + P, 0.5 + L, raster.dsm[L][P]));
		}
	}
	std::cout << "Point added" << std::endl;
	// Add faces
	for (int L = 0; L < raster.ySize-1; L++) {
		for (int P = 0; P < raster.xSize-1; P++) {
			if (pow(raster.dsm[L][P]-raster.dsm[L+1][P+1], 2) < pow(raster.dsm[L+1][P]-raster.dsm[L][P+1], 2)) {
				mesh.add_face(vertex_index[L][P], vertex_index[L+1][P+1], vertex_index[L+1][P]);
				mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P+1]);
			} else {
				mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P]);
				mesh.add_face(vertex_index[L][P+1], vertex_index[L+1][P+1], vertex_index[L+1][P]);
			}
		}
	}
	std::cout << "Faces added" << std::endl;

	save_mesh(mesh, raster, "initial-mesh.ply");

	Cost_stop_predicate stop(10);
	//SMS::Count_stop_predicate<Surface_mesh> stop(1000);
	Custom_cost cf(raster);
	Custom_placement pf(raster);
	int r = SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).get_placement(pf).visitor(My_visitor(&mesh, &raster)));
	std::cout << "\rMesh simplified                                               " << std::endl;

	save_mesh(mesh, raster, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, mesh);
}

void add_label(const Raster &raster, Surface_mesh &mesh) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool created;
	boost::tie(label, created) = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label",0);
	assert(created);

	for (auto face : mesh.faces()) {
		int face_label[LABELS.size()] = {0};
		int sum_face_label = 0;

		CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
		boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
		for (auto pixel : raster.triangle_to_pixel(mesh.point(*(vbegin++)), mesh.point(*(vbegin++)), mesh.point(*(vbegin++)))) {
			if (raster.land_cover[pixel.second][pixel.first] > -1) {
				sum_face_label++;
				face_label[raster.land_cover[pixel.second][pixel.first]]++;
			}
		}

		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		label[face] = argmax - face_label;
	}
}

void change_vertical_faces(Surface_mesh &mesh, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::unordered_map<Surface_mesh::Face_index, char> new_label;
	std::list<Surface_mesh::Face_index> remove_face;

	for (auto face : mesh.faces()) {
		new_label[face] = label[face];
		CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
		boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
		auto p0 = mesh.point(*(vbegin++));
		auto p1 = mesh.point(*(vbegin++));
		auto p2 = mesh.point(*(vbegin++));
		//TODO : use coords and not grid
		float nz = ((-p0.x() + p1.x()) * (-p0.y() + p2.y()) - (-p0.x() + p2.x()) * (-p0.y() + p1.y()));
		float surface = pow(K::Triangle_3(p0, p1, p2).squared_area(), 0.5);
		if (abs(nz)/(2*surface) < 0.5) {
			/*CGAL::Face_around_face_iterator<Surface_mesh> fbegin, fend;
			for(boost::tie(fbegin, fend) = faces_around_face(mesh.halfedge(face), mesh); fbegin != fend && new_label[face] != 4; ++fbegin) {
				if (*fbegin != boost::graph_traits<Surface_mesh>::null_face()) {
					if (label[*fbegin] == 4) {
						//Building
						new_label[face] = 4;
					} else if (label[*fbegin] == 5 && new_label[face] != 4) {
						//High vegetation
						new_label[face] = 5;
					}
					CGAL::Face_around_face_iterator<Surface_mesh> fbegin2, fend2;
					for(boost::tie(fbegin2, fend2) = faces_around_face(mesh.halfedge(*fbegin), mesh); fbegin2 != fend2 && new_label[face] != 4; ++fbegin2) {
						if (*fbegin2 != boost::graph_traits<Surface_mesh>::null_face()) {
							if (label[*fbegin2] == 4) {
								//Building
								new_label[face] = 4;
							} else if (label[*fbegin2] == 5 && new_label[face] != 4) {
								//High vegetation
								new_label[face] = 5;
							}
						}
					}
				}
			}*/

			int face_label[LABELS.size()] = {0};
			int sum_face_label = 0;

			CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
			boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
			for (auto pixel : raster.triangle_to_pixel(mesh.point(*(vbegin++)), mesh.point(*(vbegin++)), mesh.point(*(vbegin++)))) {
				if (raster.land_cover[pixel.second][pixel.first] > -1) {
					sum_face_label++;
					face_label[raster.land_cover[pixel.second][pixel.first]]++;
				}
			}

			if (face_label[4] > 0) {
				//Building
				new_label[face] = 4;
			} else if (face_label[5] > 0) {
				//High vegetation
				new_label[face] = 5;
			}

			if (new_label[face] != 4 && new_label[face] != 5) {
				new_label[face] = 0;
			}
		}
	}

	for (auto &face_label: new_label) {
		label[face_label.first] = face_label.second;
	}

}

void set_path(Surface_mesh &mesh,
				Surface_mesh::Property_map<Surface_mesh::Face_index, int> &path,
				Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> &label,
				std::vector<std::list<Surface_mesh::Face_index>> &paths,
				Surface_mesh::Face_index face,
				int path_id) {
	path[face] = path_id;
	paths[path_id].push_back(face);
	CGAL::Face_around_face_iterator<Surface_mesh> fbegin, fend;
	for(boost::tie(fbegin, fend) = faces_around_face(mesh.halfedge(face), mesh); fbegin != fend; ++fbegin) {
		if (*fbegin != boost::graph_traits<Surface_mesh>::null_face()) {
			if (label[face] == label[*fbegin] && path[*fbegin] == -1) {
				set_path(mesh, path, label, paths, *fbegin, path_id);
			}
		}
	}
}

std::vector<std::list<Surface_mesh::Face_index>> compute_path(Surface_mesh &mesh) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool created;
	boost::tie(path, created) = mesh.add_property_map<Surface_mesh::Face_index, int>("path", -1);
	assert(created);

	std::vector<std::list<Surface_mesh::Face_index>> paths;

	for (auto face: mesh.faces()) {
		if (path[face] == -1) {
			paths.push_back(std::list<Surface_mesh::Face_index>{});
			set_path(mesh, path, label, paths, face, paths.size()-1);
		}
	}

	return paths;
}

template<typename Face_handle>
CGAL::Polygon_with_holes_2<typename Face_handle::value_type::Vertex::Point::R> polygon(Face_handle face) {

	typedef typename Face_handle::value_type::Vertex::Point::R Kernel;

	CGAL::Polygon_2<Kernel> border;
	auto curr = face->outer_ccb();
	do {
		auto p = curr->target()->point();
		border.push_back(p);
	} while (++curr != face->outer_ccb());

	std::list<CGAL::Polygon_2<Kernel>> holes;
	for(auto hole = face->holes_begin(); hole != face->holes_end(); hole++) {
		auto curr = *hole;
		CGAL::Polygon_2<Kernel> hole_polygon;
		do {
			auto p = curr->target()->point();
			hole_polygon.push_back(p);
		} while (++curr != *hole);
		holes.push_back(hole_polygon);
	} 

	return CGAL::Polygon_with_holes_2<Kernel>(border, holes.begin(), holes.end());
}


std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> compute_path_polygon(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> path_polygon;

	typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
	for (int i = 0; i < paths.size(); i++) {

		Filtered_graph filtered_sm(mesh, i, path);

		if (filtered_sm.number_of_faces() > 3) {

			int lab = label[*(CGAL::faces(filtered_sm).first)];
			if (lab == 3 || lab == 8 || lab == 9) {

				Arrangement_2 arr;
				std::map<Surface_mesh::vertex_index, Arrangement_2::Vertex_handle> point_map;
				for (auto edge: CGAL::edges(filtered_sm)) {
					if (CGAL::is_border (edge, filtered_sm)) {
						auto p0 = mesh.point(CGAL::source(edge, filtered_sm));
						auto p1 = mesh.point(CGAL::target(edge, filtered_sm));
						insert_non_intersecting_curve(arr, Traits_2::X_monotone_curve_2(Exact_predicates_kernel::Point_2(p0.x(), p0.y()), Exact_predicates_kernel::Point_2(p1.x(), p1.y())));
					}
				}

				path_polygon[i] = polygon((*(arr.unbounded_face()->holes_begin()))->twin()->face());

				{ // Arrangement
					Surface_mesh skeleton;

					std::map<Arrangement_2::Vertex_handle, Surface_mesh::vertex_index> v_map;
					for (auto v = arr.vertices_begin(); v != arr.vertices_end(); v++) {
						auto p = v->point();
						auto z = raster.dsm[int(p.y())][int(p.x())];
						v_map[v] = skeleton.add_vertex(Point_3((float) p.x(), (float) p.y(), z));
					}

					for (auto he = arr.edges_begin(); he != arr.edges_end(); ++he ) {
						skeleton.add_edge(v_map[he->source()], v_map[he->target()]);
					}

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
					skeleton_name << "arr_" << lab << "_" << i << ".ply";
					std::ofstream mesh_ofile (skeleton_name.str().c_str());
					CGAL::IO::write_PLY (mesh_ofile, skeleton);
					mesh_ofile.close();
				}

				Surface_mesh part_mesh;
				CGAL::copy_face_graph(filtered_sm, part_mesh);
				std::stringstream name;
				name << "part_mesh_" << lab << "_" << i << ".ply";
				save_mesh(part_mesh, raster, name.str().c_str());

			}
		}
	}

	return path_polygon;

}


std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> compute_medial_axes(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> medial_axes;

	for (int i = 0; i < paths.size(); i++) {

		if (path_polygon.count(i) > 0) {

			int lab = label[paths.at(i).front()];

			auto poly = path_polygon.at(i);

			poly = CGAL::Polyline_simplification_2::simplify(
				poly,
				CGAL::Polyline_simplification_2::Squared_distance_cost(),
				CGAL::Polyline_simplification_2::Stop_above_cost_threshold(pow(raster.coord_distance_to_grid_distance(1.5),2))
			);

			auto iss = CGAL::create_interior_straight_skeleton_2(poly, Exact_predicates_kernel());
			medial_axes[i] = CGAL::convert_straight_skeleton_2<CGAL::Straight_skeleton_2<K>>(*iss);

			{ // Skeleton
				Surface_mesh skeleton;

				std::map<int, Surface_mesh::vertex_index> v_map;
				for (auto v = iss->vertices_begin(); v != iss->vertices_end(); v++) {
					auto p = v->point();
					auto z = raster.dsm[int(p.y())][int(p.x())];
					v_map[v->id()] = skeleton.add_vertex(Point_3((float) p.x(), (float) p.y(), z));
				}

				for (auto he = iss->halfedges_begin(); he != iss->halfedges_end(); ++he ) {
					skeleton.add_edge(v_map[he->vertex()->id()], v_map[he->opposite()->vertex()->id()]);
				}

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
				skeleton_name << "skeleton_" << lab << "_" << i << ".ply";
				std::ofstream mesh_ofile (skeleton_name.str().c_str());
				CGAL::IO::write_PLY (mesh_ofile, skeleton);
				mesh_ofile.close();
			}

			{ // Path
				Surface_mesh skeleton;

				std::map<int, Surface_mesh::vertex_index> v_map;
				for (auto v = iss->vertices_begin(); v != iss->vertices_end(); v++) {
					if (v->is_skeleton()) {
						auto p = v->point();
						auto z = raster.dsm[int(p.y())][int(p.x())];
						v_map[v->id()] = skeleton.add_vertex(Point_3((float) p.x(), (float) p.y(),z));
					}
				}

				for (auto he = iss->halfedges_begin(); he != iss->halfedges_end(); ++he ) {
					if (he->is_inner_bisector()) {
						auto v0 = v_map.find(he->vertex()->id());
						auto v1 = v_map.find(he->opposite()->vertex()->id());
						if (v0 != v_map.end() && v1 != v_map.end() && v0->second != v1->second) {
							skeleton.add_edge(v_map[he->vertex()->id()], v_map[he->opposite()->vertex()->id()]);
						}
					}
				}

				bool created;
				Surface_mesh::Property_map<Surface_mesh::Edge_index, int> edge_prop;
				boost::tie(edge_prop, created) = skeleton.add_property_map<Surface_mesh::Edge_index, int>("prop",0);
				assert(created);
				Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> red;
				boost::tie(red, created) = skeleton.add_property_map<Surface_mesh::Vertex_index, unsigned char>("red", LABELS[lab].red);
				assert(created);
				Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> green;
				boost::tie(green, created) = skeleton.add_property_map<Surface_mesh::Vertex_index, unsigned char>("green", LABELS[lab].green);
				assert(created);
				Surface_mesh::Property_map<Surface_mesh::Vertex_index, unsigned char> blue;
				boost::tie(blue, created) = skeleton.add_property_map<Surface_mesh::Vertex_index, unsigned char>("blue", LABELS[lab].blue);
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
				skeleton_name << "path_" << lab << "_" << i << ".ply";
				std::ofstream mesh_ofile (skeleton_name.str().c_str());
				CGAL::IO::write_PLY (mesh_ofile, skeleton);
				mesh_ofile.close();
			}

		}
	}

	return medial_axes;

}

std::set<std::pair<skeletonPoint,skeletonPoint>> link_paths(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> &medial_axes, const Raster &raster) {

	// Get label property
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::set<std::pair<skeletonPoint,skeletonPoint>> result;

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

						result.insert(std::make_pair(skeletonPoint(path2, v2), skeletonPoint(path1, e1, p1)));
						exit3:
							continue;
					}

				} else if (path1 == path2) {

					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Vertex_handle>, K::FT> distance_vv;
					std::map<std::pair<CGAL::Straight_skeleton_2<K>::Vertex_handle, CGAL::Straight_skeleton_2<K>::Halfedge_handle>, std::pair<K::FT, K::Point_2>> distance_vh;

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
						bool intersect = false;
						for (auto edge = path_polygon.at(path1).outer_boundary().edges_begin(); !intersect && edge != path_polygon.at(path1).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect = true;
							}
						}
						for (auto hole: path_polygon.at(path1).holes()) {
							for (auto edge = hole.edges_begin(); !intersect && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect = true;
								}
							}
						}
						if (!intersect) continue;

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
						bool intersect = false;
						for (auto edge = path_polygon.at(path1).outer_boundary().edges_begin(); !intersect && edge != path_polygon.at(path1).outer_boundary().edges_end(); edge++) {
							if (CGAL::do_intersect(*edge, segment)) {
								intersect = true;
							}
						}
						for (auto hole: path_polygon.at(path1).holes()) {
							for (auto edge = hole.edges_begin(); !intersect && edge != hole.edges_end(); edge++) {
								if (CGAL::do_intersect(*edge, segment)) {
									intersect = true;
								}
							}
						}
						if (!intersect) continue;
									
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
