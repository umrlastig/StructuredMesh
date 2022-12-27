#include "header.hpp"
#include "raster.hpp"

#include <CGAL/IO/WKT.h>
#include <CGAL/boost/graph/copy_face_graph.h>

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

float single_face_cost(const Raster &raster, const Point_3 &p0, const Point_3 &p1, const Point_3 &p2);

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
