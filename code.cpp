#include "header.hpp"
#include "raster.hpp"

#include <regex>
#include <unordered_map>

#include <ogrsf_frmts.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

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
