#include "header.hpp"
#include "raster.hpp"

#include <regex>
#include <unordered_map>

#include <ogrsf_frmts.h>
#include <CGAL/IO/WKT.h>

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
	boost::tie(label, created) = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("f:label",0);
	assert(created);

	for (auto face : mesh.faces()) {
		int face_label[LABELS.size()] = {0};
		int sum_face_label = 0;

		CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
		boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
		for (auto pixel : raster.triangle_to_pixel(mesh.point(*(vbegin++)), mesh.point(*(vbegin++)), mesh.point(*(vbegin++)))) {
			sum_face_label++;
			face_label[raster.land_cover[pixel.second][pixel.first]]++;
		}

		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		label[face] = argmax - face_label;
	}
}

void change_vertical_faces(Surface_mesh &mesh) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("f:label");
	assert(has_label);

	for (auto face : mesh.faces()) {
		if (label[face] == LABEL_RAIL) {

			CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
			boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
			auto p0 = mesh.point(*(vbegin++));
			auto p1 = mesh.point(*(vbegin++));
			auto p2 = mesh.point(*(vbegin++));

			auto n = CGAL::orthogonal_vector(p0, p1, p2);
			if (CGAL::scalar_product(n, K::Vector_3(0,0,1)) / CGAL::sqrt(n.squared_length()) < 0.98) {
				label[face] = LABEL_OTHER;
			}
		} else if (label[face] == LABEL_ROAD) {

			CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
			boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
			auto p0 = mesh.point(*(vbegin++));
			auto p1 = mesh.point(*(vbegin++));
			auto p2 = mesh.point(*(vbegin++));

			auto n = CGAL::orthogonal_vector(p0, p1, p2);
			if (CGAL::scalar_product(n, K::Vector_3(0,0,1)) / CGAL::sqrt(n.squared_length()) < 0.95) {
				label[face] = LABEL_OTHER;
			}
		}
	}

}

void compute_normal_angle_coef(Surface_mesh &mesh) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> normal_angle_coef;
	bool created_normal_angle_coef;
	boost::tie(normal_angle_coef, created_normal_angle_coef) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:n_a_coef", 1);
	assert(created_normal_angle_coef);

	for (auto face : mesh.faces()) {
		CGAL::Vertex_around_face_iterator<Surface_mesh> vbegin, vend;
		boost::tie(vbegin, vend) = vertices_around_face(mesh.halfedge(face), mesh);
		auto p0 = mesh.point(*(vbegin++));
		auto p1 = mesh.point(*(vbegin++));
		auto p2 = mesh.point(*(vbegin++));
		auto n = CGAL::orthogonal_vector (p0, p1, p2);
		normal_angle_coef[face] = sin(-CGAL::approximate_angle(n, K::Vector_3(0,0,1)) * M_PI / 180) + 1;
	}
}
