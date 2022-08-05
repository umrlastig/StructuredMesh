#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <list>
#include <vector>
#include <algorithm>

#include "ogrsf_frmts.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include "label.cpp"

typedef CGAL::Simple_cartesian<float>                       K;
typedef K::Point_2                                          Point_2;
typedef K::Point_3                                          Point_3;
typedef K::Triangle_2                                       Triangle_2;
typedef CGAL::Surface_mesh<Point_3>                         Surface_mesh;
typedef CGAL::Polygon_2<K>                                  Polygon;

namespace SMS = CGAL::Surface_mesh_simplification;

class Raster {
	public:
		std::vector<std::vector<float>> dsm;
		std::vector<std::vector<float>> dtm;
		std::vector<std::vector<char>> land_cover;
		int xSize;
		int ySize;

	private:
		OGRSpatialReference crs;
		double grid_to_crs[6] = {0,1,0,0,0,1};

		static std::pair<int,int> grid_conversion(int P, int L, double grid_to_crs[6], OGRCoordinateTransformation *crs_to_other_crs, double other_grid_to_other_crs[6]) {
			double x = grid_to_crs[0] + (0.5 + P)*grid_to_crs[1] + (0.5 + L)*grid_to_crs[2];
			double y = grid_to_crs[3] + (0.5 + P)*grid_to_crs[4] + (0.5 + L)*grid_to_crs[5];
			crs_to_other_crs->Transform(1,&x,&y);
			double fact = other_grid_to_other_crs[1]*other_grid_to_other_crs[5] - other_grid_to_other_crs[2]*other_grid_to_other_crs[4];
			x -= other_grid_to_other_crs[0];
			y -= other_grid_to_other_crs[3];
			int newP = ((int) ((other_grid_to_other_crs[5]*x - other_grid_to_other_crs[2]*y) / fact));
			int newL = ((int) ((-other_grid_to_other_crs[4]*x + other_grid_to_other_crs[1]*y) / fact));
			return std::pair<int,int>(newP, newL);
		}

	public:
		void coord_to_grid(double x, double y, float& P, float& L) const {
			double fact = grid_to_crs[1]*grid_to_crs[5] - grid_to_crs[2]*grid_to_crs[4];
			x -= grid_to_crs[0];
			y -= grid_to_crs[3];
			P = (grid_to_crs[5]*x - grid_to_crs[2]*y) / fact;
			L = (-grid_to_crs[4]*x + grid_to_crs[1]*y) / fact;
		}

		void grid_to_coord(int P, int L, double& x, double& y) const {
			x = grid_to_crs[0] + (0.5 + P)*grid_to_crs[1] + (0.5 + L)*grid_to_crs[2];
			y = grid_to_crs[3] + (0.5 + P)*grid_to_crs[4] + (0.5 + L)*grid_to_crs[5];
		}

		void grid_to_coord(float P, float L, double& x, double& y) const {
			x = grid_to_crs[0] + P*grid_to_crs[1] + L*grid_to_crs[2];
			y = grid_to_crs[3] + P*grid_to_crs[4] + L*grid_to_crs[5];
		}

		template <typename T>
		std::list<std::pair<int,int>> triangle_to_pixel(CGAL::Point_3<T> a, CGAL::Point_3<T> b, CGAL::Point_3<T> c) const {
			std::list<std::pair<int,int>> ret;
			int min_x = std::max((int) std::min({a.x(), b.x(), c.x()}), 0);
			int max_x = std::min((int) std::max({a.x(), b.x(), c.x()}), xSize-1);
			int min_y = std::max((int) std::min({a.y(), b.y(), c.y()}), 0);
			int max_y = std::min((int) std::max({a.y(), b.y(), c.y()}), ySize-1);
			Triangle_2 triangle (Point_2(a.x(), a.y()), Point_2(b.x(), b.y()), Point_2(c.x(), c.y()));
			for (int L = min_y; L <= max_y; L++) {
				for (int P = min_x; P <= max_x; P++) {
					if (triangle.bounded_side(Point_2(0.5 + P, 0.5 + L)) != CGAL::ON_UNBOUNDED_SIDE) {
						ret.push_back(std::pair<int,int>(P,L));
					}
				}
			}
			return ret;
		}

		Raster(char *dsm_path, char *dtm_path, char *land_cover_path) {
			GDALAllRegister();

			// Get DSM informations and CRS
			GDALDataset *dsm_dataset = (GDALDataset *) GDALOpen(dsm_path, GA_ReadOnly );
			if( dsm_dataset == NULL ) {
				throw std::invalid_argument(std::string("Unable to open ") + dsm_path + ".");
			}
			xSize = dsm_dataset->GetRasterBand(1)->GetXSize();
			ySize = dsm_dataset->GetRasterBand(1)->GetYSize();
			crs = *dsm_dataset->GetSpatialRef();
			if (dsm_dataset->GetGeoTransform(grid_to_crs) >= CE_Failure) {
				throw std::invalid_argument(std::string(dsm_path) + " do not contain an affine transform.");
			}
			dsm = std::vector<std::vector<float>>(ySize, std::vector<float>(xSize, 0));
			for (int L = 0; L < ySize; L++) {
				if (dsm_dataset->GetRasterBand(1)->RasterIO(GF_Read, 0, L, xSize, 1, &dsm[L][0], xSize, 1, GDT_Float32, 0, 0) >= CE_Failure) {
					throw std::invalid_argument(std::string(dsm_path) + " can't be read.");
				}
			}
			std::cout << "DSM load" << std::endl;

			// Get DTM informations and CRS transform
			dtm = std::vector<std::vector<float>>(ySize, std::vector<float>(xSize, 0));
			GDALDataset *dtm_dataset = (GDALDataset *) GDALOpen( dtm_path, GA_ReadOnly );
			if( dtm_dataset == NULL ) {
				throw std::invalid_argument(std::string("Unable to open ") + dtm_path + ".");
			}
			double dtm_grid_to_dtm_crs[6];
			if (dtm_dataset->GetGeoTransform(dtm_grid_to_dtm_crs) >= CE_Failure) {
				throw std::invalid_argument("Can't transform DTM grid to DTM CRS.");
			}
			OGRCoordinateTransformation *CRS_to_dtm_crs = OGRCreateCoordinateTransformation(
				&crs,
				dtm_dataset->GetSpatialRef());
			if (CRS_to_dtm_crs == NULL) {
				throw std::runtime_error("Can't transform DSM CRS to DTM CRS.");
			}
			for (int L = 0; L < ySize; L++) {
				for (int P = 0; P < xSize; P++) {
					std::pair<int,int> new_coord = grid_conversion(P, L, grid_to_crs, CRS_to_dtm_crs, dtm_grid_to_dtm_crs);
					if (new_coord.first >= 0 && new_coord.first < dtm_dataset->GetRasterBand(1)->GetXSize() && new_coord.second >= 0 && new_coord.second < dtm_dataset->GetRasterBand(1)->GetYSize()) {
						if (dtm_dataset->GetRasterBand(1)->RasterIO(GF_Read, new_coord.first, new_coord.second, 1, 1, &dtm[L][P], 1, 1, GDT_Float32, 0, 0) >= CE_Failure) {
							throw std::invalid_argument(std::string(dtm_path) + " can't be read.");
						}
					}
				}
			}
			std::cout << "DTM load" << std::endl;

			// Get land cover informations and CRS transform
			land_cover = std::vector<std::vector<char>>(ySize, std::vector<char>(xSize, -1));
			GDALDataset *land_cover_dataset = (GDALDataset *) GDALOpen( land_cover_path, GA_ReadOnly );
			if( dsm_dataset == NULL ) {
				throw std::invalid_argument(std::string("Unable to open ") + land_cover_path + ".");
			}

			double land_cover_grid_to_land_cover_crs[6];
			if (land_cover_dataset->GetGeoTransform(land_cover_grid_to_land_cover_crs) >= CE_Failure) {
				throw std::invalid_argument("Can't transform land cover grid to land cover CRS.");
			}
			OGRCoordinateTransformation *CRS_to_land_cover_crs = OGRCreateCoordinateTransformation(
				&crs,
				land_cover_dataset->GetSpatialRef());
			if (CRS_to_land_cover_crs == NULL) {
				throw std::runtime_error("Can't transform DSM CRS to land cover CRS.");
			}
			for (int L = 0; L < ySize; L++) {
				for (int P = 0; P < xSize; P++) {
					std::pair<int,int> new_coord = grid_conversion(P, L, grid_to_crs, CRS_to_land_cover_crs, land_cover_grid_to_land_cover_crs);
					if (new_coord.first >= 0 && new_coord.first < land_cover_dataset->GetRasterBand(1)->GetXSize() && new_coord.second >= 0 && new_coord.second < land_cover_dataset->GetRasterBand(1)->GetYSize()) {
						unsigned char value;
						if (land_cover_dataset->GetRasterBand(1)->RasterIO(GF_Read, new_coord.first, new_coord.second, 1, 1, &value, 1, 1, GDT_Byte, 0, 0) >= CE_Failure) {
							throw std::invalid_argument(std::string(land_cover_path) + " can't be read.");
						}
						land_cover[L][P] = (char) value;
					}
				}
			}
			std::cout << "Land cover load" << std::endl;
		}
};

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

int compute_LOD2(char *DSM, char *DTM, char *land_use_map, char *LOD0, char *orthophoto) {
	const Raster raster(DSM, DTM, land_use_map);

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

	SMS::Count_ratio_stop_predicate<Surface_mesh> stop(0.05);
	int r = SMS::edge_collapse(mesh, stop);
	std::cout << "Mesh simplified" << std::endl;

	typedef CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3> OutputMesh;
	OutputMesh output_mesh;
	CGAL::copy_face_graph (mesh, output_mesh);

	// Label
	OutputMesh::Property_map<OutputMesh::Face_index, unsigned char> label;
	bool created;
	boost::tie(label, created) = output_mesh.add_property_map<OutputMesh::Face_index, unsigned char>("label",0);
	assert(created);

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

	for (auto face : output_mesh.faces()) {
		int face_label[LABELS.size()] = {0};

		CGAL::Vertex_around_face_iterator<OutputMesh> vbegin, vend;
    	boost::tie(vbegin, vend) = vertices_around_face(output_mesh.halfedge(face), output_mesh);
		for (auto pixel : raster.triangle_to_pixel(output_mesh.point(*vbegin++), output_mesh.point(*(vbegin++)), output_mesh.point(*(vbegin++)))) {
			if (raster.land_cover[pixel.second][pixel.first] > -1) {
				face_label[raster.land_cover[pixel.second][pixel.first]]++;
			}
		}
		
		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		label[face] = argmax - face_label;
		red[face] = LABELS[argmax - face_label].red;
		green[face] = LABELS[argmax - face_label].green;
		blue[face] = LABELS[argmax - face_label].blue;
	}

	for(auto vertex : output_mesh.vertices()) {
		auto point = output_mesh.point(vertex);
		double x, y;
		raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
		output_mesh.point(vertex) = CGAL::Simple_cartesian<double>::Point_3(x, y, (double) point.z());
	}

	std::ofstream mesh_ofile ("mesh.ply", std::ios_base::binary);
	CGAL::IO::set_binary_mode (mesh_ofile);
	CGAL::IO::write_PLY (mesh_ofile, output_mesh);
	mesh_ofile.close();

	return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
	int opt;
	const struct option options[] = {
		{"help", no_argument, NULL, 'h'},
		{"DSM", required_argument, NULL, 's'},
		{"DTM", required_argument, NULL, 't'},
		{"land_use_map", required_argument, NULL, 'l'},
		{"LOD0", required_argument, NULL, '0'},
		{"orthophoto", required_argument, NULL, 'i'},
		{NULL, 0, 0, '\0'}
	};

	char *DSM = NULL;
	char *DTM = NULL;
	char *land_use_map = NULL;
	char *LOD0 = NULL;
	char *orthophoto = NULL;

	while ((opt = getopt_long(argc, argv, "hs:t:l:0:i:", options, NULL)) != -1) {
		switch(opt) {
			case 'h':
				std::cout << "Usage: " << argv[0] << " [OPTIONS] -s DSM -t DTM -l land_use_map" << std::endl;
				std::cout << "Build LOD2 representation from DSM, DTM and land use map." << std::endl;
				std::cout << "OPTIONS:" << std::endl;
				std::cout << " -h, --help                         Print this help anq quit." << std::endl;
				std::cout << " -s, --DSM=/file/path.tiff          DSM as TIFF file." << std::endl << std::endl;
				std::cout << " -t, --DTM=/file/path.tiff          DTM as TIFF file." << std::endl << std::endl;
				std::cout << " -l, --land_use_map=/file/path.tiff land use map as TIFF file." << std::endl << std::endl;
				std::cout << " -0, --LOD0=/file/path.shp          LOD0 as Shapefile." << std::endl << std::endl;
				std::cout << " -i, --orthophoto=/file/path.tiff   RGB orthophoto as TIFF file." << std::endl;
				return EXIT_SUCCESS;
				break;
			case 's':
				DSM = optarg;
				break;
			case 't':
				DTM = optarg;
				break;
			case 'l':
				land_use_map = optarg;
				break;
			case '0':
				LOD0 = optarg;
				break;
			case 'i':
				orthophoto = optarg;
				break;
		}
	}

	if (DSM == NULL || DTM == NULL || land_use_map == NULL) {
		std::cerr << "DSM, DTM and land_use_map are mandatory" << std::endl;
		std::cerr << "Usage: " << argv[0] << " [OPTIONS] -s DSM -t DTM -l land_use_map" << std::endl;
		std::cerr << "Use -h/--help to obtain more informations" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "DSM: " << DSM << std::endl;
	std::cout << "DTM: " << DTM << std::endl;
	std::cout << "Land use map: " << land_use_map << std::endl;
	if (LOD0 != NULL) {
		std::cout << "LOD0: " << LOD0 << std::endl;
	}
	if (orthophoto != NULL) {
		std::cout << "Orthophoto: " << orthophoto << std::endl;
	}
	std::cout << std::endl;

	return compute_LOD2(DSM, DTM, land_use_map, LOD0, orthophoto);
}
