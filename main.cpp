#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <list>
#include <vector>

#include "ogrsf_frmts.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <CGAL/draw_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2                                          Point_2;
typedef K::Point_3                                          Point_3;
typedef CGAL::Point_set_3<Point_3>                          Point_set_3;
typedef CGAL::Surface_mesh<Point_3>                         Surface_mesh;
typedef CGAL::Polygon_2<K>                                  Polygon;
typedef CGAL::Projection_traits_xy_3<K>                     Projection_traits;
typedef CGAL::Delaunay_triangulation_2<Projection_traits>   TIN;

GDALRasterBand* get_band_from_TIFF(char* path, int i = 1) {
	GDALAllRegister();

	GDALDataset* dataset;

	dataset = (GDALDataset *) GDALOpen( path, GA_ReadOnly );
	if( dataset == NULL ) {
		std::cerr << "Unable to open " << path << "." << std::endl;
		return NULL;
	}

	GDALRasterBand *band = dataset->GetRasterBand(i);

	return band;
}

std::list<Polygon> get_LOD0_from_shapefile(char* path) {
	GDALAllRegister();

	GDALDataset* dataset;

	std::list<Polygon> polygons;

	dataset = (GDALDataset*) GDALOpenEx(path, GDAL_OF_VECTOR, NULL, NULL, NULL );
	if( dataset == NULL ) {
		std::cerr << "Unable to open " << path << "." << std::endl;
		return polygons;
	}

	for( OGRLayer* layer: dataset->GetLayers() ) {
        for( const auto& feature: *layer ) {
			if ((*feature)["num_class"].GetInteger() == 4) {
				// Get only the buildings
				OGRGeometry* geometry = feature->GetGeometryRef();
				
				if (wkbFlatten(geometry->getGeometryType()) == wkbPolygon) {

					OGRPolygon * shp_polygon = geometry->toPolygon();
					Polygon cgal_polygon;
					std::istringstream wkt(shp_polygon->exportToWkt());
					CGAL::IO::read_polygon_WKT (wkt, cgal_polygon);
					polygons.push_back(cgal_polygon);

				} else if (wkbFlatten(geometry->getGeometryType()) == wkbMultiPolygon) {

					OGRMultiPolygon * shp_multi_polygon = geometry->toMultiPolygon();
					for (OGRPolygon* shp_polygon: *shp_multi_polygon) {
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

void compute_object(GDALRasterBand* land_use, std::vector<std::vector<int>>& object_id, std::vector<std::pair<std::list<std::pair<int,int>>, int>>& objects, int x, int y, int id) {
	if (x >= 0 && x < land_use->GetXSize() && y >= 0 && y < land_use->GetYSize()) { // Valid pixel
		if (object_id[x][y] == -1) { // Note labeld pixel
			// Get label
			int label;
			if (land_use->RasterIO(GF_Read, x, y, 1, 1, &label, 1, 1, GDT_Int64, 0, 0) == CE_None) {
				if (objects[id].second == -1) {
					objects[id].second = label;
				}
				if (objects[id].second == label) { // Same object th id
					object_id[x][y] = id;
					objects[id].first.push_back({x, y});
					// Check neighbors
					compute_object(land_use, object_id, objects, x-1, y, id);
					compute_object(land_use, object_id, objects, x+1, y, id);
					compute_object(land_use, object_id, objects, x, y-1, id);
					compute_object(land_use, object_id, objects, x, y+1, id);
				}
			} else {
				std::cerr << "Error while reading (" << x << ", " << y << ") value of land use map." << std::endl;
			}
		}
	}
}

int compute_LOD2(char* DSM, char* DTM, char* land_use_map, char* LOD0, char* orthophoto) {
	// Get raset info as GDALRasterBand and LOD0 as polygons
	GDALRasterBand* dsm = get_band_from_TIFF(DSM);
	std::cout << dsm->GetXSize() << "x" << dsm->GetYSize() << " DSM load." << std::endl;

	GDALRasterBand* dtm = get_band_from_TIFF(DTM);
	std::cout << dtm->GetXSize() << "x" << dtm->GetYSize() << " DTM load." << std::endl;

	GDALRasterBand* land_use = get_band_from_TIFF(land_use_map);
	std::cout << land_use->GetXSize() << "x" << land_use->GetYSize() << " land use map load." << std::endl;
	
	if (LOD0 != NULL) {
		std::list<Polygon> buildings = get_LOD0_from_shapefile(LOD0);
		std::cout << buildings.size() << " buildings load from " << LOD0 << "." << std::endl;
	}

	if (orthophoto != NULL) {
		GDALRasterBand* red = get_band_from_TIFF(orthophoto, 1);
		GDALRasterBand* green = get_band_from_TIFF(orthophoto, 2);
		GDALRasterBand* blue = get_band_from_TIFF(orthophoto, 3);
		std::cout << red->GetXSize() << "x" << red->GetYSize() << " orthophoto load." << std::endl;
	}

	// Get objects (connected components) from land use map
	std::vector<std::vector<int>> object_id(land_use->GetXSize(), std::vector<int>(land_use->GetYSize(), -1));
	std::vector<std::pair<std::list<std::pair<int,int>>, int>> objects;

	for (int x = 0; x < land_use->GetXSize(); x++) {
		for (int y = 0; y < land_use->GetYSize(); y++) {
			if (object_id[x][y] == -1) {
				int id = objects.size();
				objects.push_back({{}, -1});
				compute_object(land_use, object_id, objects, x, y, id);
				//std::cout << "id: " << id << "\nnum_points: " << objects[0].first.size();
			}
		}
	}

	int id = 0;
	for (auto object: objects) {
		if (object.first.size() >= 3) {
			Point_set_3 points;
			points.reserve(object.first.size());
			for (auto p : object.first) {
				double z;
				if (dsm->RasterIO(GF_Read, p.first, p.second, 1, 1, &z, 1, 1, GDT_Float64, 0, 0) == CE_None) {
					points.insert(Point_3(p.first, p.second, z));
				}
			}

			TIN tin (points.points().begin(), points.points().end());
			
			if (tin.number_of_faces() > 0) {
				Surface_mesh mesh;
				CGAL::copy_face_graph (tin, mesh);

				std::ofstream mesh_ofile ("mesh_" + std::to_string(id++) + ".ply", std::ios_base::binary);
				CGAL::IO::set_binary_mode (mesh_ofile);
				CGAL::IO::write_PLY (mesh_ofile, mesh);
				mesh_ofile.close();
			}
		}
	}

	return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
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

	char* DSM = NULL;
	char* DTM = NULL;
	char* land_use_map = NULL;
	char* LOD0 = NULL;
	char* orthophoto = NULL;

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
