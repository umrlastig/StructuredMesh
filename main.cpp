#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <regex>

#include "ogrsf_frmts.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include "label.cpp"

typedef CGAL::Simple_cartesian<float>                       K;
typedef K::Point_2                                          Point_2;
typedef K::Point_3                                          Point_3;
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

		OGRSpatialReference get_crs() const {
			return OGRSpatialReference(crs);
		}

		template <typename T>
		std::list<std::pair<int,int>> triangle_to_pixel(CGAL::Point_3<T> a, CGAL::Point_3<T> b, CGAL::Point_3<T> c) const {
			std::list<std::pair<int,int>> ret;
			int min_x = std::max((int) std::min({a.x(), b.x(), c.x()}), 0);
			int max_x = std::min((int) std::max({a.x(), b.x(), c.x()}), xSize-1);
			int min_y = std::max((int) std::min({a.y(), b.y(), c.y()}), 0);
			int max_y = std::min((int) std::max({a.y(), b.y(), c.y()}), ySize-1);
			K::Triangle_2 triangle (Point_2(a.x(), a.y()), Point_2(b.x(), b.y()), Point_2(c.x(), c.y()));
			for (int L = min_y; L <= max_y; L++) {
				for (int P = min_x; P <= max_x; P++) {
					if (triangle.bounded_side(Point_2(0.5 + P, 0.5 + L)) != CGAL::ON_UNBOUNDED_SIDE) {
						ret.push_back(std::pair<int,int>(P,L));
					}
				}
			}
			return ret;
		}

		void fill_holes() {
			int max_hole_size = 300;

			//Negative holes
			{
				// Find holes
				std::list<std::pair<float,std::list<std::pair<int,int>>>> holes;
				for (int L = 1; L < ySize - 1; L++) {
					for (int P = 1; P < xSize -1; P++) {
						if (land_cover[L][P] >= 1 && land_cover[L][P] <= 3 || land_cover[L][P] >= 6) {
							if (dsm[L][P] <= dsm[L+1][P] & dsm[L][P] <= dsm[L-1][P] & dsm[L][P] <= dsm[L][P+1] & dsm[L][P] <= dsm[L][P-1]) {
								holes.push_back(std::pair<float,std::list<std::pair<int,int>>>(dsm[L][P], {std::pair<int,int>(P,L)}));
							}
						}
					}
				}
				
				// Compute hole elevation
				for (auto &hole: holes) {
					bool end = false;
					std::list<std::pair<int,int>> ring ({
						std::pair<int,int>(hole.second.front().first-1,hole.second.front().second),
						std::pair<int,int>(hole.second.front().first+1,hole.second.front().second),
						std::pair<int,int>(hole.second.front().first,hole.second.front().second-1),
						std::pair<int,int>(hole.second.front().first,hole.second.front().second+1)
						});

					while(!end) {
						float min_dsm = FLT_MAX;
						std::pair<int,int> min_pixel;

						for (auto pixel: ring) {
							if (dsm[pixel.second][pixel.first] < min_dsm) {
								min_pixel = pixel;
								min_dsm = dsm[pixel.second][pixel.first];
							}
						}

						hole.first = min_dsm;
						hole.second.push_back(min_pixel);
						ring.remove(min_pixel);

						if (min_pixel.first-1 >= 0) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first-1, min_pixel.second)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first-1, min_pixel.second)) == ring.end()) {
									if (dsm[min_pixel.second][min_pixel.first-1] < hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(min_pixel.first-1, min_pixel.second));
									}
								}
							}
						}
						if (min_pixel.first+1 < xSize) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first+1, min_pixel.second)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first+1, min_pixel.second)) == ring.end()) {
									if (dsm[min_pixel.second][min_pixel.first+1] < hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(min_pixel.first+1, min_pixel.second));
									}
								}
							}
						}
						if (min_pixel.second-1 >= 0) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first, min_pixel.second-1)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first, min_pixel.second-1)) == ring.end()) {
									if (dsm[min_pixel.second-1][min_pixel.first] < hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(min_pixel.first, min_pixel.second-1));
									}
								}
							}
						}
						if (min_pixel.second+1 < ySize) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first, min_pixel.second+1)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first, min_pixel.second+1)) == ring.end()) {
									if (dsm[min_pixel.second+1][min_pixel.first] < hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(min_pixel.first, min_pixel.second+1));
									}
								}
							}
						}

						if (hole.second.size() >= max_hole_size) {
							end = true;
						}
					}

					//Fill holes
					for (auto pixel: hole.second) {
						dsm[pixel.second][pixel.first] = hole.first;
					}
				}
			}

			//Positive holes
			{
				// Find holes
				std::list<std::pair<float,std::list<std::pair<int,int>>>> holes;
				for (int L = 1; L < ySize - 1; L++) {
					for (int P = 1; P < xSize -1; P++) {
						if (land_cover[L][P] >= 1 && land_cover[L][P] <= 3 || land_cover[L][P] >= 6) {
							if (dsm[L][P] >= dsm[L+1][P] & dsm[L][P] >= dsm[L-1][P] & dsm[L][P] >= dsm[L][P+1] & dsm[L][P] >= dsm[L][P-1]) {
								holes.push_back(std::pair<float,std::list<std::pair<int,int>>>(dsm[L][P], {std::pair<int,int>(P,L)}));
							}
						}
					}
				}
				
				// Compute hole elevation
				for (auto &hole: holes) {
					bool end = false;
					std::list<std::pair<int,int>> ring ({
						std::pair<int,int>(hole.second.front().first-1,hole.second.front().second),
						std::pair<int,int>(hole.second.front().first+1,hole.second.front().second),
						std::pair<int,int>(hole.second.front().first,hole.second.front().second-1),
						std::pair<int,int>(hole.second.front().first,hole.second.front().second+1)
						});

					while(!end) {
						float max_dsm = FLT_MIN;
						std::pair<int,int> max_pixel;

						for (auto pixel: ring) {
							if (dsm[pixel.second][pixel.first] > max_dsm) {
								max_pixel = pixel;
								max_dsm = dsm[pixel.second][pixel.first];
							}
						}

						hole.first = max_dsm;
						hole.second.push_back(max_pixel);
						ring.remove(max_pixel);

						if (max_pixel.first-1 >= 0) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first-1, max_pixel.second)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first-1, max_pixel.second)) == ring.end()) {
									if (dsm[max_pixel.second][max_pixel.first-1] > hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(max_pixel.first-1, max_pixel.second));
									}
								}
							}
						}
						if (max_pixel.first+1 < xSize) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first+1, max_pixel.second)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first+1, max_pixel.second)) == ring.end()) {
									if (dsm[max_pixel.second][max_pixel.first+1] > hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(max_pixel.first+1, max_pixel.second));
									}
								}
							}
						}
						if (max_pixel.second-1 >= 0) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first, max_pixel.second-1)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first, max_pixel.second-1)) == ring.end()) {
									if (dsm[max_pixel.second-1][max_pixel.first] > hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(max_pixel.first, max_pixel.second-1));
									}
								}
							}
						}
						if (max_pixel.second+1 < ySize) {
							if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first, max_pixel.second+1)) == hole.second.end()) {
								if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first, max_pixel.second+1)) == ring.end()) {
									if (dsm[max_pixel.second+1][max_pixel.first] > hole.first) {
										end = true;
									} else {
										ring.push_back(std::pair<int,int>(max_pixel.first, max_pixel.second+1));
									}
								}
							}
						}

						if (hole.second.size() >= max_hole_size) {
							end = true;
						}
					}

					//Fill holes
					for (auto pixel: hole.second) {
						dsm[pixel.second][pixel.first] = hole.first;
					}
				}
			}
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

			fill_holes();
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

float face_cost(const Raster &raster, const Point_3 &p0, const Point_3 &p1, const Point_3 &p2) {
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
		verticality = nz;
	}

	// Eccentricity
	Point_3 S = CGAL::centroid(p0,p1,p2);
	float M = (pow(p2.x()-S.x(),2) + pow(p2.y()-S.y(),2) + pow(p2.z()-S.z(),2) + (pow(p1.x()-p0.x(),2) + pow(p1.y()-p0.y(),2) + pow(p1.z()-p0.z(),2))/3)/4;
	float N = sqrtf(pow((p2.y()-S.y())*(p1.z()-p0.z())-(p2.z()-S.z())*(p1.y()-p0.y()),2) + pow((p2.z()-S.z())*(p1.x()-p0.x())-(p2.x()-S.x())*(p1.z()-p0.z()),2) + pow((p2.x()-S.x())*(p1.y()-p0.y())-(p2.y()-S.y())*(p1.x()-p0.x()),2))/(4 * sqrtf(3));
	float eccentricity = M*M - 4*N*N;
	eccentricity = sqrtf(1-(M-eccentricity)/(M+eccentricity));

	float surface = K::Triangle_3(p0, p1, p2).squared_area ();
	return surface * (0 * least_squares + 1 * least + 1 * entropy + 0 * verticality + 0 * eccentricity);
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
			
			for (int j = 0; j < 5; j++) {
				p[j] = best_point(raster, p[j].x(), p[j].y(), profile);
			}

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

			for (int i = 0; i < 1; i++) {
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
				p[1] = best_point(raster, p[1].x(), p[1].y(), profile);
				p[3] = best_point(raster, p[1].x(), p[1].y(), profile);

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

	// Entropy
	OutputMesh::Property_map<OutputMesh::Face_index, float> quality;
	boost::tie(quality, created) = output_mesh.add_property_map<OutputMesh::Face_index, float>("quality",0);
	assert(created);

	for (auto face : output_mesh.faces()) {
		int face_label[LABELS.size()] = {0};
		int sum_face_label = 0;

		CGAL::Vertex_around_face_iterator<OutputMesh> vbegin, vend;
    	boost::tie(vbegin, vend) = vertices_around_face(output_mesh.halfedge(face), output_mesh);
		for (auto pixel : raster.triangle_to_pixel(output_mesh.point(*(vbegin++)), output_mesh.point(*(vbegin++)), output_mesh.point(*(vbegin++)))) {
			if (raster.land_cover[pixel.second][pixel.first] > -1) {
				sum_face_label++;
				face_label[raster.land_cover[pixel.second][pixel.first]]++;
			}
		}
		
		auto argmax = std::max_element(face_label, face_label+LABELS.size());
		label[face] = argmax - face_label;
		red[face] = LABELS[argmax - face_label].red;
		green[face] = LABELS[argmax - face_label].green;
		blue[face] = LABELS[argmax - face_label].blue;

		/*float entropy = 0;
		if (sum_face_label > 0) {
			for (int i = 0; i < LABELS.size(); i++) {
				if (face_label[i] > 0) {
					entropy += ((float) face_label[i])*log((float) face_label[i]);
				}
			}
			entropy = log((float) sum_face_label) - entropy/((float) sum_face_label);
		}
		quality[face] = entropy;*/
		boost::tie(vbegin, vend) = vertices_around_face(output_mesh.halfedge(face), output_mesh);
		auto pa = output_mesh.point(*(vbegin++));
		auto pb = output_mesh.point(*(vbegin++));
		auto pc = output_mesh.point(*(vbegin++));
		quality[face] = face_cost(raster, Point_3(pa.x(), pa.y(), pa.z()), Point_3(pb.x(), pb.y(), pb.z()), Point_3(pc.x(), pc.y(), pc.z()));
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
		
		void OnCollected(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<float>& cost) {
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
				std::cout << "\rCollapse: " << (initial_edge_count-current_edge_count) << "/" << initial_edge_count << " (" << ((int) (((float) (initial_edge_count-current_edge_count))/initial_edge_count*100)) << "%)" << " still " << (((float) current_edge_count) * diff.count() / (initial_edge_count-current_edge_count)) << "s" << " (" << (((float) (initial_edge_count-current_edge_count)) / diff.count()) << " op/s) - cost: " << *cost << "     " << std::flush;
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

int compute_LOD2(char *DSM, char *DTM, char *land_use_map, char *LOD0, char *orthophoto) {
	const Raster raster(DSM, DTM, land_use_map);

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

	// Cost_stop_predicate stop(10);
	SMS::Count_stop_predicate<Surface_mesh> stop(1000);
	Custom_cost cf(raster);
	Custom_placement pf(raster);
	int r = SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).get_placement(pf).visitor(My_visitor(&mesh, &raster)));
	std::cout << "\rMesh simplified                                               " << std::endl;

	save_mesh(mesh, raster, "final-mesh.ply");

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
