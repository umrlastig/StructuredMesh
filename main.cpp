#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <list>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <regex>
#include <unordered_map>

#include "ogrsf_frmts.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Polygon_mesh_processing/locate.h>

#include "label.cpp"

typedef CGAL::Simple_cartesian<float>                       K;
typedef K::Point_2                                          Point_2;
typedef K::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                         Surface_mesh;
typedef CGAL::Polygon_2<K>                                  Polygon;

namespace SMS = CGAL::Surface_mesh_simplification;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Exact_predicates_kernel;
typedef CGAL::Arr_segment_traits_2<Exact_predicates_kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement_2;

namespace PMP = CGAL::Polygon_mesh_processing;

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

		void align_land_cover_and_dsm() {
			std::vector<std::vector<char>> new_land_cover = std::vector<std::vector<char>>(ySize, std::vector<char>(xSize, -1));
			for (int L = 0; L < ySize; L++) {
				for (int P = 0; P < xSize; P++) {
					float neighboors_count[LABELS.size()] = {0};
					int n = 0;
					float moy = 0;
					float square_moy = 0;
					for (int i = -5; i <= 5; i++) {
						for (int j = -5; j <= 5; j++) {
							if (L+i >= 0 && P + j >= 0 && L+i < ySize && P+j < xSize) {
								neighboors_count[land_cover[L+i][P+j]] += 1/(0.1+abs(dsm[L][P] - dsm[L+i][P+j]));
								if (i > -3 && i < 3 && j > -3 && j < 3) {
									n++;
									moy += dsm[L+i][P+j];
									square_moy += pow(dsm[L+i][P+j], 2);
								}
							}
						}
					}
					if (pow(square_moy/n - pow(moy/n, 2), 0.5) > 0.2) {
						new_land_cover[L][P] = std::max_element(neighboors_count, neighboors_count + LABELS.size()) - neighboors_count;
					} else {
						new_land_cover[L][P] = land_cover[L][P];
					}

				}
			}

			land_cover = new_land_cover;
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

		float coord_distance_to_grid_distance(double d) const {
			float P, L;
			double x, y;
			grid_to_coord(0, 0, x, y);
			coord_to_grid(x+d/sqrt(2),y+d/sqrt(2),P,L);
			return sqrt(pow(P,2)+pow(L,2));
		}

		double grid_distance_to_coord_distance(float d) const {
			double x1, y1, x2, y2;
			grid_to_coord(0, 0, x1, y1);
			grid_to_coord((float) (d/sqrt(2)), (float) (d/sqrt(2)), x2, y2);
			return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
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

			align_land_cover_and_dsm();
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
				remove_face.push_back(face);
			}
		}
	}

	for (auto &face_label: new_label) {
		label[face_label.first] = face_label.second;
	}

	for(auto &face: remove_face) {
		CGAL::Euler::remove_face(mesh.halfedge(face), mesh);
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


std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> compute_medial_axes(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> medial_axes;

	typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
	for (int i = 0; i < paths.size(); i++) {

		Filtered_graph filtered_sm(mesh, i, path);

		if (filtered_sm.number_of_faces() > 1) {

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

				auto poly = polygon((*(arr.unbounded_face()->holes_begin()))->twin()->face());

				poly = CGAL::Polyline_simplification_2::simplify(
					poly,
					CGAL::Polyline_simplification_2::Squared_distance_cost(),
					CGAL::Polyline_simplification_2::Stop_above_cost_threshold(pow(raster.coord_distance_to_grid_distance(1.5),2))
				);

				auto iss = CGAL::create_interior_straight_skeleton_2(poly, Exact_predicates_kernel());
				medial_axes[i] = CGAL::convert_straight_skeleton_2<CGAL::Straight_skeleton_2<K>>(*iss);

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

				Surface_mesh part_mesh;
				CGAL::copy_face_graph(filtered_sm, part_mesh);
				std::stringstream name;
				name << "part_mesh_" << lab << "_" << i << ".ply";
				save_mesh(part_mesh, raster, name.str().c_str());
			}
		}
	}

	return medial_axes;

}

struct skeletonPoint {
	int path;
	Point_2 point;
	CGAL::Straight_skeleton_2<K>::Vertex_handle vertex;
	CGAL::Straight_skeleton_2<K>::Halfedge_handle halfedge;

	skeletonPoint () {}

	skeletonPoint (int _path, CGAL::Straight_skeleton_2<K>::Vertex_handle _vertex) {
		path = _path;
		vertex = _vertex;
		point = _vertex->point();
		halfedge = nullptr;
	}

	skeletonPoint (int _path, CGAL::Straight_skeleton_2<K>::Halfedge_handle _halfedge, K::Point_2 _point) {
		path = _path;
		halfedge = _halfedge;
		point = _point;
		vertex = nullptr;
	}

	friend bool operator<(const skeletonPoint& l, const skeletonPoint& r) {
		return l.point < r.point;
	}

	friend bool operator==(const skeletonPoint& l, const skeletonPoint& r) {
		return l.point == r.point;
	}
};

std::set<std::pair<skeletonPoint,skeletonPoint>> link_paths(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> &medial_axes, const Raster &raster) {

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::set<std::pair<skeletonPoint,skeletonPoint>> result;

	for (int selected_label:  {3, 8, 9}) {
		std::list<int> same_label_paths;

		for (int i = 0; i < paths.size(); i++) {
			if (label[paths[i].front()] == selected_label && medial_axes.count(i) == 1) {
				same_label_paths.push_back(i);
			}
		}

		for (int path1: same_label_paths) {
			for (int path2: same_label_paths) {
				if (path1 < path2) {
					std::map<skeletonPoint, skeletonPoint> nearests1, nearests2;

					//for each vertex of path1 compute nearest point on path2
					for (auto point1: medial_axes.at(path1)->vertex_handles()) {
						if (point1->is_skeleton()) {
							skeletonPoint source(path1, point1);
							skeletonPoint nearest;
							K::FT min_distance = -1;

							for (auto point2: medial_axes.at(path2)->vertex_handles()) {
								if (point2->is_skeleton()) {
									if (CGAL::squared_distance(point1->point(), point2->point()) < min_distance || min_distance == -1) {
										min_distance = CGAL::squared_distance(point1->point(), point2->point());
										nearest = skeletonPoint(path2, point2);
									}
								}
							}
							for (auto edge2: medial_axes.at(path2)->halfedge_handles()) {
								if (edge2->vertex()->id() < edge2->opposite()->vertex()->id()) {
									if(edge2->is_inner_bisector() && edge2->opposite()->is_inner_bisector()) {
										auto target = edge2->vertex()->point();
										auto source = edge2->opposite()->vertex()->point();
										auto proj = K::Line_2(source, target).projection(point1->point());
										if (K::Segment_2(source, target).collinear_has_on(proj)) {
											if (CGAL::squared_distance(point1->point(), proj) < min_distance) {
												min_distance = CGAL::squared_distance(point1->point(), proj);
												nearest = skeletonPoint(path2, edge2, proj);
											}
										}
									}
								}
							}
							nearests1[source] = nearest;
						}
					}

					//for each vertex of path2 compute nearest point on path1
					for (auto point2: medial_axes.at(path2)->vertex_handles()) {
						if (point2->is_skeleton()) {
							skeletonPoint source(path2, point2);
							skeletonPoint nearest;
							K::FT min_distance = -1;

							for (auto point1: medial_axes.at(path1)->vertex_handles()) {
								if (point1->is_skeleton()) {
									if (CGAL::squared_distance(point2->point(), point1->point()) < min_distance || min_distance == -1) {
										min_distance = CGAL::squared_distance(point2->point(), point1->point());
										nearest = skeletonPoint(path1, point1);
									}
								}
							}
							for (auto edge1: medial_axes.at(path1)->halfedge_handles()) {
								if (edge1->vertex()->id() < edge1->opposite()->vertex()->id()) {
									if(edge1->is_inner_bisector() && edge1->opposite()->is_inner_bisector()) {
										auto target = edge1->vertex()->point();
										auto source = edge1->opposite()->vertex()->point();
										auto proj = K::Line_2(source, target).projection(point2->point());
										if (K::Segment_2(source, target).collinear_has_on(proj)) {
											if (CGAL::squared_distance(point2->point(), proj) < min_distance) {
												min_distance = CGAL::squared_distance(point2->point(), proj);
												nearest = skeletonPoint(path1, edge1, proj);
											}
										}
									}
								}
							}
							nearests2[source] = nearest;
						}
					}

					//for projetion on path2, compute nearest point on path 1 if its a vertex
					for (auto points1: nearests1) {
						if (nearests2.count(points1.second) == 0) {
							skeletonPoint nearest = points1.first;
							K::FT min_distance = CGAL::squared_distance(points1.first.point, points1.second.point);

							for (auto point1: medial_axes.at(path1)->vertex_handles()) {
								if (point1->is_skeleton()) {
									if (CGAL::squared_distance(points1.second.point, point1->point()) < min_distance) {
										min_distance = CGAL::squared_distance(points1.second.point, point1->point());
										nearest = skeletonPoint(path1, point1);
									}
								}
							}
							for (auto edge1: medial_axes.at(path1)->halfedge_handles()) {
								if (edge1->vertex()->id() < edge1->opposite()->vertex()->id()) {
									if(edge1->is_inner_bisector() && edge1->opposite()->is_inner_bisector()) {
										auto target = edge1->vertex()->point();
										auto source = edge1->opposite()->vertex()->point();
										auto proj = K::Line_2(source, target).projection(points1.second.point);
										if (K::Segment_2(source, target).collinear_has_on(proj)) {
											if (CGAL::squared_distance(points1.second.point, proj) < min_distance) {
												min_distance = CGAL::squared_distance(points1.second.point, proj);
												nearest = skeletonPoint(path1, edge1, proj);
											}
										}
									}
								}
							}
							if (nearests1.count(nearest) == 1) {
								nearests2[points1.second] = nearest;
							}
						}
					}

					//for projetion on path1, compute nearest point on path 2 if its a vertex
					for (auto points2: nearests2) {
						if (nearests1.count(points2.second) == 0) {
							skeletonPoint nearest = points2.first;
							K::FT min_distance = CGAL::squared_distance(points2.first.point, points2.second.point);

							for (auto point2: medial_axes.at(path2)->vertex_handles()) {
								if (point2->is_skeleton()) {
									if (CGAL::squared_distance(points2.second.point, point2->point()) < min_distance) {
										min_distance = CGAL::squared_distance(points2.second.point, point2->point());
										nearest = skeletonPoint(path2, point2);
									}
								}
							}
							for (auto edge2: medial_axes.at(path2)->halfedge_handles()) {
								if (edge2->vertex()->id() < edge2->opposite()->vertex()->id()) {
									if(edge2->is_inner_bisector() && edge2->opposite()->is_inner_bisector()) {
										auto target = edge2->vertex()->point();
										auto source = edge2->opposite()->vertex()->point();
										auto proj = K::Line_2(source, target).projection(points2.second.point);
										if (K::Segment_2(source, target).collinear_has_on(proj)) {
											if (CGAL::squared_distance(points2.second.point, proj) < min_distance) {
												min_distance = CGAL::squared_distance(points2.second.point, proj);
												nearest = skeletonPoint(path2, edge2, proj);
											}
										}
									}
								}
							}
							if (nearests2.count(nearest) == 1) {
								nearests1[points2.second] = nearest;
							}
						}
					}

					for (auto points1: nearests1) {
						if (nearests2.count(points1.second) == 1 && nearests2.at(points1.second) == points1.first) {
							if (CGAL::squared_distance(points1.first.point, points1.second.point) < pow(raster.coord_distance_to_grid_distance(50),2)) {

								// slope
								float max_slope;
								if (selected_label == 8) {
									// slop for road
									max_slope = 0.3;
								} else {
									// slop for railways and water
									max_slope = 0.1;
								}
								auto z1 = raster.dsm[int(points1.first.point.y())][int(points1.first.point.x())];
								auto z2 = raster.dsm[int(points1.second.point.y())][int(points1.second.point.x())];

								if (abs(z1 - z2) / raster.grid_distance_to_coord_distance(sqrt(CGAL::squared_distance(points1.first.point, points1.second.point))) < max_slope) {

									K::Point_3 source = Point_3(points1.first.point.x(), points1.first.point.y(), z1);
									K::Vector_3 direction(source, Point_3(points1.second.point.x(), points1.second.point.y(), z2));
									float distance = sqrt(pow(direction.x(),2) + pow(direction.y(),2));

									int label_count[LABELS.size()] = {0};
									bool out = false;

									for(int t = 0; t < distance; t++) {
										K::Point_3 point = source + t/distance * direction;
										auto dsm_z = raster.dsm[int(point.y())][int(point.x())];
										auto dtm_z = raster.dtm[int(point.y())][int(point.x())];
										auto label = raster.land_cover[int(point.y())][int(point.x())];

										if (point.z() > dsm_z + 0.5) {
											out = true;
											break;
										}

										if (point.z() < dtm_z - 0.5) {
											out = true;
											break;
										}

										if (point.z() > dsm_z - 0.5) {
											label_count[label]++;
										}

									}

									if (out) {
										break;
									}

									if (label_count[1] > 0) {
										// bare ground
										break;
									}
									if (label_count[2] > 0) {
										// low vegetation
										break;
									}
									if (selected_label != 3 && label_count[3] > 0) {
										// water
										break;
									}
									if (label_count[4] > 0) {
										// building
										break;
									}
									if (label_count[6] > 0) {
										// parking
										break;
									}
									if (label_count[7] > raster.coord_distance_to_grid_distance(0.5)) {
										// pedestrian
										break;
									}
									if (selected_label == 3 && label_count[8] > 0) {
										// road
										break;
									}
									if (selected_label == 3 && label_count[9] > 0) {
										// railways
										break;
									}
									if (label_count[10] > 0) {
										// swimming pool
										break;
									}

									result.insert(std::pair<skeletonPoint,skeletonPoint>(points1.first, points1.second));

								}
							}
						}
					}

				} else if (path1 == path2) {
					std::map<skeletonPoint, skeletonPoint> nearests1;

					//for each vertex of path1 compute nearest point on path2
					for (auto point1: medial_axes.at(path1)->vertex_handles()) {
						if (point1->is_skeleton()) {
							skeletonPoint source(path1, point1);
							skeletonPoint nearest;
							K::FT min_distance = -1;
							bool nearest_point = false;

							auto he = point1->halfedge_around_vertex_begin();
							do {
								auto point2 = (*he)->opposite()->vertex();
								if (point2->is_skeleton()) {
									if (CGAL::squared_distance(point1->point(), point2->point()) < min_distance || min_distance == -1) {
										min_distance = CGAL::squared_distance(point1->point(), point2->point());
										nearest = skeletonPoint(path2, point2);
									}
								}
							} while (++he != point1->halfedge_around_vertex_begin());

							for (auto point2: medial_axes.at(path2)->vertex_handles()) {
								if (point2->is_skeleton()) {
									if (CGAL::squared_distance(point1->point(), point2->point()) < min_distance || min_distance == -1) {
										min_distance = CGAL::squared_distance(point1->point(), point2->point());
										nearest = skeletonPoint(path2, point2);
										nearest_point = true;
									}
								}
							}
							for (auto edge2: medial_axes.at(path2)->halfedge_handles()) {
								if (edge2->vertex()->id() < edge2->opposite()->vertex()->id()) {
									if(edge2->is_inner_bisector() && edge2->opposite()->is_inner_bisector()) {
										auto target = edge2->vertex()->point();
										auto source = edge2->opposite()->vertex()->point();
										auto proj = K::Line_2(source, target).projection(point1->point());
										if (K::Segment_2(source, target).collinear_has_on(proj)) {
											if (CGAL::squared_distance(point1->point(), proj) < min_distance) {
												min_distance = CGAL::squared_distance(point1->point(), proj);
												nearest = skeletonPoint(path2, edge2, proj);
												nearest_point = true;
											}
										}
									}
								}
							}

							if (nearest_point) {
								nearests1[source] = nearest;
							}
						}
					}

					for (auto points1: nearests1) {
						if (CGAL::squared_distance(points1.first.point, points1.second.point) < pow(raster.coord_distance_to_grid_distance(50),2)) {

							// slope
							float max_slope;
							if (selected_label == 8) {
								// slop for road
								max_slope = 0.3;
							} else {
								// slop for railways and water
								max_slope = 0.1;
							}
							auto z1 = raster.dsm[int(points1.first.point.y())][int(points1.first.point.x())];
							auto z2 = raster.dsm[int(points1.second.point.y())][int(points1.second.point.x())];

							if (abs(z1 - z2) / raster.grid_distance_to_coord_distance(sqrt(CGAL::squared_distance(points1.first.point, points1.second.point))) < max_slope) {

								K::Point_3 source = Point_3(points1.first.point.x(), points1.first.point.y(), z1);
								K::Vector_3 direction(source, Point_3(points1.second.point.x(), points1.second.point.y(), z2));
								float distance = sqrt(pow(direction.x(),2) + pow(direction.y(),2));

								int label_count[LABELS.size()] = {0};
								bool out = false;

								for(int t = 0; t < distance; t++) {
									K::Point_3 point = source + t/distance * direction;
									auto dsm_z = raster.dsm[int(point.y())][int(point.x())];
									auto dtm_z = raster.dtm[int(point.y())][int(point.x())];
									auto label = raster.land_cover[int(point.y())][int(point.x())];

									if (point.z() > dsm_z + 0.5) {
										out = true;
										break;
									}

									if (point.z() < dtm_z - 0.5) {
										out = true;
										break;
									}

									if (point.z() > dsm_z - 0.5) {
										label_count[label]++;
									}

								}

								if (out) {
									break;
								}

								int sum = - label_count[selected_label];
								for (int label = 0; label < LABELS.size(); label++) {
									sum += label_count[label];
								}
								if (sum == 0) {
									// no other class
									break;
								}

								if (label_count[1] > 0) {
									// bare ground
									break;
								}
								if (label_count[2] > 0) {
									// low vegetation
									break;
								}
								if (selected_label != 3 && label_count[3] > 0) {
									// water
									break;
								}
								if (label_count[4] > 0) {
									// building
									break;
								}
								if (label_count[6] > 0) {
									// parking
									break;
								}
								if (label_count[7] > raster.coord_distance_to_grid_distance(0.5)) {
									// pedestrian
									break;
								}
								if (selected_label == 3 && label_count[8] > 0) {
									// road
									break;
								}
								if (selected_label == 3 && label_count[9] > 0) {
									// railways
									break;
								}
								if (label_count[10] > 0) {
									// swimming pool
									break;
								}

								result.insert(std::pair<skeletonPoint,skeletonPoint>(points1.first, points1.second));
							}
						}
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

void add_road(
		std::list<std::pair<K::Vector_2, K::FT>> &roads,
		CGAL::Straight_skeleton_2<K>::Vertex_handle vertex,
		CGAL::Straight_skeleton_2<K>::Vertex_handle previous_vertex) {
	auto he = vertex->halfedge_around_vertex_begin();
	do {
		auto vertex2 = (*he)->opposite()->vertex();
		if (vertex2->is_skeleton() && vertex2 != previous_vertex) {
			if (abs(vertex2->time() - vertex->time()) / sqrt(CGAL::squared_distance(vertex->point(), vertex2->point())) < 0.1) {
				roads.push_back(std::pair<K::Vector_2, K::FT>(K::Vector_2(vertex2->point(), vertex->point()), vertex2->time() + vertex->time()));
			} else {
				add_road(roads, vertex2, vertex);
			}
		}
	} while (++he != vertex->halfedge_around_vertex_begin());

}

std::pair<K::FT, K::FT> road_width (std::pair<skeletonPoint,skeletonPoint> link) {
	K::Vector_2 vector(link.first.point, link.second.point);

	std::list<std::pair<K::Vector_2, K::FT>> roads1;
	std::list<std::pair<K::Vector_2, K::FT>> roads2;

	if (link.first.vertex != nullptr) {
		add_road(roads1, link.first.vertex, nullptr);
	} else {
		add_road(roads1, link.first.halfedge->vertex(), link.first.halfedge->opposite()->vertex());
		add_road(roads1, link.first.halfedge->opposite()->vertex(), link.first.halfedge->vertex());
		if (abs(link.first.halfedge->opposite()->vertex()->time() - link.first.halfedge->vertex()->time()) / sqrt(CGAL::squared_distance(link.first.halfedge->vertex()->point(), link.first.halfedge->opposite()->vertex()->point())) < 0.1) {
			roads1.push_back(std::pair<K::Vector_2, K::FT>(K::Vector_2(link.first.halfedge->opposite()->vertex()->point(), link.first.halfedge->vertex()->point()), link.first.halfedge->opposite()->vertex()->time() + link.first.halfedge->vertex()->time()));
		}
	}

	if (link.second.vertex != nullptr) {
		add_road(roads2, link.second.vertex, nullptr);
	} else {
		add_road(roads2, link.second.halfedge->vertex(), link.second.halfedge->opposite()->vertex());
		add_road(roads2, link.second.halfedge->opposite()->vertex(), link.second.halfedge->vertex());
		if (abs(link.second.halfedge->opposite()->vertex()->time() - link.second.halfedge->vertex()->time()) / sqrt(CGAL::squared_distance(link.second.halfedge->vertex()->point(), link.second.halfedge->opposite()->vertex()->point())) < 0.1) {
			roads2.push_back(std::pair<K::Vector_2, K::FT>(K::Vector_2(link.second.halfedge->opposite()->vertex()->point(), link.second.halfedge->vertex()->point()), link.second.halfedge->opposite()->vertex()->time() + link.second.halfedge->vertex()->time()));
		}
	}

	K::FT road_width1 = 0;
	K::FT sum1 = 0;
	for (auto road: roads1) {
		auto angle = abs(CGAL::scalar_product(road.first, vector) / sqrt(road.first.squared_length()));
		road_width1 += road.second * angle;
		sum1 += angle;
	}

	K::FT road_width2 = 0;
	K::FT sum2 = 0;
	for (auto road: roads2) {
		auto angle = abs(CGAL::scalar_product(road.first, vector) / sqrt(road.first.squared_length()));
		road_width2 += road.second * angle;
		sum2 += angle;
	}

	return std::pair<K::FT, K::FT>(road_width1/sum1, road_width2/sum2);

}


void bridge (std::pair<skeletonPoint,skeletonPoint> link, const Surface_mesh &mesh, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
	Filtered_graph filtered_sm1(mesh, link.first.path, path);
	Filtered_graph filtered_sm2(mesh, link.second.path, path);

	typedef CGAL::dynamic_vertex_property_t<Point_2>               Point_2_property;
	auto projection_pmap = CGAL::get(Point_2_property(), mesh);
	
	for(auto v : CGAL::vertices(mesh)) {
		const Point_3& p = mesh.point(v);
		put(projection_pmap, v, Point_2(p.x(), p.y()));
	}

	auto location1 = PMP::locate(link.first.point, filtered_sm1, CGAL::parameters::vertex_point_map(projection_pmap));
	auto location2 = PMP::locate(link.second.point, filtered_sm2, CGAL::parameters::vertex_point_map(projection_pmap));
	
	K::Vector_2 link_vector(link.first.point, link.second.point);
	auto left = link_vector.perpendicular(CGAL::CLOCKWISE);

	left /= sqrt(left.squared_length());

	auto width = road_width(link);

	Point_2 left1 = link.first.point + left*width.first/2;
	Point_2 right1 = link.first.point - left*width.first/2;
	Point_2 left2 = link.second.point + left*width.second/2;
	Point_2 right2 = link.second.point - left*width.second/2;

	auto point1 = PMP::construct_point(location1, mesh);
	auto point2 = PMP::construct_point(location2, mesh);

	{ // Skeleton
		Surface_mesh skeleton;

		auto v1 = skeleton.add_vertex(point1);
		auto v2 = skeleton.add_vertex(point2);
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left1.x(), left1.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(left2.x(), left2.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(right1.x(), right1.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(right2.x(), right2.y(), point2.z()));
		skeleton.add_edge(v1, v2);

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

		std::cout << location1.first << "\n";
		std::cout << ((int) label[location1.first]) << "\n";
		std::cout << location2.first << "\n";
		std::cout << ((int) label[location2.first]) << "\n";

		std::stringstream skeleton_name;
		skeleton_name << "bridge_" << ((int) label[location1.first]) << "_" << link.first.path << "_" << link.second.path << ".ply";
		std::ofstream mesh_ofile (skeleton_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, skeleton);
		mesh_ofile.close();
	}

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
		{"mesh", required_argument, NULL, 'M'},
		{"terrain_mesh", required_argument, NULL, 'T'},
		{NULL, 0, 0, '\0'}
	};

	char *DSM = NULL;
	char *DTM = NULL;
	char *land_use_map = NULL;
	char *LOD0 = NULL;
	char *orthophoto = NULL;
	char *MESH = NULL;
	char *TERRAIN_MESH = NULL;

	while ((opt = getopt_long(argc, argv, "hs:t:l:0:i:M:T:", options, NULL)) != -1) {
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
				std::cout << " -M, --mesh=/file/path.ply          mesh as PLY file." << std::endl;
				std::cout << " -T, --terrain_mesh=/file/path.ply  terrain mesh as PLY file." << std::endl;
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
			case 'M':
				MESH = optarg;
				break;
			case 'T':
				TERRAIN_MESH = optarg;
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

	const Raster raster(DSM, DTM, land_use_map);

	Surface_mesh terrain_mesh, mesh;
	if (MESH == NULL || TERRAIN_MESH == NULL) {
		std::tie(terrain_mesh, mesh) = compute_meshes(raster);

		std::ofstream mesh_ofile ("save_mesh.ply", std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ofile);
		CGAL::IO::write_PLY (mesh_ofile, mesh);
		mesh_ofile.close();

		mesh_ofile = std::ofstream("save_terrain_mesh.ply", std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ofile);
		CGAL::IO::write_PLY (mesh_ofile, terrain_mesh);
		mesh_ofile.close();

	} else {
		std::ifstream mesh_ifile (MESH, std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ifile);
		CGAL::IO::read_PLY (mesh_ifile, mesh);
		mesh_ifile.close();

		mesh_ifile = std::ifstream(TERRAIN_MESH, std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ifile);
		CGAL::IO::read_PLY (mesh_ifile, terrain_mesh);
		mesh_ifile.close();
		std::cout << "Mesh and terrain mesh load" << std::endl;
	}

	add_label(raster, mesh);
	change_vertical_faces(mesh, raster);
	save_mesh(mesh, raster, "final-mesh-without-facade.ply");

	std::vector<std::list<Surface_mesh::Face_index>> paths = compute_path(mesh);
	save_mesh(mesh, raster, "final-mesh-with-path.ply");

	std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> medial_axes = compute_medial_axes(mesh, paths, raster);

	std::set<std::pair<skeletonPoint,skeletonPoint>> links = link_paths(mesh, paths, medial_axes, raster);

	for (auto link: links) {
		bridge(link, mesh, raster);
	}

	return EXIT_SUCCESS;
}
