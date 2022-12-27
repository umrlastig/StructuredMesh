#ifndef RASTER_H_
#define RASTER_H_

#include "header.hpp"

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

		static std::pair<int,int> grid_conversion(int P, int L, double grid_to_crs[6], OGRCoordinateTransformation *crs_to_other_crs, double other_grid_to_other_crs[6]);

		void align_land_cover_and_dsm();

	public:
		void coord_to_grid(double x, double y, float& P, float& L) const;

		void grid_to_coord(int P, int L, double& x, double& y) const;

		void grid_to_coord(float P, float L, double& x, double& y) const;

		float coord_distance_to_grid_distance(double d) const;

		double grid_distance_to_coord_distance(float d) const;

		OGRSpatialReference get_crs() const;

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

		Raster(char *dsm_path, char *dtm_path, char *land_cover_path);
};

#endif  /* !RASTER_H_ */
