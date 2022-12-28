#ifndef BRIDGE_H_
#define BRIDGE_H_

#include "header.hpp"
#include "raster.hpp"

typedef std::pair<skeletonPoint,skeletonPoint> pathLink;

std::set<pathLink> link_paths(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> &medial_axes, const Raster &raster);

struct pathBridge {
	pathLink link;
	unsigned char label;
	int N;
	float *xl;
	float *xr;
	float *z_segment;
	float cost;

	pathBridge(pathLink link);
	~pathBridge();

};

pathBridge bridge (pathLink link, const Surface_mesh &mesh, const Raster &raster);

#endif  /* !BRIDGE_H_ */
