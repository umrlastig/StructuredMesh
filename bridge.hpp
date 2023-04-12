#ifndef BRIDGE_H_
#define BRIDGE_H_

#include "header.hpp"
#include "raster.hpp"

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <vector>
#include <map>

typedef std::pair<skeletonPoint,skeletonPoint> pathLink;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>        AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                AABB_tree;

std::set<pathLink> link_paths(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> &medial_axes, const Raster &raster);

struct pathBridge {
	pathLink link;
	unsigned char label;
	int N;
	double *xl;
	double *xr;
	double *z_segment;
	float cost;
	CGAL::Surface_mesh<CGAL::Exact_predicates_exact_constructions_kernel::Point_3> crossing_surface;

	pathBridge(pathLink link);
	pathBridge(const pathBridge& other);
	~pathBridge();

};

pathBridge bridge (pathLink link, const Surface_mesh &mesh, const AABB_tree &tree, const Raster &raster);

void close_surface_mesh(Surface_mesh &mesh);

AABB_tree index_surface_mesh(Surface_mesh &mesh);

void add_bridge_to_mesh(Surface_mesh &mesh, const std::vector<pathBridge> &bridges, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const Raster &raster);

#endif  /* !BRIDGE_H_ */
