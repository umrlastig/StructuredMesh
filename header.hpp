#ifndef HEADER_H_
#define HEADER_H_

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <list>
#include <vector>
#include <set>
#include <algorithm>
#include <cmath>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Point_set_3.h>

#include <ogr_spatialref.h>

#include "label.hpp"

typedef CGAL::Simple_cartesian<float>                       K;
typedef K::Point_2                                          Point_2;
typedef K::Point_3                                          Point_3;
typedef CGAL::Surface_mesh<Point_3>                         Surface_mesh;
typedef CGAL::Polygon_2<K>                                  Polygon;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Exact_predicates_kernel;

typedef Exact_predicates_kernel Point_set_kernel;
typedef CGAL::Point_set_3<Point_set_kernel::Point_3> Point_set;

/// A point on a skeleton
struct skeletonPoint {
	/// The skeleton id
	int path;
	/// The point coordinate
	Point_2 point;
	/// A handle to the vertex if the point is a vertex
	CGAL::Straight_skeleton_2<K>::Vertex_handle vertex;
	/// Or a handle to the edge otherwise
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
		if (l.path == r.path) return l.point < r.point;
		return l.path < r.path;
	}

	friend bool operator==(const skeletonPoint& l, const skeletonPoint& r) {
		return (l.path == r.path && l.point == r.point);
	}
};

// contains information on mesh location
struct Surface_mesh_info {
	const OGRSpatialReference crs;
	const double x_0;
	const double y_0;

	Surface_mesh_info();
	Surface_mesh_info(OGRSpatialReference crs, double x_0, double y_0);

	void save_mesh (const Surface_mesh &mesh, const char *filename) const;
};

#endif  /* !HEADER_H_ */
