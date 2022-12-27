#ifndef HEADER_H_
#define HEADER_H_

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
#include <limits>

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
#include <CGAL/Delaunay_triangulation_2.h>
#include <Eigen/SparseQR>

#include "label.hpp"

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

#endif  /* !HEADER_H_ */