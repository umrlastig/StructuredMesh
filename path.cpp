#include "header.hpp"
#include "raster.hpp"

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Straight_skeleton_converter_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Arr_segment_traits_2<Exact_predicates_kernel> Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                       Arrangement_2;

void save_mesh(const Surface_mesh &mesh, const Raster &raster, const char *filename);

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


std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> compute_path_polygon(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> path_polygon;

	typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
	for (int i = 0; i < paths.size(); i++) {

		Filtered_graph filtered_sm(mesh, i, path);

		if (filtered_sm.number_of_faces() > 3) {

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

				path_polygon[i] = polygon((*(arr.unbounded_face()->holes_begin()))->twin()->face());

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

				Surface_mesh part_mesh;
				CGAL::copy_face_graph(filtered_sm, part_mesh);
				std::stringstream name;
				name << "part_mesh_" << lab << "_" << i << ".ply";
				save_mesh(part_mesh, raster, name.str().c_str());

			}
		}
	}

	return path_polygon;

}


std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> compute_medial_axes(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const Raster &raster) {
	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
	assert(has_label);

	std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> medial_axes;

	for (int i = 0; i < paths.size(); i++) {

		if (path_polygon.count(i) > 0) {

			int lab = label[paths.at(i).front()];

			auto poly = path_polygon.at(i);

			poly = CGAL::Polyline_simplification_2::simplify(
				poly,
				CGAL::Polyline_simplification_2::Squared_distance_cost(),
				CGAL::Polyline_simplification_2::Stop_above_cost_threshold(pow(raster.coord_distance_to_grid_distance(1.5),2))
			);

			auto iss = CGAL::create_interior_straight_skeleton_2(poly, Exact_predicates_kernel());
			medial_axes[i] = CGAL::convert_straight_skeleton_2<CGAL::Straight_skeleton_2<K>>(*iss);

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

		}
	}

	return medial_axes;

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
	if (roads1.size() > 0) {
		K::FT sum1 = 0;
		for (auto road: roads1) {
			auto angle = abs(CGAL::scalar_product(road.first, vector) / sqrt(road.first.squared_length())) + 0.01;
			road_width1 += road.second * angle;
			sum1 += angle;
		}
		road_width1 /= sum1;
	} else if (link.first.vertex != nullptr) {
		road_width1 = 2*link.first.vertex->time();
	} else {
		road_width1 = link.first.halfedge->vertex()->time() + link.first.halfedge->opposite()->vertex()->time();
	}

	K::FT road_width2 = 0;
	if (roads2.size() > 0) {
		K::FT sum2 = 0;
		for (auto road: roads2) {
			auto angle = abs(CGAL::scalar_product(road.first, vector) / sqrt(road.first.squared_length())) + 0.01;
			road_width2 += road.second * angle;
			sum2 += angle;
		}
		road_width2 /= sum2;
	} else if (link.second.vertex != nullptr) {
		road_width2 = 2*link.second.vertex->time();
	} else {
		road_width2 = link.second.halfedge->vertex()->time() + link.second.halfedge->opposite()->vertex()->time();
	}

	return std::pair<K::FT, K::FT>(road_width1, road_width2);

}


std::pair<Surface_mesh::Face_index, Point_2> point_on_path_border(const Surface_mesh &mesh, Surface_mesh::Face_index face, std::list<K::Segment_2> segments) {
	assert(face != Surface_mesh::null_face());

	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = mesh.property_map<Surface_mesh::Face_index, int>("path");
	assert(has_path);

	typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
	Filtered_graph filtered_mesh(mesh, path[face], path);

	// find first intersecting edge
	K::Segment_2 segment;
	auto he = CGAL::halfedge(face, filtered_mesh);
	while (!segments.empty()) {
		segment = segments.front();
		segments.pop_front();
		do {
			auto source = mesh.point(CGAL::source(he, filtered_mesh));
			auto target = mesh.point(CGAL::target(he, filtered_mesh));
			if (CGAL::do_intersect(segment, K::Segment_2(Point_2(source.x(), source.y()), Point_2(target.x(), target.y())))) {
				goto goto1_point_on_path_border;
			}
			he = CGAL::next(he, mesh);
		} while (he != CGAL::halfedge(face, filtered_mesh));
	}
	return std::make_pair(face, segment.target());
	goto1_point_on_path_border:
	face = CGAL::face(CGAL::opposite(he, filtered_mesh), filtered_mesh);
	while(face != Surface_mesh::null_face()) {
		he = CGAL::next(CGAL::opposite(he, filtered_mesh), filtered_mesh);
		auto source = mesh.point(CGAL::source(he, filtered_mesh));
		auto target = mesh.point(CGAL::target(he, filtered_mesh));
		if (CGAL::do_intersect(segment, K::Segment_2(Point_2(source.x(), source.y()), Point_2(target.x(), target.y())))) {
			face = CGAL::face(CGAL::opposite(he, filtered_mesh), filtered_mesh);
			continue;
		}
		he = CGAL::next(he, filtered_mesh);
		source = mesh.point(CGAL::source(he, filtered_mesh));
		target = mesh.point(CGAL::target(he, filtered_mesh));
		if (CGAL::do_intersect(segment, K::Segment_2(Point_2(source.x(), source.y()), Point_2(target.x(), target.y())))) {
			face = CGAL::face(CGAL::opposite(he, filtered_mesh), filtered_mesh);
			continue;
		}

		he = CGAL::halfedge(face, filtered_mesh);
		while (!segments.empty()) {
			segment = segments.front();
			segments.pop_front();
			do {
				auto source = mesh.point(CGAL::source(he, filtered_mesh));
				auto target = mesh.point(CGAL::target(he, filtered_mesh));
				if (CGAL::do_intersect(segment, K::Segment_2(Point_2(source.x(), source.y()), Point_2(target.x(), target.y())))) {
					goto goto2_point_on_path_border;
				}
				he = CGAL::next(he, mesh);
			} while (he != CGAL::halfedge(face, filtered_mesh));
		}
		return std::make_pair(face, segment.target());
		goto2_point_on_path_border:
		face = CGAL::face(CGAL::opposite(he, filtered_mesh), filtered_mesh);
	}
	face = CGAL::face(he, filtered_mesh);
	auto source = mesh.point(CGAL::source(he, filtered_mesh));
	auto target = mesh.point(CGAL::target(he, filtered_mesh));
	auto result = CGAL::intersection(segment, K::Segment_2(Point_2(source.x(), source.y()), Point_2(target.x(), target.y())));
	assert(result);
	if (const K::Segment_2* s = boost::get<K::Segment_2>(&*result)) {
		if (CGAL::squared_distance(Point_2(source.x(), source.y()), segment.source()) < CGAL::squared_distance(Point_2(target.x(), target.y()), segment.source())) {
			return std::make_pair(face, Point_2(source.x(), source.y()));
		} else {
			return std::make_pair(face, Point_2(target.x(), target.y()));
		}
	} else {
		const Point_2* p = boost::get<Point_2 >(&*result);
		return std::make_pair(face, *p);
	}
}