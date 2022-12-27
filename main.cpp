#include "header.hpp"
#include "raster.hpp"

#include <CGAL/Polygon_mesh_processing/locate.h>
#include <Eigen/SparseQR>

namespace PMP = CGAL::Polygon_mesh_processing;

void save_mesh(const Surface_mesh &mesh, const Raster &raster, const char *filename);

std::tuple<Surface_mesh, Surface_mesh> compute_meshes(const Raster &raster);

void add_label(const Raster &raster, Surface_mesh &mesh);

void change_vertical_faces(Surface_mesh &mesh, const Raster &raster);

std::vector<std::list<Surface_mesh::Face_index>> compute_path(Surface_mesh &mesh);

std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> compute_path_polygon(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const Raster &raster);


std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> compute_medial_axes(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const Raster &raster);


std::set<std::pair<skeletonPoint,skeletonPoint>> link_paths(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> &medial_axes, const Raster &raster);

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


void bridge (std::pair<skeletonPoint,skeletonPoint> link, const Surface_mesh &mesh, const Raster &raster) {
	std::cout << "Bridge " << link.first.path << " (" << link.first.point << ") -> " << link.second.path << " (" << link.second.point << ")\n";


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
	auto point1 = PMP::construct_point(location1, mesh);
	auto point2 = PMP::construct_point(location2, mesh);

	unsigned char bridge_label = label[location1.first];

	K::Vector_2 link_vector(link.first.point, link.second.point);
	float length = sqrt(link_vector.squared_length());
	auto l = link_vector / length;
	auto n = l.perpendicular(CGAL::COUNTERCLOCKWISE);

	auto width = road_width(link);

	Point_2 left1 = link.first.point - n*width.first/2;
	Point_2 right1 = link.first.point + n*width.first/2;
	Point_2 left2 = link.second.point - n*width.second/2;
	Point_2 right2 = link.second.point + n*width.second/2;

	auto left1_border = point_on_path_border(mesh, location1.first, {K::Segment_2(link.first.point, left1)});
	auto right1_border = point_on_path_border(mesh, location1.first, {K::Segment_2(link.first.point, right1)});
	auto left2_border = point_on_path_border(mesh, location2.first, {K::Segment_2(link.second.point, left2)});
	auto right2_border = point_on_path_border(mesh, location2.first, {K::Segment_2(link.second.point, right2)});

	float dl0 = sqrt(CGAL::squared_distance(link.first.point, left1_border.second));
	float dr0 = sqrt(CGAL::squared_distance(link.first.point, right1_border.second));
	float dlN = sqrt(CGAL::squared_distance(link.second.point, left2_border.second));
	float drN = sqrt(CGAL::squared_distance(link.second.point, right2_border.second));

	int N = int(sqrt(link_vector.squared_length ()));
	float xl[N+1];
	float xr[N+1];
	float z_segment[N+1];
	for (int i = 0; i <= N; i++) {
		xl[i] = dl0 + ((float) i)/N*(dlN-dl0);
		xr[i] = dr0 + ((float) i)/N*(drN-dr0);
		z_segment[i] = point1.z() + (point2.z() - point1.z())*((float) i)/N;
	}

	float tunnel_height = 3; // in meter

	for(int loop = 0; loop < 10; loop++) {

		std::cout << "Solve surface\n";

		// Surface solver
		{
			std::vector<Eigen::Triplet<float>> tripletList;
			tripletList.reserve(10*N+2);
			std::vector<float> temp_b;
			temp_b.reserve(9*N+2);

			// regularity of the surface
			float alpha = 10;
			for (int i = 0; i < N; i++) {
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i, alpha));
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i+1, -alpha));
				temp_b.push_back(0);
			}

			// attachment to DSM data
			float beta = 1;
			for (int i = 0; i <= N; i++) {
				int j = std::max({0, ((int) (- xl[i]))});
				while(j <= xr[i]) {
					auto p = link.first.point + ((float) i)/N*link_vector + j * n;
					if (p.y() >= 0 && p.y() < raster.ySize && p.x() >= 0 && p.y() < raster.xSize) {
						if (z_segment[i] > raster.dsm[p.y()][p.x()] - tunnel_height / 2) {
							tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i, beta));
							temp_b.push_back(raster.dsm[p.y()][p.x()]);
						} else if (z_segment[i] > raster.dsm[p.y()][p.x()] - tunnel_height) {
							tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i, beta));
							temp_b.push_back(raster.dsm[p.y()][p.x()] - tunnel_height);
						}
					}
					j++;
				}
				j = std::max({1, ((int) (- xr[i]))});
				while(j <= xl[i]) {
					auto p = link.first.point + ((float) i)/N*link_vector - j * n;
					if (p.y() >= 0 && p.y() < raster.ySize && p.x() >= 0 && p.y() < raster.xSize) {
						if (z_segment[i] > raster.dsm[p.y()][p.x()] - tunnel_height / 2) {
							tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i, beta));
							temp_b.push_back(raster.dsm[p.y()][p.x()]);
						} else if (z_segment[i] > raster.dsm[p.y()][p.x()] - tunnel_height) {
							tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i, beta));
							temp_b.push_back(raster.dsm[p.y()][p.x()] - tunnel_height);
						}
					}
					j++;
				}
			}

			// border
			float zeta = 10;
			tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), 0, zeta));
			temp_b.push_back(zeta * point1.z());
			tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), N, zeta));
			temp_b.push_back(zeta * point2.z());

			// solving
			Eigen::SparseMatrix<float> A(temp_b.size(), N+1);
			A.setFromTriplets(tripletList.begin(), tripletList.end());
			Eigen::VectorXf b = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(temp_b.data(), temp_b.size());

			Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>> solver;
			solver.compute(A);
			if(solver.info()!=Eigen::Success) {
				std::cout << "decomposition failed\n";
				return;
			}
			Eigen::VectorXf x = solver.solve(b);
			if(solver.info()!=Eigen::Success) {
				std::cout << "solving failed\n";
				return;
			}
			for (int i = 0; i <= N; i++) {
				z_segment[i] = x[i];
			}
		}

		std::cout << "Solve contour\n";

		// Contour solver
		{
			std::vector<Eigen::Triplet<float>> tripletList;
			tripletList.reserve(2*N + N+1 + 4 + 4);
			std::vector<float> temp_b;
			temp_b.reserve(2*N + N+1 + 4 + 4);

			//regularity of the contour
			float gamma = 5;
			for (int j = 0; j < N; j++) {
				// x^l_j − x^l_{j+1}
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), j,gamma));
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), j+1,-gamma));
				temp_b.push_back(0);
			}
			for (int j = N+1; j < 2*N+1; j++) {
				// x^r_{j-N-1} − x^r_{j+1-N-1}
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), j,gamma));
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), j+1,-gamma));
				temp_b.push_back(0);
			}

			//width of the reconstructed surface
			float delta = 10;
			for (int j = 0; j <= N; j++) {
				// x^l_j − x^l_{j+1}
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), j,delta)); //x^l_j
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), j+N+1,delta)); //x^r_j
				temp_b.push_back(delta * (width.first + j*(width.second - width.first)/N));
			}

			//centering of the surface on the link vertices
			float epsilon = 1;
			tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), 0,epsilon)); //x^l_0
			tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), N+1,-epsilon)); //x^r_0
			temp_b.push_back(0);
			tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), N,epsilon)); //x^l_N
			tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), 2*N+1,-epsilon)); //x^r_N
			temp_b.push_back(0);

			//keeping surface where cost is low
			float theta = 1;
			float iota  = 10;
			// left contour
			for (int i = 0; i <= N; i++) {
				float j = xl[i];
				float j_min = j;
				float v_min;
				auto p = link.first.point + ((float) i)/N*link_vector + j * n;
				if (p.y() >= 0 && p.y() < raster.ySize && p.x() >= 0 && p.y() < raster.xSize) {
					if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height) {
						continue;
					} else if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height/2) {
						v_min = pow(z_segment[i] - raster.dsm[p.y()][p.x()] + tunnel_height, 2);
					} else {
						v_min = pow(z_segment[i] - raster.dsm[p.y()][p.x()], 2);
						if (raster.land_cover[p.y()][p.x()] != bridge_label) {
							if ((raster.land_cover[p.y()][p.x()] != 8 && raster.land_cover[p.y()][p.x()] != 9) || (bridge_label != 8 && bridge_label != 9)) {
								v_min += theta;
							}
						}
					}
				} else {
					continue;
				}

				j--;
				while(j > -xr[i]) {
					auto p = link.first.point + ((float) i)/N*link_vector + j * n;
					if (p.y() >= 0 && p.y() < raster.ySize && p.x() >= 0 && p.y() < raster.xSize) {
						if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height) {
							j_min = j;
							break;
						} else if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height/2) {
							float v = pow(z_segment[i] - raster.dsm[p.y()][p.x()] + tunnel_height, 2);
							if (v < v_min) {
								j_min = j;
								v_min = v;
							}
						} else {
							float v = pow(z_segment[i] - raster.dsm[p.y()][p.x()], 2);
							if (raster.land_cover[p.y()][p.x()] != bridge_label) {
								if ((raster.land_cover[p.y()][p.x()] != 8 && raster.land_cover[p.y()][p.x()] != 9) || (bridge_label != 8 && bridge_label != 9)) {
									v += theta;
								}
							}
							if (v < v_min) {
								j_min = j;
								v_min = v;
							}
						}
					} else {
						break;
					}
					j--;
				}
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i, iota));
				temp_b.push_back(iota * j_min);
			}
			// right contour
			for (int i = 0; i <= N; i++) {
				float j = xr[i];
				float j_min = j;
				float v_min;
				auto p = link.first.point + ((float) i)/N*link_vector + j * n;
				if (p.y() >= 0 && p.y() < raster.ySize && p.x() >= 0 && p.y() < raster.xSize) {
					if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height) {
						continue;
					} else if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height/2) {
						v_min = pow(z_segment[i] - raster.dsm[p.y()][p.x()] + tunnel_height, 2);
					} else {
						v_min = pow(z_segment[i] - raster.dsm[p.y()][p.x()], 2);
						if (raster.land_cover[p.y()][p.x()] != bridge_label) {
							if ((raster.land_cover[p.y()][p.x()] != 8 && raster.land_cover[p.y()][p.x()] != 9) || (bridge_label != 8 && bridge_label != 9)) {
								v_min += theta;
							}
						}
					}
				} else {
					continue;
				}

				j--;
				while(j > -xl[i]) {
					auto p = link.first.point + ((float) i)/N*link_vector + j * n;
					if (p.y() >= 0 && p.y() < raster.ySize && p.x() >= 0 && p.y() < raster.xSize) {
						if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height) {
							j_min = j;
							break;
						} else if (z_segment[i] <= raster.dsm[p.y()][p.x()] - tunnel_height/2) {
							float v = pow(z_segment[i] - raster.dsm[p.y()][p.x()] + tunnel_height, 2);
							if (v < v_min) {
								j_min = j;
								v_min = v;
							}
						} else {
							float v = pow(z_segment[i] - raster.dsm[p.y()][p.x()], 2);
							if (raster.land_cover[p.y()][p.x()] != bridge_label) {
								if ((raster.land_cover[p.y()][p.x()] != 8 && raster.land_cover[p.y()][p.x()] != 9) || (bridge_label != 8 && bridge_label != 9)) {
									v += theta;
								}
							}
							if (v < v_min) {
								j_min = j;
								v_min = v;
							}
						}
					} else {
						break;
					}
					j--;
				}
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), i + N+1, iota));
				temp_b.push_back(iota * j_min);
			}

			// fix border
			float eta = 100;
			if (xl[0] > dl0) {
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), 0, eta)); //x^l_0
				temp_b.push_back(eta*dl0);
			}
			if (xl[N] > dlN) {
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), N, eta)); //x^l_N
				temp_b.push_back(eta*dlN);
			}
			if (xr[0] > dr0) {
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), N+1, eta)); //x^r_0
				temp_b.push_back(eta*dr0);
			}
			if (xr[N] > drN) {
				tripletList.push_back(Eigen::Triplet<float>(temp_b.size(), 2*N+1, eta)); //x^r_N
				temp_b.push_back(eta*drN);
			}

			// solving
			Eigen::SparseMatrix<float> A(temp_b.size(),2*N+2);
			A.setFromTriplets(tripletList.begin(), tripletList.end());
			Eigen::VectorXf b = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(temp_b.data(), temp_b.size());
			
			Eigen::SparseQR<Eigen::SparseMatrix<float>, Eigen::COLAMDOrdering<int>> solver;
			solver.compute(A);
			if(solver.info()!=Eigen::Success) {
				std::cout << "decomposition failed\n";
				return;
			}
			Eigen::VectorXf x = solver.solve(b);
			if(solver.info()!=Eigen::Success) {
				std::cout << "solving failed\n";
				return;
			}
			for (int j = 0; j < N; j++) {
				// x^l_j
				xl[j] = x[j];
				// x^r_j
				xr[j] = x[j+N+1];
			}

		}
	}
	
	if (xl[0] > dl0) xl[0] = dl0;
	if (xl[N] > dlN) xl[N] = dlN;
	if (xr[0] > dr0) xr[0] = dr0;
	if (xr[N] > drN) xr[N] = drN;

	std::cout << "Save\n";

	{ // Surface

		Surface_mesh bridge_mesh;
		Surface_mesh::Vertex_index Xl[N+1];
		Surface_mesh::Vertex_index Xr[N+1];

		// Add points
		for (int i = 0; i <= N; i++) {
			auto p1 = link.first.point + ((float) i)/N*link_vector - xl[i] * n;
			auto p2 = link.first.point + ((float) i)/N*link_vector + xr[i] * n;
			Xl[i] = bridge_mesh.add_vertex(Point_3(p1.x(), p1.y(), z_segment[i]));
			Xr[i] = bridge_mesh.add_vertex(Point_3(p2.x(), p2.y(), z_segment[i]));
		}

		// Add faces
		for (int i = 0; i < N; i++) {
			bridge_mesh.add_face(Xl[i], Xr[i], Xr[i+1]);
			bridge_mesh.add_face(Xl[i], Xr[i+1], Xl[i+1]);
		}

		double min_x, min_y;
		raster.grid_to_coord(0, 0, min_x, min_y);

		for(auto vertex : bridge_mesh.vertices()) {
			auto point = bridge_mesh.point(vertex);
			double x, y;
			raster.grid_to_coord((float) point.x(), (float) point.y(), x, y);
			bridge_mesh.point(vertex) = Point_3(x-min_x, y-min_y, point.z());
		}

		std::stringstream bridge_mesh_name;
		bridge_mesh_name << "bridge_mesh_" << ((int) bridge_label) << "_" << link.first.path << "_" << link.second.path << "_" << link.first.point << "_" << link.second.point << ".ply";
		std::ofstream mesh_ofile (bridge_mesh_name.str().c_str());
		CGAL::IO::write_PLY (mesh_ofile, bridge_mesh);
		mesh_ofile.close();

		std::cout << bridge_mesh_name.str() << " saved\n";

	}


	{ // Skeleton
		Surface_mesh skeleton;

		auto v1 = skeleton.add_vertex(point1);
		auto v2 = skeleton.add_vertex(point2);
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left1_border.second.x(), left1_border.second.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(left2_border.second.x(), left2_border.second.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(right1_border.second.x(), right1_border.second.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(right2_border.second.x(), right2_border.second.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left1.x(), left1.y(), point1.z()));
		v2 = skeleton.add_vertex(Point_3(right1.x(), right1.y(), point1.z()));
		skeleton.add_edge(v1, v2);

		v1 = skeleton.add_vertex(Point_3(left2.x(), left2.y(), point2.z()));
		v2 = skeleton.add_vertex(Point_3(right2.x(), right2.y(), point2.z()));
		skeleton.add_edge(v1, v2);

		/*auto init = skeleton.add_vertex(Point_3(Xr[0].x(), Xr[0].y(), point1.z()));
		auto temp1 = CGAL::SM_Vertex_index(init);
		for (int i = 1; i <= N; i++) {
			auto temp2 = skeleton.add_vertex(Point_3(Xr[i].x(), Xr[i].y(), point1.z()));
			skeleton.add_edge(temp1, temp2);
			temp1 = temp2;
		}
		for (int i = N; i >= 0; i--) {
			auto temp2 = skeleton.add_vertex(Point_3(Xl[i].x(), Xl[i].y(), point1.z()));
			skeleton.add_edge(temp1, temp2);
			temp1 = temp2;
		}
		skeleton.add_edge(temp1, init);*/

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
		skeleton_name << "bridge_" << ((int) bridge_label) << "_" << link.first.path << "_" << link.second.path << "_" << link.first.point << "_" << link.second.point << ".ply";
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

	std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> path_polygon = compute_path_polygon(mesh, paths, raster);
	std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> medial_axes = compute_medial_axes(mesh, paths, path_polygon, raster);

	std::set<std::pair<skeletonPoint,skeletonPoint>> links = link_paths(mesh, paths, path_polygon, medial_axes, raster);

	for (auto link: links) {
		if (link.first.path == 2 && link.second.path == 75 || link.first.path == 53 && link.second.path == 94) {
		bridge(link, mesh, raster);
		}
	}

	return EXIT_SUCCESS;
}
