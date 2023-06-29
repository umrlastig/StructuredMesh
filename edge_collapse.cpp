#include "header.hpp"
#include "raster.hpp"

#include <chrono>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

namespace SMS = CGAL::Surface_mesh_simplification;

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

struct My_visitor : SMS::Edge_collapse_visitor_base<Surface_mesh> {
	private:
		int i_collecte = 0;
		Surface_mesh &mesh;
		const Surface_mesh_info &mesh_info;
		std::chrono::time_point<std::chrono::system_clock> start_collecte;
		std::chrono::time_point<std::chrono::system_clock> start_collapse;
		bool output[5] = {false};

	public:
		My_visitor(const Surface_mesh &mesh, const Surface_mesh_info &mesh_info) : mesh(mesh), mesh_info(mesh_info) {}

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
					mesh_info.save_mesh(mesh, "mesh-1e-4.ply");
				} else if(*cost > 0 && !output[1]) {
					output[1] = true;
					mesh_info.save_mesh(mesh, "mesh-0.ply");
				}
				if(current_edge_count <= 100000 && !output[2]) {
					output[2] = true;
					mesh_info.save_mesh(mesh, "mesh-100000.ply");
				} else if(current_edge_count <= 10000 && !output[3]) {
					output[3] = true;
					mesh_info.save_mesh(mesh, "mesh-10000.ply");
				} else if(current_edge_count <= 5000 && !output[4]) {
					output[4] = true;
					mesh_info.save_mesh(mesh, "mesh-5000.ply");
				} else if(current_edge_count <= 1000000 && !output[5]) {
					output[5] = true;
					mesh_info.save_mesh(mesh, "mesh-1000000.ply");
				}
			}

		}

};

std::tuple<Surface_mesh, Surface_mesh> compute_meshes(const Raster &raster, const Surface_mesh_info &mesh_info) {

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

	mesh_info.save_mesh(terrain_mesh, "initial-terrain-mesh.ply");

	SMS::edge_collapse(terrain_mesh, Cost_stop_predicate(10));
	std::cout << "Terrain mesh simplified" << std::endl;

	mesh_info.save_mesh(terrain_mesh, "terrain-mesh.ply");

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

	mesh_info.save_mesh(mesh, "initial-mesh.ply");

	Cost_stop_predicate stop(10);
	//SMS::Count_stop_predicate<Surface_mesh> stop(1000);
	Custom_cost cf(raster);
	Custom_placement pf(raster);
	int r = SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).get_placement(pf).visitor(My_visitor(mesh, mesh_info)));
	std::cout << "\rMesh simplified                                               " << std::endl;

	mesh_info.save_mesh(mesh, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, mesh);
}