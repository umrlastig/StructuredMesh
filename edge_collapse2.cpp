#include "header.hpp"
#include "raster.hpp"

#include <chrono>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_profile.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_matrix.h>
#include <CGAL/Eigen_vector.h>
#include <CGAL/Eigen_svd.h>
#endif

typedef CGAL::Eigen_svd::Vector Eigen_vector;
typedef CGAL::Eigen_svd::Matrix Eigen_matrix;

namespace SMS = CGAL::Surface_mesh_simplification;

void save_mesh(const Surface_mesh &mesh, const Raster &raster, const char *filename);

std::pair <Eigen_matrix, Eigen_vector> volume_preservation (const SMS::Edge_profile<Surface_mesh>& profile) {

	K::Vector_3 n (CGAL::NULL_VECTOR);
	K::FT det = 0;

	SMS::Edge_profile<Surface_mesh>::Triangle_vector triangles = profile.triangles();
	for (auto triangle: triangles) {
		if (!((triangle.v0 == profile.v0() && (triangle.v1 == profile.v1() || triangle.v2 == profile.v1()))
			|| (triangle.v1 == profile.v0() && (triangle.v0 == profile.v1() || triangle.v2 == profile.v1()))
			|| (triangle.v2 == profile.v0() && (triangle.v0 == profile.v1() || triangle.v1 == profile.v1()))
		)) {
			Point_3 p0 = get(profile.vertex_point_map(),triangle.v0);
			Point_3 p1 = get(profile.vertex_point_map(),triangle.v1);
			Point_3 p2 = get(profile.vertex_point_map(),triangle.v2);

			if (!CGAL::collinear (p0, p1, p2)) {
				/*std::cout << "triangle: " << p0 << ", " << p1 << ", " << p2 << "\n";
				std::cout << "normal: " << CGAL::normal(p0, p1, p2) << "\n";
				std::cout << "determinant: " << CGAL::determinant(K::Vector_3(Point_3(CGAL::ORIGIN), p0), K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p2)) << "\n";*/
				
				n += CGAL::normal(p0, p1, p2);
				det += CGAL::determinant(K::Vector_3(Point_3(CGAL::ORIGIN), p0), K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p2));
			}
		} else {
			//std::cout << "no way" << "\n";
		}
	}

	Eigen_matrix A(1,3);
	A.set(0, 0, n.x());
	A.set(0, 1, n.y());
	A.set(0, 2, n.z());

	Eigen_vector B(1);
	B.set(0,det);

	A *= 1.732; // sqrt(3)
	B *= 1.732;

	std::cout << "A: " << A << "\n";
	std::cout << "B: " << B << "\n";

	return std::make_pair(A,B);
}

std::pair <Eigen_matrix, Eigen_vector> boundary_preservation (const SMS::Edge_profile<Surface_mesh>& profile, K::Vector_3 * ee1 = nullptr, K::Vector_3 * ee2 = nullptr) {
	K::Vector_3 e1 (CGAL::NULL_VECTOR);
	K::Vector_3 e2 (CGAL::NULL_VECTOR);

	for (auto edge: profile.border_edges()) {
		Point_3 p0 = get(profile.vertex_point_map(), profile.surface_mesh().source(edge));
		Point_3 p1 = get(profile.vertex_point_map(), profile.surface_mesh().target(edge));

		//std::cout << "edge: " << p0 << ", " << p1 << "\n";

		e1 += K::Vector_3(p0, p1);
		e2 += CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0));

		/*std::cout << "e1: " << K::Vector_3(p0, p1) << "\n";
		std::cout << "e2: " << CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0)) << "\n";*/
	}

	if (ee1 != nullptr) *ee1 = e1;
	if (ee1 != nullptr) *ee1 = e1;

	K::Vector_3 e3 = CGAL::cross_product(e1, e2);

	/*std::cout << "E1: " << e1 << "\n";
	std::cout << "E2: " << e2 << "\n";
	std::cout << "E3: " << e3 << "\n";*/

	K::Vector_3 r1 = CGAL::scalar_product (e1, e1) * e3;
	K::FT r2 = CGAL::scalar_product (e3, e3);

	K::Vector_3 r3 = CGAL::cross_product(e1, e3);

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();

	Eigen_matrix A(2,3);
	A.set(0, 0, r1.x());
	A.set(0, 1, r1.y());
	A.set(0, 2, r1.z());
	A.set(1, 0, r3.x());
	A.set(1, 1, r3.y());
	A.set(1, 2, r3.z());

	Eigen_vector B(2);
	B.set(0, r2);
	B.set(1, 0);

	A *= 1.225*squared_length; // sqrt(3/2)
	B *= 1.225*squared_length;

	std::cout << "A: " << A << "\n";
	std::cout << "B: " << B << "\n";

	return std::make_pair(A,B);
}

std::pair <Eigen_matrix, Eigen_vector> volume_optimisation (const SMS::Edge_profile<Surface_mesh>& profile, K::FT *constant = nullptr) {
	Eigen_matrix A(3,3);
	K::Vector_3 b (CGAL::NULL_VECTOR);

	if (constant != nullptr) *constant = 0;

	SMS::Edge_profile<Surface_mesh>::Triangle_vector triangles = profile.triangles();
	for (auto triangle: triangles) {
		if (!((triangle.v0 == profile.v0() && (triangle.v1 == profile.v1() || triangle.v2 == profile.v1()))
			|| (triangle.v1 == profile.v0() && (triangle.v0 == profile.v1() || triangle.v2 == profile.v1()))
			|| (triangle.v2 == profile.v0() && (triangle.v0 == profile.v1() || triangle.v1 == profile.v1()))
		)) {
			Point_3 p0 = get(profile.vertex_point_map(),triangle.v0);
			Point_3 p1 = get(profile.vertex_point_map(),triangle.v1);
			Point_3 p2 = get(profile.vertex_point_map(),triangle.v2);
			
			if (!CGAL::collinear (p0, p1, p2)) {
				K::Vector_3 n = CGAL::normal(p0, p1, p2);
				Eigen_vector nd (3);
				nd.set(0, n.x());
				nd.set(1, n.y());
				nd.set(2, n.z());

				A += nd * nd.transpose();
				K::FT det = CGAL::determinant(K::Vector_3(Point_3(CGAL::ORIGIN), p0), K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p2));
				b += det * n;
				if (constant != nullptr) *constant += det*det;
			}
		}
	}

	Eigen_vector B(3);
	B.set(0, b.x());
	B.set(1, b.y());
	B.set(2, b.z());

	A /= 18;
	B /= 18;
	if (constant != nullptr) *constant /= 18*2;

	std::cout << "A: " << A << "\n";
	std::cout << "B: " << B << "\n";

	return std::make_pair(A,B);
}

std::pair <Eigen_matrix, Eigen_vector> boundary_optimization (const SMS::Edge_profile<Surface_mesh>& profile, K::FT *constant = nullptr) {
	Eigen_matrix A(3,3);
	K::Vector_3 b (CGAL::NULL_VECTOR);

	if (constant != nullptr) *constant = 0;

	for (auto edge: profile.border_edges()) {
		Point_3 p0 = get(profile.vertex_point_map(), profile.surface_mesh().source(edge));
		Point_3 p1 = get(profile.vertex_point_map(), profile.surface_mesh().target(edge));

		K::Vector_3 e1 (p0, p1);
		K::Vector_3 e2 = CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), p1), K::Vector_3(Point_3(CGAL::ORIGIN), p0));

		Eigen_matrix e1_cross (3,3);
		e1_cross.set(1,0, e1.z());
		e1_cross.set(2,0, -e1.y());
		e1_cross.set(0,1, -e1.z());
		e1_cross.set(2,1, e1.x());
		e1_cross.set(0,2, e1.y());
		e1_cross.set(1,2, -e1.x());

		A += e1_cross.transpose() * e1_cross;
		b += - CGAL::cross_product(e1, e2);
		if (constant != nullptr) *constant += CGAL::scalar_product(e2,e2);
	}

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();

	Eigen_vector B(3);

	B.set(0, b.x());
	B.set(1, b.y());
	B.set(2, b.z());

	A *= squared_length/2;
	B *= squared_length/2;
	if (constant != nullptr) *constant *= squared_length/(2*2);

	std::cout << "A: " << A << "\n";
	std::cout << "B: " << B << "\n";

	return std::make_pair(A,B);
}

std::pair <Eigen_matrix, Eigen_vector> triangle_shape_optimization (const SMS::Edge_profile<Surface_mesh>& profile, K::FT *constant = nullptr) {
	int a;
	K::Vector_3 b (CGAL::NULL_VECTOR);

	if (constant != nullptr) *constant = 0;

	for (auto v: profile.link()) {
		K::Vector_3 vi (K::Point_3(CGAL::ORIGIN), get(profile.vertex_point_map(), v));

		a += 1;
		b += vi;

		if (constant != nullptr) *constant += CGAL::scalar_product(vi,vi);
	}

	auto squared_length = K::Vector_3(profile.p0(), profile.p1()).squared_length();

	Eigen_matrix A(3,3);
	Eigen_vector B(3);

	A.set(0,0,a);
	A.set(1,1,a);
	A.set(2,2,a);

	B.set(0, b.x());
	B.set(1, b.y());
	B.set(2, b.z());

	A *= 2*squared_length*squared_length;
	B *= 2*squared_length*squared_length;
	if (constant != nullptr) *constant *= squared_length*squared_length;

	std::cout << "A: " << A << "\n";
	std::cout << "B: " << B << "\n";

	return std::make_pair(A,B);
}

struct LindstromTurk_param {
	float volume_preservation;
	float boundary_preservation;
	float volume_optimisation;
	float boundary_optimization;
	float triangle_shape_optimization;

	LindstromTurk_param(float volume_preservation,
						float boundary_preservation,
						float volume_optimisation,
						float boundary_optimization,
						float triangle_shape_optimization) :
						volume_preservation(volume_preservation),
						boundary_preservation(boundary_preservation),
						volume_optimisation(volume_optimisation),
						boundary_optimization(boundary_optimization),
						triangle_shape_optimization(triangle_shape_optimization) {}
};

class Custom_placement {
	const Raster &raster;
	LindstromTurk_param params;

	public:
		Custom_placement (const Raster &raster) : raster(raster), params(LindstromTurk_param(0,0,1,0,0.1)) {}

		boost::optional<SMS::Edge_profile<Surface_mesh>::Point> operator()(const SMS::Edge_profile<Surface_mesh>& profile) const {
			typedef boost::optional<SMS::Edge_profile<Surface_mesh>::Point> result_type;

			Eigen_vector B(12);
			Eigen_matrix A(12, 3);

			auto r1 = volume_preservation(profile);
			for(std::size_t i=0; i<1; ++i) {
				for(std::size_t j = 0; j < 3; ++j) {
					A.set(i, j, r1.first(i,j)*params.volume_preservation);
				}
				B.set(i, r1.second(i)*params.volume_preservation);
			}

			auto r2 = boundary_preservation(profile);
			for(std::size_t i=0; i<2; ++i) {
				for(std::size_t j = 0; j < 3; ++j) {
					A.set(i+1, j, r2.first(i,j)*params.boundary_preservation);
				}
				B.set(i+1, r2.second(i)*params.boundary_preservation);
			}

			auto r3 = volume_optimisation(profile);
			for(std::size_t i=0; i<3; ++i) {
				for(std::size_t j = 0; j < 3; ++j) {
					A.set(i+3, j, r3.first(i,j)*params.volume_optimisation);
				}
				B.set(i+3, r3.second(i)*params.volume_optimisation);
			}

			auto r4 = boundary_optimization(profile);
			for(std::size_t i=0; i<3; ++i) {
				for(std::size_t j = 0; j < 3; ++j) {
					A.set(i+6, j, r4.first(i,j)*params.boundary_optimization);
				}
				B.set(i+6, r4.second(i)*params.boundary_optimization);
			}

			auto r5 = triangle_shape_optimization(profile);
			for(std::size_t i=0; i<3; ++i) {
				for(std::size_t j = 0; j < 3; ++j) {
					A.set(i+9, j, r5.first(i,j)*params.triangle_shape_optimization);
				}
				B.set(i+9, r5.second(i)*params.triangle_shape_optimization);
			}

			/*std::cout << "------------\n";

			std::cout << r3.first << "\n";
			std::cout << r3.second << "\n";*/

			std::cout << A << "\n";
			std::cout << B << "\n";

			std::cout << "---------------------\n";

			// Solve AX=B
			auto C = B;
			CGAL::Eigen_svd::solve(A, B);
			
			// Print result
			Point_3 placement(B.vector()[0], B.vector()[1], B.vector()[2]);

			std::cout << "p0: " << profile.p0() << ", p1: " << profile.p1() << ", p: " << placement << "\n";

			std::cout << A*B - C << "\n";

			std::cout << "---------------------\n";

			std::cout << "test---------------------\n";

			// Solve AX=B
			auto r = CGAL::midpoint(profile.p0(), profile.p1());

			if (CGAL::squared_distance(placement, r) > 1) std::cout << "ALERT : " << CGAL::squared_distance(placement, r) << "\n";

			/*B.set(0, r.x());
			B.set(1, r.y());
			B.set(2, r.z());
			
			// Print result
			placement = Point_3(B.vector()[0], B.vector()[1], B.vector()[2]);

			std::cout << "p0: " << profile.p0() << ", p1: " << profile.p1() << ", p: " << placement << "\n";


			std::cout << A*B - C << "\n";*/

			std::cout << "---------------------\n";

			return result_type(placement);
		}
};

class Custom_cost {
	const Raster &raster;
	LindstromTurk_param params;

	public:
		Custom_cost (const Raster &raster) : raster(raster), params(LindstromTurk_param(0,0,1,0,0.1)) {}

		boost::optional<SMS::Edge_profile<Surface_mesh>::FT> operator()(const SMS::Edge_profile<Surface_mesh>& profile, const boost::optional<SMS::Edge_profile<Surface_mesh>::Point>& placement) const {
			typedef boost::optional<SMS::Edge_profile<Surface_mesh>::FT> result_type;

			if (placement) {
				K::FT cost = 0;

				Eigen_matrix v(3,1);
				v.set(0,0,placement->x());
				v.set(1,0,placement->y());
				v.set(2,0,placement->z());

				auto r1 = volume_preservation(profile);
				cost += pow((r1.first * v)(0,0) / r1.second(0,0), 2)*params.volume_preservation;

				K::Vector_3 e1, e2;
				auto r2 = boundary_preservation(profile, &e1, &e2);
				cost += (CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), *placement), e1) + e2).squared_length()/4*params.boundary_preservation;

				K::FT constant3;
				auto r3 = volume_optimisation(profile, &constant3);
				cost += ((v.transpose()*r3.first*v/2-r3.second*v)(0,0)+constant3)*params.volume_optimisation;

				K::FT constant4;
				auto r4 = boundary_optimization(profile, &constant4);
				cost += ((v.transpose()*r4.first*v/2-r4.second*v)(0,0)+constant4)*params.boundary_optimization;

				K::FT constant5;
				auto r5 = triangle_shape_optimization(profile, &constant5);
				cost += ((v.transpose()*r5.first*v/2-r5.second*v)(0,0)+constant5)*params.triangle_shape_optimization;

				std::cout << "---------------------\n";

				std::cout << "volume_preservation: " << pow((r1.first * v)(0,0) / r1.second(0,0), 2)*params.volume_preservation << "\n";
				std::cout << "boundary_preservation: " << (CGAL::cross_product(K::Vector_3(Point_3(CGAL::ORIGIN), *placement), e1) + e2).squared_length()/4*params.boundary_preservation << "\n";
				std::cout << "volume_optimisation: " << ((v.transpose()*r3.first*v/2-r3.second*v)(0,0)+constant3)*params.volume_optimisation << "\n";
				std::cout << "boundary_optimization: " << ((v.transpose()*r4.first*v/2-r4.second*v)(0,0)+constant4)*params.boundary_optimization << "\n";
				std::cout << "triangle_shape_optimization: " << ((v.transpose()*r5.first*v/2-r5.second*v)(0,0)+constant5)*params.triangle_shape_optimization << "\n";

				std::cout << "---------------------\n";

				return result_type(cost);
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

	//Cost_stop_predicate stop(10);
	SMS::Count_stop_predicate<Surface_mesh> stop(1000);
	Custom_cost cf(raster);
	Custom_placement pf(raster);
	int r = SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).get_placement(pf).visitor(My_visitor(&mesh, &raster)));
	std::cout << "\rMesh simplified                                               " << std::endl;

	save_mesh(mesh, raster, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, mesh);
}
