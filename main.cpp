#include "header.hpp"
#include "raster.hpp"
#include "edge_collapse.hpp"

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <getopt.h>
#include <cstdlib>

void add_label(const Raster &raster, Surface_mesh &mesh);

void change_vertical_faces(Surface_mesh &mesh);

void compute_normal_angle_coef(Surface_mesh &mesh);

std::vector<std::list<Surface_mesh::Face_index>> compute_path(Surface_mesh &mesh);

std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> compute_path_polygon(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const Surface_mesh_info &mesh_info);

std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> compute_medial_axes(const Surface_mesh &mesh, const std::vector<std::list<Surface_mesh::Face_index>> &paths, const std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> &path_polygon, const Surface_mesh_info &mesh_info);

#include "bridge.hpp"

Surface_mesh_info::Surface_mesh_info() : x_0(0), y_0(0) {}
Surface_mesh_info::Surface_mesh_info(OGRSpatialReference crs, double x_0, double y_0) : crs(crs), x_0(x_0), y_0(y_0) {}

void Surface_mesh_info::save_mesh(const Surface_mesh &mesh, const char *filename) const {
	Surface_mesh output_mesh (mesh);

	Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> label;
	bool has_label;
	boost::tie(label, has_label) = output_mesh.property_map<Surface_mesh::Face_index, unsigned char>("f:label");
	if (has_label) {
		// Color
		bool created;
		Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> red;
		Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> green;
		Surface_mesh::Property_map<Surface_mesh::Face_index, unsigned char> blue;
		boost::tie(red, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("red",0);
		assert(created);
		boost::tie(green, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("green",0);
		assert(created);
		boost::tie(blue, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("blue",0);
		assert(created);

		for (auto face : output_mesh.faces()) {
			red[face] = LABELS.at(label[face]).red;
			green[face] = LABELS.at(label[face]).green;
			blue[face] = LABELS.at(label[face]).blue;
		}
	}

	Surface_mesh::Property_map<Surface_mesh::Face_index, int> path;
	bool has_path;
	boost::tie(path, has_path) = output_mesh.property_map<Surface_mesh::Face_index, int>("path");
	if (has_path) {
		// Entropy
		bool created;
		Surface_mesh::Property_map<Surface_mesh::Face_index, float> quality;
		boost::tie(quality, created) = output_mesh.add_property_map<Surface_mesh::Face_index, float>("quality",0);
		assert(created);
		
		for (auto face : output_mesh.faces()) {
			quality[face] = path[face];
		}
	}

	char *temp;
	/*const char *options_wkt[] = { "MULTILINE=NO", "FORMAT=WKT2", NULL };
	crs.exportToWkt(&temp, options_wkt);*/
	crs.exportToProj4(&temp); // WKT format is too long for MeshLab
	std::string crs_as_string(temp);
	CPLFree(temp);

	std::ofstream mesh_ofile (filename, std::ios_base::binary);
	CGAL::IO::set_binary_mode (mesh_ofile);
	CGAL::IO::write_PLY (mesh_ofile, output_mesh, "crs " + crs_as_string + "\nx0 " + std::to_string(x_0) + "\ny0 " + std::to_string(y_0) );
	mesh_ofile.close();
}

std::tuple<Surface_mesh, std::tuple<Surface_mesh, Point_set>> compute_meshes(const Raster &raster, const Surface_mesh_info &mesh_info) {

	std::cout << "Terrain mesh" << std::endl;
	Surface_mesh terrain_mesh;
	
	double x, y;
	
	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> terrain_vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			raster.grid_to_coord(P, L, x, y);
			terrain_vertex_index[L][P] = terrain_mesh.add_vertex(Point_3(x - mesh_info.x_0, y - mesh_info.y_0, raster.dtm[L][P]));
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

	// Return mesh if coords are in reverse order
	double x_0, y_0, x_1, y_1;	
	raster.grid_to_coord(0, 0, x_0, y_0);
	raster.grid_to_coord(1, 1, x_1, y_1);
	if ((x_1-x_0)*(y_1-y_0) < 0) {
		CGAL::Polygon_mesh_processing::reverse_face_orientations(terrain_mesh); 	
	}

	mesh_info.save_mesh(terrain_mesh, "initial-terrain-mesh.ply");

	SMS::edge_collapse(terrain_mesh, Cost_stop_predicate(10));
	std::cout << "Terrain mesh simplified" << std::endl;

	mesh_info.save_mesh(terrain_mesh, "terrain-mesh.ply");

	std::cout << "Surface mesh" << std::endl;
	Surface_mesh mesh;

	Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
	bool created_point_in_face;
	boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points", std::list<Point_set::Index>());
	assert(created_point_in_face);

	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
	bool created_face_costs;
	boost::tie(face_costs, created_face_costs) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:cost", 0);
	assert(created_face_costs);

	Point_set point_cloud;
	bool created_point_label;
	Point_set::Property_map<unsigned char> point_cloud_label;
	boost::tie (point_cloud_label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", LABEL_OTHER);
	assert(created_point_label);

	// Add points
	std::vector<std::vector<Surface_mesh::Vertex_index>> vertex_index(raster.ySize, std::vector<Surface_mesh::Vertex_index>(raster.xSize, Surface_mesh::Vertex_index()));
	std::vector<std::vector<Point_set::Index>> point_index(raster.ySize, std::vector<Point_set::Index>(raster.xSize, Point_set::Index()));
	for (int L = 0; L < raster.ySize; L++) {
		for (int P = 0; P < raster.xSize; P++) {
			raster.grid_to_coord(P, L, x, y);
			vertex_index[L][P] = mesh.add_vertex(Point_3(x - mesh_info.x_0, y - mesh_info.y_0, raster.dsm[L][P]));
			point_index[L][P] = *(point_cloud.insert(Point_set::Point_3(x - mesh_info.x_0, y - mesh_info.y_0, raster.dsm[L][P])));
			point_cloud_label[point_index[L][P]] = raster.land_cover[L][P];
		}
	}
	std::cout << "Point added" << std::endl;

	// Add faces
	for (int L = 0; L < raster.ySize-1; L++) {
		for (int P = 0; P < raster.xSize-1; P++) {
			if (pow(raster.dsm[L][P]-raster.dsm[L+1][P+1], 2) < pow(raster.dsm[L+1][P]-raster.dsm[L][P+1], 2)) {
				auto f1 = mesh.add_face(vertex_index[L][P], vertex_index[L+1][P+1], vertex_index[L+1][P]);
				face_costs[f1] = 0;
				auto f2 = mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P+1]);
				face_costs[f2] = 0;
				point_in_face[f1].push_back(point_index[L][P]);
				if (L == raster.ySize - 2) point_in_face[f1].push_back(point_index[L+1][P]);
				if (P == raster.xSize - 2) point_in_face[f2].push_back(point_index[L][P+1]);
				if (L == raster.ySize - 2 && P == raster.xSize - 2) point_in_face[f2].push_back(point_index[L+1][P+1]);
			} else {
				auto f1 = mesh.add_face(vertex_index[L][P], vertex_index[L][P+1], vertex_index[L+1][P]);
				face_costs[f1] = 0;
				auto f2 = mesh.add_face(vertex_index[L][P+1], vertex_index[L+1][P+1], vertex_index[L+1][P]);
				face_costs[f2] = 0;
				point_in_face[f1].push_back(point_index[L][P]);
				if (L == raster.ySize - 2) point_in_face[f1].push_back(point_index[L+1][P]);
				if (P == raster.xSize - 2) point_in_face[f2].push_back(point_index[L][P+1]);
				if (L == raster.ySize - 2 && P == raster.xSize - 2) point_in_face[f2].push_back(point_index[L+1][P+1]);
			}
		}
	}
	std::cout << "Faces added" << std::endl;

	float alpha = 2, beta = 1, gamma = 0.01;

	K::FT mean_point_per_area = get_mean_point_per_area(mesh, point_cloud);
	K::FT min_point_per_area = mean_point_per_area / 2;

	for(auto face: mesh.faces()) {
		K::FT min_point = 0;
		if (min_point_per_area > 0) {
			auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
			K::FT face_area = CGAL::sqrt(K::Triangle_3(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++)).squared_area());
			min_point = min_point_per_area * face_area;
		}	
		if (point_in_face[face].size() > 0) {
			if (point_in_face[face].size() > min_point) {
				int face_label[LABELS.size()] = {0};
				for (auto point: point_in_face[face]) {
					face_label[point_cloud_label[point]]++;
				}
				auto argmax = std::max_element(face_label, face_label+LABELS.size());
				face_costs[face] += beta * (point_in_face[face].size() - *argmax);
			} else { // We don't have information about this face.
				face_costs[face] += beta * point_in_face[face].size();
			}
		}
	}
	std::cout << "Faces cost" << std::endl;

	// Return mesh if coords are in reverse order
	if ((x_1-x_0)*(y_1-y_0) < 0) {
		CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh); 	
	}

	Cost_stop_predicate stop(5);
	//SMS::Count_stop_predicate<Surface_mesh> stop(50);
	const LindstromTurk_param params (10,1,10,1,0.000001,1,0.01);
	Custom_placement pf(params, mesh, point_cloud);
	Custom_cost cf(params, alpha, beta, gamma, 0.01, min_point_per_area, mesh, point_cloud);
	My_visitor mv (params, alpha, beta, gamma, min_point_per_area, mesh, mesh_info, point_cloud);
	SMS::Bounded_normal_change_filter<> filter;
	SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(mv));

	mesh_info.save_mesh(mesh, "final-mesh.ply");

	return std::make_tuple(terrain_mesh, std::make_tuple(mesh, point_cloud));
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
		{"point_cloud", required_argument, NULL, 'P'},
		{NULL, 0, 0, '\0'}
	};

	char *DSM = NULL;
	char *DTM = NULL;
	char *land_use_map = NULL;
	char *LOD0 = NULL;
	char *orthophoto = NULL;
	char *MESH = NULL;
	char *POINT_CLOUD = NULL;

	while ((opt = getopt_long(argc, argv, "hs:t:l:0:i:M:P:", options, NULL)) != -1) {
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
				std::cout << " -P, --point_cloud=/file/path.ply   point cloud as PLY file." << std::endl;
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
			case 'P':
				POINT_CLOUD = optarg;
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

	std::srand(38401);

	const Raster raster(DSM, DTM, land_use_map);

	double min_x, min_y;
	raster.grid_to_coord(0, 0, min_x, min_y);
	Surface_mesh_info mesh_info(raster.get_crs(), min_x, min_y);

	Surface_mesh terrain_mesh, mesh;
	Point_set point_cloud;
	if (MESH == NULL) {
		std::tuple<Surface_mesh&, Point_set&> nested = std::tie(mesh, point_cloud);
		std::tie(terrain_mesh, nested) = compute_meshes(raster, mesh_info);

		std::ofstream mesh_ofile ("save_mesh.ply", std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ofile);
		CGAL::IO::write_PLY (mesh_ofile, mesh);
		mesh_ofile.close();

		mesh_ofile = std::ofstream("save_pointcloud.ply", std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ofile);
		CGAL::IO::write_PLY (mesh_ofile, point_cloud);
		mesh_ofile.close();

		mesh_ofile = std::ofstream("save_terrain_mesh.ply", std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ofile);
		CGAL::IO::write_PLY (mesh_ofile, terrain_mesh);
		mesh_ofile.close();

		return EXIT_SUCCESS;

	} else {
		std::ifstream mesh_ifile (MESH, std::ios_base::binary);
		CGAL::IO::set_binary_mode (mesh_ifile);
		CGAL::IO::read_PLY (mesh_ifile, mesh);
		mesh_ifile.close();
		std::cout << "Mesh load" << std::endl;

		if (POINT_CLOUD != NULL) {
			std::ifstream mesh_ifile (POINT_CLOUD, std::ios_base::binary);
			CGAL::IO::set_binary_mode (mesh_ifile);
			CGAL::IO::read_PLY (mesh_ifile, point_cloud);
			mesh_ifile.close();
			std::cout << "Point_cloud load" << std::endl;
		} else {
			point_cloud = compute_point_cloud(mesh);
			std::cout << "Point_cloud compute" << std::endl;
		}

		associate_mesh_point_cloud(mesh, point_cloud);
		std::cout << "Mesh and point cloud associate" << std::endl;
	}

	change_vertical_faces(mesh);
	mesh_info.save_mesh(mesh, "final-mesh-without-facade.ply");
	std::cout << "Label set for vertical face" << std::endl;

	compute_normal_angle_coef(mesh);
	mesh_info.save_mesh(mesh, "final-mesh-with-normal-angle-ceof.ply");
	std::cout << "Normal angle coef computed" << std::endl;

	std::vector<std::list<Surface_mesh::Face_index>> paths = compute_path(mesh);
	mesh_info.save_mesh(mesh, "final-mesh-with-path.ply");
	std::cout << "Path computed" << std::endl;

	std::map<int, CGAL::Polygon_with_holes_2<Exact_predicates_kernel>> path_polygon = compute_path_polygon(mesh, paths, mesh_info);
	std::map<int, boost::shared_ptr<CGAL::Straight_skeleton_2<K>>> medial_axes = compute_medial_axes(mesh, paths, path_polygon, mesh_info);
	std::cout << "Medial axes computed" << std::endl;

	std::set<pathLink> links = link_paths(mesh, paths, path_polygon, medial_axes, mesh_info);
	std::cout << "Links computed" << std::endl;

	close_surface_mesh(mesh);

	mesh_info.save_mesh(mesh, "final-closed-mesh-with-path.ply");
	std::cout << "Mesh waterthighted" << std::endl;

	AABB_tree tree = index_surface_mesh(mesh);

	std::vector<pathBridge> bridges_to_add;
	int i = 0;
	std::cout << "Computing " << links.size() << " bridges" << std::endl;
	for (auto link: links) {
		std::cout << "\rBridge " << i++ << "/" << links.size() << "               ";
		std::cout.flush();
		pathBridge bridge_result = bridge(link, mesh, tree, mesh_info);
		if (bridge_result.cost < 10) {
			bridges_to_add.push_back(bridge_result);
		}
	}
	std::cout << "\rBridges computed               " << std::endl;

	add_bridge_to_mesh(mesh, bridges_to_add, path_polygon, mesh_info);

	mesh_info.save_mesh(mesh, "final-closed-mesh-with-path-and-bridges.ply");
	std::cout << "Bridges added to mesh" << std::endl;

	return EXIT_SUCCESS;
}
