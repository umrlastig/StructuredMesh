#include "header.hpp"
#include "raster.hpp"

#include <getopt.h>
#include <cstdlib>

std::tuple<Surface_mesh, Surface_mesh> compute_meshes(const Raster &raster, const Surface_mesh_info &mesh_info);

void add_label(const Raster &raster, Surface_mesh &mesh);

void change_vertical_faces(Surface_mesh &mesh, const Raster &raster);

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
	boost::tie(label, has_label) = output_mesh.property_map<Surface_mesh::Face_index, unsigned char>("label");
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
			red[face] = LABELS[label[face]].red;
			green[face] = LABELS[label[face]].green;
			blue[face] = LABELS[label[face]].blue;
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

	std::srand(38401);

	const Raster raster(DSM, DTM, land_use_map);

	double min_x, min_y;
	raster.grid_to_coord(0, 0, min_x, min_y);
	Surface_mesh_info mesh_info(raster.get_crs(), min_x, min_y);

	Surface_mesh terrain_mesh, mesh;
	if (MESH == NULL || TERRAIN_MESH == NULL) {
		std::tie(terrain_mesh, mesh) = compute_meshes(raster, mesh_info);

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

	/*change_vertical_faces(mesh, raster);*/ // Need label information from point cloud
	mesh_info.save_mesh(mesh, "final-mesh-without-facade.ply");
	std::cout << "Label set for vertical face" << std::endl;

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
		if (bridge_result.cost < 50) {
			bridges_to_add.push_back(bridge_result);
		}
	}
	std::cout << "\rBridges computed               " << std::endl;

	add_bridge_to_mesh(mesh, bridges_to_add, path_polygon, mesh_info);

	mesh_info.save_mesh(mesh, "final-closed-mesh-with-path-and-bridges.ply");
	std::cout << "Bridges added to mesh" << std::endl;

	return EXIT_SUCCESS;
}
