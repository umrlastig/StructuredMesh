#include "header.hpp"
#include "edge_collapse.hpp"

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <getopt.h>
#include <cstdlib>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>        AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                AABB_tree;

Surface_mesh_info::Surface_mesh_info() : x_0(0), y_0(0) {}

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

	std::ofstream mesh_ofile (filename, std::ios_base::binary);
	CGAL::IO::set_binary_mode (mesh_ofile);
	CGAL::IO::write_PLY (mesh_ofile, output_mesh);
	mesh_ofile.close();
}

int main(int argc, char **argv) {
	int opt;
	const struct option options[] = {
		{"help", no_argument, NULL, 'h'},
		{"mesh", required_argument, NULL, 'm'},
		{"point_cloud", required_argument, NULL, 'p'},
		{NULL, 0, 0, '\0'}
	};

	char *mesh_file = NULL;
	char *point_cloud_file = NULL;

	while ((opt = getopt_long(argc, argv, "hm:p:", options, NULL)) != -1) {
		switch(opt) {
			case 'h':
				std::cout << "Usage: " << argv[0] << " [OPTIONS] -m mesh" << std::endl;
				std::cout << "Simplify mesh." << std::endl;
				std::cout << "OPTIONS:" << std::endl;
				std::cout << " -h, --help                         Print this help anq quit." << std::endl;
				std::cout << " -m, --mesh=/file/path.ply          mesh as PLY file." << std::endl;
				std::cout << " -p, --point_cloud=/file/path.ply   point cloud as PLY file." << std::endl;
				return EXIT_SUCCESS;
				break;
			case 'm':
				mesh_file = optarg;
				break;
			case 'p':
				point_cloud_file = optarg;
				break;
		}
	}

	if (mesh_file == NULL) {
		std::cerr << "mesh is mandatory" << std::endl;
		std::cerr << "Usage: " << argv[0] << " [OPTIONS] -m mesh" << std::endl;
		std::cerr << "Use -h/--help to obtain more informations" << std::endl;
		return EXIT_FAILURE;
	}

	std::cout << "Mesh: " << mesh_file << std::endl;
	if (point_cloud_file != NULL) std::cout << "Point cloud: " << point_cloud_file << std::endl;

	// Load mesh
	Surface_mesh mesh;
	Surface_mesh_info mesh_info;
	if (!CGAL::IO::read_polygon_mesh(mesh_file, mesh)) {
		throw std::invalid_argument(std::string("Unable to open ") + mesh_file + ".");
	}

	// Load point cloud
	Point_set point_cloud;
	if (point_cloud_file != NULL) {
		if (!CGAL::IO::read_point_set(point_cloud_file, point_cloud)) {
			throw std::invalid_argument(std::string("Unable to open ") + point_cloud_file + ".");
		}
	}

	// Create point_in_face
	std::map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	AABB_tree mesh_tree;
	PMP::build_AABB_tree(mesh, mesh_tree);
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;
	for(auto ph: point_cloud) {
		auto location = PMP::locate_with_AABB_tree(type_converter(point_cloud.point(ph)), mesh_tree, mesh);
		point_in_face[location.first].push_back(ph);
	}

	// Get or create created_point_label
	bool created_point_label;
	Point_set::Property_map<unsigned char> label;
	if (point_cloud.has_property_map<unsigned char>("p:label")) {
		boost::tie (label, created_point_label) = point_cloud.property_map<unsigned char>("p:label");
		std::cout << "Point cloud has label" << std::endl;
	} else {
		boost::tie (label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", LABEL_OTHER);
		std::cout << "Point cloud has no label" << std::endl;
	}
	assert(created_point_label);
	
	bool created_point_isborder;
	Point_set::Property_map<bool> isborder;
	boost::tie (isborder, created_point_isborder) = point_cloud.add_property_map<bool>("p:isborder", false);
	assert(created_point_isborder);

	// Set isBorder
	int N = 10;
	Point_tree point_tree(point_cloud.begin(), point_cloud.end(), Point_tree::Splitter(), TreeTraits(point_cloud.point_map()));
	Neighbor_search::Distance tr_dist(point_cloud.point_map());
	for (auto point: point_cloud) {
		Neighbor_search search(point_tree, point_cloud.point(point), N, 0, true, tr_dist, false);
		unsigned char l = label[search.begin()->first];
		for(auto it = search.begin() + 1; it != search.end() && !isborder[point]; ++it) {
			if (label[it->first] != l) isborder[point] = true;
		}
	}
	std::cout << "Border points found" << std::endl;

	auto created_label = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label", LABEL_OTHER);
	assert(created_label.second);

	add_label(mesh, point_cloud, point_in_face);

	mesh_info.save_mesh(mesh, "initial-mesh.ply");

	return EXIT_SUCCESS;
}
