#include "header.hpp"
#include "edge_collapse.hpp"

#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

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
	std::ofstream myfile;
	myfile.open ("commande_line.txt");
	std::copy( argv+1, argv+argc, std::ostream_iterator<const char*>( myfile, " " ) );
	myfile << "\n";

	int opt;
	int option_index = 0;
	const struct option options[] = {
		{"help", no_argument, NULL, 'h'},
		{"mesh", required_argument, NULL, 'm'},
		{"point_cloud", required_argument, NULL, 'p'},
		{"l1", required_argument, NULL, 0},
		{"l2", required_argument, NULL, 0},
		{"l3", required_argument, NULL, 0},
		{"l4", required_argument, NULL, 0},
		{"l5", required_argument, NULL, 0},
		{"l6", required_argument, NULL, 0},
		{"l7", required_argument, NULL, 0},
		{"c1", required_argument, NULL, 0},
		{"c2", required_argument, NULL, 0},
		{"c3", required_argument, NULL, 0},
		{"c4", required_argument, NULL, 0},
		{"cs", required_argument, NULL, 0},
		{"ns", required_argument, NULL, 0},
		{"baseline", required_argument, NULL, 0},
		{NULL, 0, 0, '\0'}
	};

	char *mesh_file = NULL;
	char *point_cloud_file = NULL;
	float l1=10, l2=1, l3=10, l4=1, l5=0.001, l6=1, l7=1, c1=1, c2=1, c3=1, c4=1, cs=1;
	int ns = 0;
	int baseline = -1;

	while ((opt = getopt_long(argc, argv, "hm:p:", options, &option_index)) != -1) {
		switch(opt) {
			case 0:
				switch(option_index) {
					case 3:
						l1 = atof(optarg);
						break;
					case 4:
						l2 = atof(optarg);
						break;
					case 5:
						l3 = atof(optarg);
						break;
					case 6:
						l4 = atof(optarg);
						break;
					case 7:
						l5 = atof(optarg);
						break;
					case 8:
						l6 = atof(optarg);
						break;
					case 9:
						l7 = atof(optarg);
						break;
					case 10:
						c1 = atof(optarg);
						break;
					case 11:
						c2 = atof(optarg);
						break;
					case 12:
						c3 = atof(optarg);
						break;
					case 13:
						c4 = atof(optarg);
						break;
					case 14:
						cs = atof(optarg);
						break;
					case 15:
						ns = atoi(optarg);
						break;
					case 16:
						baseline = atoi(optarg);
						break;
				}
				break;
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
	std::cout << "l1=" << l1 << ", l2=" << l2 << ", l3=" << l3 << ", l4=" << l4 << ", l5=" << l5 << ", l6=" << l6 << ", l7=" << l7 << ", c1=" << c1 << ", c2=" << c2 << ", c3=" << c3 << ", c4=" << c4 << "\n";
	if (ns > 0) {
		std::cout << "ns=" << ns << "\n";
	} else {
		std::cout << "cs=" << cs << "\n";
	}
	myfile << "Mesh: " << mesh_file << std::endl;
	if (point_cloud_file != NULL) myfile << "Point cloud: " << point_cloud_file << std::endl;
	myfile << "l1=" << l1 << ", l2=" << l2 << ", l3=" << l3 << ", l4=" << l4 << ", l5=" << l5 << ", l6=" << l6 << ", l7=" << l7 << ", c1=" << c1 << ", c2=" << c2 << ", c3=" << c3 << ", c4=" << c4 << "\n";
	if (ns > 0) {
		myfile << "ns=" << ns << "\n";
	} else {
		myfile << "cs=" << cs << "\n";
	}
	myfile.close();

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

	// Create point_in_face
	Surface_mesh::Property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>> point_in_face;
	bool created_point_in_face;
	boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::vector<Point_set::Index>>("f:points", std::vector<Point_set::Index>());
	assert(created_point_in_face);

	Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
	bool created_face_costs;
	boost::tie(face_costs, created_face_costs) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:points", 0);
	assert(created_face_costs);

	AABB_tree mesh_tree;
	PMP::build_AABB_tree(mesh, mesh_tree);
	CGAL::Cartesian_converter<Exact_predicates_kernel,K> type_converter;
	for(auto ph: point_cloud) {
		auto p = type_converter(point_cloud.point(ph));
		auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
		point_in_face[location.first].push_back(ph);
		face_costs[location.first] += c1 * CGAL::squared_distance(p, PMP::construct_point(location, mesh));
	}
	for(auto face: mesh.faces()) {
		if (point_in_face[face].size() > 0) {
			int face_label[LABELS.size()] = {0};
			int sum_face_label = 0;
			for (auto ph: point_in_face[face]) {
				sum_face_label++;
				face_label[label[ph]]++;
			}
			auto argmax = std::max_element(face_label, face_label+LABELS.size());
			face_costs[face] += c2 * (sum_face_label - *argmax);
		}
	}
	
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
	{
		bool created_point_isborder;
		Point_set::Property_map<int> isborder2;
		boost::tie (isborder2, created_point_isborder) = point_cloud.add_property_map<int>("p:isborder2", false);
		assert(created_point_isborder);
		for (auto point: point_cloud) {
			isborder2[point] = isborder[point];
		}
	}

	{
		bool created;
		Point_set::Property_map<unsigned char> red;
		Point_set::Property_map<unsigned char> green;
		Point_set::Property_map<unsigned char> blue;
		boost::tie(red, created) = point_cloud.add_property_map<unsigned char>("red",0);
		assert(created);
		boost::tie(green, created) = point_cloud.add_property_map<unsigned char>("green",0);
		assert(created);
		boost::tie(blue, created) = point_cloud.add_property_map<unsigned char>("blue",0);
		assert(created);

		for (auto point: point_cloud) {
			red[point] = LABELS[label[point]].red;
			green[point] = LABELS[label[point]].green;
			blue[point] = LABELS[label[point]].blue;
		}
	}

	CGAL::IO::write_point_set("pc_with_border.ply", point_cloud);
	std::cout << "Border points found" << std::endl;

	auto created_label = mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("label", LABEL_OTHER);
	assert(created_label.second);

	int min_surface = 10;
	K::FT mean_point_per_area = add_label(mesh, point_cloud, min_surface);
	mesh_info.save_mesh(mesh, "initial-mesh.ply");
	
	if(baseline >= 0) {
		SMS::Bounded_normal_change_filter<> filter;
		if (baseline == 0) {
			std::cout << "LindstromTurk" << std::endl;
			if (ns > 0) {
				SMS::Count_stop_predicate<Surface_mesh> stop(ns);
				SMS::edge_collapse(mesh, stop, CGAL::parameters::filter(filter));
			} else {
				Cost_stop_predicate stop(cs);
				SMS::edge_collapse(mesh, stop, CGAL::parameters::filter(filter));
			}
		} else {
			std::cout << "GarlandHeckbert_plane" << std::endl;
			SMS::GarlandHeckbert_policies<Surface_mesh, K> gh_policies(mesh);
			auto gh_cost = gh_policies.get_cost();
			auto gh_placement = gh_policies.get_placement();
			if (ns > 0) {
				SMS::Count_stop_predicate<Surface_mesh> stop(ns);
				SMS::edge_collapse(mesh, stop, CGAL::parameters::filter(filter)
																.get_cost(gh_cost)
																.get_placement(gh_placement));
			} else {
				Cost_stop_predicate stop(cs);
				SMS::edge_collapse(mesh, stop, CGAL::parameters::filter(filter)
																.get_cost(gh_cost)
																.get_placement(gh_placement));
			}
		}
		std::cout << "\rMesh simplified" << std::endl;
		mesh_info.save_mesh(mesh, "final-mesh.ply");
		return EXIT_SUCCESS;
	}

	Surface_mesh::Property_map<Surface_mesh::Halfedge_index, K::FT> placement_costs;
	bool created_placement_costs;
	boost::tie(placement_costs, created_placement_costs) = mesh.add_property_map<Surface_mesh::Halfedge_index, K::FT>("h:p_cost", 0);
	assert(created_placement_costs);
	
	const LindstromTurk_param params (l1,l2,l3,l4,l5,l6,l7);
	Custom_placement pf(params, mesh, point_cloud);
	Custom_cost cf(c1, c2, c3, c4, min_surface, mean_point_per_area, mesh, point_cloud);
	My_visitor mv(c1, c2, min_surface, mean_point_per_area, mesh, mesh_info, point_cloud);
	SMS::Bounded_normal_change_filter<> filter;
	if (ns > 0) {
		SMS::Count_stop_predicate<Surface_mesh> stop(ns);
		SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(mv));
	} else {
		Cost_stop_predicate stop(cs);
		SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(mv));
	}
	std::cout << "\rMesh simplified                                               " << std::endl;

	// add_label(mesh, point_cloud, min_surface, mean_point_per_area);
	mesh_info.save_mesh(mesh, "final-mesh.ply");

	return EXIT_SUCCESS;
}
