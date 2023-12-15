#include "header.hpp"
#include "edge_collapse.hpp"

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>

#include <getopt.h>
#include <cstdlib>

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
		boost::tie(green, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("green",0);
		boost::tie(blue, created) = output_mesh.add_property_map<Surface_mesh::Face_index, unsigned char>("blue",0);

		for (auto face : output_mesh.faces()) {
			red[face] = LABELS.at(label[face]).red;
			green[face] = LABELS.at(label[face]).green;
			blue[face] = LABELS.at(label[face]).blue;
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
	float l1=10, l2=1, l3=10, l4=1, l5=0.00001, l6=1, l7=0.01, c1=1, c2=1, c3=0.01, c4=1, cs=1;
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
	myfile << "Mesh=" << mesh_file << std::endl;
	if (point_cloud_file != NULL) {
		myfile << "Point cloud=" << point_cloud_file << std::endl;
	} else {
		myfile << "Point cloud=" << std::endl;
	}
	myfile << "l1=" << l1 << "\n";
	myfile << "l2=" << l2 << "\n";
	myfile << "l3=" << l3 << "\n";
	myfile << "l4=" << l4 << "\n";
	myfile << "l5=" << l5 << "\n";
	myfile << "l6=" << l6 << "\n";
	myfile << "l7=" << l7 << "\n";
	myfile << "c1=" << c1 << "\n";
	myfile << "c2=" << c2 << "\n";
	myfile << "c3=" << c3 << "\n";
	myfile << "c4=" << c4 << "\n";
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
	boost::tie (label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", LABEL_OTHER);
	if (created_point_label) {
		std::cout << "Point cloud has no label" << std::endl;

		if (point_cloud.has_property_map<int>("label")) {
			bool has_point_int_label;
			Point_set::Property_map<int> int_label;
			boost::tie (int_label, has_point_int_label) = point_cloud.property_map<int>("label");
			assert(has_point_int_label);
			for (auto ph: point_cloud) {
				int value = int_label[ph];
				if (value >= LABELS.size()) value = LABEL_OTHER;
				if (value < 0) value = LABEL_OTHER;
				label[ph] = static_cast<unsigned char>(value);
			}
			std::cout << "Point cloud label from int_label property" << std::endl;
		}

	} else {
		std::cout << "Point cloud has label" << std::endl;
	}

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

	K::FT mean_point_per_area = get_mean_point_per_area(mesh, point_cloud);
	K::FT min_point_per_area = mean_point_per_area / 10;
	
	const LindstromTurk_param params (l1,l2,l3,l4,l5,l6,l7);
	Custom_placement pf(params, mesh, point_cloud);
	Custom_cost cf(params, c1, c2, c3, c4, min_point_per_area, mesh, point_cloud);
	My_visitor mv(params, c1, c2, c3, min_point_per_area, mesh, mesh_info, point_cloud);
	SMS::Bounded_normal_change_filter<> filter;
	if (ns > 0) {
		SMS::Count_stop_predicate<Surface_mesh> stop(ns);
		SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(mv));
	} else {
		Cost_stop_predicate stop(cs);
		SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cf).filter(filter).get_placement(pf).visitor(mv));
	}

	mesh_info.save_mesh(mesh, "final-mesh.ply");

	return EXIT_SUCCESS;
}
