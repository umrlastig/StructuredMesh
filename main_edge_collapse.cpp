#include "header.hpp"
#include "edge_collapse.hpp"

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_filter.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <getopt.h>
#include <cstdlib>
#include <random>
#include <filesystem>

typedef CGAL::Search_traits_3<Point_set_kernel>             Traits_base;
typedef CGAL::Search_traits_adapter<Point_set::Index, Point_set::Point_map, Traits_base> TreeTraits;
typedef CGAL::Orthogonal_k_neighbor_search<TreeTraits>      Neighbor_search;
typedef Neighbor_search::Tree                               Point_tree;

typedef CGAL::AABB_face_graph_triangle_primitive<Surface_mesh> AABB_face_graph_primitive;
typedef CGAL::AABB_traits<K, AABB_face_graph_primitive>        AABB_face_graph_traits;
typedef CGAL::AABB_tree<AABB_face_graph_traits>                AABB_tree;

namespace PMP = CGAL::Polygon_mesh_processing;

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
		{"subsample", required_argument, NULL, 0},
		{"baseline", required_argument, NULL, 0},
		{"min_point_factor", required_argument, NULL, 0},
		{"no_subdivide", no_argument, NULL, 0},
		{"no_direct_search", no_argument, NULL, 0},
		{"no_border_point", no_argument, NULL, 0},
		{"save_step_mesh", no_argument, NULL, 0},
		{"next_mesh", required_argument, NULL, 0},
		{NULL, 0, 0, '\0'}
	};

	char *mesh_file = NULL;
	char *point_cloud_file = NULL;
	float l1=10, l2=1, l3=10, l4=1, l5=0.00001, l6=1, l7=0.01, c1=2, c2=1, c3=0.01, c4=0.01, cs=1;
	int ns = 0;
	float subsample = -1;
	int baseline = -1;
	float min_point_factor = 10;
	bool subdivide = true, direct_search = true, border_point = true, step_mesh = false;
	char *next_mesh = NULL;

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
						subsample = atof(optarg);
						break;
					case 17:
						baseline = atoi(optarg);
						break;
					case 18:
						min_point_factor = atof(optarg);
						break;
					case 19:
						subdivide = false;
						break;
					case 20:
						direct_search = false;
						break;
					case 21:
						border_point = false;
						break;
					case 22:
						step_mesh = true;
						break;
					case 23:
						next_mesh = optarg;
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
	if (subsample != -1) {
		if (subsample <= 0) subsample = -1;
		if (subsample >= 1) {
			subsample = int(subsample);
			if (min_point_factor > 0) {
				std::cout << "subsample >= 1 and min_point_factor > 0 are incompatible, min_point_factor is set to 0.\n";
				min_point_factor = 0;
			}
		}
		std::cout << "subsample=" << subsample << "\n";
	}
	std::cout << "min_point_factor=" << min_point_factor << "\n";
	std::cout << "subdivide=" << int(subdivide) << "\n";
	std::cout << "direct_search=" << int(direct_search) << "\n";
	std::cout << "border_point=" << int(border_point) << "\n";
	
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
	myfile << "subsample=" << subsample << "\n";
	myfile << "min_point_factor=" << min_point_factor << "\n";
	myfile << "subdivide=" << int(subdivide) << "\n";
	myfile << "direct_search=" << int(direct_search) << "\n";
	myfile << "border_point=" << int(border_point) << "\n";
	myfile.close();

	if (next_mesh != nullptr) {
		if (!std::filesystem::exists(next_mesh)) {
			std::filesystem::create_directory(next_mesh);
		}
	}

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
	Point_set::Property_map<unsigned char> point_cloud_label;
	boost::tie (point_cloud_label, created_point_label) = point_cloud.add_property_map<unsigned char>("p:label", LABEL_OTHER);
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
				point_cloud_label[ph] = static_cast<unsigned char>(value);
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

	Ablation_study ablation (subdivide, direct_search, border_point,step_mesh);
	ablation.ground_truth_point_cloud = point_cloud;
	ablation.ground_truth_surface_mesh = mesh;

	K::FT min_point_per_area;
	if (min_point_factor > 0) {
		K::FT mean_point_per_area = get_mean_point_per_area(mesh, point_cloud);
		min_point_per_area = mean_point_per_area / min_point_factor;
	} else {
		min_point_per_area = 0;
	}

	if (subsample >= 0) {
		//subsample point cloud

		if (subsample >= 1) {

			// Keep only {subsample:int} points per face

			// Create point_in _face
			bool created_point_in_face;
			Surface_mesh::Property_map<Surface_mesh::Face_index, std::list<Point_set::Index>> point_in_face;
			boost::tie(point_in_face, created_point_in_face) = mesh.add_property_map<Surface_mesh::Face_index, std::list<Point_set::Index>>("f:points", std::list<Point_set::Index>());

			// Create face_costs
			bool created_face_costs;
			Surface_mesh::Property_map<Surface_mesh::Face_index, K::FT> face_costs;
			boost::tie(face_costs, created_face_costs) = mesh.add_property_map<Surface_mesh::Face_index, K::FT>("f:cost", 0);

			CGAL::Cartesian_converter<Point_set_kernel,K> type_converter;
			float alpha = c1;
			float beta = c2;

			//get border point
			bool created_point_isborder;
			Point_set::Property_map<bool> isborder;
			boost::tie (isborder, created_point_isborder) = point_cloud.add_property_map<bool>("p:isborder", false);
			if (border_point && created_point_isborder) {
				int N = 10;
				Point_tree point_tree(point_cloud.begin(), point_cloud.end(), Point_tree::Splitter(), TreeTraits(point_cloud.point_map()));
				Neighbor_search::Distance tr_dist(point_cloud.point_map());
				for (auto point: point_cloud) {
					Neighbor_search search(point_tree, point_cloud.point(point), N, 0, true, tr_dist, false);
					unsigned char l = point_cloud_label[search.begin()->first];
					for(auto it = search.begin() + 1; it != search.end() && !isborder[point]; ++it) {
						if (point_cloud_label[it->first] != l) isborder[point] = true;
					}
				}
				std::cout << "Border points computed" << std::endl;
			}

			// drop points in face to keep only subsample points per face
			std::random_device rd;
			std::mt19937 g(rd());
			for (auto face: mesh.faces()) {
				std::vector<Point_set::Index> temp(point_in_face[face].begin(), point_in_face[face].end());
				std::shuffle(temp.begin(), temp.end(), g);

				point_in_face[face].clear();
				for(std::size_t i = 0; i < temp.size(); i++) {
					if (i < subsample) {
						point_in_face[face].push_back(temp[i]);
					} else if (isborder[temp[i]]) {
						point_in_face[face].push_back(temp[i]);
					}
				}
			}

			if(created_point_in_face) {
				AABB_tree mesh_tree;
				PMP::build_AABB_tree(mesh, mesh_tree);
				for(auto ph: point_cloud) {
					auto p = type_converter(point_cloud.point(ph));
					auto location = PMP::locate_with_AABB_tree(p, mesh_tree, mesh);
					point_in_face[location.first].push_back(ph);
					if (created_face_costs && alpha > 0) face_costs[location.first] += alpha * CGAL::squared_distance(p, PMP::construct_point(location, mesh));
				}
			} else if (created_face_costs && alpha > 0) {
				for (auto face: mesh.faces()) {
					auto r = mesh.vertices_around_face(mesh.halfedge(face)).begin();
					K::Triangle_3 triangle_face(mesh.point(*r++), mesh.point(*r++), mesh.point(*r++));
					for (auto ph: point_in_face[face]) {
						auto p = type_converter(point_cloud.point(ph));
						face_costs[face] += alpha * CGAL::squared_distance(p, triangle_face);
					}
				}
			}
			if (created_face_costs && beta > 0) {
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
//if (face.idx() == 58591) std::cerr << "face_costs[58521] = " << face_costs[face] << "\n";
				}
			}

			std::cout << "Points subsampled with " << subsample << " points per face." << std::endl;

		} else {

			// Keep only {subsample:proportion} of initial points.
			auto iterator_to_first_to_remove = CGAL::random_simplify_point_set (point_cloud, (1 - subsample) * 100);
			point_cloud.remove(iterator_to_first_to_remove, point_cloud.end());
			point_cloud.collect_garbage();

			min_point_per_area *= subsample;

			std::cout << "Points subsampled with " << subsample << " points of original point cloud." << std::endl;

		}

	}
	
	const LindstromTurk_param params (l1,l2,l3,l4,l5,l6,l7);
	Custom_placement pf(params, mesh, point_cloud, ablation);
	Custom_cost cf(params, c1, c2, c3, c4, min_point_per_area, mesh, point_cloud, next_mesh);
	My_visitor mv(params, c1, c2, c3, min_point_per_area, mesh, mesh_info, point_cloud, ablation);
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
