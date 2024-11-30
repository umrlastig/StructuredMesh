# Compute LOD2
Create a Structured mesh from a mesh and a labeled point cloud.

# Requirements (test version)
- CGAL (5.4.5)
- GDAL (3.5.0)
- Eigen (3.4.0)

# Compilation
```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

# Usage
Usage: `./structured-mesh` [OPTIONS] -m MESH -p POINT_CLOUD

OPTIONS:
- `-h`, `--help`: Print this help anq quit.
- `-m`, `--mesh=/file/path.ply`: mesh as PLY file.
- `-p`, `--point_cloud=/file/path.ply`: point cloud as PLY file.

See main_edge_collapse.cpp for other options.