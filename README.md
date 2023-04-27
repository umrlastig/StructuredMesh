# Compute LOD2
Create LOD2 model from DSM, DTM, land use map, LOD0 building footprint and orthophoto.

# Requirements (test version)
- CGAL (5.4.1)
- GDAL (3.5.0)
- Eigen (3.4.0)
- Ceres (2.1.0)

# Compilation
```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ make
```

# Usage
Usage: `./compute-LOD2` [OPTIONS] -s DSM -t DTM -l land_use_map

OPTIONS:
- `-h`, `--help`: Print this help anq quit.
- `-s`, `--DSM=/file/path.tiff`: DSM as TIFF file.
- `-t`, `--DTM=/file/path.tiff`: DTM as TIFF file.
- `-l`, `--land_use_map=/file/path.tiff`: land use map as TIFF file.
- `-0`, `--LOD0=/file/path.shp`: LOD0 as Shapefile.
- `-i`, `--orthophoto=/file/path.tiff`: RGB orthophoto as TIFF file.
