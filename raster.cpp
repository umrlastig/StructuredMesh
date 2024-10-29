#include "raster.hpp"

#include <gdal_priv.h>

std::pair<int,int> Raster::grid_conversion(int P, int L, double grid_to_crs[6], OGRCoordinateTransformation *crs_to_other_crs, double other_grid_to_other_crs[6]) {
    double x = grid_to_crs[0] + (0.5 + P)*grid_to_crs[1] + (0.5 + L)*grid_to_crs[2];
    double y = grid_to_crs[3] + (0.5 + P)*grid_to_crs[4] + (0.5 + L)*grid_to_crs[5];
    crs_to_other_crs->Transform(1,&x,&y);
    double fact = other_grid_to_other_crs[1]*other_grid_to_other_crs[5] - other_grid_to_other_crs[2]*other_grid_to_other_crs[4];
    x -= other_grid_to_other_crs[0];
    y -= other_grid_to_other_crs[3];
    int newP = ((int) ((other_grid_to_other_crs[5]*x - other_grid_to_other_crs[2]*y) / fact));
    int newL = ((int) ((-other_grid_to_other_crs[4]*x + other_grid_to_other_crs[1]*y) / fact));
    return std::pair<int,int>(newP, newL);
}

void Raster::align_land_cover_and_dsm() {
    std::vector<std::vector<unsigned char>> new_land_cover = std::vector<std::vector<unsigned char>>(ySize, std::vector<unsigned char>(xSize, LABEL_UNKNOWN));
    for (int L = 0; L < ySize; L++) {
        for (int P = 0; P < xSize; P++) {
            float neighboors_count[LABELS.size()] = {0};
            int n = 0;
            float moy = 0;
            float square_moy = 0;
            for (int i = -5; i <= 5; i++) {
                for (int j = -5; j <= 5; j++) {
                    if (L+i >= 0 && P + j >= 0 && L+i < ySize && P+j < xSize) {
                        neighboors_count[land_cover[L+i][P+j]] += 1/(0.1+abs(dsm[L][P] - dsm[L+i][P+j]));
                        if (i > -3 && i < 3 && j > -3 && j < 3) {
                            n++;
                            moy += dsm[L+i][P+j];
                            square_moy += pow(dsm[L+i][P+j], 2);
                        }
                    }
                }
            }
            if (pow(square_moy/n - pow(moy/n, 2), 0.5) > 1) {
                new_land_cover[L][P] = std::max_element(neighboors_count, neighboors_count + LABELS.size()) - neighboors_count;
            } else {
                new_land_cover[L][P] = land_cover[L][P];
            }

        }
    }

    land_cover = new_land_cover;
}

void Raster::coord_to_grid(double x, double y, float& P, float& L) const {
    double fact = grid_to_crs[1]*grid_to_crs[5] - grid_to_crs[2]*grid_to_crs[4];
    x -= grid_to_crs[0];
    y -= grid_to_crs[3];
    P = (grid_to_crs[5]*x - grid_to_crs[2]*y) / fact;
    L = (-grid_to_crs[4]*x + grid_to_crs[1]*y) / fact;
}

void Raster::grid_to_coord(int P, int L, double& x, double& y) const {
    x = grid_to_crs[0] + (0.5 + P)*grid_to_crs[1] + (0.5 + L)*grid_to_crs[2];
    y = grid_to_crs[3] + (0.5 + P)*grid_to_crs[4] + (0.5 + L)*grid_to_crs[5];
}

void Raster::grid_to_coord(K::FT P_, K::FT L_, double& x, double& y) const {
    double P = CGAL::to_double(P_);
    double L = CGAL::to_double(L_);

    x = grid_to_crs[0] + P*grid_to_crs[1] + L*grid_to_crs[2];
    y = grid_to_crs[3] + P*grid_to_crs[4] + L*grid_to_crs[5];
}

float Raster::coord_distance_to_grid_distance(double d) const {
    float P, L;
    double x, y;
    grid_to_coord(0, 0, x, y);
    coord_to_grid(x+d/sqrt(2),y+d/sqrt(2),P,L);
    return sqrt(pow(P,2)+pow(L,2));
}

double Raster::grid_distance_to_coord_distance(K::FT d_) const {
    double d = CGAL::to_double(d_);
    double x1, y1, x2, y2;
    grid_to_coord(0, 0, x1, y1);
    grid_to_coord((float) (d/sqrt(2)), (float) (d/sqrt(2)), x2, y2);
    return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
}

OGRSpatialReference Raster::get_crs() const {
    return OGRSpatialReference(crs);
}

void Raster::fill_holes() {
    unsigned int max_hole_size = 300;

    //Negative holes
    {
        // Find holes
        std::list<std::pair<float,std::list<std::pair<int,int>>>> holes;
        for (int L = 1; L < ySize - 1; L++) {
            for (int P = 1; P < xSize -1; P++) {
                if (land_cover[L][P] == LABEL_RAIL || land_cover[L][P] == LABEL_ROAD || land_cover[L][P] == LABEL_WATER) {
                    if ((dsm[L][P] <= dsm[L+1][P]) && (dsm[L][P] <= dsm[L-1][P]) && (dsm[L][P] <= dsm[L][P+1]) && (dsm[L][P] <= dsm[L][P-1])) {
                        holes.push_back(std::pair<float,std::list<std::pair<int,int>>>(dsm[L][P], {std::pair<int,int>(P,L)}));
                    }
                }
            }
        }

        // Compute hole elevation
        for (auto &hole: holes) {
            bool end = false;
            std::list<std::pair<int,int>> ring ({
                std::pair<int,int>(hole.second.front().first-1,hole.second.front().second),
                std::pair<int,int>(hole.second.front().first+1,hole.second.front().second),
                std::pair<int,int>(hole.second.front().first,hole.second.front().second-1),
                std::pair<int,int>(hole.second.front().first,hole.second.front().second+1)
                });

            while(!end) {
                float min_dsm = FLT_MAX;
                std::pair<int,int> min_pixel;

                for (auto pixel: ring) {
                    if (dsm[pixel.second][pixel.first] < min_dsm) {
                        min_pixel = pixel;
                        min_dsm = dsm[pixel.second][pixel.first];
                    }
                }

                hole.first = min_dsm;
                hole.second.push_back(min_pixel);
                ring.remove(min_pixel);

                if (min_pixel.first-1 >= 0) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first-1, min_pixel.second)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first-1, min_pixel.second)) == ring.end()) {
                            if (dsm[min_pixel.second][min_pixel.first-1] < hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(min_pixel.first-1, min_pixel.second));
                            }
                        }
                    }
                }
                if (min_pixel.first+1 < xSize) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first+1, min_pixel.second)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first+1, min_pixel.second)) == ring.end()) {
                            if (dsm[min_pixel.second][min_pixel.first+1] < hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(min_pixel.first+1, min_pixel.second));
                            }
                        }
                    }
                }
                if (min_pixel.second-1 >= 0) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first, min_pixel.second-1)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first, min_pixel.second-1)) == ring.end()) {
                            if (dsm[min_pixel.second-1][min_pixel.first] < hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(min_pixel.first, min_pixel.second-1));
                            }
                        }
                    }
                }
                if (min_pixel.second+1 < ySize) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(min_pixel.first, min_pixel.second+1)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(min_pixel.first, min_pixel.second+1)) == ring.end()) {
                            if (dsm[min_pixel.second+1][min_pixel.first] < hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(min_pixel.first, min_pixel.second+1));
                            }
                        }
                    }
                }

                if (hole.second.size() >= max_hole_size) {
                    end = true;
                }
            }

            //Fill holes
            for (auto pixel: hole.second) {
                dsm[pixel.second][pixel.first] = hole.first;
            }
        }
    }

    //Positive holes
    {
        // Find holes
        std::list<std::pair<float,std::list<std::pair<int,int>>>> holes;
        for (int L = 1; L < ySize - 1; L++) {
            for (int P = 1; P < xSize -1; P++) {
                if (land_cover[L][P] == LABEL_RAIL || land_cover[L][P] == LABEL_ROAD || land_cover[L][P] == LABEL_WATER) {
                    if ((dsm[L][P] >= dsm[L+1][P]) && (dsm[L][P] >= dsm[L-1][P]) && (dsm[L][P] >= dsm[L][P+1]) && (dsm[L][P] >= dsm[L][P-1])) {
                        holes.push_back(std::pair<float,std::list<std::pair<int,int>>>(dsm[L][P], {std::pair<int,int>(P,L)}));
                    }
                }
            }
        }

        // Compute hole elevation
        for (auto &hole: holes) {
            bool end = false;
            std::list<std::pair<int,int>> ring ({
                std::pair<int,int>(hole.second.front().first-1,hole.second.front().second),
                std::pair<int,int>(hole.second.front().first+1,hole.second.front().second),
                std::pair<int,int>(hole.second.front().first,hole.second.front().second-1),
                std::pair<int,int>(hole.second.front().first,hole.second.front().second+1)
                });

            while(!end) {
                float max_dsm = FLT_MIN;
                std::pair<int,int> max_pixel;

                for (auto pixel: ring) {
                    if (dsm[pixel.second][pixel.first] > max_dsm) {
                        max_pixel = pixel;
                        max_dsm = dsm[pixel.second][pixel.first];
                    }
                }

                hole.first = max_dsm;
                hole.second.push_back(max_pixel);
                ring.remove(max_pixel);

                if (max_pixel.first-1 >= 0) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first-1, max_pixel.second)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first-1, max_pixel.second)) == ring.end()) {
                            if (dsm[max_pixel.second][max_pixel.first-1] > hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(max_pixel.first-1, max_pixel.second));
                            }
                        }
                    }
                }
                if (max_pixel.first+1 < xSize) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first+1, max_pixel.second)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first+1, max_pixel.second)) == ring.end()) {
                            if (dsm[max_pixel.second][max_pixel.first+1] > hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(max_pixel.first+1, max_pixel.second));
                            }
                        }
                    }
                }
                if (max_pixel.second-1 >= 0) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first, max_pixel.second-1)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first, max_pixel.second-1)) == ring.end()) {
                            if (dsm[max_pixel.second-1][max_pixel.first] > hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(max_pixel.first, max_pixel.second-1));
                            }
                        }
                    }
                }
                if (max_pixel.second+1 < ySize) {
                    if (std::find(hole.second.begin(), hole.second.end(), std::pair<int,int>(max_pixel.first, max_pixel.second+1)) == hole.second.end()) {
                        if (std::find(ring.begin(), ring.end(), std::pair<int,int>(max_pixel.first, max_pixel.second+1)) == ring.end()) {
                            if (dsm[max_pixel.second+1][max_pixel.first] > hole.first) {
                                end = true;
                            } else {
                                ring.push_back(std::pair<int,int>(max_pixel.first, max_pixel.second+1));
                            }
                        }
                    }
                }

                if (hole.second.size() >= max_hole_size) {
                    end = true;
                }
            }

            //Fill holes
            for (auto pixel: hole.second) {
                dsm[pixel.second][pixel.first] = hole.first;
            }
        }
    }
}

Raster::Raster(char *dsm_path, char *dtm_path, char *land_cover_path) {
    GDALAllRegister();

    // Get DSM informations and CRS
    GDALDataset *dsm_dataset = (GDALDataset *) GDALOpen(dsm_path, GA_ReadOnly );
    if( dsm_dataset == NULL ) {
        throw std::invalid_argument(std::string("Unable to open ") + dsm_path + ".");
    }
    xSize = dsm_dataset->GetRasterBand(1)->GetXSize();
    ySize = dsm_dataset->GetRasterBand(1)->GetYSize();
    crs = *dsm_dataset->GetSpatialRef();
    if (dsm_dataset->GetGeoTransform(grid_to_crs) >= CE_Failure) {
        throw std::invalid_argument(std::string(dsm_path) + " do not contain an affine transform.");
    }
    dsm = std::vector<std::vector<float>>(ySize, std::vector<float>(xSize, 0));
    for (int L = 0; L < ySize; L++) {
        if (dsm_dataset->GetRasterBand(1)->RasterIO(GF_Read, 0, L, xSize, 1, &dsm[L][0], xSize, 1, GDT_Float32, 0, 0) >= CE_Failure) {
            throw std::invalid_argument(std::string(dsm_path) + " can't be read.");
        }
    }
    std::cout << "DSM load" << std::endl;

    // Get DTM informations and CRS transform
    dtm = std::vector<std::vector<float>>(ySize, std::vector<float>(xSize, 0));
    GDALDataset *dtm_dataset = (GDALDataset *) GDALOpen( dtm_path, GA_ReadOnly );
    if( dtm_dataset == NULL ) {
        throw std::invalid_argument(std::string("Unable to open ") + dtm_path + ".");
    }
    double dtm_grid_to_dtm_crs[6];
    if (dtm_dataset->GetGeoTransform(dtm_grid_to_dtm_crs) >= CE_Failure) {
        throw std::invalid_argument("Can't transform DTM grid to DTM CRS.");
    }
    OGRCoordinateTransformation *CRS_to_dtm_crs = OGRCreateCoordinateTransformation(
        &crs,
        dtm_dataset->GetSpatialRef());
    if (CRS_to_dtm_crs == NULL) {
        throw std::runtime_error("Can't transform DSM CRS to DTM CRS.");
    }
    for (int L = 0; L < ySize; L++) {
        for (int P = 0; P < xSize; P++) {
            std::pair<int,int> new_coord = grid_conversion(P, L, grid_to_crs, CRS_to_dtm_crs, dtm_grid_to_dtm_crs);
            if (new_coord.first >= 0 && new_coord.first < dtm_dataset->GetRasterBand(1)->GetXSize() && new_coord.second >= 0 && new_coord.second < dtm_dataset->GetRasterBand(1)->GetYSize()) {
                if (dtm_dataset->GetRasterBand(1)->RasterIO(GF_Read, new_coord.first, new_coord.second, 1, 1, &dtm[L][P], 1, 1, GDT_Float32, 0, 0) >= CE_Failure) {
                    throw std::invalid_argument(std::string(dtm_path) + " can't be read.");
                }
            }
        }
    }
    std::cout << "DTM load" << std::endl;

    // Get land cover informations and CRS transform
    land_cover = std::vector<std::vector<unsigned char>>(ySize, std::vector<unsigned char>(xSize, LABEL_UNKNOWN));
    GDALDataset *land_cover_dataset = (GDALDataset *) GDALOpen( land_cover_path, GA_ReadOnly );
    if( dsm_dataset == NULL ) {
        throw std::invalid_argument(std::string("Unable to open ") + land_cover_path + ".");
    }

    double land_cover_grid_to_land_cover_crs[6];
    if (land_cover_dataset->GetGeoTransform(land_cover_grid_to_land_cover_crs) >= CE_Failure) {
        throw std::invalid_argument("Can't transform land cover grid to land cover CRS.");
    }
    OGRCoordinateTransformation *CRS_to_land_cover_crs = OGRCreateCoordinateTransformation(
        &crs,
        land_cover_dataset->GetSpatialRef());
    if (CRS_to_land_cover_crs == NULL) {
        throw std::runtime_error("Can't transform DSM CRS to land cover CRS.");
    }
    for (int L = 0; L < ySize; L++) {
        for (int P = 0; P < xSize; P++) {
            std::pair<int,int> new_coord = grid_conversion(P, L, grid_to_crs, CRS_to_land_cover_crs, land_cover_grid_to_land_cover_crs);
            if (new_coord.first >= 0 && new_coord.first < land_cover_dataset->GetRasterBand(1)->GetXSize() && new_coord.second >= 0 && new_coord.second < land_cover_dataset->GetRasterBand(1)->GetYSize()) {
                unsigned char value;
                if (land_cover_dataset->GetRasterBand(1)->RasterIO(GF_Read, new_coord.first, new_coord.second, 1, 1, &value, 1, 1, GDT_Byte, 0, 0) >= CE_Failure) {
                    throw std::invalid_argument(std::string(land_cover_path) + " can't be read.");
                }
                if (value >= LABELS.size()) value = LABEL_UNKNOWN;
                land_cover[L][P] = value;
            }
        }
    }
    std::cout << "Land cover load" << std::endl;

    align_land_cover_and_dsm();

    fill_holes();
}
