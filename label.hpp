#ifndef LABEL_H_
#define LABEL_H_

#include <string>
#include <array>

//Digitanie COS43 labels

struct Label {
    unsigned char value;
    std::string label;
    unsigned char red;
    unsigned char green;
    unsigned char blue;
    Label(unsigned char value, std::string label, unsigned char red, unsigned char green, unsigned char blue): value(value), label(label), red(red), green(green), blue(blue) {}
};

static std::array<Label, 45> LABELS {
    Label(0, "other", 255, 255, 255),
    Label(1, "construction_site", 140, 90, 100),
    Label(2, "bare_ground", 71, 58, 17),
    Label(3, "bare_parking", 160, 130, 105),
    Label(4, "bare_road", 160, 119, 34),
    Label(5, "bare_pedestrian", 230, 167, 31),
    Label(6, "sand", 234, 220, 90),
    Label(7, "snow", 240, 239, 220),
    Label(8, "field", 235, 255, 6),
    Label(9, "sport_vegetation", 190, 215, 165),
    Label(10, "grassland", 140, 240, 118),
    Label(11, "aquaculture", 11, 222, 189),
    Label(12, "hedge", 119, 211, 0),
    Label(13, "shrub", 113, 184, 48),
    Label(14, "vegetation", 0, 210, 50),
    Label(15, "arboriculture", 120, 155, 100),
    Label(16, "tree", 0, 120, 15),
    Label(17, "spiney", 30, 143, 100),
    Label(18, "forest", 0, 70, 0),
    Label(19, "winter_high_vegetation", 90, 180, 170),
    Label(20, "ice", 59, 173, 240),
    Label(21, "river", 30, 145, 246),
    Label(22, "pond", 0, 75, 190),
    Label(23, "sea", 0, 55, 105),
    Label(24, "swimmingpool", 100, 130, 255),
    Label(25, "bridge", 160, 160, 246),
    Label(26, "boat", 130, 41, 244),
    Label(27, "railways", 75, 20, 132),
    Label(28, "road", 141, 91, 210),
    Label(29, "private_road", 205, 140, 242),
    Label(30, "central_reservation", 163, 127, 180),
    Label(31, "parking", 170, 60, 160),
    Label(32, "pedestrian", 190, 38, 194),
    Label(33, "yard", 255, 175, 200),
    Label(34, "sport", 255, 115, 180),
    Label(35, "cemetery", 125, 15, 70),
    Label(36, "impervious", 219, 20, 123),
    Label(37, "terrace", 200, 20, 79),
    Label(38, "container", 255, 65, 75),
    Label(39, "storage_tank", 195, 85, 0),
    Label(40, "greenhouse", 255, 150, 85),
    Label(41, "building", 240, 0, 0),
    Label(42, "high_building", 127, 1, 0),
    Label(43, "pipeline", 35, 85, 85),
    Label(44, "unknown", 250, 250, 150)
};

const unsigned char LABEL_OTHER = 44;

const unsigned char LABEL_RAIL = 27;
const unsigned char LABEL_WATER = 21;
const unsigned char LABEL_ROAD = 28;

#endif  /* !LABEL_H_ */
