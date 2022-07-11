#include <string>
#include <array>

//Digitanie COS11 labels

struct Label {
    unsigned char value;
    std::string label;
    unsigned char red;
    unsigned char green;
    unsigned char blue;
    Label(int value, std::string label, unsigned char red, unsigned char green, unsigned char blue): value(value), label(label), red(red), green(green), blue(blue) {}
};

std::array<Label, 11> LABELS {
    Label(0, "other", 0, 0, 0),
    Label(1, "bare ground", 146, 109, 39),
    Label(2, "low vegetation", 131, 160, 110),
    Label(3, "water", 40, 67, 96),
    Label(4, "building", 179, 130, 102),
    Label(5, "high vegetation", 134, 206, 108),
    Label(6, "parking", 188, 183, 174),
    Label(7, "pedestrian", 180, 71, 53),
    Label(8, "road", 80, 89, 85),
    Label(9, "railways", 170, 106, 96),
    Label(10, "swimming pool", 105, 185, 211)
};
