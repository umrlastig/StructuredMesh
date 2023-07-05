#ifndef LABEL_H_
#define LABEL_H_

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

static std::array<Label, 12> LABELS {
    Label(0, "other", 255, 255, 255),
    Label(1, "bare ground", 100, 50, 0),
    Label(2, "low vegetation", 0, 250, 50),
    Label(3, "water", 0, 50, 250),
    Label(4, "building", 250, 50, 50),
    Label(5, "high vegetation", 0, 100, 50),
    Label(6, "parking", 200, 200, 200),
    Label(7, "pedestrian", 200, 150, 50),
    Label(8, "road", 100, 100, 100),
    Label(9, "railways", 200, 100, 200),
    Label(10, "swimming pool", 50, 150, 250),
    Label(11, "rail crossing", 250, 150, 0)
};

const int LABEL_OTHER = 0;

#endif  /* !LABEL_H_ */
