#include "colors.h"
namespace colors
{
const float                     whitecolor[3] = {255, 255, 255};
const float                     redcolor[3]   = {1.0, 0.0, 0.0};
const float                     greencolor[3] = {0.0, 1.0, 0.0};
const float                     blackcolor[3] = {0.0, 0.0, 0.0};
const float                     bluecolor[3]  = {0.0, 0.0, 1.0};
const std::vector<const float*> colors        = {redcolor, greencolor, bluecolor};
const std::vector<const float*> colorsPoints  = {blackcolor, greencolor, bluecolor};
}