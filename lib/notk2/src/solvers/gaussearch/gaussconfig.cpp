#include "gaussconfig.h"

namespace notk
{
bool GaussSearchConfig::get_use_random()
{
    return use_random;
}
size_t GaussSearchConfig::get_max_preprocess_iterations()
{
    return maximal_preprocess_iterations;
}
}