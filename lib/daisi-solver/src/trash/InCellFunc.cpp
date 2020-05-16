#include "ParticleGridInterface2d.h"

template class ParticleGridInterface2d<float, std::vector<float>>;
template class ParticleGridInterface2d<double, std::vector<double>>;

template <class PointType, class DataContainer>
bool ParticleGridInterface2d<PointType, DataContainer>::InCell(PointType x1, PointType x2,
                                                               unsigned int basePoint,
                                                               unsigned int levelHigh)
{
    if (x1 >= x1Array[basePoint] && x1 <= x1Array[basePoint + 1] && x2 >= x2Array[basePoint] &&
        x2 <= x2Array[levelHigh])
        return true;
    return false;
};
template <class PointType, class DataContainer>

/*bool <PointType, DataContainer>::InCell1(PointType x1, PointType x2, )
{
        PointType tmp1 = x1Array[p0];
        PointType tmp2 = x1Array[p1];
        PointType tmp3 = x1Array[p2];
        PointType tmp4 = x1Array[basePoint];
        if (x1 >= x1Array[p0] && x1 <= x1Array[p1] && x2 >= x2Array[p0] && x2 <= x2Array[p2] &&
std::abs(x1 - x1Array[basePoint])<gridSteps<PointType>::H1 / 2 && std::abs(x2 -
x2Array[basePoint])<gridSteps<PointType>::H2 / 2) return true; return false;
};*/