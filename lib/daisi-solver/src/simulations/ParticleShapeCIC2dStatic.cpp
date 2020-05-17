#include "ParticleShape2d.h"

template class ParticleShapeCIC2dStatic<float>;
template class ParticleShapeCIC2dStatic<double>;

template <class PointType>
PointType ParticleShapeCIC2dStatic<PointType>::CellArea()
{
    return gridSteps<PointType>::H1 * gridSteps<PointType>::H2;
};
template <class PointType>
bool ParticleShapeCIC2dStatic<PointType>::InCell(PointType x1, PointType x2,
                                                 const std::vector<PointType>& x1Array,
                                                 const std::vector<PointType>& x2Array,
                                                 unsigned int basePoint, unsigned int levelHigh)
{
    if (levelHigh == -1)
        return false;
    if (x1 >= x1Array[basePoint] && x1 <= x1Array[basePoint + 1] && x2 >= x2Array[basePoint] &&
        x2 <= x2Array[levelHigh])
        return true;
    return false;
};

template <class PointType>
bool ParticleShapeCIC2dStatic<PointType>::InCell(
    PointType x1, PointType x2, PointType x3, const std::vector<PointType>& x1Array,
    const std::vector<PointType>& x2Array, const std::vector<PointType>& x3Array,
    unsigned int basePoint, const std::vector<int>& levelHigh, const std::vector<int>& levelZ)
{
    if (levelHigh[basePoint] == -1)
        return false;
    if (levelZ[basePoint] == -1)
        return false;

    if (x1 >= x1Array[basePoint] && x1 <= x1Array[basePoint + 1] && x2 >= x2Array[basePoint] &&
        x2 <= x2Array[levelHigh[basePoint]] && x3 >= x3Array[basePoint] &&
        x3 <= x3Array[levelZ[basePoint]])
        return true;
    return false;
};

template <class PointType>
bool ParticleShapeCIC2dStatic<PointType>::InCellWithEps(PointType x1, PointType x2,
                                                        const std::vector<PointType>& x1Array,
                                                        const std::vector<PointType>& x2Array,
                                                        unsigned int                  basePoint,
                                                        unsigned int                  levelHigh)
{
    if (levelHigh == -1)
        return false;

    PointType epsx1 = std::abs(1e-2 * (x1Array[basePoint + 1] - x1Array[basePoint]));
    PointType epsx2 = std::abs(1e-2 * (x2Array[levelHigh] - x2Array[basePoint]));

    if (levelHigh == -1)
        return false;
    if ((x1 >= x1Array[basePoint] || std::abs(x1 - x1Array[basePoint]) < epsx1) &&
        (x1 <= x1Array[basePoint + 1] || std::abs(x1 - x1Array[basePoint + 1]) < epsx1) &&
        (x2 >= x2Array[basePoint] || std::abs(x2 - x2Array[basePoint]) < epsx2) &&
        (x2 <= x2Array[levelHigh] || std::abs(x2 - x2Array[levelHigh]) < epsx2))
        return true;
    return false;
};
template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::Wcalculate(PointType x1, PointType x2,
                                                     const std::vector<PointType>& x1Array,
                                                     const std::vector<PointType>& x2Array,
                                                     std::vector<PointType>&       W,
                                                     unsigned int basePoint, unsigned int levelHigh)
{
    PointType H1   = x1Array[basePoint + 1] - x1Array[basePoint];
    PointType H2   = x2Array[levelHigh] - x2Array[basePoint];
    PointType tt   = x1Array[basePoint + 1];
    PointType wx11 = (x1Array[basePoint + 1] - x1) / H1;
    PointType wx12 = 1 - wx11;
    PointType wx21 = (x2Array[levelHigh] - x2) / H2;
    PointType wx22 = 1 - wx21;
    W[0]           = wx11 * wx21;
    W[1]           = wx12 * wx21;
    W[2]           = wx11 * wx22;
    W[3]           = wx12 * wx22;
};

template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::Wcalculate3d(
    PointType x1, PointType x2, PointType x3, const std::vector<PointType>& x1Array,
    const std::vector<PointType>& x2Array, const std::vector<PointType>& x3Array,
    std::vector<PointType>& W, unsigned int basePoint, const std::vector<int>& levelHigh,
    const std::vector<int>& levelZ)
{
    PointType H1 = x1Array[basePoint + 1] - x1Array[basePoint];
    PointType H2 = x2Array[levelHigh[basePoint]] - x2Array[basePoint];
    PointType H3 = x3Array[levelZ[basePoint]] - x3Array[basePoint];

    PointType wx11 = (x1Array[basePoint + 1] - x1) / H1;
    PointType wx12 = 1 - wx11;

    PointType wx21 = (x2Array[levelHigh[basePoint]] - x2) / H2;
    PointType wx22 = 1 - wx21;

    PointType wx31 = (x3Array[levelZ[basePoint]] - x3) / H3;
    PointType wx32 = 1 - wx31;

    W[0] = wx11 * wx21 * wx31;
    W[1] = wx12 * wx21 * wx31;
    W[2] = wx11 * wx22 * wx31;
    W[3] = wx12 * wx22 * wx31;

    W[4] = wx11 * wx21 * wx32;
    W[5] = wx12 * wx21 * wx32;
    W[6] = wx11 * wx22 * wx32;
    W[7] = wx12 * wx22 * wx32;
};

template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::WcalculatePolar(PointType x1, PointType x2,
                                                          const std::vector<PointType>& x1Array,
                                                          const std::vector<PointType>& x2Array,
                                                          std::vector<PointType>&       W,
                                                          unsigned int                  basePoint,
                                                          unsigned int                  levelHigh)
{

    PointType vol1 =
        (x1 * x1 - x1Array[basePoint] * x1Array[basePoint]) * (x2 - x2Array[basePoint]) / 2;

    PointType vol2 =
        (x1 * x1 - x1Array[basePoint] * x1Array[basePoint]) * (x2Array[levelHigh] - x2) / 2;

    PointType vol3 =
        (x1Array[basePoint + 1] * x1Array[basePoint + 1] - x1 * x1) * (x2 - x2Array[basePoint]) / 2;

    PointType vol4 =
        (x1Array[basePoint + 1] * x1Array[basePoint + 1] - x1 * x1) * (x2Array[levelHigh] - x2) / 2;

    PointType V = vol1 + vol2 + vol3 + vol4;

    W[0] = vol4 / V;
    W[1] = vol2 / V;
    W[2] = vol3 / V;
    W[3] = vol1 / V;
};

template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::ValueInterpolate(const std::vector<PointType>& W,
                                                           const std::vector<PointType>& ValueArray,
                                                           PointType&                    result,
                                                           unsigned int                  basePoint,
                                                           unsigned int                  levelHigh)
{
    result = W[0] * ValueArray[basePoint] + W[1] * ValueArray[basePoint + 1] +
             W[2] * ValueArray[levelHigh] + W[3] * ValueArray[levelHigh + 1];
};

template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::ValueInterpolate3d(
    const std::vector<PointType>& W, const std::vector<PointType>& ValueArray, PointType& result,
    unsigned int basePoint, const std::vector<int>& levelHigh, const std::vector<int>& levelZ)
{
    result = W[0] * ValueArray[basePoint] + W[1] * ValueArray[basePoint + 1] +
             W[2] * ValueArray[levelHigh[basePoint]] + W[3] * ValueArray[levelHigh[basePoint] + 1] +
             W[4] * ValueArray[levelZ[basePoint]] + W[5] * ValueArray[levelZ[basePoint] + 1] +
             W[6] * ValueArray[levelZ[levelHigh[basePoint]]] +
             W[7] * ValueArray[levelZ[levelHigh[basePoint]] + 1];
};

template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::ChargeCalculate(const std::vector<PointType>& W,
                                                          PointType               particleCharge,
                                                          std::vector<PointType>& rho,
                                                          unsigned int            basePoint,
                                                          unsigned int            levelHigh)
{
    rho[basePoint]     = rho[basePoint] + W[0] * particleCharge;
    rho[basePoint + 1] = rho[basePoint + 1] + W[1] * particleCharge;
    rho[levelHigh]     = rho[levelHigh] + W[2] * particleCharge;
    rho[levelHigh + 1] = rho[levelHigh + 1] + W[3] * particleCharge;
};
template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::ChargeZeroing(){};
template <class PointType>
void ParticleShapeCIC2dStatic<PointType>::CellVolume(const std::vector<PointType>& rArray,
                                                     const std::vector<PointType>& zArray,
                                                     unsigned int basePoint, unsigned int levelHigh,
                                                     PointType& vol)
{
    PointType H1 = rArray[basePoint + 1] - rArray[basePoint];
    PointType H2 = zArray[levelHigh] - zArray[basePoint];
    vol          = 1;
    vol          = PointType(PI() * gridSteps<PointType>::H2 *
                    ((rArray[basePoint] + H1 / 2) * (rArray[basePoint] + H1 / 2) -
                     (rArray[basePoint] - H1 / 2) * (rArray[basePoint] - H1 / 2)));
};
