#ifndef PARTICLESHAPE2D_H
#define PARTICLESHAPE2D_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <common_tools/constants.h>

template <class PointType>
class gridSteps
{
  public:
    static PointType H1;
    static PointType H2;
};
template <class PointType>
PointType gridSteps<PointType>::H2;
template <class PointType>
PointType gridSteps<PointType>::H1;

template <class PointType>
class ParticleShapeCIC2dStatic
{
  public:
    PointType   CellArea();
    static bool InCell(PointType x1, PointType x2, PointType x3,
                       const std::vector<PointType>& x1Array, const std::vector<PointType>& x2Array,
                       const std::vector<PointType>& x3Array, unsigned int basePoint,
                       const std::vector<int>& levelHigh, const std::vector<int>& levelZ);
    static bool InCellWithEps(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                              const std::vector<PointType>& x2Array, unsigned int basePoint,
                              unsigned int levelHigh);
    static bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                       const std::vector<PointType>& x2Array, unsigned int basePoint,
                       unsigned int levelHigh);
    static void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                           const std::vector<PointType>& x2Array, std::vector<PointType>& W,
                           unsigned int basePoint, unsigned int levelHigh);
    static void Wcalculate3d(PointType x1, PointType x2, PointType x3,
                             const std::vector<PointType>& x1Array,
                             const std::vector<PointType>& x2Array,
                             const std::vector<PointType>& x3Array, std::vector<PointType>& W,
                             unsigned int basePoint, const std::vector<int>& levelHigh,
                             const std::vector<int>& levelZ);
    static void WcalculatePolar(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                                const std::vector<PointType>& x2Array, std::vector<PointType>& W,
                                unsigned int basePoint, unsigned int levelHigh);

    static void ValueInterpolate(const std::vector<PointType>& W,
                                 const std::vector<PointType>& ValueArray, PointType& result,
                                 unsigned int basePoint, unsigned int levelHigh);
    static void ValueInterpolate3d(const std::vector<PointType>& W,
                                   const std::vector<PointType>& ValueArray, PointType& result,
                                   unsigned int basePoint, const std::vector<int>& levelHigh,
                                   const std::vector<int>& levelZ);
    static void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                                std::vector<PointType>& rho, unsigned int basePoint,
                                unsigned int levelHigh);
    static void ChargeZeroing();
    static void CellVolume(const std::vector<PointType>& rArray,
                           const std::vector<PointType>& zArray, unsigned int basePoint,
                           unsigned int levelHigh, PointType& vol);
};

#endif
