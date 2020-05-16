#ifndef PARTICLESHAPE2DCIC_H
#define PARTICLESHAPE2DCIC_H

#include "armadillo"

#include <memory>
#include <vector>

#include <Dmath.h>

template <class PointType>
class BoundaryContainer2d;
template <class PointType>
class BoundaryContainer3d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class GridData2d;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;
template <class PointType>
class GridData3dpolar;
template <class PointType>
class NearCathodeVolume;
template <class PointType>
class ParticleShape2dCIC
{
  public:
    std::vector<int> levelHigh;
    std::vector<int> levelZ;

    std::vector<std::vector<unsigned int>> emCells;

    std::vector<char>      isBoundary;
    std::vector<PointType> cellVolume;

    const PointType* x1Array;
    const PointType* x2Array;
    const PointType* x3Array;

    unsigned int cellNumber;

    long long GetMemorySize()
    {
        long long result = 0;
        result           = result + Dmath::vectorsize(emCells);
        result           = result + Dmath::vectorsize(levelHigh);
        result           = result + Dmath::vectorsize(levelZ);
        result           = result + Dmath::vectorsize(isBoundary);
        result           = result + Dmath::vectorsize(cellVolume);
        return result;
    };

    ParticleShape2dCIC()
    {
        x1Array = NULL;
        x2Array = NULL;
        x3Array = NULL;
    };

    int
    InitEmCells(const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumeS,
                PointType r1, PointType z1, int index, int flagInit);

    void init(int size, const PointType* x1ArrayIn, const PointType* x2ArrayIn);
    void init(int size, const PointType* x1ArrayIn, const PointType* x2ArrayIn,
              const PointType* x3ArrayIn);

    void AddCell(const std::shared_ptr<GridData2daxs<PointType>>& gridData,
                 unsigned int levelHighIn, const std::vector<int>& flagOut,
                 std::vector<int>                                      pointsIn,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);
    void AddCell(const std::shared_ptr<GridData2d<PointType>>& gridData, unsigned int levelHighIn,
                 const std::vector<int>& flagOut, std::vector<int> pointsIn,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);
    void AddCell(const std::shared_ptr<GridData2dpolar<PointType>>& gridData,
                 unsigned int levelHighIn, const std::vector<int>& flagOut,
                 std::vector<int>                                      pointsIn,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);
    void AddCell(const std::shared_ptr<GridData3d<PointType>>& gridData, unsigned int levelHighIn,
                 unsigned int levelZIn, const std::vector<int>& flagOut, std::vector<int> pointsIn,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);

    int InCell(PointType x1, PointType x2);
    int InCell(PointType x1, PointType x2, PointType x3);
    bool InCell(int index, PointType x1, PointType x2);
    bool InCell(int index, PointType x1, PointType x2, PointType x3);

    int InCellWithEps(PointType x1, PointType x2);
    int InCellWithEps(PointType x1, PointType x2, PointType x3);
    bool InCellWithEps(int index, PointType x1, PointType x2);
    bool InCellWithEps(int index, PointType x1, PointType x2, PointType x3);

    void Wcalculate(unsigned int* basePoint, PointType* x1, PointType* x2, PointType (*W)[9],
                    int i1, int i2, PointType* p1, PointType* p2, PointType* tmp1, PointType* tmp2,
                    PointType* gamma, double& size);

    void Wcalculate(int index, PointType x1, PointType x2, PointType* W);
    void WcalculatePolar(int index, PointType x1, PointType x2, PointType* W);

    void WcalculatePolar(unsigned int* basePoint, PointType* x1, PointType* x2, PointType (*W)[9],
                         int i1, int i2, PointType* p1, PointType* p2, PointType* tmp1,
                         PointType* tmp2, PointType* gamma, double& size);
    void WcalculatePolar(unsigned int* basePoint, PointType* x1, PointType* x2, PointType (*W)[9],
                         int i1, int i2);

    void Wcalculate(unsigned int* basePoint, PointType* x1, PointType* x2, PointType* x3,
                    PointType (*W)[9], int i1, int i2, PointType* p1, PointType* p2, PointType* p3,
                    PointType* tmp1, PointType* tmp2, PointType* gamma, double& size);

    void Wcalculate(unsigned int* basePoint, PointType* x1, PointType* x2, PointType* x3,
                    PointType (*W)[9], int i1, int i2);
    void Wcalculate(unsigned int* basePoint, PointType* x1, PointType* x2, PointType (*W)[9],
                    int i1, int i2);

    void Wcalculate(int basePoint, PointType x1, PointType x2, PointType x3, PointType* W);

    void ValueInterpolate(int index, const PointType* W, const std::vector<PointType>& ValueArray,
                          PointType& result);
    void ValueInterpolate3d(int index, const PointType* W, const std::vector<PointType>& ValueArray,
                            PointType& result);

    void ChargeCalculate(int index, const PointType* W, PointType particleCharge, PointType* rho);
    void ChargeCalculate3d(int index, const PointType* W, PointType particleCharge, PointType* rho);

    void ChargeCalculate(
        unsigned int* isPeriodical, int i1, int i2,
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        int emType, arma::mat& Icoef, PointType* kI, unsigned char* flag,
        const std::vector<unsigned int>& emissionCells, unsigned int* startCell, PointType* r1,
        PointType* z1, unsigned int* index1, PointType (*W1)[9], PointType* r2, PointType* z2,
        unsigned int* index2, PointType (*W2)[9], PointType* I, PointType dt, PointType* rho,
        PointType* cartesianX1, PointType* cartesianX2, int flagType, PointType* cartesianX12,
        PointType* cartesianX22);
    void ChargeCalculate(
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        int emType, arma::mat& Icoef, PointType kI, unsigned char flag,
        const std::vector<unsigned int>& emissionCells, unsigned int startCell, PointType r1,
        PointType z1, int index1, const PointType* W1, PointType r2, PointType z2, int index2,
        const PointType* W2, PointType I, PointType dt, PointType* rho, PointType cartesianX1,
        PointType cartesianX2, int flagType, PointType cartesianX12, PointType cartesianX22);

    PointType CellVolume(int index);
    PointType GetH1(int index);
    PointType GetH2(int index);
};
#endif