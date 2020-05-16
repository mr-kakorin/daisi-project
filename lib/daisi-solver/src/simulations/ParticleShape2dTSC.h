#ifndef PARTICLESHAPE2DTSC_H
#define PARTICLESHAPE2DTSC_H
#include <memory>
#include <vector>

template <class PointType> class BoundaryContainer2d;
template <class PointType> class BoundaryContainer3d;
template <class PointType> class GridData2daxs;
template <class PointType> class GridData2d;
template <class PointType> class GridData2dpolar;
template <class PointType> class GridData3d;
template <class PointType> class GridData3dpolar;
template <class PointType> class NearCathodeVolume;
template <class PointType> class ParticleShape2dTSC
{
  public:
    std::vector<unsigned int>     levelHigh;
    std::vector<unsigned int>     levelLow;
    std::vector<char>             cellType;
    std::vector<char>             isBoundary;
    std::vector<std::vector<int>> points;
    std::vector<PointType>        cellVolume;

    const PointType* x1Array;
    const PointType* x2Array;
    unsigned int     cellNumber;

    ParticleShape2dTSC(){

    };

    void init(int size, const PointType* x1ArrayIn, const PointType* x2ArrayIn);

    void AddCell(const std::shared_ptr<GridData2d<PointType>>& gridData, unsigned int levelLowIn, unsigned int levelHighIn,
                 char cellTypeIn, std::vector<int> pointsIn, const std::vector<int>& flagOut,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);

    void AddCell(const std::shared_ptr<GridData2daxs<PointType>>& gridData, unsigned int levelLowIn, unsigned int levelHighIn,
                 char cellTypeIn, std::vector<int> pointsIn, const std::vector<int>& flagOut,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);

    void AddCell(const std::shared_ptr<GridData2dpolar<PointType>>& gridData, unsigned int levelLowIn,
                 unsigned int levelHighIn, char cellTypeIn, std::vector<int> pointsIn, const std::vector<int>& flagOut,
                 const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
                 const std::vector<int>& boundarypoints, int problemType);

    bool InCell(int index, PointType x1, PointType x2);
    int InCell(PointType x1, PointType x2);

    int InCellWithEps(PointType x1, PointType x2)
    {
        return 0;
    };
    bool InCellWithEps(int index, PointType x1, PointType x2)
    {
        return true;
    };

    bool InCell0(int index, PointType x1, PointType x2);
    bool InCell1(int index, PointType x1, PointType x2);
    bool InCell2(int index, PointType x1, PointType x2);
    bool InCell3(int index, PointType x1, PointType x2);
    bool InCell4(int index, PointType x1, PointType x2);

    void WcalculatePolar(int index, PointType x1, PointType x2, PointType* W);

    void Wcalculate(int index, PointType x1, PointType x2, PointType* W);
    void Wcalculate0(int index, PointType x1, PointType x2, PointType* W);
    void Wcalculate1(int index, PointType x1, PointType x2, PointType* W);
    void Wcalculate2(int index, PointType x1, PointType x2, PointType* W);
    void Wcalculate3(int index, PointType x1, PointType x2, PointType* W);
    void Wcalculate4(int index, PointType x1, PointType x2, PointType* W);
    void Wcalculate0Vector(int* index, PointType* x1, PointType* x2, PointType (*W)[9], int size);

    void ValueInterpolate(int index, const PointType* W, const std::vector<PointType>& ValueArray, PointType& result);
    void ValueInterpolate0(int index, const PointType* W, const std::vector<PointType>& ValueArray, PointType& result);
    void ValueInterpolate1(int index, const PointType* W, const std::vector<PointType>& ValueArray, PointType& result);
    void ValueInterpolate2(int index, const PointType* W, const std::vector<PointType>& ValueArray, PointType& result);
    void ValueInterpolate3(int index, const PointType* W, const std::vector<PointType>& ValueArray, PointType& result);
    void ValueInterpolate4(int index, const PointType* W, const std::vector<PointType>& ValueArray, PointType& result);

    void ChargeCalculate(int index, const PointType* W, PointType particleCharge, PointType* rho);
    void ChargeCalculate0(int index, const PointType* W, PointType particleCharge, PointType* rho);
    void ChargeCalculate1(int index, const PointType* W, PointType particleCharge, PointType* rho);
    void ChargeCalculate2(int index, const PointType* W, PointType particleCharge, PointType* rho);
    void ChargeCalculate3(int index, const PointType* W, PointType particleCharge, PointType* rho);
    void ChargeCalculate4(int index, const PointType* W, PointType particleCharge, PointType* rho);

    void ChargeCalculate(PointType r1, PointType z1, int index1, const PointType* W1, PointType r2, PointType z2,
                         int index2, const PointType* W2, PointType dr, PointType dz, PointType I, PointType dt,
                         PointType* rho);

    void      ChargeZeroing();
    PointType CellArea();
    PointType CellVolume(int index);

    PointType GetH1(int index);
    PointType GetH2(int index);
};
#endif