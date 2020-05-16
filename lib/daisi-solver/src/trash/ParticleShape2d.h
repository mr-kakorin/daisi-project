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
class IParticleShape2d
{
    friend class boost::serialization::access;

  public:
    unsigned int basePoint;
    unsigned int levelLow;
    unsigned int levelHigh;
    IParticleShape2d(){};
    IParticleShape2d(unsigned int basePointIn, unsigned int levelLowIn, unsigned int levelHighIn)
    {
        basePoint = basePointIn;
        levelLow  = levelLowIn;
        levelHigh = levelHighIn;
    };
    virtual bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                        const std::vector<PointType>& x2Array) = 0;
    virtual void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                            const std::vector<PointType>& x2Array, std::vector<PointType>& W) = 0;
    virtual void ValueInterpolate(const std::vector<PointType>& W,
                                  const std::vector<PointType>& ValueArray, PointType& result) = 0;
    virtual void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                                 std::vector<PointType>& rho)          = 0;
    virtual void      ChargeZeroing()                                  = 0;
    virtual PointType CellArea()                                       = 0;
    virtual PointType CellVolume(const std::vector<PointType>& rArray) = 0;

  public:
    IParticleShape2d& operator=(const IParticleShape2d& right)
    {
        if (this == &right)
        {
            return *this;
        }
        basePoint = right.basePoint;
        levelLow  = right.levelLow;
        levelHigh = right.levelHigh;
        return *this;
    }
};
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
template <class PointType>
class ParticleShapeCIC2d : public IParticleShape2d<PointType>
{
    friend class boost::serialization::access;

  public:
    PointType CellArea()
    {
        return gridSteps<PointType>::H1 * gridSteps<PointType>::H2;
    };
    bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                const std::vector<PointType>& x2Array)
    {
        if (x1 >= x1Array[basePoint] && x1 <= x1Array[basePoint + 1] && x2 >= x2Array[basePoint] &&
            x2 <= x2Array[levelHigh])
            return true;
        return false;
    };
    void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                    const std::vector<PointType>& x2Array, std::vector<PointType>& W)
    {
        PointType wx11 = (x1Array[basePoint + 1] - x1) / gridSteps<PointType>::H1;
        PointType wx12 = 1 - wx11;
        PointType wx21 = (x2Array[levelHigh] - x2) / gridSteps<PointType>::H2;
        PointType wx22 = 1 - wx21;
        W[0]           = wx11 * wx21;
        W[1]           = wx12 * wx21;
        W[2]           = wx11 * wx22;
        W[3]           = wx12 * wx22;
    };
    void ValueInterpolate(const std::vector<PointType>& W, const std::vector<PointType>& ValueArray,
                          PointType& result)
    {
        result = W[0] * ValueArray[basePoint] + W[1] * ValueArray[basePoint + 1] +
                 W[2] * ValueArray[levelHigh] + W[3] * ValueArray[levelHigh + 1];
    };
    void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                         std::vector<PointType>& rho)
    {
        rho[basePoint]     = rho[basePoint] + W[0] * particleCharge;
        rho[basePoint + 1] = rho[basePoint + 1] + W[1] * particleCharge;
        rho[levelHigh]     = rho[levelHigh] + W[2] * particleCharge;
        rho[levelHigh + 1] = rho[levelHigh + 1] + W[3] * particleCharge;
    };
    void      ChargeZeroing(){};
    PointType CellVolume(const std::vector<PointType>& rArray)
    {
        return PointType(commtools::PI() * gridSteps<PointType>::H2 *
                         ((rArray[basePoint] + gridSteps<PointType>::H1 / 2) *
                              (rArray[basePoint] + gridSteps<PointType>::H1 / 2) -
                          (rArray[basePoint] - gridSteps<PointType>::H1 / 2) *
                              (rArray[basePoint] - gridSteps<PointType>::H1 / 2)));
    };
    ParticleShapeCIC2d(){};
    ParticleShapeCIC2d(unsigned int basePointIn, unsigned int levelLowIn, unsigned int levelHighIn);
};
template <class PointType>
ParticleShapeCIC2d<PointType>::ParticleShapeCIC2d(unsigned int basePointIn, unsigned int levelLowIn,
                                                  unsigned int levelHighIn)
    : IParticleShape2d<PointType>(basePointIn, levelLowIn, levelHighIn){

      };
template <class PointType>
class ParticleShapeBoundaryCIC2d : public IParticleShape2d<PointType>
{
    unsigned int p0;
    unsigned int p1;
    unsigned int p2;
    unsigned int p3;
    friend class boost::serialization::access;

  public:
    PointType CellArea()
    {
        return gridSteps<PointType>::H1 * gridSteps<PointType>::H2 / 4;
    };
    PointType CellVolume(const std::vector<PointType>& rArray)
    {
        if (basePoint == p0 || basePoint == p2)
            return PointType(0.5 * commtools::PI() * gridSteps<PointType>::H2 *
                             ((rArray[basePoint] + gridSteps<PointType>::H1 / 2) *
                                  (rArray[basePoint] + gridSteps<PointType>::H1 / 2) -
                              (rArray[basePoint]) * (rArray[basePoint])));
        if (basePoint == p1 || basePoint == p3)
            return PointType(0.5 * commtools::PI() * gridSteps<PointType>::H2 *
                             ((rArray[basePoint]) * (rArray[basePoint]) -
                              (rArray[basePoint] - gridSteps<PointType>::H1 / 2) *
                                  (rArray[basePoint] - gridSteps<PointType>::H1 / 2)));
        return 0;
    };
    bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                const std::vector<PointType>& x2Array)
    {
        PointType tmp1 = x1Array[p0];
        PointType tmp2 = x1Array[p1];
        PointType tmp3 = x1Array[p2];
        PointType tmp4 = x1Array[basePoint];
        if (x1 >= x1Array[p0] && x1 <= x1Array[p1] && x2 >= x2Array[p0] && x2 <= x2Array[p2] &&
            std::abs(x1 - x1Array[basePoint]) < gridSteps<PointType>::H1 / 2 &&
            std::abs(x2 - x2Array[basePoint]) < gridSteps<PointType>::H2 / 2)
            return true;
        return false;
    };
    void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                    const std::vector<PointType>& x2Array, std::vector<PointType>& W)
    {

        PointType XX = (x2Array[basePoint] - x2) / gridSteps<PointType>::H2;

        PointType wx21;
        PointType wx22;
        if (basePoint == p1 || basePoint == p0)
        {
            wx21 = PointType(0.75 - XX * XX);
            wx22 =
                PointType(0.5 * (0.5 + XX) * (0.5 + XX)) + PointType(0.5 * (0.5 - XX) * (0.5 - XX));
        };

        if (basePoint == p2 || basePoint == p3)
        {
            wx22 = PointType(0.75 - XX * XX);
            wx21 =
                PointType(0.5 * (0.5 + XX) * (0.5 + XX)) + PointType(0.5 * (0.5 - XX) * (0.5 - XX));
        };

        PointType x = std::abs(x1 - x1Array[basePoint]) / gridSteps<PointType>::H1;

        // PointType wx12 = (x1Array[p1] - x1) / gridSteps<PointType>::H1;

        PointType wx11 =
            8 * x * x * x - 8.000047376600016 * x * x + 3.000011137807037 * x + 0.000000313816663;

        /*if (std::abs(XX1 - 1) > 0.00000001)
        {
                wx11 = wx11 + (0.5*(0.5 - (1 - XX1))*(0.5 - (1 - XX1)));
        };*/

        PointType wx12 = 1 - wx11;

        //	PointType wx11 = (x1Array[p1] - x1) / gridSteps<PointType>::H1;
        //	PointType wx12 = 1 - wx11;
        // PointType wx21 = (x2Array[p2] - x2) / gridSteps<PointType>::H2;
        //	PointType wx22 = 1 - wx21;
        W[0] = wx11 * wx21;
        W[1] = wx12 * wx21;
        W[2] = wx11 * wx22;
        W[3] = wx12 * wx22;
    };
    void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                         std::vector<PointType>& rho)
    {
        rho[p0] = rho[p0] + W[0] * particleCharge;
        rho[p1] = rho[p1] + W[1] * particleCharge;
        rho[p2] = rho[p2] + W[2] * particleCharge;
        rho[p3] = rho[p3] + W[3] * particleCharge;
    };
    void ChargeZeroing(){};
    void ValueInterpolate(const std::vector<PointType>& W, const std::vector<PointType>& ValueArray,
                          PointType& result)
    {
        result = W[0] * ValueArray[p0] + W[1] * ValueArray[p1] + W[2] * ValueArray[p2] +
                 W[3] * ValueArray[p3];
    };
    ParticleShapeBoundaryCIC2d(){};
    ParticleShapeBoundaryCIC2d(unsigned int base, unsigned int p0In, unsigned int p1In,
                               unsigned int p2In, unsigned int p3In)
    {
        basePoint = base;
        p0        = p0In;
        p1        = p1In;
        p2        = p2In;
        p3        = p3In;
    };
};
template <class PointType>
class ParticleShapeBoundary3CellCIC2d : public IParticleShape2d<PointType>
{
    std::vector<ParticleShapeBoundaryCIC2d<PointType>> cells;
    friend class boost::serialization::access;

  public:
    PointType CellArea()
    {
        return 3 * gridSteps<PointType>::H1 * gridSteps<PointType>::H2 / 4;
    };
    bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                const std::vector<PointType>& x2Array)
    {
        for (int i = 0; i < cells.size(); i++)
        {
            if (cells[i].InCell(x1, x2, x1Array, x2Array))
                return true;
        }
        return false;
    };
    void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                    const std::vector<PointType>& x2Array, std::vector<PointType>& W)
    {
        for (int i = 0; i < cells.size(); i++)
        {
            if (cells[i].InCell(x1, x2, x1Array, x2Array))
            {
                cells[i].Wcalculate(x1, x2, x1Array, x2Array, W);
                W[4] = float(i);
                return;
            }
        }
    };
    void ValueInterpolate(const std::vector<PointType>& W, const std::vector<PointType>& ValueArray,
                          PointType& result)
    {
        cells[int(W[4])].ValueInterpolate(W, ValueArray, result);
    };
    void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                         std::vector<PointType>& rho)
    {
        cells[int(W[4])].ChargeCalculate(W, particleCharge, rho);
    };
    PointType CellVolume(const std::vector<PointType>& rArray)
    {
        return cells[0].CellVolume(rArray) + cells[1].CellVolume(rArray) +
               cells[2].CellVolume(rArray);
    };
    void ChargeZeroing(){};
    ParticleShapeBoundary3CellCIC2d(){};
    void AddCell(ParticleShapeBoundaryCIC2d<PointType> in)
    {
        cells.push_back(in);
    }
};
template <class PointType>
class ParticleShapeBoundaryCICx1TSCx2_2d : public IParticleShape2d<PointType>
{
    unsigned int p0;
    unsigned int p1;
    unsigned int p2;
    unsigned int p3;
    unsigned int p4;
    unsigned int p5;
    friend class boost::serialization::access;

  public:
    PointType CellArea()
    {
        return 2 * gridSteps<PointType>::H1 * gridSteps<PointType>::H2 / 4;
    };
    bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                const std::vector<PointType>& x2Array)
    {
        if (std::abs(x1 - x1Array[basePoint]) < gridSteps<PointType>::H1 / 2 &&
            std::abs(x2 - x2Array[basePoint]) < gridSteps<PointType>::H2 / 2 && x1 <= x1Array[p1] &&
            x1 >= x1Array[p0])
            return true;
        return false;
    };
    void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                    const std::vector<PointType>& x2Array, std::vector<PointType>& W)
    {

        PointType XX = (x1 - x1Array[p0]) / gridSteps<PointType>::H1;

        PointType wx11 = (x1Array[p1] - x1) / gridSteps<PointType>::H1;

        if (std::abs(XX - 1) > 0.00000001)
        {
            wx11 = wx11 + (0.5 * (0.5 - (1 - XX)) * (0.5 - (1 - XX)));
        };

        PointType wx12 = 1 - wx11;

        PointType X = (x2Array[p2] - x2) / gridSteps<PointType>::H2;

        PointType wx2_1 = PointType(0.5 * (0.5 + X) * (0.5 + X));
        PointType wx21  = PointType(0.5 * (0.5 - X) * (0.5 - X));
        PointType wx20  = PointType(1 - wx2_1 - wx21);
        PointType wx20T = PointType(0.75 - X * X);

        W[0] = wx11 * wx2_1;
        W[1] = wx12 * wx2_1;
        W[2] = wx11 * wx20;
        W[3] = wx12 * wx20;
        W[4] = wx11 * wx21;
        W[5] = wx12 * wx21;
    };
    void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                         std::vector<PointType>& rho)
    {
        rho[p0] = rho[p0] + W[0] * particleCharge;
        rho[p1] = rho[p1] + W[1] * particleCharge;
        rho[p2] = rho[p2] + W[2] * particleCharge;
        rho[p3] = rho[p3] + W[3] * particleCharge;
        rho[p4] = rho[p4] + W[4] * particleCharge;
        rho[p5] = rho[p5] + W[5] * particleCharge;
    };
    PointType CellVolume(const std::vector<PointType>& rArray)
    {

        if (basePoint == p2)
            return PointType(commtools::PI() * gridSteps<PointType>::H2 *
                             ((rArray[basePoint] + gridSteps<PointType>::H1 / 2) *
                                  (rArray[basePoint] + gridSteps<PointType>::H1 / 2) -
                              (rArray[basePoint]) * (rArray[basePoint])));
        if (basePoint == p3)
            return PointType(commtools::PI() * gridSteps<PointType>::H2 *
                             ((rArray[basePoint]) * (rArray[basePoint]) -
                              (rArray[basePoint] - gridSteps<PointType>::H1 / 2) *
                                  (rArray[basePoint] - gridSteps<PointType>::H1 / 2)));
        return 0;
    };

    void ChargeZeroing(){};

    void ValueInterpolate(const std::vector<PointType>& W, const std::vector<PointType>& ValueArray,
                          PointType& result)
    {
        result = W[0] * ValueArray[p0] + W[1] * ValueArray[p1] + W[2] * ValueArray[p2] +
                 W[3] * ValueArray[p3] + W[4] * ValueArray[p4] + W[5] * ValueArray[p5];
    };
    ParticleShapeBoundaryCICx1TSCx2_2d(){};
    ParticleShapeBoundaryCICx1TSCx2_2d(unsigned int basePointIn, unsigned int p0In,
                                       unsigned int p1In, unsigned int p2In, unsigned int p3In,
                                       unsigned int p4In, unsigned int p5In)
    {
        basePoint = basePointIn;
        p0        = p0In;
        p1        = p1In;
        p2        = p2In;
        p3        = p3In;
        p4        = p4In;
        p5        = p5In;
    };
};
template <class PointType>
class ParticleShapeBoundaryCICx2TSCx1_2d : public IParticleShape2d<PointType>

{
    unsigned int p0;
    unsigned int p1;
    unsigned int p2;
    unsigned int p3;
    unsigned int p4;
    unsigned int p5;
    friend class boost::serialization::access;

  public:
    PointType CellArea()
    {
        return 3 * gridSteps<PointType>::H1 * gridSteps<PointType>::H2 / 4;
    };
    bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                const std::vector<PointType>& x2Array)
    {
        if (std::abs(x1 - x1Array[basePoint]) < gridSteps<PointType>::H1 / 2 &&
            std::abs(x2 - x2Array[basePoint]) < gridSteps<PointType>::H2 / 2 && x2 <= x2Array[p3] &&
            x2 >= x2Array[p0])
            return true;
        return false;
    };
    void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                    const std::vector<PointType>& x2Array, std::vector<PointType>& W)
    {
        PointType XX = (x2Array[basePoint] - x2) / gridSteps<PointType>::H2;

        PointType wx21;
        PointType wx22;
        if (basePoint == p1)
        {
            wx21 = PointType(0.75 - XX * XX);
            wx22 =
                PointType(0.5 * (0.5 + XX) * (0.5 + XX)) + PointType(0.5 * (0.5 - XX) * (0.5 - XX));
        };

        if (basePoint == p4)
        {
            wx22 = PointType(0.75 - XX * XX);
            wx21 =
                PointType(0.5 * (0.5 + XX) * (0.5 + XX)) + PointType(0.5 * (0.5 - XX) * (0.5 - XX));
        };
        //	PointType wx21 = (x2Array[p3] - x2) / gridSteps<PointType>::H2;
        //	PointType wx22 = 1 - wx21;

        PointType X = (x1Array[p1] - x1) / gridSteps<PointType>::H1;

        PointType wx1_1 = PointType(0.5 * (0.5 + X) * (0.5 + X));
        PointType wx11  = PointType(0.5 * (0.5 - X) * (0.5 - X));
        PointType wx10  = PointType(1 - wx1_1 - wx11);
        PointType wx10T = PointType(0.75 - X * X);

        W[0] = wx21 * wx1_1;
        W[1] = wx21 * wx10;
        W[2] = wx21 * wx11;
        W[3] = wx22 * wx1_1;
        W[4] = wx22 * wx10;
        W[5] = wx22 * wx11;
    };
    void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                         std::vector<PointType>& rho)
    {
        rho[p0] = rho[p0] + W[0] * particleCharge;
        rho[p1] = rho[p1] + W[1] * particleCharge;
        rho[p2] = rho[p2] + W[2] * particleCharge;
        rho[p3] = rho[p3] + W[3] * particleCharge;
        rho[p4] = rho[p4] + W[4] * particleCharge;
        rho[p5] = rho[p5] + W[5] * particleCharge;
    };
    void      ChargeZeroing(){};
    PointType CellVolume(const std::vector<PointType>& rArray)
    {
        return PointType(0.5 * commtools::PI() * gridSteps<PointType>::H2 *
                         ((rArray[basePoint] + gridSteps<PointType>::H1 / 2) *
                              (rArray[basePoint] + gridSteps<PointType>::H1 / 2) -
                          (rArray[basePoint] - gridSteps<PointType>::H1 / 2) *
                              (rArray[basePoint] - gridSteps<PointType>::H1 / 2)));
    };

    void ValueInterpolate(const std::vector<PointType>& W, const std::vector<PointType>& ValueArray,
                          PointType& result)
    {
        result = W[0] * ValueArray[p0] + W[1] * ValueArray[p1] + W[2] * ValueArray[p2] +
                 W[3] * ValueArray[p3] + W[4] * ValueArray[p4] + W[5] * ValueArray[p5];
    };
    ParticleShapeBoundaryCICx2TSCx1_2d(){};
    ParticleShapeBoundaryCICx2TSCx1_2d(unsigned int basePointIn, unsigned int p0In,
                                       unsigned int p1In, unsigned int p2In, unsigned int p3In,
                                       unsigned int p4In, unsigned int p5In)
    {
        basePoint = basePointIn;
        p0        = p0In;
        p1        = p1In;
        p2        = p2In;
        p3        = p3In;
        p4        = p4In;
        p5        = p5In;
    };
};

template <class PointType>
class ParticleShapeTSC2d : public IParticleShape2d<PointType>
{
    friend class boost::serialization::access;

  public:
    PointType CellArea()
    {
        return gridSteps<PointType>::H1 * gridSteps<PointType>::H2;
    };
    bool InCell(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                const std::vector<PointType>& x2Array)
    {
        if (std::abs(x1 - x1Array[basePoint]) < gridSteps<PointType>::H1 / 2 &&
            std::abs(x2 - x2Array[basePoint]) < gridSteps<PointType>::H2 / 2)
            return true;
        return false;
    };
    void Wcalculate(PointType x1, PointType x2, const std::vector<PointType>& x1Array,
                    const std::vector<PointType>& x2Array, std::vector<PointType>& W)
    {
        PointType X1 = (x1Array[basePoint] - x1) / gridSteps<PointType>::H1;

        PointType wx1_1 = PointType(0.5 * (0.5 + X1) * (0.5 + X1));
        PointType wx11  = PointType(0.5 * (0.5 - X1) * (0.5 - X1));
        PointType wx10  = PointType(1 - wx1_1 - wx11);
        PointType wx10T = PointType(0.75 - X1 * X1);

        PointType X2 = (x2Array[basePoint] - x2) / gridSteps<PointType>::H2;

        PointType wx2_1 = PointType(0.5 * (0.5 + X2) * (0.5 + X2));
        PointType wx21  = PointType(0.5 * (0.5 - X2) * (0.5 - X2));
        PointType wx20  = PointType(1 - wx2_1 - wx21);
        PointType wx20T = PointType(0.75 - X2 * X2);

        W[0] = wx2_1 * wx1_1;
        W[1] = wx2_1 * wx10;
        W[2] = wx2_1 * wx11;
        W[3] = wx20 * wx1_1;
        W[4] = wx20 * wx10;
        W[5] = wx20 * wx11;
        W[6] = wx21 * wx1_1;
        W[7] = wx21 * wx10;
        W[8] = wx21 * wx11;
    };
    void ChargeCalculate(const std::vector<PointType>& W, PointType particleCharge,
                         std::vector<PointType>& rho)
    {
        rho[levelLow - 1]  = rho[levelLow - 1] + W[0] * particleCharge;
        rho[levelLow]      = rho[levelLow] + W[1] * particleCharge;
        rho[levelLow + 1]  = rho[levelLow + 1] + W[2] * particleCharge;
        rho[basePoint - 1] = rho[basePoint - 1] + W[3] * particleCharge;
        rho[basePoint]     = rho[basePoint] + W[4] * particleCharge;
        rho[basePoint + 1] = rho[basePoint + 1] + W[5] * particleCharge;
        rho[levelHigh - 1] = rho[levelHigh - 1] + W[6] * particleCharge;
        rho[levelHigh]     = rho[levelHigh] + W[7] * particleCharge;
        rho[levelHigh + 1] = rho[levelHigh + 1] + W[8] * particleCharge;
    };
    void      ChargeZeroing(){};
    PointType CellVolume(const std::vector<PointType>& rArray)
    {
        return PointType(commtools::PI() * gridSteps<PointType>::H2 *
                         ((rArray[basePoint] + gridSteps<PointType>::H1 / 2) *
                              (rArray[basePoint] + gridSteps<PointType>::H1 / 2) -
                          (rArray[basePoint] - gridSteps<PointType>::H1 / 2) *
                              (rArray[basePoint] - gridSteps<PointType>::H1 / 2)));
    };
    void ValueInterpolate(const std::vector<PointType>& W, const std::vector<PointType>& ValueArray,
                          PointType& result)
    {
        result = W[0] * ValueArray[levelLow - 1] + W[1] * ValueArray[levelLow] +
                 W[2] * ValueArray[levelLow + 1] + W[3] * ValueArray[basePoint - 1] +
                 W[4] * ValueArray[basePoint] + W[5] * ValueArray[basePoint + 1] +
                 W[6] * ValueArray[levelHigh - 1] + W[7] * ValueArray[levelHigh] +
                 W[8] * ValueArray[levelHigh + 1];
    };
    ParticleShapeTSC2d(unsigned int basePointIn, unsigned int levelLowIn, unsigned int levelHighIn);
};
template <class PointType>
ParticleShapeTSC2d<PointType>::ParticleShapeTSC2d(unsigned int basePointIn, unsigned int levelLowIn,
                                                  unsigned int levelHighIn)
    : IParticleShape2d<PointType>(basePointIn, levelLowIn, levelHighIn){};
#endif