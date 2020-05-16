#ifndef ParticleGridInterface_H
#define ParticleGridInterface_H
#include "armadillo"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

template <class PointType>
class ElectrodeCurrent;
class BoundaryConditions;
template <class PointType>
class BoundaryContainer2d;
template <class PointType>
class BoundaryContainer3d;
template <class PointType>
class IParticleShape2d;
template <class PointType>
class NearCathodeVolume;

namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
};
namespace Dmath
{
class imat;
};
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
class MeshContainer2d;
template <class PointType>
class Particles3dcil;
template <class PointType>
class Particles2dpolar;
template <class PointType>
class Particles2d;
template <class PointType>
class Particles3d;
template <class PointType>
class ParticleShape2dTSC;
template <class PointType>
class ParticleShape2dCIC;

template <class PointType>
class particlesFields;

template <class PointType>
class ParticleGridInterface
{
    friend class boost::serialization::access;
    friend MeshContainer2d<PointType>;
    friend Particles2d<PointType>;
    friend Particles2dpolar<PointType>;
    friend Particles3dcil<PointType>;

  private:
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& priorityParticleShapeType;
        ar& epsilon;
        ar& fabsorpFrac;
        ar& removeParticlesFlag;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& priorityParticleShapeType;
        ar& epsilon;
        ar& fabsorpFrac;
        ar& removeParticlesFlag;
    }
    double fabsorpFrac;
    double epsilon;

    std::shared_ptr<ParticleShape2dTSC<PointType>> TSCcellArray;
    std::shared_ptr<ParticleShape2dCIC<PointType>> CICcellArray;

    int                                    removeParticlesFlag;
    int                                    priorityParticleShapeType;
    std::vector<int>                       CICArray;
    std::vector<bool>                      isBoundary;
    std::vector<char>                      cellTypes;
    std::vector<std::vector<unsigned int>> searchIndexes;
    std::vector<std::vector<unsigned int>> searchIndexesLarge;
    std::vector<std::vector<unsigned int>> searchIndexesLarge1;
    std::vector<std::vector<unsigned int>> searchIndexesSpecial;
    std::vector<std::vector<unsigned int>> startCellNumbers;

    std::vector<std::vector<PointType>> tmpData1;
    std::vector<std::vector<PointType>> tmpData2;
    std::vector<std::vector<PointType>> gamma;

    const std::shared_ptr<BoundaryConditions> boundaryConditions;

    void AddCell(unsigned short particleShapeType, unsigned int basePoint, unsigned int levelLow,
                 unsigned int levelHigh);
    void AddCell(unsigned int basePoint, unsigned int p1, unsigned int p2, unsigned int p3,
                 unsigned int p4);
    void AddCell(unsigned short particleShapeType, unsigned int basePoint, unsigned int p1,
                 unsigned int p2, unsigned int p3, unsigned int p4, unsigned int p5,
                 unsigned int p6);
    void AddCell(unsigned int basePoint, const std::vector<int>& v);
    int InCell(PointType x1, PointType x2, const std::vector<unsigned int>& searchIndexes);
    int InCell(PointType x1, PointType x2, PointType x3,
               const std::vector<unsigned int>& searchIndexes);

    int InCellWithEps(PointType x1, PointType x2, const std::vector<unsigned int>& searchIndexes);
    int InCellWithEps(PointType x1, PointType x2, PointType x3,
                      const std::vector<unsigned int>& searchIndexes);

  public:
    void
    SearchBoundariesCells(std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries);
    void SearchBoundariesCells(
        std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries){};

    long long GetMemorySize();
    void init(int nflow, const std::shared_ptr<GridData3d<PointType>>& gridData,
              Dmath::imat templNumbIn, Dmath::imat flagMatrix,
              const std::vector<int>&                               boundarypoints,
              const std::shared_ptr<BoundaryContainer3d<PointType>> domainBoundary, int size,
              int numThreads, int blockSize);
    std::vector<unsigned int>& CheckParticlesBoundaries(
        const std::shared_ptr<BoundaryConditions>&                          boundaryConditions,
        std::vector<unsigned int>&                                          result,
        const std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>&    conductorList,
        const std::shared_ptr<Particles3d<PointType>>&                      state1,
        const std::shared_ptr<Particles3d<PointType>>& state2, PointType dt, PointType charge,
        PointType mass, int i1, int i2, int thread);
    void
    InitEmCells(const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
                int emType, const std::shared_ptr<Particles3d<PointType>>& state1, double dH,
                int flag, int flagInit);
    void Particles2GridPTI(
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        int emType, arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<Particles3d<PointType>>& state1, PointType (*W1)[9],
        const std::shared_ptr<Particles3d<PointType>>& state2, PointType (*W2)[9], PointType* rho,
        PointType dt, int flagType, int i1, int i2, int thread);
    std::vector<double> GetParameters();
    void SetParameters(std::vector<double> in);
    ParticleGridInterface();
    void SearchStartCellsEmission(
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<Particles3d<PointType>>>&         particles);

    template <class gridDatatype>
    void init(int nflow, const std::shared_ptr<gridDatatype>& gridData, Dmath::imat templNumbIn,
              Dmath::imat flagMatrix, const std::vector<int>& boundaryPoints,
              const std::shared_ptr<BoundaryContainer2d<PointType>>& domainBoundary, int size,
              int numThreads, int blockSize);

    void init(int nflow, const std::shared_ptr<GridData3d<PointType>>& gridData,
              Dmath::imat templNumbIn, Dmath::imat flagMatrix,
              const std::vector<int>&                                boundarypoints,
              const std::shared_ptr<BoundaryContainer2d<PointType>>& domainBoundary, int size,
              int numThreads, int blockSize);

    PointType GetH1(int index);
    PointType GetH2(int index);

    std::shared_ptr<Dmath::imat>              templNumb;
    std::vector<IParticleShape2d<PointType>*> cellArray;

    template <class particlesType>
    std::vector<unsigned int>& CheckParticlesBoundaries(
        const std::shared_ptr<BoundaryConditions>&                          boundaryConditions,
        std::vector<unsigned int>&                                          result,
        const std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
        const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>&    conductorList,
        const std::shared_ptr<particlesType>& state1, const std::shared_ptr<particlesType>& state2,
        PointType dt, PointType charge, PointType mass, int i1, int i2, int thread);

    template <class particlesType>
    void ApplyBoundaryCondition(
        const std::shared_ptr<BoundaryConditions>& boundaryConditions, int flag,
        std::vector<unsigned int>& empty, const std::string& conditionType,
        const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& conductorList,
        std::vector<double> conditionProperties, particlesType state1, particlesType state2,
        const std::vector<unsigned int>&           particlesNumbers,
        const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
        const std::vector<DGeo::Point<PointType>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, PointType dt, PointType charge,
        PointType mass, int i1, int thread);

    template <class particlesType>
    void ApplyDefaultCondition(
        const std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& conductorList,
        particlesType state1, particlesType state2,
        const std::vector<unsigned int>&           particlesNumbers,
        const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
        const std::vector<DGeo::Point<PointType>>& intersectionPoints,
        const std::vector<int>& intersectionBoundaryNumber, PointType dt, PointType charge,
        PointType mass, int thread);

    void Grid2Particles(const std::shared_ptr<Particles2d<PointType>>&      particles,
                        particlesFields<PointType>&                         particlesTmp,
                        const std::shared_ptr<const GridData2d<PointType>>& gridData,
                        PointType (*W)[9], int i1, int i2, int step, int recalculate);
    void Grid2Particles(const std::shared_ptr<Particles3dcil<PointType>>& particles,
                        particlesFields<PointType>&                       particlesTmp,
                        const std::shared_ptr<GridData2daxs<PointType>>&  gridData,
                        PointType (*W)[9], int i1, int i2, int step, int recalculate);
    void Grid2Particles(const std::shared_ptr<Particles2dpolar<PointType>>& particles,
                        particlesFields<PointType>&                         particlesTmp,
                        const std::shared_ptr<GridData2dpolar<PointType>>&  gridData,
                        PointType (*W)[9], int i1, int i2, int step, int recalculate);
    void Grid2Particles(const std::shared_ptr<Particles3d<PointType>>& particles,
                        particlesFields<PointType>&                    particlesTmp,
                        const std::shared_ptr<GridData3d<PointType>>&  gridData, PointType (*W)[9],
                        int i1, int i2, int step, int recalculate);

    void Particles2Grid(const std::shared_ptr<Particles2d<PointType>>& particles, PointType* rho,
                        PointType (*W)[9], int i1, int i2);
    void Particles2Grid(const std::shared_ptr<Particles3dcil<PointType>>& particles, PointType* rho,
                        PointType (*W)[9], int i1, int i2);
    void Particles2Grid(const std::shared_ptr<Particles2dpolar<PointType>>& particles,
                        PointType* rho, PointType (*W)[9], int i1, int i2);
    void Particles2Grid(const std::shared_ptr<Particles3d<PointType>>& particles, PointType* rho,
                        PointType (*W)[9], int i1, int i2);

    template <class particlesType>
    void Particles2GridPTI(
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        int emType, arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<particlesType>& state1, PointType (*W1)[9],
        const std::shared_ptr<particlesType>& state2, PointType (*W2)[9], PointType* rho,
        PointType dt, int flagType, int i1, int i2, int thread);

    void Particles2GridPTI(
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        int emType, arma::mat& Icoef, const std::vector<unsigned int>& emissionCells,
        const std::shared_ptr<Particles2dpolar<PointType>>& state1, PointType (*W1)[9],
        const std::shared_ptr<Particles2dpolar<PointType>>& state2, PointType (*W2)[9],
        PointType* rho, PointType dt, int flagType, int i1, int i2, int thread);

    template <class particlesType>
    void
    InitEmCells(const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
                int emType, const std::shared_ptr<particlesType>& state1, double dH, int flag,
                int flagInit);

    void Wcalculate(const std::shared_ptr<Particles3dcil<PointType>>& particles, PointType (*W)[9],
                    int i1, int i2, int thread, int flagStepEstimate);
    void Wcalculate(const std::shared_ptr<Particles2d<PointType>>& particles, PointType (*W)[9],
                    int i1, int i2, int thread, int flagStepEstimate);
    void Wcalculate(const std::shared_ptr<Particles2dpolar<PointType>>& particles,
                    PointType (*W)[9], int i1, int i2, int thread, int flagStepEstimate);
    void Wcalculate(const std::shared_ptr<Particles3d<PointType>>& particles, PointType (*W)[9],
                    int i1, int i2, int thread, int flagStepEstimate);

    template <class particlesType>
    std::vector<unsigned int> InCell(const std::shared_ptr<particlesType>& particles,
                                     std::vector<unsigned int> EmptyPlaces, int i1, int i2,
                                     int thread, int flow);

    std::vector<unsigned int> InCell(const std::shared_ptr<Particles3d<PointType>>& particles,
                                     std::vector<unsigned int> EmptyPlaces, int i1, int i2,
                                     int thread, int flow);

    template <class PointType1>
    void axsPolar(Particles3dcil<PointType1>*                     particles,
                  const std::shared_ptr<GridData2dpolar<double>>& gridData,
                  std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread);

    template <class PointType1>
    void axsPolar(Particles3d<PointType1>*                        particles,
                  const std::shared_ptr<GridData2dpolar<double>>& gridData,
                  std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread){

    };

    template <class PointType1>
    void axsPolar(Particles2dpolar<PointType1>*                   particles,
                  const std::shared_ptr<GridData2dpolar<double>>& gridData,
                  std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread){

    };

    template <class PointType1>
    void axsPolar(Particles2d<PointType1>*                        particles,
                  const std::shared_ptr<GridData2dpolar<double>>& gridData,
                  std::vector<unsigned int> EmptyPlaces, int i1, int i2, int thread){

    };

    std::vector<unsigned int>
    InCellWithEps(const std::shared_ptr<Particles3d<PointType>>& particles, int i1, int i2,
                  int thread, int flow);

    template <class particlesType>
    void SearchStartCells(int flow, const std::shared_ptr<particlesType>& particles);

    void SearchStartCells(int flow, const std::shared_ptr<Particles3d<PointType>>& particles);

    template <class particlesType>
    void SearchStartCellsEmission(
        const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
        std::vector<std::shared_ptr<particlesType>>&                  particles);

    void Charge2Density(PointType* rho, std::vector<int>& nonZeros);

    int InCellWithEps(PointType x1, PointType x2);

    int InCell(PointType x1, PointType x2);

    template <class particlesType>
    std::vector<unsigned int> InCellWithEps(const std::shared_ptr<particlesType>& particles, int i1,
                                            int i2, int thread, int flow);

    bool InCell0(PointType x1, PointType x2);
};

#endif