#ifndef PARTICLEGRIDINTERFACE2D_H
#define PARTICLEGRIDINTERFACE2D_H
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "FlagStringsSolver.h"
#include "GridData.h"
#include "Particle.h"
#include "ParticleShape2d.h"
#include "ParticleShape2dCIC.h"
#include "ParticleShape2dTSC.h"
#include "armadillo"
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
template <class PointType> class MeshContainer2d;
template <class PointType> class Particles3dcil;
template <class PointType> class Particles2dpolar;
template <class PointType> class Particles2d;
template <class PointType> class ParticleGridInterface2d
{
    friend class boost::serialization::access;
    friend MeshContainer2d<PointType>;
    friend Particles2d<PointType>;
    friend Particles2dpolar<PointType>;
    friend Particles3dcil<PointType>;

  private:
    int iter;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const
    {
        ar& priorityParticleShapeType;
    }
    template <class Archive> void load(Archive& ar, const unsigned int)
    {
        ar& priorityParticleShapeType;
    }

    ParticleShape2dTSC<PointType> TSCcellArray;
    ParticleShape2dCIC<PointType> CICcellArray;

    int               priorityParticleShapeType;
    std::vector<int>  CICArray;
    std::vector<bool> isBoundary;
    int               s;
    std::vector<int>  list1;
    std::vector<char> cellTypes;
    std::vector<int>  searchIndexes;

    std::vector<std::vector<int>> searchIndexesAll;

    BoundaryConditions*           boundaryConditions;
    std::vector<std::vector<int>> list;
    const PointType*              x1Array;
    const PointType*              x2Array;

    std::vector<int> indexes;
    std::vector<int> cells;

    void AddCell(unsigned short particleShapeType, unsigned int basePoint, unsigned int levelLow,
                 unsigned int levelHigh);
    void AddCell(unsigned int basePoint, unsigned int p1, unsigned int p2, unsigned int p3, unsigned int p4);
    void AddCell(unsigned short particleShapeType, unsigned int basePoint, unsigned int p1, unsigned int p2,
                 unsigned int p3, unsigned int p4, unsigned int p5, unsigned int p6);
    void AddCell(unsigned int basePoint, const std::vector<int>& v);
    int InCell(PointType x1, PointType x2, const std::vector<int>& searchIndexes);
    int InCellWithEps(PointType x1, PointType x2, const std::vector<int>& searchIndexes);

  public:
    template <class gridDatatype>
    void init(gridDatatype gridData, Dmath::imat templNumbIn, Dmath::imat flagMatrix,
              const std::vector<int>& boundaryPoints, const BoundaryContainer2d<PointType>& domainBoundary, int size);

    PointType GetH1(int index);
    PointType GetH2(int index);

    Dmath::imat                               templNumb;

    template <class particlesType>
    std::vector<unsigned int> CheckParticlesBoundaries(const std::vector<BoundaryContainer2d<PointType>>& boundaries,
                                                       std::vector<ElectrodeCurrent<PointType>>&          conductorList,
                                                       particlesType state1, particlesType state2, PointType dt,
                                                       PointType charge, PointType mass);

    template <class particlesType>
    void ApplyBoundaryCondition(std::vector<unsigned int>& empty, const std::string& conditionType,
                                std::vector<ElectrodeCurrent<PointType>>& conductorList,
                                std::vector<double> conditionProperties, particlesType state1, particlesType state2,
                                const std::vector<unsigned int>&           particlesNumbers,
                                const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
                                const std::vector<DGeo::Point<PointType>>& intersectionPoints,
                                const std::vector<int>& intersectionBoundaryNumber, PointType dt, PointType charge,
                                PointType mass);

    template <class particlesType>
    void ApplyDefaultCondition(std::vector<ElectrodeCurrent<PointType>>& conductorList, particlesType state1,
                               particlesType state2, const std::vector<unsigned int>& particlesNumbers,
                               const std::vector<DGeo::Edge<PointType>>&  intersectionEdges,
                               const std::vector<DGeo::Point<PointType>>& intersectionPoints,
                               const std::vector<int>& intersectionBoundaryNumber, PointType dt, PointType charge,
                               PointType mass);

    void Grid2Particles(Particles2d<PointType>* particles, const GridData2d<PointType>* const gridData,
                        PointType (*W)[9]);
    void Grid2Particles(Particles3dcil<PointType>* particles, const GridData2daxs<PointType>* const gridData,
                        PointType (*W)[9]);
    void Grid2Particles(Particles2dpolar<PointType>* particles, const GridData2dpolar<PointType>* const gridData,
                        PointType (*W)[9]);
    void Particles2Grid(const Particles2d<PointType>* const particles, GridData2d<PointType>* gridData,
                        PointType (*W)[9]);
    void Particles2Grid(const Particles3dcil<PointType>* const particles, GridData2daxs<PointType>* gridData,
                        PointType (*W)[9]);
    void Particles2Grid(const Particles2dpolar<PointType>* const particles, GridData2dpolar<PointType>* gridData,
                        PointType (*W)[9]);

    template <class particlesType>
    void Particles2GridPTI(const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
                           arma::mat& Icoef, const std::vector<unsigned int>& emissionCells, particlesType state1,
                           PointType (*W1)[9], particlesType state2, PointType (*W2)[9], PointType* rho, PointType dt);

    template <class particlesType>
    void InitEmCells(const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
                     particlesType state1, double dH);

    void Wcalculate(Particles3dcil<PointType>* particles, PointType (*W)[9]);
    void Wcalculate(Particles2d<PointType>* particles, PointType (*W)[9]);
    void Wcalculate(Particles2dpolar<PointType>* particles, PointType (*W)[9]);

    template <class particlesType> std::vector<unsigned int> InCell(particlesType particles);

    template <class particlesType> void SearchStartCells(particlesType particles);

    template <class particlesType>
    void SearchStartCellsEmission(const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes,
                                  particlesType                                                 particles);

    void Charge2Density(PointType* rho);

    int InCellWithEps(PointType x1, PointType x2);

    int InCell(PointType x1, PointType x2);

    template <class particlesType> std::vector<unsigned int> InCellWithEps(particlesType particles);

    bool InCell0(PointType x1, PointType x2);

    int GetPriorityParticleShapeType()
    {
        return priorityParticleShapeType;
    };
    void SetPriorityParticleShapeType(int in)
    {
        priorityParticleShapeType = in;
    };
    ParticleGridInterface2d()
    {
        iter                      = 0;
        priorityParticleShapeType = 0;
    };
    void initBoundaries(BoundaryConditions* boundaryConditions)
    {
        list1 = boundaryConditions->GetDefaultConditionsList();
        s     = boundaryConditions->PropertyConditionListSize();

        list.clear();
        for (int i = 0; i < s; i++)
        {
            list.push_back(boundaryConditions->GetPropertyConditionsBoundariesList(i));
        }
        this->boundaryConditions = boundaryConditions;
    };
};

#endif