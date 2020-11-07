#ifndef PARTICLESOURCE_H
#define PARTICLESOURCE_H
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

template <class PointType>
class GridData2dpolar;
template <class PointType>
class BoundaryContainer2d;

namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
}
template <class PointType>
class ElectrodeCurrent;
template <class PointType>
class ParticleSourceEdge
{
  public:
    friend class boost::serialization::access;
    std::shared_ptr<DGeo::Edge<PointType>> extractingEdge;
    double                                 curveLength;
    double                                 currentDensity;
    double                                 maximalCurrentDensity;
    double                                 accumulatedCurrentDensity;
    double                                 normalX;
    double                                 normalY;
    double                                 alphaNormal;
    double                                 flagNormal;
    int                                    cellNumber;
    double                                 E;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    ParticleSourceEdge(DGeo::Edge<PointType>& Edge, double curveLengthPrev, int cell, double X,
                       double Y);
    ParticleSourceEdge();
};
template <class PointType>
class ParticleSource2d
{
    friend class boost::serialization::access;
    bool EdgeCmp(DGeo::Edge<PointType>& Edge1, DGeo::Edge<PointType>& Edge2);
    std::shared_ptr<ElectrodeCurrent<PointType>> electrode;

  public:
    double sourceCurrent;
    void SetFlowCurrent(double res);
    double                 ErAverage;
    int                    lastIndex;
    std::vector<PointType> polinom;

    std::vector<ParticleSourceEdge<PointType>> sourceSurface;

    void SetZCoordinate(double Z);
    void resetSearch();
    void setErAverage(double in);
    double getErAverage();

    std::vector<DGeo::Point<double>> GetSufacePoints();

    std::vector<PointType> GetParticle(PointType L1, PointType L2, int flag);

	bool GetParticleOptimized(PointType L1, PointType L2, int flag, PointType*const& out, double const timestep);

    double length();

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;

    template <class Archive>
    void load(Archive& ar, const unsigned int);

    double GetEmissionCurrent(int flag);

    template <class gridDataType>
    void
    InitEmissionBoundary(std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
                         const std::shared_ptr<gridDataType>&                         grid,
                         std::vector<double> parametersIn, std::string& errorMsg);

    void
    InitEmissionBoundary(std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
                         const std::shared_ptr<GridData2dpolar<PointType>>&           gridData,
                         std::vector<double> parametersIn, std::string& errorMsg);
    void InitEmissionBoundary(const std::shared_ptr<ElectrodeCurrent<PointType>>& electrode,
                              const std::shared_ptr<GridData2dpolar<PointType>>&  gridData,
                              std::vector<double> parametersIn, std::string& errorMsg);

    // void InitEmissionBoundary(std::vector<BoundaryContainer2d <PointType>>
    // boundaryIn,const std::shared_ptr<GridData2daxs<PointType>>& grid);

    std::vector<std::vector<float>> GetCurrentDensityDistribution();

    std::vector<std::vector<double>> GetEmitterField();
};

#endif