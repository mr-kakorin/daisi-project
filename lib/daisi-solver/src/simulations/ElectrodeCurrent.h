#ifndef ELECTRODECURRENT_H
#define ELECTRODECURRENT_H
#include <boost/archive/binary_iarchive.hpp>
#include <deque>
#include <vector>

template <class PointType>
class BoundaryContainer2d;
template <class PointType>
class BoundaryContainer3d;
namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
}

template <class PointType>
class ElectrodeCurrent
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    PointType totalCurrent;
    int       n;
    int       mapSize;

  public:
    std::vector<DGeo::Edge<PointType>>  ElectrodeEdges;
    std::deque<std::vector<PointType>>  CurrentAtPointsTimes;
    std::deque<std::vector<PointType>>  EnergyAtPointsTimes;
    std::deque<std::vector<PointType>>  CurrentAtPointsLongTimes;
    std::deque<std::vector<PointType>>  EnergyAtPointsLongTimes;
    std::vector<std::vector<PointType>> CurrentAtPointsThreaded;
    std::vector<std::vector<PointType>> EnergyAtPointsThreaded;

    std::vector<PointType> averageIrradiatedPowerDensity;
    std::vector<PointType> averageIrradiatedCurrentDensity;
    std::vector<PointType> averageCollectedCurrentDensity;
    std::vector<PointType> CurrentAtPointsAverage;
    std::vector<PointType> CurrentAtPointsAverageCollected;
    std::vector<PointType> EnergyAtPointsAverage;
    std::vector<PointType> averageCollectedCurrentDensitySim;

    std::vector<double> parameters;
    //	template <class PointType1>
    // std::vector<std::vector<PointType1>> GetElectrodeValue(int flag);

    std::vector<std::vector<float>> GetElectrodeValueF(int flag);
    std::vector<std::vector<double>> GetElectrodeValueD(int flag);

    std::vector<int> boundaryNumbers;
    double           GetCurrent();
    void SetCurrent(PointType cur);
    std::vector<std::vector<PointType>> PowerAtPointsTmp;
    void InitParallel(int numThreads);
    void SetBoundaryList(std::vector<int>                                              list,
                         std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>>& boundaries,
                         std::string&                                                  errorMsg);
    void SetBoundaryList(std::vector<int>                                              list,
                         std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                         std::string&                                                  errorMsg);
    void AddCharge(PointType charge, const DGeo::Point<PointType>& point, PointType dt,
                   PointType energy, PointType chargeParticle, int thread);
    void ResetCurrent();
    void PowerAndCurrentsCalculate(int k, double t, double dt, int type);
    void ResetPower();
    ElectrodeCurrent();
    std::vector<double> GetParameteres();
    void SetParameteres(const std::vector<double>& input);
};
#endif