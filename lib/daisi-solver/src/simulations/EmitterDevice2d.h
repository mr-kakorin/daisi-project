#ifndef EMITTERDEVICE2D_H
#define EMITTERDEVICE2D_H
#include "EmitterDeviceBase.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

template <class PointType>
class EmitterDevice2d : public EmitterDeviceBase<PointType>
{
    friend class boost::serialization::access;

  private:
    int    nParticlesXY; //����� ������������ ������ � ������ ����� ������ � ������ �����
    double phiXY;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    void GenerateSyncParticle(const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
                              PointType                                           restMass);
    void GenerateSyncParticle(const std::shared_ptr<Particles2d<PointType>>& particlesData,
                              PointType                                      restMass);
    std::vector<double> GetAdditionalSourceInf();
    void SetAdditionalSourceInf(std::vector<double> inf);
    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                  error,
                           std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaryIn,
                           const std::shared_ptr<GridData2d<PointType>>&                 grid);
    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                  error,
                           std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaryIn,
                           const std::shared_ptr<GridData2dpolar<PointType>>&            grid);
    void PreliminaryGeneration(const std::shared_ptr<Particles2dpolar<PointType>>& particlesDataZ0,
                               PointType                                           restMass);
    void PreliminaryGeneration(const std::shared_ptr<Particles2d<PointType>>& particlesDataZ0,
                               PointType                                      restMass);
    void GenerateParticles(std::vector<int> indesexPerThread, int thread, int numThreads,
                           std::vector<unsigned int>&                     EmptyPlaces,
                           const std::shared_ptr<Particles2d<PointType>>& particlesData,
                           PointType restMass, double charge, int flagClear, double dt,
                           int stepNumber, const std::shared_ptr<GridData2d<PointType>>& grid,
                           int flagLocate, int flagDistrib);
    void GenerateParticlesOld(std::vector<int> indesexPerThread, int thread, int numThreads,
                              std::vector<unsigned int>&                     EmptyPlaces,
                              const std::shared_ptr<Particles2d<PointType>>& particlesData,
                              PointType restMass, double charge, int flagClear, double dt,
                              int stepNumber, const std::shared_ptr<GridData2d<PointType>>& grid,
                              int flagLocate, int flagDistrib);
    void GenerateParticles(std::vector<int> indesexPerThread, int thread, int numThreads,
                           std::vector<unsigned int>&                          EmptyPlaces,
                           const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
                           PointType restMass, double charge, int flagClear, double dt,
                           int stepNumber, const std::shared_ptr<GridData2dpolar<PointType>>& grid,
                           int flagLocate, int flagDistrib);

    void GenerateParticlesLinac(int flagTest, int thread, int numThreads,
                                std::vector<unsigned int>&                     EmptyPlaces,
                                const std::shared_ptr<Particles2d<PointType>>& particlesData,
                                PointType restMass, short chargeSign, int flagClear, double dt,
                                int stepNumber, const std::shared_ptr<GridData2d<PointType>>& grid,
                                int flagLocate, int flagDistrib);
    void GenerateParticlesLinac(int flagTest, int thread, int numThreads,
                                std::vector<unsigned int>&                          EmptyPlaces,
                                const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
                                PointType restMass, short chargeSign, int flagClear, double dt,
                                int                                                stepNumber,
                                const std::shared_ptr<GridData2dpolar<PointType>>& grid,
                                int flagLocate, int flagDistrib);

    double GetEmissionCurrent();
    EmitterDevice2d(int DistributionStyleIn);
    EmitterDevice2d();
    std::vector<std::vector<float>> GetCurrentDensityDistribution();
};
#endif