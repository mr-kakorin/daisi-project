#ifndef EMITTERDEVICE2DAXS_H
#define EMITTERDEVICE2DAXS_H
#include "EmitterDeviceBase.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

template <class PointType>
class EmitterDevice2daxs : public EmitterDeviceBase<PointType>
{
    friend class boost::serialization::access;

  private:
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<EmitterDeviceBase<PointType>>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<EmitterDeviceBase<PointType>>(*this);
    }

  public:
    void PreliminaryGeneration(const std::shared_ptr<Particles3dcil<PointType>>& particlesDataZ0,
                               PointType                                         restMass);

    void GenerateSyncParticle(const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
                              PointType                                         restMass);

    std::vector<double> GetAdditionalSourceInf();
    void SetAdditionalSourceInf(std::vector<double> inf);

    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                 error,
                           std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
                           const std::shared_ptr<GridData2daxs<PointType>>&             grid);
    void GenerateParticles(std::vector<int> indesexPerThread, int thread, int numThreads,
                           std::vector<unsigned int>&                        EmptyPlaces,
                           const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
                           PointType restMass, double charge, int flagClear, double dt,
                           int stepNumber, const std::shared_ptr<GridData2daxs<PointType>>& grid,
                           int flagLocate, int flagDistrib);
    void GenerateParticlesOld(std::vector<int> indesexPerThread, int thread, int numThreads,
                              std::vector<unsigned int>&                        EmptyPlaces,
                              const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
                              PointType restMass, double charge, int flagClear, double dt,
                              int stepNumber, const std::shared_ptr<GridData2daxs<PointType>>& grid,
                              int flagLocate, int flagDistrib);
    void GenerateParticlesLinac(int flagTest, int thread, int numThreads,
                                std::vector<unsigned int>&                        EmptyPlaces,
                                const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
                                PointType restMass, short chargeSign, int flagClear, double dt,
                                int                                              stepNumber,
                                const std::shared_ptr<GridData2daxs<PointType>>& grid,
                                int flagLocate, int flagDistrib);

    double GetEmissionCurrent();
    EmitterDevice2daxs(int DistributionStyleIn);
    EmitterDevice2daxs(){

    };
    std::vector<std::vector<float>> GetCurrentDensityDistribution();
};

#endif