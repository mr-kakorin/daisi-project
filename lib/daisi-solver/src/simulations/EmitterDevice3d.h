#ifndef EMITTERDEVICE3D_H
#define EMITTERDEVICE3D_H
#include "EmitterDeviceBase.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

template <class PointType>
class EmitterDevice3d : public EmitterDeviceBase<PointType>
{
    friend class boost::serialization::access;

  private:
    int                 nParticlesXY; //����� ������������ ������ � ������ ����� ������ � ������ �����
    int                 nParticlesZ;
    int                 ndZ;
    double              phiXY;
    double              Z1;
    double              Z2;
    std::vector<double> zArray;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

    std::vector<std::shared_ptr<ParticleSource2d<PointType>>> particleSources;

  public:
    void GenerateSyncParticle(const std::shared_ptr<Particles3d<PointType>>& particlesData,
                              PointType                                      restMass);
    void SetFlowCurrent(double res);
    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                 error,
                           std::vector<std::shared_ptr<BoundaryContainer3d<PointType>>> boundaryIn,
                           const std::shared_ptr<GridData3d<PointType>>&                grid);

    double getErAverage();

    double GetEmissionCurrent();

    double GetSourceSize();

    std::vector<std::shared_ptr<ParticleSource2d<PointType>>>& GetParticleSource();

    std::vector<double> GetAdditionalSourceInf();

    void SetAdditionalSourceInf(std::vector<double> inf);

    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                 error,
                           std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
                           const std::shared_ptr<GridData2d<PointType>>&                grid);

    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                 error,
                           std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
                           const std::shared_ptr<GridData2dpolar<PointType>>&           grid);

    void SetBoundariesList(std::vector<int> in, std::vector<double> parametersIn,
                           std::string&                                                 error,
                           std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>> boundaryIn,
                           const std::shared_ptr<GridData3d<PointType>>&                grid);

    void GenerateParticles(std::vector<int> indesexPerThread, int thread, int numThreads,
                           std::vector<unsigned int>&                     EmptyPlaces,
                           const std::shared_ptr<Particles3d<PointType>>& particlesData,
                           PointType restMass, double charge, int flagClear, double dt,
                           int stepNumber, const std::shared_ptr<GridData3d<PointType>>& grid,
                           int flagLocate, int flagDistrib);

    EmitterDevice3d(){

    };
    EmitterDevice3d(int DistributionStyleIn);

    std::vector<std::shared_ptr<ParticleSource2d<PointType>>>& GetParticleSources();

    std::vector<std::vector<float>> GetCurrentDensityDistribution();

    void GenerateParticlesLinac(int flagTest, int thread, int numThreads,
                                std::vector<unsigned int>&                     EmptyPlaces,
                                const std::shared_ptr<Particles3d<PointType>>& particlesData,
                                PointType restMass, short chargeSign, int flagClear, double dt,
                                int stepNumber, const std::shared_ptr<GridData3d<PointType>>& grid,
                                int flagLocate, int flagDistrib);
    void PreliminaryGeneration(const std::shared_ptr<Particles3d<PointType>>& particlesDataZ0,
                               PointType                                      restMass);
};
#endif