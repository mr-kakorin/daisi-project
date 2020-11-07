#ifndef EMITTERDEVICEINTERFACE_H
#define EMITTERDEVICEINTERFACE_H
#include "GridData.h"
#include <algorithm>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <vector>
#include <daisi-solver/ModelInterface.h>

template <class PointType>
class BoundaryContainer2d;
template <class PointType>
class BoundaryContainer3d;
template <class PointType>
class GridData2d;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;
template <class PointType>
class GridData3dpolar;
template <class PointType>
class Particles2dpolar;
template <class PointType>
class Particles3d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class Particles3dcil;
template <class PointType>
class Particles2d;

template <class PointType>
class EmitterDeviceBaseGPU;
template <class PointType>
class ParticleSource2d;

template <class PointType>
class ElectrodeCurrent;

template <class PointType>
class myvector;

using EnergyDistribution = double(*)();

template <class PointType>
class EmitterDeviceBase
{
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    friend class boost::serialization::access;
    friend class EmitterDeviceBaseGPU<PointType>;

    // Linac style
    // emission type  NumbersParams[0]
    // emit period  NumbersParams[1]
    // particles on bunch NumbersParams[2]

    // Current  DistribParams[0]
    // Energy  DistribParams[1]
    // Energy spread DistribParams[2]
    // bunch frequency DistribParams[3]

    // Source style
    // emission type  NumbersParams[0]
    // emit period  NumbersParams[1]
    // particles on Energy NumbersParams[2]
    // particles on Emitter NumbersParams[3]

    // Energy  DistribParams[0]

    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    std::shared_ptr<ElectrodeCurrent<PointType>> Assignedelectrode;

  protected:
    std::vector<int>    NumbersParams;
    std::vector<int>    flagsParams;
    std::vector<double> DistribParams;

    std::vector<double>                          startPoint;
    std::vector<double>                          endPoint;
    std::shared_ptr<ParticleSource2d<PointType>> particleSource;
    std::vector<int>                             boundaryList;
    std::vector<double>                          emitterInitParam;

    double                           emissionCurrent;
    std::vector<std::vector<double>> Data;

  public:
	Daisi::Emission::EnergyDistributionType             energy_distribution{ Daisi::Emission::EnergyDistributionType_Bimodal };
    std::vector<double>                                 GetEmitterInitParameters();
    const std::shared_ptr<ElectrodeCurrent<PointType>>& GetAssignedElectrode();
    void SetGetAssignedElectrode(const std::shared_ptr<ElectrodeCurrent<PointType>>& in);
    std::vector<double>              GetDistribParams();
    double                           GetLambda();
    int                              GetMaxParticles();
    int                              GetNumbersOfParticlesGeneration();
    int                              GetEmitPeriod();
    std::vector<std::vector<double>> GetAllParameters();
    void SetAllParameters(const std::vector<std::vector<double>>& In);

    int  GetEmissionType();
    void SetFlowCurrent(double res);
    double                                                    GetSourceSize();
    int                                                       GetnParticlesEmitter();
    int                                                       GetnParticlesBunch();
    double                                                    getErAverage();
    std::vector<std::vector<double>>                       GetEmitterField();
    std::vector<unsigned int>                                 newIndexes;
    std::shared_ptr<ParticleSource2d<PointType>>              GetParticleSource();
    std::vector<std::shared_ptr<ParticleSource2d<PointType>>> GetParticleSources();
    std::vector<int>                                          GetBoundariesList();
    EmitterDeviceBase(int DistributionStyleIn);
    void GenerateEllipses(double restMass);
    void GetSliceIndexes(std::vector<int>& sliceIndexesParallel, int flagTest, double t1, double t2,
                         double& phMin, double& phMax, int numThreads, int thread);
    EmitterDeviceBase() = default;
    void SetDirectionPoints(std::vector<double> sP, std::vector<double> eP); //
    std::vector<std::vector<double>> GetDirectionPoints();                   //
};
#endif