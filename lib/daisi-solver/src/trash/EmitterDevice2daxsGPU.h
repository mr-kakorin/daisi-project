#ifndef EMITTERDEVICE2DAXSGPU_H
#define EMITTERDEVICE2DAXSGPU_H
#include "EmitterDevice2daxs.h"
#include "EmitterDeviceBaseGPU.h"
#include "ParticleGPU.h"
#include "ParticleSourceGPU.h"
template <class PointType> class EmitterDevice2daxsGPU : public EmitterDeviceBase2dGPU<PointType>
{
    friend class boost::serialization::access;

  private:
    int    nParticlesRZ;
    int    nParticlesPhiR;
    double phiRZ;
    double phiPhiR;
    double KphPhiR;

  public:
    EmitterDevice2daxsGPU(){

    };
    template <class PointType>
    EmitterDevice2daxsGPU(EmitterDevice2daxs<PointType>& obj)
        : EmitterDeviceBase2dGPU<PointType>(obj.GetParticlesNumber(), obj.GetDistributionParameters(),
                                            obj.GetEmissionType(), obj.GetEmissionCurrent(), obj.GetParticleSource())
    {
        std::vector<int>    res  = obj.GetParticlesNumber();
        std::vector<double> res1 = obj.GetDistributionParameters();
        nParticlesRZ             = res[2];
        nParticlesPhiR           = res[3];
        phiRZ                    = res1[2];
        phiPhiR                  = res1[3];
        KphPhiR                  = res1[4];
    };
    /*	int GetNumbersOfParticlesGeneration()
            {
                    return nParticlesEmitter*nParticlesEnergy*nParticlesRZ*nParticlesPhiR;
            };
            void GenerateParticles(const std::vector<unsigned int>& EmptyPlaces, Particles3dcil<PointType>&
       particlesData, PointType restMass, int flagClear); void InitEmissionBoundary(std::vector<BoundaryContainer2d
       <PointType>> boundaryIn); double GetEmissionCurrent(); void SetEmissionCurrent(double current);
            EmitterDevice2daxs();
            std::vector<int> GetParticlesNumber();
            void SetParticlesNumber(std::vector<int> in);
            void SetDistributionParameters(std::vector<double> in);
            std::vector<double> GetDistributionParameters();
            std::vector<std::vector<float>> GetCurrentDensityDistribution();*/
};

#endif