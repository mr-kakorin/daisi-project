#ifndef EMITTERDEVICE2DGPU_H
#define EMITTERDEVICE2DGPU_H
#include "EmitterDevice2d.h"
#include "EmitterDeviceBaseGPU.h"
#include "ParticleGPU.h"
#include "ParticleSourceGPU.h"
template <class PointType> class EmitterDevice2dGPU : public EmitterDeviceBase2dGPU<PointType>
{
    friend class boost::serialization::access;

  private:
    int    nParticlesXY; //число эмиттируемых частиц с разным углом вылета в каждой точке
    double phiXY;

  public:
    EmitterDevice2dGPU(){

    };
    template <class PointType>
    EmitterDevice2dGPU(EmitterDevice2d<PointType>& obj)
        : EmitterDeviceBase2dGPU<PointType>(obj.GetParticlesNumber(), obj.GetDistributionParameters(),
                                            obj.GetEmissionType(), obj.GetEmissionCurrent(), obj.GetParticleSource())
    {
        particleSource           = *obj.GetParticleSource();
        std::vector<int>    res  = obj.GetParticlesNumber();
        std::vector<double> res1 = obj.GetDistributionParameters();
        nParticlesXY             = res[2];
        phiXY                    = res1[2];
    };

    /*	int GetNumbersOfParticlesGeneration()
            {
                    return nParticlesEmitter*nParticlesEnergy*nParticlesXY;
            };
            void GenerateParticles(const std::vector<unsigned int>& EmptyPlaces, Particles2d<PointType>& particlesData,
       PointType restMass, int flagClear); void InitEmissionBoundary(std::vector<BoundaryContainer2d <PointType>>
       boundaryIn); double GetEmissionCurrent(); void SetEmissionCurrent(double current); EmitterDevice2d();
            std::vector<int> GetParticlesNumber();
            void SetParticlesNumber(std::vector<int> in);
            void SetDistributionParameters(std::vector<double> in);
            std::vector<double> GetDistributionParameters();
            std::vector<std::vector<float>> GetCurrentDensityDistribution();*/
};
#endif