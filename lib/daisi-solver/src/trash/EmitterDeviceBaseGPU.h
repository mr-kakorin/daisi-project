#ifndef EMITTERDEVICEINTERFACEGPU_H
#define EMITTERDEVICEINTERFACEGPU_H
#include "BoundaryContainer2d.h"
#include "ParticleSourceGPU.h"
#include <algorithm>
#include <vector>
template <class PointType> class EmitterDeviceBase2dGPU
{
    friend class boost::serialization::access;

  protected:
    ParticleSource2dGPU<PointType> particleSource;
    int                            emissionType;
    int                            distributionEnergyType; // тип энергетического распределения
    int                            generateType;           // тип распределения (случайное/неслучайное)
    int                            nParticlesEnergy;       // число эмиттируемых частиц с различной энергией
    double                         energyAverage;          // средняя энергия эмиттируемых частиц
    double                         energySpread;           // энергетический разброc эмиттируемых частиц
    int                            nParticlesEmitter;      //число эмиттируемых частиц с разным положением на эмиттере
    double                         totalCurrent;
    std::vector<int>               boundaryList;

  public:
    EmitterDeviceBase2dGPU()
    {
    }
    EmitterDeviceBase2dGPU(std::vector<int> pNumbers, std::vector<double> param, int emissionTypeIn,
                           double totalCurrentIn, ParticleSource2d<PointType>* obj)
    {
        particleSource = *obj;
        emissionType   = emissionTypeIn;

        totalCurrent = totalCurrentIn;

        generateType           = 0;
        distributionEnergyType = 0;
        nParticlesEnergy       = pNumbers[0];
        nParticlesEmitter      = pNumbers[1];

        energyAverage = param[0];
        energySpread  = param[1];
    };
    /*	std::vector<unsigned int> newIndexes;

            virtual std::vector<std::vector<float>> GetCurrentDensityDistribution() = 0;

            std::vector<int> GetBoundariesList()
            {
                    return boundaryList;
            };
            void SetBoundariesList(std::vector<int>  in, std::vector< BoundaryContainer2d <PointType>> boundaryIn)
            {
                    boundaryList = in;
                    InitEmissionBoundary(boundaryIn);
            };
            virtual void InitEmissionBoundary(std::vector<BoundaryContainer2d <PointType>> boundaryIn) = 0;
            EmitterDeviceInterface()
            {
                    emissionType = 0;
                    distributionEnergyType = 0;
                    generateType = 0;
                    nParticlesEnergy = 1;
                    nParticlesEmitter = 1;
                    energyAverage = 0;
                    energySpread = 0;
            };
            void SetEmissionType(int type)
            {
                    emissionType = type;
            };
            int GetEmissionType()
            {
                    return emissionType;
            };*/
};
#endif