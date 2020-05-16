#ifndef PARTICLESFLOWGPU_H
#define PARTICLESFLOWGPU_H
#include "BoundaryConditions.h"
#include "EmitterDeviceBase.h"
#include "Particle.h"
#include "ParticleGPU.h"
#include "ParticlesFlow.h"
#include "boost/archive/text_iarchive.hpp"
#include "boost/archive/text_oarchive.hpp"
template <class ParticlesDataType, class EmitterDeviceType, class PointType> class ParticlesFlowGPU
{
    int               ParticlesType;
    double            mass;
    double            charge;
    ParticlesDataType DynamicsData;
    EmitterDeviceType EmitterDeviceStatus;

  public:
    BoundaryConditions boundaryConditions;
    ParticlesFlowGPU(){

    };

    template <class ParticlesDataTypeCPU, class EmitterDeviceTypeCPU>
    ParticlesFlowGPU(ParticlesFlow<ParticlesDataTypeCPU, EmitterDeviceTypeCPU, PointType>& obj)
    {
        ParticlesType       = obj.GetParticlesType();
        mass                = obj.GetMass();
        charge              = obj.GetCharge();
        boundaryConditions  = obj.boundaryConditions;
        DynamicsData        = *obj.GetDynamicsData();
        EmitterDeviceStatus = obj.EmitterDeviceStatus;
    };

    /*
    BoundaryConditions* GetboundaryConditions()
    {
            return &boundaryConditions;
    };
    


    PointType GetMass()
    {
            return mass;
    };
    PointType GetCharge()
    {
            return charge;
    };
    PointType GetAlpha()commtools::LIGHT_VELOCITY()commtools::LIGHT_VELOCITY()
    {
            return charge / (mass*Dconst::LIGTH_VELOCITY*Dconst::LIGTH_VELOCITY);
    };
    //EmitterDeviceType*  EmitterDeviceStatus;
    EmitterDeviceType*  EmitterDeviceStatus;

    EmitterDeviceType* GetEmitterDevice()
    {
            return EmitterDeviceStatus;
    };

    {
            return &boundaryConditions;
    };

    ParticlesDataType* GetDynamicsData(){
            return &DynamicsData;
    }
    std::vector<float> GetData(int Index)
    {
            return DynamicsData.GetData(Index, mass, charge);
    };
    std::vector<void*> GetData()
    {
            return DynamicsData.GetData();
    };
    ParticlesFlow()
    {
            EmitterDeviceStatus = new EmitterDeviceType();
    };
    ParticlesFlow(int ParticlesTypeIn, int MassNumber, int ChargeNumber, std::vector<int> boundariesList) //�����������
    ��� �����
    {
            EmitterDeviceStatus = new EmitterDeviceType();
            ParticlesType = ParticlesTypeIn;
            mass = commtools::PROTON_MASS()*MassNumber;
            charge = commtools::ELECTRON_CHARGE()*ChargeNumber;
            boundaryConditions.SetDefaultConditionsList(boundariesList);
    };
    ParticlesFlow(int ParticlesTypeIn, std::vector<int> boundariesList) //����������� ����������
    {
            EmitterDeviceStatus = new EmitterDeviceType();
            ParticlesType = ParticlesTypeIn;
            mass = commtools::ELECTRON_MASS();
            charge = commtools::ELECTRON_CHARGE();
            boundaryConditions.SetDefaultConditionsList(boundariesList);
    };
    std::vector<double> GetFlowProperties()
    {
            std::vector<double> properties;
            properties.push_back(double(ParticlesType));
            properties.push_back(mass);
            properties.push_back(charge);
            return properties;
    };
    std::vector<unsigned int> GetNewIndexes(std::vector<unsigned int> EmptyPlaces)
    {
            int particlesNumber = EmitterDeviceStatus->GetNumbersOfParticlesGeneration();
            std::vector<unsigned int> result(particlesNumber);


            if (particlesNumber<EmptyPlaces.size())
            {
                    for (int i = 0; i < particlesNumber; i++)
                            result[i] = EmptyPlaces[i];
                    return result;
            }

            int k = 0;

            for (int i = 0; i < EmptyPlaces.size(); i++)
                    result[i] = EmptyPlaces[i];

            for (int i = int(EmptyPlaces.size()); i < particlesNumber; i++)
            {
                    result[i] = DynamicsData.NParticles + k;
                    k++;
            }
            return result;
    };
    void GenerateParticles(std::vector<unsigned int> EmptyPlaces, int flagClear)
    {
            EmitterDeviceStatus->GenerateParticles(EmptyPlaces, DynamicsData, mass, flagClear);
    };*/
};
#endif