#ifndef DEVICEGPU_H
#define DEVICEGPU_H

#include "BoundaryConditionsGPU.h"
#include "EmitterDevice2dGPU.h"
#include "EmitterDevice2daxsGPU.h"
//#include "Device.h"
#include "GridDataGPU.h"
#include "ParticleGPU.h"
#include "ParticlesFlowGPU.h"
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
class DeviceStatus;
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType>
class DeviceStatusGPU
{
    //	std::vector<ParticlesFlow<ParticlesDataType, ParticleSourceType, EmitterDeviceType,
    // PointType>*>
    // particlesFlowStatus;

    GPUTypes::vector<ParticlesFlowGPU<ParticlesDataType, EmitterDeviceType, PointType>>
                                            particlesFlowStatus;
    GridDataType                            gridData;
    BoundaryConditionsGPU                   boundaryConditions;
    GPUTypes::vector<BoundaryContainerType> boundaries;

  public:
    DeviceStatusGPU(){

    };
    template <class ParticlesDataTypeCPU, class EmitterDeviceTypeCPU, class GridDataTypeCPU,
              class PointTypeCPU, class BoundaryContainerTypeCPU, class MeshContainerTypeCPU>
    DeviceStatusGPU(
        DeviceStatus<ParticlesDataTypeCPU, EmitterDeviceTypeCPU, GridDataTypeCPU, PointTypeCPU,
                     BoundaryContainerTypeCPU, MeshContainerTypeCPU>& DeviceStatusCPU)
    {
        boundaryConditions  = DeviceStatusCPU.boundaryConditions;
        gridData            = DeviceStatusCPU.gridData;
        particlesFlowStatus = DeviceStatusCPU.GetFlows();
        boundaries          = DeviceStatusCPU.boundaries;
    };

    DeviceStatusGPU* AllocateOnGPU()
    {
        DeviceStatusGPU<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                        BoundaryContainerType>* result;
        cudaMalloc((void**)&result, sizeof(*this));
        cudaMemcpy(result, this, sizeof(*this), cudaMemcpyHostToDevice);
        return result;
    };

    /*GridDataType* GetGridData(){
            return &gridData;
    };

    BoundaryConditions* GetboundaryConditions()
    {
            return &boundaryConditions;
    };

    ParticlesFlow<ParticlesDataType, ParticleSourceType, EmitterDeviceType, PointType>* GetFlow(int
    flowNumber)
    {
            return particlesFlowStatus[flowNumber];
    };

    void GenerateParticles(std::vector<unsigned int> EmptyPlaces, int flagClear)
    {
            for (int i = 0; i < particlesFlowStatus.size(); i++)
            {
                    particlesFlowStatus[i]->GenerateParticles(EmptyPlaces, flagClear);
            };

    };*/
};

#endif