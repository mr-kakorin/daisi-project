#include "Device.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "DataTypes.h"
#include "ElectrodeCurrent.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice2daxs.h"
#include "Geom.h"
#include "GridData.h"
#include "MeshContainer2d.h"
#include "Particle.h"
#include "ParticlesFlow.h"
#include "VTKIncludeSolver.h"
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>

template class DeviceStatus<Particles3dcil<float>, EmitterDevice2daxs<float>, GridData2daxs<float>,
                            float, BoundaryContainer2d<float>, MeshContainer2d<float>>;

template class DeviceStatus<Particles3dcil<double>, EmitterDevice2daxs<double>,
                            GridData2daxs<double>, double, BoundaryContainer2d<double>,
                            MeshContainer2d<double>>;

template class DeviceStatus<Particles2d<float>, EmitterDevice2d<float>, GridData2d<float>, float,
                            BoundaryContainer2d<float>, MeshContainer2d<float>>;

template class DeviceStatus<Particles2d<double>, EmitterDevice2d<double>, GridData2d<double>,
                            double, BoundaryContainer2d<double>, MeshContainer2d<double>>;

template class DeviceStatus<Particles2dpolar<float>, EmitterDevice2d<float>, GridData2dpolar<float>,
                            float, BoundaryContainer2d<float>, MeshContainer2d<float>>;
template class DeviceStatus<Particles2dpolar<double>, EmitterDevice2d<double>,
                            GridData2dpolar<double>, double, BoundaryContainer2d<double>,
                            MeshContainer2d<double>>;

template class DeviceStatus<Particles3d<double>, EmitterDevice3d<double>, GridData3d<double>,
                            double, BoundaryContainer2d<double>, MeshContainer3d<double>>;
template class DeviceStatus<Particles3d<float>, EmitterDevice3d<float>, GridData3d<float>, float,
                            BoundaryContainer2d<float>, MeshContainer3d<float>>;

template class DeviceStatus<Particles3d<double>, EmitterDevice3d<double>, GridData3d<double>,
                            double, BoundaryContainer3d<double>, MeshContainer3d<double>>;

template void
device2daxsfloat::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                             const unsigned int file_version);
template void
device2daxsdouble::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);
template void
device2dfloat::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                          const unsigned int file_version);
template void
device2ddouble::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                           const unsigned int file_version);
template void
device2dpolarfloat::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);
template void
device2dpolardouble::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                const unsigned int file_version);
template void
device3dExtrdouble::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                               const unsigned int file_version);
template void
device3dExtrfloat::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);
template void
device3ddouble::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                           const unsigned int file_version);

template void
device2daxsfloat::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                             const unsigned int file_version);
template void
device2daxsdouble::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                              const unsigned int file_version);
template void
device2dfloat::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                          const unsigned int file_version);
template void
device2ddouble::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                           const unsigned int file_version);
template void
device2dpolarfloat::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                               const unsigned int file_version);
template void
device2dpolardouble::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                                const unsigned int file_version);
template void
device3dExtrdouble::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                               const unsigned int file_version);
template void
device3dExtrfloat::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                              const unsigned int file_version);
template void
device3ddouble::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                           const unsigned int file_version);

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
template <class Archive>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::save(Archive& ar,
                                                                  const unsigned int) const
{
    //
    ar& particlesFlowStatus;
    ar& boundaryConditions;
    ar& mesh;
    ar& boundaries;
    ar& domainBoundaries;
    ar& conductorList;
    ar& domainBoundary;
    ar& boundariesForFlows;
    ar& globalFieldConditions;
    ar& gridData;
    ar& isFieldSimulated;
    ar& isFieldSolverInit;
}
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
template <class Archive>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::load(Archive& ar, const unsigned int)
{
    // ar & gridData;
    ar& particlesFlowStatus;
    ar& boundaryConditions;
    ar& mesh;
    ar& boundaries;
    ar& domainBoundaries;
    ar& conductorList;
    ar& domainBoundary;
    ar& boundariesForFlows;
    ar& globalFieldConditions;
    ar& gridData;
    ar& isFieldSimulated;
    ar& isFieldSolverInit;
    isFieldSolverInit = 0;
}

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
bool DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::isTimeVarying()
{

    for (int i = 0; i < boundaryConditions->GetNumberProperties(); i++)
    {
        if (boundaryConditions->GetPotentialFrequency(i))
            return true;
    };
    return false;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::initParallel(int numThreads,
                                                                          int flagRestart,
                                                                          int memorySize,
                                                                          int solverType,
                                                                          int blockSize)
{
    std::vector<unsigned int> t;
    gridData->InitParallel(numThreads);
    for (int i = 0; i < conductorList.size(); i++)
    {
        conductorList[i]->InitParallel(numThreads);
        if (flagRestart)
            conductorList[i]->ResetPower();
    };
    for (int i = 0; i < particlesFlowStatus.size(); i++)
    {
        particlesFlowStatus[i]->InitParallelTmpBlocks(numThreads, blockSize, solverType);
        particlesFlowStatus[i]->GetDynamicsDataStart()->clear();
        if (flagRestart)
        {
            particlesFlowStatus[i]->clearParallel();
        };
        particlesFlowStatus[i]->InitParallel(memorySize, numThreads, solverType);
    }
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<double> DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                                 BoundaryContainerType, MeshContainerType>::GetMasses()
{
    std::vector<double> result;
    for (int i = 0; i < particlesFlowStatus.size(); i++)
        result.push_back(particlesFlowStatus[i]->GetMass());

    return result;
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<double> DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                                 BoundaryContainerType, MeshContainerType>::GetCharges()
{
    std::vector<double> result;
    for (int i = 0; i < particlesFlowStatus.size(); i++)
        result.push_back(particlesFlowStatus[i]->GetCharge());

    return result;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType,
                  MeshContainerType>::SetglobalFieldConditions(const std::vector<double>& in)
{
    if (globalFieldConditions != in)
        isFieldSimulated = 0;

    globalFieldConditions = in;
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<double>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetglobalFieldConditions()
{
    return globalFieldConditions;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
long long DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                       BoundaryContainerType, MeshContainerType>::GetMemorySize()
{
    volatile long long result = 0;
    for (int i = 0; i < particlesFlowStatus.size(); i++)
        result = result + particlesFlowStatus[i]->GetMemorySize();

    return result;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::ChangeGeom(double p1, double p2)
{
    isFieldSimulated = 0;
    for (int i = 0; i < boundariesForFlows.size(); i++)
    {
        for (int j = 0; j < boundariesForFlows[i]->EdgesData.size(); j++)
        {
            if (std::abs(boundariesForFlows[i]->EdgesData[j].point1.x - p1) < 1e-7)
                boundariesForFlows[i]->EdgesData[j].point1.x = p2;
            if (std::abs(boundariesForFlows[i]->EdgesData[j].point1.y - p1) < 1e-7)
                boundariesForFlows[i]->EdgesData[j].point1.y = p2;

            if (std::abs(boundariesForFlows[i]->EdgesData[j].point2.x - p1) < 1e-7)
                boundariesForFlows[i]->EdgesData[j].point2.x = p2;
            if (std::abs(boundariesForFlows[i]->EdgesData[j].point2.y - p1) < 1e-7)
                boundariesForFlows[i]->EdgesData[j].point2.y = p2;
        }
        boundariesForFlows[i]->ConvertBoundary2VTKUnstructuredGrid();
    }

    for (int i = 0; i < boundaries.size(); i++)
    {
        for (int j = 0; j < boundaries[i]->EdgesData.size(); j++)
        {
            if (std::abs(boundaries[i]->EdgesData[j].point1.x - p1) < 1e-7)
                boundaries[i]->EdgesData[j].point1.x = p2;
            if (std::abs(boundaries[i]->EdgesData[j].point1.y - p1) < 1e-7)
                boundaries[i]->EdgesData[j].point1.y = p2;

            if (std::abs(boundaries[i]->EdgesData[j].point2.x - p1) < 1e-7)
                boundaries[i]->EdgesData[j].point2.x = p2;
            if (std::abs(boundaries[i]->EdgesData[j].point2.y - p1) < 1e-7)
                boundaries[i]->EdgesData[j].point2.y = p2;
        }
        boundaries[i]->ConvertBoundary2VTKUnstructuredGrid();
    }
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType,
                  MeshContainerType>::AddSetOfPotentials(std::string filename)
{
    isFieldSimulated = 0;

    FILE* fp = fopen(filename.c_str(), "r");

    std::string      filenameTmp;
    char             ss[250];
    std::vector<int> boundariesTmp;
    int              tmp;

    int i = 0;
    while (fgets(ss, 250, fp))
    {
        filenameTmp = ss;

        if (filenameTmp.back() == '\n')
            filenameTmp.pop_back();

        boundaryConditions->AddPropertyCondition(std::string("potential"), 0);
        boundaryConditions->SetConditionPropertiesFromFile(i, filenameTmp);

        boundariesTmp.clear();

        fgets(ss, 250, fp);
        int j = 0;

        int jend;

        if (ss[strlen(ss) - 1] == '\n')
            jend = strlen(ss) - 2;
        else
            jend = strlen(ss) - 1;

        jend = strlen(ss);
        while (j < jend)
        {
            int rees = sscanf(ss + j, "%d", tmp);
            if (rees > 0)
                boundariesTmp.push_back(tmp);
            j = j + 2;
        }

        boundaryConditions->SetPropertyConditionsBoundariesList(i, boundariesTmp);

        i++;
    };
    fclose(fp);
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::ResetFlows()
{
    particlesFlowStatus.clear();
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::DeleteFlow(int number)
{
    particlesFlowStatus.erase(particlesFlowStatus.begin() + number);
}

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetElectrodes()
{
    return conductorList;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<std::shared_ptr<EmitterDeviceType>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetEmittersVector()
{
    std::vector<std::shared_ptr<EmitterDeviceType>> result(particlesFlowStatus.size());
    for (int i = 0; i < particlesFlowStatus.size(); i++)
    {
        result[i] = particlesFlowStatus[i]->GetEmitterDevice();
    };
    return result;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType,
                  MeshContainerType>::SetEmitterBoundariesList(int flowNumber, std::vector<int> in,
                                                               std::vector<double> params,
                                                               std::string&        errorMsg)
{
    std::vector<std::shared_ptr<BoundaryContainerType>> bound;

    if (GetFlow(flowNumber)->GetDistributionStyle() == 0)
    {
        for (int i = 0; i < in.size(); i++)
            bound.push_back(boundaries[in[i]]);
    }
    if (GetFlow(flowNumber)->GetDistributionStyle() == 1)
    {
        std::shared_ptr<BoundaryContainerType> boundary =
            std::shared_ptr<BoundaryContainerType>(new BoundaryContainerType());

        boundary->EdgesData = conductorList[in[0]]->ElectrodeEdges;
        if (!boundary->EdgesData.size())
            return;
        bound.push_back(boundary);
        GetFlow(flowNumber)->GetEmitterDevice()->SetGetAssignedElectrode(conductorList[in[0]]);
    }
    GetFlow(flowNumber)
        ->GetEmitterDevice()
        ->SetBoundariesList(in, params, errorMsg, bound, gridData);
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::AddBoundary(std::string filename,
                                                                         std::string& errorMsg)
{
    std::shared_ptr<BoundaryContainerType> tmp;
    try
    {
        tmp = std::shared_ptr<BoundaryContainerType>(new BoundaryContainerType(filename, errorMsg));
        if (errorMsg.size())
            throw 1;
    }
    catch (int& a)
    {
        return;
    }
    boundaries.push_back(tmp);

    boundariesForFlows.push_back(
        std::shared_ptr<BoundaryContainerType>(new BoundaryContainerType(filename, errorMsg)));
    // boundariesForFlows.back().RemoveFineEdges(0.0001);

    float k = 0.001;

    float kloc = k;
    while (1)
    {
        boundariesForFlows.back()->RemoveFineEdges(kloc);
        if (boundariesForFlows.back()->EdgesData.size() < 0.2 * boundaries.back()->EdgesData.size())
        {
            *boundariesForFlows.back() = *boundaries.back();
            kloc                       = kloc / 2;
        }
        else
            break;
    }
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::AddBoundaries(std::string filename,
                                                                           std::string& errorMsg)
{
    FILE* fp = fopen(filename.c_str(), "r");

    std::string filenameTmp;
    char        ss[250];

    int i = 0;
    while (fgets(ss, 250, fp))
    {
        filenameTmp = ss;

        if (filenameTmp.back() == '\n')
            filenameTmp.pop_back();

        std::string errorMsgLoc;
        AddBoundary(filenameTmp, errorMsgLoc);
        if (errorMsgLoc.size())
        {
            errorMsg = errorMsgLoc + " Error string in config is " + std::to_string(i);
            break;
        };
        int i = 0;
    };

    fclose(fp);
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<int> DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                              BoundaryContainerType, MeshContainerType>::GetBoundariesList()
{
    std::vector<int> allboundariesList;
    for (int i = 0; i < boundaries.size(); i++)
        allboundariesList.push_back(i);
    return allboundariesList;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<std::vector<int>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetConductorsList()
{
    std::vector<std::vector<int>> result;
    for (int i = 0; i < conductorList.size(); i++)
        result.push_back(conductorList[i]->boundaryNumbers);
    return result;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType,
                  MeshContainerType>::SetConductorsList(std::vector<int> list, int number, double l,
                                                        std::string& errorMsg)
{
    conductorList[number]->SetBoundaryList(list, boundaries, errorMsg);
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::AddConductor()
{
    conductorList.push_back(
        std::shared_ptr<ElectrodeCurrent<PointType>>(new ElectrodeCurrent<PointType>()));
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<int> DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                              BoundaryContainerType, MeshContainerType>::GetDomainBoundaryList()
{
    return domainBoundaries;
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType,
                  MeshContainerType>::SetDomainBoundaryList(const std::vector<int>& in)
{
    domainBoundaries = in;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::MergeDomainBoundary()
{
    *domainBoundary = *boundaries[domainBoundaries[0]];
    for (int i = 1; i < domainBoundaries.size(); i++)
        domainBoundary->Merge(boundaries[domainBoundaries[i]]);
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::shared_ptr<BoundaryContainerType>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetDomainBoundary()
{
    return domainBoundary;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<std::shared_ptr<BoundaryContainerType>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetboundariesForFlows()
{
    return boundariesForFlows;
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<std::shared_ptr<BoundaryContainerType>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::Getboundaries()
{
    return boundaries;
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::shared_ptr<MeshContainerType>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::Getmesh()
{
    return mesh;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
vtkSmartPointer<vtkUnstructuredGrid>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetDomainBoundaryVTKUnstructuredGrid()
{
    if (domainBoundaries.size() == 0)
        return NULL;
    return GetDomainBoundary()->GetBoundaryVTKUnstructuredGrid();
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::shared_ptr<GridDataType>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetGridData()
{
    return gridData;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::shared_ptr<BoundaryConditions>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetboundaryConditions()
{
    return boundaryConditions;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::vector<std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetFlows()
{
    return particlesFlowStatus;
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::GetFlow(int flowNumber)
{
    return particlesFlowStatus[flowNumber];
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
int DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                 BoundaryContainerType, MeshContainerType>::GetNumberParticlesFlows()
{
    return int(particlesFlowStatus.size());
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::AddFlow(int ParticleType,
                                                                     int DistributionStyle)
{
    particlesFlowStatus.push_back(
        std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>(
            new ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>(
                ParticleType, DistributionStyle, GetBoundariesList())));
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::AddFlow(int ParticleType,
                                                                     int    DistributionStyle,
                                                                     double massNumber,
                                                                     double chargeNumber)
{
    particlesFlowStatus.push_back(
        std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>(
            new ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>(
                ParticleType, DistributionStyle, massNumber, chargeNumber, GetBoundariesList())));
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType, BoundaryContainerType,
             MeshContainerType>::DeviceStatus()
{
    gridData           = std::shared_ptr<GridDataType>(new GridDataType());
    mesh               = std::shared_ptr<MeshContainerType>(new MeshContainerType());
    domainBoundary     = std::shared_ptr<BoundaryContainerType>(new BoundaryContainerType());
    boundaryConditions = std::shared_ptr<BoundaryConditions>(new BoundaryConditions());
    globalFieldConditions.resize(4);
    isFieldSolverInit = 0;
    isFieldSimulated  = 0;
};
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType,
                  MeshContainerType>::GenerateParticles(std::vector<unsigned int> EmptyPlaces,
                                                        int                       flagClear){
    /*	gridData->densityReset();
    for (int i = 0; i < particlesFlowStatus.size(); i++)
    {
    particlesFlowStatus[i].GenerateParticles(EmptyPlaces, flagClear, 0, gridData);
    };*/
};

template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
void DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType, PointType,
                  BoundaryContainerType, MeshContainerType>::SyncronizeThreads(int step, double dt)
{
    gridData->Summrho();

    for (int i = 0; i < conductorList.size(); i++)
        conductorList[i]->PowerAndCurrentsCalculate(step, GetFlow(0)->GetDynamicsData(0)->Time /
                                                              (LIGHT_VELOCITY()),
                                                    dt, gridData->GetType());
};
