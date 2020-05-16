#ifndef DEVICE_H
#define DEVICE_H

#include "VTKIncludeSolver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>

template <class ParticlesDataType, class EmitterDeviceType, class PointType>
class ParticlesFlow;
class GridDataType;
class MeshContainerType;
template <class PointType>
class ElectrodeCurrent;
class BoundaryConditions;
template <class ParticlesDataType, class EmitterDeviceType, class GridDataType, class PointType,
          class BoundaryContainerType, class MeshContainerType>
class DeviceStatus
{
    std::vector<std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>>
                                                              particlesFlowStatus;
    std::shared_ptr<GridDataType>                             gridData;
    std::vector<std::shared_ptr<BoundaryContainerType>>       boundaries;
    std::vector<std::shared_ptr<BoundaryContainerType>>       boundariesForFlows;
    std::shared_ptr<BoundaryContainerType>                    domainBoundary;
    std::vector<int>                                          domainBoundaries;
    std::shared_ptr<MeshContainerType>                        mesh;
    std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>> conductorList;
    std::shared_ptr<BoundaryConditions>                       boundaryConditions;
    std::vector<double>                                       globalFieldConditions;

    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    int isFieldSimulated;
    int isFieldSolverInit;

    void SetglobalFieldConditions(const std::vector<double>& in);
    std::vector<double> GetglobalFieldConditions();
    void                MergeDomainBoundary();
    std::vector<double> GetMasses();
    std::vector<double> GetCharges();
    long long           GetMemorySize();
    void ChangeGeom(double p1, double p2);
    void AddSetOfPotentials(std::string filename);
    void ResetFlows();
    void DeleteFlow(int number);
    std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>> GetElectrodes();
    std::vector<std::shared_ptr<EmitterDeviceType>>           GetEmittersVector();
    void SetEmitterBoundariesList(int flowNumber, std::vector<int> in, std::vector<double> params,
                                  std::string& errorMsg);
    void AddBoundary(std::string filename, std::string& errorMsg);
    void AddBoundaries(std::string filename, std::string& errorMsg);
    std::vector<int>              GetBoundariesList();
    std::vector<std::vector<int>> GetConductorsList();
    void SetConductorsList(std::vector<int> list, int number, double l, std::string& errorMsg);
    void AddConductor();
    bool isTimeVarying();

    std::vector<int> GetDomainBoundaryList();
    void SetDomainBoundaryList(const std::vector<int>& in);
    std::shared_ptr<BoundaryContainerType>              GetDomainBoundary();
    std::vector<std::shared_ptr<BoundaryContainerType>> GetboundariesForFlows();
    std::vector<std::shared_ptr<BoundaryContainerType>> Getboundaries();
    std::shared_ptr<MeshContainerType>                  Getmesh();
    vtkSmartPointer<vtkUnstructuredGrid>                GetDomainBoundaryVTKUnstructuredGrid();
    std::shared_ptr<GridDataType>                       GetGridData();
    std::shared_ptr<BoundaryConditions>                 GetboundaryConditions();

    std::vector<std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>>
    GetFlows();

    std::shared_ptr<ParticlesFlow<ParticlesDataType, EmitterDeviceType, PointType>>
    GetFlow(int flowNumber);

    int  GetNumberParticlesFlows();
    void AddFlow(int ParticleType, int DistributionStyle);
    void AddFlow(int ParticleType, int DistributionStyle, double massNumber, double chargeNumber);
    DeviceStatus();
    void GenerateParticles(std::vector<unsigned int> EmptyPlaces, int flagClear);
    void initParallel(int numThreads, int flagRestart, int memorySize, int solverType,
                      int blockSize);
    void SyncronizeThreads(int step, double dt);
};
#endif