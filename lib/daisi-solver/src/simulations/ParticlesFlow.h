#ifndef PARTICLESFLOW_H
#define PARTICLESFLOW_H
#include <boost/archive/binary_iarchive.hpp>
//#include "ParticleGPU.h"
class BoundaryConditions;
template <class ParticlesDataType, class EmitterDeviceType, class PointType>
class ParticlesFlow
{
    friend class boost::serialization::access;
    int                                                          ParticlesType;
    int                                                          DistributionStyle;
    double                                                       mass;
    double                                                       charge;
    double                                                       maxTime;
    std::shared_ptr<ParticlesDataType>                           newParticlesZ0;
    std::shared_ptr<ParticlesDataType>                           tmpData;
    std::shared_ptr<ParticlesDataType>                           DynamicsData;
    std::vector<std::shared_ptr<ParticlesDataType>>              DynamicsDataParallel;
    std::vector<std::shared_ptr<ParticlesDataType>>              DynamicsDataParallelTmp;
    std::vector<std::vector<std::shared_ptr<ParticlesDataType>>> Emittances;
    // std::vector<ParticlesDataType> newParticlesParallel;
    std::vector<double>                                 EmittancesList;
    std::vector<std::vector<int>>                       indesexPerThread;
    std::vector<std::vector<std::vector<unsigned int>>> writeIndexes;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    std::shared_ptr<BoundaryConditions> boundaryConditions;
    std::shared_ptr<EmitterDeviceType>  EmitterDeviceStatus;

    void GetEmittanceData(std::vector<std::vector<float>>& data, int flagNumber, int emFlag);

    void clearParallel();
    void AddEmittancesData(int thread);

    void MergeEmittancesData();

    void GetRmsEmittances(std::vector<float>& emittance);

    std::vector<double> GetEmittancesList();
    void AddEmittance(double param);
    void checkEmittances(int thread, int i1, int i2);

    long long GetMemorySize();
    int       GetNumberOfThreads();
    int       GetNumberOfParticles(); // total number of particles in DynamicsDataParallel

    void GenerateSyncParticle();

    void CalculateFlowCurrent();
    int  GetMaxParticles();

    void DivideForThreads(int numThreads);
    void InitParallel(int size, int numThreads, int flagMethod);

    bool isConfigure();

    void InitParallelTmpBlocks(int numThreads, int blockSize, int solverType);
    std::vector<double> GetFlowMCNumbers();

    void SetFlowMCNumbers(std::vector<double> numbers);

    int       GetParticlesType();
    int       GetDistributionStyle();
    PointType GetMass();
    PointType GetCharge();
    PointType GetAlpha();
    // EmitterDeviceType*  EmitterDeviceStatus;

    std::shared_ptr<EmitterDeviceType> GetEmitterDevice();

    std::shared_ptr<BoundaryConditions> GetboundaryConditions();
    std::shared_ptr<ParticlesDataType>  GetDynamicsData();

    std::shared_ptr<ParticlesDataType> GetDynamicsDataStart();

    std::shared_ptr<ParticlesDataType> GetDynamicsDataTmp();
    std::shared_ptr<ParticlesDataType> GetDynamicsData(int thread);
    std::shared_ptr<ParticlesDataType> GetDynamicsDataTmp(int thread);

    std::vector<std::shared_ptr<ParticlesDataType>> GetDynamicsDataParallelArray();

    std::vector<float> GetData(int Index);

    std::vector<void*> GetData();

    ParticlesFlow();

    void ReserveMemory(int size);

    ParticlesFlow(int ParticlesTypeIn, int DistributionStyleIn, double MassNumber,
                  double ChargeNumber,
                  std::vector<int> boundariesList); //����������� ��� �����

    bool CheckTimeLimit(int thread);

    ParticlesFlow(int ParticlesTypeIn, int DistributionStyleIn,
                  std::vector<int> boundariesList); //����������� ����������

    std::vector<double> GetFlowProperties();

    std::vector<unsigned int> GetNewIndexes(std::vector<unsigned int> EmptyPlaces);

    template <class GridDataType>
    void GenerateParticles(std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                           int stepNumber, const std::shared_ptr<GridDataType>& grid,
                           int flagLocate, int flagDistrib);

    template <class GridDataType>
    void GenerateParticlesThreaded(int thread, int numThreads,
                                   std::vector<unsigned int>& EmptyPlaces, int flagClear, double dt,
                                   int stepNumber, const std::shared_ptr<GridDataType>& grid,
                                   int flagLocate, int flagDistrib);

    template <class GridDataType>
    void GenerateParticlesThreadedTest(int thread, int numThreads,
                                       std::vector<unsigned int>& EmptyPlaces, int flagClear,
                                       double dt, int stepNumber,
                                       const std::shared_ptr<GridDataType>& grid, int flagLocate,
                                       int flagDistrib);

    /*template <class GridDataType>
    void GenerateParticlesVector(double dt, int stepNumber)
    {
            newParticles = ParticlesDataType();
            EmitterDeviceStatus->GenerateParticles(std::vector <unsigned int > {}, newParticles,
    mass,
    Dmath::sign(charge), 1, dt, stepNumber);
    };*/
    void setNewParticles(const std::vector<unsigned int>& EmptyPlaces);
    void CopyDynamicsDataToTmp(int i1, int i2);
    void CopyDynamicsDataToTmpThreaded(int thread, int i1, int i2);
};
#endif