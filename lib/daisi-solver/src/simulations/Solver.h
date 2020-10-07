#ifndef SOLVER_H
#define SOLVER_H
#include <boost/serialization/vector.hpp>
#include <mutex>
#include <thread>

template <class PointType>
class ParticleGridInterface;
template <class PointType>
class ParticlesMover;
template <class PointType>
class MeshGenerator;
template <class PointType>
class PoissonSolver;
template <class PointType>
class MagneticFieldSolver;
class SimulationData;
class DynamicsData;
template <class PointType>
class EmissionCurrentSolverPIC;
template <class PointType>
class EmissionCurrentSolverPTI;
template <class PointType>
class particlesFields;
/*class freq
{
        std::vector <double> Lcell;
        std::vector <double> kList;
        std::vector <double> Rleft;
        std::vector <double> Rright;
        std::vector <double> frequency;
        std::vector<double> LDifarray;
        void InterpFreq(double Linput, double Rleftinput, double Rrightinput, double frequencyInpt,
double LgapEstimate,
std::vector<double>& out);

public:
        freq(std::string parFile);
        void GetCellParameters(double Linput, double Rleftinput, double frequencyInpt, double
LgapEstimate,
std::vector<double>& out);
};*/

template <class PointType>
class Solver
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    std::vector<double> relaxations;
    std::vector<float>  errors;
    std::vector<double> parameters;
    std::vector<double> outTimes;
    int                 recalculate;

    std::shared_ptr<EmissionCurrentSolverPIC<PointType>> emissionCurrentSolverPIC;
    std::shared_ptr<EmissionCurrentSolverPTI<PointType>> emissionCurrentSolverPTI;
    std::shared_ptr<PoissonSolver<PointType>>            fieldSolver;
    std::shared_ptr<MeshGenerator<PointType>>            meshGenerator;
    std::shared_ptr<ParticleGridInterface<PointType>>    particleGridInterface;
    std::shared_ptr<ParticlesMover<PointType>>           particleMover;
    std::shared_ptr<MagneticFieldSolver<PointType>>      magneticFieldSolver;

    int                                                 restartPossible;
    std::vector<unsigned int>                           v;
    std::vector<std::vector<int>>                       blockIndexes;
    std::vector<std::vector<std::vector<unsigned int>>> startCellNumbersGridTmp;
    std::vector<std::vector<std::vector<unsigned int>>> startCellNumbersGridTmp1;
    std::vector<std::vector<std::vector<unsigned int>>> startCellNumbersTmp;
    std::vector<std::vector<std::vector<unsigned int>>> startCellNumbersTmp1;
    std::vector<PointType (*)[9]>                       W;
    std::vector<PointType (*)[9]>                       Wtmp;
    std::vector<particlesFields<PointType>>             fields;
    std::vector<unsigned int>                           saveIndexes;
    std::vector<unsigned int>                           writeIndexes;
    std::vector<std::vector<unsigned int>>              EmptyPlaces;
    std::vector<std::vector<unsigned int>>              EmptyPlaces1;
    std::vector<std::vector<std::vector<unsigned int>>> EmptyPlacesPIC;
    std::vector<std::vector<std::vector<int>>>          flagSave;
    std::vector<PointType>                              alpas;
    volatile int                                        memorySize;
    double                                              outTime;
    int                                                 saveParam;
    double                                              tracesSaveProbability;
    int                                                 numThreads;
    int                                                 blockSize;
    std::vector<int>                                    nonZeros;
    std::vector<double>                                 outTimesLoc;
    std::vector<std::vector<unsigned int>>              startCells;
    std::vector<int>                                    boundaryPoints;
    double                                              maxtime;

    // std::vector<PointType(*)[9]> W;

    //	std::vector<PointType(*)[9]> Wtmp;

    void SetSolverGeneralParameters(const std::vector<std::vector<double>>& par);

    void InitLocalVariables(int nflows, int fieldSize);
    void
    InitOutputData(std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
                   int nflows, int SpaceSize, std::vector<double> masses,
                   std::vector<double> charges, int size);

    std::vector<std::vector<double>> GetParametersPIC();
    std::vector<std::vector<double>> GetParametersPTI();
    void SetParametersPIC(const std::vector<std::vector<double>>& par);
    void SetParametersPTI(const std::vector<std::vector<double>>& par);

    std::vector<std::string> GetVisNames(int solver);
    Solver();

    template <class deviceType>
    void
    SimulateCPUPIC(double& progress, double& progressLoc, std::vector<std::string>& status,
                   bool& flagAbort, const std::shared_ptr<deviceType>& deviceStatus,
                   std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
                   std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData,
                   int flagRestart, std::string& errorMsg);

    // template <class deviceType>
    // void
    // SimulateCPUPTI(double& progress, bool& flagAbort,
    //                const std::shared_ptr<deviceType>& deviceStatus,
    //                std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&
    //                outputData,
    //                std::mutex& plotMutex, const std::shared_ptr<SimulationData>& simulationData);

    template <class deviceType>
    void TimeStepEstimate(std::vector<double>&               result,
                          const std::shared_ptr<deviceType>& deviceStatus, int numThreads,
                          double& progress);

    template <class deviceType>
    void
    doIterartion(bool& flagAbort, const std::vector<unsigned int>& emissionCells, double tolerance,
                 const std::shared_ptr<deviceType>& deviceStatus, std::mutex& plotMutex,
                 const std::shared_ptr<SimulationData>& simulationData,
                 std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData,
                 int flagSaveData, int& iteration, int saveParam, int numberSaveTraces,
                 int numThreads, int blockSize);

    template <class deviceType>
    std::vector<float> estimateError(const std::shared_ptr<deviceType>& deviceStatus, float& maxEr);

    template <class deviceType>
    void ErrorEstimateEvent(const std::shared_ptr<deviceType>& deviceStatus);

    template <class deviceTypeGPU>
    void
    SimulateGPU(double& progress, bool& flagAbort, deviceTypeGPU& deviceStatus,
                std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& outputData){
        /*	solver.SimulateGPU(progress, timeDilatation, deviceStatus, outputData);
        InitFieldSolver();

        std::vector<int> boundaryPoints = fieldSolver->getNearBoundarypoints_all();

        meshGenerator->GetMesh()->Convert2ParticleGridInterface(&particleGridInterface,
        boundaryPoints);



        int size = 10000000;
        int threadsPerBlock = BLOCK_SIZE;
        int blocksPerGrid = size / BLOCK_SIZE + 1;*/

        //	SimulationGPU <DeviceStatusGPU<ParticlesDataTypeGPU, ParticleSourceTypeGPU,
        // EmitterDeviceTypeGPU,
        // GridDataTypeGPU, PointType>, GPUTypes::vector<BoundaryContainerTypeGPU>> <<
        // <blocksPerGrid, threadsPerBlock
        // >>
        //>(deviceStatusDevice, deviceStatus->boundariesDevice); 	start(blocksPerGrid,
        // threadsPerBlock,
        // deviceStatusDevice, deviceStatus->boundariesDevice); 	test(); 		test
        // <<<blocksPerGrid,
        // threadsPerBlock>>>();
    }

    template <class deviceType>
    void InitFieldSolver(const std::shared_ptr<deviceType>& deviceStatus, double& progress)
    {
        //	if (flagInitSolver == 0)
        //	{
        //	deviceStatus->Getmesh()->Convert2GridData(deviceStatus->GetGridData());
        fieldSolver->InitSolver(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                                deviceStatus->Getboundaries(),
                                deviceStatus->GetboundaryConditions(), progress);
        //		flagInitSolver = 1;
        //	}
    }
    template <class deviceType>
    void InitSolvers(const std::shared_ptr<deviceType>& deviceStatus, int flagRestart,
                     double& progressLoc, std::vector<std::string>& status);

    template <class deviceType>
    void CheckConfiguration(const std::shared_ptr<deviceType>& deviceStatus, int flagRestart,
                            std::string& errorMsg);
};

#endif