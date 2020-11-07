#ifndef PICTemplate_H
#define PICTemplate_H

#include "ModelInterface.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <mutex>
template <class PointType>
class Solver;
class deviceStatus;
class BoundaryConditions;

template <class DeviceType, class PointType>
class ModelTemplate : public ModelInterface
{
    friend class boost::serialization::access;

  public:
    virtual ~ModelTemplate() override = default;
	std::vector<std::vector<double>> get_particles_positions() override;
	std::vector<std::vector<double>> get_particles_moments() override;
	std::vector<double> get_particles_charge() override;
	std::vector<std::vector<std::pair<double,double>>> get_mesh() override;
	std::vector<double> get_volume_charge() override;
    int  i;
	virtual void setEnergyDistribution( Daisi::Emission::EnergyDistributionType const type ) override {
		auto dev_em_statuses = deviceStatus->GetEmittersVector();
		for( auto& x: dev_em_statuses )
			x->energy_distribution = type;
	};

    void SetglobalFieldConditions(const std::vector<double>& in);
    std::vector<double>                GetglobalFieldConditions();
    std::vector<std::string>           PlotNames;
    std::shared_ptr<Solver<PointType>> solver;
    std::shared_ptr<DeviceType>        deviceStatus;

    // DeviceStatusType deviceStatus;

    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>> outputData;

    std::vector<std::shared_ptr<lineplot>> lineplots;
    std::mutex                             plotMutex;

    std::vector<std::vector<double>> GetSolverParametersPIC();
    std::vector<std::vector<double>> GetSolverParametersPTI();
    std::vector<std::vector<double>> GetSolverGeneralParameters();
    std::vector<std::string>         flagsPIC;
    std::vector<std::string>         flagsPTI;

    std::vector<int>    GetNumberParticlesFlowsTypes();
    std::vector<double> GetEmitterInitParameters(int currentFlow);
    std::vector<double> GetElectrodeParametersList(int number);
    void setElectrodeParametersList(int number, std::vector<double>& input);
    void SetSolverParametersPIC(const std::vector<std::vector<double>>& par);
    void SetSolverParametersPTI(const std::vector<std::vector<double>>& par);
    void SetSolverGeneralParameters(const std::vector<std::vector<double>>& par);
    std::vector<double>      GetFieldSolverParameters();
    std::vector<std::string> GetConditionPropertiesNames(std::string flag1, int flag2, int i);
    void SetAllEmitterParameters(int flow, const std::vector<std::vector<double>>& params);
    std::vector<std::shared_ptr<SimulationData>> simulationData;
    void GenerateMesh(std::string meshParam, double&, std::string& errorMsg);
    void AddFlow(int ParticleType, int DistributionStyle);
    void AddFlow(int, int, double, double);
    void ChangeGeom(double p1, double p2);
    std::vector<DynamicsData*>       GetDynamicsDataLinac();
    std::vector<std::vector<double>> GetLinacApproximationParameters();
    void SetLinacApproximationParameters(std::vector<std::vector<double>> in);
    std::vector<double> GetLinacDynSimParameters();
    void GetEmittanceData(std::vector<std::vector<float>>& data, int flowNumber, int flagNumber,
                          int emFlag);
    std::vector<double> GetEmittancesList(int flow);
    void AddEmittance(int flow, double param);
    std::vector<double> GetAdditionalSourceInf(int flowNumber);
    void SetAdditionalSourceInf(int flowNumber, std::vector<double> inf);
    void SetConditionPropertiesFromFile(std::string file, int number);
    void GetPlotXYBoundarySpecial(std::vector<std::vector<float>>& Xarr,
                                  std::vector<std::vector<float>>& Yarr);
    std::vector<std::vector<double>> GetEmitterField(int flowNumber);
    void ErrorEstimate(double& progress);
    std::vector<double> GetElectrodesCurrents();
    std::vector<float>  Get_Errors();
    std::vector<double> GetFlowMCNumbers(int flow);
    void SetFlowMCNumbers(int flow, std::vector<double> numbers);
    std::vector<std::vector<float>> GetElectrodeValue(int conductor, int flag);

    void DeleteFlow(int number);
    std::vector<std::vector<int>> GetConductorsList();
    void SetConductorsList(std::vector<int> list, int number, double l, std::string& errorMsg);
    void                         AddConductor();
    void                         InitFieldSolver();
    vtkSmartPointer<vtkPolyData> GetVTKBoundaryPoints();
    std::vector<std::string> GetVisNames(int i);
    std::vector<std::shared_ptr<SimulationData>> GetSimulationData();
    std::vector<std::string>                     GetPlotNames();
    void SetPlotNames(std::vector<std::string> vector);
    void GetGridData(void* Array[1], int& size, int& sizeElement, std::string flag,
                     int PlotTypeFlag);
    void GetGridData(void* Array[1], int& size, int& sizeElement, int flag, int PlotTypeFlag);
    void GetPlot(int plotNumber, std::vector<float>& Xdata, std::vector<std::vector<float>>& Ydata);
    std::vector<std::shared_ptr<lineplot>> GetPlotsVector();
    void addPlot(std::vector<std::string> flag, double x1, double y1, double x2, double y2,
                 int PlotTypeFlag);
    void addPlot(std::vector<std::string> flag, double x1, double y1, double z1, double x2,
                 double y2, double z2, int PlotTypeFlag);
    float GetDeviceCurrentTime();
    void SaveData(std::string fileName);
    void LoadData(std::string fileName);
    void GetPlotXYBoundary(int boundary, std::vector<float>& Xarr, std::vector<float>& Yarr);
    void GetPlotBoundaryRotate(int boundary, std::vector<std::vector<float>>& Xarr,
                               std::vector<std::vector<float>>& Yarr);
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>& GetDynamicsData();
    std::shared_ptr<BoundaryConditions> ParseBoundaryArgument(std::string flag1, int flag2);
    void SimulateGPU(double& progress, bool& flagAbort);
    void SimulateCPUPIC(double& progress, double& progressLoc, std::vector<std::string>& status,
                        bool& flagAbort, int flagRestart, std::string& errorMsg);
    void SimulateCPUPTI(double& progress, bool& flagAbort);
    std::vector<std::vector<double>> GetSolverEmissionModelParameters();
    void SetSolverEmissionModelParameters(std::vector<std::vector<double>> in);
    void SetFieldSolverParameters(std::vector<double> in);
    void FieldSimulate(double t, double& progress, double& progressLoc,
                       std::vector<std::string>& status);
    void SetAccParams(std::vector<double> in, std::string filename);
    void GenerateGeometry(std::vector<double> in, double& progress, bool& flagAbort,
                          std::string meshParamFile, std::string projectFolder);
    std::vector<double> GetAcceleratorParams();
    void GetParticlesCloud(int flag, int flowNumber, std::vector<std::vector<void*>>& pointArray,
                           std::vector<int>& sizeArray, int& sizeElement);
	void GetMomentumsParticlesCloud(int flag, int flowNumber, std::vector<std::vector<void*>>& pointArray,
	                       std::vector<int>& sizeArray, int& sizeElement);
    std::vector<std::vector<float>> GetCurrentDensityDistribution(int flowNumber);
    void GenerateParticles(int flowNumber, int flagClear);
    double GetEmissionCurrent(int flowNumber);
    void SetEmissionCurrent(int flowNumber, double current);
    int                                  GetNumberParticlesFlows();
    vtkSmartPointer<vtkUnstructuredGrid> GetMeshBoundaryVTKGrid();
    vtkUnstructuredGrid*                 GetVTKGrid();
    vtkUnstructuredGrid* GetVTKGrid(int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
                                    std::string name);
    std::vector<double> GetMeshParam();
    void SetMeshParam(std::vector<double> in);
    std::vector<int> GetMeshBoundariesList();
    void SetMeshBoundariesList(std::vector<int> in);
    std::vector<int> GetBoundariesList();
    std::vector<int> GetEmitterBoundariesList(int flowNumber);
    void SetEmitterBoundariesList(int flowNumber, std::vector<int> in, std::vector<double> params,
                                  std::string& errorMsg);
    void SetConditionProperties(std::string flag1, int flag2, int i,
                                std::vector<std::vector<double>> cond);
    std::vector<double> GetConditionProperties(std::string flag1, int flag2, int i);
    void SetDefaultConditionsList(std::string flag1, int flag2, const std::vector<int>& in);
    void SetPropertyConditionsBoundariesList(std::string flag1, int flag2, int i,
                                             const std::vector<int>& in);
    std::vector<int> GetPropertyConditionsBoundariesList(std::string flag1, int flag2, int i);
    int GetPropertyConditionTypeFlag(std::string flag1, int flag2, int i);
    std::vector<double> GetPropertyConditionManualRestictions(std::string flag1, int flag2, int i);
    void SetPropertyConditionManualRestictions(std::string flag1, int flag2, int i,
                                               std::vector<double> params);
    void AddPropertyCondition(std::string flag1, int flag2, std::string type, int boundaryTypeFlag);
    int GetNumberPropertyConditions(std::string flag1, int flag2);
    std::vector<int> GetDefaultConditionsList(std::string flag1, int flag2);
    std::string GetConditionPropertyType(std::string flag1, int flag2, int i);
    void AddBoundary(std::string filename, std::string& error);
    void AddSetOfPotentials(std::string filename);
    void AddBoundaries(std::string filename, std::string& error);
    int                              GetNumberBoundaries();
    std::vector<std::vector<double>> GetAllEmitterParameters(int currentFlow);
    vtkSmartPointer<vtkUnstructuredGrid> GetBoundaryVTKUnstructuredGrid(int number);
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    void SetDirectionPoints(int flowNumber, std::vector<double> sP,
                            std::vector<double> eP);                     /////////////
    std::vector<std::vector<double>> GetDirectionPoints(int flowNumber); ////////////////////
    void SerTest();
    ModelTemplate();
    std::vector<double> GetFlowProperties(int flowNumber);
    void SetBoundary(std::string InputFileName);
    void WriteMesh2VTK(std::string InputFileName);
    void WriteBoundary2VTK(std::string InputFileName);
    void SetFlowData(int FlowNumber);
};

#endif