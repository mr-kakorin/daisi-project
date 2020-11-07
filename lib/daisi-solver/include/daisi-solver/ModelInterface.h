#ifndef PICInterface_H
#define PICInterface_H

#include <memory>
#include <thread>
#include <vector>

class DynamicsData;
class SimulationData;
class lineplot;
class vtkUnstructuredGrid;
class vtkFloatArray;
class vtkPolyData;

template <class T>
class vtkSmartPointer;

namespace Daisi::Emission
{
	enum EnergyDistributionType
	{
		EnergyDistributionType_Bimodal,
		EnergyDistributionType_Maxwell
	};
}

class ModelInterface
{
  public:
    virtual ~ModelInterface() = default;
    virtual std::vector<std::vector<double>> get_particles_positions() = 0;
    virtual std::vector<std::vector<double>> get_particles_moments() = 0;
    virtual std::vector<double> get_particles_charge() = 0;
    virtual std::vector<double> GetElectrodeParametersList(int number) = 0;
    virtual std::vector<std::vector<std::pair<double,double>>> get_mesh() = 0;
    virtual std::vector<double> get_volume_charge() = 0;
    virtual void setElectrodeParametersList(int number, std::vector<double>& input) = 0;
	virtual void setEnergyDistribution( Daisi::Emission::EnergyDistributionType const type ) {};
    virtual void SetglobalFieldConditions(const std::vector<double>& in)             = 0;
    virtual std::vector<double>              GetglobalFieldConditions()              = 0;
    virtual std::vector<std::vector<double>> GetSolverParametersPTI()                = 0;
    virtual std::vector<std::vector<double>> GetSolverParametersPIC()                = 0;
    virtual std::vector<std::vector<double>> GetSolverGeneralParameters()            = 0;
    virtual std::vector<double>              GetFieldSolverParameters()              = 0;
    virtual std::vector<double> GetEmitterInitParameters(int currentFlow)            = 0;
    virtual void SetSolverParametersPTI(const std::vector<std::vector<double>>&)     = 0;
    virtual void SetSolverParametersPIC(const std::vector<std::vector<double>>&)     = 0;
    virtual void SetSolverGeneralParameters(const std::vector<std::vector<double>>&) = 0;

    virtual void                SaveData(std::string)           = 0;
    virtual void                LoadData(std::string)           = 0;
    virtual void                AddSetOfPotentials(std::string) = 0;
    virtual std::vector<double> GetEmittancesList(int flow)     = 0;
    virtual void AddEmittance(int flow, double param) = 0;
    virtual float GetDeviceCurrentTime() = 0;
    virtual void GenerateParticles(int Flow, int flagClear)       = 0;
    virtual void AddFlow(int ParticleType, int DistributionStyle) = 0;
    virtual void AddFlow(int, int, double, double) = 0;
    virtual void SerTest() = 0;
    virtual void AddBoundary(std::string, std::string& errorMsg)   = 0;
    virtual void AddBoundaries(std::string, std::string& errorMsg) = 0;
    virtual void SetAllEmitterParameters(int                                     flow,
                                         const std::vector<std::vector<double>>& params) = 0;
    virtual std::vector<int> GetNumberParticlesFlowsTypes()                              = 0;

    virtual std::vector<double> GetLinacDynSimParameters() = 0;

    virtual std::vector<std::vector<double>> GetAllEmitterParameters(int currentFlow) = 0;

    virtual void GetEmittanceData(std::vector<std::vector<float>>& data, int flowNumber,
                                  int flagNumber, int emFlag) = 0;

    virtual int GetNumberBoundaries() = 0;

    virtual std::vector<double> GetMeshParam()                    = 0;
    virtual void                SetMeshParam(std::vector<double>) = 0;
    virtual void GenerateMesh(std::string meshParam, double&, std::string& errorMsg) = 0;
    virtual std::vector<int> GetMeshBoundariesList()                 = 0;
    virtual void             SetMeshBoundariesList(std::vector<int>) = 0;
    virtual void ChangeGeom(double p1, double p2) = 0;

    virtual void SetConditionPropertiesFromFile(std::string file, int number)       = 0;
    virtual std::vector<int> GetDefaultConditionsList(std::string flag1, int flag2) = 0;
    virtual void SetDefaultConditionsList(std::string flag1, int flag2,
                                          const std::vector<int>&) = 0;
    virtual std::vector<int> GetPropertyConditionsBoundariesList(std::string flag1, int flag2,
                                                                 int) = 0;
    virtual int GetPropertyConditionTypeFlag(std::string flag1, int flag2, int) = 0;
    virtual std::vector<double> GetPropertyConditionManualRestictions(std::string flag1, int flag2,
                                                                      int) = 0;
    virtual void SetPropertyConditionManualRestictions(std::string flag1, int flag2, int,
                                                       std::vector<double> params) = 0;

    virtual void SetPropertyConditionsBoundariesList(std::string flag1, int flag2, int,
                                                     const std::vector<int>&) = 0;
    virtual void AddPropertyCondition(std::string flag1, int flag2, std::string type,
                                      int boundaryTypeFlag) = 0;
    virtual int GetNumberPropertyConditions(std::string flag1, int flag2) = 0;
    virtual std::string GetConditionPropertyType(std::string flag1, int flag2, int i) = 0;

    virtual vtkSmartPointer<vtkUnstructuredGrid> GetBoundaryVTKUnstructuredGrid(int) = 0;
    virtual vtkSmartPointer<vtkUnstructuredGrid> GetMeshBoundaryVTKGrid()            = 0;

    virtual std::vector<std::vector<double>> GetEmitterField(int flowNumber) = 0;

	virtual std::vector<std::vector<float>> GetEmitterFieldFloat(int flowNumber)
	{
		std::vector<std::vector<float>> resultfloat;
		std::vector<std::vector<double>>resultdouble = GetEmitterField( flowNumber );
		resultfloat.resize( resultdouble.size() );
		for (int i=0; i < resultdouble.size(); ++i)
		{
			resultfloat[i].resize( resultdouble[i].size() );
			for (int j=0; j<resultdouble[i].size(); ++j)
				resultfloat[i][j] = static_cast<float>( resultdouble[i][j] );
		}
		return resultfloat;
	}
    virtual vtkUnstructuredGrid* GetVTKGrid()                               = 0;
    virtual vtkUnstructuredGrid* GetVTKGrid(int flag, double param,
                                            vtkSmartPointer<vtkFloatArray>& vtkData,
                                            std::string                     name) = 0;

    virtual int                              GetNumberParticlesFlows() = 0;
    virtual std::vector<double>              GetFlowProperties(int)    = 0;
    virtual std::vector<std::vector<double>> GetDirectionPoints(int)   = 0;             //
    virtual void SetDirectionPoints(int, std::vector<double>, std::vector<double>) = 0; //

    virtual std::vector<int> GetEmitterBoundariesList(int flowNumber) = 0;
    virtual void SetEmitterBoundariesList(int flowNumber, std::vector<int> in,
                                          std::vector<double> params, std::string& errorMsg) = 0;
    virtual std::vector<double> GetAdditionalSourceInf(int flowNumber) = 0;
    virtual void SetAdditionalSourceInf(int flowNumber, std::vector<double> inf) = 0;

    virtual std::vector<std::vector<double>> GetSolverEmissionModelParameters() = 0;

    virtual void SetSolverEmissionModelParameters(std::vector<std::vector<double>> in) = 0;

    virtual double GetEmissionCurrent(int flowNumber) = 0;
    virtual void SetEmissionCurrent(int flowNumber, double current) = 0;

    virtual void SetFieldSolverParameters(std::vector<double> in) = 0;

    virtual void GetParticlesCloud(int flag, int flowNumber,
                                   std::vector<std::vector<void*>>& pointArray,
                                   std::vector<int>& sizeArray, int& sizeElement) = 0;
	virtual void GetMomentumsParticlesCloud(int flag, int flowNumber,
	                               std::vector<std::vector<void*>>& pointArray,
	                               std::vector<int>& sizeArray, int& sizeElement) = 0;
    virtual void GetGridData(void* Array[1], int& size, int& sizeElement, std::string flag,
                             int PlotTypeFlag) = 0;
    virtual void GetGridData(void* Array[1], int& size, int& sizeElement, int flag,
                             int PlotTypeFlag)                                            = 0;
    virtual std::vector<int>                GetBoundariesList()                           = 0;
    virtual std::vector<std::vector<float>> GetCurrentDensityDistribution(int flowNumber) = 0;
    virtual void FieldSimulate(double t, double& progress, double& progressLoc,
                               std::vector<std::string>& status) = 0;
    virtual void SimulateCPUPIC(double& progress, double& progressLoc,
                                std::vector<std::string>& status, bool& flagAbort, int flagRestart,
                                std::string& errorMsg) = 0;
    virtual void SimulateCPUPTI(double& progress, bool& flagAbort) = 0;
    virtual void SimulateGPU(double& progress, bool& flagAbort)    = 0;

    virtual void GenerateGeometry(std::vector<double> in, double& progress, bool& flagAbort,
                                  std::string meshParamFile, std::string projectFolder) = 0;
    virtual void SetAccParams(std::vector<double> in, std::string filename)             = 0;
    virtual std::vector<double> GetAcceleratorParams() = 0;

    virtual std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&
    GetDynamicsData() = 0;

    virtual std::vector<DynamicsData*> GetDynamicsDataLinac() = 0;

    virtual void GetPlotXYBoundary(int boundary, std::vector<float>& Xarr,
                                   std::vector<float>& Yarr) = 0;
    virtual void GetPlotBoundaryRotate(int boundary, std::vector<std::vector<float>>& Xarr,
                                       std::vector<std::vector<float>>& Yarr) = 0;
    virtual void GetPlotXYBoundarySpecial(std::vector<std::vector<float>>& Xarr,
                                          std::vector<std::vector<float>>& Yarr) = 0;
    virtual std::vector<std::shared_ptr<lineplot>> GetPlotsVector()              = 0;
    virtual void addPlot(std::vector<std::string> flag, double x1, double y1, double x2, double y2,
                         int PlotTypeFlag) = 0;
    virtual void addPlot(std::vector<std::string> flag, double x1, double y1, double z1, double x2,
                         double y2, double z2, int PlotTypeFlag) = 0;

    virtual void GetPlot(int plotNumber, std::vector<float>& Xdata,
                         std::vector<std::vector<float>>& Ydata)             = 0;
    virtual std::vector<std::string> GetPlotNames()                          = 0;
    virtual void SetPlotNames(std::vector<std::string> vector)               = 0;
    virtual std::vector<std::shared_ptr<SimulationData>> GetSimulationData() = 0;
    virtual std::vector<std::string> GetVisNames(int i)                      = 0;

    virtual vtkSmartPointer<vtkPolyData> GetVTKBoundaryPoints() = 0;
    virtual void                         InitFieldSolver()      = 0;

    virtual std::vector<std::vector<int>> GetConductorsList() = 0;
    virtual void SetConductorsList(std::vector<int> list, int number, double l,
                                   std::string& errorMsg) = 0;
    virtual void AddConductor()                           = 0;
    virtual void DeleteFlow(int number)                   = 0;

    virtual std::vector<std::vector<float>> GetElectrodeValue(int conductor, int flag) = 0;

    virtual std::vector<double> GetFlowMCNumbers(int flow) = 0;
    virtual void SetFlowMCNumbers(int flow, std::vector<double> numbers) = 0;
    virtual std::vector<double> GetElectrodesCurrents() = 0;
    virtual std::vector<float>  Get_Errors()            = 0;
    virtual void ErrorEstimate(double& progress)        = 0;

    virtual void SetConditionProperties(std::string flag1, int flag2, int i,
                                        std::vector<std::vector<double>> cond) = 0;
    virtual std::vector<double> GetConditionProperties(std::string flag1, int flag2, int i) = 0;
    virtual std::vector<std::string> GetConditionPropertiesNames(std::string flag1, int flag2,
                                                                 int i) = 0;
};
#endif