#include "ModelTemplate.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "DataTypes.h"
#include "ElectrodeCurrent.h"
#include "EmissionCurrentSolver.h"
#include "FlagStringsSolver.h"
#include "Geom.h"
#include "GridData.h"
#include "MagneticFieldSolver.h"
#include "MeshGenerator.h"
#include "Particle.h"
#include "ParticleGridInterface.h"
#include "ParticlesFlow.h"
#include "ParticlesMover.h"
#include "PoissonSolver.h"
#include "Results.h"
#include "Solver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include "MeshContainer2d.h"

template <class DeviceType, class PointType>
std::vector<std::vector<std::pair<double,double>>> ModelTemplate<DeviceType, PointType>::get_mesh()
{
	std::vector<std::vector<DGeo::Point<PointType>>> mesh = deviceStatus->Getmesh()->meshData;
	std::vector<std::vector<std::pair<double,double>>> result;
	result.resize( mesh.size() );
	for ( int i=0; i< mesh.size(); ++i)
	{
		result[i].resize(mesh[i].size());
		for(int j=0;j<mesh[i].size();++j)
			result[i][j] = { mesh[i][j].x, mesh[i][j].y};
	}
	return result;

}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::get_volume_charge()
{

}

template <class DeviceType, class PointType>
std::vector<std::vector<double>> ModelTemplate<DeviceType, PointType>::get_particles_positions()
{
	auto const& vec = deviceStatus->GetFlow( 0 )->GetDynamicsData()->positions;
	std::vector<std::vector<double>> result;
	result.resize( vec.size() );
	for(int i=0; i < vec.size(); ++i)
	{
		result[i].resize(vec[i].size());
		for (int j =0; j< vec[i].size();++j)
			result[i][j] = vec[i][j];
	}
	return result;
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>> ModelTemplate<DeviceType, PointType>::get_particles_moments()
{
	auto const& vec = deviceStatus->GetFlow( 0 )->GetDynamicsData()->momentums;
	std::vector<std::vector<double>> result;
	result.resize( vec.size() );
	for(int i=0; i < vec.size(); ++i)
	{
		result[i].resize(vec[i].size());
		for (int j =0; j< vec[i].size();++j)
			result[i][j] = vec[i][j];
	}
	return result;
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::get_particles_charge()
{
	auto const& vec = deviceStatus->GetFlow( 0 )->GetDynamicsData()->q;
	std::vector<double> result;
	result.resize( vec.size() );
	for(int i=0; i < vec.size(); ++i)
	{
			result[i] = vec[i];
	}
	return result;
}

template <class DeviceType, class PointType>
std::vector<int> ModelTemplate<DeviceType, PointType>::GetNumberParticlesFlowsTypes()
{
    std::vector<int> result;

    for (int j = 0; j < deviceStatus->GetNumberParticlesFlows(); j++)
        result.push_back(deviceStatus->GetFlow(j)->GetDistributionStyle());

    return result;
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetElectrodeParametersList(int number)
{
    return deviceStatus->GetElectrodes()[number]->GetParameteres();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::setElectrodeParametersList(int number,
                                                                      std::vector<double>& input)
{
    deviceStatus->GetElectrodes()[number]->SetParameteres(input);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetglobalFieldConditions(const std::vector<double>& in)
{
    deviceStatus->SetglobalFieldConditions(in);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetglobalFieldConditions()
{
    return deviceStatus->GetglobalFieldConditions();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetSolverParametersPTI(
    const std::vector<std::vector<double>>& par)
{
    solver->SetParametersPTI(par);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetSolverParametersPIC(
    const std::vector<std::vector<double>>& par)
{
    solver->SetParametersPIC(par);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetSolverGeneralParameters(
    const std::vector<std::vector<double>>& par)
{
    solver->SetSolverGeneralParameters(par);
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>> ModelTemplate<DeviceType, PointType>::GetSolverParametersPTI()
{
    return solver->GetParametersPTI();
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>> ModelTemplate<DeviceType, PointType>::GetSolverParametersPIC()
{
    return solver->GetParametersPIC();
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>> ModelTemplate<DeviceType, PointType>::GetSolverGeneralParameters()
{
    std::vector<std::vector<double>> result;

    std::vector<double> result1 = solver->particleGridInterface->GetParameters();

    result.push_back(result1);

    std::vector<double> result2 = solver->particleMover->GetParameters();
    result.push_back(result2);

    result.push_back(solver->parameters);

    return result;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GenerateMesh(std::string meshParam, double& progress,
                                                        std::string& errorMsg)
{
    if (!deviceStatus->GetDomainBoundaryList().size())
    {
        errorMsg = "Domain boundaries list is emplty!";
        return;
    }

    deviceStatus->MergeDomainBoundary();
    solver->meshGenerator->MeshGenerate(meshParam, progress, deviceStatus->GetDomainBoundary(),
                                        deviceStatus->Getmesh(), 0, errorMsg);

    if (errorMsg.size() == 0)
    {
        deviceStatus->Getmesh()->Convert2GridData(deviceStatus->GetGridData());
        deviceStatus->isFieldSimulated  = 0;
        deviceStatus->isFieldSolverInit = 0;
    }
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddFlow(int ParticleType, int DistributionStyle)
{
    int                      n     = GetNumberParticlesFlows();
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + std::to_string(n);
    // simulationData.AddFlags(names);
    deviceStatus->AddFlow(ParticleType, DistributionStyle);

    if (0 == DistributionStyle || 1 == DistributionStyle)
    {
        solver->emissionCurrentSolverPIC->addFlow(n);
        solver->emissionCurrentSolverPTI->addFlow(n);
    }
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddFlow(int ParticleType, int DistributionStyle,
                                                   double massNumber, double chargeNumber)
{
    int                      n     = GetNumberParticlesFlows();
    std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
    for (int i   = 0; i < names.size(); i++)
        names[i] = names[i] + std::to_string(n);
    // simulationData.AddFlags(names);
    deviceStatus->AddFlow(ParticleType, DistributionStyle, massNumber, chargeNumber);

    if (0 == DistributionStyle || 1 == DistributionStyle)
    {
        solver->emissionCurrentSolverPIC->addFlow(n);
        solver->emissionCurrentSolverPTI->addFlow(n);
    }
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetEmitterInitParameters(int currentFlow)
{
    return deviceStatus->GetFlow(currentFlow)->GetEmitterDevice()->GetEmitterInitParameters();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::ChangeGeom(double p1, double p2)
{
    deviceStatus->ChangeGeom(p1, p2);
}

template <class DeviceType, class PointType>
std::vector<DynamicsData*> ModelTemplate<DeviceType, PointType>::GetDynamicsDataLinac()
{
    std::vector<DynamicsData*> t;
    return t;
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>>
ModelTemplate<DeviceType, PointType>::GetLinacApproximationParameters()
{
    std::vector<std::vector<double>> t;
    return t;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetLinacApproximationParameters(
    std::vector<std::vector<double>> in){

}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetLinacDynSimParameters()
{
    std::vector<double> t;
    return t;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetEmittanceData(std::vector<std::vector<float>>& data,
                                                            int flowNumber, int flagNumber,
                                                            int emFlag)
{
    deviceStatus->GetFlow(flowNumber)->GetEmittanceData(data, flagNumber, emFlag);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetEmittancesList(int flow)
{
    return deviceStatus->GetFlow(flow)->GetEmittancesList();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddEmittance(int flow, double param)
{
    deviceStatus->GetFlow(flow)->AddEmittance(param);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetAdditionalSourceInf(int flowNumber)
{
    return deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->GetAdditionalSourceInf();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetAdditionalSourceInf(int flowNumber,
                                                                  std::vector<double> inf)
{
    deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->SetAdditionalSourceInf(inf);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetConditionPropertiesFromFile(std::string file,
                                                                          int number)
{
    deviceStatus->GetboundaryConditions()->SetConditionPropertiesFromFile(number, file);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetPlotXYBoundarySpecial(
    std::vector<std::vector<float>>& Xarr, std::vector<std::vector<float>>& Yarr)
{
    for (int i = 0; i < deviceStatus->GetboundaryConditions()->PropertyConditionListSize(); i++)
    {
        for (int j = 0;
             j <
             deviceStatus->GetboundaryConditions()->GetPropertyConditionsBoundariesList(i).size();
             j++)
        {
            int k =
                deviceStatus->GetboundaryConditions()->GetPropertyConditionsBoundariesList(i)[j];
            deviceStatus->GetboundariesForFlows()[k]->GetPlotRotate1(Xarr, Yarr, 108);
        }
    }
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>>
ModelTemplate<DeviceType, PointType>::GetEmitterField(int flowNumber)
{
    return deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->GetEmitterField();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::ErrorEstimate(double& progress)
{
    solver->ErrorEstimateEvent(deviceStatus);
    progress = 1;
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetElectrodesCurrents()
{
    std::vector<double> result;
    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
        result.push_back(deviceStatus->GetElectrodes()[i]->GetCurrent());

    return result;
}

template <class DeviceType, class PointType>
std::vector<float> ModelTemplate<DeviceType, PointType>::Get_Errors()
{
    return solver->errors;
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetFlowMCNumbers(int flow)
{
    return deviceStatus->GetFlow(flow)->GetFlowMCNumbers();
}
template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetFlowMCNumbers(int flow, std::vector<double> numbers)
{
    deviceStatus->GetFlow(flow)->SetFlowMCNumbers(numbers);
}

template <class DeviceType, class PointType>
std::vector<std::vector<float>>
ModelTemplate<DeviceType, PointType>::GetElectrodeValue(int conductor, int flag)
{
    return deviceStatus->GetElectrodes()[conductor]->GetElectrodeValueF(flag);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::DeleteFlow(int number)
{
    deviceStatus->DeleteFlow(number);
}

template <class DeviceType, class PointType>
std::vector<std::vector<int>> ModelTemplate<DeviceType, PointType>::GetConductorsList()
{
    return deviceStatus->GetConductorsList();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetConductorsList(std::vector<int> list, int number,
                                                             double l, std::string& errorMsg)
{
    deviceStatus->SetConductorsList(list, number, l, errorMsg);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddConductor()
{
    deviceStatus->AddConductor();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::InitFieldSolver()
{
    double progressLoc;
    solver->InitFieldSolver(deviceStatus, progressLoc);
    deviceStatus->isFieldSolverInit = 1;
}

template <class DeviceType, class PointType>
vtkSmartPointer<vtkPolyData> ModelTemplate<DeviceType, PointType>::GetVTKBoundaryPoints()
{
    return deviceStatus->Getmesh()->GetVTKBoundaryPoints();
}

template <class DeviceType, class PointType>
std::vector<std::string> ModelTemplate<DeviceType, PointType>::GetVisNames(int i)
{
    int nflow       = deviceStatus->GetNumberParticlesFlows();
    int nElectrodes = deviceStatus->GetElectrodes().size();

    std::vector<std::string> result;

    for (int j = 0; j < nflow; j++)
    {
        if (deviceStatus->GetFlow(j)->GetDistributionStyle() < 5)
        {
            for (int i = 0; i < flagStringsSolver::simulationDataNamesFlowPIC.size(); i++)
                result.push_back(flagStringsSolver::simulationDataNamesFlowPIC[i] + " of flow " +
                                 std::to_string(j));
        }
        else
        {
            for (int i = 0; i < flagStringsSolver::simulationDataNamesFlowPICAccel.size(); i++)
                result.push_back(flagStringsSolver::simulationDataNamesFlowPICAccel[i] +
                                 " of flow " + std::to_string(j));
        }
    }

    for (int i = 0; i < flagStringsSolver::simulationDataNamesElectrode.size(); i++)
    {
        for (int j = 0; j < nElectrodes; j++)
            result.push_back(flagStringsSolver::simulationDataNamesElectrode[i] + " of electrode " +
                             std::to_string(j));
    }

    return result;
}

template <class DeviceType, class PointType>
std::vector<std::shared_ptr<SimulationData>>
ModelTemplate<DeviceType, PointType>::GetSimulationData()
{
    return simulationData;
}

template <class DeviceType, class PointType>
std::vector<std::string> ModelTemplate<DeviceType, PointType>::GetPlotNames()
{
    return PlotNames;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetPlotNames(std::vector<std::string> vector)
{
    PlotNames = vector;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetGridData(void* Array[1], int& size, int& sizeElement,
                                                       std::string flag, int PlotTypeFlag)
{
    deviceStatus->GetGridData()->GetData(Array, size, sizeElement, flag, PlotTypeFlag);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetGridData(void* Array[1], int& size, int& sizeElement,
                                                       int flag, int PlotTypeFlag)
{
    deviceStatus->GetGridData()->GetDataIntFlag(Array, size, sizeElement, flag, PlotTypeFlag);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetPlot(int plotNumber, std::vector<float>& Xdata,
                                                   std::vector<std::vector<float>>& Ydata)
{
    int                             npoints = 1000;
    std::vector<DGeo::Edge<double>> resArr  = lineplots[plotNumber]->line.resize(npoints);

    float length = 0;
    for (int i = 0; i < npoints; i++)
    {
        Xdata.push_back(length);
        length = length + resArr[i].length();
    }
    Xdata.push_back(length);

    plotMutex.lock();
    for (int j = 0; j < lineplots[plotNumber]->flag.size(); j++)
    {
        std::vector<float> Ydatatmp;
        for (int i = 0; i < npoints; i++)
        {
            Ydatatmp.push_back(deviceStatus->GetGridData()->interpolatePoint(
                resArr[i].point1.x, resArr[i].point1.y, resArr[i].point1.z,
                lineplots[plotNumber]->flag[j], lineplots[plotNumber]->PlotTypeFlag));
        }
        Ydatatmp.push_back(deviceStatus->GetGridData()->interpolatePoint(
            resArr[npoints - 1].point2.x, resArr[npoints - 1].point2.y,
            resArr[npoints - 1].point2.z, lineplots[plotNumber]->flag[j],
            lineplots[plotNumber]->PlotTypeFlag));
        Ydata.push_back(Ydatatmp);
    }
    plotMutex.unlock();
}

template <class DeviceType, class PointType>
std::vector<std::shared_ptr<lineplot>> ModelTemplate<DeviceType, PointType>::GetPlotsVector()
{
    return lineplots;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::addPlot(std::vector<std::string> flag, double x1,
                                                   double y1, double x2, double y2,
                                                   int PlotTypeFlag)
{
    lineplots.push_back(
        std::shared_ptr<lineplot>(new lineplot(flag, x1, y1, x2, y2, PlotTypeFlag)));
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::addPlot(std::vector<std::string> flag, double x1,
                                                   double y1, double z1, double x2, double y2,
                                                   double z2, int PlotTypeFlag)
{
    lineplots.push_back(
        std::shared_ptr<lineplot>(new lineplot(flag, x1, y1, z1, x2, y2, z2, PlotTypeFlag)));
}

template <class DeviceType, class PointType>
float ModelTemplate<DeviceType, PointType>::GetDeviceCurrentTime()
{
    if (deviceStatus->GetNumberParticlesFlows() != 0)
        return deviceStatus->GetFlow(0)->GetDynamicsData()->Time / (LIGHT_VELOCITY());
    return 0;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SaveData(std::string fileName)
{
    std::ofstream                   ofs(fileName.c_str(), std::ios::out | std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    oa << outputData;
    oa << simulationData;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::LoadData(std::string fileName)
{
    std::ifstream                   ifs(fileName.c_str(), std::ios::in | std::ios::binary);
    boost::archive::binary_iarchive ia(ifs);
    ia >> outputData;
    ia >> simulationData;
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetPlotXYBoundary(int boundary, std::vector<float>& Xarr,
                                                             std::vector<float>& Yarr)
{
    deviceStatus->Getboundaries()[boundary]->GetPlotXY(Xarr, Yarr);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetPlotBoundaryRotate(
    int boundary, std::vector<std::vector<float>>& Xarr, std::vector<std::vector<float>>& Yarr)
{
    deviceStatus->Getboundaries()[boundary]->GetPlotRotate(Xarr, Yarr, 100);
}

template <class DeviceType, class PointType>
std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&
ModelTemplate<DeviceType, PointType>::GetDynamicsData()
{
    return outputData;
}

template <class DeviceType, class PointType>
std::shared_ptr<BoundaryConditions>
ModelTemplate<DeviceType, PointType>::ParseBoundaryArgument(std::string flag1, int flag2)
{
    if (flag1 == flagStringsSolver::poisson)
        return deviceStatus->GetboundaryConditions();
    if (flag1 == flagStringsSolver::flowBoundaryList)
        return deviceStatus->GetFlow(flag2)->GetboundaryConditions();
    return deviceStatus->GetboundaryConditions();
}

// void InCell(ParticleGridInterfaceType ParticleGridInterface, )
template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SimulateGPU(double& progress, bool& flagAbort){
    // DeviceTypeGPU* deviceStatusHost = new DeviceTypeGPU(deviceStatus);
    // DeviceTypeGPU* deviceStatusDevice = deviceStatusHost->AllocateOnGPU();
    //	solver.SimulateGPU(progress, timeDilatation, deviceStatusDevice, outputData);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SimulateCPUPIC(double& progress, double& progressLoc,
                                                          std::vector<std::string>& status,
                                                          bool& flagAbort, int flagRestart,
                                                          std::string& errorMsg)
{
    int              nflow       = deviceStatus->GetNumberParticlesFlows();
    int              nElectrodes = deviceStatus->GetElectrodes().size();
    std::vector<int> flowsStyles;

    for (int j = 0; j < nflow; j++)
        flowsStyles.push_back(deviceStatus->GetFlow(j)->GetDistributionStyle());

    simulationData.push_back(
        std::shared_ptr<SimulationData>(new SimulationData(flowsStyles, nElectrodes, 0)));

    solver->SimulateCPUPIC(progress, progressLoc, status, flagAbort, deviceStatus, outputData,
                           plotMutex, simulationData.back(), flagRestart, errorMsg);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SimulateCPUPTI(double& progress, bool& flagAbort){
    // simulationData.push_back(std::shared_ptr<SimulationData>(new SimulationData()));
    // solver->SimulateCPUPTI(progress, flagAbort, deviceStatus, outputData, plotMutex,
    //                        simulationData.back());
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>>
ModelTemplate<DeviceType, PointType>::GetSolverEmissionModelParameters()
{
    return solver->emissionCurrentSolverPIC->GetParameters();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetSolverEmissionModelParameters(
    std::vector<std::vector<double>> in)
{
    solver->emissionCurrentSolverPIC->SetParameters(in);
    solver->emissionCurrentSolverPTI->SetParameters(in);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetFieldSolverParameters(std::vector<double> in)
{
    solver->fieldSolver->SetParameters(in);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetFieldSolverParameters()
{
    return solver->fieldSolver->GetParameters();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::FieldSimulate(double t, double& progress,
                                                         double&                   progressLoc,
                                                         std::vector<std::string>& status)
{
    // deviceStatus->GetGridData()->densityReset();
    status.push_back(".......................................");
    status.push_back("External field simulations process started");
    progressLoc = 0.0;
    progress    = 0.0;

    if (!deviceStatus->isFieldSolverInit)
    {
        status.push_back("Field solver initialization");
        solver->InitFieldSolver(deviceStatus, progressLoc);
        progressLoc = 1.0;
    }
    progress = 0.5;
    status.push_back("Solving linear system");

    progressLoc = 0.0;

    solver->fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                                       deviceStatus->GetboundaryConditions(), t, progressLoc);
    deviceStatus->GetGridData()->ApplyTimeDepending(deviceStatus->GetglobalFieldConditions(), t);
    deviceStatus->isFieldSimulated  = 1;
    deviceStatus->isFieldSolverInit = 1;
    progress                        = 1.0;
    status.push_back("Done");
    status.push_back(".......................................");
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetAccParams(std::vector<double> in,
                                                        std::string filename)
{
    //	deviceStatus->SetAcceleratorParams(in, filename);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GenerateGeometry(std::vector<double> in,
                                                            double& progress, bool& flagAbort,
                                                            std::string meshParamFile,
                                                            std::string projectFolder){
    // deviceStatus->SetAcceleratorParams(in);
    //	solver->GenerateGeometry(progress, flagAbort, deviceStatus, outputData, plotMutex,
    // simulationData,
    // meshParamFile, projectFolder);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetAcceleratorParams()
{
    return {};
    //	return deviceStatus->GetAcceleratorParams();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetParticlesCloud(
    int flag, int flowNumber, std::vector<std::vector<void*>>& pointArray,
    std::vector<int>& sizeArray, int& sizeElement)
{
    int n = deviceStatus->GetFlow(flowNumber)->GetNumberOfThreads();

    if (n == 0)
    {
        pointArray.resize(1);
        sizeArray.resize(1);
        deviceStatus->GetFlow(flowNumber)
            ->GetDynamicsData()
            ->GetParticlesCloud(flag, pointArray[0], sizeArray[0], sizeElement);
        return;
    }
    pointArray.resize(n);
    sizeArray.resize(n);
    for (int i = 0; i < n; i++)
        deviceStatus->GetFlow(flowNumber)
            ->GetDynamicsData(i)
            ->GetParticlesCloud(flag, pointArray[i], sizeArray[i], sizeElement);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GetMomentumsParticlesCloud(
		int flag, int flowNumber, std::vector<std::vector<void*>>& pointArray,
		std::vector<int>& sizeArray, int& sizeElement)
{
	int n = deviceStatus->GetFlow(flowNumber)->GetNumberOfThreads();

	if (n == 0)
	{
		pointArray.resize(1);
		sizeArray.resize(1);
		deviceStatus->GetFlow(flowNumber)
				->GetDynamicsData()
				->GetMomentumsParticlesCloud(flag, pointArray[0], sizeArray[0], sizeElement);
		return;
	}
	pointArray.resize(n);
	sizeArray.resize(n);
	for (int i = 0; i < n; i++)
		deviceStatus->GetFlow(flowNumber)
				->GetDynamicsData(i)
				->GetMomentumsParticlesCloud(flag, pointArray[i], sizeArray[i], sizeElement);
}

template <class DeviceType, class PointType>
std::vector<std::vector<float>>
ModelTemplate<DeviceType, PointType>::GetCurrentDensityDistribution(int flowNumber)
{
    return deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->GetCurrentDensityDistribution();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::GenerateParticles(int flowNumber, int flagClear)
{
    std::vector<unsigned int> t;
    //	deviceStatus->GetGridData()->densityReset();
    deviceStatus->GetFlow(flowNumber)
        ->GenerateParticles(t, flagClear, 0, 0, deviceStatus->GetGridData(), 0, 0);
}

template <class DeviceType, class PointType>
double ModelTemplate<DeviceType, PointType>::GetEmissionCurrent(int flowNumber)
{
    return deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->GetEmissionCurrent();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetEmissionCurrent(int flowNumber, double current)
{
    solver->emissionCurrentSolverPIC->SetEmissionCurrent(
        deviceStatus->GetFlow(flowNumber)->GetEmitterDevice(), current);
}

template <class DeviceType, class PointType>
int ModelTemplate<DeviceType, PointType>::GetNumberParticlesFlows()
{
    return deviceStatus->GetNumberParticlesFlows();
}

template <class DeviceType, class PointType>
vtkSmartPointer<vtkUnstructuredGrid> ModelTemplate<DeviceType, PointType>::GetMeshBoundaryVTKGrid()
{
    return deviceStatus->GetDomainBoundaryVTKUnstructuredGrid();
}

template <class DeviceType, class PointType>
vtkUnstructuredGrid* ModelTemplate<DeviceType, PointType>::GetVTKGrid()
{
    return deviceStatus->Getmesh()->GetVTKGrid();
}

template <class DeviceType, class PointType>
vtkUnstructuredGrid* ModelTemplate<DeviceType, PointType>::GetVTKGrid(
    int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData, std::string name)
{
    return deviceStatus->Getmesh()->GetVTKGrid(flag, param, vtkData, deviceStatus->GetGridData(),
                                               name);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetMeshParam()
{
    return solver->meshGenerator->GetMeshParam();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetMeshParam(std::vector<double> in){
    //		solver.meshGenerator->SetMeshParam(in);
}

template <class DeviceType, class PointType>
std::vector<int> ModelTemplate<DeviceType, PointType>::GetMeshBoundariesList()
{
    return deviceStatus->GetDomainBoundaryList();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetMeshBoundariesList(std::vector<int> in)
{
    deviceStatus->SetDomainBoundaryList(in);
    /*SetMeshBoundariesList( in);
    std::vector<BoundaryContainerType>  bound;
    for (int i = 0; i < in.size(); i++)
    {
    bound.push_back(deviceStatus->Getboundaries()[in[i]]);
    };
    solver.meshGenerator->SetBoundaryList(in, bound);*/
}

template <class DeviceType, class PointType>
std::vector<int> ModelTemplate<DeviceType, PointType>::GetBoundariesList()
{
    return deviceStatus->GetBoundariesList();
}

template <class DeviceType, class PointType>
std::vector<int> ModelTemplate<DeviceType, PointType>::GetEmitterBoundariesList(int flowNumber)
{
    return deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->GetBoundariesList();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetEmitterBoundariesList(int flowNumber,
                                                                    std::vector<int>    in,
                                                                    std::vector<double> params,
                                                                    std::string&        errorMsg)
{
    deviceStatus->SetEmitterBoundariesList(flowNumber, in, params, errorMsg);
    solver->emissionCurrentSolverPIC->SetEmissionCurrent(
        deviceStatus->GetFlow(flowNumber)->GetEmitterDevice(), 1.0);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetConditionProperties(
    std::string flag1, int flag2, int i, std::vector<std::vector<double>> cond)
{

    int TypeFlag = GetPropertyConditionTypeFlag(flag1, flag2, i);

    if (TypeFlag == 0)
    {
        if (flag1 == flagStringsSolver::poisson)
        {
            if (ParseBoundaryArgument(flag1, flag2)->GetConditionPropertiesSimple(i) != cond[0])
            {
                deviceStatus->isFieldSimulated  = 0;
                deviceStatus->isFieldSolverInit = 0;
            }
        }
        ParseBoundaryArgument(flag1, flag2)->SetConditionProperties(i, cond[0]);
    }
    else
    {
        ParseBoundaryArgument(flag1, flag2)->SetPropertyConditionManualRestictions(i, cond[0]);
        ParseBoundaryArgument(flag1, flag2)->SetConditionProperties(i, cond[1]);
    }
}
template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetConditionProperties(std::string flag1,
                                                                                 int flag2, int i)
{
    return ParseBoundaryArgument(flag1, flag2)->GetConditionPropertiesSimple(i);
}

template <class DeviceType, class PointType>
std::vector<std::string>
ModelTemplate<DeviceType, PointType>::GetConditionPropertiesNames(std::string flag1, int flag2,
                                                                  int i)
{
    return ParseBoundaryArgument(flag1, flag2)->GetConditionPropertiesSimpleNames(i);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetDefaultConditionsList(std::string flag1, int flag2,
                                                                    const std::vector<int>& in)
{
    ParseBoundaryArgument(flag1, flag2)->SetDefaultConditionsList(in);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetPropertyConditionsBoundariesList(
    std::string flag1, int flag2, int i, const std::vector<int>& in)
{
    if (flag1 == flagStringsSolver::poisson)
    {
        if (ParseBoundaryArgument(flag1, flag2)->GetPropertyConditionsBoundariesList(i) != in)
        {
            deviceStatus->isFieldSimulated  = 0;
            deviceStatus->isFieldSolverInit = 0;
        }
    }
    ParseBoundaryArgument(flag1, flag2)->SetPropertyConditionsBoundariesList(i, in);
}

template <class DeviceType, class PointType>
std::vector<int>
ModelTemplate<DeviceType, PointType>::GetPropertyConditionsBoundariesList(std::string flag1,
                                                                          int flag2, int i)
{
    return ParseBoundaryArgument(flag1, flag2)->GetPropertyConditionsBoundariesList(i);
}

template <class DeviceType, class PointType>
int ModelTemplate<DeviceType, PointType>::GetPropertyConditionTypeFlag(std::string flag1, int flag2,
                                                                       int i)
{
    return ParseBoundaryArgument(flag1, flag2)->GetPropertyConditionTypeFlag(i);
}

template <class DeviceType, class PointType>
std::vector<double>
ModelTemplate<DeviceType, PointType>::GetPropertyConditionManualRestictions(std::string flag1,
                                                                            int flag2, int i)
{
    return ParseBoundaryArgument(flag1, flag2)->GetPropertyConditionManualRestictions(i);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetPropertyConditionManualRestictions(
    std::string flag1, int flag2, int i, std::vector<double> params)
{
    ParseBoundaryArgument(flag1, flag2)->SetPropertyConditionManualRestictions(i, params);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddPropertyCondition(std::string flag1, int flag2,
                                                                std::string type,
                                                                int         boundaryTypeFlag)
{
    ParseBoundaryArgument(flag1, flag2)->AddPropertyCondition(type, boundaryTypeFlag);
}

template <class DeviceType, class PointType>
int ModelTemplate<DeviceType, PointType>::GetNumberPropertyConditions(std::string flag1, int flag2)
{
    int k = ParseBoundaryArgument(flag1, flag2)->GetNumberProperties();
    return k;
}

template <class DeviceType, class PointType>
std::vector<int> ModelTemplate<DeviceType, PointType>::GetDefaultConditionsList(std::string flag1,
                                                                                int flag2)
{
    return ParseBoundaryArgument(flag1, flag2)->GetDefaultConditionsList();
}

template <class DeviceType, class PointType>
std::string ModelTemplate<DeviceType, PointType>::GetConditionPropertyType(std::string flag1,
                                                                           int flag2, int i)
{
    return ParseBoundaryArgument(flag1, flag2)->GetConditionPropertyType(i);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddBoundary(std::string filename, std::string& errorMsg)
{
    deviceStatus->AddBoundary(filename, errorMsg);
    deviceStatus->GetboundaryConditions()->AddDefaultConditionsList(
        int(deviceStatus->Getboundaries().size()) - 1);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddSetOfPotentials(std::string filename)
{
    deviceStatus->AddSetOfPotentials(filename);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::AddBoundaries(std::string filename,
                                                         std::string& errorMsg)
{
    int i1 = deviceStatus->Getboundaries().size();

    deviceStatus->AddBoundaries(filename, errorMsg);

    int i2 = deviceStatus->Getboundaries().size();

    for (int i = i1; i < i2; i++)
        deviceStatus->GetboundaryConditions()->AddDefaultConditionsList(i);
}

template <class DeviceType, class PointType>
int ModelTemplate<DeviceType, PointType>::GetNumberBoundaries()
{
    return int(deviceStatus->Getboundaries().size());
}

template <class DeviceType, class PointType>
vtkSmartPointer<vtkUnstructuredGrid>
ModelTemplate<DeviceType, PointType>::GetBoundaryVTKUnstructuredGrid(int number)
{
    return deviceStatus->Getboundaries()[number]->GetBoundaryVTKUnstructuredGrid();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetAllEmitterParameters(
    int flow, const std::vector<std::vector<double>>& params)
{
    deviceStatus->GetFlow(flow)->GetEmitterDevice()->SetAllParameters(params);
    /*std::vector<int> numbers;
    numbers.push_back(params[0][0]);
    for(int i=0;i<params[2].size();i++)
            numbers.push_back(params[2][i]);

    SetParticlesNumber(flow, numbers);

    SetEmissionCurrent(flow, params[1][0]);
    SetDistributionParameters(flow, params[3]);*/
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>>
ModelTemplate<DeviceType, PointType>::GetAllEmitterParameters(int currentFlow)
{
    return deviceStatus->GetFlow(currentFlow)->GetEmitterDevice()->GetAllParameters();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetDirectionPoints(int flowNumber,
                                                              std::vector<double> sP,
                                                              std::vector<double> eP) /////////////
{
    deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->SetDirectionPoints(sP, eP);
}

template <class DeviceType, class PointType>
std::vector<std::vector<double>>
ModelTemplate<DeviceType, PointType>::GetDirectionPoints(int flowNumber) ////////////////////
{
    return deviceStatus->GetFlow(flowNumber)->GetEmitterDevice()->GetDirectionPoints();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SerTest()
{
    i++;
}

template <class DeviceType, class PointType>
ModelTemplate<DeviceType, PointType>::ModelTemplate()
{

    deviceStatus = std::shared_ptr<DeviceType>(new DeviceType());

    solver = std::shared_ptr<Solver<PointType>>(new Solver<PointType>());

    //	solver = SolverType();
    /*std::vector<std::string> PlotNames;
    DeviceType deviceStatus;
    SolverType solver;
    //DeviceStatusType deviceStatus;

    std::vector<std::vector<DynamicsData*>>  outputData;
    std::vector<lineplot> lineplots;
    std::mutex plotMutex;


    SimulationData simulationData;*/

    //		deviceStatus = DeviceStatus<ParticlesDataType, EmitterDeviceType, GridDataType,
    // PointType>(); 		solver.meshGenerator = MeshGeneratorType();
    //		particleGridInterface = ParticleGridInterfaceType();
    // devicestatus = &DeviceStatus(1);
}

template <class DeviceType, class PointType>
std::vector<double> ModelTemplate<DeviceType, PointType>::GetFlowProperties(int flowNumber)
{
    return deviceStatus->GetFlow(flowNumber)->GetFlowProperties();
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetBoundary(std::string InputFileName){
    //	solver->meshGenerator = MeshGeneratorType(InputFileName);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::WriteMesh2VTK(std::string InputFileName)
{
    solver->meshGenerator->WriteMesh2VTK(InputFileName);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::WriteBoundary2VTK(std::string InputFileName)
{
    solver->meshGenerator->WriteBoundary2VTK(InputFileName);
}

template <class DeviceType, class PointType>
void ModelTemplate<DeviceType, PointType>::SetFlowData(int FlowNumber){

}

template <class DeviceType, class PointType>
template <class Archive>
void ModelTemplate<DeviceType, PointType>::save(Archive& ar, const unsigned int) const
{
    ar& solver;
    ar& deviceStatus;
    ar& lineplots;
    ar& PlotNames;
    // ar & simulationData;
}
template <class DeviceType, class PointType>
template <class Archive>
void ModelTemplate<DeviceType, PointType>::load(Archive& ar, const unsigned int)
{
    ar& solver;
    ar& deviceStatus;
    ar& lineplots;
    ar& PlotNames;
    // ar & simulationData;
    // deviceStatus->Getmesh()->Convert2GridData(deviceStatus->GetGridData());
}

template void ModelTemplate<device2daxsfloat, float>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void ModelTemplate<device2daxsdouble, double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device2daxsfloat, float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device2daxsdouble, double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);

template void ModelTemplate<device2dfloat, float>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void ModelTemplate<device2ddouble, double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device2dfloat, float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device2ddouble, double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);

template void ModelTemplate<device2dpolarfloat, float>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void ModelTemplate<device2dpolardouble, double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device2dpolarfloat, float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device2dpolardouble, double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);

template void ModelTemplate<device3dExtrfloat, float>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void ModelTemplate<device3dExtrdouble, double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device3dExtrfloat, float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;
template void ModelTemplate<device3dExtrdouble, double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);

template class ModelTemplate<device2daxsfloat, float>;
template class ModelTemplate<device2daxsdouble, double>;
template class ModelTemplate<device2dfloat, float>;
template class ModelTemplate<device2ddouble, double>;
template class ModelTemplate<device2dpolarfloat, float>;
template class ModelTemplate<device2dpolardouble, double>;
template class ModelTemplate<device3dExtrfloat, float>;
template class ModelTemplate<device3dExtrdouble, double>;
