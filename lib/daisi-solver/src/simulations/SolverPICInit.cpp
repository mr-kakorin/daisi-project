#include "DataTypes.h"
#include "ElectrodeCurrent.h"
#include "EmissionCurrentSolver.h"
#include "Particle.h"
#include "ParticleGridInterface.h"
#include "ParticlesFlow.h"
#include "ParticlesMover.h"
#include "PoissonSolver.h"
#include "Results.h"
#include "Solver.h"
#include <Constants.h>
#include "BoundaryConditions.h"

template
class Solver<float>;

template
class Solver<double>;

template<class PointType>
void Solver<PointType>::InitLocalVariables(int nflows, int fieldSize) {
  outTime = parameters[1];
  saveParam = parameters[0];
  tracesSaveProbability = parameters[1];
  numThreads = parameters[2];
  blockSize = int(parameters[3]);
  EmptyPlaces.resize(numThreads);
  outTime = outTime * (1e-9) * LIGHT_VELOCITY();
  alpas.resize(nflows);
  memorySize = 10000000 / numThreads;
  nonZeros.reserve(1000);
  EmptyPlacesPIC.clear();
  EmptyPlacesPIC.resize(nflows);
  flagSave.clear();
  flagSave.resize(nflows);
  outTimesLoc.resize(outTimes.size());

  W.resize(numThreads);
  fields.resize(numThreads);
  blockIndexes.resize(numThreads);

  recalculate = fieldSolver->RecalculateParameter;

  for (int i = 0; i < numThreads; i++) {
    W[i] = new PointType[blockSize * 2][9];
    fields[i].Init(nflows, fieldSize, blockSize * 2);
  };

  for (int i = 0; i < nflows; i++) {
    EmptyPlacesPIC[i].resize(numThreads);
    flagSave[i].resize(numThreads);
    for (int thread = 0; thread < numThreads; thread++) {
      flagSave[i][thread].resize(outTimes.size());
      for (int j = 0; j < outTimes.size(); j++)
        flagSave[i][thread][j] = 0;

      EmptyPlacesPIC[i][thread].reserve(blockSize * 2);
    }
  }

  for (int thread = 0; thread < numThreads; thread++)
    blockIndexes[thread].resize(2 * blockSize);
};

template<class PointType>
void Solver<PointType>::InitOutputData(
    std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>> &outputData, int nflows,
    int SpaceSize, std::vector<double> masses, std::vector<double> charges, int size) {
  outputData.resize(outputData.size() + 1);
  outputData.back().resize(outTimes.size());
  for (int j = 0; j < outTimes.size(); j++) {
    outTimesLoc[j] = outTimes[j] * (1e-9) * LIGHT_VELOCITY();
    outputData.back()[j].resize(nflows);
    for (int i = 0; i < nflows; i++) {
      outputData.back()[j][i] = std::shared_ptr<DynamicsData>(new DynamicsData());
      outputData.back()[j][i]->Init(SpaceSize, masses[i], charges[i], size, numThreads, 0);
    }
  }
};

template void Solver<double>::InitSolvers<device2daxsdouble>(
    const std::shared_ptr<device2daxsdouble> &deviceStatus, int flagRestart, double &progressLoc,
    std::vector<std::string> &status);

template void
Solver<float>::InitSolvers<device2daxsfloat>(const std::shared_ptr<device2daxsfloat> &deviceStatus,
                                             int flagRestart, double &progressLoc,
                                             std::vector<std::string> &status);

template void
Solver<double>::InitSolvers<device2ddouble>(const std::shared_ptr<device2ddouble> &deviceStatus,
                                            int flagRestart, double &progressLoc,
                                            std::vector<std::string> &status);

template void
Solver<float>::InitSolvers<device2dfloat>(const std::shared_ptr<device2dfloat> &deviceStatus,
                                          int flagRestart, double &progressLoc,
                                          std::vector<std::string> &status);

template void Solver<double>::InitSolvers<device2dpolardouble>(
    const std::shared_ptr<device2dpolardouble> &deviceStatus, int flagRestart, double &progressLoc,
    std::vector<std::string> &status);

template void Solver<float>::InitSolvers<device2dpolarfloat>(
    const std::shared_ptr<device2dpolarfloat> &deviceStatus, int flagRestart, double &progressLoc,
    std::vector<std::string> &status);

template void
Solver<double>::InitSolvers<device3ddouble>(const std::shared_ptr<device3ddouble> &deviceStatus,
                                            int flagRestart, double &progressLoc,
                                            std::vector<std::string> &status);

template void Solver<float>::InitSolvers<device3dExtrfloat>(
    const std::shared_ptr<device3dExtrfloat> &deviceStatus, int flagRestart, double &progressLoc,
    std::vector<std::string> &status);

template void Solver<double>::InitSolvers<device3dExtrdouble>(
    const std::shared_ptr<device3dExtrdouble> &deviceStatus, int flagRestart, double &progressLoc,
    std::vector<std::string> &status);

template<class PointType>
template<class deviceType>
void Solver<PointType>::InitSolvers(const std::shared_ptr<deviceType> &deviceStatus,
                                    int flagRestart, double &progressLoc,
                                    std::vector<std::string> &status) {
  particleMover->InitParallel(numThreads);
  particleMover->flagInit = 0;
  progressLoc = 0;
  if (!deviceStatus->isFieldSolverInit) {
    status.push_back("Field solver initialization");
    InitFieldSolver(deviceStatus, progressLoc);
    deviceStatus->isFieldSolverInit = 1;
    progressLoc = 1;
  }
  progressLoc = 0;
//  auto pp = deviceStatus->GetboundaryConditions();
//  auto pppp = pp->GetConditionPropertiesSimple(1);
//  std::cout << pp->GetConditionPropertiesSimple(0)[0] << std::endl;
//  std::cout << pp->GetConditionPropertiesSimple(1)[0] << std::endl;
//  std::cout << pp->GetConditionPropertiesSimple(2)[0] << std::endl;
  boundaryPoints = fieldSolver->getNearBoundarypoints_all();

  status.push_back("Initializations of another solvers");

  particleGridInterface->init(
      deviceStatus->GetNumberParticlesFlows(), deviceStatus->GetGridData(),
      deviceStatus->Getmesh()->templNumb, deviceStatus->Getmesh()->flagMatrix, boundaryPoints,
      deviceStatus->GetDomainBoundary(), memorySize, numThreads, blockSize);

  emissionCurrentSolverPIC->reset();
  std::vector<unsigned int> t;

  for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++) {
    emissionCurrentSolverPIC->SetValueOnSource(
        deviceStatus->GetFlow(i)->GetEmitterDevice()->GetParticleSources()[0], {1.0}, i, 2);
    deviceStatus->GetFlow(i)->GenerateParticles(t, 0, 0, 0, deviceStatus->GetGridData(), 0, 0);
    particleGridInterface->SearchStartCells(i,
                                            deviceStatus->GetFlow(i)->GetDynamicsDataStart());
    //		particleGridInterface->SearchBoundariesCells(deviceStatus->GetboundariesForFlows());
    emissionCurrentSolverPIC->init(
        i, deviceStatus->GetGridData(),
        deviceStatus->GetFlow(i)->GetEmitterDevice()->GetEmissionType(),
        deviceStatus->GetFlow(i)->GetEmitterDevice()->GetParticleSource(),
        particleGridInterface, deviceStatus->GetGridData()->GetType());
    alpas[i] = deviceStatus->GetFlow(i)->GetAlpha();
  }
  particleMover->init(alpas, 1);

  maxtime = 0;

  for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++) {
    emissionCurrentSolverPIC->SetValueOnSource(
        deviceStatus->GetFlow(i)->GetEmitterDevice()->GetParticleSources()[0], {1.0}, i, 2);
    if (deviceStatus->GetFlow(i)->GetFlowProperties()[4] > maxtime)
      maxtime = deviceStatus->GetFlow(i)->GetFlowProperties()[4];
  }
  progressLoc = 1;
  auto qq = deviceStatus->GetboundaryConditions();
  auto ff = deviceStatus->GetglobalFieldConditions();
  if (!deviceStatus->isFieldSimulated) {
    status.push_back("External field simulation");
    auto pp = deviceStatus->GetboundaryConditions();
    fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
                               deviceStatus->GetboundaryConditions(),
                               deviceStatus->GetFlow(0)->GetDynamicsData(0)->Time /
                               (LIGHT_VELOCITY()),
                               progressLoc);
    deviceStatus->isFieldSimulated = 1;
    progressLoc = 1;
  }
  progressLoc = 0;

  if (flagRestart) {
    status.push_back("Timestep estimation");

    deviceStatus->GetGridData()->densityReset();
    std::vector<double> dt;

    if (particleMover->params[1] == 2) {
      dt.resize(deviceStatus->GetNumberParticlesFlows());
      for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
        dt[i] = LIGHT_VELOCITY() / (deviceStatus->GetglobalFieldConditions()[0] *
                                               particleMover->GetParameters()[0]);
    } else
      TimeStepEstimate(dt, deviceStatus, numThreads, progressLoc);

    volatile double minDt = dt[0];

    for (int i = 1; i < dt.size(); i++) {
      if (dt[i] < minDt)
        minDt = dt[i];
    };
    for (int i = 1; i < dt.size(); i++)
      dt[i] = minDt;

    particleMover->SetTimeSteps(dt);

    for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++)
      deviceStatus->GetFlow(i)->clearParallel();

    for (int i = 0; i < deviceStatus->GetElectrodes().size(); i++)
      deviceStatus->GetElectrodes()[i]->ResetPower();
  }

  progressLoc = 1;

  particleMover->flagInit = 1;
};

template void Solver<double>::CheckConfiguration<device2daxsdouble>(
    const std::shared_ptr<device2daxsdouble> &deviceStatus, int flagRestart, std::string &errorMsg);

template void Solver<float>::CheckConfiguration<device2daxsfloat>(
    const std::shared_ptr<device2daxsfloat> &deviceStatus, int flagRestart, std::string &errorMsg);

template void Solver<double>::CheckConfiguration<device2ddouble>(
    const std::shared_ptr<device2ddouble> &deviceStatus, int flagRestart, std::string &errorMsg);

template void
Solver<float>::CheckConfiguration<device2dfloat>(const std::shared_ptr<device2dfloat> &deviceStatus,
                                                 int flagRestart, std::string &errorMsg);

template void Solver<double>::CheckConfiguration<device2dpolardouble>(
    const std::shared_ptr<device2dpolardouble> &deviceStatus, int flagRestart,
    std::string &errorMsg);

template void Solver<float>::CheckConfiguration<device2dpolarfloat>(
    const std::shared_ptr<device2dpolarfloat> &deviceStatus, int flagRestart,
    std::string &errorMsg);

template void Solver<double>::CheckConfiguration<device3ddouble>(
    const std::shared_ptr<device3ddouble> &deviceStatus, int flagRestart, std::string &errorMsg);

template void Solver<float>::CheckConfiguration<device3dExtrfloat>(
    const std::shared_ptr<device3dExtrfloat> &deviceStatus, int flagRestart, std::string &errorMsg);

template void Solver<double>::CheckConfiguration<device3dExtrdouble>(
    const std::shared_ptr<device3dExtrdouble> &deviceStatus, int flagRestart,
    std::string &errorMsg);

template<class PointType>
template<class deviceType>
void Solver<PointType>::CheckConfiguration(const std::shared_ptr<deviceType> &deviceStatus,
                                           int flagRestart, std::string &errorMsg) {
  errorMsg.clear();

  if (!restartPossible && !flagRestart)
    errorMsg =
        errorMsg + "It is impossible to continue simulations because it will first start\n";

  for (int i = 0; i < parameters.size(); i++) {
    if (parameters[i] <= 0) {
      errorMsg = errorMsg + "Some input parameter is <=0\n";
      break;
    };
  }
  if (!flagRestart && !deviceStatus->isFieldSolverInit)
    errorMsg =
        errorMsg +
        "It is impossible to continue simulations because of field or mesh was been changed\n";

  if (!deviceStatus->Getmesh()->GetVTKGrid())
    errorMsg = errorMsg + "Mesh is not generated\n";

  if (!deviceStatus->GetNumberParticlesFlows())
    errorMsg = errorMsg + "There are no flows\n";

  for (int i = 0; i < deviceStatus->GetNumberParticlesFlows(); i++) {
    if (!deviceStatus->GetFlow(i)->isConfigure())
      errorMsg = errorMsg + "flow " + std::to_string(i) + " is not completely configured\n";
  };

  if (particleMover->params[1] == 2 && std::abs(deviceStatus->GetglobalFieldConditions()[0]) < 1e-9)
    errorMsg =
        errorMsg +
        "Unable to use this timestep estimation method. Global frequency is very small\n";
};
