#include "DTLSolver.h"
#include "../base/AccelFlow.h"
#include "AccelToolsGeneral.h"
#include "DTLDevice.h"
#include "DTLStrings.h"
#include "DTLTools.h"
#include "Results.h"
#include "Tools.h"

DTLSolver::DTLSolver(const std::string& dataFolderIn) : AccelSolverBase(dataFolderIn, 0)
{
/*     this->AddSolver(DTLSolverGenerateElectrodes, std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}},
              std::vector<std::string>{"Flow number"}, std::vector<std::vector<std::string>>{});

    this->AddSolver(DTLSolverCalculateAcceptance, std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}},
              std::vector<std::string>{"Flow number", "Number of particles", "Number of particles", "dX/dZ max range",
                                       "dY/dZ max range"},
              std::vector<std::vector<std::string>>{});

    this->AddSolver(DTLSolverSimulateEnvelopes, std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}},
              std::vector<std::string>{"Number of saved particles", "Timestep parameter"},
              std::vector<std::vector<std::string>>{});

    this->AddSolver(DTLSolverSimulateBeam, std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}},
              std::vector<std::string>{"Number of saved particles", "Timestep parameter"},
              std::vector<std::vector<std::string>>{});

    this->AddSolver(DTLSolverExport,
              std::vector<std::vector<std::string>>{{}, {"DTL.mdl"}, {}, {}, {"file with DAISI model"}, {}},
              std::vector<std::string>{}, std::vector<std::vector<std::string>>{});

    this->AddSolver(DTLOptimization, std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}},
              std::vector<std::string>{"Number of threads"}, std::vector<std::vector<std::string>>{}); */
};
void DTLSolver::GenerateDTLForFlow(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                                   std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                   std::string&                                       errorMsg){

};
void DTLSolver::SimulateDynamicsEnvelopes(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                                          std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                          std::string&                                       errorMsg){

};
void DTLSolver::CalculateAcceptance(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                                    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                    std::string&                                       errorMsg){

};

void DTLSolver::SimulateDynamicsBeam(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                                     std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                     std::string&                                       errorMsg){

};
void DTLSolver::Export(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                       std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg)
{

/*     Model2daxsdouble* currentModel = new Model2daxsdouble();

    ///////Enter here code for boundaries and boundary contitions///////////

    ///////////////////////////////////////////////////////////////////////
    std::ofstream                   ofs(outpuFilePaths[DTLSolverExport][0].c_str(), std::ios::out | std::ios::binary);
    boost::archive::binary_oarchive oa(ofs);
    int                             problemType   = 2;
    int                             precisionType = 0;
    std::string                     version       = "1.1.0";
    oa << version;
    oa << 2;
    oa << precisionType;
    oa << *currentModel; */
};
void DTLSolver::Optimization(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                             std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg){

};
