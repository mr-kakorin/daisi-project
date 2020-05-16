#include "DTLModel.h"
#include "DTLDevice.h"
#include "DTLSolver.h"
#include "DTLStrings.h"
#include "DTLTools.h"

void DTLModel::Solve(const std::string& solverName, double& progress, bool& flagAbort, std::string& errorMsg,
                     std::vector<std::string>& status)
{
    errorMsg.clear();
    device->checkInputParameters(errorMsg);
    solver->checkInputParameters(dataFolder, errorMsg, solverName);
    if (errorMsg.size())
        return;

    if (!device->GetLinacFlows().size())
    {
        errorMsg = "There are no flows";
        return;
    };

    if (solverName == DTLSolverGenerateElectrodes)
        solver->GenerateDTLForFlow(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == DTLSolverSimulateEnvelopes)
        solver->SimulateDynamicsEnvelopes(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == DTLSolverCalculateAcceptance)
        solver->CalculateAcceptance(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == DTLSolverSimulateBeam)
        solver->SimulateDynamicsBeam(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == DTLSolverExport)
        solver->Export(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == DTLOptimization)
        solver->Optimization(device, progress, flagAbort, outputData, errorMsg);

    progress = 1.0;
};
void DTLModel::GetAccelElemetsDescription(std::vector<std::vector<double>>& props, std::vector<std::string>& names){};

template void DTLModel::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                                   const unsigned int               file_version);
template void DTLModel::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                                   const unsigned int               file_version);

template <class Archive> void DTLModel::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelModelTemplate>(*this);
};
template <class Archive> void DTLModel::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelModelTemplate>(*this);
};
