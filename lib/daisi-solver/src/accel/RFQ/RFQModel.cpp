#include "RFQModel.h"
#include "RFQDevice.h"
#include "RFQSolver.h"
#include "RFQStrings.h"
#include "RFQTools.h"

void RFQModel::Solve(const std::string& solverName, double& progress, bool& flagAbort,
                     std::string& errorMsg, std::vector<std::string>& status)
{
    if (solverName == GenerateRFQForFlowName)
        solver->GenerateRFQForFlow(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SimRFQForFlowName)
        solver->SimulateDynamics(device, progress, flagAbort, outputData, errorMsg, status);

    if (solverName == ExportRFQForFlowName)
        solver->ExportToLidos(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == MatcherOptName)
        solver->MatcherOpt(device, progress, flagAbort, outputData, errorMsg, status);

    if (solverName == OptName)
        solver->Opt(device, progress, flagAbort, outputData, errorMsg, status);

    /*	switch ()
            {
            //case 0: solver->DynamicsSimulationTwiss(device, progress, flagAbort, outputData,
       errorMsg); break;
            case 1: solver->GenerateRFQForFlow(device, progress, flagAbort, outputData, errorMsg);
       break;
            case 2: solver->SimulateDynamics(device, progress, flagAbort, outputData, errorMsg);
       break;
            };*/
    progress = 1.0;
};
void RFQModel::GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                          std::vector<std::string>&         names){};

template void
RFQModel::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                     const unsigned int               file_version);
template void
RFQModel::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                     const unsigned int               file_version);

template void RFQModel::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);

template void
RFQModel::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                const unsigned int file_version) const;

template <class Archive>
void RFQModel::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelModelTemplate>(*this);
};
template <class Archive>
void RFQModel::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelModelTemplate>(*this);
};
