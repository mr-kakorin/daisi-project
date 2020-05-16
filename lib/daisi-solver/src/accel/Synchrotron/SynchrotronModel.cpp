#include "SynchrotronModel.h"
#include "SynchrotronDevice.h"
#include "SynchrotronSolver.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
/*SynchrotronModel::SynchrotronModel(const std::string& dataFolderIn)
{
        AccelModelTemplate(dataFolderIn);
};*/

void SynchrotronModel::Solve(const std::string& solverName, double& progress, bool& flagAbort,
                             std::string& errorMsg, std::vector<std::string>& status)
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
    device->InitSequenceWithErrors();
    if (!device->GetOpticElementsSequence()->length())
    {
        errorMsg = "Optic elements sequence is not defined";
        return;
    }

	if (solverName == SynchrotronSolversNameClosedOrbit)
		solver->ClosedOrbitCalculation(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversNameCenterOfMass)
        solver->DynamicsSimulationCenterMass(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversNameTwiss)
        solver->DynamicsSimulationTwiss(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversBeam)
        solver->DynamicsSimulationBeam(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversMADX)
        solver->DynamicsSimulationTwissMadx(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversSVD)
        solver->OrbitCorrectionSVD(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversMikado)
        solver->OrbitCorrectionMikado(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversOrbitCalc)
        solver->OrbitCalculation(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversAcceptance)
        solver->AcceptanceCalculation(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversCorrMatrix)
        solver->CorrMatrixCalc(device, progress, flagAbort, outputData, errorMsg);

    if (solverName == SynchrotronSolversOptimization)
        solver->SynchrotronOptimization(device, progress, flagAbort, outputData, errorMsg, status);

    if (solverName == SynchrotronSolversTolerancesNameReconstruction)
        solver->SynchrotronSolversTolerancesReconstruction(device, progress, flagAbort, outputData,
                                                           errorMsg, status);

    if (solverName == SynchrotronSolversCorrEffName)
        solver->SynchrotronSolversCorrEff(device, progress, flagAbort, outputData, errorMsg,
                                          status);

	if (solverName == SynchrotronSolversMagnetShuffleName)
		solver->SynchrotronSolversMagnetShuffle(device, progress, flagAbort, outputData, errorMsg,
			status);

    progress = 1.0;
};

void SynchrotronModel::GetAccelElemetsDescription(std::vector<std::vector<double>>& props,
                                                  std::vector<std::string>&         names)
{
    props.clear();
    names.clear();
    props.resize(3);
    for (int i = 0; i < device->GetOpticElementsSequence()->length(); i++)
    {
        if (device->GetOpticElementsSequence()->GetType(i) != "DRIFT")
        {
            props[0].push_back(device->GetOpticElementsSequence()->getParameter(i, "at"));
            props[1].push_back(device->GetOpticElementsSequence()->getParameter(i, "L"));
            names.push_back(device->GetOpticElementsSequence()->GetType(i));
            // props[2].push_back(device->GetOpticElementsSequence()->getParameter(i, "L"));
        };
    };
};
template void
SynchrotronModel::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                             const unsigned int file_version);
template void
SynchrotronModel::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                             const unsigned int file_version);

template void
SynchrotronModel::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                        const unsigned int file_version);

template void
SynchrotronModel::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                        const unsigned int file_version) const;

template <class Archive>
void SynchrotronModel::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelModelTemplate>(*this);
}
template <class Archive>
void SynchrotronModel::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelModelTemplate>(*this);
    device->saveCorrectors(dataFolder);
    device->saveErrors(dataFolder);
}