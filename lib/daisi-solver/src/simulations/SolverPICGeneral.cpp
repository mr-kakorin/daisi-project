#include "DataTypes.h"
#include "ElectrodeCurrent.h"
#include "EmissionCurrentSolver.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice3d.h"
#include "FlagStringsSolver.h"
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
#include <iostream>


template <class PointType>
template <class Archive>
void Solver<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& outTimes;
    ar& relaxations;
    ar& parameters;
    ar& emissionCurrentSolverPIC;
    ar& emissionCurrentSolverPTI;
    ar& particleGridInterface;
    ar& particleMover;
    ar& fieldSolver;
    ar& meshGenerator;
}

template <class PointType>
template <class Archive>
void Solver<PointType>::load(Archive& ar, const unsigned int)
{
    ar& outTimes;
    ar& relaxations;
    ar& parameters;
    ar& emissionCurrentSolverPIC;
    ar& emissionCurrentSolverPTI;
    ar& particleGridInterface;
    ar& particleMover;
    ar& fieldSolver;
    ar& meshGenerator;
    restartPossible = 0;
}

template <class PointType>
std::vector<std::vector<double>> Solver<PointType>::GetParametersPIC()
{
    std::vector<std::vector<double>> result;
    result.push_back(outTimes);
    return result;
}

template <class PointType>
void Solver<PointType>::SetParametersPIC(const std::vector<std::vector<double>>& par)
{
    outTimes = par[0];
}

template <class PointType>
std::vector<std::string> Solver<PointType>::GetVisNames(int solver)
{
    if (solver == 0)
        return flagStringsSolver::simulationDataNamesBasePIC;
    else
        return flagStringsSolver::simulationDataNamesBasePTI;

    //	return  flagStringsSolver::simulationDataNames;
}

template <class PointType>
void Solver<PointType>::SetParametersPTI(const std::vector<std::vector<double>>& par)
{
    relaxations = par[0];
}

template <class PointType>
std::vector<std::vector<double>> Solver<PointType>::GetParametersPTI()
{
    std::vector<std::vector<double>> result;
    result.push_back(relaxations);
    return result;
}

template <class PointType>
void Solver<PointType>::SetSolverGeneralParameters(const std::vector<std::vector<double>>& par)
{
    particleGridInterface->SetParameters(par[0]);
    particleMover->SetParameters(par[1]);
    parameters = par[2];
}

template <class PointType>
Solver<PointType>::Solver()
{
    relaxations.resize(4);
    relaxations[0] = 0.08;
    relaxations[1] = 0.3;
    relaxations[2] = 0.3;
    relaxations[3] = 0.3;

    outTimes.resize(1);
    outTimes[0] = 0;

    emissionCurrentSolverPIC = std::shared_ptr<EmissionCurrentSolverPIC<PointType>>(
        new EmissionCurrentSolverPIC<PointType>());
    emissionCurrentSolverPTI = std::shared_ptr<EmissionCurrentSolverPTI<PointType>>(
        new EmissionCurrentSolverPTI<PointType>());

    parameters.resize(5);
    parameters[0] = 10;
    parameters[1] = 0.1;
    parameters[2] = 2;
    parameters[3] = 5000;
    parameters[4] = 1;
    particleGridInterface =
        std::shared_ptr<ParticleGridInterface<PointType>>(new ParticleGridInterface<PointType>());
    particleMover = std::shared_ptr<ParticlesMover<PointType>>(new ParticlesMover<PointType>());
    fieldSolver   = std::shared_ptr<PoissonSolver<PointType>>(new PoissonSolver<PointType>());
    meshGenerator = std::shared_ptr<MeshGenerator<PointType>>(new MeshGenerator<PointType>());
    magneticFieldSolver =
        std::shared_ptr<MagneticFieldSolver<PointType>>(new MagneticFieldSolver<PointType>());
}

template class Solver<float>;
template class Solver<double>;

template void
Solver<float>::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                          const unsigned int file_version);
template void
Solver<double>::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                           const unsigned int file_version);

template void
Solver<double>::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                           const unsigned int file_version);
template void
Solver<float>::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                          const unsigned int file_version);
