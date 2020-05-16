#include <common_tools/constants.h>

#include "../base/AccelFlow.h"

#include "Results.h"
#include "SynchrotronDevice.h"
#include "SynchrotronSolver.h"
#include "SynchrotronStrings.h"
#include "SynchrotronTools.h"
#include "Tools.h"

template void
SynchrotronSolver::serialize<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);
template void
SynchrotronSolver::serialize<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                              const unsigned int file_version);

template void
SynchrotronSolver::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                         const unsigned int file_version);

template void
SynchrotronSolver::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                         const unsigned int file_version) const;

template <class Archive>
void SynchrotronSolver::save(Archive& ar, const unsigned int) const
{
    ar& boost::serialization::base_object<AccelSolverBase>(*this);
};
template <class Archive>
void SynchrotronSolver::load(Archive& ar, const unsigned int)
{
    ar& boost::serialization::base_object<AccelSolverBase>(*this);
};

SynchrotronSolver::SynchrotronSolver(const std::string& dataFolderIn)
    : AccelSolverBase(dataFolderIn, 0)
{
    std::vector<std::string>              empty;
    std::vector<std::vector<std::string>> empty1;
    std::vector<std::string>              faik = {"dfdf"};

    AddSolver(SynchrotronSolversNameClosedOrbit,
              std::vector<std::vector<std::string>>{{"", ""},
                                                    {},
                                                    {"*.json", "*.json"},
                                                    {"*.json", "*.json"},
                                                    {},
                                                    {"Tolerances", "NOTK config"}},
              empty, empty1);

    AddSolver(
        SynchrotronSolversNameCenterOfMass,
        std::vector<std::vector<std::string>>{{""}, {}, {"*.json"}, {"*.json"}, {}, {"Tolerances"}},
        empty, empty1);

    AddSolver(SynchrotronSolversNameTwiss,
              std::vector<std::vector<std::string>>{
                  {}, {"twiss_daisi.txt"}, {}, {}, {"file with twiss parameters"}, {}},
              empty, empty1);

    std::vector<std::string> tt1 = {"Number of visualized traces",
                                    "Number of saved traces to file"};

    AddSolver(SynchrotronSolversBeam,
              {{""},
               {"trackone_daisi.txt"},
               {"*.json"},
               {"*.json"},
               {"file with dynamics"},
               {"Tolerances"}},
              tt1, empty1);

    std::vector<std::vector<std::string>> tt2 = {{"Use alignment err", "Not Use alignment err"}};

    AddSolver(SynchrotronSolversMADX,
              {{},
               {"twiss_madx.txt", "trackone_madx.txt"},
               {},
               {},
               {"file with twiss parameters", "file with dynamics "},
               {}},
              faik, tt2);

    AddSolver(SynchrotronSolversCorrMatrix,
              std::vector<std::vector<std::string>>{
                  {}, {"CorrMatrix.dat"}, {}, {}, {"File with correction matrix"}, {}},
              empty, empty1);

    AddSolver(SynchrotronSolversSVD,
              {{""},
               {"SVD_results.txt"},
               {"*.dat"},
               {"*.dat"},
               {"File with correction results"},
               {"Correction matrix"}},
              empty, empty1);

    std::vector<std::string> tt3 = {"Number of correctors (X)", "Number of correctors (Y)"};

    AddSolver(SynchrotronSolversMikado,
              std::vector<std::vector<std::string>>{{""},
                                                    {"Mikado_results.txt"},
                                                    {"*.dat"},
                                                    {"*.dat"},
                                                    {"File with correction results"},
                                                    {"Correction matrix"}},
              tt3, empty1);

    AddSolver(SynchrotronSolversOrbitCalc,
              std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}}, empty, empty1);

    std::vector<std::string> tt4 = {"Number of particles", "dY/dZ search range",
                                    "dX/dZ search range", "Relative channel radius"};

    AddSolver(SynchrotronSolversAcceptance,
              std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}}, tt4, empty1);

    std::vector<std::string> t5 = {"Number of threads", "Maximal number of iterations",
                                   "Population size", "Parameters search range (Relative)"};

    std::vector<std::vector<std::string>> tt5 = {{"Particle swarm and annealing sim", "Genetic"},
                                                 {"Maximal", "RMS"}};

    AddSolver(SynchrotronSolversOptimization,
              std::vector<std::vector<std::string>>{{}, {}, {}, {}, {}, {}}, t5, tt5);

    /*AddSolver(SynchrotronSolversTolerancesNameReconstruction,
              std::vector<std::vector<std::string>>{ {""}, { "Tolerances.dat" }, {
       "*.txt"
       }, {
       "*.txt" }, {"File with recinstructed tolerances"}, { "Weights" }},
              std::vector<std::string>{"Number of threads", "Maximal number of
       iterations",
                                       "Population size", "Relative error"},
              std::vector<std::vector<std::string>>{
                  {"Particle swarm and annealing sim", "Genetic"},
                  {"Vector fitness", "Matrix fitness"},
                  {"Monitors default", "Monitors everywhere", "Monitors
       enumeration"}});*/

    std::vector<std::string> t6 = {"dip. comp. err. sigma start", "dip. comp. err. sigma end",
                                   "number of steps", "number of iterations"};

    AddSolver(SynchrotronSolversCorrEffName,
              std::vector<std::vector<std::string>>{
                  {""}, {}, {"*.madx"}, {"*.madx"}, {}, {"MAD-X config"}},
              t6, empty1);

    std::vector<std::vector<std::string>> tt7{{"Vector fitness", "Matrix fitness"}};

    AddSolver(SynchrotronSolversTolerancesNameReconstruction,
              std::vector<std::vector<std::string>>{
                  {"", "", "", ""},
                  {"Comparison.dat", "Tolerances.json", "Rec_tolerances.json"},
                  {"*.dat", "*.json", "*.dat", "*.dat"},
                  {"*.dat", "*.json"},
                  {"File with tolerances comparison", "File with generated tolerances",
                   "File with reconstructed tolerances"},
                  {"Monitiors", "NOTK config", "Initial data", "Errors"}},
              faik, tt7);

    std::vector<std::string> tt8{"Number of iterations"};

    AddSolver(SynchrotronSolversMagnetShuffleName,
              std::vector<std::vector<std::string>>{{""},
                                                    {"BestShuffes.dat"},
                                                    {"*.txt"},
                                                    {"*.txt"},
                                                    {"File with best shuffes"},
                                                    {"Lenghts"}},
              tt8, empty1);

    BaseInit(dataFolderIn);
};
