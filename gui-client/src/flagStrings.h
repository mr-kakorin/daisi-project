#ifndef FLAGSTRINGS_H
#define FLAGSTRINGS_H
#include <iostream>
#include <vector>
namespace flagStrings
{
const static std::string model               = "model";
const static std::string glDef               = "glDef";
const static std::string device              = "device";
const static std::string boundaries          = "boundaries";
const static std::string boundariesList      = "boundariesList";
const static std::string potentialList       = "potentialList";
const static std::string poisson             = "poisson";
const static std::string neumann             = "neumann";
const static std::string globalField         = "globalField";
const static std::string mesh                = "mesh";
const static std::string flows               = "flows";
const static std::string flowList            = "flowList";
const static std::string emitterList         = "emitterList";
const static std::string beamState           = "beamState";
const static std::string vizuaization        = "vizuaization";
const static std::string fieldSolver         = "fieldSolver";
const static std::string solver              = "solver";
const static std::string flowBoundary        = "flowBoundary";
const static std::string flowBoundaryList    = "flowBoundaryList";
const static std::string fabsopbtion          = "fabsopbtion";
const static std::string simulatePIC         = "simulatePIC";
const static std::string simulatePTI         = "simulatePTI";
const static std::string results             = "results";
const static std::string devicestate         = "devicestate";
const static std::string plots2d             = "plots2d";
const static std::string lineplotsList       = "lineplotsList";
const static std::string lineplots           = "lineplots";
const static std::string conductors          = "conductors";
const static std::string conductorsList      = "conductorsList";
const static std::string solverSettings      = "solverSettings";
const static std::string estimateErrors      = "estimateErrors";
const static std::string generateGeometry    = "generateGeometry";
const static std::string solverEmission      = "solverEmission";
const static std::string rfq                 = "rfq";
const static std::string flowsAccel          = "flowsAccel";
const static std::string flowListAccel       = " flowListAccel";
const static std::string RFQCavityParameters = "RFQCavityParameters";
const static std::string rfqSolverParams     = "rfqSolverParams";
const static std::string RFQOpt              = "RFQOpt";
const static std::string NuclAccelParams     = "NuclAccelParams";
const static std::string AccelSolvers        = "AccelSolvers";

const static std::string AccelParams     = "AccelParams";
const static std::string AccelParamsCalc = "AccelParamsCalc";

const static std::vector<std::string> radioDistributionStyleNames = {"Source-style", "Linac-style"};
const static std::vector<std::string> radioEmissionTypeNames      = {"Space-charge Limited", "Constant Density"};
const static std::vector<std::string> problemTypesNames = {"1d",       "2d Cartesian", "2d Cilindrical axisymmetric",
                                                           "2d Polar", "3d Cartesian", "3d Cilindrical"};
const static std::vector<std::string> precisionNames         = {"double", "float"};
const static std::vector<std::string> radioParticleTypeNames = {"Electon", "Ion"};
const static std::vector<std::string> solverTypeNames = {"Particle tracking iterative", "Particle-in-cell Poisson",
                                                         "Particle - in - cell FDTD"};
const static std::vector<std::string> ParticlesNumberNames2d    = {"Energy distribution", "Along emitter", "XY plane"};
const static std::vector<std::string> ParticlesNumberNames2daxs = {"Energy distribution", "Along emitter", "RZ plane",
                                                                   "RPhi plane"};
const static std::vector<std::string> FlowParametersNames = {"Beam cloud", "Current density distribution"};
const static std::vector<std::string> ShapeTypeNames      = {"CIC", "TSC"};
const static std::vector<std::string> MoverTypeNames      = {"Leap-Frog", "Boris"};

const static std::vector<std::string> ParticlesNumberDistributionParameters2d = {
    "Average energy, eV", "Energy spread, eV", "XY plane angle spread, deg"};
const static std::vector<std::string> ParticlesNumberDistributionParameters2daxs = {
    "Average energy, eV", "Energy spread, eV", "RZ plane angle spread, deg", "RPhi plane angle spread, deg",
    "Distribution parameter"};
const static std::vector<std::string> PlotFlags2daxs = {"Er, V/m", "Ez, V/m", "V, V", "Bphi, T",
                                                        "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags2d = {"Ex, V/m", "Ey, V/m", "V, V", "Charge density, cl/m^3"};

const static std::vector<std::string> flowBoundaryTypeNames = {"Metal Absorbtion", "Reflection"};

const static std::vector<std::string> ResultsNames2daxs = {"RZ plane", "R(t)", "Z(z)", "XY plane", "Currents"};
};
#endif