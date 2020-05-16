#ifndef FLAGSTRINGSD_H
#define FLAGSTRINGSD_H
#include <QString>
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

const static std::vector<QString> brouseFlags{"NuclotronOpt", "NuclotronLine"};
const static std::vector<QString> FilesExtensions = {"Optics description, (*.opt)",
                                                     "Sequence description, (*.line)"};

const static std::vector<std::string> fieldParameters = {
    std::string("potential offset, V"), std::string("potential amplitude, V"),
    std::string("frequency, Hz"), std::string("initial phase, rad")};
const static std::vector<std::string> fieldParametersGlobal = {std::string("potential offset, V"),
                                                               std::string("frequency, Hz"),
                                                               std::string("initial phase, rad")};
const static std::vector<std::string> RFQMainParameters = {"Frequency (Hz)", "Channel Radius"};
const static std::vector<std::string> RFQSectionsNames  = {
    "Shaper section", "Gentle buncher", "Former section", "Accelerator", "Macther", "Output"};
const static std::vector<std::string> RFQSectionsNamesFields = {
    "Number of cells", "Output phase", "Output Modulation", "Appr. degr. ph_s", "Appr. degr. mod"};
const static std::vector<std::string> RFQMatcherNames = {"Number of cells", "MaximalRadius"};

const static std::vector<std::string> LinacSimParameters = {
    "Steps per period", "Z grid steps",      "R grid steps",      "Rel Grid length",
    "Channel aperture", "Recalculate param", "Colounb force flag"};

const static std::vector<std::string> RFQControlsNames  = {"Modulation", "Sync phases"};
const static std::vector<std::string> RFQFlowParameters = {
    "Intervane Voltage", "Mass Number",        "Charge Number",   "XdX Emittance",
    "YdY Emittance",     "Average energy",     "Momentum spread", "Current per impulse",
    "N particles Long",  "N particles Transv", "X Max",           "dX Max",
    "Out energy",        "Rel channel radius"};
const static std::vector<std::string> RFQCavityParametersNames = {
    "Cells lengths", "Minimal radii", "Matcher radii", "All radii", "Energy", "Acceleration Eff"};

const static std::vector<std::string> LinacParamsNames = {"Frequency (Hz)",
                                                          "Initial Energy (eV)",
                                                          "Minimum (Lgap/L) Factor",
                                                          "Outer radius of first DT, m",
                                                          "Inner radius of first DT, m",
                                                          "DT rounding radius, m",
                                                          "Resonator radius, m",
                                                          "Voltage between DT, V",
                                                          "Particle mass number",
                                                          "Particle charge number",
                                                          "Number of periods"};

const static std::vector<std::string> PlotType = {"Total", "External", "Space-Charge"};

const static std::vector<std::string> NamesOfNumbersParams2daxs = {
    "Emission type", "Emit period", "Energy distribution", "Along emitter"};
const static std::vector<std::string> NamesOfDistribParamsParams2daxs = {"Temperature, eV"};

const static std::vector<std::string> NamesOfNumbersParamsLinac2daxs = {
    "Emission type", "Emit period", "Particles per bunch", "Number of bunches"};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac2daxs = {
    "Average bunch current, A",
    "Average energy, eV ",
    "Momentum spread, %",
    "Bunches frequency, Hz",
    "Initial Phase, deg",
    "Left phase border, deg",
    "Right phase border, deg",
    "dW/dPh, rad",
    "X, m",
    "dX, mrad",
    "X/dX, rad",
    "Y, m",
    "dY, mrad",
    "Y/dY, rad"};

const static std::vector<std::string> NamesOfNumbersParams2d       = {};
const static std::vector<std::string> NamesOfDistribParamsParams2d = {};

const static std::vector<std::string> NamesOfNumbersParamsLinac2d       = {};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac2d = {};

const static std::vector<std::string> NamesOfNumbersParams2dpolar       = {};
const static std::vector<std::string> NamesOfDistribParamsParams2dpolar = {};

const static std::vector<std::string> NamesOfNumbersParamsLinac2dpolar       = {};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac2dpolar = {};

const static std::vector<std::string> NamesOfNumbersParams3dfromExtrusion       = {};
const static std::vector<std::string> NamesOfDistribParamsParams3dfromExtrusion = {};

const static std::vector<std::string> NamesOfNumbersParamsLinac3dfromExtrusion = {
    "Emission type", "Emit period", "Energy distribution", "Particles per step"};
const static std::vector<std::string> NamesOfDistribParamsParamsLinac3dfromExtrusion = {
    "Average bunch current, A", "Average energy, eV ", "Beam radius, m"};

const static std::vector<std::string> NuclotronMainNames = {"Elements description",
                                                            "Elements sequence", "Accel length"};

const static std::vector<std::string> radioDistributionStyleNames = {"Source-style", "Linac-style"};
const static std::vector<std::string> radioEmissionTypeNames      = {"Space-charge Limited",
                                                                "Constant Density"};
const static std::vector<std::string> problemTypesNames = {"1d",
                                                           "2d Cartesian",
                                                           "2d Cilindrical axisymmetric",
                                                           "2d Polar",
                                                           "3d Cartesian from extrusion",
                                                           "3d Cartesian",
                                                           "3d Cilindrical",
                                                           "RFQ Design",
                                                           "DTL Design",
                                                           "Nuclotron"};
const static std::vector<std::string> precisionNames         = {"double", "float"};
const static std::vector<std::string> radioParticleTypeNames = {"Electon", "Ion"};
const static std::vector<std::string> solverTypeNames        = {
    "Particle tracking iterative", "Particle-in-cell Poisson", "Particle - in - cell FDTD"};

const static std::vector<std::string> FlowParametersNames = {
    "Beam cloud", "Current density distribution", "Electric field on emitter"};
const static std::vector<std::string> ShapeTypeNames = {"CIC", "TSC"};
const static std::vector<std::string> MoverTypeNames = {"Constant automatic step",
                                                        "Adaptive (iterative only)"};

const static std::vector<std::string> DeviceTypeNames = {"CPU", "GPU"};

const static std::vector<std::string> PlotFlags2daxs = {
    "Er, V/m", "Ez, V/m", "Enorm, V/m", "V, V", "Bphi, T", "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags2d = {"Ex, V/m", "Ey, V/m", "Enorm, V/m", "V, V",
                                                     "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags2dpolar = {"Er, V/m", "Ephi, V/m", "Enorm, V/m",
                                                          "V, V", "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags3d = {
    "Ex, V/m",    "Ey, V/m", "Ez, V/m",
    "Enorm, V/m", "Bx, V/m", "By, V/m",
    "Bz, V/m",    "V, V",    "Charge density, cl/m^3"};
const static std::vector<std::string> PlotFlags3dPlane = {"YZ", "XZ", "XY"};

const static std::vector<std::string> flowBoundaryTypeNames = {
    "Metal Absorbtion", "Reflection", "Backscattering", "Absolut Absorbtion"};

const static std::vector<std::string> flowBoundaryTypeFlagsNames = {"Boundary", "Manual"};

const static std::vector<std::string> ResultsNames2daxs = {"RZ plane",
                                                           "Energy (Z)",
                                                           "Energy",
                                                           "Energy phi",
                                                           "XY plane",
                                                           "RZ plane Rotating",
                                                           "RZ plane Rotating 3 Trace",
                                                           "Charge(t)",
                                                           "R(Z)",
                                                           "R(t)",
                                                           "Z(t)",
                                                           "Phi(t)"};
const static std::vector<std::string> ResultsNames2dpolar = {
    "XY plane", "X(t)", "Y(t)", "Currents", "Energy", "XY plane rotate", "Energy phi"};
const static std::vector<std::string> ResultsNames2d = {"XY plane", "X(t)", "Y(t)", "Currents",
                                                        "Energy"};
const static std::vector<std::string> ResultsNames3d = {"X(t)", "Y(t)", "Z(t)",  "X(z)",
                                                        "Y(z)", "X(Y)", "Energy"};

const static std::vector<std::string> RFQResults = {"Energy(z)",
                                                    "Phase(t)",
                                                    "Rx(z)",
                                                    "Ry(z)",
                                                    "R(z)",
                                                    "XdX Ellipses",
                                                    "YdY Ellipses",
                                                    "Channel and Beam",
                                                    "Acceleration and Transmission",
                                                    "defocusing",
                                                    "Velocity"};
const static std::vector<std::string> VisNames     = {"Beam Cloud"};
const static std::vector<std::string> VisFlowNames = {"Current Density",
                                                      "Electric field on emitter"};
const static std::vector<std::string> VisElectrodeNames = {"Power Density", "Energy Density"};

const static std::vector<std::string> SpChargeTypeNames = {"Gauss's law", "Child-Langmuir las"};
};
#endif