#ifndef FLAGSTRINGS_H
#define FLAGSTRINGS_H

#include "FlagStringsSolver.h"
#include <iostream>
#include <vector>

namespace flagStringsSolver
{
std::vector<std::string> PlotFlags2daxs = {"Er, V/m", "Ez, V/m", "Enorm, V/m",
                                           "V, V",    "Bphi, T", "Charge density, cl/m^3"};
std::vector<std::string> PlotFlags2d = {"Ex, V/m", "Ey, V/m", "Enorm, V/m", "V, V",
                                        "Charge density, cl/m^3"};
std::vector<std::string> PlotFlags2dpolar = {"Er, V/m", "Ephi, V/m", "Enorm, V/m", "V, V",
                                             "Charge density, cl/m^3"};
std::vector<std::string> PlotFlags3d = {"Ex, V/m",    "Ey, V/m", "Ez, V/m",
                                        "Enorm, V/m", "Bx, V/m", "By, V/m",
                                        "Bz, V/m",    "V, V",    "Charge density, cl/m^3"};

std::vector<std::string> flowBoundaryTypeNames = {"Transparency", "Reflection", "Backscattering"};

std::vector<std::string> ResultsNames2daxs = {"RZ plane", "R(t)", "Z(z)", "XY plane", "Currents"};
std::string              flowBoundaryList  = "flowBoundaryList";
std::string              poisson           = "poisson";

std::vector<std::string> simulationDataNamesFlowPIC = {
    "N particles(t); of flow", "I(t), A; of flow", "Er_av(t), V/m; of flow"};

std::vector<std::string> simulationDataNamesFlowPICAccel = {
    "N particles(t); of flow", "XdX RMS Emittance, pi*cm*mrad", "YdY RMS Emittance, pi*cm*mrad",
    "dWPh RMS Emittance, kEv*ns"};
std::vector<std::string> simulationDataNamesElectrode = {"Collected current"};
std::vector<std::string> simulationDataNamesFlowPTI   = {"I(iter)", "Er_av(iter)"};

std::vector<std::string> simulationDataNamesBasePTI = {"Error(iter)"};
std::vector<std::string> simulationDataNamesBasePIC = {"Error(t)"};
}
#endif