#ifndef FLAGSTRINGSD_H
#define FLAGSTRINGSD_H
#include <QString>
#include <iostream>
#include <vector>
namespace flagStrings
{
const static std::vector<std::string> PlotType       = {"External", "Space-Charge"};
const static std::vector<std::string> ShapeTypeNames = {"CIC", "TSC"};
const static std::vector<std::string> MoverTypeNames = {
    "Constant step using CFL", "Adaptive using CFL (iterative only)", "Constant using wave length"};

const static std::vector<std::string> flowBoundaryTypeNames      = {"Transparency", "Reflection", "Backscattering"};
const static std::vector<std::string> flowBoundaryTypeFlagsNames = {"Boundary", "Analytical"};
const static std::vector<std::string> DeviceTypeNames            = {"CPU", "GPU"};

const static std::vector<std::string> VisNames          = {"Beam Cloud"};
const static std::vector<std::string> VisFlowNames      = {"Current Density", "Electric field on emitter"};
const static std::vector<std::string> VisElectrodeNames = {"Average power density", "Average current density"};

const static std::vector<std::string> SpChargeTypeNames = {"Child-Langmuir law", "Gauss's law"};
};
#endif