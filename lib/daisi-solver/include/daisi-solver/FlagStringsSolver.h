#ifndef FLAGSTRINGS_H
#define FLAGSTRINGS_H

#include <iostream>
#include <vector>

namespace flagStringsSolver
{
extern std::vector<std::string> PlotFlags2daxs;
extern std::vector<std::string> PlotFlags2d;
extern std::vector<std::string> PlotFlags2dpolar;
extern std::vector<std::string> PlotFlags3d;

extern std::vector<std::string> flowBoundaryTypeNames;

extern std::vector<std::string> ResultsNames2daxs;
extern std::string              flowBoundaryList;
extern std::string              poisson;

extern std::vector<std::string> simulationDataNamesFlowPIC;

extern std::vector<std::string> simulationDataNamesFlowPICAccel;
extern std::vector<std::string> simulationDataNamesElectrode;
extern std::vector<std::string> simulationDataNamesFlowPTI;

extern std::vector<std::string> simulationDataNamesBasePTI;
extern std::vector<std::string> simulationDataNamesBasePIC;
};
#endif