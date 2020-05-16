#include "RFQStrings.h"

std::string GenerateRFQForFlowName = "Generate Electrodes";
std::string SimRFQForFlowName      = "Simulate dynamics";
std::string ExportRFQForFlowName   = "Export to LIDOS";
std::string MatcherOptName         = "Matcher optimization";
std::string OptName                = "Full optimization";

std::vector<std::string> RFQSolversNames = {"Load from paramteq", "Generate Electrodes",
                                            "Dynamics Sim", "Optimization", "Export to LIDOS"};
std::vector<std::vector<std::string>> RFQSolverParameteresNames = {
    {}, {"Flow number"}, {"Flow number", "Number of visualized traces"}, {"Number of threads"}, {}};
std::vector<std::vector<std::string>> RFQOutputFiles     = {{}, {}, {}, {}, {"paramteq.dat"}};
std::vector<std::vector<std::string>> RFQInputFiles      = {{""}, {}, {}, {}, {}};
std::vector<std::vector<std::string>> RFQInputFilesProps = {{"*.dat"}, {}, {}, {}, {}};
std::vector<std::vector<std::string>> RFQInputFilesExt   = {{"*.dat"}, {}, {}, {}, {}};
std::vector<std::string>              RFQMainParameters  = {"Frequency, Hz", "Channel Radius",
                                              "Number of sections"};
std::vector<std::string> RFQResults = {"Energy(z)",
                                       "Phase(t)",
                                       "Rx(z)",
                                       "Ry(z)",
                                       "XdX Ellipses",
                                       "YdY Ellipses",
                                       "Acceleration and Transmission"};

std::vector<std::string> RFQSectionNames0 = {"Number of cells", "Maximal Radius"};
std::vector<std::string> RFQSectionNames  = {"Number of cells", "Output phase", "Appr. deg. phase",
                                            "Output modulation", "Appr. deg. mod"};