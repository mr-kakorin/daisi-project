#ifndef RFQStrings_H
#define RFQStrings_H

#include <string>
#include <vector>

extern std::string GenerateRFQForFlowName;
extern std::string SimRFQForFlowName;
extern std::string ExportRFQForFlowName;
extern std::string MatcherOptName;
extern std::string OptName;

extern std::vector<std::string>              RFQSolversNames;
extern std::vector<std::vector<std::string>> RFQSolverParameteresNames;
extern std::vector<std::vector<std::string>> RFQOutputFiles;
extern std::vector<std::vector<std::string>> RFQInputFiles;
extern std::vector<std::vector<std::string>> RFQInputFilesProps;
extern std::vector<std::vector<std::string>> RFQInputFilesExt;
extern std::vector<std::string>              RFQMainParameters;
extern std::vector<std::string>              RFQResults;

extern std::vector<std::string> RFQSectionNames0;
extern std::vector<std::string> RFQSectionNames;

#endif