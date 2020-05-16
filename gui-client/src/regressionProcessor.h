#ifndef regProc
#define regProc
#include "map"
#include "vector"
#include <string>

void loadXls(std::string& errorMessage, std::map<std::string, std::vector<std::vector<float>>>& Oilwells);
void processOilwells(std::string& errorMessage, double& progress, std::vector<std::vector<float>>& results,
                     std::vector<std::vector<float>>& resultsLineX, std::vector<std::vector<float>>& resultsLineY,
                     std::vector<std::vector<float>>& eps);
#endif
