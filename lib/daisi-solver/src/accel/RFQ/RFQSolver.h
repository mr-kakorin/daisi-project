#ifndef RFQSolver_H
#define RFQSolver_H

#include <armadillo>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "../base/AccelSolverBase.h"

class SimulationDataAccel;
class RFQDevice;
class RFQSolver : public AccelSolverBase
{
  public:
    RFQSolver(const std::string& dataFolderIn);
    RFQSolver(){};
    void GenerateRFQForFlow(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                            std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                            std::string&                                       errorMsg);
    void SimulateDynamics(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                          std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                          std::string& errorMsg, std::vector<std::string>& status);
    void ExportToLidos(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                       std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                       std::string&                                       errorMsg);
    void MatcherOpt(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
                    std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                    std::string& errorMsg, std::vector<std::string>& status);
    void Opt(std::shared_ptr<RFQDevice> device, double& progress, bool& flagAbort,
             std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
             std::vector<std::string>& status);
};

#endif