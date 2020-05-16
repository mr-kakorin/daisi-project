#ifndef DTLSolver_H
#define DTLSolver_H
#include "AccelSolverBase.h"
#include <armadillo>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
class SimulationDataAccel;
class DTLDevice;
class DTLSolver : public AccelSolverBase
{
  public:
    DTLSolver(const std::string& dataFolderIn);
    DTLSolver(){};
    void GenerateDTLForFlow(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                            std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg);
    void SimulateDynamicsEnvelopes(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                                   std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                   std::string&                                       errorMsg);
    void CalculateAcceptance(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                             std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg);
    void SimulateDynamicsBeam(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                              std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg);
    void Export(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg);
    void Optimization(std::shared_ptr<DTLDevice> device, double& progress, bool& flagAbort,
                      std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg);
};

#endif