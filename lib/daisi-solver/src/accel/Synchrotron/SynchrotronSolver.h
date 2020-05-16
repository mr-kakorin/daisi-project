#ifndef SynchrotronSolver_H
#define SynchrotronSolver_H

#include <armadillo>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

#include "../base/AccelSolverBase.h"

class SynchrotronDevice;
class OpticElementsSequence;
class SimulationDataAccel;
class SynchrotronSolver : public AccelSolverBase
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive> void save(Archive& ar, const unsigned int) const;
    template <class Archive> void load(Archive& ar, const unsigned int);
    void Mikado(int flag, std::shared_ptr<SynchrotronDevice>& device, arma::mat& R,
                std::vector<int>& k1, int N, arma::mat& dx0, const arma::mat& x0);

  public:
    SynchrotronSolver(const std::string& dataFolderIn);
    SynchrotronSolver(){};
    void DynamicsSimulationTwiss(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                                 bool&                                              flagAbort,
                                 std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                 std::string& errorMsg, bool flagSave = true);
    void DynamicsSimulationTwissMadx(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                                     bool&                                              flagAbort,
                                     std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                     std::string&                                       errorMsg);
    void DynamicsSimulationBeam(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                                bool&                                              flagAbort,
                                std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                std::string&                                       errorMsg);
    void DynamicsSimulationCenterMass(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                                      bool&                                              flagAbort,
                                      std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                      std::string&                                       errorMsg);

	void ClosedOrbitCalculation(std::shared_ptr<SynchrotronDevice>& device, double& progress,
		bool&                                              flagAbort,
		std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
		std::string&                                       errorMsg);

    void OrbitCorrectionSVD(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                            bool&                                              flagAbort,
                            std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                            std::string&                                       errorMsg);
    void OrbitCorrectionMikado(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                               bool&                                              flagAbort,
                               std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                               std::string&                                       errorMsg);
    void OrbitCalculation(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                          bool&                                              flagAbort,
                          std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                          std::string&                                       errorMsg);
    void AcceptanceCalculation(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                               bool&                                              flagAbort,
                               std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                               std::string&                                       errorMsg);
    void CorrMatrixCalc(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                        bool&                                              flagAbort,
                        std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                        std::string&                                       errorMsg);
    void SynchrotronOptimization(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                                 bool&                                              flagAbort,
                                 std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                 std::string& errorMsg, std::vector<std::string>& status);
    void SynchrotronSolversTolerancesReconstruction(
        std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
        std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
        std::vector<std::string>& status);
    void SynchrotronSolversCorrEff(std::shared_ptr<SynchrotronDevice>& device, double& progress,
                                   bool&                                              flagAbort,
                                   std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
                                   std::string& errorMsg, std::vector<std::string>& status);
    void SynchrotronSolversCorrEffStep(
        std::shared_ptr<SynchrotronDevice>& device, double& progress, bool& flagAbort,
        std::vector<std::shared_ptr<SimulationDataAccel>>& outputData, std::string& errorMsg,
        std::vector<std::string>& status, double sigma);

	void SynchrotronSolversMagnetShuffle(std::shared_ptr<SynchrotronDevice>& device, double& progress,
		bool&                                              flagAbort,
		std::vector<std::shared_ptr<SimulationDataAccel>>& outputData,
		std::string&                                       errorMsg, std::vector<std::string>& status);
};
void mtrans(const std::shared_ptr<OpticElementsSequence>& seq, double rorbit, double lambdak,
            double Brho, std::vector<arma::mat>& Tx, std::vector<arma::mat>& Ty,
            std::vector<arma::mat>& sTx, std::vector<arma::mat>& sTy);
#endif