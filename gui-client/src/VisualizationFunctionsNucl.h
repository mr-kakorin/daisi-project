#ifndef VisualizationFunctionsNucl_H
#define VisualizationFunctionsNucl_H
#include <memory>
#include <string>
#include <vector>
class vtkComponent;
class ModelInterface;
class AccelModelInterface;
class SimulationDataAccel;

void ShowSimulationNuclTwiss(vtkComponent* input, AccelModelInterface* currentModel,
                             std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                             std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                             float& ymin);
void ShowSimulationNuclBeam(vtkComponent* input, AccelModelInterface* currentModel,
                            std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                            std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                            float& ymin);
void ShowNuclElements(vtkComponent* input, AccelModelInterface* currentModel,
                      std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                      std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin);
void ShowNuclCorrOrbit(vtkComponent* input, AccelModelInterface* currentModel,
                       std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                       std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin);
void ShowSimulationNuclCM(vtkComponent* input, AccelModelInterface* currentModel,
                          std::shared_ptr<SimulationDataAccel> simData, const std::string& flag, std::vector<int> props,
                          std::vector<double> plotProperties, int& flagSeq, float& ymax, float& ymin);
void ShowSimulationNuclOrbit(vtkComponent* input, AccelModelInterface* currentModel,
                             std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                             std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                             float& ymin);
void ShowSimulationNuclAcceptance(vtkComponent* input, AccelModelInterface* currentModel,
                                  std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                                  std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
                                  float& ymin);
void SynchrotronRec(vtkComponent* input, AccelModelInterface* currentModel,
	std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
	std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
	float& ymin);
void SynchrotronEst(vtkComponent* input, AccelModelInterface* currentModel,
	std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
	std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
	float& ymin);
void SynchrotronShuffle(vtkComponent* input, AccelModelInterface* currentModel,
	std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
	std::vector<int> props, std::vector<double> plotProperties, int& flagSeq, float& ymax,
	float& ymin);

#endif