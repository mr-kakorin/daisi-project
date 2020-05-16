#ifndef VisualizationFunctionsAccel_H
#define VisualizationFunctionsAccel_H
#include <memory>
#include <vector>
class vtkComponent;
class ModelInterface;
class AccelModelInterface;
class SimulationDataAccel;
void ShowLinacSequence(vtkComponent* input, AccelModelInterface* currentModel, int flag);
void ShowAccelCalcParameters(vtkComponent* input, AccelModelInterface* currentModel, int flag);
void ShowSimulationNuclDataPlot(vtkComponent* input, AccelModelInterface* currentModel,
                                std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                                std::vector<int> props, std::vector<double> plotProperties);
void ShowSimulationAccelDataPlot(vtkComponent* input, AccelModelInterface* currentModel,
                                 std::shared_ptr<SimulationDataAccel> simData, const std::string& flag,
                                 std::vector<int> props, std::vector<double> plotProperties);
#endif