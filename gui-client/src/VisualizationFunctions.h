#ifndef VisualizationFunctions_H
#define VisualizationFunctions_H
#include <memory>
#include <string>
#include <vector>
class vtkComponent;
class ModelInterface;
class AccelModelInterface;
class DynamicsData;
class lineplot;
void SimpleShowPoints(vtkComponent* input, const std::vector<float>& dataX, const std::vector<float>& dataY,
                      const std::vector<std::vector<float>>& dataXLine,
                      const std::vector<std::vector<float>>& dataYLine);
void ShowRFQRFQCavityParametersPlots(vtkComponent* input, ModelInterface* currentModel, int flag);
void ShowLinacControls(vtkComponent* input, ModelInterface* currentModel, int flag);
void ShowBoundaries(vtkComponent* input, ModelInterface* currentModel, const std::vector<int>& list1,
                    const float* color1, const std::vector<int>& list2, const float* color2);
void HighLigthBoundary(vtkComponent* input, ModelInterface* currentModel, int n, const float* color);
void ShowMesh(vtkComponent* input, ModelInterface* currentModel, const float* color);
void ShowParticlesCloud(vtkComponent* input, ModelInterface* currentModel, int flowNumber, int precisionType,
                        int problemType);
void ShowPlotXY(vtkComponent* input, std::vector<float> X, std::vector<float> Y, std::string xName, std::string yName);
void ShowCurrentDensity(vtkComponent* input, ModelInterface* currentModel, int flowNumber);
void ShowEmitterField(vtkComponent* input, ModelInterface* currentModel, int flowNumber);
void ShowLinacResultsTraces(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                            int flag, int Distrtype, std::vector<int> props, int numberOfTraces);
void ShowLinacResultsChar(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                          int flag, int Distrtype, std::vector<int> props, int numberOfTraces);
void ShowPositionsWithGeometryPlane(vtkComponent* input, ModelInterface* currentModel,
                                    std::shared_ptr<DynamicsData> data, int flowNumber, std::string xName,
                                    std::string yName, int problemType, int numberOfTraces, int flag, int Distrtype,
                                    std::vector<int> props);
void Show3PositionsWithGeometryPlane(vtkComponent* input, ModelInterface* currentModel,
                                     std::shared_ptr<DynamicsData> data, int flowNumber, std::string xName,
                                     std::string yName, int problemType, int numberOfTraces, int flag, int Distrtype,
                                     std::vector<int> props);
void ShowRT(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data, int flowNumber,
            int numberOfTraces, int Distrtype, std::vector<int> props);
void ShowXY(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data, int flowNumber,
            int numberOfTraces, int Distrtype, std::vector<int> props);
void ShowEnergyEphiPolar(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                         int flowNumber, int numberOfTraces, int Distrtype, std::vector<int> props);
void ShowEnergyPolar(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                     int flowNumber, int numberOfTraces, int Distrtype, std::vector<int> props);
void ShowAllCloud(vtkComponent* input, ModelInterface* currentModel, int precisionType, int problemType);
void ShowPlotLine(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<lineplot>& plotData);
void ShowLinePlot(vtkComponent* input, ModelInterface* currentModel, int plotNumber,
                  std::shared_ptr<lineplot>& plotData);
void Show2dPlot(vtkComponent* input, ModelInterface* currentModel, std::string flag, int precisionType,
                int PlotTypeFlag);
void Show2dPlot3d(vtkComponent* input, ModelInterface* currentModel, std::string flag, int flag1, double param,
                  int precisionType, int PlotTypeFlag);
void ShowSimulationDataPlot(vtkComponent* input, ModelInterface* currentModel, int nSimData, int plotType, int nflow,
                            int nplot, const std::string& yName);
void ShowValueAlongConductor(vtkComponent* input, ModelInterface* currentModel, int conductor, int flag);
void ShowErrors(vtkComponent* input, ModelInterface* currentModel);
void ShowPositionsWithGeometryPlaneSpecial(vtkComponent* input, ModelInterface* currentModel,
                                           std::shared_ptr<DynamicsData> data, int flowNumber, std::string xName,
                                           std::string yName, int problemType, int numberOfTraces, int Distrtype,
                                           std::vector<int> props);
void ShowEmittanceDataPlot(vtkComponent* input, ModelInterface* currentModel, int flowNumber, int flagNumber,
                           int emFlag, int nCharts);
void Show3dTrace(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data, int flowNumber,
                 int numberOfTraces, int flag, int Distrtype, std::vector<int> props);
void ShowTracesT(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                 std::string xName, std::string yName, int position, int flowNumber, int numberOfTraces, int Distrtype,
                 std::vector<int> props, int problemType);
void ShowEnergy(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                std::string xName, int flowNumber, int numberOfTraces, int Distrtype, std::vector<int> props,
                int axsFlag);
void SimpleShowChart(vtkComponent* input, const std::vector<float>& data);
#endif