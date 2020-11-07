#ifndef VTKCOMPONENT_H
#define VTKCOMPONENT_H

#include <functional>
#include <thread>
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkPNGWriter.h>

#include <vtkImageCanvasSource2D.h>
#include <vtkImageCast.h>
#include "colors.h"

class QVTKWidget;
class vtkRenderer;
class vtkActor;
class vtkDataSetMapper;
class vtkRenderWindow;
class vtkRenderWindowInteractor;
class vtkContextView;
class vtkContextActor;
class vtkTable;
class vtkChartXY;
class vtkChartXYZ;
class vtkScalarBarActor;
class vtkCubeAxesActor;
class vtkDataSet;
class vtkFloatArray;
class vtkComponent
{
    vtkSmartPointer<vtkRenderer>                   renderer;
    std::vector<vtkSmartPointer<vtkActor>>         actorArray;
    std::vector<vtkSmartPointer<vtkDataSetMapper>> mapperArray;
    vtkSmartPointer<vtkRenderWindow>               renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor>     renderWindowInteractor;
    vtkSmartPointer<vtkContextView>                view;
    vtkSmartPointer<vtkContextActor>               contextActor;

    std::vector<vtkSmartPointer<vtkTable>> tables;
    vtkSmartPointer<vtkChartXY>            chart;
    vtkSmartPointer<vtkChartXYZ>           chart3d;
    vtkSmartPointer<vtkScalarBarActor>     scalarBar;
    vtkSmartPointer<vtkCubeAxesActor>      cubeAxesActor;
    int                                    flagReset;
    double*                                position;
    double*                                focalpoint;
    int                                    viewType;
    int                                    irender;
    int                                    flagScalarBar;
    int                                    flagcubeAxesActor;
    std::function<void(vtkComponent*)>     dataSetter;
    std::vector<float>                     XdataCurrent;
    std::vector<std::vector<float>>        YdataCurrent;

  public:
    QVTKWidget* graphicsWidget;

    void SaveData2File(std::string name);
    void SaveData2EpsFile(std::string name);

    void CheckType(int viewTypeIn);
    vtkComponent();
    void refresh(int time);

    template <class _Fn, class... _Args> void setDataFunction(_Fn&& _Fx, _Args&&... _Ax)
    {
        dataSetter = std::bind(std::forward<_Fn>(_Fx), std::placeholders::_1, std::forward<_Args>(_Ax)...);
    };
    void setWidget(QVTKWidget* graphicsWidgetIn);
    void addVisualizationDataBoundary(vtkDataSet* input, const float color[3] = colors::whitecolor, int LineWidth = 2);
    void addVisualizationDataPlotsAutoY(const std::vector<std::vector<float>>& Xdata,
                                        const std::vector<std::vector<float>>& Ydata, const float* color, float width,
                                        std::string xName, std::string yName, int n, int Distrtype,
                                        std::vector<int> props, float& ymax, float& ymin, int startPoint = 0,
                                        int endPoint = -1);
	void addVisualizationDataPlotsAutoY(const std::vector<std::vector<double>>& Xdata,
	                                    const std::vector<std::vector<double>>& Ydata, const float* color, float width,
	                                    std::string xName, std::string yName, int n, int Distrtype,
	                                    std::vector<int> props, float& ymax, float& ymin, int startPoint = 0,
	                                    int endPoint = -1);
	void addVisualizationDataPlots(const std::vector<std::vector<double>>& Xdata,
	                                             const std::vector<std::vector<double>>& Ydata, const float* color,
	                                             float width, std::string xName, std::string yName, int n, int Distrtype,
	                                             std::vector<int> props, float yMin, float yMax);
    void addVisualizationDataPlots(const std::vector<std::vector<float>>& Xdata,
                                   const std::vector<std::vector<float>>& Ydata, const float* color, float width,
                                   std::string xName, std::string yName, int n, int Distrtype, std::vector<int> props,
                                   float yMin, float yMax);
    void addVisualizationDataFillRect(const std::vector<float>& Xdata, const std::vector<float>& Ydata,
                                      const float* color, std::string legend = "", float opacity = 0.0);
    void addVisualizationDataPlot(const std::vector<float>& Xdata, const std::vector<float>& Ydata, const float* color,
                                  float width, std::string xName, std::string yName, int flag, int style,
                                  int Markerwidth, std::string legend = "", float opacity = 0.0,
                                  bool flagFixAxis = true);
    void addVisualizationDataPlot(const std::vector<float>& Xdata, const std::vector<float>& Ydata, float width,
                                  std::string xName, std::string yName, std::string Name);

    template <class Datatype>
    void addVisualizationDataPlot(const std::vector<Datatype>& Xdata, const std::vector<Datatype>& Ydata, float width,
                                  std::string xName, std::string yName, std::string Name, float& ymax, float& ymin,
                                  bool flagFixAxis = true, int startPoint = 0, int endPoint = -1);
    void addVisualizationDataPlot(const std::vector<float>& Xdata, const std::vector<std::vector<float>>& Ydata,
                                  std::vector<const float*> colors, std::vector<float> width, std::string xName,
                                  std::string yName, std::vector<std::string> legends);
    void addVisualizationDataCloud(vtkDataSet* input, float, const float* color = &colors::blackcolor[0]);
    void addVisualizationDataMesh(vtkDataSet* input, int flagEdge = 0);
    void addVisualizationDataMesh(vtkDataSet* input, vtkSmartPointer<vtkFloatArray> vtkData, double minVal,
                                  double maxVal, std::string title);
    void setAxis();
    void removeVisualizationData(int n);
    void clear();
    void AddSurfacePlot();
    void GridData();
    void AddText(const std::vector<std::string>& input, std::vector<const float*> colors, int pos);
};
#endif