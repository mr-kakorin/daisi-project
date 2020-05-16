#ifndef VTKCOMPONENT_H
#define VTKCOMPONENT_H
#include "VTKInclude.h"
#include "colors.h"
#include <functional>
#include <thread>
int flagReset;
class vtkComponent
{

    vtkSmartPointer<vtkRenderer>                   renderer;
    std::vector<vtkSmartPointer<vtkActor>>         actorArray;
    std::vector<vtkSmartPointer<vtkDataSetMapper>> mapperArray;
    vtkSmartPointer<vtkRenderWindow>               renderWindow;
    vtkSmartPointer<vtkRenderWindowInteractor>     renderWindowInteractor;
    std::thread                                    vizualization_thread;
    vtkSmartPointer<vtkMutexLock>                  mutex;
    double*                                        position;
    double*                                        focalpoint;
    int                                            viewType;
    int                                            irender;
    void CheckType(int viewTypeIn)
    {
        if (viewTypeIn != viewType || (viewTypeIn == 1 && viewType == 1))
        {
            renderWindow->RemoveRenderer(renderer);
            renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
            renderWindow->AddRenderer(renderer);
            graphicsWidget->SetRenderWindow(renderWindow);
        }
        viewType = viewTypeIn;
    };

  public:
    QVTKWidget* graphicsWidget;
    vtkComponent()
    {
        viewType     = 0;
        irender      = 0;
        renderer     = vtkSmartPointer<vtkRenderer>::New();
        renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        mutex                  = vtkSmartPointer<vtkMutexLock>::New();
        position               = new double();
        focalpoint             = new double();
    };
    std::function<void(vtkComponent*)> dataSetter;
    void refresh(int time)
    {
        clear();
        render(time);
    };
    void render(int time)
    {
        int i = 0;
        if (time == 0)
        {
            std::this_thread::sleep_for(time);
            mutex->Lock();
            dataSetter(this);
            QMetaObject::invokeMethod(graphicsWidget, "update", Qt::QueuedConnection);
            irender = 1;
            mutex->Unlock();
            flagReset = 1;
        }
        else
        {
            while (1)
            {
                std::this_thread::sleep_for(time);
                mutex->Lock();
                dataSetter(this);
                QMetaObject::invokeMethod(graphicsWidget, "update", Qt::QueuedConnection);
                irender = 1;
                mutex->Unlock();
                if (i == 0)
                    flagReset = 1;
                i++;
            };
        }
    };
    void reset()
    {
        renderer->ResetCamera();
    };
    template <class _Fn, class... _Args> void setDataFunction(_Fn&& _Fx, _Args&&... _Ax)
    {
        dataSetter = std::bind(std::forward<_Fn>(_Fx), std::placeholders::_1, std::forward<_Args>(_Ax)...);
    };
    void setWidget(QVTKWidget* graphicsWidgetIn)
    {
        graphicsWidget = graphicsWidgetIn;
        graphicsWidget->SetRenderWindow(renderWindow);
    };
    void addVizualizationDataBoundary(vtkDataSet* input, const float color[3] = colors::whitecolor, int LineWidth = 2)
    {
        CheckType(0);
        mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
        mapperArray.back()->SetInputData(input);
        actorArray.push_back(vtkSmartPointer<vtkActor>::New());
        actorArray.back()->SetMapper(mapperArray.back());
        actorArray.back()->GetProperty()->SetColor(color[0], color[1], color[2]);
        actorArray.back()->GetProperty()->SetLineWidth(LineWidth);
        mutex->Lock();
        renderer->AddActor(actorArray.back());
        mutex->Unlock();
    };
    void addVizualizationDataPlot(const std::vector<float>& Xdata, const std::vector<std::vector<float>>& Ydata,
                                  std::vector<const float*> colors, std::vector<float> width, std::string xName,
                                  std::string yName, std::vector<std::string> legends)
    {
        CheckType(1);
        clear();

        int nPlots = int(Ydata.size());

        vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

        vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("X Axis");
        table->AddColumn(arrX);

        for (int i = 0; i < nPlots; i++)
        {
            vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
            arrY->SetName(legends[i].c_str());
            table->AddColumn(arrY);
        };

        int numPoints = int(Xdata.size());
        table->SetNumberOfRows(numPoints);

        for (int i = 0; i < numPoints; i++)
        {
            table->SetValue(i, 0, Xdata[i]);
            for (int j = 0; j < nPlots; j++)
            {
                table->SetValue(i, j + 1, Ydata[j][i]);
            }
        }

        vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
        view->GetRenderer()->SetBackground(1.0, 1.0, 1.0);

        vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
        view->GetScene()->AddItem(chart);
        vtkPlot* line = chart->AddPlot(vtkChart::LINE);

        for (int i = 0; i < nPlots; i++)
        {
            line->SetInputData(table, 0, i + 1);
            line->SetColor(colors[i][0], colors[i][1], colors[i][2]);
            line->SetWidth(width[i]);
        }

        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
        chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());

        view->SetRenderWindow(renderWindow);
    }
    void addVizualizationDataCloud(vtkDataSet* input)
    {
        CheckType(0);

        mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
        mapperArray.back()->SetInputData(input);
        actorArray.push_back(vtkSmartPointer<vtkActor>::New());
        actorArray.back()->SetMapper(mapperArray.back());
        actorArray.back()->GetProperty()->SetPointSize(0.5);
        renderer->AddActor(actorArray.back());
    };
    void addVizualizationDataMesh(vtkDataSet* input)
    {
        CheckType(0);
        mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
        mapperArray.back()->SetInputData(input);
        actorArray.push_back(vtkSmartPointer<vtkActor>::New());
        actorArray.back()->SetMapper(mapperArray.back());
        actorArray.back()->GetProperty()->SetColor(0.0, 1.0, 0.0);
        actorArray.back()->GetProperty()->EdgeVisibilityOn();
        renderer->AddActor(actorArray.back());
    };
    void removeVizualizationData(int n)
    {
        renderer->RemoveActor(actorArray[n]);
        mapperArray.erase(mapperArray.begin() + n);
        actorArray.erase(actorArray.begin() + n);
    };
    void clear()
    {
        /*	for (int i = 0; i < actorArray.size(); i++)
                {
                        renderer->RemoveActor(actorArray[i]);
                };*/
        actorArray.clear();
        mapperArray.clear();

        /*renderer = vtkSmartPointer<vtkRenderer>::New();
        renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        graphicsWidget->SetRenderWindow(renderWindow);*/
    };
};
#endif