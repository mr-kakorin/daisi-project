#include "vtkComponent.h"

#include "VTKInclude.h"

#include "vtkBrush.h"
#include "vtkContext2D.h"
#include "vtkContextActor.h"
#include "vtkContextItem.h"
#include "vtkContextScene.h"
#include "vtkObjectFactory.h"
#include "vtkOpenGLContextDevice2D.h"
#include "vtkPen.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkStdString.h"
#include "vtkTextProperty.h"
// double background[] = { 128.0/255,128.0 / 255,128.0 / 255 };
// double backgroundLinePlot[] = { 128.0 / 255,128.0 / 255,128.0 / 255 };

double backgroundLinePlot[] = {255, 255, 255};
double background[]         = {255, 255, 255};

int labelFont = 22;
int titleFont = 22;

/*class APIDiagram : public vtkPlot
{
public:
        static APIDiagram *New();
        vtkTypeMacro(APIDiagram, vtkContextItem);
        // Paint event for the chart, called whenever the chart needs to be drawn
        virtual bool Paint(vtkContext2D *paintera);
};

vtkStandardNewMacro(APIDiagram);
// This function draws our API diagram
bool APIDiagram::Paint(vtkContext2D *painter)
{
        chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(1.1*ymax);

        // Drawing a hard wired diagram 800x600 as a demonstration of the 2D API
        painter->GetTextProp()->SetColor(0.0, 0.0, 0.0);
        painter->GetTextProp()->SetFontSize(24);
        painter->GetPen()->SetColor(0, 0, 0);
        painter->GetTextProp()->Get
        painter->GetBrush()->SetColor(100, 255, 100);
        painter->DrawRect(100, 50, 200, 100);
        painter->DrawString(200, 100, "OpenGL");

        painter->GetBrush()->SetColor(255, 100, 0);
        painter->DrawRect(300, 50, 200, 100);
        painter->DrawString(400, 100, "Others?");

        painter->GetBrush()->SetColor(100, 0, 255);
        painter->DrawRect(500, 50, 200, 100);
        painter->DrawString(600, 100, "Others?");

        painter->GetBrush()->SetColor(180, 180, 255);
        painter->DrawRect(100, 150, 600, 100);
        painter->DrawString(400, 200, "2D API");

        painter->GetBrush()->SetColor(255, 255, 180);
        painter->DrawRect(100, 250, 600, 200);
        painter->DrawString(400, 400, "Canvas API");

        painter->GetBrush()->SetColor(180, 255, 180);
        painter->DrawRect(100, 250, 300, 100);
        painter->DrawString(250, 300, "Point Mark");

        painter->GetBrush()->SetColor(255, 255, 255);
        painter->DrawRect(100, 450, 600, 100);
        painter->DrawString(400, 500, "Canvas View");

        return true;
};
*/

void vtkComponent::SaveData2EpsFile(std::string name)
{
    vtkSmartPointer<vtkGL2PSExporter> vtext;
    vtext = vtkSmartPointer<vtkGL2PSExporter>::New();
    // vtext = vtkGL2PSExporter();
    vtext->SetFileFormatToEPS();
    vtext->SetLineWidthFactor(1.5);
    vtext->SetRenderWindow(renderWindow);
    vtext->SetFilePrefix(name.c_str());
    vtext->Write();
};

void vtkComponent::AddText(const std::vector<std::string>& input, std::vector<const float*> colors, int pos)
{
    for (int i = 0; i < input.size(); i++)
    {
        vtkSmartPointer<vtkTextActor> textActor = vtkSmartPointer<vtkTextActor>::New();
        textActor->SetInput(input[i].c_str());
        textActor->SetPosition(pos * i, 0);
        textActor->GetTextProperty()->SetFontSize(12);
        textActor->GetTextProperty()->SetColor(colors[i][0], colors[i][1], colors[i][2]);
        view->GetRenderer()->AddActor2D(textActor);
    }
}

void vtkComponent::SaveData2File(std::string name)
{
    FILE* fp = fopen(name.c_str(), "w");
    for (int i = 0; i < XdataCurrent.size(); i++)
    {
        fprintf(fp, "%.5f\t", XdataCurrent[i]);
        for (int j = 0; j < YdataCurrent.size(); j++)
            fprintf(fp, "%.5f\t", YdataCurrent[j][i]);
        fprintf(fp, "\n");
    };
    fclose(fp);
}
void vtkComponent::CheckType(int viewTypeIn)
{

    if (viewTypeIn != viewType || (viewTypeIn == 1) || (viewTypeIn == 2))
    {
        renderWindow->RemoveRenderer(renderer);
        renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
        renderWindow->AddRenderer(renderer);
        renderer->SetBackground(background);
        graphicsWidget->SetRenderWindow(renderWindow);
    }
    viewType = viewTypeIn;

    if (viewTypeIn == 1)
    {
        view = vtkSmartPointer<vtkContextView>::New();
        view->GetRenderer()->SetBackground(backgroundLinePlot);
        chart = vtkSmartPointer<vtkChartXY>::New();
        view->GetScene()->AddItem(chart);
        view->SetRenderWindow(renderWindow);
        return;
    }
    if (viewTypeIn == 2)
    {
        view    = vtkSmartPointer<vtkContextView>::New();
        chart3d = vtkSmartPointer<vtkChartXYZ>::New();
        view->GetRenderer()->SetBackground(backgroundLinePlot);
        view->GetScene()->AddItem(chart3d.GetPointer());
        view->SetRenderWindow(renderWindow);
        view->GetRenderWindow()->SetMultiSamples(0);

        return;
    }
};
vtkComponent::vtkComponent()
{
    flagReset         = 0;
    viewType          = 0;
    irender           = 0;
    flagScalarBar     = 0;
    flagcubeAxesActor = 0;
    renderer          = vtkSmartPointer<vtkRenderer>::New();
    renderWindow      = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    renderer->SetBackground(background);
    renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    contextActor           = vtkSmartPointer<vtkContextActor>::New();
    view                   = vtkSmartPointer<vtkContextView>::New();
    view->GetRenderer()->SetBackground(backgroundLinePlot);
    chart         = vtkSmartPointer<vtkChartXY>::New();
    chart3d       = vtkSmartPointer<vtkChartXYZ>::New();
    scalarBar     = vtkSmartPointer<vtkScalarBarActor>::New();
    cubeAxesActor = vtkSmartPointer<vtkCubeAxesActor>::New();
    position      = new double();
    focalpoint    = new double();
};
void vtkComponent::refresh(int time)
{
    dataSetter(this);
    graphicsWidget->update();
    irender = 1;
    flagReset++;
    if (flagReset == 1)
        renderer->ResetCamera();
};
void vtkComponent::setWidget(QVTKWidget* graphicsWidgetIn)
{
    graphicsWidget = graphicsWidgetIn;
    graphicsWidget->SetRenderWindow(renderWindow);
};
void vtkComponent::addVisualizationDataBoundary(vtkDataSet* input, const float color[3], int LineWidth)
{
    mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
    mapperArray.back()->SetInputData(input);
    actorArray.push_back(vtkSmartPointer<vtkActor>::New());
    actorArray.back()->SetMapper(mapperArray.back());
    actorArray.back()->GetProperty()->SetColor(color[0], color[1], color[2]);
    actorArray.back()->GetProperty()->SetLineWidth(LineWidth);
    renderer->AddActor(actorArray.back());
};
void vtkComponent::addVisualizationDataPlotsAutoY(const std::vector<std::vector<float>>& Xdata,
                                                  const std::vector<std::vector<float>>& Ydata, const float* color,
                                                  float width, std::string xName, std::string yName, int n,
                                                  int Distrtype, std::vector<int> props, float& ymax, float& ymin,
                                                  int startPoint, int endPoint)
{

    xName                   = xName + std::string("\n ");
    yName                   = yName + std::string("\n ");
    int              nPlots = int(Ydata.size());
    std::vector<int> dashed;
    std::vector<int> indexes;
    ymax = -1e20;
    ymin = 1e20;
    if (0 == n)
    {
        indexes.resize(nPlots);
        for (int i     = 0; i < nPlots; i++)
            indexes[i] = i;

        dashed = {1, 2, 10, 11};
        dashed = {};
    }
    else
    {
        if (Distrtype == 0)
        {
            indexes.resize(nPlots);
            for (int i     = 0; i < nPlots; i++)
                indexes[i] = i;

            for (int i = 0; i < nPlots - n; i++)
            {
                int k = rand() % indexes.size();
                indexes.erase(indexes.begin() + k);
            }
        }
        if (Distrtype == 1)
        {
            std::vector<int> indexesTmp;

            indexesTmp.resize(nPlots);
            for (int i        = 0; i < nPlots; i++)
                indexesTmp[i] = i;

            std::sort(indexesTmp.begin(), indexesTmp.end(), [&](int i, int j) { return Xdata[i][0] > Xdata[j][0]; });

            int dn = nPlots / n;

            for (int i = 0; i < n; i++)
                indexes.push_back(indexesTmp[dn * i]);
        }
    }

    //	view->GetScene()->AddItem(chart);

    for (int j = 0; j < indexes.size(); j++)
    {
        int i = indexes[j];

        int endpointLoc = int(Xdata[i].size());

        int flag = 0;
        if (endPoint != -1 && endPoint < endpointLoc)
            endpointLoc = endPoint;

        if (endpointLoc - startPoint < 2)
            flag = 1;
        // continue;
        if (endpointLoc - startPoint == 0)
            continue;

        vtkPlot* line;

        if (flag == 0)
            line = chart->AddPlot(vtkChart::LINE);

        if (flag == 1)
        {
            line = chart->AddPlot(vtkChart::POINTS);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(width);
        }

        tables.push_back(vtkSmartPointer<vtkTable>::New());
        vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("X Axis");
        tables.back()->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName("Y axis");
        tables.back()->AddColumn(arrY);

        /*	if (n != 0)
                {
                int nrow = numPoints / 20;
                tables.back()->SetNumberOfRows(nrow);

                int kk = 0;
                for (int jj = 0; jj < numPoints; jj++)
                {
                if (jj % 20 == 0 && kk<nrow)
                {
                tables.back()->SetValue(kk, 0, Xdata[i][jj]);
                tables.back()->SetValue(kk, 1, Ydata[i][jj]);
                kk++;
                }
                }
                }*/

        //		if (n == 0)
        //		{

        tables.back()->SetNumberOfRows(endpointLoc - startPoint);

        for (int jj = startPoint; jj < endpointLoc; jj++)
        {
            tables.back()->SetValue(jj, 0, Xdata[i][jj]);
            tables.back()->SetValue(jj, 1, Ydata[i][jj]);
            if (Ydata[i][jj] > ymax)
                ymax = Ydata[i][jj];
            if (Ydata[i][jj] < ymin)
                ymin = Ydata[i][jj];
        }
        //	}

        line->SetInputData(tables.back(), 0, 1);
        line->SetColor(color[0], color[1], color[2]);
        line->SetWidth(width);

        for (int k = 0; k < dashed.size(); k++)
        {
            if (dashed[k] == j)
                line->GetPen()->SetLineType(2);
        }
    }

    /*if (0 == n)
    {
            chart->GetAxis(vtkAxis::BOTTOM)->SetUnscaledMaximum(0.85);
            chart->GetAxis(vtkAxis::BOTTOM)->SetUnscaledMinimum(-0.02);
            chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::FIXED);

            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(0.22);
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(-0.02);
            chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
    }*/

    double xmax = chart->GetAxis(vtkAxis::BOTTOM)->GetUnscaledMaximum();
    double xmin = chart->GetAxis(vtkAxis::BOTTOM)->GetUnscaledMinimum();

    // chart->GetAxis(vtkAxis::BOTTOM)->SetNumberOfTicks(props[1]);

    vtkSmartPointer<vtkDoubleArray> ar = vtkSmartPointer<vtkDoubleArray>::New();

    if (props.size() == 0)
    {
        return;
    }

    if (props[1] != -1)
    {
        double dx = (xmax - xmin) / props[1];

        dx = floor(dx * 10) / 10;

        if (xmin < 0 && xmin + dx > 0)
            xmin = 0;

        for (int i = 0; i < props[1] + 1; i++)
            ar->InsertNextValue(double(xmin + double(i) * dx));

        // chart->GetAxis(vtkAxis::BOTTOM)->SetTickPositions(ar);

        chart->GetAxis(vtkAxis::BOTTOM)->SetGridVisible(false);
        // ar->RemoveFirstTuple();
    }
    if (props[2] != -1)
    {
        chart->GetAxis(vtkAxis::LEFT)->SetNumberOfTicks(props[2]);
        chart->GetAxis(vtkAxis::LEFT)->SetGridVisible(false);
    }

    /*chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);


    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);*/

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::AUTO);
    chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::AUTO);

    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);

    // chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::);
    //	ymax = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMaximum();
    //	ymin = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMinimum();
}

void vtkComponent::addVisualizationDataPlots(const std::vector<std::vector<float>>& Xdata,
                                             const std::vector<std::vector<float>>& Ydata, const float* color,
                                             float width, std::string xName, std::string yName, int n, int Distrtype,
                                             std::vector<int> props, float yMin, float yMax)
{
    xName                   = xName + std::string("\n ");
    yName                   = yName + std::string("\n ");
    int              nPlots = int(Ydata.size());
    std::vector<int> dashed;
    std::vector<int> indexes;

    if (0 == n)
    {
        indexes.resize(nPlots);
        for (int i     = 0; i < nPlots; i++)
            indexes[i] = i;

        dashed = {1, 2, 10, 11};
    }
    else
    {
        if (Distrtype == 0)
        {
            indexes.resize(nPlots);
            for (int i     = 0; i < nPlots; i++)
                indexes[i] = i;

            for (int i = 0; i < nPlots - n; i++)
            {
                int k = rand() % indexes.size();
                indexes.erase(indexes.begin() + k);
            }
        }
        if (Distrtype == 1)
        {
            std::vector<int> indexesTmp;

            indexesTmp.resize(nPlots);
            for (int i        = 0; i < nPlots; i++)
                indexesTmp[i] = i;

            std::sort(indexesTmp.begin(), indexesTmp.end(), [&](int i, int j) { return Xdata[i][0] > Xdata[j][0]; });

            int dn = nPlots / n;

            for (int i = 0; i < n; i++)
                indexes.push_back(indexesTmp[dn * i]);
        }
    }

    //	view->GetScene()->AddItem(chart);
    int flag = 0;
    for (int j = 0; j < indexes.size(); j++)
    {
        int i = indexes[j];

        int numPoints = int(Xdata[i].size());
        if (numPoints < 2)
            flag = 1;
        // continue;
        if (numPoints == 0)
            continue;

        vtkPlot* line;

        if (flag == 0)
            line = chart->AddPlot(vtkChart::LINE);

        if (flag == 1)
        {
            line = chart->AddPlot(vtkChart::POINTS);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(width);
        }

        tables.push_back(vtkSmartPointer<vtkTable>::New());
        vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("X Axis");
        tables.back()->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName("Y axis");
        tables.back()->AddColumn(arrY);

        /*	if (n != 0)
        {
        int nrow = numPoints / 20;
        tables.back()->SetNumberOfRows(nrow);

        int kk = 0;
        for (int jj = 0; jj < numPoints; jj++)
        {
        if (jj % 20 == 0 && kk<nrow)
        {
        tables.back()->SetValue(kk, 0, Xdata[i][jj]);
        tables.back()->SetValue(kk, 1, Ydata[i][jj]);
        kk++;
        }
        }
        }*/

        //		if (n == 0)
        //		{
        tables.back()->SetNumberOfRows(numPoints);

        for (int jj = 0; jj < numPoints; jj++)
        {
            tables.back()->SetValue(jj, 0, Xdata[i][jj]);
            tables.back()->SetValue(jj, 1, Ydata[i][jj]);
        }
        //	}

        line->SetInputData(tables.back(), 0, 1);
        line->SetColor(color[0], color[1], color[2]);
        line->SetWidth(width);

        for (int k = 0; k < dashed.size(); k++)
        {
            if (dashed[k] == j)
                line->GetPen()->SetLineType(2);
        }
    }

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());

    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);

    // chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetBold(5);
    // chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetBold(10);

    chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(yMax);
    chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(yMin);
    chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);

    vtkSmartPointer<vtkStringArray> labelArray = vtkSmartPointer<vtkStringArray>::New();
}

void vtkComponent::addVisualizationDataFillRect(const std::vector<float>& Xdata, const std::vector<float>& Ydata,
                                                const float* color, std::string legend, float opacity){
    /*CheckType(1);
    clear();


    XdataCurrent = Xdata;
    YdataCurrent.clear();
    YdataCurrent.push_back(Ydata);

    int nPlots = 1;
    ///vtkSmartPointer<vtkContext2D> painter = vtkSmartPointer<vtkContext2D>::New();


    vtkPlot *line = chart->AddPlot(vtkChart::LINE);
    //painter->DrawRect(100, 50, 200, 100);
    APIDiagram * diag = new APIDiagram();
    chart->AddItem(diag);
    //chart->Paint(painter);
    //contextActor->GetScene()->AddItem(view);
    tables.push_back(vtkSmartPointer<vtkTable>::New());
    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("X Axis");
    tables.back()->AddColumn(arrX);

    vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
    arrY->SetName("Y axis");
    tables.back()->AddColumn(arrY);
    arrY->SetName(legend.c_str());

    int numPoints = int(Xdata.size());

    tables.back()->SetNumberOfRows(numPoints);

    for (int j = 0; j < numPoints; j++)
    {
            tables.back()->SetValue(j, 0, Xdata[j]);
            tables.back()->SetValue(j, 1, Ydata[j]);
    }

    line->SetInputData(tables.back(), 0, 1);
    line->SetColor(color[0], color[1], color[2]);
    line->GetPen()->SetOpacityF(opacity);
    line->SetWidth(140);
    contextActor->GetScene()->AddItem(line);

    renderer->AddActor(contextActor);

    //	view->GetScene()->AddItem(chart);
    /*float ymax = -1e20;
    float ymin = 1e20;
    vtkPlot *line;

    for (int i = 0; i < nPlots; i++)
    {
            line = chart->AddPlot(vtkChart::FILL_RECT);


            int numPoints = int(Xdata.size());
            if (numPoints < 2)
                    continue;

            tables.push_back(vtkSmartPointer<vtkTable>::New());
            vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
            arrX->SetName("X Axis");
            tables.back()->AddColumn(arrX);

            vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
            arrY->SetName("Y axis");
            tables.back()->AddColumn(arrY);
            arrY->SetName(legend.c_str());


            tables.back()->SetNumberOfRows(numPoints);

            for (int j = 0; j < numPoints; j++)
            {
                    tables.back()->SetValue(j, 0, Xdata[j]);
                    tables.back()->SetValue(j, 1, Ydata[j]);

                    if (Ydata[j] < ymin)
                            ymin = Ydata[j];

                    if (Ydata[j] > ymax)
                            ymax = Ydata[j];
            }

            line->SetInputData(tables.back(), 0, 1);
            line->SetColor(color[0], color[1], color[2]);
            line->GetPen()->SetOpacityF(opacity);
            //line->set
    };*/
};
void vtkComponent::addVisualizationDataPlot(const std::vector<float>& Xdata, const std::vector<float>& Ydata,
                                            const float* color, float width, std::string xName, std::string yName,
                                            int flag, int style, int Markerwidth, std::string legend, float opacity,
                                            bool flagFixAxis)
{
    if (xName.size() != 0)
    {
        xName = xName + std::string("\n ");
        chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    };

    if (yName.size() != 0)
    {
        yName = yName + std::string("\n ");
        chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    };

    XdataCurrent = Xdata;
    YdataCurrent.clear();
    YdataCurrent.push_back(Ydata);

    int nPlots = 1;

    //	view->GetScene()->AddItem(chart);
    float ymaxCur = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMaximum();
    float yminCur = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMinimum();
    float ymax    = -1e20;
    float ymin    = 1e20;
    for (int i = 0; i < nPlots; i++)
    {
        int numPoints = int(Xdata.size());
        if (numPoints < 2)
            flag = 1;
        // continue;
        if (numPoints == 0)
            continue;

        vtkPlot* line;

        if (flag == 0)
        {
            line = chart->AddPlot(vtkChart::LINE);
            line->SetWidth(width);

            if (Markerwidth != 0)
            {
                vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
                vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(Markerwidth);
            }
        }
        if (flag == 1)
        {
            line = chart->AddPlot(vtkChart::POINTS);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(width);
        }
        if (flag == 2)
        {
            line = chart->AddPlot(vtkChart::BAR);
            //	vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            //	vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(width);
        }

        tables.push_back(vtkSmartPointer<vtkTable>::New());
        vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("X Axis");
        tables.back()->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName("Y axis");
        tables.back()->AddColumn(arrY);
        arrY->SetName(legend.c_str());

        tables.back()->SetNumberOfRows(numPoints);

        for (int j = 0; j < numPoints; j++)
        {
            tables.back()->SetValue(j, 0, Xdata[j]);
            tables.back()->SetValue(j, 1, Ydata[j]);

            if (Ydata[j] < ymin)
                ymin = Ydata[j];

            if (Ydata[j] > ymax)
                ymax = Ydata[j];
        }

        line->SetInputData(tables.back(), 0, 1);
        line->SetColor(color[0], color[1], color[2]);
        // line->SetWidth(width);
        // line->GetPen()->SetLineType(style);
        // line->GetPen()->SetOpacityF(opacity);
    }

    float k1, k2;

    if (ymaxCur > ymax)
    {
        ymax = ymaxCur;
        k1   = 1;
    }
    else
        k1 = 1.1;

    if (yminCur < ymin)
    {
        ymin = yminCur;
        k2   = 1;
    }
    else
        k2 = 1.1;

    if (flagFixAxis)
    {
        if (ymax > 0)
        {
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(k1 * ymax);
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(std::min(float(0.0), float(k2 * ymin)));
            chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
        }

        if (ymax < 0)
        {
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(std::max(float(0.0), float(k1 * ymax)));
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(k2 * ymin);
            chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
        }
    }
    else
    {
        chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::AUTO);
        chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::AUTO);
    };

    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);
}
template void vtkComponent::addVisualizationDataPlot<float>(const std::vector<float>& Xdata,
                                                            const std::vector<float>& Ydata, float width,
                                                            std::string xName, std::string yName, std::string Name,
                                                            float& ymax, float& ymin, bool flagFixAxis, int startPoint,
                                                            int endPoint);
// template void vtkComponent::addVisualizationDataPlot<double>(const std::vector<double>& Xdata, const
// std::vector<double>& Ydata, float width, std::string xName, std::string yName, std::string Name, float& ymax, float&
// ymin, bool flagFixAxis);

template <class Datatype>
void vtkComponent::addVisualizationDataPlot(const std::vector<Datatype>& Xdata, const std::vector<Datatype>& Ydata,
                                            float width, std::string xName, std::string yName, std::string Name,
                                            float& ymax, float& ymin, bool flagFixAxis, int startPoint, int endPoint)
{

    xName = xName + std::string("\n ");
    yName = yName + std::string("\n ");

    XdataCurrent = Xdata;
    YdataCurrent.clear();

    YdataCurrent.push_back(Ydata);

    int nPlots = 1;

    //	view->GetScene()->AddItem(chart);
    //	float ymaxCur = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMaximum();
    //	float yminCur = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMinimum();
    float ymaxCur = 0;
    float yminCur = 0;
    // ymax = -1e20;
    // ymin = 1e20;

    int flag = 0;
    for (int i = 0; i < nPlots; i++)
    {
        int endpointLoc = int(std::min(Xdata.size(), Ydata.size()));

        tables.push_back(vtkSmartPointer<vtkTable>::New());
        vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("X Axis");
        tables.back()->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(Name.c_str());
        tables.back()->AddColumn(arrY);

        if (endPoint != -1)
            endpointLoc = endPoint;

        int numPoints = endpointLoc - startPoint;

        if (numPoints < 2)
            flag = 1;
        // continue;
        if (numPoints == 0)
            continue;

        vtkPlot* line;

        if (flag == 0)
            line = chart->AddPlot(vtkChart::LINE);
        if (flag == 1)
        {
            line = chart->AddPlot(vtkChart::POINTS);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(width);
        }

        tables.back()->SetNumberOfRows(numPoints);

        for (int j = startPoint; j < numPoints; j++)
        {
            tables.back()->SetValue(j, 0, Xdata[j]);
            tables.back()->SetValue(j, 1, Ydata[j]);

            if (Ydata[j] < ymin)
                ymin = Ydata[j];

            if (Ydata[j] > ymax)
                ymax = Ydata[j];
        }

        line->SetInputData(tables.back(), 0, 1);
        line->SetWidth(width);
    }

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());

    float k1, k2;

    if (ymaxCur > ymax)
    {
        ymax = ymaxCur;
        k1   = 1;
    }
    else
        k1 = 1.1;

    if (yminCur < ymin)
    {
        ymin = yminCur;
        k2   = 1;
    }
    else
        k2 = 1.1;

    if (flagFixAxis)
    {
        if (ymax > 0)
        {
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(k1 * ymax);
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(std::min(float(0.0), float(k2 * ymin)));
            chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
        }

        if (ymax < 0)
        {
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(std::max(float(0.0), float(k1 * ymax)));
            chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(k2 * ymin);
            chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
        }
    }
    else
    {
        chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::AUTO);
        chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::AUTO);
    };

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);
}

void vtkComponent::addVisualizationDataPlot(const std::vector<float>& Xdata, const std::vector<float>& Ydata,
                                            float width, std::string xName, std::string yName, std::string Name)
{

    xName = xName + std::string("\n ");
    yName = yName + std::string("\n ");

    XdataCurrent = Xdata;
    YdataCurrent.clear();

    YdataCurrent.push_back(Ydata);

    int nPlots = 1;
    int flag   = 0;
    for (int i = 0; i < nPlots; i++)
    {
        int numPoints = int(std::min(Xdata.size(), Ydata.size()));
        if (numPoints < 2)
            flag = 1;
        if (numPoints == 0)
            continue;
        //

        vtkPlot* line;

        if (flag == 0)
            line = chart->AddPlot(vtkChart::LINE);

        if (flag == 1)
        {
            line = chart->AddPlot(vtkChart::POINTS);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
            vtkPlotPoints::SafeDownCast(line)->SetMarkerSize(width);
        }

        tables.push_back(vtkSmartPointer<vtkTable>::New());
        vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
        arrX->SetName("X Axis");
        tables.back()->AddColumn(arrX);

        vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(Name.c_str());
        tables.back()->AddColumn(arrY);

        tables.back()->SetNumberOfRows(numPoints);

        for (int j = 0; j < numPoints; j++)
        {
            tables.back()->SetValue(j, 0, Xdata[j]);
            tables.back()->SetValue(j, 1, Ydata[j]);
        }

        line->SetInputData(tables.back(), 0, 1);
        line->SetWidth(width);
    }

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::AUTO);
    chart->GetAxis(vtkAxis::BOTTOM)->SetBehavior(vtkAxis::AUTO);
}

void vtkComponent::addVisualizationDataPlot(const std::vector<float>&              Xdata,
                                            const std::vector<std::vector<float>>& Ydata,
                                            std::vector<const float*> colors, std::vector<float> width,
                                            std::string xName, std::string yName, std::vector<std::string> legends)
{
    if (!Xdata.size())
        return;
    xName = xName + std::string("\n ");
    yName = yName + std::string("\n ");

    XdataCurrent = Xdata;
    YdataCurrent = Ydata;

    int nPlots = int(Ydata.size());

    float ymax = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMaximum();
    float ymin = chart->GetAxis(vtkAxis::LEFT)->GetUnscaledMinimum();

    tables.push_back(vtkSmartPointer<vtkTable>::New());

    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("X Axis");
    tables.back()->AddColumn(arrX);

    for (int i = 0; i < nPlots; i++)
    {
        vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
        arrY->SetName(legends[i].c_str());
        tables.back()->AddColumn(arrY);
    };

    int numPoints = int(Xdata.size());
    tables.back()->SetNumberOfRows(numPoints);

    for (int i = 0; i < numPoints; i++)
    {
        tables.back()->SetValue(i, 0, Xdata[i]);
        for (int j = 0; j < nPlots; j++)
        {
            tables.back()->SetValue(i, j + 1, Ydata[j][i]);

            if (Ydata[j][i] < ymin)
                ymin = Ydata[j][i];

            if (Ydata[j][i] > ymax)
                ymax = Ydata[j][i];
        }
    }

    int flag = 0;
    if (numPoints < 2)
        flag = 1;
    // continue;
    if (numPoints == 0)
        return;

    vtkPlot* line;

    if (flag == 0)
        line = chart->AddPlot(vtkChart::LINE);

    if (flag == 1)
    {
        line = chart->AddPlot(vtkChart::POINTS);
        vtkPlotPoints::SafeDownCast(line)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
    }

    for (int i = 0; i < nPlots; i++)
    {
        line->SetInputData(tables.back(), 0, i + 1);
        line->SetColor(colors[i][0], colors[i][1], colors[i][2]);
        line->SetWidth(width[i]);
    }

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());

    if (ymax > 0)
    {
        chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(1.1 * ymax);
        chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(std::min(float(0.0), float(1.1 * ymin)));
        chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
    }

    if (ymax < 0)
    {
        chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMaximum(std::max(float(0.0), float(1.1 * ymax)));
        chart->GetAxis(vtkAxis::LEFT)->SetUnscaledMinimum(1.1 * ymin);
        chart->GetAxis(vtkAxis::LEFT)->SetBehavior(vtkAxis::FIXED);
    }

    chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xName.c_str());
    chart->GetAxis(vtkAxis::BOTTOM)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::BOTTOM)->GetLabelProperties()->SetFontSize(labelFont);

    chart->GetAxis(vtkAxis::LEFT)->SetTitle(yName.c_str());
    chart->GetAxis(vtkAxis::LEFT)->GetTitleProperties()->SetFontSize(titleFont);
    chart->GetAxis(vtkAxis::LEFT)->GetLabelProperties()->SetFontSize(labelFont);
}
void vtkComponent::addVisualizationDataCloud(vtkDataSet* input, float size, const float* color)
{
    mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
    mapperArray.back()->SetInputData(input);
    actorArray.push_back(vtkSmartPointer<vtkActor>::New());
    actorArray.back()->SetMapper(mapperArray.back());
    actorArray.back()->GetProperty()->SetPointSize(size);
    actorArray.back()->GetProperty()->SetColor(color[0], color[1], color[2]);
    renderer->AddActor(actorArray.back());
};
void vtkComponent::addVisualizationDataMesh(vtkDataSet* input, int flagEdge)
{
    mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
    mapperArray.back()->SetInputData(input);
    actorArray.push_back(vtkSmartPointer<vtkActor>::New());
    actorArray.back()->SetMapper(mapperArray.back());
    actorArray.back()->GetProperty()->SetColor(background);
    if (flagEdge == 0)
        actorArray.back()->GetProperty()->EdgeVisibilityOn();
    renderer->AddActor(actorArray.back());
};
void vtkComponent::addVisualizationDataMesh(vtkDataSet* input, vtkSmartPointer<vtkFloatArray> vtkData, double minVal,
                                            double maxVal, std::string title)
{
    flagScalarBar = 1;
    mapperArray.push_back(vtkSmartPointer<vtkDataSetMapper>::New());
    mapperArray.back()->SetInputData(input);
    actorArray.push_back(vtkSmartPointer<vtkActor>::New());
    actorArray.back()->SetMapper(mapperArray.back());
    actorArray.back()->GetProperty()->SetColor(0.0, 1.0, 0.0);

    renderer->AddActor(actorArray.back());

    //	actorArray.back()->GetProperty()->SetDiffuseColor()

    vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
    poly->DeepCopy(input);
    poly->GetPointData()->SetScalars(vtkData);

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

    mapper->SetInputData(poly);
    mapper->ScalarVisibilityOn();
    mapper->SetScalarModeToUsePointData();
    mapper->SetColorModeToMapScalars();

    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle(title.c_str());
    scalarBar->GetTitleTextProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1], colors::blackcolor[2]);
    scalarBar->SetNumberOfLabels(10);
    scalarBar->GetLabelTextProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1], colors::blackcolor[2]);

    vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
    hueLut->SetTableRange(minVal, maxVal);
    hueLut->Build();

    mapper->SetLookupTable(hueLut);
    scalarBar->SetLookupTable(hueLut);

    renderer->AddActor2D(scalarBar);
};

void vtkComponent::removeVisualizationData(int n)
{
    renderer->RemoveActor(actorArray[n]);
    mapperArray.erase(mapperArray.begin() + n);
    actorArray.erase(actorArray.begin() + n);
};
void vtkComponent::clear()
{
    for (int i = 0; i < actorArray.size(); i++)
    {
        renderer->RemoveActor(actorArray[i]);
    };

    if (flagScalarBar)
        renderer->RemoveActor2D(scalarBar);

    if (flagcubeAxesActor)
        renderer->RemoveActor(cubeAxesActor);

    flagScalarBar     = 0;
    flagcubeAxesActor = 0;
    //	view->GetScene()->RemoveItem(chart);
    actorArray.clear();
    mapperArray.clear();

    if (chart->GetNumberOfPlots() != 0)
    {
        view->GetScene()->RemoveItem(chart);
        //	chart->ClearPlots();
    }
    if (chart3d.GetPointer()->GetNumberOfItems() != 0)
    {
        view->GetScene()->RemoveItem(chart3d);
        //	chart->ClearPlots();
    }
    tables.clear();
};
void vtkComponent::setAxis()
{
    flagcubeAxesActor = 1;

    double bounds[6];
    double boundsMax[6] = {0, 0, 0, 0, 0, 0};

    for (int i = 0; i < mapperArray.size(); i++)
    {
        mapperArray[i]->GetBounds(bounds);
        for (int j = 0; j < 6; j++)
        {
            if (bounds[j] > boundsMax[j])
                boundsMax[j] = bounds[j];
        };
    }

    for (int i = 0; i < 3; i++)
    {
        if (boundsMax[2 * i] < 0)
            boundsMax[2 * i] = boundsMax[2 * i] * 1.1;
        else
            boundsMax[2 * i] = boundsMax[2 * i] * 0.9;

        if (boundsMax[2 * i + 1] < 0)
            boundsMax[2 * i + 1] = boundsMax[2 * i + 1] * 0.9;
        else
            boundsMax[2 * i + 1] = boundsMax[2 * i + 1] * 1.1;
    }

    cubeAxesActor->SetBounds(boundsMax);
    cubeAxesActor->SetCamera(renderer->GetActiveCamera());
    cubeAxesActor->GetTitleTextProperty(0)->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);
    cubeAxesActor->GetTitleTextProperty(0)->SetFontSize(20.0);
    cubeAxesActor->GetLabelTextProperty(0)->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);
    cubeAxesActor->GetLabelTextProperty(0)->SetFontSize(60.0);

    cubeAxesActor->GetTitleTextProperty(1)->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);
    cubeAxesActor->GetTitleTextProperty(1)->SetFontSize(30.0);
    cubeAxesActor->GetLabelTextProperty(1)->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);

    cubeAxesActor->GetTitleTextProperty(2)->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);
    cubeAxesActor->GetLabelTextProperty(2)->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);

    cubeAxesActor->GetXAxesLinesProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);
    cubeAxesActor->GetYAxesLinesProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);
    cubeAxesActor->GetZAxesLinesProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1],
                                                     colors::blackcolor[2]);

    cubeAxesActor->SetFlyModeToStaticTriad();

    cubeAxesActor->XAxisMinorTickVisibilityOff();
    cubeAxesActor->YAxisMinorTickVisibilityOff();
    cubeAxesActor->ZAxisMinorTickVisibilityOff();

    renderer->AddActor(cubeAxesActor);
};

void vtkComponent::AddSurfacePlot()
{
    tables.push_back(vtkSmartPointer<vtkTable>::New());
    float numPoints = 70;
    float inc       = 9.424778 / (numPoints - 1);
    for (float i = 0; i < numPoints; ++i)
    {
        vtkSmartPointer<vtkFloatArray> arr = vtkSmartPointer<vtkFloatArray>::New();
        tables.back()->AddColumn(arr.GetPointer());
    }
    tables.back()->SetNumberOfRows(numPoints);
    for (float i = 0; i < numPoints; ++i)
    {
        float x = i * inc;
        for (float j = 0; j < numPoints; ++j)
        {
            float y = j * inc;
            tables.back()->SetValue(i, j, std::sin(sqrt(x * x + y * y)));
        }
    }

    vtkSmartPointer<vtkPlotSurface> plot = vtkSmartPointer<vtkPlotSurface>::New();

    chart3d->SetGeometry(vtkRectf(200.0, 200.0, 300, 300));
    plot->SetXRange(0, 10.0);
    plot->SetYRange(0, 10.0);

    //	plot->SetInputData(tables.back().GetPointer());
    //	plot->
    chart3d.GetPointer()->AddPlot(plot.GetPointer());
};