#include "VisualizationFunctions.h"
#include "VTKInclude.h"
#include "colors.h"
#include "vtkComponent.h"

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

const double PI                  = 3.141592653;
const double LIGHT_VELOCITY      = 299792458.0;
const double ELECTRON_CHARGE     = -1.60217656e-19;
const double ELECTRON_MASS       = 9.10938215e-31;
const double PROTON_MASS         = 1.67262178e-27;
const double VACUUM_PERMITTIVITY = 8.85418781e-12;
const double VACUUM_PERMEABILITY = 1.25663706e-6;

void SimpleShowPoints(vtkComponent* input, const std::vector<float>& dataX, const std::vector<float>& dataY,
                      const std::vector<std::vector<float>>& dataXLine,
                      const std::vector<std::vector<float>>& dataYLine)
{
    input->CheckType(1);
    input->clear();
    input->addVisualizationDataPlot(dataX, dataY, colors::blackcolor, 4.0, "X", "Y", 1, 1, 0, "distribution", 1.0,
                                    false);

    for (int i = 0; i < dataYLine.size(); i++)
    {
        input->addVisualizationDataPlot(dataXLine[i], dataYLine[i], colors::colors[i], 2, "X", "Y", 0, 1, 0, "Line",
                                        1.0, false);
        input->addVisualizationDataPlot(dataXLine[i], dataYLine[i], colors::colors[i], 2, "X", "Y", 0, 1, 0, "Line",
                                        1.0, false);
    }
};

void SimpleShowChart(vtkComponent* input, const std::vector<float>& data)
{
    input->CheckType(1);
    input->clear();
    std::vector<float> XChart;
    std::vector<float> YChart;
    int                nCharts = 50;
    if (data.size() == 0)
        return;

    volatile float wmin = data[0];
    volatile float wmax = data[0];

    for (int i = 0; i < data.size(); i++)
    {
        if (data[i] > wmax)
            wmax = data[i];
        if (data[i] < wmin)
            wmin = data[i];
    };

    if (std::abs(wmax - wmin) < 1e-7)
        return;

    if (wmin < 0)
        wmin = wmin * 1.1;
    else
        wmin = wmin * 0.9;

    if (wmax > 0)
        wmax = wmax * 1.1;
    else
        wmax = wmax * 0.9;

    float d = (wmax - wmin) / nCharts;
    XChart.resize(nCharts + 1);
    YChart.resize(nCharts + 1);

    XChart[0]     = wmin;
    XChart.back() = wmax;

    for (int i = 1; i < nCharts; i++)
        XChart[i] = wmin + i * d;

    for (int i = 1; i < nCharts; i++)
    {
        for (int j = 0; j < data.size(); j++)
        {
            if (data[j] > XChart[i] && data[j] < XChart[i + 1])
                YChart[i]++;
        }
    }

    input->addVisualizationDataPlot(XChart, YChart, colors::blackcolor, 5, "", "", 2, 1, 0, "distribution", 1.0, false);
}
void ShowLinacResultsTraces(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                            int flag, int Distrtype, std::vector<int> props, int numberOfTraces)
{
   /* float ymax;
    float ymin;
    input->CheckType(1);
    input->clear();
    std::string                     yName, xName;
    std::vector<std::vector<double>> tmp;
    std::vector<double>              Zav(data->data[1][0].size());
    std::vector<double>              beta = data->dataAdd[4];
    float                           gamma;
    float                           en;
    switch (flag)
    {
    case 0:
        yName = "Energy, Ev";
        xName = "Z, m";
        tmp   = data->data[3];

        en = LIGHT_VELOCITY * LIGHT_VELOCITY * data->mass;

        for (int i = 0; i < data->data[0].size(); i++)
        {
            std::vector<float> Xtmp;
            for (int j = 0; j < data->data[0][i].size(); j++)
            {
                float V   = tmp[i][j] * LIGHT_VELOCITY;
                tmp[i][j] = data->mass * V * V / (std::abs(ELECTRON_CHARGE) * 2);
            }
        }

        // input->addVisualizationDataPlots(data->data[0], tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
        // Distrtype, props);
        input->addVisualizationDataPlotsAutoY(data->data[0], tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
                                              Distrtype, props, ymax, ymin);
        break;
    case 1:
        yName = "Phase, rad";
        xName = "T, ns";

        tmp = data->data[0];

        for (int k = 0; k < data->data[0][0].size(); k++)
        {
            Zav[k] = data->data[0][0][k];
            /*int s = 0;
            for (int i = 0; i < data->TimeArray.size(); i++)
            {
                    if (k < data->data[0][i].size())
                    {
                            Zav[k] = Zav[k] + data->data[0][i][k];
                            s++;
                    }
            }
            Zav[k] = Zav[k] / s;

            for (int i = 0; i < data->TimeArray.size(); i++)
            {
                tmp[i][k] = tmp[i][k] - Zav[k];
                tmp[i][k] = 2 * PI * tmp[i][k] / (beta[k] * data->lambda);
            }
        }

        input->addVisualizationDataPlots(data->TimeArray, tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
                                         Distrtype, props, -3.5, 3.5);
        break;
    case 2:
        yName = "Rx, m";
        xName = "Z, m";
        tmp.resize(data->data[1].size());

        for (int i = 0; i < data->data[1].size(); i++)
        {
            tmp[i].resize(data->data[1][i].size());
            for (int j = 0; j < data->data[1][i].size(); j++)
                tmp[i][j] = sqrt(data->data[1][i][j]);
        }

        input->addVisualizationDataPlotsAutoY(data->data[0], tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
                                              Distrtype, props, ymax, ymin);
 //       input->addVisualizationDataPlot(data->dataAdd[6], data->dataAdd[7], 2, xName, yName, std::string("Channel"),
//                                        ymax, ymin);

        break;
    case 3:

        yName = "Ry, m";
        xName = "Z, m";
        tmp.resize(data->data[2].size());
        for (int i = 0; i < data->data[2].size(); i++)
        {
            tmp[i].resize(data->data[2][i].size());
            for (int j = 0; j < data->data[2][i].size(); j++)
                tmp[i][j] = sqrt(data->data[2][i][j]);
        }

        input->addVisualizationDataPlotsAutoY(data->data[0], tmp, colors::blackcolor, 1, xName, yName, numberOfTraces,
                                              Distrtype, props, ymax, ymin);
        input->addVisualizationDataPlot(data->dataAdd[6], data->dataAdd[7], 2, xName, yName, std::string("Channel"),
                                        ymax, ymin);

        break;

    case 5:
    {

        yName = "X, m";
        xName = "dX, rad";

        std::vector<std::vector<float>> Xar(numberOfTraces);
        std::vector<std::vector<float>> Yar1(numberOfTraces);
        std::vector<std::vector<float>> Yar2(numberOfTraces);

        for (int k = 0; k < numberOfTraces; k++)
        {

            int    tr = rand() % data->dataAdd[8].size();
            double V  = data->data[3][0].back();

            double em = data->emittanceX * 1e-6;

            double betaX  = data->dataAdd[8][tr] / (em);
            double AlphaX = -data->dataAdd[9][tr] / em;
            double gammaX = data->dataAdd[10][tr] / (em);

            double xM  = sqrt(data->dataAdd[8][tr]);
            double dxM = sqrt(data->dataAdd[10][tr]);
            int    n   = 500;
            double dx  = 2 * xM / n;
            double xc  = -1.5 * xM;
            double dxC;
            while (1)
            {
                if (std::abs(xM - std::abs(xc)) < 0.1 * xM)
                    dxC = 0.04 * dx;
                else
                    dxC = dx;

                Xar[k].push_back(xc);

                double b = 2 * AlphaX * xc;
                double d = b * b - 4 * betaX * (gammaX * xc * xc - em);

                Yar1[k].push_back((1.0 / V) * (-b + sqrt(d)) / (2 * betaX));
                Yar2[k].push_back((1.0 / V) * (-b - sqrt(d)) / (2 * betaX));

                xc = xc + dxC;

                if (xc > 1.5 * xM)
                    break;
            }
        }
        input->addVisualizationDataPlotsAutoY(Xar, Yar1, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(Xar, Yar2, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
    }
    break;

    case 6:
    {
        yName = "Y, m";
        xName = "dY, rad";

        std::vector<std::vector<float>> Xar(numberOfTraces);
        std::vector<std::vector<float>> Yar1(numberOfTraces);
        std::vector<std::vector<float>> Yar2(numberOfTraces);

        for (int k = 0; k < numberOfTraces; k++)
        {

            int    tr = rand() % data->dataAdd[8].size();
            double V  = data->data[3][0].back();

            double em = data->emittanceY * 1e-6;

            double betaX  = data->dataAdd[11][tr] / (em);
            double AlphaX = -data->dataAdd[12][tr] / em;
            double gammaX = data->dataAdd[13][tr] / (em);

            double xM  = sqrt(data->dataAdd[11][tr]);
            double dxM = sqrt(data->dataAdd[13][tr]);
            int    n   = 500;
            double dx  = 2 * xM / n;
            double xc  = -1.5 * xM;
            double dxC;
            while (1)
            {
                if (std::abs(xM - std::abs(xc)) < 0.1 * xM)
                    dxC = 0.04 * dx;
                else
                    dxC = dx;

                Xar[k].push_back(xc);

                double b = 2 * AlphaX * xc;
                double d = b * b - 4 * betaX * (gammaX * xc * xc - em);

                Yar1[k].push_back((1.0 / V) * (-b + sqrt(d)) / (2 * betaX));
                Yar2[k].push_back((1.0 / V) * (-b - sqrt(d)) / (2 * betaX));

                xc = xc + dxC;

                if (xc > 1.5 * xM)
                    break;
            }
        }
        input->addVisualizationDataPlotsAutoY(Xar, Yar1, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(Xar, Yar2, colors::blackcolor, 1, xName, yName, numberOfTraces, Distrtype,
                                              props, ymax, ymin);
    }
    break;
    }*/
};

void ShowLinacResultsChar(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                          int flag, int Distrtype, std::vector<int> props, int numberOfTraces)
{
/*    input->CheckType(1);
    input->clear();
    std::string yName, xName;
    float       ymax;
    float       ymin;
    float       gamma;
    float       en;
    switch (flag)
    {
    case 7:
        yName = "R chan";
        xName = "T, ns";
//        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[0], 2, xName, yName, std::string("R channel"),
 //                                       ymax, ymin);
        yName = "R beam";
        xName = "T, ns";
//        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[1], 2, xName, yName, std::string("R beam"),
  //                                      ymax, ymin);
        break;
    case 8:
        yName = "Transmission";
        xName = "T, ns";
//        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[2], 2, xName, yName,
//                                        std::string("Transmission"), ymax, ymin);
        yName = "Acceleration";
        xName = "T, ns";
//        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[3], 2, xName, yName,
 //                                       std::string("Acceleration"), ymax, ymin);
        break;

    case 9:
        yName = "Defocusing";
        xName = "T, ns";
        input->addVisualizationDataPlot(data->TimeArrayAdd, data->dataAdd[5], 2, xName, yName,
                                        std::string("Defocusing"), ymax, ymin);
        break;

    case 10:
        yName                    = "Energy, Ev";
        xName                    = "Z, m";
        std::vector<double> tmp   = data->data[3][0];
        std::vector<double> tmpT  = data->TimeArray[0];
        std::vector<double> tmpT1 = data->dataAdd[14];

        en = LIGHT_VELOCITY * LIGHT_VELOCITY * data->mass;

        for (int i = 0; i < tmp.size(); i++)
            tmp[i] = tmp[i] * LIGHT_VELOCITY;

        float omega = 2 * PI * 108e6;

        for (int i = 0; i < tmp.size(); i++)
        {
            tmpT[i] = tmpT[i] * 1e-9 * omega;
            while (tmpT[i] > PI)
                tmpT[i] = tmpT[i] - 2 * PI;
        }

        std::vector<float> cells;
        for (int i = 0; i < data->dataAdd[14].size(); i++)
        {
            cells.push_back(i);

            while (tmpT1[i] > PI / 2)
                tmpT1[i] = tmpT1[i] - PI;
        }
        // input->addVisualizationDataPlot(data->data[0][0], tmpT, 2, xName, yName, std::string("Velocity"));
        //input->addVisualizationDataPlot(cells, tmpT1, 2, xName, yName, std::string("Velocity"), ymax, ymin);

        // input->addVisualizationDataPlot(data->TimeArray[0], tmp, 2, xName, yName, std::string("Velocity"));
    }*/
};

void ShowRFQRFQCavityParametersPlots(vtkComponent* input, ModelInterface* currentModel, int flag)
{
    input->CheckType(1);
    input->clear();
    std::vector<double> data;
    float               ymax;
    float               ymin;
    // data = currentModel->GetRFQCavityParameters(flag);

    std::vector<float> dataX(data.size());
    std::vector<float> dataY(data.size());

    for (int i = 0; i < data.size(); i++)
    {
        dataX[i] = i + 1;
        dataY[i] = data[i];
    };

    std::string yName;

    if (flag == 0)
        yName = "Cells lengths, m";

    if (flag == 1)
        yName = "Regular part channel minimal radii, m";

    if (flag == 2)
        yName = "Mather profile radii, m";

    if (flag == 3)
        yName = "Channel minimal radii, m";

    if (flag == 3)
        yName = "Average Energy, eV";

    std::string xName  = "Cell number";
    std::string contor = "Parameter";

    input->addVisualizationDataPlot(dataX, dataY, 2, xName, yName, contor, ymax, ymin);

    if (flag == 0)
    {
        std::vector<float> dataY1 = dataY;
        dataY1[0]                 = 0;
        for (int i = 0; i < data.size() - 1; i++)
            dataY1[i + 1] = dataY1[i + 1] + dataY1[i];

        for (int i = 0; i < dataY1.size(); i++)
            dataY1[i] = dataY1[i] / 100;

        std::string contor1 = "Length";

        input->addVisualizationDataPlot(dataX, dataY1, 1, xName, yName, contor1, ymax, ymin);
    };
}

void ShowLinacControls(vtkComponent* input, ModelInterface* currentModel, int flag)
{
    input->CheckType(1);
    input->clear();
    float               ymax;
    float               ymin;
    std::vector<double> data;

    //	data = currentModel->GetControl(flag);

    std::vector<float> dataX(data.size());
    std::vector<float> dataY(data.size());

    for (int i = 0; i < data.size(); i++)
    {
        dataX[i] = i;
        dataY[i] = data[i];
    };

    std::string yName  = "control";
    std::string xName  = "cell number";
    std::string contor = "control";

    input->addVisualizationDataPlot(dataX, dataY, 2, xName, yName, contor, ymax, ymin);
};

void ShowBoundaries(vtkComponent* input, ModelInterface* currentModel, const std::vector<int>& list1,
                    const float* color1, const std::vector<int>& list2, const float* color2)
{
    // �������� ��� ��������� �������. ��������� ������� �� ������ highLigthNumbers
    input->CheckType(0);
    input->clear();

    for (int i = 0; i < list1.size(); i++)
        input->addVisualizationDataBoundary(currentModel->GetBoundaryVTKUnstructuredGrid(list1[i]), color1);

    for (int i = 0; i < list2.size(); i++)
        input->addVisualizationDataBoundary(currentModel->GetBoundaryVTKUnstructuredGrid(list2[i]), color2);

    // input->setAxis();
};
;

void HighLigthBoundary(vtkComponent* input, ModelInterface* currentModel, int n, const float* color)
{
    // �������� ��� ��������� ������� �� ����� ��� ������������ �����. ��������� ������� n
    // input->removeVisualizationData(n);
    input->addVisualizationDataBoundary(currentModel->GetBoundaryVTKUnstructuredGrid(n), color);
};

void ShowMesh(vtkComponent* input, ModelInterface* currentModel, const float* color)
{
    input->CheckType(0);
    input->clear();
    input->addVisualizationDataBoundary(currentModel->GetMeshBoundaryVTKGrid(), color, 3);
    input->setAxis();

    input->addVisualizationDataMesh(currentModel->GetVTKGrid());
    vtkSmartPointer<vtkPolyData> vertex = currentModel->GetVTKBoundaryPoints();

    input->addVisualizationDataCloud(vertex, 2.0, colors::redcolor);
};
void ShowParticlesCloud(vtkComponent* input, ModelInterface* currentModel, int flowNumber, int precisionType,
                        int problemType)
{
    input->CheckType(0);
    input->clear();
    std::vector<int>                arraySize;
    int                             elemSize;
    std::vector<std::vector<void*>> pointArray;
    currentModel->GetParticlesCloud(0, flowNumber, pointArray, arraySize, elemSize);

    vtkSmartPointer<vtkPoints>    points   = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType                     pid;
    switch (precisionType)
    {
    case 0:
    {
        std::vector<double> tmp = {0, 0, 0};

        for (int thread = 0; thread < pointArray.size(); thread++)
        {
            for (int i = 0; i < arraySize[thread]; i++)
            {
                for (int k = 0; k < pointArray[thread].size(); k++)
                    tmp[k] = *((double*)((char*)pointArray[thread][k] + i * elemSize));

                if (problemType == 2)
                    pid = points->InsertNextPoint(tmp[0], tmp[1], 0);
                else
                    pid = points->InsertNextPoint(tmp[0], tmp[1], tmp[2]);
                vertices->InsertNextCell(1, &pid);
            }
        }
        break;
    }
    case 1:
    {
        std::vector<float> tmpF = {0, 0, 0};
        for (int thread = 0; thread < pointArray.size(); thread++)
        {
            for (int i = 0; i < arraySize[thread]; i++)
            {
                for (int k = 0; k < pointArray[thread].size(); k++)
                    tmpF[k] = *((float*)((char*)pointArray[thread][0] + i * elemSize));

                if (problemType == 2)
                    pid = points->InsertNextPoint(tmpF[0], tmpF[1], 0);
                else
                    pid = points->InsertNextPoint(tmpF[0], tmpF[1], tmpF[2]);
                vertices->InsertNextCell(1, &pid);
            }
        }
        break;
    }
    }

    vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
    pointCloud->SetPoints(points);
    pointCloud->SetVerts(vertices);

    input->addVisualizationDataCloud(pointCloud, 0.5);

    for (int i = 0; i < currentModel->GetNumberBoundaries(); i++)
        input->addVisualizationDataBoundary(currentModel->GetBoundaryVTKUnstructuredGrid(i), colors::redcolor);
    input->setAxis();
};
void ShowAllCloud(vtkComponent* input, ModelInterface* currentModel, int precisionType, int problemType)
{
    input->CheckType(0);
    input->clear();
    int elemSize;

    std::vector<int>                arraySize;
    std::vector<std::vector<void*>> pointArray;

    for (int i = 0; i < currentModel->GetNumberBoundaries(); i++)
        input->addVisualizationDataBoundary(currentModel->GetBoundaryVTKUnstructuredGrid(i), colors::redcolor);

    input->setAxis();

    for (int k = 0; k < currentModel->GetNumberParticlesFlows(); k++)
    {
        currentModel->GetParticlesCloud(0, k, pointArray, arraySize, elemSize);

        vtkSmartPointer<vtkPoints>    points   = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
        vtkIdType                     pid;
        switch (precisionType)
        {
        case 0:
        {
            std::vector<double> tmp = {0, 0, 0};

            for (int thread = 0; thread < pointArray.size(); thread++)
            {
                for (int i = 0; i < arraySize[thread]; i++)
                {
                    for (int k = 0; k < pointArray[thread].size(); k++)
                        tmp[k] = *((double*)((char*)pointArray[thread][k] + i * elemSize));

                    if (problemType == 2)
                        pid = points->InsertNextPoint(tmp[0], tmp[1], 0);
                    else
                        pid = points->InsertNextPoint(tmp[0], tmp[1], tmp[2]);
                    vertices->InsertNextCell(1, &pid);
                }
            }
            break;
        }
        case 1:
        {
            std::vector<float> tmpF = {0, 0, 0};
            for (int thread = 0; thread < pointArray.size(); thread++)
            {
                for (int i = 0; i < arraySize[thread]; i++)
                {
                    for (int k = 0; k < pointArray[thread].size(); k++)
                        tmpF[k] = *((float*)((char*)pointArray[thread][0] + i * elemSize));

                    if (problemType == 2)
                        pid = points->InsertNextPoint(tmpF[0], tmpF[1], 0);
                    else
                        pid = points->InsertNextPoint(tmpF[0], tmpF[1], tmpF[2]);
                    vertices->InsertNextCell(1, &pid);
                }
            }
            break;
        }
        }

        vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
        pointCloud->SetPoints(points);
        pointCloud->SetVerts(vertices);

        if (k > colors::colorsPoints.size())
        {
            float randcolor[3];
            for (int i = 0; i < 3; i++)
                randcolor[i] = float(rand() % 255) / 255.0;
            input->addVisualizationDataCloud(pointCloud, 0.5, randcolor);
        }
        else
            input->addVisualizationDataCloud(pointCloud, 0.5, colors::colorsPoints[k]);
    }
    // input->addVisualizationDataBoundary(currentModel->GetMeshBoundaryVTKGrid(), colors::redcolor);
};

void ShowEmittanceDataPlot(vtkComponent* input, ModelInterface* currentModel, int flowNumber, int flagNumber,
                           int emFlag, int nCharts)
{
    input->CheckType(1);
    input->clear();

    std::vector<std::vector<float>> data;

    std::string xName, yName;

    if (emFlag <= 3)
    {
        currentModel->GetEmittanceData(data, flowNumber, flagNumber, emFlag);
        if (!data.size())
            return;
        switch (emFlag)
        {
        case 0:
            xName = "time, ns";
            yName = "Delta Pz, %";
            break;
        case 1:
            xName = "X, m";
            yName = "dX/dZ, rad";
            break;
        case 2:
            xName = "Y, m";
            yName = "dY/dZ, rad";
            break;
        case 3:
            xName = "X, m";
            yName = "Y, m";
            break;
        }

        input->addVisualizationDataPlot(data[0], data[1], colors::blackcolor, 4.0, xName, yName, 1, 1, 0,
                                        "distribution", 1.0, false);
    }
    else
    {
        volatile int N = 1;
        currentModel->GetEmittanceData(data, flowNumber, flagNumber, 0);
        if (!data.size())
            return;
        switch (emFlag)
        {
        case 4:
            xName = "Delta W, %";
            yName = "Number of particles, total  " + std::to_string(data[0].size()) + "  particles";
            N     = 1;
            break;
        case 5:
            xName = "Phi, deg";
            yName = "Number of particles, total  " + std::to_string(data[0].size()) + "  particles";
            N     = 0;
            break;
        }

        std::vector<float> XChart;
        std::vector<float> YChart;

        if (data[N].size() == 0)
            return;

        volatile float wmin = data[N][0];
        volatile float wmax = data[N][0];

        for (int i = 0; i < data[N].size(); i++)
        {
            if (data[N][i] > wmax)
                wmax = data[N][i];
            if (data[N][i] < wmin)
                wmin = data[N][i];
        };

        if (std::abs(wmax - wmin) < 1e-7)
            return;

        if (wmin < 0)
            wmin = wmin * 1.1;
        else
            wmin = wmin * 0.9;

        if (wmax > 0)
            wmax = wmax * 1.1;
        else
            wmax = wmax * 0.9;

        float d = (wmax - wmin) / nCharts;
        XChart.resize(nCharts + 1);
        YChart.resize(nCharts + 1);

        XChart[0]     = wmin;
        XChart.back() = wmax;

        for (int i = 1; i < nCharts; i++)
            XChart[i] = wmin + i * d;

        for (int i = 1; i < nCharts; i++)
        {
            for (int j = 0; j < data[N].size(); j++)
            {
                if (data[N][j] > XChart[i] && data[N][j] < XChart[i + 1])
                    YChart[i]++;
            }
        }

        input->addVisualizationDataPlot(XChart, YChart, colors::blackcolor, 5, xName, yName, 2, 1, 0, "distribution",
                                        1.0, false);
    }
};

void ShowPlotXY(vtkComponent* input, std::vector<float> X, std::vector<float> Y, std::string xName, std::string yName,
                int numberOfTraces, int Distrtype, std::vector<int> props)
{
    input->CheckType(1);
    input->clear();

    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();

    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName(xName.c_str());
    table->AddColumn(arrX);

    vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
    arrY->SetName(yName.c_str());
    table->AddColumn(arrY);

    table->SetNumberOfRows(X.size());
    for (int i = 0; i < X.size(); ++i)
    {
        table->SetValue(i, 0, X[i]);
        table->SetValue(i, 1, Y[i]);
    }

    // input->addVisualizationDataPlot(table);
};
void ShowCurrentDensity(vtkComponent* input, ModelInterface* currentModel, int flowNumber)
{
    input->CheckType(1);
    input->clear();

    std::vector<std::vector<float>> tmp   = currentModel->GetCurrentDensityDistribution(flowNumber);
    std::string                     xName = "Emitter Length";
    std::string                     yName = "Current Density";
    std::vector<const float*>       vc;
    vc.push_back(colors::redcolor);
    std::vector<std::vector<float>> Y;
    Y.push_back(tmp[1]);
    input->addVisualizationDataPlot(tmp[0], Y, vc, std::vector<float>{3.0}, xName, yName,
                                    std::vector<std::string>{"Emitter current"});
};
void ShowEmitterField(vtkComponent* input, ModelInterface* currentModel, int flowNumber)
{
    input->CheckType(1);
    input->clear();

    std::vector<std::vector<float>> tmp   = currentModel->GetEmitterFieldFloat(flowNumber);
    std::string                     xName = "Emitter Length";
    std::string                     yName = "Electric field, V/m";
    std::vector<const float*>       vc;
    vc.push_back(colors::redcolor);
    std::vector<std::vector<float>> Y;
    Y.push_back(tmp[1]);
    input->addVisualizationDataPlot(tmp[0], Y, vc, std::vector<float>{3.0}, xName, yName,
                                    std::vector<std::string>{"Electric field, V/m"});
};

void ShowPositionsWithGeometryPlaneSpecial(vtkComponent* input, ModelInterface* currentModel,
                                           std::shared_ptr<DynamicsData> data, int flowNumber, std::string xName,
                                           std::string yName, int problemType, int numberOfTraces, int Distrtype,
                                           std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    float ymax, ymin;

    std::vector<int> list = currentModel->GetBoundariesList();

    std::vector<std::vector<float>> XdataBoundary;
    std::vector<std::vector<float>> YdataBoundary;

    currentModel->GetPlotXYBoundarySpecial(XdataBoundary, YdataBoundary);

    input->addVisualizationDataPlotsAutoY(data->data[0], data->data[1], colors::blackcolor, 1.5, xName, yName,
                                          numberOfTraces, Distrtype, props, ymax, ymin);
    input->addVisualizationDataPlotsAutoY(XdataBoundary, YdataBoundary, colors::blackcolor, 3, xName, yName, 0,
                                          Distrtype, props, ymax, ymin);
};

void Show3PositionsWithGeometryPlane(vtkComponent* input, ModelInterface* currentModel,
                                     std::shared_ptr<DynamicsData> data, int flowNumber, std::string xName,
                                     std::string yName, int problemType, int numberOfTraces, int flag, int Distrtype,
                                     std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    std::vector<int> list = currentModel->GetBoundariesList();
    float            ymax, ymin;

    std::vector<std::vector<float>> XdataBoundary;
    std::vector<std::vector<float>> YdataBoundary;

    for (int i = 0; i < list.size(); i++)
    {
        std::vector<float> Xtmp;
        std::vector<float> Ytmp;
        currentModel->GetPlotXYBoundary(list[i], Xtmp, Ytmp);
        XdataBoundary.push_back(Xtmp);
        YdataBoundary.push_back(Ytmp);
    };

    if (flag == 0)
    {
        input->addVisualizationDataPlotsAutoY(data->data[0], data->data[1], colors::blackcolor, 1.5, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(XdataBoundary, YdataBoundary, colors::blackcolor, 3, xName, yName, 0,
                                              Distrtype, props, ymax, ymin);
    }
    if (flag == 1)
    {
        for (int i = 0; i < data->data[1].size(); i++)
        {
            if (data->data[1][i][0] < 0.025 && data->data[1][i][0] > 0.01 &&
                data->TimeArray[i].back() - data->TimeArray[i][0] > 7)
            {
               // input->addVisualizationDataPlot(data->data[1][i], data->data[0][i], colors::blackcolor, 1.5, yName,
                //                                xName, 0, 1, 0);
                break;
            }
        }

        for (int i = 0; i < data->data[1].size(); i++)
        {

            if (data->data[1][i][0] < 0.26 && data->data[1][i][0] > 0.25 &&
                data->TimeArray[i].back() - data->TimeArray[i][0] > 7)
            {
              //  input->addVisualizationDataPlot(data->data[1][i], data->data[0][i], colors::blackcolor, 1.5, yName,
            //                                  xName, 0, 1, 0);
                //	input->addVisualizationDataPlot(data->data[1][i], data->data[0][i], colors::redcolor, 1.5,
                // yName, xName, 0, 2, 0);
                break;
            }
        }

        for (int i = 0; i < data->data[1].size(); i++)
        {
            // if (data->data[1][i][0]>0.475 && data->TimeArray[i].back() - data->TimeArray[i][0]>20)
            if (data->data[1][i][0] > 0.475 && data->TimeArray[i].back() - data->TimeArray[i][0] > 7)
            {
              //  input->addVisualizationDataPlot(data->data[1][i], data->data[0][i], colors::blackcolor, 1.5, yName,
             //                                   xName, 0, 1, 0);
                //	input->addVisualizationDataPlot(data->data[1][i], data->data[0][i], colors::bluecolor, 2.5,
                // yName, xName, 0, 5, 0);
                break;
            }
        }

        input->addVisualizationDataPlotsAutoY(YdataBoundary, XdataBoundary, colors::blackcolor, 3, yName, xName, 0,
                                              Distrtype, props, ymax, ymin);
    }
};

void ShowPositionsWithGeometryPlane(vtkComponent* input, ModelInterface* currentModel,
                                    std::shared_ptr<DynamicsData> data, int flowNumber, std::string xName,
                                    std::string yName, int problemType, int numberOfTraces, int flag, int Distrtype,
                                    std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    float ymax, ymin;

    std::vector<int> list = currentModel->GetBoundariesList();

    std::vector<std::vector<float>> XdataBoundary;
    std::vector<std::vector<float>> YdataBoundary;

    for (int i = 0; i < list.size(); i++)
    {
        std::vector<float> Xtmp;
        std::vector<float> Ytmp;
        currentModel->GetPlotXYBoundary(list[i], Xtmp, Ytmp);
        XdataBoundary.push_back(Xtmp);
        YdataBoundary.push_back(Ytmp);
    };

    if (flag == 0)
    {
        input->addVisualizationDataPlotsAutoY(data->data[0], data->data[1], colors::blackcolor, 1.5, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(XdataBoundary, YdataBoundary, colors::blackcolor, 3, xName, yName, 0,
                                              Distrtype, props, ymax, ymin);
    }
    if (flag == 1)
    {
        input->addVisualizationDataPlotsAutoY(data->data[1], data->data[0], colors::blackcolor, 1.5, yName, xName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        input->addVisualizationDataPlotsAutoY(YdataBoundary, XdataBoundary, colors::blackcolor, 3, yName, xName, 0,
                                              Distrtype, props, ymax, ymin);
    }

    /*	if (problemType == 2)
            {
                    input->addVisualizationDataPlots(data->data[1], data->data[0], colors::blackcolor, 1, xName, yName);
                    input->addVisualizationDataPlots(YdataBoundary, XdataBoundary, colors::redcolor, 2, xName, yName);
            }
            else
            {
                    input->addVisualizationDataPlots(data->data[0], data->data[1], colors::blackcolor, 1, xName, yName);
                    input->addVisualizationDataPlots(XdataBoundary, YdataBoundary, colors::redcolor, 2, xName, yName);
            }*/
};
void ShowRT(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data, int flowNumber,
            int numberOfTraces, int Distrtype, std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    float       ymax, ymin;
    std::string yName = "R, m";
    std::string xName = "T, ns";
    input->addVisualizationDataPlotsAutoY(data->TimeArray, data->data[0], colors::blackcolor, 1, xName, yName,
                                          numberOfTraces, Distrtype, props, ymax, ymin);
};

void Show3dTrace(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data, int flowNumber,
                 int numberOfTraces, int flag, int Distrtype, std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    std::string yName, xName;
    float       ymax, ymin;
    switch (flag)
    {
    case 0:
        yName = "X, m";
        xName = "T, ns";
        input->addVisualizationDataPlotsAutoY(data->TimeArray, data->data[0], colors::blackcolor, 1, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        break;
    case 1:
        yName = "Y, m";
        xName = "T, ns";
        input->addVisualizationDataPlotsAutoY(data->TimeArray, data->data[1], colors::blackcolor, 1, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        break;
    case 2:
        yName = "Z, m";
        xName = "T, ns";
        input->addVisualizationDataPlotsAutoY(data->TimeArray, data->data[2], colors::blackcolor, 1, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        break;
    case 3:
        yName = "X, m";
        xName = "Z, m";
        input->addVisualizationDataPlotsAutoY(data->data[2], data->data[0], colors::blackcolor, 1, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        break;
    case 4:
        yName = "Y, m";
        xName = "Z, m";
        input->addVisualizationDataPlotsAutoY(data->data[2], data->data[1], colors::blackcolor, 1, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        break;
    case 5:
        yName = "X, m";
        xName = "Y, m";
        input->addVisualizationDataPlotsAutoY(data->data[1], data->data[0], colors::blackcolor, 1, xName, yName,
                                              numberOfTraces, Distrtype, props, ymax, ymin);
        break;
    }
}

void ShowTracesT(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                 std::string xName, std::string yName, int position, int flowNumber, int numberOfTraces, int Distrtype,
                 std::vector<int> props, int problemType)
{
    input->CheckType(1);
    input->clear();
    float ymax, ymin;
    if (xName == "time, ns")
        input->addVisualizationDataPlotsAutoY(data->TimeArray, data->data[position], colors::blackcolor, 1, xName,
                                              yName, numberOfTraces, Distrtype, props, ymax, ymin);
    if (xName == "Z, m")
    {
        if (problemType == 2)
            input->addVisualizationDataPlotsAutoY(data->data[1], data->data[position], colors::blackcolor, 1, xName,
                                                  yName, numberOfTraces, Distrtype, props, ymax, ymin);
        else
            input->addVisualizationDataPlotsAutoY(data->data[2], data->data[position], colors::blackcolor, 1, xName,
                                                  yName, numberOfTraces, Distrtype, props, ymax, ymin);
    }
}

void ShowEnergy(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                std::string xName, int flowNumber, int numberOfTraces, int Distrtype, std::vector<int> props,
                int axsFlag)
{
    input->CheckType(1);
    input->clear();
    std::string         yName = "Energy, eV";
    std::vector<double> prop  = currentModel->GetFlowProperties(flowNumber);

    std::vector<std::vector<float>> Ydata;
    float                           ymax, ymin;

    double gamma;
    double v;
    numberOfTraces = std::min(int(numberOfTraces), int(data->data[0].size()));

    float en = LIGHT_VELOCITY * LIGHT_VELOCITY * prop[2];

    for (int i = 0; i < data->data[0].size(); i++)
    {
        std::vector<float> Xtmp;
        for (int j = 0; j < data->data[0][i].size(); j++)
        {
            if (!axsFlag)
                gamma = sqrt(1 + data->data[3][i][j] * data->data[3][i][j] + data->data[4][i][j] * data->data[4][i][j] +
                             data->data[5][i][j] * data->data[5][i][j]);
            else
                gamma = sqrt(1 + data->data[3][i][j] * data->data[3][i][j] + data->data[4][i][j] * data->data[4][i][j] +
                             (data->data[5][i][j] / data->data[0][i][j]) * (data->data[5][i][j] / data->data[0][i][j]));

            if (gamma < 1.002)
            {
                if (!axsFlag)
                    v = LIGHT_VELOCITY *
                        sqrt(data->data[3][i][j] * data->data[3][i][j] + data->data[4][i][j] * data->data[4][i][j] +
                             data->data[5][i][j] * data->data[5][i][j]);
                else
                    v = LIGHT_VELOCITY *
                        sqrt(data->data[3][i][j] * data->data[3][i][j] + data->data[4][i][j] * data->data[4][i][j] +
                             (data->data[5][i][j] / data->data[0][i][j]) * (data->data[5][i][j] / data->data[0][i][j]));

                Xtmp.push_back(-prop[2] * v * v / (2 * ELECTRON_CHARGE));
            }
            else
                Xtmp.push_back(-(gamma - 1) * en / ELECTRON_CHARGE);
        }
        Ydata.push_back(Xtmp);
    }

    //if (xName == "time, ns")
    //    input->addVisualizationDataPlotsAutoY(data->TimeArray, Ydata, colors::blackcolor, 1, xName, yName,
     //                                         numberOfTraces, Distrtype, props, ymax, ymin);
    //if (xName == "Z, m")
    //    input->addVisualizationDataPlotsAutoY(data->data[2], Ydata, colors::blackcolor, 1, xName, yName, numberOfTraces,
     //                                         Distrtype, props, ymax, ymin);
}

void ShowErrors(vtkComponent* input, ModelInterface* currentModel)
{
    input->CheckType(1);
    input->clear();
    float                            ymax;
    float                            ymin;
    std::vector<float>               errors = currentModel->Get_Errors();
    std::vector<std::vector<double>> param;
    //	 = currentModel->GetSolverParameters();
    std::vector<float> Xdata;

    for (int i = 0; i < errors.size(); i++)
        Xdata.push_back(i * param[1][2]);

    std::string yName = "error";
    std::string xName = "T, ns";
    std::string flag  = "error";

    input->addVisualizationDataPlot(Xdata, errors, 2, xName, yName, flag, ymax, ymin);
};
void ShowXY(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data, int flowNumber,
            int numberOfTraces, int Distrtype, std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    std::string      xName = "X, m";
    std::string      yName = "Y, m";
    float            ymax, ymin;
    std::vector<int> list = currentModel->GetBoundariesList();

    std::vector<std::vector<float>> XdataBoundary;
    std::vector<std::vector<float>> YdataBoundary;

    std::vector<std::vector<float>> X;

    std::vector<int> indexes;

    numberOfTraces = std::min(int(numberOfTraces), int(data->data[0].size()));

    if (0 == numberOfTraces)
    {
        indexes.resize(data->data[0].size());
        for (int i = 0; i < data->data[0].size(); i++)
            indexes[i] = i;
    }
    else
    {
        indexes.resize(data->data[0].size());
        for (int i = 0; i < data->data[0].size(); i++)
            indexes[i] = i;

        for (int i = 0; i < data->data[0].size() - numberOfTraces; i++)
        {
            int k = rand() % indexes.size();
            indexes.erase(indexes.begin() + k);
        }
    }

    std::vector<std::vector<float>> Xtmp(indexes.size());
    std::vector<std::vector<float>> Ytmp(indexes.size());

    for (int i = 0; i < indexes.size(); i++)
    {

        for (int j = 0; j < data->data[0][indexes[i]].size(); j++)
        {
            Xtmp[i].push_back(data->data[0][indexes[i]][j] * std::sin(data->data[4][indexes[i]][j]));
            Ytmp[i].push_back(data->data[0][indexes[i]][j] * std::cos(data->data[4][indexes[i]][j]));
        }
        //	if (Xtmp.size()>1)
        //	input->addVisualizationDataPlot(Xtmp, Ytmp, colors::blackcolor, 1, xName, yName, 0, 1, 0);
    }
    input->addVisualizationDataPlotsAutoY(Xtmp, Ytmp, colors::blackcolor, 2, xName, yName, numberOfTraces, Distrtype,
                                          props, ymax, ymin);

    for (int i = 0; i < list.size(); i++)
    {
        std::vector<std::vector<float>> XdataBoundary;
        std::vector<std::vector<float>> YdataBoundary;
        currentModel->GetPlotBoundaryRotate(list[i], XdataBoundary, YdataBoundary);
        input->addVisualizationDataPlotsAutoY(YdataBoundary, XdataBoundary, colors::redcolor, 2, xName, yName,
                                              YdataBoundary.size(), Distrtype, props, ymax, ymin);
    };
};
void ShowPlotLine(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<lineplot>& plotData)
{
#ifdef SIM

    input->CheckType(0);
    input->clear();

    std::vector<int> list = currentModel->GetBoundariesList();
    for (int i = 0; i < list.size(); i++)
        input->addVisualizationDataBoundary(currentModel->GetBoundaryVTKUnstructuredGrid(list[i]), colors::redcolor);

    input->addVisualizationDataBoundary(plotData->line.GetVTKGrid(), colors::blackcolor);
#endif
};
void ShowLinePlot(vtkComponent* input, ModelInterface* currentModel, int plotNumber,
                  std::shared_ptr<lineplot>& plotData)
{
#ifdef SIM

    input->CheckType(1);
    input->clear();
    float                           ymax;
    float                           ymin;
    std::vector<float>              Xdata;
    std::vector<std::vector<float>> Ydata;
    currentModel->GetPlot(plotNumber, Xdata, Ydata);

    std::string yName;
    for (int i = 0; i < plotData->flag.size(); i++)
    {
        yName = yName + std::string(" ") + plotData->flag[i];
    }

    std::string xName = "L, m";

    for (int i = 0; i < plotData->flag.size(); i++)
    {
        input->addVisualizationDataPlot(Xdata, Ydata[i], 2, xName, yName, plotData->flag[i], ymax, ymin);
    }
#endif
}
void Show2dPlot(vtkComponent* input, ModelInterface* currentModel, std::string flag, int precisionType,
                int PlotTypeFlag)
{
    /*input->CheckType(2);
    input->clear();
    input->AddSurfacePlot();*/
    input->CheckType(0);
    input->clear();

    void* Array[1];
    int   size;
    int   elemSize;

    if (flag == std::string("Enorm, V/m"))
    {
        void* ArrayTmp1[1];
        void* ArrayTmp2[1];
        currentModel->GetGridData(ArrayTmp1, size, elemSize, 0, PlotTypeFlag);
        currentModel->GetGridData(ArrayTmp2, size, elemSize, 1, PlotTypeFlag);
        Array[0] = malloc(size * elemSize);

        switch (precisionType)
        {
        case 0:
            for (int i = 0; i < size; i++)
            {
                double valtmp1                               = *((double*)((char*)ArrayTmp1[0] + i * elemSize));
                double valtmp2                               = *((double*)((char*)ArrayTmp2[0] + i * elemSize));
                *((double*)((char*)Array[0] + i * elemSize)) = sqrt(valtmp1 * valtmp1 + valtmp2 * valtmp2);
            }
            break;
        case 1:
            for (int i = 0; i < size; i++)
            {
                float valtmp1                               = *((float*)((char*)ArrayTmp1[0] + i * elemSize));
                float valtmp2                               = *((float*)((char*)ArrayTmp2[0] + i * elemSize));
                *((float*)((char*)Array[0] + i * elemSize)) = sqrt(valtmp1 * valtmp1 + valtmp2 * valtmp2);
            }
            break;
        }
    }
    else
        currentModel->GetGridData(Array, size, elemSize, flag, PlotTypeFlag);

    vtkSmartPointer<vtkFloatArray> vtkData = vtkSmartPointer<vtkFloatArray>::New();

    vtkData->SetName(flag.c_str());

    vtkUnstructuredGrid*                 VTKgrid    = currentModel->GetVTKGrid();
    vtkSmartPointer<vtkUnstructuredGrid> newVTKgrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    newVTKgrid->ShallowCopy(VTKgrid);

    double minVal = *((double*)((char*)Array[0]));
    double maxVal = *((double*)((char*)Array[0]));
    double valtmp;
    for (int i = 1; i < size; i++)
    {
        valtmp = *((double*)((char*)Array[0] + i * elemSize));
        if (valtmp < minVal && !std::isnan(valtmp) && !std::isinf(valtmp))
            minVal = valtmp;
        if (valtmp > maxVal && !std::isnan(valtmp) && !std::isinf(valtmp))
            maxVal = valtmp;
    }
    double newmaxVal = maxVal - minVal;
    switch (precisionType)
    {
    case 0:
        for (int i = 0; i < size; i++)
        {
            valtmp = (*((double*)((char*)Array[0] + i * elemSize)) - minVal) / newmaxVal;
            vtkData->InsertNextValue(valtmp);
        }
        break;
    case 1:
        for (int i = 0; i < size; i++)
        {
            valtmp = *((float*)((char*)Array[0] + i * elemSize)) - minVal / newmaxVal;
            vtkData->InsertNextValue(valtmp);
        }
        break;
    }

    input->CheckType(0);
    input->clear();
    newVTKgrid->GetPointData()->SetScalars(vtkData);

    input->addVisualizationDataMesh(newVTKgrid, vtkData, minVal, maxVal, flag);
};

void Show2dPlot3d(vtkComponent* input, ModelInterface* currentModel, std::string flag, int flag1, double param,
                  int precisionType, int PlotTypeFlag)
{
    vtkSmartPointer<vtkFloatArray> vtkData = vtkSmartPointer<vtkFloatArray>::New();

    input->CheckType(0);
    input->clear();

    vtkData->SetName(flag.c_str());

    void* Array[1];
    int   size;
    int   elemSize;

    vtkUnstructuredGrid* VTKgrid = currentModel->GetVTKGrid(flag1, param, vtkData, flag);

    float valuesRange[2];
    vtkData->GetValueRange(valuesRange);
    double range = valuesRange[1] - valuesRange[0];
    float  tmp;
    if (range != 0)
    {
        for (int i = 0; i < vtkData->GetNumberOfTuples(); i++)
        {
            tmp = vtkData->GetValue(i);
            tmp = (tmp - valuesRange[0]) / range;
            vtkData->SetValue(i, tmp);
        };
    }

    input->addVisualizationDataBoundary(currentModel->GetMeshBoundaryVTKGrid(), colors::blackcolor, 3);
    VTKgrid->GetPointData()->SetScalars(vtkData);
    input->addVisualizationDataMesh(VTKgrid, vtkData, valuesRange[0], valuesRange[1], flag);
};

void ShowSimulationDataPlot(vtkComponent* input, ModelInterface* currentModel, int nSimData, int plotType, int nflow,
                            int nplot, const std::string& yName)
{
    float ymax;
    float ymin;
    input->CheckType(1);
    input->clear();

    std::string xName = "Time, ns";
    if (plotType == 0)
        input->addVisualizationDataPlot(currentModel->GetSimulationData()[nSimData]->XData,
                                        currentModel->GetSimulationData()[nSimData]->YDataFlow[nflow][nplot], 2, xName,
                                        yName, yName);

    if (plotType == 1)
        input->addVisualizationDataPlot(currentModel->GetSimulationData()[nSimData]->XData,
                                        currentModel->GetSimulationData()[nSimData]->YDataElectrode[nplot][0], 2, xName,
                                        yName, yName);
};

void ShowEnergyPolar(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                     int flowNumber, int numberOfTraces, int Distrtype, std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    std::string         yName = "Energy, eV";
    std::string         xName = "T, ns";
    std::vector<double> prop  = currentModel->GetFlowProperties(flowNumber);
    float               ymax, ymin;

    std::vector<std::vector<float>> Ydata;

    float gamma;

    numberOfTraces = std::min(int(numberOfTraces), int(data->data[0].size()));

    float en = LIGHT_VELOCITY * LIGHT_VELOCITY * prop[2];

    for (int i = 0; i < data->data[0].size(); i++)
    {
        std::vector<float> Xtmp;
        for (int j = 0; j < data->data[0][i].size(); j++)
        {
            gamma = sqrt(1 + data->data[2][i][j] * data->data[2][i][j] +
                         (data->data[3][i][j] / data->data[5][i][j]) * (data->data[3][i][j] / data->data[5][i][j]));
            Xtmp.push_back(-(gamma - 1) * en / ELECTRON_CHARGE);
        }
        Ydata.push_back(Xtmp);
    }

   // input->addVisualizationDataPlotsAutoY(data->TimeArray, Ydata, colors::blackcolor, 1, xName, yName, numberOfTraces,
  //                                        Distrtype, props, ymax, ymin);
}

void ShowEnergyEphiPolar(vtkComponent* input, ModelInterface* currentModel, std::shared_ptr<DynamicsData> data,
                         int flowNumber, int numberOfTraces, int Distrtype, std::vector<int> props)
{
    input->CheckType(1);
    input->clear();
    std::string         yName = "Energy, eV";
    std::string         xName = "T, ns";
    std::vector<double> prop  = currentModel->GetFlowProperties(flowNumber);
    numberOfTraces            = std::min(int(numberOfTraces), int(data->data[0].size()));

    std::vector<std::vector<float>> Ydata;
    float                           ymax, ymin;

    float gamma;

    float en = LIGHT_VELOCITY * LIGHT_VELOCITY * prop[2];

    for (int i = 0; i < data->data[0].size(); i++)
    {
        std::vector<float> Xtmp;
        for (int j = 0; j < data->data[0][i].size(); j++)
        {
            gamma = sqrt(1 + (data->data[3][i][j] / data->data[5][i][j]) * (data->data[3][i][j] / data->data[5][i][j]));
            Xtmp.push_back(-(gamma - 1) * en / ELECTRON_CHARGE);
        }
        Ydata.push_back(Xtmp);
    }

  //  input->addVisualizationDataPlots(data->TimeArray, Ydata, colors::blackcolor, 1, xName, yName, numberOfTraces,
   //                                  Distrtype, props, ymax, ymin);
}

void ShowValueAlongConductor(vtkComponent* input, ModelInterface* currentModel, int conductor, int flag)
{
    input->CheckType(1);
    input->clear();
    std::string xName = "conductor lenght";
    std::string yName;

    std::vector<std::vector<float>> tmp = currentModel->GetElectrodeValue(conductor, flag);
    if (flag == 0)
    {
        yName = "Av. irr. power density, W/(cm^2)";
        for (int i = 0; i < tmp[1].size(); i++)
            tmp[1][i] = tmp[1][i] * 1e-4;
    }

    if (flag == 1)
    {
        yName = "Av. irr. current density, A/(cm^2)";
        for (int i = 0; i < tmp[1].size(); i++)
            tmp[1][i] = tmp[1][i] * 1e-4;
    }

    if (flag == 2)
        yName = "Av. collected Current density, A/m";

    std::vector<const float*> vc;
    vc.push_back(colors::redcolor);
    // input->addVisualizationDataPlot(tmp[0], Y, vc, std::vector <float> {3.0}, xName, yName, std::vector < std::string
    // > {"power density"});
    input->addVisualizationDataPlot(tmp[0], tmp[1], 3.0, xName, yName, std::string{"Density along conductor"});
};