#include <QtWidgets/QSplitter>
#include <QVTKWidget.h>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QListWidgetItem>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "FlagStringsD.h"
#include "GeneralTools.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"


std::vector<std::string> simulationDataNamesFlowPIC   = {"N particles(t); of flow",       "I(t), A; of flow",
                                                       "Er_av(t), V/m; of flow",        "XdX RMS Emittance, pi*cm*mrad",
                                                       "YdY RMS Emittance, pi*cm*mrad", "dWPh RMS Emittance, kEv*ns"};
std::vector<std::string> simulationDataNamesElectrode = {"Collected current"};

void Daizy::RefreshGraphics(std::vector<vtkComponent*> VTKArray)
{

    for (int i = 0; i < VTKArray.size(); i++)
    {
        VTKArray[i]->refresh(0);
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(8000));
}
void Daizy::RefreshGraphics1(std::vector<QVTKWidget*> graphicsArray, QGridLayout* rigthWidgetGrid)
{

    std::vector<vtkComponent*> VTKArray;
    for (int i = 0; i < 3; i++)
    {
        VTKArray.push_back(new vtkComponent());
        VTKArray.back()->setWidget(graphicsArray.back());
    }

    for (int i = 0; i < VTKArray.size(); i++)
    {
        VTKArray[i]->refresh(0);
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
}

void Daizy::PrepareVisualisation()
{
    std::vector<std::string>               NowPlotNames = currentProject->currentModel->GetPlotNames();
    std::vector<std::shared_ptr<lineplot>> plotsData    = currentProject->currentModel->GetPlotsVector();

    VTKArray.clear();

    std::vector<std::string> names;
    for (int i = 0; i < NowPlotNames.size(); i++)
    {
        VTKArray.push_back(new vtkComponent());

        int s = currentProject->currentModel->GetSimulationData().size();

        for (int k = 0; k < simulationDataNamesFlowPIC.size(); k++)
        {
            if (simulationDataNamesFlowPIC[k] == NowPlotNames[i].substr(0, simulationDataNamesFlowPIC[k].size()))
            {
                int nflow = std::atoi(&NowPlotNames[i][NowPlotNames[i].size() - 1]);
                VTKArray.back()->setDataFunction(ShowSimulationDataPlot, currentProject->currentModel, s, 0, nflow, k,
                                                 NowPlotNames[i]);
                break;
            };
        }

        for (int k = 0; k < simulationDataNamesElectrode.size(); k++)
        {
            if (simulationDataNamesElectrode[k] == NowPlotNames[i].substr(0, simulationDataNamesElectrode[k].size()))
            {
                int nEl = std::atoi(&NowPlotNames[i][NowPlotNames[i].size() - 1]);
                VTKArray.back()->setDataFunction(ShowSimulationDataPlot, currentProject->currentModel, s, 1, k, nEl,
                                                 NowPlotNames[i]);
                break;
            };
        }

        for (int k = 0; k < plotsData.size(); k++)
        {
            if (NowPlotNames[i] == std::string("plot") + QString::number(k).toStdString())
            {
                VTKArray.back()->setDataFunction(ShowLinePlot, currentProject->currentModel, k, plotsData[k]);
                break;
            }
        };

        if (NowPlotNames[i] == flagStrings::VisNames[0])
            VTKArray.back()->setDataFunction(ShowAllCloud, currentProject->currentModel, currentProject->precisionType,
                                             (currentProject)->problemType);

        for (int k = 0; k < currentProject->currentModel->GetNumberParticlesFlows(); k++)
        {
            for (int j = 0; j < flagStrings::VisFlowNames.size(); j++)
            {
                if (NowPlotNames[i] == flagStrings::VisFlowNames[j] + QString::number(k).toStdString())
                {
                    if (j == 0)
                    {
                        VTKArray.back()->setDataFunction(ShowCurrentDensity, currentProject->currentModel, k);
                        break;
                    }
                    if (j == 1)
                    {
                        VTKArray.back()->setDataFunction(ShowEmitterField, currentProject->currentModel, k);
                        break;
                    }
                }
            }
        }

        for (int k = 0; k < currentProject->currentModel->GetConductorsList().size(); k++)
        {
            for (int j = 0; j < flagStrings::VisElectrodeNames.size(); j++)
            {
                if (NowPlotNames[i] == flagStrings::VisElectrodeNames[j] + QString::number(k).toStdString())
                {

                    VTKArray.back()->setDataFunction(ShowValueAlongConductor, currentProject->currentModel, k, j);
                    break;
                }
            }
        }
    };

    clearLayout(rigthWidgetGrid);

    graphicsArray.clear();
    int s  = VTKArray.size();
    int s1 = round(sqrt(s) + 0.49);

    int ind;
    for (int i = 0; i < s1; i++)
    {
        for (int j = 0; j < s1; j++)
        {
            ind = i * s1 + j;
            if (ind < s)
            {
                graphicsArray.push_back(new QVTKWidget());
                rigthWidgetGrid->addWidget(graphicsArray.back(), i, j);
                VTKArray[ind]->setWidget(graphicsArray.back());
            }
        }
    }
}

int  num = 0;
void Daizy::updateGrViz()
{
    num++;
    for (int i = 0; i < VTKArray.size(); i++)
    {
        VTKArray[i]->refresh(0);
    };

    // clearLayout(SolutionParamsGrid);

    /*	std::shared_ptr<SimulationData> simulationData = currentProject->currentModel->GetSimulationData()[0];

            int k = 0;
            int s = simulationData->YDataFlow[k].size()-1;

            if (s == -1)
                    return;

            int s1 = currentProject->currentModel->GetNumberParticlesFlows();

            for (int i = 0; i < s1; i++)
            {
                    for (int j = 0; j < simulationData->dataFlags.size() / s1; j++)
                    {
                            QLabel* label1 = new QLabel();
                            label1->setText(simulationData->dataFlags[k].c_str());
                            SolutionParamsGrid->addWidget(label1, i, 2*j);

                            QLabel* label2= new QLabel();
                            float r;
                            if (simulationData->YDataFlow[k].size() != 0)
                            //	r = simulationData->YDataFlow[k].back()
                                    ;

                            else
                                    r = 0;
                            label2->setText(QString::number(r));
                            SolutionParamsGrid->addWidget(label2, i, 2 * j+1);
                            k++;
                    }
            }
            std::vector<double> currents = currentProject->currentModel->GetElectrodesCurrents();
            for (int i = 0; i < currentProject->currentModel->GetConductorsList().size(); i++)
            {
                    QLabel* label1 = new QLabel();
                    label1->setText(QString("Electrode ") + QString::number(i) + QString(" current"));
                    SolutionParamsGrid->addWidget(label1, i+s1, 0);

                    QLabel* label2 = new QLabel();
                    label2->setText(QString::number(currents[i]));
                    SolutionParamsGrid->addWidget(label2, i + s1, 1);
            }*/
}