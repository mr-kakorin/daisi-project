#include <QVTKWidget.h>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QTreeWidgetItem>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "GeneralTools.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "flagStringsTree.h"
#include "vtkComponent.h"

void Daizy::currItemClicked(QTreeWidgetItem* item, int)
{
    MyTreeItem* item1 = dynamic_cast<MyTreeItem*>(item);
    currentClickItem  = item1;
    currItemClickedSimSave(item1);

    middleWidgetAggregator->currentClickItem = item1;
    middleWidgetAggregator->clear();

    prevFlag                = item1->flag;
    prevFlagInt             = item1->flag1;
    prevflagSearchBoundary  = item1->flagSearchBoundary;
    prevflagSearchBoundaryI = item1->flagSearchBoundaryI;

    if (item1->flag != flagStrings::results && lastFlag == flagStrings::results)
    {
        /*VTKArray.clear();
        graphicsArray.clear();
        graphicsArray.push_back(new QVTKWidget());
        VTKArray.push_back(new vtkComponent());
        VTKArray[0]->setWidget(graphicsArray[0]);
        clearLayout(rigthWidgetGrid);
        rigthWidgetGrid->addWidget(graphicsArray.back(), 0, 0);*/
    };
    stopTimer();

    lastFlag = item1->flag;

    if (item1->flag == flagStrings::model || item1->flag == flagStrings::glDef)
    {
        middleWidgetAggregator->showSummary();
        return;
    }
    if (item1->flag == flagStrings::boundaries)
    {
        middleWidgetAggregator->ShowAddBoundary();
        return;
    }

    if (item1->flag == flagStrings::boundariesList)
    {
        std::vector<int> allBoundariesList = currentProject->currentModel->GetBoundariesList();
        allBoundariesList.erase(allBoundariesList.begin() + item1->flag1);
        VTKArray[0]->setDataFunction(ShowBoundaries, currentProject->currentModel, allBoundariesList,
                                     colors::blackcolor, std::vector<int>{item1->flag1}, colors::redcolor);
        VTKArray[0]->refresh(0);
        // VTKArray[0]->refresh();
        return;
    }

    if (item1->flag == flagStrings::neumann)
    {
        middleWidgetAggregator->ShowDefaultCondMenu(item1->flag1, QString("Neumann conditions"), flagStrings::poisson);
        return;
    }

    if (item1->flag == flagStrings::fabsopbtion)
    {
        middleWidgetAggregator->ShowDefaultCondMenu(item1->flag1, QString("Absorption conditions"),
                                                    flagStrings::flowBoundaryList);
        return;
    }

    if (item1->flag == flagStrings::conductors)
    {
        middleWidgetAggregator->AddConductor();
        return;
    }

    if (item1->flag == flagStrings::flowBoundary)
    {
        middleWidgetAggregator->AddFlowCondition();
        return;
    }
    if (item1->flag == flagStrings::poisson)
    {
        middleWidgetAggregator->AddCondition();
        return;
    }
    if (item1->flag == flagStrings::globalField)
    {
        middleWidgetAggregator->ShowGlobalFieldMenu();
        return;
    };
    if (item1->flag == flagStrings::potentialList)
    {
        middleWidgetAggregator->PotentialFieldMenu(item1->flag1);
        return;
    }

    if (item1->flag == flagStrings::conductorsList)
    {
        middleWidgetAggregator->ShowConductorSelectionMenu(item1->flag1);
        return;
    }

    if (item1->flag == flagStrings::flowBoundaryList)
    {
        middleWidgetAggregator->currentConditioin = item1->flag1;
        middleWidgetAggregator->ShowSelectionMenu(item1->flag1, item1->flagSearchBoundary, item1->flagSearchBoundaryI);
        return;
        /*if (item1->flag3 == flagStrings::flowBoundaryTypeNames[0])
        {
                middleWidgetAggregator->ShowSelectionMenu(item1->flag1, item1->flagSearchBoundary,
        item1->flagSearchBoundaryI, { std::string("transparency value") }); return;
        }
        if (item1->flag3 == flagStrings::flowBoundaryTypeNames[1])
        {
                middleWidgetAggregator->ShowSelectionMenu(item1->flag1, item1->flagSearchBoundary,
        item1->flagSearchBoundaryI); return;
        }
        if (item1->flag3 == flagStrings::flowBoundaryTypeNames[2])
        {
                middleWidgetAggregator->ShowSelectionMenu(item1->flag1, item1->flagSearchBoundary,
        item1->flagSearchBoundaryI, { std::string("alpha fraction"), std::string("beta fraction") }); return;
        }*/
    }

    if (item1->flag == flagStrings::mesh)
    {
        middleWidgetAggregator->ShowMeshMenu();
        VTKArray[0]->setDataFunction(ShowMesh, currentProject->currentModel, colors::blackcolor);
        VTKArray[0]->refresh(0);
        return;
    }
    if (item1->flag == flagStrings::flows)
    {
        middleWidgetAggregator->ShowAddFlowMenu();
        return;
    }

    if (item1->flag == flagStrings::flowList)
    {
        middleWidgetAggregator->ShowFlowSummary(item1->flag1);
        return;
    }

    if (item1->flag == flagStrings::emitterList)
    {
        middleWidgetAggregator->ShowFlowEmitterProperty(item1->flag1);

        // VTKArray[0]->setDataFunction(ShowBoundaries, currentProject->currentModel, List2Vector(leftList),
        // colors::blackcolor, List2Vector(rigthList), colors::redcolor);  VTKArray[0]->refresh(0);

        /*std::vector<double> prop = currentProject->currentModel->GetFlowProperties(item1->flag1);
        int distributionStyle = (int)prop[1];
        if (1 == distributionStyle)
                return;
        VTKArray[0]->setDataFunction(ShowBoundaries, currentProject->currentModel, List2Vector(leftList),
        colors::blackcolor, List2Vector(rigthList), colors::redcolor); VTKArray[0]->refresh(0);*/

        return;
    }
    if (item1->flag == flagStrings::beamState)
    {
        currentFlow = item1->flag1;
        middleWidgetAggregator->ShowFlowState(item1->flag1);
        return;
    }
    if (item1->flag == flagStrings::Visuaization)
    {
        middleWidgetAggregator->ShowAddPlotMenu();
        return;
    }

    if (item1->flag == flagStrings::solver)
    {
        return;
    }
    if (item1->flag == flagStrings::fieldSolver)
    {
        middleWidgetAggregator->ShowFieldSolverMenu();
        return;
    }
    if (item1->flag == flagStrings::simulatePIC)
    {
        middleWidgetAggregator->ShowSolverMenuPIC();
        return;
    }

    if (item1->flag == flagStrings::simulatePTI)
    {
        middleWidgetAggregator->ShowSolverMenuPTI();
        return;
    }

    if (item1->flag == flagStrings::solverSettings)
    {
        middleWidgetAggregator->ShowSolverSettings();
        return;
    }

    if (item1->flag == flagStrings::results)
    {
        ShowResultsMenu();
        return;
    }
    if (item1->flag == flagStrings::lineplots)
    {
        middleWidgetAggregator->ShowAddLinePlotMenu();
        return;
    }
    if (item1->flag == flagStrings::lineplotsList)
    {
        currentplotNumber = item1->flag1;
        middleWidgetAggregator->ShowLinePlotMenu(currentplotNumber);
        return;
    }

    if (item1->flag == flagStrings::plots2d)
    {
        currentplotNumber = item1->flag1;
        middleWidgetAggregator->ShowPlot2dMenu();
        return;
    }

    if (item1->flag == flagStrings::solverEmission)
    {
        middleWidgetAggregator->ShowEmissionModelSolverSettings();
        return;
    };
};

void Daizy::currItemClickedSimSave(MyTreeItem* item1)
{
    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;
    std::string                           errorMessage;

    bool ok = middleWidgetAggregator->FetchParameters(parameters1, parameters2, parameters3);

    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Incorrect non-digit input");
        return;
    }

    if (prevFlag == flagStrings::globalField)
        currentProject->currentModel->SetglobalFieldConditions(parameters1[0]);

    if (prevFlag == flagStrings::potentialList || prevFlag == flagStrings::flowBoundaryList)
        currentProject->currentModel->SetConditionProperties(prevflagSearchBoundary, prevflagSearchBoundaryI,
                                                             prevFlagInt, parameters1);

    if (prevFlag == flagStrings::flowList)
        currentProject->currentModel->SetFlowMCNumbers(prevFlagInt, parameters1[1]);

    if (prevFlag == flagStrings::emitterList)
    {
        std::vector<double> prop              = currentProject->currentModel->GetFlowProperties(prevFlagInt);
        int                 distributionStyle = (int)prop[1];

        if (distributionStyle == 5)
        {
            if (parameters1.size() > 1)
                currentProject->currentModel->SetAllEmitterParameters(
                    prevFlagInt, std::vector<std::vector<double>>(parameters1.begin(), parameters1.end()));
        }
        else
        {
            if (parameters1.size() > 1)
                currentProject->currentModel->SetAllEmitterParameters(
                    prevFlagInt, std::vector<std::vector<double>>(parameters1.begin() + 1, parameters1.end()));
        }
    }
    if (prevFlag == flagStrings::conductorsList)
        middleWidgetAggregator->ApplyConductorBoundaries(prevFlagInt);

    if (prevFlag == flagStrings::fieldSolver)
        currentProject->currentModel->SetFieldSolverParameters(parameters1[0]);

    if (prevFlag == flagStrings::solverEmission)
        currentProject->currentModel->SetSolverEmissionModelParameters(parameters1);

    if (prevFlag == flagStrings::solverSettings)
        currentProject->currentModel->SetSolverGeneralParameters(parameters1);

    if (prevFlag == flagStrings::simulatePIC)
        middleWidgetAggregator->ApplySolverSettingsPIC();
};