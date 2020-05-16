#include <QtWidgets/QMessageBox>
#include <QtWidgets/QTreeWidgetItem>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "flagStringsTree.h"

void Daizy::currItemClickedAccel(QTreeWidgetItem* item, int)
{
    MyTreeItem* item1 = dynamic_cast<MyTreeItem*>(item);
    currentClickItem  = item1;
    currItemClickedAccelSave(item1);

    middleWidgetAggregator->currentClickItem = item1;
    middleWidgetAggregator->clear();
    prevFlag       = item1->flag;
    prevFlagInt    = item1->flag1;
    prevFlagString = item1->flag3;
    if (item1->flag == flagStrings::AccelSolvers)
    {
        middleWidgetAggregator->showAccelSolverMenu(item1->flag3);
        currentsolver = item1->flag3;
        return;
    }

    if (item1->flag == flagStrings::AccelParams)
    {
        middleWidgetAggregator->showAccelParameters(0);
        return;
    }

    if (item1->flag == flagStrings::AccelParamsCalc)
    {
        middleWidgetAggregator->showAccelParametersCalc(0);
        return;
    }

    if (item1->flag == flagStrings::model || item1->flag == flagStrings::glDef)
    {
        middleWidgetAggregator->showSummary();
        return;
    }

    if (item1->flag == flagStrings::flowsAccel)
    {
        middleWidgetAggregator->ShowAddFlowAccelMenu();
        return;
    }

    if (item1->flag == flagStrings::flowListAccel)
    {
        middleWidgetAggregator->ShowFlowAccel(item1->flag1);
        return;
    }

    if (item1->flag == flagStrings::results)
    {
        ShowResultsMenu();
        return;
    }
};
void Daizy::currItemClickedAccelSave(MyTreeItem* item1)
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

    if (prevFlag == flagStrings::AccelParams)
    {
        currentProject->accelModel->SetMainAccelParameters(parameters1[0]);
        return;
    }

    if (prevFlag == flagStrings::flowListAccel)
    {
        currentProject->accelModel->SetParametersAccelFlow(parameters1[0], prevFlagInt);
        return;
    }

    if (prevFlag == flagStrings::AccelSolvers)
    {
        std::vector<double> parameters22;
        for (int i = 0; i < parameters1.size() - 1; i++)
            parameters22.push_back(parameters1[1 + i][0]);

        currentProject->accelModel->SetSolverAllParameters(prevFlagString, parameters3[0], parameters2[0],
                                                           parameters1[0], parameters22, errorMessage);
        return;
    }
};
