#include "GroupBoxWithItems.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "VisualizationFunctionsAccel.h"
#include "vtkComponent.h"

void MiddleWidget::InitTrees(int numberOfResultPlots)
{
    TreeList.resize(numberOfResultPlots);
    for (int i = 0; i < numberOfResultPlots; i++)
    {
        TreeList[i] = new QTreeWidget();
        TreeList[i]->setHeaderLabel("Simulations results");
        middleWidgetGrid->addWidget(TreeList[i], i, 0);
    }
    itemvector.resize(numberOfResultPlots);
};

void MiddleWidget::showAccelResuts(int numberOfResultPlots)
{
    clear();
    InitTrees(numberOfResultPlots);

    for (int k = 0; k < numberOfResultPlots; k++)
    {
        for (int t = 0; t < (*currentProject)->accelModel->GetSimulationData().size(); t++)
        {
            MyTreeItem* res = new MyTreeItem();
            res->setText(0, (*currentProject)->accelModel->GetSimulationData()[t]->tag.c_str());
            TreeList[k]->addTopLevelItem(res);
            res->flag1 = -1;
            for (int i = 0; i < (*currentProject)->accelModel->GetSimulationData()[t]->dataFlags.size(); i++)
            {
                itemvector[k].push_back(new MyTreeItem());
                itemvector[k].back()->flag = (*currentProject)->accelModel->GetSimulationData()[t]->dataFlags[i];
                itemvector[k].back()->setText(0, itemvector[k].back()->flag.c_str());
                itemvector[k].back()->flag1 = t;
                itemvector[k].back()->flag2 = k;
                res->addChild(itemvector[k].back());
                /*itemvector.push_back(new MyTreeItem());
                itemvector.back()->flag = (*currentProject)->accelModel->GetSimulationData()[t]->dataFlags[i];
                itemvector.back()->setText(0, itemvector.back()->flag.c_str());
                itemvector.back()->flag1 = t;
                itemvector.back()->flag = k;
                res->addChild(itemvector.back());*/
            };
        }
        connect(TreeList[k], SIGNAL(itemClicked(QTreeWidgetItem*, int)), this,
                SLOT(listShowResultClickAccel(QTreeWidgetItem*, int)));
    }

    groupBoxes.push_back(new GroupBoxWithItems("Main parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), numberOfResultPlots, 0);
    groupBoxes.back()->Create({"Number of particles", "L end (-1 means Lmax)"}, std::vector<int>{100, -1});
};

void MiddleWidget::listShowResultClickAccel(QTreeWidgetItem* item, int)
{

    std::vector<int>    props = {-1, -1, -1};
    std::vector<double> props1;
    bool                ok = groupBoxes[0]->GetParameters(props1);

    if (!ok)
        QMessageBox::critical(this, "Daisi error", "Incorrect input data, unable to plot result");

    MyTreeItem* item1 = dynamic_cast<MyTreeItem*>(item);

    if (item1->flag == "")
        return;

    if ((*currentProject)->problemType == 9)
        (*VTKArray)[item1->flag2]->setDataFunction(ShowSimulationNuclDataPlot, (*currentProject)->accelModel,
                                                   (*currentProject)->accelModel->GetSimulationData()[item1->flag1],
                                                   item1->flag, props, props1);
    else
        (*VTKArray)[item1->flag2]->setDataFunction(ShowSimulationAccelDataPlot, (*currentProject)->accelModel,
                                                   (*currentProject)->accelModel->GetSimulationData()[item1->flag1],
                                                   item1->flag, props, props1);

    (*VTKArray)[item1->flag2]->refresh(0);
    return;
}