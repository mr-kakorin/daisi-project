#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "flagStringsTree.h"
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QTreeWidgetItem>

//#include "flagStrings.h"
std::vector<std::string> namesProblems = {"RFQ", "DTL", "Synchrotron"};
void                     Daizy::ShowProjectTreeAccel()
{
    leftWidgetTree->clear();

    leftWidgetTree->setContextMenuPolicy(Qt::CustomContextMenu);
    //	connect(leftWidgetTree, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this,
    // SLOT(currItemClickedAccel(QTreeWidgetItem*, int)));
    QObject::connect(leftWidgetTree, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this,
                     SLOT(currItemClickedAccel(QTreeWidgetItem*, int)));

    MyTreeItem* header = new MyTreeItem();
    header->setText(0, currentProject->projectName.c_str());
    leftWidgetTree->setHeaderItem(header);

    MyTreeItem* model = new MyTreeItem();
    model->setText(0, currentProject->currentModelFileName.c_str());
    model->flag = flagStrings::model;
    leftWidgetTree->addTopLevelItem(model);

    MyTreeItem* Accel = new MyTreeItem();
    Accel->setText(0, namesProblems[currentProject->problemType - 7].c_str());
    model->addChild(Accel);

    MyTreeItem* AccelParams = new MyTreeItem();
    AccelParams->setText(0, "Accelerator parameters");
    AccelParams->flag = flagStrings::AccelParams;
    Accel->addChild(AccelParams);

    MyTreeItem* AccelParamsCalc = new MyTreeItem();
    AccelParamsCalc->setText(0, "Calculated accelerator parameters");
    AccelParamsCalc->flag = flagStrings::AccelParamsCalc;
    Accel->addChild(AccelParamsCalc);

    MyTreeItem* Flows = new MyTreeItem();
    Flows->setText(0, "Particles flows");
    Flows->flag = flagStrings::flowsAccel;
    Accel->addChild(Flows);

    std::vector<MyTreeItem*> FlowsVector;

    for (int i = 0; i < currentProject->accelModel->GetNumberOfAccelFlows(); i++)
    {
        FlowsVector.push_back(new MyTreeItem());
        FlowsVector.back()->setText(0, (std::string("Flow") + std::to_string(i)).c_str());
        FlowsVector.back()->flag  = flagStrings::flowListAccel;
        FlowsVector.back()->flag1 = i;
        Flows->addChild(FlowsVector.back());
    }

    MyTreeItem* Solvers = new MyTreeItem();
    Solvers->setText(0, "Solvers");
    // Solvers->flag = flagStrings::flowsAccel;
    Accel->addChild(Solvers);

    std::vector<std::string> sNames = currentProject->accelModel->GetSolversNames();

    for (int i = 0; i < sNames.size(); i++)
    {
        MyTreeItem* Solver = new MyTreeItem();
        Solver->setText(0, sNames[i].c_str());
        Solver->flag  = flagStrings::AccelSolvers;
        Solver->flag1 = i;
        Solver->flag3 = sNames[i];
        Solvers->addChild(Solver);
    }

    MyTreeItem* Results = new MyTreeItem();
    Results->setText(0, "Results");
    Results->flag = flagStrings::results;
    leftWidgetTree->addTopLevelItem(Results);
};