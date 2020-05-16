#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "FlagStringsD.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "flagStringsTree.h"

void Daizy::SimulateIterative()
{
    currItemClickedSimSave(NULL);

    PrepareVisualisation();

    int device = 0;

    connect(timerViz, SIGNAL(timeout()), this, SLOT(updateGrViz()));

    timerViz->start(3000);

    if (work_thread.joinable())
        work_thread.join();

    progress  = 0;
    flagAbort = true;

    work_thread = std::thread(&ModelInterface::SimulateCPUPTI, currentProject->currentModel, std::ref(progress),
                              std::ref(flagAbort));
};
void Daizy::SimulatePIC()
{

    //	Visualization();
    currItemClickedSimSave(NULL);

    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;
    std::string                           errorMessage;

    bool ok = middleWidgetAggregator->FetchParameters(parameters1, parameters2, parameters3);

    PrepareVisualisation();

    int device = 0;

    connect(timerViz, SIGNAL(timeout()), this, SLOT(updateGrViz()));

    timerViz->start(1000 / parameters1[0][0]);

    if (work_thread.joinable())
        work_thread.join();

    progress  = 0;
    flagAbort = true;

    work_thread = std::thread(&ModelInterface::SimulateCPUPIC, currentProject->currentModel, std::ref(progress),
                              std::ref(progressLoc), std::ref(status), std::ref(flagAbort), 1, std::ref(errorMsg));
}
void Daizy::SimulatePICContinue()
{

    //	Visualization();
    currItemClickedSimSave(NULL);
    PrepareVisualisation();

    int device = 0;

    connect(timerViz, SIGNAL(timeout()), this, SLOT(updateGrViz()));

    timerViz->start(3000);

    if (work_thread.joinable())
        work_thread.join();

    progress  = 0;
    flagAbort = true;

    work_thread = std::thread(&ModelInterface::SimulateCPUPIC, currentProject->currentModel, std::ref(progress),
                              std::ref(progressLoc), std::ref(status), std::ref(flagAbort), 0, std::ref(errorMsg));
}

void Daizy::InitFieldSolver()
{
    currentProject->currentModel->InitFieldSolver();
};
void Daizy::FieldSimulate()
{
    currItemClickedSimSave(NULL);
    if (work_thread.joinable())
        work_thread.join();

    work_thread = std::thread(&ModelInterface::FieldSimulate, currentProject->currentModel, 0, std::ref(progress),
                              std::ref(progressLoc), std::ref(status));
    // currentProject->currentModel->FieldSimulate(0, progress, progressLoc,	status);
};

void Daizy::AddPotentialEvent()
{
    currentProject->currentModel->AddPropertyCondition(flagStrings::poisson, 0, std::string("potential"), 0);

    MyTreeItem* newItem = new MyTreeItem();
    int         i       = currentProject->currentModel->GetNumberPropertyConditions(flagStrings::poisson, 0);
    currentProject->currentModel->SetConditionProperties(currentClickItem->flagSearchBoundary, 0, i - 1,
                                                         {{0, 0, 0, 0}});

    newItem->setText(0, (std::string("potential") + std::to_string(i - 1)).c_str());
    newItem->flag               = flagStrings::potentialList;
    newItem->flag1              = i - 1;
    newItem->flagSearchBoundary = flagStrings::poisson;

    currentClickItem->addChild(newItem);
    currentClickItem->setExpanded(true);
    newItem->setSelected(true);
};

void Daizy::AddConductorEvent()
{
    currentProject->currentModel->AddConductor();

    MyTreeItem*                   newItem = new MyTreeItem();
    std::vector<std::vector<int>> List    = currentProject->currentModel->GetConductorsList();
    newItem->setText(0, QString("conductor") + QString::number(List.size() - 1));
    newItem->flag  = flagStrings::conductorsList;
    newItem->flag1 = List.size() - 1;

    currentClickItem->addChild(newItem);
    currentClickItem->setExpanded(true);
    newItem->setSelected(true);
};

void Daizy::AddFlowConditionEvent()
{

    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;

    bool ok = middleWidgetAggregator->FetchParameters(parameters1, parameters2, parameters3);

    int boundaryType     = parameters1[0][0];
    int boundaryTypeFlag = parameters1[1][0];

    currentProject->currentModel->AddPropertyCondition(flagStrings::flowBoundaryList, currentClickItem->flag1,
                                                       flagStrings::flowBoundaryTypeNames[boundaryType].c_str(),
                                                       boundaryTypeFlag);

    MyTreeItem* newItem = new MyTreeItem();
    int         p       = currentProject->currentModel->GetNumberPropertyConditions(flagStrings::flowBoundaryList,
                                                                      currentClickItem->flag1);

    currentProject->currentModel->SetConditionProperties(flagStrings::flowBoundaryList, currentClickItem->flag1, p - 1,
                                                         {{0, 0, 0, 0}, {0, 0, 0, 0}});

    newItem->setText(0, (flagStrings::flowBoundaryTypeNames[boundaryType].c_str() + std::to_string(p - 1)).c_str());
    newItem->flag                = flagStrings::flowBoundaryList;
    newItem->flag1               = p - 1;
    newItem->flag3               = flagStrings::flowBoundaryTypeNames[boundaryType];
    newItem->flagSearchBoundaryI = currentClickItem->flag1;
    newItem->flagSearchBoundary  = flagStrings::flowBoundaryList;
    currentClickItem->addChild(newItem);
    currentClickItem->setExpanded(true);
    newItem->setSelected(true);
};

void Daizy::GenerateMesh()
{
    // clearLayout(rigthWidgetGrid);

    //	work_thread.join();
    if (work_thread.joinable())
        work_thread.join();

    work_thread = std::thread(&ModelInterface::GenerateMesh, currentProject->currentModel,
                              currentProject->meshParamFileName, std::ref(progress), std::ref(errorMsg));

    // currentProject->currentModel->GenerateMesh(progress);
    // ShowMesh();
};

void Daizy::AddFlow()
{
    if (!currentProject->currentModel->GetVTKGrid())
    {
        QMessageBox::warning(this, "Daisi warning", "Mesh should be generated before flow adding");
        return;
    };

    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;

    bool ok = middleWidgetAggregator->FetchParameters(parameters1, parameters2, parameters3);

    if (!ok)
    {
        QMessageBox::critical(this, "Daisi error", "Incorrect input");
        return;
    };

    currentProject->currentModel->AddFlow(1, parameters1[1][0], parameters1[2][0], parameters1[2][1]);

    MyTreeItem* newItem = new MyTreeItem();
    int         i       = currentProject->currentModel->GetNumberParticlesFlows();
    newItem->setText(0,
                     QApplication::translate("MainWindow", (std::string("Flow") + std::to_string(i - 1)).c_str(), 0));
    newItem->flag  = flagStrings::flowList;
    newItem->flag1 = i - 1;
    currentClickItem->addChild(newItem);
    currentClickItem->setExpanded(true);

    if (currentProject->problemType >= 7)
        ShowProjectTreeAccel();
    else
        ShowProjectTreeSim();
};
