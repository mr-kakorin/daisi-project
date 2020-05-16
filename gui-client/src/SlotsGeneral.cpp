#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "Menu.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#ifdef NUCL
#include "namesNucl.h"
#else
#include "names.h"
#endif // NUCL
void Daizy::MyConnect()
{
    connect(middleWidgetAggregator, SIGNAL(CreateNewProject()), this, SLOT(CreateNewProjectEvent()));
    connect(middleWidgetAggregator, SIGNAL(AddFlowAccel()), this, SLOT(AddFlowAccel()));
    connect(middleWidgetAggregator, SIGNAL(AccelSolve()), this, SLOT(AccelSolve()));
    connect(middleWidgetAggregator, SIGNAL(ShowProjectTreeSim()), this, SLOT(ShowProjectTreeSim()));
    connect(middleWidgetAggregator, SIGNAL(ShowProjectTreeAccel()), this, SLOT(ShowProjectTreeAccel()));
    connect(middleWidgetAggregator, SIGNAL(AddPotentialEvent()), this, SLOT(AddPotentialEvent()));
    connect(middleWidgetAggregator, SIGNAL(AddConductorEvent()), this, SLOT(AddConductorEvent()));
    connect(middleWidgetAggregator, SIGNAL(AddFlowConditionEvent()), this, SLOT(AddFlowConditionEvent()));
    connect(middleWidgetAggregator, SIGNAL(GenerateMesh()), this, SLOT(GenerateMesh()));
    connect(middleWidgetAggregator, SIGNAL(AddFlow()), this, SLOT(AddFlow()));
    connect(middleWidgetAggregator, SIGNAL(DefaultBuilder()), this, SLOT(DefaultBuilder()));
    connect(middleWidgetAggregator, SIGNAL(FieldSimulate()), this, SLOT(FieldSimulate()));
    connect(middleWidgetAggregator, SIGNAL(ResetFlags()), this, SLOT(ResetFlags()));

    connect(middleWidgetAggregator, SIGNAL(SimulateIterative()), this, SLOT(SimulateIterative()));
    connect(middleWidgetAggregator, SIGNAL(SimulatePIC()), this, SLOT(SimulatePIC()));
    connect(middleWidgetAggregator, SIGNAL(SimulatePICContinue()), this, SLOT(SimulatePICContinue()));

    connect(mainMenu, SIGNAL(ShowProjectTreeAccel()), this, SLOT(ShowProjectTreeAccel()));
    connect(mainMenu, SIGNAL(ShowProjectTreeSim()), this, SLOT(ShowProjectTreeSim()));
    connect(mainMenu, SIGNAL(ResetFlags()), this, SLOT(ResetFlags()));
    connect(mainMenu, SIGNAL(saveParameters()), this, SLOT(saveParameters()));
    connect(mainMenu, SIGNAL(saveParameters()), this, SLOT(saveParameters()));
    connect(mainMenu, SIGNAL(showSummary()), middleWidgetAggregator, SLOT(showSummary()));
    connect(mainMenu, SIGNAL(clear()), middleWidgetAggregator, SLOT(clear()));
    connect(mainMenu, SIGNAL(ShowCreateNewProject()), middleWidgetAggregator, SLOT(ShowCreateNewProject()));
};

void Daizy::DefaultBuilder()
{
    leftWidgetTree->clear();
    QTreeWidgetItem* qtreewidgetitem = new QTreeWidgetItem();
    qtreewidgetitem->setText(0, "ModelBuilder");
    leftWidgetTree->setHeaderItem(qtreewidgetitem);
};

void Daizy::CreateNewProjectEvent()
{
    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;

    middleWidgetAggregator->FetchParameters(parameters1, parameters2, parameters3);

    currentProject = new Dproject::project(int(parameters1[0][0]), int(parameters1[1][0]), 1, parameters2[0][1],
                                           parameters2[0][0], version);

    if (currentProject->problemType >= 7)
        ShowProjectTreeAccel();
    else
        ShowProjectTreeSim();

    middleWidgetAggregator->clear();
    middleWidgetAggregator->showSummary();
};

void Daizy::saveParameters()
{
    if (!currentProject)
        return;
    if (currentProject->problemType >= 7)
        currItemClickedAccelSave(currentClickItem);
    else
        currItemClickedSimSave(currentClickItem);
};
void Daizy::ResetFlags()
{
    stopTimer();
    prevFlag                = "";
    prevFlagInt             = 0;
    currentClickItem        = NULL;
    prevflagSearchBoundary  = "";
    prevflagSearchBoundaryI = 0;
};

void Daizy::VisualThr(QListWidgetItem* item)
{
    QMetaObject::invokeMethod(this, "Visual", Qt::QueuedConnection, Q_ARG(QListWidgetItem*, item));
}
