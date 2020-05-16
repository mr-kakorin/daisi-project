#include <QtWidgets/QMessageBox>
#include <QtWidgets/QTreeWidgetItem>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "colors.h"
#include "flagStringsTree.h"
#include "vtkComponent.h"


void Daizy::ShowProjectTreeSim()
{
    leftWidgetTree->clear();

    // MyTreeItem * header = new MyTreeItem();
    // header->setText(0, currentProject->projectName.c_str());
    // leftWidgetTree->setHeaderItem(header);

    MyTreeItem* model = new MyTreeItem();
    model->setText(0, currentProject->currentModelFileName.c_str());
    model->flag = flagStrings::model;
    leftWidgetTree->addTopLevelItem(model);

    MyTreeItem* Device = new MyTreeItem();
    Device->setText(0, "Device");
    Device->flag = flagStrings::device;
    model->addChild(Device);

    MyTreeItem* Boundaries = new MyTreeItem();
    Boundaries->setText(0, "Boundaries");
    Boundaries->flag = flagStrings::boundaries;
    Device->addChild(Boundaries);

    std::vector<MyTreeItem*> BoundariesVector;
    for (int i = 0; i < currentProject->currentModel->GetNumberBoundaries(); i++)
    {
        BoundariesVector.push_back(new MyTreeItem());
        BoundariesVector.back()->setText(0, (std::string("boundary") + std::to_string(i)).c_str());
        BoundariesVector.back()->flag  = flagStrings::boundariesList;
        BoundariesVector.back()->flag1 = i;
        Boundaries->addChild(BoundariesVector.back());
    };

    MyTreeItem* Poisson = new MyTreeItem();
    Poisson->setText(0, "Electric fields");
    Poisson->flag = flagStrings::poisson;
    Device->addChild(Poisson);

    MyTreeItem* Ampere = new MyTreeItem();
    Ampere->setText(0, "Magnetic fields");
    Ampere->flag = flagStrings::ampere;
    Device->addChild(Ampere);

    MyTreeItem* Insulation = new MyTreeItem();
    Insulation->setText(0, "Magnetic insulation");
    Insulation->flag               = flagStrings::insulation;
    Insulation->flagSearchBoundary = flagStrings::poisson;
    Ampere->addChild(Insulation);

    MyTreeItem* Neumann = new MyTreeItem();
    Neumann->setText(0, "Zero charge");
    Neumann->flag               = flagStrings::neumann;
    Neumann->flagSearchBoundary = flagStrings::poisson;
    Poisson->addChild(Neumann);

    MyTreeItem* globalField = new MyTreeItem();
    globalField->setText(0, "Global parameters");
    globalField->flag = flagStrings::globalField;
    Poisson->addChild(globalField);

    std::vector<MyTreeItem*> PotentialVector;
    for (int i = 0; i < currentProject->currentModel->GetNumberPropertyConditions(flagStrings::poisson, 0); i++)
    {
        PotentialVector.push_back(new MyTreeItem());
        PotentialVector.back()->setText(
            0, (currentProject->currentModel->GetConditionPropertyType(flagStrings::poisson, 0, i) + std::to_string(i))
                   .c_str());
        PotentialVector.back()->flag               = flagStrings::potentialList;
        PotentialVector.back()->flag1              = i;
        PotentialVector.back()->flagSearchBoundary = flagStrings::poisson;
        Poisson->addChild(PotentialVector.back());
    };

    MyTreeItem* Conductors = new MyTreeItem();
    Conductors->setText(0, "Conductors");
    Conductors->flag = flagStrings::conductors;
    Device->addChild(Conductors);

    std::vector<std::vector<int>> conductorsList = currentProject->currentModel->GetConductorsList();
    std::vector<MyTreeItem*>      conductorsVector;

    for (int i = 0; i < conductorsList.size(); i++)
    {
        conductorsVector.push_back(new MyTreeItem());
        conductorsVector.back()->setText(0, (std::string("conductor") + std::to_string(i)).c_str());
        conductorsVector.back()->flag  = flagStrings::conductorsList;
        conductorsVector.back()->flag1 = i;
        Conductors->addChild(conductorsVector.back());
    };

    MyTreeItem* Mesh = new MyTreeItem();
    Mesh->setText(0, "Mesh");
    Mesh->flag = flagStrings::mesh;
    Device->addChild(Mesh);

    MyTreeItem* Flows = new MyTreeItem();
    Flows->setText(0, "Particles flows");
    Flows->flag = flagStrings::flows;
    Device->addChild(Flows);

    std::vector<MyTreeItem*> FlowsVector;
    std::vector<MyTreeItem*> EmitterVector;
    std::vector<MyTreeItem*> BeamStateVector;

    for (int i = 0; i < currentProject->currentModel->GetNumberParticlesFlows(); i++)
    {
        FlowsVector.push_back(new MyTreeItem());
        FlowsVector.back()->setText(0, (std::string("Flow") + std::to_string(i)).c_str());
        FlowsVector.back()->flag  = flagStrings::flowList;
        FlowsVector.back()->flag1 = i;
        Flows->addChild(FlowsVector.back());

        EmitterVector.push_back(new MyTreeItem());
        EmitterVector.back()->setText(0, (std::string("Emitter") + std::to_string(i)).c_str());
        EmitterVector.back()->flag  = flagStrings::emitterList;
        EmitterVector.back()->flag1 = i;
        FlowsVector.back()->addChild(EmitterVector.back());

        BeamStateVector.push_back(new MyTreeItem());
        BeamStateVector.back()->setText(0, (std::string("Flow State") + std::to_string(i)).c_str());
        BeamStateVector.back()->flag  = flagStrings::beamState;
        BeamStateVector.back()->flag1 = i;
        FlowsVector.back()->addChild(BeamStateVector.back());

        MyTreeItem* FlowBoundary = new MyTreeItem();
        FlowBoundary->setText(0, "Flow boundary conditions");
        FlowBoundary->flag  = flagStrings::flowBoundary;
        FlowBoundary->flag1 = i;
        FlowsVector.back()->addChild(FlowBoundary);

        MyTreeItem* Absopbtion = new MyTreeItem();
        Absopbtion->setText(0, "Absolut Absopbtion");
        Absopbtion->flag                = flagStrings::fabsopbtion;
        Absopbtion->flag1               = i;
        Absopbtion->flagSearchBoundary  = flagStrings::flowBoundaryList;
        Absopbtion->flagSearchBoundaryI = i;
        FlowBoundary->addChild(Absopbtion);

        std::vector<MyTreeItem*> BoundaryVector;
        for (int k = 0; k < currentProject->currentModel->GetNumberPropertyConditions(flagStrings::flowBoundaryList, i);
             k++)
        {
            BoundaryVector.push_back(new MyTreeItem());
            BoundaryVector.back()->flag3 =
                currentProject->currentModel->GetConditionPropertyType(flagStrings::flowBoundaryList, i, k);
            BoundaryVector.back()->setText(
                0, (currentProject->currentModel->GetConditionPropertyType(flagStrings::flowBoundaryList, i, k) +
                    std::to_string(k))
                       .c_str());
            BoundaryVector.back()->flag                = flagStrings::flowBoundaryList;
            BoundaryVector.back()->flagSearchBoundaryI = i;
            BoundaryVector.back()->flag1               = k;
            BoundaryVector.back()->flagSearchBoundary  = flagStrings::flowBoundaryList;
            FlowBoundary->addChild(BoundaryVector.back());
        };
    };

    MyTreeItem* State = new MyTreeItem();
    State->setText(0, "Device state");
    State->flag = flagStrings::devicestate;
    Device->addChild(State);

    MyTreeItem* lineplots = new MyTreeItem();
    lineplots->setText(0, "line plots");
    lineplots->flag = flagStrings::lineplots;
    State->addChild(lineplots);

    MyTreeItem* plots2d = new MyTreeItem();
    plots2d->setText(0, "2d plots");
    plots2d->flag = flagStrings::plots2d;
    State->addChild(plots2d);

    std::vector<std::shared_ptr<lineplot>> plotsData = currentProject->currentModel->GetPlotsVector();

    std::vector<MyTreeItem*> plotsVector;

    for (int i = 0; i < plotsData.size(); i++)
    {
        plotsVector.push_back(new MyTreeItem());
        plotsVector.back()->flag = flagStrings::lineplotsList;
        plotsVector.back()->setText(0, QString("plot") + QString::number(i));
        plotsVector.back()->flag1 = i;
        lineplots->addChild(plotsVector.back());
    }

    leftWidgetTree->setContextMenuPolicy(Qt::CustomContextMenu);

    connect(leftWidgetTree, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this,
            SLOT(currItemClicked(QTreeWidgetItem*, int)));

    MyTreeItem* Solver = new MyTreeItem();
    Solver->setText(0, "Solver");
    Solver->flag = flagStrings::solver;
    model->addChild(Solver);

    MyTreeItem* FieldSolver = new MyTreeItem();
    FieldSolver->setText(0, "Electric field Solver");
    FieldSolver->flag = flagStrings::fieldSolver;
    Solver->addChild(FieldSolver);

    MyTreeItem* SolverSettingsEmission = new MyTreeItem();
    SolverSettingsEmission->setText(0, "Emission model");
    SolverSettingsEmission->flag = flagStrings::solverEmission;
    Solver->addChild(SolverSettingsEmission);

    MyTreeItem* SolverSettings = new MyTreeItem();
    SolverSettings->setText(0, "General solvers Settings");
    SolverSettings->flag = flagStrings::solverSettings;
    Solver->addChild(SolverSettings);

    MyTreeItem* SimulatePic = new MyTreeItem();
    SimulatePic->setText(0, "PIC solver");
    SimulatePic->flag = flagStrings::simulatePIC;
    Solver->addChild(SimulatePic);

    /*MyTreeItem * SimulatePTI = new MyTreeItem();
    SimulatePTI->setText(0, "Iterative solver");
    SimulatePTI->flag = flagStrings::simulatePTI;
    Solver->addChild(SimulatePTI);*/

    MyTreeItem* Visuaization = new MyTreeItem();
    Visuaization->setText(0, "Visuaization elements");
    Visuaization->flag = flagStrings::Visuaization;
    Solver->addChild(Visuaization);

    MyTreeItem* Results = new MyTreeItem();
    Results->setText(0, "Results");
    Results->flag = flagStrings::results;
    model->addChild(Results);
};
