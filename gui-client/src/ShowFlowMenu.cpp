#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QListWidgetItem>
#include <QtWidgets/QPushButton>

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "flagStringsFlow.h"
#include "vtkComponent.h"
//#include "flagStrings.h"
const static std::vector<std::string> FlowParametersNames = {"Beam cloud", "Current density distribution",
                                                             "Electric field on emitter"};

void MiddleWidget::ShowAddFlowMenu()
{

    groupBoxes.push_back(new GroupBoxWithItems("Particles type"));
    groupBoxes.back()->Create({}, {}, flagStrings::radioParticleTypeNames, 1);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->SetRadioButtonsSignals(std::bind(&MiddleWidget::ChangeFlowDistrStyle, this));

    groupBoxes.push_back(new GroupBoxWithItems("Emission type"));
    groupBoxes.back()->Create({}, {}, flagStrings::radioDistributionStyleNames, 0, {1, 1, 0, 0, 0, 1});
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Mass and Charge"));
    groupBoxes.back()->Create({"Mass/Proton mass", "Charge/Proton charge"}, std::vector<double>{1, 1});
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 2, 0);

    mainMiddleWidgetButton = new QPushButton(tr("Add"));
    mainMiddleWidgetButton->setMaximumSize(QSize(70, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 3, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(AddFlow()));

    AddLastFictiveBox(4);
};

void MiddleWidget::ShowFlowSummary(int i)
{
    currentFlow              = i;
    std::vector<double> MC   = (*currentProject)->currentModel->GetFlowMCNumbers(i);
    std::vector<double> prop = (*currentProject)->currentModel->GetFlowProperties(i);

    groupBoxes.push_back(new GroupBoxWithItems("Current properties"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->Create({"Emission type: ", "Particle mass: ", "Particle charge: "},
                              {flagStrings::radioDistributionStyleNames[prop[1]],
                               QString::number(prop[2]).toStdString(), QString::number(prop[3]).toStdString()},
                              "str");

    groupBoxes.push_back(new GroupBoxWithItems("Change properties"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);
    groupBoxes.back()->Create({"Mass number", "Charge number", "Maximal time, ns"},
                              std::vector<double>{MC[0], MC[1], prop[4]});

    AddLastFictiveBox(2);
};

void MiddleWidget::ShowFlowEmitterProperty(int currentFlowIn)
{
    std::vector<double> prop              = (*currentProject)->currentModel->GetFlowProperties(currentFlowIn);
    int                 distributionStyle = (int)prop[1];
    switch (distributionStyle)
    {
    case 0:
        ShowFlowEmitterPropertyPlasma(currentFlowIn);
        break;
    case 1:
        ShowFlowEmitterPropertyEvaporation(currentFlowIn);
        break;
    case 2:
        ShowFlowEmitterPropertyField(currentFlowIn);
        break;
    case 3:
        ShowFlowEmitterPropertyThermionic(currentFlowIn);
        break;
    case 4:
        ShowFlowEmitterPropertyPhotoemission(currentFlowIn);
        break;
    case 5:
        ShowFlowEmitterPropertyAccelerator(currentFlowIn);
        break;
    }
}

/*void Daizy::SetDirectionPoints()//
{
        bool ok = true;
        std::vector<double> startPoint(3);
        for (int i = 0; i < 3; i++)
        {
                startPoint[i]=plotPoint1[i]->toDoubleMy(&ok);
        }
        std::vector<double> endPoint(3);
        for (int i = 0; i < 3; i++)
        {
                endPoint[i]=plotPoint2[i]->toDoubleMy(&ok);
        }
        (*currentProject)->currentModel->SetDirectionPoints(currentFlow,startPoint,endPoint);
        if (!ok)
        {
                QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
        }
}*/
void MiddleWidget::SetEmissionParametersLinac()
{
    /*bool ok = true;

    int n = int(ParticleNumberTextEditVector.size());
    std::vector<int > paramNumber(20);
    paramNumber[0] = 1;
    for (int j = 1; j < n+1; j++)
    {
            paramNumber[j]=ParticleNumberTextEditVector[j-1]->toPlainText().toInt();
    };
    (*currentProject)->currentModel->SetParticlesNumber(currentFlow, paramNumber);

    int n1 = int(ParticleDistributionParametersTextEditVector.size());
    std::vector<double > paramDistrib(20);
    for (int j = 0; j < n1; j++)
    {
            paramDistrib[j] = ParticleDistributionParametersTextEditVector[j]->toDoubleMy(&ok);
    };
    (*currentProject)->currentModel->SetDistributionParameters(currentFlow, paramDistrib);*/
}

void MiddleWidget::ApplyEmitterBoundary()
{
    /*if ((*currentProject)->problemType == 4)
    {
            std::vector<double> z;
            z.push_back(MiddleWidgetTextEdit5->toDoubleMy(&ok));
            z.push_back(MiddleWidgetTextEdit6->toDoubleMy(&ok));
            z.push_back(MiddleWidgetTextEdit4->toDoubleMy(&ok));
            (*currentProject)->currentModel->SetAdditionalSourceInf(currentFlow, z);
    }*/
    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;
    std::string                           errorMessage;

    bool ok = FetchParameters(parameters1, parameters2, parameters3);

    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Incorrect non-digit input");
        return;
    }

    std::vector<std::vector<int>> list1;
    std::vector<std::vector<int>> list2;
    GetListsData(list1, list2);
    std::string error;

    if (list2[0].size())
    {
        (*currentProject)->currentModel->SetEmitterBoundariesList(currentFlow, list2[0], parameters1[0], error);
        clear();
        ShowFlowEmitterProperty(currentFlow);
        if (error.size())
            QMessageBox::critical(this, "Daisi error", error.c_str());
    }
    else
        QMessageBox::critical(this, "Daisi error", "Error! The boundary list is empty!");
};
void MiddleWidget::ShowFlowState(int flow)
{
    bool ok                 = true;
    currentFlow             = flow;
    QListWidget*     List   = new QListWidget();
    std::vector<int> styles = (*currentProject)->currentModel->GetNumberParticlesFlowsTypes();

    std::vector<std::string> FlowParametersNamesCurrent;
    if (styles[flow] == 5)
        FlowParametersNamesCurrent =
            std::vector<std::string>(FlowParametersNames.begin(), FlowParametersNames.begin() + 1);
    else
        FlowParametersNamesCurrent = FlowParametersNames;

    for (int i = 0; i < FlowParametersNamesCurrent.size(); i++)
    {
        QListWidgetItem* newitem = new QListWidgetItem();
        newitem->setText(FlowParametersNamesCurrent[i].c_str());
        List->addItem(newitem);
    };

    middleWidgetGrid->addWidget(List, 0, 0);

    connect(List, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listShowFlowStateClick(QListWidgetItem*)));

    if (styles[flow] < 5)
        return;

    QTreeWidget* TreeList = new QTreeWidget();
    TreeList->setHeaderLabel("Emittances characteristics");

    std::vector<double>                   Em = (*currentProject)->currentModel->GetEmittancesList(flow);
    std::vector<std::vector<MyTreeItem*>> itemvector(Em.size());

    std::vector<std::string> currentNames = {"Longitudinal", "XdX", "YdY", "XY", "W Spectrum", "Phi Spectrum"};

    for (int k = 0; k < Em.size(); k++)
    {
        MyTreeItem* EmRes = new MyTreeItem();
        EmRes->setText(0, QString::number(Em[k]));
        TreeList->addTopLevelItem(EmRes);
        EmRes->flag1 = -1;

        for (int i = 0; i < currentNames.size(); i++)
        {
            itemvector[k].push_back(new MyTreeItem());
            itemvector[k].back()->setText(0, currentNames[i].c_str());
            itemvector[k].back()->flag1 = flow;
            itemvector[k].back()->flag2 = k;
            itemvector[k].back()->flag4 = i;

            EmRes->addChild(itemvector[k].back());
        };
    }

    middleWidgetGrid->addWidget(TreeList, 1, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Emittance monitor properties"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 2, 0);
    groupBoxes.back()->Create({"Position, for new monitor Z = "}, std::vector<double>{0.01});

    mainMiddleWidgetButton = new QPushButton(tr("Add emittance monitor"));
    mainMiddleWidgetButton->setMaximumSize(QSize(140, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 3, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Histograms properties"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 4, 0);
    groupBoxes.back()->Create({"Number of colums in histograms"}, std::vector<double>{20});

    connect(TreeList, SIGNAL(itemClicked(QTreeWidgetItem*, int)), this,
            SLOT(listShowResultClickEmittances(QTreeWidgetItem*, int)));

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(AddEmittance()));
};
void MiddleWidget::AddEmittance()
{
    std::vector<double> p;
    bool                ok = groupBoxes[0]->GetParameters(p);

    if (!ok)
    {
        QMessageBox::critical(this, "Daisi warning", "Incorrect input");
        return;
    };

    (*currentProject)->currentModel->AddEmittance(currentFlow, p[0]);
};
void MiddleWidget::listShowResultClickEmittances(QTreeWidgetItem* item, int)
{
    std::vector<double> p;
    bool                ok = groupBoxes[1]->GetParameters(p);

    if (!ok)
    {
        QMessageBox::critical(this, "Daisi warning", "Incorrect input");
        return;
    };

    MyTreeItem* item1 = dynamic_cast<MyTreeItem*>(item);
    if (item1->flag1 == -1)
        return;
    (*VTKArray)[0]->setDataFunction(ShowEmittanceDataPlot, (*currentProject)->currentModel, item1->flag1, item1->flag2,
                                    item1->flag4, p[0]);
    (*VTKArray)[0]->refresh(0);
    return;
}

void MiddleWidget::listShowFlowStateClick(QListWidgetItem* item)
{
    //	Visualization_thread = std::thread(&Daizy::VisualThr, this, item);

    std::string name = item->text().toStdString();
    if (name == FlowParametersNames[0])
    {
        (*VTKArray)[0]->setDataFunction(ShowParticlesCloud, (*currentProject)->currentModel, currentFlow,
                                        (*currentProject)->precisionType, (*currentProject)->problemType);
        (*VTKArray)[0]->refresh(0);
        return;
    };
    if (name == FlowParametersNames[1])
    {

        (*VTKArray)[0]->setDataFunction(ShowCurrentDensity, (*currentProject)->currentModel, currentFlow);
        (*VTKArray)[0]->refresh(0);
        return;
    };
    if (name == FlowParametersNames[2])
    {

        (*VTKArray)[0]->setDataFunction(ShowEmitterField, (*currentProject)->currentModel, currentFlow);
        (*VTKArray)[0]->refresh(0);
        return;
    };
};

void MiddleWidget::ChangeFlowDistrStyle()
{
    std::vector<double> param;
    groupBoxes[0]->GetParameters(param);
    if (param[0] == 0)
        groupBoxes[2]->SetParameters({5.44616975512e-04, -1});
    else
        groupBoxes[2]->SetParameters({1, 1});
};