#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "GroupBoxWithTextItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "flagStringsFlow.h"
#include "vtkComponent.h"
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QListWidgetItem>
#include <QtWidgets/QPushButton>

void MiddleWidget::ShowFlowEmitterPropertyAccelerator(int currentFlowIn) //
{
    int  pos = 0;
    bool ok  = true;

    std::vector<std::vector<double>> allParameters =
        (*currentProject)->currentModel->GetAllEmitterParameters(currentFlowIn);

    groupBoxes.push_back(new GroupBoxWithItems("Distribution"));
    groupBoxes.back()->Create({}, {}, {"Gauss", "Uniform"}, allParameters[0][0]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);

    pos++;

    groupBoxes.push_back(new GroupBoxWithItems("Number of particles"));
    groupBoxes.back()->Create(flagStrings::NamesOfNumbersParamsLinac, allParameters[1]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    ;
    pos++;

    groupBoxes.push_back(new GroupBoxWithItems("Distribution Parameters"));
    groupBoxes.back()->Create(flagStrings::NamesOfDistribParamsParamsLinac, allParameters[2]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    pos++;

    AddLastFictiveBox(pos);
}

void MiddleWidget::ShowFlowEmitterPropertyPlasma(int currentFlowIn)
{
    currentFlow = currentFlowIn;

    std::vector<int> v2 = (*currentProject)->currentModel->GetBoundariesList();

    std::vector<int> v1 = (*currentProject)->currentModel->GetEmitterBoundariesList(currentFlow);

    vectorSubtraction(v2, v1);

    int pos = 0;

    ListsSelections.push_back(new ListSelectionMenu("Emitter boundary"));
    ListsSelections.back()->Create(v2, v1, NULL,
                                   std::bind(&MiddleWidget::HighLightBoundary, this, std::placeholders::_1));
    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 0, 0);
    pos++;

    groupBoxes.push_back(new GroupBoxWithItems("Initialization parameters"));
    groupBoxes.back()->Create({"Number of discretization steps", "Position of supposed flow propagation point, x1",
                               "Position of supposed flow propagation point, x2",
                               "Position of supposed flow propagation point, x2"},
                              (*currentProject)->currentModel->GetEmitterInitParameters(currentFlow));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);

    pos++;
    QPushButton* mainMiddleWidgetButton = new QPushButton(tr("Create emitter"));
    mainMiddleWidgetButton->setMaximumSize(QSize(120, 30));

    /*if ((*currentProject)->problemType == 4)
    {
    std::vector<double> z = (*currentProject)->currentModel->GetAdditionalSourceInf(currentFlow);
    QLabel* labelZ1 = new QLabel();
    labelZ1->setText(QString("Z start"));
    middleWidgetGrid->addWidget(labelZ1, pos, 0);

    MiddleWidgetTextEdit5 = new MyQTextEdit();
    MiddleWidgetTextEdit5->setMaximumSize(QSize(16777215, 25));
    MiddleWidgetTextEdit5->setText(QString::number(z[0]));

    middleWidgetGrid->addWidget(MiddleWidgetTextEdit5, pos, 1);
    pos++;

    QLabel* labelZ2 = new QLabel();
    labelZ2->setText(QString("Z end"));
    middleWidgetGrid->addWidget(labelZ2, pos, 0);

    MiddleWidgetTextEdit6 = new MyQTextEdit();
    MiddleWidgetTextEdit6->setMaximumSize(QSize(16777215, 25));
    MiddleWidgetTextEdit6->setText(QString::number(z[1]));

    middleWidgetGrid->addWidget(MiddleWidgetTextEdit6, pos, 1);
    pos++;

    QLabel* labelZ3 = new QLabel();
    labelZ3->setText(QString("Number of particles along Z axis"));
    middleWidgetGrid->addWidget(labelZ3, pos, 0);

    MiddleWidgetTextEdit4 = new MyQTextEdit();
    MiddleWidgetTextEdit4->setMaximumSize(QSize(16777215, 25));
    MiddleWidgetTextEdit4->setText(QString::number(z[2]));

    middleWidgetGrid->addWidget(MiddleWidgetTextEdit4, pos, 1);
    pos++;
    };*/

    middleWidgetGrid->addWidget(mainMiddleWidgetButton, pos, 0);
    pos++;

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(ApplyEmitterBoundary()));

    if (!v1.size())
    {
        QLabel* labelM = new QLabel();
        labelM->setText("Emitter should be created before the further flow configuring");
        middleWidgetGrid->addWidget(labelM, pos, 0);
        return;
        pos++;
    };

    std::vector<std::vector<double>> allParameters =
        (*currentProject)->currentModel->GetAllEmitterParameters(currentFlow);

    double             current            = (*currentProject)->currentModel->GetEmissionCurrent(currentFlow);
    GroupBoxWithItems* MainparametersCalc = new GroupBoxWithItems("Calculated parameters");
    middleWidgetGrid->addWidget(MainparametersCalc->GetPointer(), pos, 0);
    MainparametersCalc->Create(std::vector<std::string>({"Emission current, A"}), std::vector<double>({current}), "labels");
    pos++;

    /*int current = (*currentProject)->currentModel->GetEmissionCurrent(currentFlow);
    groupBoxWithTextItems.push_back(new GroupBoxWithTextItems("Emission current"));
    groupBoxWithTextItems.back()->Create({ "Emission current, A" }, std::vector<std::string>{std::to_string(current) });
    middleWidgetGrid->addWidget(groupBoxWithTextItems.back()->GetPointer(), pos, 0);
    pos++;*/

    // groupBoxes.back()->SetRadioButtonsSignals(std::bind(&MiddleWidget::ChangeEmissionCurrentStyle, this));

    pos++;
    // ChangeEmissionCurrentStyle();

    groupBoxes.push_back(new GroupBoxWithItems("Distribution along emitter"));
    groupBoxes.back()->Create({}, {}, {"Quiet", "Random"}, allParameters[0][0]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);

    pos++;

    const static std::vector<std::string> NamesOfNumbersParams = {"Emit each timestep", "At each position",
                                                                  "Along emitter"};
    const static std::vector<std::string> NamesOfDistribParamsParams = {"Plasma temperature, eV",
                                                                        "Plasma density, cm^-3"};

    groupBoxes.push_back(new GroupBoxWithItems("Number of particles"));
    groupBoxes.back()->Create(NamesOfNumbersParams, allParameters[1]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    pos++;

    groupBoxes.push_back(new GroupBoxWithItems("Distribution Parameters"));
    groupBoxes.back()->Create(NamesOfDistribParamsParams, allParameters[2]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    pos++;

    AddLastFictiveBox(pos);
};

void MiddleWidget::ShowFlowEmitterPropertyEvaporation(int currentFlowIn)
{
    currentFlow = currentFlowIn;

    std::vector<std::vector<int>> vv2 = (*currentProject)->currentModel->GetConductorsList();
    std::vector<int>              v2;
    for (int i = 0; i < vv2.size(); i++)
        v2.push_back(i);
    std::vector<int> v1 = (*currentProject)->currentModel->GetEmitterBoundariesList(currentFlow);

    vectorSubtraction(v2, v1);

    int pos = 0;

    ListsSelections.push_back(new ListSelectionMenu("Select electrode for emission boundary"));
    ListsSelections.back()->Create(v2, v1, NULL, NULL);
    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 0, 0);

    pos++;

    groupBoxes.push_back(new GroupBoxWithItems("Initialization parameters"));
    groupBoxes.back()->Create({"Number of discretization steps", "Position of supposed flow propagation point, x1",
                               "Position of supposed flow propagation point, x2",
                               "Position of supposed flow propagation point, x2"},
                              (*currentProject)->currentModel->GetEmitterInitParameters(currentFlow));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);

    pos++;

    QPushButton* mainMiddleWidgetButton = new QPushButton(tr("Create emitter"));
    mainMiddleWidgetButton->setMaximumSize(QSize(120, 30));

    middleWidgetGrid->addWidget(mainMiddleWidgetButton, pos, 0);
    pos++;

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(ApplyEmitterBoundary()));

    if (!v1.size())
    {
        QLabel* labelM = new QLabel();
        labelM->setText("Emitter should be created before the further flow configuring");
        middleWidgetGrid->addWidget(labelM, pos, 0);
        return;
        pos++;
    };

    std::vector<std::vector<double>> allParameters =
        (*currentProject)->currentModel->GetAllEmitterParameters(currentFlow);

    double             current            = (*currentProject)->currentModel->GetEmissionCurrent(currentFlow);
    GroupBoxWithItems* MainparametersCalc = new GroupBoxWithItems("Calculated parameters");
    middleWidgetGrid->addWidget(MainparametersCalc->GetPointer(), pos, 0);
    MainparametersCalc->Create(std::vector<std::string>({"Emission current, A"}), std::vector<double>({current}), "labels");
    pos++;

    /*int current = (*currentProject)->currentModel->GetEmissionCurrent(currentFlow);
    groupBoxWithTextItems.push_back(new GroupBoxWithTextItems("Emission current"));
    groupBoxWithTextItems.back()->Create({ "Emission current, A" }, std::vector<std::string>{std::to_string(current) });
    middleWidgetGrid->addWidget(groupBoxWithTextItems.back()->GetPointer(), pos, 0);
    pos++;*/

    // groupBoxes.back()->SetRadioButtonsSignals(std::bind(&MiddleWidget::ChangeEmissionCurrentStyle, this));

    pos++;
    // ChangeEmissionCurrentStyle();

    groupBoxes.push_back(new GroupBoxWithItems("Distribution along emitter"));
    groupBoxes.back()->Create({}, {}, {"Quiet", "Random"}, allParameters[0][0]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);

    pos++;

    const static std::vector<std::string> NamesOfNumbersParams = {"Emit each timestep", "At each position",
                                                                  "Along emitter"};
    // const static std::vector<std::string> NamesOfDistribParamsParams = { "Thermal conductivity"};
    const static std::vector<std::string> NamesOfDistribParamsParams = {"Temperature of emitted particles"};
    groupBoxes.push_back(new GroupBoxWithItems("Number of particles"));
    groupBoxes.back()->Create(NamesOfNumbersParams, allParameters[1]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    pos++;

    groupBoxes.push_back(new GroupBoxWithItems("Distribution Parameters"));
    groupBoxes.back()->Create(NamesOfDistribParamsParams, allParameters[2]);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
    pos++;

    AddLastFictiveBox(pos);
};
void MiddleWidget::ShowFlowEmitterPropertyField(int currentFlowIn){

};
void MiddleWidget::ShowFlowEmitterPropertyThermionic(int currentFlowIn){

};
void MiddleWidget::ShowFlowEmitterPropertyPhotoemission(int currentFlowIn){

};
void MiddleWidget::ChangeEmissionCurrentStyle()
{
    std::vector<double> param;
    groupBoxes[1]->GetParameters(param);
    if (param[0] == 0)
        groupBoxes[0]->SetTextBoxEnable(0, false);
    if (param[0] == 1)
        groupBoxes[0]->SetTextBoxEnable(0, true);
};