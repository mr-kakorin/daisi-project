#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GroupBoxWithItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "VisualizationFunctions.h"
#include "flagStringsTree.h"
#include "vtkComponent.h"
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>


const static std::vector<std::string> fieldParametersGlobal = {std::string("frequency, Hz"),
                                                               std::string("initial phase, rad")};

void MiddleWidget::ShowGlobalFieldMenu()
{
    std::vector<double> param = (*currentProject)->currentModel->GetglobalFieldConditions();
    groupBoxes.push_back(new GroupBoxWithItems("Global field parameters"));
    groupBoxes.back()->Create(fieldParametersGlobal, param);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    AddLastFictiveBox(1);
};

void MiddleWidget::PotentialFieldMenu(int currentConditioinIn)
{
    currentConditioin = currentConditioinIn;
    ShowSelectionMenu(currentConditioin, flagStrings::poisson, 0);

    /*QGroupBox *groupBox = new QGroupBox(tr("Set from file"));
    QVBoxLayout *vbox = new QVBoxLayout;
    groupBox->setMaximumHeight(150);
    QLabel* label = new QLabel();
    label->setObjectName(QStringLiteral("File path"));
    label->setText(QApplication::translate("MainWindow", "File path", 0));

    vbox->addWidget(label);

    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    vbox->addWidget(MiddleWidgetTextEdit1);

    QPushButton *pushBrouse = new QPushButton(tr("Browse"));
    pushBrouse->setMaximumSize(QSize(60, 40));
    vbox->addWidget(pushBrouse);

    QPushButton* set = new QPushButton(tr("Set"));
    set->setMaximumSize(QSize(50, 30));
    vbox->addWidget(set);


    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    middleWidgetGrid->addWidget((QWidget*)groupBox, 4, 0);

    QGroupBox * groupBox5 = new QGroupBox();
    groupBox5->setStyleSheet("border:0;");
    middleWidgetGrid->addWidget(groupBox5, 5, 0);
    setLayout(middleWidgetGrid);

    connect(set, SIGNAL(clicked()), this, SLOT(SetConditionFromFile()));
    connect(pushBrouse, SIGNAL(clicked()), this, SLOT(fileBoundaryBrouse()));*/
};

void MiddleWidget::AddCondition()
{
    /*QGroupBox *groupBox = new QGroupBox(tr("Add potential"));
    QVBoxLayout *vbox = new QVBoxLayout;

    mainMiddleWidgetButton = new QPushButton(tr("Add potential"));
    mainMiddleWidgetButton->setMaximumSize(QSize(80, 30));
    vbox->addWidget(mainMiddleWidgetButton);


    QLabel* label = new QLabel(groupBox);
    label->setObjectName(QStringLiteral("File path"));
    label->setText(QApplication::translate("MainWindow", "File path", 0));

    vbox->addWidget(label);

    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    vbox->addWidget(MiddleWidgetTextEdit1);

    QPushButton *pushBrouse = new QPushButton(tr("Browse"));
    pushBrouse->setMaximumSize(QSize(60, 40));
    vbox->addWidget(pushBrouse);

    MiddleWidgetButton1 = new QPushButton(tr("Create set of potentials"));
    MiddleWidgetButton1->setMaximumSize(QSize(150, 30));
    vbox->addWidget(MiddleWidgetButton1);


    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    middleWidgetGrid->addWidget((QWidget*)groupBox, 2, 0);

    connect(pushBrouse, SIGNAL(clicked()), this, SLOT(fileBoundaryBrouse()));

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(AddPotentialEvent()));
    connect(MiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(AddSetOfPotentialsEvent()));


    setLayout(middleWidgetGrid);*/

    mainMiddleWidgetButton = new QPushButton("Add potential");
    mainMiddleWidgetButton->setMaximumSize(QSize(100, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 0, 0);
    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(AddPotentialEvent()));
    AddLastFictiveBox(1);
};

void MiddleWidget::applyglobalFieldConditions(){
    //	(*currentProject)->currentModel->SetglobalFieldConditions(param);
};

/*void Daizy::AddSetOfPotentialsEvent()
{
(*currentProject)->currentModel->AddSetOfPotentials(MiddleWidgetTextEdit1->toPlainText().toStdString());

ShowProjectTreeSim();
};*/

/*void Daizy::SetConditionFromFile()
{
std::vector<int> tmp;
for (int i = 0; i < rigthList->count(); i++)
{
tmp.push_back(rigthList->item(i)->text().toInt());
};
(*currentProject)->currentModel->SetPropertyConditionsBoundariesList(currentClickItem->flagSearchBoundary,
currentClickItem->flagSearchBoundaryI, currentConditioin, tmp);

tmp.clear();
for (int i = 0; i < leftList->count(); i++)
{
tmp.push_back(leftList->item(i)->text().toInt());
};
(*currentProject)->currentModel->SetDefaultConditionsList(currentClickItem->flagSearchBoundary,
currentClickItem->flagSearchBoundaryI, tmp);


(*currentProject)->currentModel->SetConditionPropertiesFromFile(MiddleWidgetTextEdit1->toPlainText().toStdString(),
currentConditioin);
};*/