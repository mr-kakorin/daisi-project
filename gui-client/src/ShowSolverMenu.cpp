#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "FlagStringsD.h"
#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "GroupBoxWithTextItems.h"
#include "MiddleWidget.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"
#include <QTimer>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>

void MiddleWidget::ShowEmissionModelSolverSettings()
{

    std::vector<std::vector<double>> param = (*currentProject)->currentModel->GetSolverEmissionModelParameters();

    if (param[1].size() == 0)
        return;

    groupBoxes.push_back(new GroupBoxWithItems("Space-charge limited emission model"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->Create({}, {}, flagStrings::SpChargeTypeNames, param[0][0], {0, 0});

    std::vector<std::string> strs = {"Emission cell length", "Emission cell height"};
    std::vector<double>      vals(2);
    for (int i = 0; i < param[1].size(); i++)
    {
        groupBoxes.push_back(
            new GroupBoxWithItems("Flow " + QString::number(param[1][i]) + " emission model parameters"));
        middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), i + 1, 0);
        vals[0] = param[2][i];
        vals[1] = param[3][i];
        groupBoxes.back()->Create(strs, vals);
    }

    AddLastFictiveBox(param[1].size() + 1);
}

void MiddleWidget::ApplySolverEmissionModelSettings(){
    /*bool ok = true;

    std::vector<std::vector<double>> param(4);

    for (int i = 0; i<2; i++)
    {
            if (radioSpChargeType[i]->isChecked())
            {
                    param[0].push_back(i);
                    break;
            }
    };


    for (int i = 0; i < MiddleWidgetTextEditVector.size()/2; i++)
    {
            param[2].push_back(MiddleWidgetTextEditVector[i * 2]->toDoubleMy(&ok));
            param[3].push_back(MiddleWidgetTextEditVector[i * 2 + 1]->toDoubleMy(&ok));
    }
    if (!ok)
    {
            QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }
    (*currentProject)->currentModel->SetSolverEmissionModelParameters(param);*/
};
void MiddleWidget::ShowSolverMenuPTI()
{
    std::vector<std::string>         names = {"omegaE 1", "omegaB 1", "omegaE 2", "omegaB 2"};
    std::vector<std::vector<double>> param = (*currentProject)->currentModel->GetSolverParametersPTI();

    groupBoxes.push_back(new GroupBoxWithItems("Relaxation parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->Create(names, param[0]);

    mainMiddleWidgetButton = new QPushButton(tr("Simulate dynamics"));
    mainMiddleWidgetButton->setMaximumSize(QSize(100, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 1, 0);
    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(SimulateIterative()));

    AddLastFictiveBox(2);
};

void MiddleWidget::ShowSolverMenuPIC()
{
    clearLayout(middleWidgetGrid);
    std::vector<double> param = (*currentProject)->currentModel->GetSolverParametersPIC()[0];

    QString str = "";
    for (int i = 0; i < param.size(); i++)
        str = str + QString::number(param[i]) + ";";

    groupBoxWithTextItems.push_back(new GroupBoxWithTextItems("Output time moments"));
    middleWidgetGrid->addWidget(groupBoxWithTextItems.back()->GetPointer(), 0, 0);
    groupBoxWithTextItems.back()->Create({""}, {str.toStdString()});

    groupBoxes.push_back(new GroupBoxWithItems("Visualization parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);
    groupBoxes.back()->Create({"Frame per second"}, std::vector<double>{5});

    mainMiddleWidgetButton = new QPushButton("Restart PIC simulations");
    mainMiddleWidgetButton->setMaximumSize(QSize(150, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 2, 0);
    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(SimulatePIC()));

    /*mainMiddleWidgetButton1 = new QPushButton("Continue PIC simulations");
    mainMiddleWidgetButton1->setMaximumSize(QSize(150, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton1, 2, 0);
    connect(mainMiddleWidgetButton1, SIGNAL(clicked()), this, SIGNAL(SimulatePICContinue()));*/

    AddLastFictiveBox(3);
};

void MiddleWidget::ShowFieldSolverMenu()
{
    bool ok = true;

    clearLayout(middleWidgetGrid);

    std::vector<double> param = (*currentProject)->currentModel->GetFieldSolverParameters();

    groupBoxes.push_back(new GroupBoxWithItems("Field solver parameters"));

    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);

    std::vector<std::string> paramFieldSolver = {"Omega External",       "Tolerance External", "Omega charge",
                                                 "Tolerance charge",     "Tolerance geometry", "Tolerance cut cells",
                                                 "Space-charge  margin", "Recalculate param"};

    // groupBoxes.back()->Create(paramFieldSolver, { "SOR" , "BiCGStab" }, { "Time depending flag","Space charge
    // flag","Fast init flag" }, param);

    groupBoxes.back()->Create(paramFieldSolver, param);

    mainMiddleWidgetButton = new QPushButton("Simulate fields");
    mainMiddleWidgetButton->setMaximumSize(QSize(150, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 1, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(FieldSimulate()));

    AddLastFictiveBox(2);
};

void MiddleWidget::ShowSolverSettings()
{

    std::vector<std::vector<double>> param = (*currentProject)->currentModel->GetSolverGeneralParameters();

    groupBoxes.push_back(new GroupBoxWithItems("Particle - mesh"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->Create({"boundaries tol", "Absorption fraction"}, flagStrings::ShapeTypeNames,
                              {"Remove non-CFL particles"}, param[0], {0, 0});

    groupBoxes.push_back(new GroupBoxWithItems("Time step estimation method"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);
    groupBoxes.back()->Create({"Steps per cell or per wave length"}, param[1], flagStrings::MoverTypeNames,
                              param[1].back(), {1, 0, 1});

    groupBoxes.push_back(new GroupBoxWithItems("Another settings"));
    std::vector<std::string> anotherNames = {"Save each step", "Trace save probability", "Number of threads",
                                             "Max particles per block"};
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 2, 0);
    groupBoxes.back()->Create(anotherNames, param[2]);

    AddLastFictiveBox(3);
};

void MiddleWidget::ApplySolverSettingsPTI(){
    /*bool ok = true;
    std::vector<std::vector<double>>  param;
    std::vector<double> p1;
    for (int i = 0; i < 4; i++)
    {
            p1.push_back(MiddleWidgetTextEditVector1[i]->toDoubleMy(&ok));
    };
    param.push_back(p1);
    (*currentProject)->currentModel->SetSolverPTIParameters(param);
    if (!ok)
    {
            QMessageBox::warning(this, "Daisi warning", "Non-digit input is converted to zero");
    }*/
};
void MiddleWidget::ApplySolverSettingsPIC()
{
    std::vector<double>      p1;
    std::string              error;
    std::vector<std::string> pp1;
    groupBoxWithTextItems[0]->GetParameters(pp1);

    char* pp = &pp1[0][0];
    strsplit(pp, ";", p1, error);

    if (error.size() != 0)
    {
        QMessageBox::warning(this, "Daisi warning", error.c_str());
        return;
    }

    std::vector<std::vector<double>> param;
    param.push_back(p1);
    (*currentProject)->currentModel->SetSolverParametersPIC(param);
};

/*void Daizy::ShowEstErrors()
{

        clearLayout(middleWidgetGrid);

        VTKArray[0]->setDataFunction(ShowErrors, (*currentProject)->currentModel);
        VTKArray[0]->refresh(0);

        mainMiddleWidgetButton = new QPushButton(tr("Error Estimate"));
        mainMiddleWidgetButton->setMaximumSize(QSize(100, 30));
        middleWidgetGrid->addWidget(mainMiddleWidgetButton, 0, 0);


        QGroupBox * groupBox1 = new QGroupBox();

        groupBox1->setStyleSheet("border:0;");

        middleWidgetGrid->addWidget(groupBox1, 1, 0);

        connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(ErrorEstimate()));
}*/

/*void Daizy::ShowGenerateGeometryMenu()
{
        std::vector<double> param = (*currentProject)->currentModel->GetAcceleratorParams();
        clearLayout(middleWidgetGrid);

        QGridLayout *vbox5 = new QGridLayout;
        QGroupBox *groupBox5 = new QGroupBox(tr("Accelerator parameters"));
        groupBox5->setMaximumHeight(500);
        groupBox5->setLayout(vbox5);

        middleWidgetGrid->addWidget((QWidget*)groupBox5, 0, 0);

        int s = flagStrings::LinacParamsNames.size();

        MiddleWidgetTextEditVector.resize(s);

        for (int i = 0; i < s; i++)
        {
                QLabel* label52 = new QLabel();
                label52->setText(flagStrings::LinacParamsNames[i].c_str());
                vbox5->addWidget(label52, i, 0);

                MiddleWidgetTextEditVector[i] = new MyQTextEdit();
                MiddleWidgetTextEditVector[i]->setMaximumSize(QSize(16777215, 25));
                MiddleWidgetTextEditVector[i]->setText(QString::number(param[i]));
                vbox5->addWidget(MiddleWidgetTextEditVector[i], i, 1);
        };


        QGridLayout *vbox11 = new QGridLayout;
        QGroupBox *groupBox11 = new QGroupBox(tr("Resonator type"));
        groupBox11->setMaximumHeight(150);
        groupBox11->setLayout(vbox11);

        middleWidgetGrid->addWidget((QWidget*)groupBox11, 1, 0);
        radioShapeType[0] = new QRadioButton("beta*lambda");
        radioShapeType[1] = new QRadioButton("beta*lambda/2");
        vbox11->addWidget(radioShapeType[0], 1, 0);
        vbox11->addWidget(radioShapeType[1], 2, 0);

        radioShapeType[int(param[s])]->setChecked(true);


        QGroupBox *groupBox1 = new QGroupBox(tr("Sync phase file"));
        QVBoxLayout *vbox1 = new QVBoxLayout;
        groupBox1->setMaximumHeight(100);


        MiddleWidgetTextEdit1 = new MyQTextEdit();
        MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

        vbox1->addWidget(MiddleWidgetTextEdit1);

        QPushButton *pushBrouse = new QPushButton(tr("Browse"));
        pushBrouse->setMaximumSize(QSize(60, 40));
        vbox1->addWidget(pushBrouse);

        groupBox1->setLayout(vbox1);

        middleWidgetGrid->addWidget((QWidget*)groupBox1, 2, 0);


        mainMiddleWidgetButton = new QPushButton(tr("Save Parameters"));
        mainMiddleWidgetButton->setMaximumSize(QSize(150, 30));
        middleWidgetGrid->addWidget(mainMiddleWidgetButton, 3, 0);

        QPushButton* mainMiddleWidgetButton1 = new QPushButton(tr("Apply (Generate Geometry)"));
        mainMiddleWidgetButton1->setMaximumSize(QSize(150, 30));
        middleWidgetGrid->addWidget(mainMiddleWidgetButton1, 4, 0);
        


        QGroupBox * groupBox2 = new QGroupBox();

        groupBox2->setStyleSheet("border:0;");

        middleWidgetGrid->addWidget(groupBox2, 5, 0);

        connect(pushBrouse, SIGNAL(clicked()), this, SLOT(fileBoundaryBrouse()));
        connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(SaveParams()));
        connect(mainMiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(GenerateGeometry()));

}*/