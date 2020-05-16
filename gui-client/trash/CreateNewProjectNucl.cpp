#include "DaisiClient.h"
#include "FlagStringsD.h"

void Daizy::ShowCreateNewProject()
{

    DefaultBuilder();

    clearLayout(middleWidgetGrid);

    QVBoxLayout* vbox     = new QVBoxLayout;
    QGroupBox*   groupBox = new QGroupBox(tr("Problem type"));

    groupBox->setMaximumHeight(300);
    middleWidgetGrid->addWidget((QWidget*)groupBox, 0, 0);
    for (int i = 0; i < flagStrings::problemTypesNames.size(); i++)
    {
        radioProblemType[i] = new QRadioButton(tr(flagStrings::problemTypesNames[i].c_str()));
        vbox->addWidget(radioProblemType[i]);
        radioProblemType[i]->setEnabled(false);
    }
    radioProblemType[9]->setChecked(true);

    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    QVBoxLayout* vbox1     = new QVBoxLayout;
    QGroupBox*   groupBox1 = new QGroupBox(tr("Precision"));
    groupBox1->setMaximumHeight(80);
    middleWidgetGrid->addWidget((QWidget*)groupBox1, 1, 0);

    radioPrecisionType[0] = new QRadioButton(tr(flagStrings::precisionNames[0].c_str()));
    radioPrecisionType[1] = new QRadioButton(tr(flagStrings::precisionNames[1].c_str()));

    radioPrecisionType[0]->setEnabled(false);
    radioPrecisionType[1]->setEnabled(false);

    radioPrecisionType[0]->setChecked(true);
    vbox1->addWidget(radioPrecisionType[0]);
    vbox1->addWidget(radioPrecisionType[1]);
    vbox1->addStretch(1);
    groupBox1->setLayout(vbox1);

    QVBoxLayout* vbox11     = new QVBoxLayout;
    QGroupBox*   groupBox11 = new QGroupBox(tr("Solver type"));
    groupBox11->setMaximumHeight(120);
    middleWidgetGrid->addWidget((QWidget*)groupBox11, 2, 0);

    radioSolverType[0] = new QRadioButton(tr(flagStrings::solverTypeNames[0].c_str()));
    radioSolverType[1] = new QRadioButton(tr(flagStrings::solverTypeNames[1].c_str()));
    radioSolverType[2] = new QRadioButton(tr(flagStrings::solverTypeNames[2].c_str()));

    radioSolverType[0]->setEnabled(false);
    radioSolverType[1]->setEnabled(false);
    radioSolverType[2]->setEnabled(false);

    // radioSolverType[0]->setChecked(true);
    vbox11->addWidget(radioSolverType[0]);
    vbox11->addWidget(radioSolverType[1]);
    vbox11->addWidget(radioSolverType[2]);
    vbox11->addStretch(1);
    groupBox11->setLayout(vbox11);

    QVBoxLayout* vbox2     = new QVBoxLayout;
    QGroupBox*   groupBox2 = new QGroupBox(tr("Project name"));
    groupBox2->setMaximumHeight(150);
    middleWidgetGrid->addWidget(groupBox2, 3, 0);
    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setObjectName(QStringLiteral("textEdit"));
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    QLabel* label = new QLabel(groupBox2);
    label->setObjectName(QStringLiteral("Project path"));
    label->setText(QApplication::translate("MainWindow", "Project path", 0));
    QPushButton* push = new QPushButton(tr("Browse"));
    push->setMaximumSize(QSize(60, 40));

    MiddleWidgetTextEdit2 = new MyQTextEdit();
    MiddleWidgetTextEdit2->setObjectName(QStringLiteral("textEdit"));
    MiddleWidgetTextEdit2->setMaximumSize(QSize(16777215, 25));
    MiddleWidgetTextEdit2->setText(QStringLiteral("..\\Projects"));

    vbox2->addWidget(MiddleWidgetTextEdit1);
    vbox2->addWidget(label);
    vbox2->addWidget(MiddleWidgetTextEdit2);
    vbox2->addWidget(push);

    vbox2->addStretch(1);
    groupBox2->setLayout(vbox2);

    mainMiddleWidgetButton = new QPushButton(tr("Create"));
    mainMiddleWidgetButton->setMaximumSize(QSize(50, 30));
    middleWidgetGrid->addWidget((QWidget*)mainMiddleWidgetButton, 4, 0);
    setLayout(middleWidgetGrid);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(CreateNewProjectEvent()));
};
void Daizy::CreateNewProjectEvent()
{
    for (int i = 0; i < 2; i++)
    {
        if (radioPrecisionType[i]->isChecked())
        {
            precisionType = i;
            break;
        }
    };
    for (int i = 0; i < flagStrings::problemTypesNames.size(); i++)
    {
        if (radioProblemType[i]->isChecked())
        {
            problemType = i;
            break;
        }
    };
    for (int i = 0; i < 3; i++)
    {
        if (radioSolverType[i]->isChecked())
        {
            solverType = i;
            break;
        }
    };
    currentProject = new Dproject::project(problemType, precisionType, solverType,
                                           MiddleWidgetTextEdit1->toPlainText().toStdString(),
                                           MiddleWidgetTextEdit2->toPlainText().toStdString());
    ShowProjectTreeAccel();
    showSummary();
};
