#include "GroupBoxWithItems.h"
#include "MiddleWidget.h"

#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>

int  currentFlow;
void MiddleWidget::ShowAddFlowAccelMenu()
{
    mainMiddleWidgetButton = new QPushButton(tr("Add"));
    mainMiddleWidgetButton->setMaximumSize(QSize(70, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 0, 0);

    AddLastFictiveBox(1);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(AddFlowAccel()));
};
void MiddleWidget::ShowFlowAccel(int flow)
{

    std::vector<double>      params;
    std::vector<std::string> keys;

    (*currentProject)->accelModel->GetParametersAccelFlow(keys, params, flow);

    groupBoxes.push_back(new GroupBoxWithItems("Flow parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    std::vector<std::string> radioButtonsNames = {"Gauss initial distr", "Uniform initial distr"};
    groupBoxes.back()->Create(std::vector<std::string>(keys.begin(), keys.end() - 1),
                              std::vector<double>(params.begin(), params.end() - 1), radioButtonsNames, params.back());

    AddLastFictiveBox(1);
};