#include "DaisiClient.h"
#include "FlagStringsD.h"
#include "GroupBoxWithItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>

void MiddleWidget::AddFlowCondition()
{

    groupBoxes.push_back(new GroupBoxWithItems("Boundary condition type"));
    groupBoxes.back()->Create({}, {}, flagStrings::flowBoundaryTypeNames, 1);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Restriction type"));
    groupBoxes.back()->Create({}, {}, flagStrings::flowBoundaryTypeFlagsNames, 0);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);

    mainMiddleWidgetButton = new QPushButton(tr("Add condition"));
    mainMiddleWidgetButton->setMaximumSize(QSize(200, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 2, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(AddFlowConditionEvent()));
    AddLastFictiveBox(3);
};
