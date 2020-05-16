#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GeneralTools.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>

void MiddleWidget::ShowMeshMenu()
{
    std::vector<int> v2 = (*currentProject)->currentModel->GetBoundariesList();

    std::vector<int> v1 = (*currentProject)->currentModel->GetMeshBoundariesList();

    vectorSubtraction(v2, v1);

    ListsSelections.push_back(new ListSelectionMenu("Computation domain boundary"));
    ListsSelections.back()->Create(v2, v1, std::bind(&MiddleWidget::SetMeshBoundary, this),
                                   std::bind(&MiddleWidget::HighLightBoundary, this, std::placeholders::_1));
    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 0, 0);

    mainMiddleWidgetButton = new QPushButton(tr("Generate"));
    mainMiddleWidgetButton->setMaximumSize(QSize(70, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 1, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(GenerateMesh()));

    AddLastFictiveBox(2);
};
void MiddleWidget::SetMeshBoundary()
{
    std::vector<std::vector<int>> list1;
    std::vector<std::vector<int>> list2;
    GetListsData(list1, list2);

    (*currentProject)->currentModel->SetMeshBoundariesList(list2[0]);
}