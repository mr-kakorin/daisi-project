#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"
std::vector<int>    Prevlist2;
std::vector<double> Prevlist3;
//#include "flagStrings.h"
std::vector<QString> names = {"Average power density", "Average irradiated current density",
                              "Average collected current density"};
void MiddleWidget::ShowConductorSelectionMenu(int number)
{

    QListWidget* List1 = new QListWidget();

    QListWidgetItem* newitem = new QListWidgetItem();
    newitem->setText(names[0]);
    List1->addItem(newitem);

    QListWidgetItem* newitem1 = new QListWidgetItem();
    newitem1->setText(names[1]);
    List1->addItem(newitem1);

    QListWidgetItem* newitem2 = new QListWidgetItem();
    newitem2->setText(names[2]);
    List1->addItem(newitem2);

    middleWidgetGrid->addWidget(List1, 0, 0);

    connect(List1, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listShowConductorStateClick(QListWidgetItem*)));

    std::vector<std::vector<int>> List = (*currentProject)->currentModel->GetConductorsList();
    ;
    ;

    std::vector<int> v2 = (*currentProject)->currentModel->GetBoundariesList();

    std::vector<int> v1 = List[number];

    vectorSubtraction(v2, v1);

    /*for (int i = 0; i < v1.size(); i++)
    v2[v1[i]] = -1;

    int s = int(v2.size());

    for (int i = 0; i < s; i++)
    {
    if (v2[i] == -1)
    {
    v2.erase(v2.begin() + i);
    i--;
    s--;
    }
    }*/

    Prevlist2 = v1;

    ListsSelections.push_back(new ListSelectionMenu("Conductor boundary"));
    // ListsSelections.back()->Create(v2, v1, std::bind(&MiddleWidget::ApplyConductorBoundaries, this),
    // std::bind(&MiddleWidget::HighLightBoundary, this, std::placeholders::_1));
    ListsSelections.back()->Create(v2, v1, NULL,
                                   std::bind(&MiddleWidget::HighLightBoundary, this, std::placeholders::_1));

    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 1, 0);

    std::vector<double> params = (*currentProject)->currentModel->GetElectrodeParametersList(number);
    Prevlist3                  = params;

    groupBoxes.push_back(new GroupBoxWithItems("Parameters"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 2, 0);
    groupBoxes.back()->Create({"Discretization length, m", "Averaging time, ns"}, params);

    (*VTKArray)[0]->setDataFunction(ShowValueAlongConductor, (*currentProject)->currentModel, number, 0);
    (*VTKArray)[0]->refresh(0);

    AddLastFictiveBox(3);
}

void MiddleWidget::listShowConductorStateClick(QListWidgetItem* item)
{
    //	Visualization_thread = std::thread(&Daizy::VisualThr, this, item);

    int i = 0;
    for (i = 0; i < names.size(); i++)
    {
        if (item->text() == names[i])
            break;
    }

    (*VTKArray)[0]->setDataFunction(ShowValueAlongConductor, (*currentProject)->currentModel, currentClickItem->flag1,
                                    i);
    (*VTKArray)[0]->refresh(0);
};

void MiddleWidget::ApplyConductorBoundaries(int n)
{
    std::vector<std::vector<int>> list1;
    std::vector<std::vector<int>> list2;
    GetListsData(list1, list2);
    std::vector<double> parameters;
    groupBoxes[0]->GetParameters(parameters);
    std::string error;
    if (Prevlist2 != list2[0] || Prevlist3 != parameters)
    {
        (*currentProject)->currentModel->setElectrodeParametersList(n, parameters);
        (*currentProject)->currentModel->SetConductorsList(list2[0], n, parameters[0], error);
    }
    if (error.size())
        QMessageBox::critical(this, "Daisi error", error.c_str());
};

void MiddleWidget::AddConductor()
{

    mainMiddleWidgetButton = new QPushButton("Add conductor");
    mainMiddleWidgetButton->setMaximumSize(QSize(80, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 0, 0);
    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(AddConductorEvent()));
    AddLastFictiveBox(1);
};