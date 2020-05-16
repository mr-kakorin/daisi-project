#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"

#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QPushButton>
//#include "flagStrings.h"
void MiddleWidget::ShowDefaultCondMenu(int currentConditioin, QString name, std::string flag)
{
    std::vector<int> allBoundariesList = (*currentProject)->currentModel->GetBoundariesList();
    std::vector<int> DefaultConditionsList =
        (*currentProject)->currentModel->GetDefaultConditionsList(flag, currentConditioin);
    vectorSubtraction(allBoundariesList, DefaultConditionsList);

    ListsSelections.push_back(new ListSelectionMenu("Boundaries"));
    ListsSelections.back()->Create(allBoundariesList, DefaultConditionsList,
                                   std::bind(&MiddleWidget::ApplyDefaultBoundaries, this),
                                   std::bind(&MiddleWidget::HighLightBoundary, this, std::placeholders::_1));
    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 0, 0);

    (*VTKArray)[0]->setDataFunction(ShowBoundaries, (*currentProject)->currentModel, DefaultConditionsList,
                                    colors::redcolor, allBoundariesList, colors::blackcolor);
    (*VTKArray)[0]->refresh(0);

    AddLastFictiveBox(1);
};

void MiddleWidget::ShowSelectionMenu(int conditionNumber, std::string flag1, int flag2)
{

    std::vector<int> v3 = (*currentProject)->currentModel->GetBoundariesList();
    std::vector<int> v2 = (*currentProject)
                              ->currentModel->GetDefaultConditionsList(currentClickItem->flagSearchBoundary,
                                                                       currentClickItem->flagSearchBoundaryI);

    std::vector<int> v1 =
        (*currentProject)->currentModel->GetPropertyConditionsBoundariesList(flag1, flag2, conditionNumber);

    vectorSubtraction(v3, v1);

    int TypeFlag = (*currentProject)->currentModel->GetPropertyConditionTypeFlag(flag1, flag2, conditionNumber);

    int pos = 0;
    if (TypeFlag == 0)
    {

        ListsSelections.push_back(new ListSelectionMenu("Boundaries"));
        ListsSelections.back()->Create(v2, v1, std::bind(&MiddleWidget::applyBoundary, this),
                                       std::bind(&MiddleWidget::HighLightBoundary, this, std::placeholders::_1));
        middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), pos, 0);
        // listExec = std::bind(&Daizy::applyBoundary, this);
        pos++;

        std::vector<std::vector<int>> list1;
        std::vector<std::vector<int>> list2;

        GetListsData(list1, list2);

        (*VTKArray)[0]->setDataFunction(ShowBoundaries, (*currentProject)->currentModel,
                                        (*currentProject)->currentModel->GetBoundariesList(), colors::blackcolor,
                                        list2[0], colors::redcolor);
        //(*VTKArray)[0]->setDataFunction(ShowBoundaries, (*currentProject)->currentModel, v3, colors::blackcolor,
        // list2[0], colors::redcolor);

        (*VTKArray)[0]->refresh(0);
    }
    else
    {
        std::vector<double> manualRestictions =
            (*currentProject)->currentModel->GetPropertyConditionManualRestictions(flag1, flag2, conditionNumber);
        std::vector<std::string> radioBoundaryRestiction;

        switch ((*currentProject)->problemType)
        {
        case 1:
            radioBoundaryRestiction = {"X", "Y"};
            break;
        case 2:
            radioBoundaryRestiction = {"R", "Z"};
            break;
        case 3:
            radioBoundaryRestiction = {"R", "ZPhi"};
            break;
        case 4:
            radioBoundaryRestiction = {"X", "Y"};
            break;
        };

        groupBoxes.push_back(new GroupBoxWithItems("Position of boundary condition"));
        middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
        groupBoxes.back()->Create({"Plane position"}, radioBoundaryRestiction, {},
                                  std::vector<double>(manualRestictions.begin() + 1, manualRestictions.end()));
        pos++;
    }
    /*for (int i = 0; i < radioBoundaryRestiction.size(); i++)
    {
            vbox1->addWidget(radioBoundaryRestiction[i]);
            connect(radioBoundaryRestiction[i], SIGNAL(clicked()), this, SLOT(applyBoundary()));
    };

    middleWidgetGrid->addWidget(groupBox1, pos, 0);
    pos++;


    radioBoundaryRestiction[int(manualRestictions[0])]->setChecked(true);


    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));
    middleWidgetGrid->addWidget(MiddleWidgetTextEdit1, pos, 0);
    pos++;
    MiddleWidgetTextEdit1->setText(QString::number(manualRestictions[2]));
    connect(MiddleWidgetTextEdit1, SIGNAL(textChanged()), this, SLOT(applyBoundary()));



    QLabel* label = new QLabel();
    label->setText("Attach to conductor");
    middleWidgetGrid->addWidget(label, pos, 0);
    pos++;

    MiddleWidgetTextEdit2 = new MyQTextEdit();
    MiddleWidgetTextEdit2->setMaximumSize(QSize(16777215, 25));
    MiddleWidgetTextEdit2->setText(QString::number(manualRestictions[1]));

    middleWidgetGrid->addWidget(MiddleWidgetTextEdit2, pos, 0);
    connect(MiddleWidgetTextEdit2, SIGNAL(textChanged()), this, SLOT(applyBoundary()));

    pos++;

};*/
    std::vector<std::string> names =
        (*currentProject)->currentModel->GetConditionPropertiesNames(flag1, flag2, conditionNumber);

    if (names.size())
    {
        std::vector<double> props =
            (*currentProject)->currentModel->GetConditionProperties(flag1, flag2, conditionNumber);

        groupBoxes.push_back(new GroupBoxWithItems("Parameters"));
        middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);
        groupBoxes.back()->Create(names, props);
        pos++;
    }

    AddLastFictiveBox(pos);
};
void MiddleWidget::applyBoundary()
{

    bool ok = true;

    std::vector<std::vector<int>> list1;
    std::vector<std::vector<int>> list2;

    GetListsData(list1, list2);

    /*	if (MiddleWidgetTextEditVector[0]->toPlainText().toStdString() != std::string("no"))
            {
                    std::vector<double> props;
                    for (int i = 0; i < MiddleWidgetTextEditVector.size(); i++)
                            props.push_back(MiddleWidgetTextEditVector[i]->toDoubleMy(&ok));

                    (*currentProject)->currentModel->SetConditionProperties(currentClickItem->flagSearchBoundary,
       currentClickItem->flagSearchBoundaryI, currentConditioin, props);
            }*/

    int TypeFlag =
        (*currentProject)
            ->currentModel->GetPropertyConditionTypeFlag(currentClickItem->flagSearchBoundary,
                                                         currentClickItem->flagSearchBoundaryI, currentConditioin);

    if (TypeFlag == 0)
    {
        (*currentProject)
            ->currentModel->SetPropertyConditionsBoundariesList(currentClickItem->flagSearchBoundary,
                                                                currentClickItem->flagSearchBoundaryI,
                                                                currentConditioin, list2[0]);
        (*currentProject)
            ->currentModel->SetDefaultConditionsList(currentClickItem->flagSearchBoundary,
                                                     currentClickItem->flagSearchBoundaryI, list1[0]);
    }
    else
    {
        /*	std::vector<double> props(3);

                for (int i = 0; i<radioBoundaryRestiction.size(); i++)
                {
                        if (radioBoundaryRestiction[i]->isChecked())
                        {
                                props[0] = double(i);
                                break;
                        }
                };

                props[1] = MiddleWidgetTextEdit2->toDoubleMy(&ok);
                props[2] = MiddleWidgetTextEdit1->toDoubleMy(&ok);

                (*currentProject)->currentModel->SetPropertyConditionManualRestictions(currentClickItem->flagSearchBoundary,
           currentClickItem->flagSearchBoundaryI, currentConditioin, props);
                */
    };
};

void MiddleWidget::ApplyDefaultBoundaries()
{
    std::vector<std::vector<int>> list1;
    std::vector<std::vector<int>> list2;
    GetListsData(list1, list2);
    (*currentProject)
        ->currentModel->SetDefaultConditionsList(currentClickItem->flagSearchBoundary,
                                                 currentClickItem->flagSearchBoundaryI, list2[0]);
};

void MiddleWidget::HighLightBoundary(int n)
{
    std::vector<std::vector<int>> list1;
    std::vector<std::vector<int>> list2;

    GetListsData(list1, list2);
    (*VTKArray)[0]->setDataFunction(ShowBoundaries, (*currentProject)->currentModel,
                                    (*currentProject)->currentModel->GetBoundariesList(), colors::blackcolor, list2[0],
                                    colors::redcolor);
    (*VTKArray)[0]->refresh(0);
    (*VTKArray)[0]->setDataFunction(HighLigthBoundary, (*currentProject)->currentModel, n, colors::greencolor);
    (*VTKArray)[0]->refresh(0);
};
