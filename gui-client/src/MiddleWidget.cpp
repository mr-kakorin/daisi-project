#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "MiddleWidget.h"
#include "FilesBrouses.h"
#include "GeneralTools.h"
#include "GroupBoxWithItems.h"
#include "GroupBoxWithTextItems.h"
#include "ListSelectionMenu.h"
#include "MyTreeItem.h"
#include "regressionProcessor.h"
#include "vtkComponent.h"

#include <QFileDialog>
#include <QRadioButton>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListWidgetItem>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QVBoxLayout>

#ifdef NUCL
#define flagProblems 2;
#endif

#ifdef SIM
#define flagProblems 1;
#ifdef NUCL
#define flagProblems 2;
#endif
#endif

const static std::vector<std::string> problemTypesNames = {"1d",
                                                           "2d Cartesian",
                                                           "2d Cylindrical axisymmetric",
                                                           "2d Polar",
                                                           "3d Cartesian from extrusion",
                                                           "3d Cartesian",
                                                           "3d Cylindrical",
                                                           "RFQ Design",
                                                           "DTL Design",
                                                           "Dynamics in synchrotron optics"};
const static std::vector<std::string> precisionNames = {"double", "float"};

void MiddleWidget::InitTrees(int numberOfResultPlots)
{
    TreeList.resize(numberOfResultPlots);
    for (int i = 0; i < numberOfResultPlots; i++)
    {
        TreeList[i] = new QTreeWidget();
        TreeList[i]->setHeaderLabel("Simulations results");
        middleWidgetGrid->addWidget(TreeList[i], i, 0);
    }
    itemvector.resize(numberOfResultPlots);
};

QGridLayout* MiddleWidget::GetMiddleWidgetGrid()
{
    return middleWidgetGrid;
};

void MiddleWidget::GetListsData(std::vector<std::vector<int>>& list1, std::vector<std::vector<int>>& list2)
{
    list1.resize(ListsSelections.size());
    list2.resize(ListsSelections.size());

    for (int i = 0; i < ListsSelections.size(); i++)
    {
        ListsSelections[i]->GetListsData(list1[i], list2[i]);
    };
};
void MiddleWidget::GetListsData(std::vector<std::vector<std::string>>& list1,
                                std::vector<std::vector<std::string>>& list2)
{
    list1.resize(ListsSelections.size());
    list2.resize(ListsSelections.size());

    for (int i = 0; i < ListsSelections.size(); i++)
    {
        ListsSelections[i]->GetListsData(list1[i], list2[i]);
    };
};

void MiddleWidget::clear()
{
    for (int i = 0; i < groupBoxes.size(); i++)
        delete groupBoxes[i];

    for (int i = 0; i < filesBrouses.size(); i++)
        delete filesBrouses[i];

    for (int i = 0; i < groupBoxWithTextItems.size(); i++)
        delete groupBoxWithTextItems[i];

    for (int i = 0; i < TreeList.size(); i++)
        delete TreeList[i];

    for (int i = 0; i < ListsSelections.size(); i++)
        delete ListsSelections[i];

    for (int i = 0; i < itemvector.size(); i++)
    {
        //	for (int j = 0; j < itemvector[i].size(); j++)
        //		delete itemvector[i][j];
        itemvector[i].clear();
    }

    TreeList.clear();
    itemvector.clear();
    groupBoxes.clear();
    filesBrouses.clear();
    ListsSelections.clear();
    groupBoxWithTextItems.clear();

    clearLayout(middleWidgetGrid);
}
void MiddleWidget::AddLastFictiveBox(int pos)
{
    groupBoxFictive = new QGroupBox();
    groupBoxFictive->setStyleSheet("border:0;");
    middleWidgetGrid->addWidget(groupBoxFictive, pos, 0);
};

bool MiddleWidget::FetchParameters(std::vector<std::vector<double>>&      parameters1,
                                   std::vector<std::vector<std::string>>& parameters2,
                                   std::vector<std::vector<std::string>>& parameters3)
{
    bool ok = true;
    parameters1.resize(groupBoxes.size());

    if (!groupBoxes.size())
        parameters1.resize(1);

    for (int i = 0; i < groupBoxes.size(); i++)
    {
        if (!groupBoxes[i]->GetParameters(parameters1[i]))
            ok = false;
    }

    parameters2.resize(groupBoxWithTextItems.size());

    if (!groupBoxWithTextItems.size())
        parameters2.resize(1);

    for (int i = 0; i < groupBoxWithTextItems.size(); i++)
        groupBoxWithTextItems[i]->GetParameters(parameters2[i]);

    parameters3.resize(filesBrouses.size());

    if (!filesBrouses.size())
        parameters3.resize(1);

    for (int i = 0; i < filesBrouses.size(); i++)
        filesBrouses[i]->GetParameters(parameters3[i]);

    return ok;
};

MiddleWidget::MiddleWidget(){

};
MiddleWidget::MiddleWidget(QWidget* widget, Dproject::project** currentProjectIn,
                           std::vector<vtkComponent*>* VTKArrayIn)
{
    middleWidgetGrid = new QGridLayout(widget);
    currentProject   = currentProjectIn;
    VTKArray         = VTKArrayIn;
    groupBoxFictive  = new QGroupBox();
    groupBoxFictive->setStyleSheet("border:0;");
};

void MiddleWidget::ShowCreateNewProject()
{
    clear();
    emit DefaultBuilder();
    emit ResetFlags();

    int              flag = flagProblems;
    std::vector<int> flags;

    if (0 == flag)
        flags = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    if (1 == flag)
        flags = {0, 1, 1, 0, 0, 0, 0, 0, 0, 0};

    if (2 == flag)
        flags = {0, 1, 1, 0, 0, 0, 0, 1, 1, 1};

    groupBoxes.push_back(new GroupBoxWithItems("Problem type"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->Create({}, {problemTypesNames}, {}, {9}, flags);

    groupBoxes.push_back(new GroupBoxWithItems("Precision type"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);
    groupBoxes.back()->Create({}, {precisionNames}, {}, {0});

    groupBoxWithTextItems.push_back(new GroupBoxWithTextItems("Path and name"));
    middleWidgetGrid->addWidget(groupBoxWithTextItems.back()->GetPointer(), 2, 0);
    groupBoxWithTextItems.back()->Create({"Project path", "Project name"}, {"../Projects", ""});

    mainMiddleWidgetButton = new QPushButton("Create");
    mainMiddleWidgetButton->setMaximumSize(QSize(50, 30));
    middleWidgetGrid->addWidget((QWidget*)mainMiddleWidgetButton, 3, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SIGNAL(CreateNewProject()));
    AddLastFictiveBox(4);
};

void MiddleWidget::showSummary()
{

    int n1 = (*currentProject)->problemType;
    int n2 = (*currentProject)->precisionType;

    groupBoxes.push_back(new GroupBoxWithItems("Project properties"));
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 0, 0);
    groupBoxes.back()->Create({"Problem type", "Precision type"}, {problemTypesNames[n1], precisionNames[n2]}, "str");
    AddLastFictiveBox(1);
};