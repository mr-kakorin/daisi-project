#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "FlagStringsD.h"
#include "GroupBoxWithItems.h"
#include "ListSelectionMenu.h"
#include "MiddleWidget.h"
#include "MyTreeItem.h"
#include "VisualizationFunctions.h"
#include "flagStringsResults.h"
#include "flagStringsTree.h"
#include "vtkComponent.h"

#include <QVTKWidget.h>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QListWidgetItem>

void MiddleWidget::ShowAddPlotMenu()
{

    std::vector<std::shared_ptr<lineplot>> plotsData = (*currentProject)->currentModel->GetPlotsVector();

    std::vector<std::string> PlotNames;
    std::vector<std::string> NowPlotNames = (*currentProject)->currentModel->GetPlotNames();

    std::vector<std::string> SimNames = (*currentProject)->currentModel->GetVisNames(0);

    for (int i = 0; i < plotsData.size(); i++)
        PlotNames.push_back(std::string("plot") + QString::number(i).toStdString());

    for (int i = 0; i < SimNames.size(); i++)
        PlotNames.push_back(SimNames[i]);

    for (int i = 0; i < flagStrings::VisNames.size(); i++)
        PlotNames.push_back(flagStrings::VisNames[i]);

    for (int i = 0; i < (*currentProject)->currentModel->GetNumberParticlesFlows(); i++)
    {
        if ((*currentProject)->currentModel->GetNumberParticlesFlowsTypes()[i] == 5)
            continue;

        for (int j = 0; j < flagStrings::VisFlowNames.size(); j++)
            PlotNames.push_back(flagStrings::VisFlowNames[j] + QString::number(i).toStdString());
    }

    for (int i = 0; i < (*currentProject)->currentModel->GetConductorsList().size(); i++)
    {
        for (int j = 0; j < flagStrings::VisElectrodeNames.size(); j++)
            PlotNames.push_back(flagStrings::VisElectrodeNames[j] + QString::number(i).toStdString());
    }

    for (int i = 0; i < NowPlotNames.size(); i++)
    {
        int s = plotsData.size();
        for (int j = 0; j < s; j++)
        {
            if (NowPlotNames[i] == PlotNames[j])
            {
                PlotNames.erase(PlotNames.begin() + j);
                j--;
                s--;
            };
        };
    };

    ListsSelections.push_back(new ListSelectionMenu("Boundaries"));
    ListsSelections.back()->Create(PlotNames, NowPlotNames, std::bind(&MiddleWidget::ChangeVisualizationPlots, this),
                                   NULL);
    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 0, 0);
};
void MiddleWidget::ChangeVisualizationPlots()
{
    std::vector<std::vector<std::string>> list1;
    std::vector<std::vector<std::string>> list2;
    GetListsData(list1, list2);
    (*currentProject)->currentModel->SetPlotNames(list2[0]);
};

void MiddleWidget::ShowAddLinePlotMenu()
{

    std::vector<std::string> Names;
    if ((*currentProject)->problemType == 1)
        Names = flagStrings::PlotFlags2d;

    if ((*currentProject)->problemType == 2)
        Names = flagStrings::PlotFlags2daxs;

    if ((*currentProject)->problemType == 3)
        Names = flagStrings::PlotFlags2dpolar;

    if ((*currentProject)->problemType == 4)
        Names = flagStrings::PlotFlags3d;

    ListsSelections.push_back(new ListSelectionMenu("Grid value"));
    ListsSelections.back()->Create(Names, std::vector<std::string>(), NULL, NULL);
    middleWidgetGrid->addWidget(ListsSelections.back()->GetPointer(), 0, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Point 1"));
    groupBoxes.back()->Create({"x1", "x2", "x3"}, std::vector<double>{0, 0, 0});
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 1, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Point 2"));
    groupBoxes.back()->Create({"x1", "x2", "x3"}, std::vector<double>{0, 0, 0});
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 2, 0);

    groupBoxes.push_back(new GroupBoxWithItems("Grid value type"));
    groupBoxes.back()->Create({}, {}, flagStrings::PlotType, 0);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), 3, 0);

    mainMiddleWidgetButton = new QPushButton(tr("Add"));
    mainMiddleWidgetButton->setMaximumSize(QSize(70, 30));
    middleWidgetGrid->addWidget(mainMiddleWidgetButton, 4, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(AddPlot()));
};
void MiddleWidget::AddPlot()
{
    std::vector<std::vector<double>>      parameters1;
    std::vector<std::vector<std::string>> parameters2;
    std::vector<std::vector<std::string>> parameters3;
    std::string                           errorMessage;

    bool ok = FetchParameters(parameters1, parameters2, parameters3);

    if (!ok)
    {
        QMessageBox::warning(this, "Daisi warning", "Incorrect non-digit input");
        return;
    }

    std::vector<std::vector<std::string>> list1;
    std::vector<std::vector<std::string>> list2;
    GetListsData(list1, list2);

    int PlotTypeFlag = parameters1[2][0];

    std::vector<std::string> flag = list2[0];
    if ((*currentProject)->problemType == 1 || (*currentProject)->problemType == 2 ||
        (*currentProject)->problemType == 3)
        (*currentProject)
            ->currentModel->addPlot(flag, parameters1[0][0], parameters1[0][1], parameters1[1][0], parameters1[1][1],
                                    PlotTypeFlag);
    if ((*currentProject)->problemType == 4)
        (*currentProject)
            ->currentModel->addPlot(flag, parameters1[0][0], parameters1[0][1], parameters1[0][2], parameters1[1][0],
                                    parameters1[1][1], parameters1[1][2], PlotTypeFlag);

    std::vector<std::shared_ptr<lineplot>> plotsData = (*currentProject)->currentModel->GetPlotsVector();

    std::vector<MyTreeItem*> plotsVector;

    MyTreeItem* newItem = new MyTreeItem();
    int         i       = plotsData.size();
    newItem->setText(0, QString("plot") + QString::number(i - 1));
    newItem->flag  = flagStrings::lineplotsList;
    newItem->flag1 = i - 1;

    currentClickItem->addChild(newItem);
    currentClickItem->setExpanded(true);
    newItem->setSelected(true);
};
void MiddleWidget::ShowLinePlotMenu(int currentplotNumberIn)
{

    QListWidget* List = new QListWidget();

    QListWidgetItem* newitem1 = new QListWidgetItem();
    newitem1->setText("Line");
    List->addItem(newitem1);

    QListWidgetItem* newitem2 = new QListWidgetItem();
    newitem2->setText("Plot along line");
    List->addItem(newitem2);

    middleWidgetGrid->addWidget(List, 0, 0);

    connect(List, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listShowPlotClick(QListWidgetItem*)));

    currentplotNumber = currentplotNumberIn;

    std::vector<std::shared_ptr<lineplot>> plotsData = (*currentProject)->currentModel->GetPlotsVector();

    (*VTKArray)[0]->setDataFunction(ShowLinePlot, (*currentProject)->currentModel, currentplotNumber,
                                    plotsData[currentplotNumber]);
    (*VTKArray)[0]->refresh(0);
};
void MiddleWidget::ShowPlot2dMenu()
{

    std::vector<std::string> Names;
    if ((*currentProject)->problemType == 1)
        Names = flagStrings::PlotFlags2d;

    if ((*currentProject)->problemType == 2)
        Names = flagStrings::PlotFlags2daxs;

    if ((*currentProject)->problemType == 3)
        Names = flagStrings::PlotFlags2dpolar;

    if ((*currentProject)->problemType == 4)
        Names = flagStrings::PlotFlags3d;

    QListWidget* List = new QListWidget();

    for (int i = 0; i < Names.size(); i++)
    {
        QListWidgetItem* newitem2 = new QListWidgetItem();
        newitem2->setText(Names[i].c_str());
        List->addItem(newitem2);
    }

    middleWidgetGrid->addWidget(List, 0, 0);
    int pos = 1;

    /*if ((*currentProject)->problemType == 4)
    {
            QGridLayout *vbox = new QGridLayout;
            QGroupBox *groupBox = new QGroupBox(tr("2d cut plane"));

            groupBox->setMaximumHeight(150);
            middleWidgetGrid->addWidget((QWidget*)groupBox, 2, 0);
            ListPlanes[0] = new QRadioButton(tr(flagStrings::PlotFlags3dPlane[0].c_str()));
            ListPlanes[1] = new QRadioButton(tr(flagStrings::PlotFlags3dPlane[1].c_str()));
            ListPlanes[2] = new QRadioButton(tr(flagStrings::PlotFlags3dPlane[2].c_str()));

            vbox->addWidget(ListPlanes[0], 0, 0);
            vbox->addWidget(ListPlanes[1], 1, 0);
            vbox->addWidget(ListPlanes[2], 2, 0);

            groupBox->setLayout(vbox);

            ListPlanes[0]->setChecked(true);

            middleWidgetGrid->addWidget(groupBox, 1, 0);


            MiddleWidgetTextEdit1 = new MyQTextEdit();
            MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));
            MiddleWidgetTextEdit1->setText(QString::number(0));
            middleWidgetGrid->addWidget(MiddleWidgetTextEdit1, 2, 0);

            pos = 3;
    }*/

    groupBoxes.push_back(new GroupBoxWithItems("Grid value type"));
    groupBoxes.back()->Create({}, {}, flagStrings::PlotType, 0);
    middleWidgetGrid->addWidget(groupBoxes.back()->GetPointer(), pos, 0);

    connect(List, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listShowPlot2dClick(QListWidgetItem*)));
};
void MiddleWidget::listShowPlot2dClick(QListWidgetItem* item)
{
    std::vector<double> paramf;
    groupBoxes.back()->GetParameters(paramf);

    int PlotTypeFlag = paramf[0];
    if ((*currentProject)->problemType == 4)
    {
        /*std::string name = item->text().toStdString();

        int planeType;
        for (int i = 0; i<3; i++)
        {
        if (ListPlanes[i]->isChecked())
        {
        planeType = i;
        break;
        }
        };
        double param = MiddleWidgetTextEdit1->toDoubleMy(&ok);
        (*VTKArray)[0]->setDataFunction(Show2dPlot3d, (*currentProject)->currentModel, name, planeType, param,
        (*currentProject)->precisionType, PlotTypeFlag);
        (*VTKArray)[0]->refresh(0);*/
    }
    else
    {
        std::string name = item->text().toStdString();
        (*VTKArray)[0]->setDataFunction(Show2dPlot, (*currentProject)->currentModel, name,
                                        (*currentProject)->precisionType, PlotTypeFlag);
        (*VTKArray)[0]->refresh(0);
    }
};
void MiddleWidget::listShowPlotClick(QListWidgetItem* item)
{
    //	Visualization_thread = std::thread(&Daizy::VisualThr, this, item);

    std::string name = item->text().toStdString();

    std::vector<std::shared_ptr<lineplot>> plotsData = (*currentProject)->currentModel->GetPlotsVector();

    if (name == std::string("Line"))
    {
        (*VTKArray)[0]->setDataFunction(ShowPlotLine, (*currentProject)->currentModel, plotsData[currentplotNumber]);
        (*VTKArray)[0]->refresh(0);
        return;
    };
    if (name == std::string("Plot along line"))
    {

        (*VTKArray)[0]->setDataFunction(ShowLinePlot, (*currentProject)->currentModel, currentplotNumber,
                                        plotsData[currentplotNumber]);
        (*VTKArray)[0]->refresh(0);
        return;
    };
};
