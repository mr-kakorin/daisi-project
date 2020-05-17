#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "FlagStringsD.h"
#include "GeneralTools.h"
#include "MiddleWidget.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"
#include <memory>

#include <QVTKWidget.h>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSplitter>
#include <QtWidgets/QTreeWidgetItem>
#include <QtWidgets/QVBoxLayout>

void Daizy::PrepareRightWidget()
{
    for (int i = 1; i < numberOfResultPlots; i++)
    {
        graphicsArray.push_back(new QVTKWidget());
        VTKArray.push_back(new vtkComponent());
        VTKArray.back()->setWidget(graphicsArray.back());
    }

    int s1 = round(sqrt(numberOfResultPlots) + 0.49);

    int ind;
    int flagBreak = 0;
    for (int i = 0; i < s1; i++)
    {
        for (int j = 0; j < s1; j++)
        {
            ind = i * s1 + j;
            if (ind > 0 && ind < numberOfResultPlots)
            {
                rigthWidgetGrid->addWidget(graphicsArray[ind], j, i);
            }
        }
    }
};
void Daizy::ShowResultsMenu()
{

    PrepareRightWidget();

    if (currentProject->problemType < 7)
        middleWidgetAggregator->showSimulationsResuts(numberOfResultPlots);

    /*QGridLayout *grid = new QGridLayout;
    QGroupBox *groupBox1 = new QGroupBox(tr("Font and Tics"));
    groupBox1->setMaximumHeight(3000);
    groupBox1->setLayout(grid);

    middleWidgetGrid->addWidget(groupBox1, numberOfResultPlots + 1, 0);

    std::vector<std::string> names = { "Font Size", "X axis Tics", "Y axis Tics" };
    std::vector<int> values = { 25, -1, -1 };

    MiddleWidgetTextEditVector.clear();
    int i = 0;
    for (i = 0; i < names.size(); i++)
    {
            QLabel* l = new QLabel();
            l->setText(names[i].c_str());
            grid->addWidget(l, i, 0);

            MiddleWidgetTextEditVector.push_back(new MyQTextEdit());
            MiddleWidgetTextEditVector.back()->setMaximumSize(QSize(16777215, 25));
            MiddleWidgetTextEditVector.back()->setText(QString::number(values[i]));
            grid->addWidget(MiddleWidgetTextEditVector.back(), i, 1);
    }*/

    MiddleWidgetButton1 = new QPushButton("Add Plot");
    MiddleWidgetButton2 = new QPushButton("Remove Plot");

    middleWidgetAggregator->GetMiddleWidgetGrid()->addWidget(MiddleWidgetButton1, numberOfResultPlots + 1, 0);
    middleWidgetAggregator->GetMiddleWidgetGrid()->addWidget(MiddleWidgetButton2, numberOfResultPlots + 2, 0);

    connect(MiddleWidgetButton1, SIGNAL(clicked()), this, SLOT(AddResultPlot()));
    connect(MiddleWidgetButton2, SIGNAL(clicked()), this, SLOT(RemoveResultPlot()));
};

void Daizy::AddResultPlot()
{
    numberOfResultPlots++;
    graphicsArray.push_back(new QVTKWidget());
    VTKArray.push_back(new vtkComponent());
    VTKArray.back()->setWidget(graphicsArray.back());

    int s1 = round(sqrt(numberOfResultPlots) + 0.49);

    int ind, i, j;
    int flagBreak = 0;
    for (i = 0; i < s1; i++)
    {
        for (j = 0; j < s1; j++)
        {
            ind = i * s1 + j;
            if (ind == numberOfResultPlots - 1)
            {
                flagBreak = 1;
                break;
            }
        }
        if (flagBreak)
            break;
    }

    rigthWidgetGrid->addWidget(graphicsArray.back(), j, i);
    ShowResultsMenu();
};
void Daizy::RemoveResultPlot()
{
    if (numberOfResultPlots == 1)
        return;

    // for(int i=0;i<numberOfResultPlots;i++)
    //	rigthWidgetGrid->removeWidget(graphicsArray[i]);

    clearLayout(rigthWidgetGrid);

    VTKArray.pop_back();

    numberOfResultPlots--;
    graphicsArray.clear();

    for (int i = 0; i < VTKArray.size(); i++)
        delete VTKArray[i];

    VTKArray.clear();
    for (int i = 0; i < numberOfResultPlots; i++)
    {
        VTKArray.push_back(new vtkComponent());
        graphicsArray.push_back(new QVTKWidget());
        VTKArray[i]->setWidget(graphicsArray[i]);
    }

    int s1 = round(sqrt(numberOfResultPlots) + 0.49);

    int ind;
    for (int i = 0; i < s1; i++)
    {
        for (int j = 0; j < s1; j++)
        {
            ind = i * s1 + j;
            if (ind < numberOfResultPlots)
            {
                rigthWidgetGrid->addWidget(graphicsArray[ind], j, i);
            }
        }
    }

    //	for (int i = 0; i < numberOfResultPlots; i++)
    //	VTKArray[i]->refresh(0);

    ShowResultsMenu();
};