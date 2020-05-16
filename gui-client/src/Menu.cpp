#include <daisi-solver/Results.h>
#include <daisi-solver/project.h>

#include "DaisiClient.h"
#include "VisualizationFunctions.h"
#include "vtkComponent.h"
#include <QtWidgets>

#include "DaisiClient.h"
#include "GeneralTools.h"
#include "Menu.h"
#include "MyTreeItem.h"
#include "vtkComponent.h"
#include <QtWidgets/QAction>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>

#include "names.h"

void MainMenu::MenuShow(QMenuBar* menuBar, Daizy* mainWindowIn, Dproject::project** currentProjectIn,
                        std::vector<vtkComponent*>* VTKArrayIn, QGridLayout* middleWidgetGridIn)
{
    middleWidgetGrid = middleWidgetGridIn;
    currentProject   = currentProjectIn;
    VTKArray         = VTKArrayIn;
    mainWindow       = mainWindowIn;
    menuFile         = new QMenu(menuBar);
    menuFile->setTitle("Project");
    menuBar->addAction(menuFile->menuAction());

    menuModel = new QMenu(menuBar);
    menuModel->setTitle("Model");
    menuBar->addAction(menuModel->menuAction());

    actionSave_Model = new QAction(mainWindow);
    menuModel->addAction(actionSave_Model);
    actionSave_Model->setText("Save Model");
    connect(actionSave_Model, SIGNAL(triggered()), this, SLOT(SaveModelEvent()));

    actionLoad_Model = new QAction(mainWindow);
    menuModel->addAction(actionLoad_Model);
    connect(actionLoad_Model, SIGNAL(triggered()), this, SLOT(LoadModelEvent()));
    actionLoad_Model->setText("Load Model");

    actionNew_Project = new QAction(mainWindow);
    menuFile->addAction(actionNew_Project);
    connect(actionNew_Project, SIGNAL(triggered()), this, SIGNAL(ShowCreateNewProject()));
    actionNew_Project->setText("New Project");

    actionSave_Project = new QAction(mainWindow);
    menuFile->addAction(actionSave_Project);
    connect(actionSave_Project, SIGNAL(triggered()), this, SLOT(SaveProjectEvent()));
    actionSave_Project->setText("Save Project");

    actionLoad_Project = new QAction(mainWindow);
    menuFile->addAction(actionLoad_Project);
    connect(actionLoad_Project, SIGNAL(triggered()), this, SLOT(LoadProjectEvent()));
    actionLoad_Project->setText("Load Project");

    menuData = new QMenu(menuBar);
    menuData->setTitle("Simulations data");
    menuBar->addAction(menuData->menuAction());

    actionSave_Data = new QAction(mainWindow);
    menuData->addAction(actionSave_Data);
    actionSave_Data->setText("Save Simulations data");
    connect(actionSave_Data, SIGNAL(triggered()), this, SLOT(SaveDataEvent()));

    actionLoad_Data = new QAction(mainWindow);
    menuData->addAction(actionLoad_Data);
    connect(actionLoad_Data, SIGNAL(triggered()), this, SLOT(LoadDataEvent()));
    actionLoad_Data->setText("Load Simulations data");

    menuExport = new QMenu(menuBar);
    menuExport->setTitle("Export data");
    menuBar->addAction(menuExport->menuAction());

    actionExport_Data = new QAction(mainWindow);
    menuExport->addAction(actionExport_Data);
    connect(actionExport_Data, SIGNAL(triggered()), this, SLOT(ExportDataEvent()));
    actionExport_Data->setText("Export plot data to text");

    actionExport_DataEps = new QAction(mainWindow);
    menuExport->addAction(actionExport_DataEps);
    connect(actionExport_DataEps, SIGNAL(triggered()), this, SLOT(ExportDataEventEps()));
    actionExport_DataEps->setText("Export plot data to eps");

    about = new QMenu(menuBar);
    about->setTitle("Help");
    menuBar->addAction(about->menuAction());

    actionAbout = new QAction(mainWindow);
    about->addAction(actionAbout);
    connect(actionAbout, SIGNAL(triggered()), this, SLOT(AboutEvent()));
    actionAbout->setText("About DAISI");

    actionDevelopers = new QAction(mainWindow);
    about->addAction(actionDevelopers);
    connect(actionDevelopers, SIGNAL(triggered()), this, SLOT(DevelopersEvent()));
    actionDevelopers->setText("Developers");

    actionLicense = new QAction(mainWindow);
    about->addAction(actionLicense);
    connect(actionLicense, SIGNAL(triggered()), this, SLOT(LicenseEvent()));
    actionLicense->setText("License");
};
void MainMenu::AboutEvent()
{
    QMessageBox::about(mainWindow, "About", aboutStr0.c_str());
};
void MainMenu::LicenseEvent()
{

    QMessageBox::about(mainWindow, "License", licStr.c_str());
};
void MainMenu::DevelopersEvent()
{
    QMessageBox::about(mainWindow, "Developers", developers.c_str());
};

void MainMenu::ExportDataEvent()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };
    emit clear();
    emit ResetFlags();

    QGroupBox*   groupBox = new QGroupBox(tr("Save Data"));
    QVBoxLayout* vbox     = new QVBoxLayout;

    QLabel* label = new QLabel(groupBox);
    label->setObjectName(QStringLiteral("File name"));
    label->setText("File name");

    vbox->addWidget(label);

    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setObjectName(QStringLiteral("textEdit"));
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    vbox->addWidget(MiddleWidgetTextEdit1);

    mainMiddleWidgetButton = new QPushButton(tr("Save Data"));
    mainMiddleWidgetButton->setMaximumSize(QSize(120, 30));
    ;
    vbox->addWidget(mainMiddleWidgetButton);

    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    middleWidgetGrid->addWidget((QWidget*)groupBox, 2, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(ExportDataFileEvent()));
    // setLayout(middleWidgetGrid);
};

void MainMenu::ExportDataEventEps()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };
    emit clear();
    emit ResetFlags();

    QGroupBox*   groupBox = new QGroupBox(tr("Save Data"));
    QVBoxLayout* vbox     = new QVBoxLayout;

    QLabel* label = new QLabel(groupBox);
    label->setObjectName(QStringLiteral("File name"));
    label->setText("File name");

    vbox->addWidget(label);

    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setObjectName(QStringLiteral("textEdit"));
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    vbox->addWidget(MiddleWidgetTextEdit1);

    mainMiddleWidgetButton = new QPushButton(tr("Save Data To Eps"));
    mainMiddleWidgetButton->setMaximumSize(QSize(120, 30));
    vbox->addWidget(mainMiddleWidgetButton);

    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    middleWidgetGrid->addWidget((QWidget*)groupBox, 2, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(ExportDataFileEventEps()));
    // setLayout(middleWidgetGrid);
};

void MainMenu::ExportDataFileEvent()
{
    std::string name = (*currentProject)->projectFolder + std::string("/dataFiles/") +
                       MiddleWidgetTextEdit1->toPlainText().toStdString() + std::string(".txt");
    (*VTKArray)[0]->SaveData2File(name);
};

void MainMenu::ExportDataFileEventEps()
{
    std::string name = (*currentProject)->projectFolder + std::string("/dataFiles/") +
                       MiddleWidgetTextEdit1->toPlainText().toStdString();
    (*VTKArray)[0]->SaveData2EpsFile(name);
};

void MainMenu::SaveDataEvent()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };

    emit clear();
    emit ResetFlags();

    QGroupBox*   groupBox = new QGroupBox(tr("Save data"));
    QVBoxLayout* vbox     = new QVBoxLayout;

    QLabel* label = new QLabel(groupBox);
    label->setObjectName(QStringLiteral("File name"));
    label->setText("File name");

    vbox->addWidget(label);

    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setObjectName(QStringLiteral("textEdit"));
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    vbox->addWidget(MiddleWidgetTextEdit1);

    mainMiddleWidgetButton = new QPushButton(tr("Save data"));
    mainMiddleWidgetButton->setMaximumSize(QSize(120, 30));
    vbox->addWidget(mainMiddleWidgetButton);

    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    middleWidgetGrid->addWidget((QWidget*)groupBox, 2, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(SaveDataFileEvent()));
};
void MainMenu::SaveProjectEvent()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };
    emit saveParameters();
    (*currentProject)->SaveProject(version);
};
void MainMenu::LoadProjectEvent()
{
    emit    ResetFlags();
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open File"), "..//", tr("Daizy projects (*.dproj)"));
    // delete (*currentProject);

    if (fileName.size() == 0)
        return;

    if (!(*currentProject))
        (*currentProject) = new Dproject::project();

    std::string errorMsgOpen;
    (*currentProject)->LoadProject(fileName.toStdString(), version, errorMsgOpen);

    if (errorMsgOpen.size())
    {
        QMessageBox::critical(this, "Daisi error", errorMsgOpen.c_str());
        return;
    };

    if ((*currentProject)->problemType >= 7)
        emit ShowProjectTreeAccel();
    else
        emit ShowProjectTreeSim();

    emit clear();
    emit showSummary();
};
void MainMenu::LoadDataEvent()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };

    QString fileName = QFileDialog::getOpenFileName(
        this, tr("Open File"), ((*currentProject)->projectFolder + std::string("/dataFiles")).c_str(),
        tr("Simulations data (*.dyn)"));

    if (fileName.size() == 0)
        return;

    (*currentProject)->LoadData(fileName.toStdString());
};

void MainMenu::LoadModelEvent()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };

    emit ResetFlags();

    QString fileName = QFileDialog::getOpenFileName(
        this, tr("Open File"), ((*currentProject)->projectFolder + std::string("/modelFiles")).c_str(),
        tr("Daisi models (*.mdl)"));

    if (!fileName.size())
        return;

    std::string errorMsgOpen;

    (*currentProject)->ModelLoad(fileName.toStdString(), version, errorMsgOpen);

    if (errorMsgOpen.size())
    {
        QMessageBox::critical(this, "Daisi error", errorMsgOpen.c_str());
        return;
    };

    if ((*currentProject)->problemType >= 7)
        emit ShowProjectTreeAccel();
    else
        emit ShowProjectTreeSim();
};
void MainMenu::SaveModelEvent()
{
    if (!(*currentProject))
    {
        QMessageBox::warning(mainWindow, "Daisi warning", "No project load");
        return;
    };
    emit saveParameters();

    emit         ResetFlags();
    emit         clear();
    QGroupBox*   groupBox = new QGroupBox(tr("Save model"));
    QVBoxLayout* vbox     = new QVBoxLayout;

    QLabel* label = new QLabel(groupBox);
    label->setObjectName(QStringLiteral("File name"));
    label->setText(QApplication::translate("MainWindow", "File name", 0));

    vbox->addWidget(label);

    MiddleWidgetTextEdit1 = new MyQTextEdit();
    MiddleWidgetTextEdit1->setObjectName(QStringLiteral("textEdit"));
    MiddleWidgetTextEdit1->setMaximumSize(QSize(16777215, 25));

    vbox->addWidget(MiddleWidgetTextEdit1);

    mainMiddleWidgetButton = new QPushButton(tr("Save model"));
    mainMiddleWidgetButton->setMaximumSize(QSize(120, 30));
    vbox->addWidget(mainMiddleWidgetButton);

    vbox->addStretch(1);
    groupBox->setLayout(vbox);

    middleWidgetGrid->addWidget((QWidget*)groupBox, 2, 0);

    connect(mainMiddleWidgetButton, SIGNAL(clicked()), this, SLOT(SaveModelFileEvent()));
};
void MainMenu::SaveModelFileEvent()
{
    (*currentProject)
        ->SaveCurrentModel((MiddleWidgetTextEdit1->toPlainText().toStdString() + std::string(".mdl")), version);
};
void MainMenu::SaveDataFileEvent()
{
    (*currentProject)->SaveData((MiddleWidgetTextEdit1->toPlainText().toStdString() + std::string(".dyn")));
};