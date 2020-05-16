#ifndef Menu_H
#define Menu_H
#include <QtWidgets>
class QMenu;
class QAction;
class Daizy;
namespace Dproject
{
class project;
}
class vtkComponent;
class MyQTextEdit;
class QGridLayout;
class QPushButton;
class MiddleWidget;
class MainMenu : public QWidget
{
    Q_OBJECT

    Dproject::project**         currentProject;
    std::vector<vtkComponent*>* VTKArray;

    QMenu* menuFile;
    QMenu* menuModel;
    QMenu* menuData;
    QMenu* menuExport;
    QMenu* menuExportEps;
    QMenu* about;
    QMenu* contextMenu;
    QMenu* menuRClick;

    QAction* actionNew_Project;
    QAction* actionSave_Project;
    QAction* actionLoad_Project;
    QAction* actionSave_Model;
    QAction* actionLoad_Model;
    QAction* actionSave_Data;
    QAction* actionExport_Data;
    QAction* actionExport_DataEps;

    QAction* actionLoad_Data;
    QAction* actionAbout;
    QAction* actionDevelopers;
    QAction* actionLicense;
    QAction* actionLicense1;

    Daizy*       mainWindow;
    QAction*     actionMenuRClick;
    MyQTextEdit* MiddleWidgetTextEdit1;
    QGridLayout* middleWidgetGrid;
    QPushButton* mainMiddleWidgetButton;

  public:
    void MenuShow(QMenuBar* menuBar, Daizy* mainWindow, Dproject::project** currentProjectIn,
                  std::vector<vtkComponent*>* VTKArrayIn, QGridLayout* middleWidgetGridIn);
  public slots:
    void AboutEvent();
    void DevelopersEvent();
    void LicenseEvent();
    void ExportDataFileEvent();
    void ExportDataEvent();
    void ExportDataFileEventEps();
    void ExportDataEventEps();
    void SaveDataFileEvent();
    void SaveDataEvent();
    void LoadDataEvent();
    void SaveProjectEvent();
    void LoadProjectEvent();
    void LoadModelEvent();
    void SaveModelEvent();
    void SaveModelFileEvent();
  signals:
    void ShowProjectTreeAccel();
    void ShowProjectTreeSim();
    void saveParameters();
    void ResetFlags();
    void ShowCreateNewProject();
    void clear();
    void showSummary();
};
#endif