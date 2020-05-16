#ifndef DAIZY_H
#define DAIZY_H
#include <QMainWindow>
#include <functional>
#include <iostream>
#include <thread>

/////////forward declarations////////
class QDialogButtonBox;
class QLabel;
class QProgressBar;
class QPushButton;
class QAction;
class QGroupBox;
class vtkComponent;
class MainMenu;
class QSplitter;
class MyTreeItem;
class QVBoxLayout;
class QTreeWidget;
class QGridLayout;
class MiddleWidget;
class QHBoxLayout;
class QListWidgetItem;
class QTreeWidgetItem;
class QVTKWidget;
class QTextEdit;
class QScrollBar;

namespace Dproject
{
class project;
};

/////////Main class////////
class Daizy : public QMainWindow
{
    Q_OBJECT
  private:
    /////////////Project/////////////

    Dproject::project* currentProject;

    /////////////some flags/////////////
    std::string              prevFlag;
    int                      prevFlagInt;
    std::string              prevFlagString;
    std::string              prevflagSearchBoundary;
    int                      prevflagSearchBoundaryI;
    bool                     flagAbort11;
    std::string              lastFlag;
    std::string              version;
    bool                     flagAbort;
    int                      numberOfResultPlots;
    int                      currentAccelSection;
    int                      currentplotNumber;
    std::string              currentsolver;
    MyTreeItem*              currentClickItem;
    bool                     flagExit;
    std::string              errorMsg;
    double                   progress;
    double                   progressLoc;
    std::vector<std::string> status;
    int                      problemType;
    int                      precisionType;
    int                      solverType;
    int                      currentConditioin;
    int                      currentFlow;

    /////////////Threads////////////////
    std::thread              work_thread;
    std::thread              Visualization_thread;
    std::thread              progressBarThread;
    std::vector<std::thread> viz_threads;

    ////////////Time and progress/////////
    QTimer*       timerViz;
    QLabel*       labelTime;
    QProgressBar* progressBar;
    QProgressBar* progressBarLocal;
    QTextEdit*    monitor;
    QScrollBar*   scrollBar;

    //////////Visualization////////////////

    std::vector<vtkComponent*> VTKArray;
    std::vector<QVTKWidget*>   graphicsArray;

    ///////////Some GUI obejects//////////
    QSplitter* splitter;
    QSplitter* splitterStatus;

    QWidget*     centralWidget;
    QPushButton* buttonAbort;
    QVBoxLayout* rightLayout;
    QMenuBar*    menuBar;
    QPushButton* MiddleWidgetButton1;
    QPushButton* MiddleWidgetButton2;
    QWidget*     leftWidget;
    QTreeWidget* leftWidgetTree;
    QGridLayout* leftWidgetGrid;
    QWidget*     rigthWidget;
    QHBoxLayout* horizontalLayout;
    QGridLayout* SolutionParamsGrid;
    QGroupBox*   SolutionParamsBox;
    QGridLayout* middleWidgetGrid;
    QGridLayout* rigthWidgetGrid;
    QWidget*     middleWidget;
    ///////////Aggregation GUI obejects//////////
    MainMenu*     mainMenu;
    MiddleWidget* middleWidgetAggregator;

    ///////General///////////
    void MyConnect();
    void DefaultBuilder();
    void setProgress();
    void VisualThr(QListWidgetItem* item);
    void ShowResultsMenu();

    /////////Visualization/////////
    void RefreshGraphics1(std::vector<QVTKWidget*> graphicsArray, QGridLayout* rigthWidgetGrid);
    void PrepareRightWidget();
    void RefreshGraphics(std::vector<vtkComponent*> VTKArray);
    void PrepareVisualisation();

    //////project tree////////

  public:
    Daizy(QWidget* parent = 0);
    ~Daizy();
    void retranslateUi();
  public slots:

    //////////General///////////
    void stopTimer();
    void ResetFlags();
    void ShowMB();
    void CreateNewProjectEvent();
    void updateGrViz();
    void updateGr();
    void saveParameters();
    void AddResultPlot();
    void RemoveResultPlot();
    void AbortSimulations();
    void ShowProjectTreeAccel();
    void ShowProjectTreeSim();
    void currItemClickedAccelSave(MyTreeItem* item1);
    void currItemClickedSimSave(MyTreeItem* item1);

    ///////////////Accelerator/////////////
    void AccelSolve();
    void AddFlowAccel();
    void currItemClickedAccel(QTreeWidgetItem*, int);

    ///////////////Simulations/////////////
    void AddConductorEvent();
    void InitFieldSolver();
    void SimulateIterative();
    void SimulatePIC();
    void SimulatePICContinue();

    void currItemClicked(QTreeWidgetItem*, int);
    void AddPotentialEvent();
    void GenerateMesh();
    void AddFlow();
    void FieldSimulate();
    void AddFlowConditionEvent();

  signals:
    void showMBSignal();
};

#endif // DAIZY_H
