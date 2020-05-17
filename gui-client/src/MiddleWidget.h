#ifndef MiddleWidget_H
#define MiddleWidget_H
#include <QtWidgets/QWidget>

/////////forward declarations////////
namespace Dproject
{
class project;
};
class vtkComponent;
class GroupBoxWithItems;
class FilesBrouses;
class GroupBoxWithTextItems;
class ListSelectionMenu;
class MyTreeItem;
class QGroupBox;
class QTreeWidget;
class QListWidgetItem;
class QPushButton;
class QGridLayout;
class QTreeWidgetItem;

class MiddleWidget : public QWidget
{
    Q_OBJECT
    std::vector<GroupBoxWithItems*>       groupBoxes;
    std::vector<FilesBrouses*>            filesBrouses;
    std::vector<QTreeWidget*>             TreeList;
    std::vector<ListSelectionMenu*>       ListsSelections;
    std::vector<std::vector<MyTreeItem*>> itemvector;
    std::vector<GroupBoxWithTextItems*>   groupBoxWithTextItems;
    QGridLayout*                          middleWidgetGrid;
    QPushButton*                          mainMiddleWidgetButton;
    QPushButton*                          mainMiddleWidgetButton1;
    Dproject::project**                   currentProject;
    std::vector<vtkComponent*>*           VTKArray;
    QGroupBox*                            groupBoxFictive;
    int                                   currentFlow;
    int                                   currentplotNumber;

  public:
    int currentConditioin;
    /////////////////General///////////////
    QGridLayout* GetMiddleWidgetGrid();
    MyTreeItem*  currentClickItem;
    MiddleWidget();
    MiddleWidget(QWidget* widget, Dproject::project** currentProjectIn, std::vector<vtkComponent*>* VTKArrayIn);
    void InitTrees(int n);
    void AddLastFictiveBox(int pos);
    bool FetchParameters(std::vector<std::vector<double>>&      parameters1,
                         std::vector<std::vector<std::string>>& parameters2,
                         std::vector<std::vector<std::string>>& parameters3);
    void GetListsData(std::vector<std::vector<int>>& list1, std::vector<std::vector<int>>& list2);
    void GetListsData(std::vector<std::vector<std::string>>& list1, std::vector<std::vector<std::string>>& list2);
    /////////////////Accelerator///////////////
    void ShowDefaultCondMenu(int currentConditioin, QString name, std::string flag);

    ///////////////Simulations///////////////
    void showSimulationsResuts(int numberOfResultPlots);
    void ShowAddBoundary();
    void AddCondition();
    void PotentialFieldMenu(int currentConditioinIn);
    void ShowSelectionMenu(int conditionNumber, std::string flag1, int flag2);
    void HighLightBoundary(int n);
    void applyBoundary();
    void ShowGlobalFieldMenu();
    void ApplyDefaultBoundaries();
    void applyglobalFieldConditions();
    void AddConductor();
    void ShowConductorSelectionMenu(int number);
    void ApplyConductorBoundaries(int n);
    void ShowMeshMenu();
    void SetMeshBoundary();

    /// Flow///
    void ShowAddFlowMenu();
    void AddFlowCondition();
    void ShowFlowEmitterProperty(int);
    void ShowFlowEmitterPropertyPlasma(int);
    void ShowFlowEmitterPropertyEvaporation(int);
    void ShowFlowEmitterPropertyField(int);
    void ShowFlowEmitterPropertyThermionic(int);
    void ShowFlowEmitterPropertyPhotoemission(int);
    void ShowFlowEmitterPropertyAccelerator(int);
    void SetEmissionParametersLinac();
    void ShowFlowSummary(int i);
    void ShowFlowState(int);
    void ChangeEmissionCurrentStyle();

    /////Visualizations///
    void ShowAddPlotMenu();
    void ShowAddLinePlotMenu();
    void ShowLinePlotMenu(int currentplotNumber);
    void ShowPlot2dMenu();

    ////Solver///
    void ShowEmissionModelSolverSettings();
    void ApplySolverEmissionModelSettings();
    void ShowSolverMenuPIC();
    void ShowSolverMenuPTI();
    void ShowFieldSolverMenu();
    void ShowSolverSettings();
    void ApplySolverSettingsPTI();
    void ApplySolverSettingsPIC();

  public slots:
    ///////////////General///////////////////
    void ShowCreateNewProject();
    void clear();
    void showSummary();

    //////////////Simulations///////////////
    void listShowResultClick(QTreeWidgetItem* item, int);
    int AddBoundariesEvent(const std::string& input, std::string& errorMsg);
    int AddBoundaryEvent(const std::string& input, std::string& errorMsg);
    void listShowConductorStateClick(QListWidgetItem* item);
    void AddEmittance();
    void listShowResultClickEmittances(QTreeWidgetItem* item, int);
    void listShowFlowStateClick(QListWidgetItem*);
    void ChangeFlowDistrStyle();
    void ChangeGeom();
    void ApplyEmitterBoundary();
    void AddPlot();
    void listShowPlotClick(QListWidgetItem* item);
    void listShowPlot2dClick(QListWidgetItem* item);
    void ChangeVisualizationPlots();

  signals:
    void DefaultBuilder();
    void CreateNewProject();
    void ShowProjectTreeSim();
    void AddPotentialEvent();
    void AddConductorEvent();
    void AddFlowConditionEvent();
    void GenerateMesh();
    void AddFlow();
    void FieldSimulate();
    void SimulateIterative();
    void SimulatePIC();
    void SimulatePICContinue();
    void ResetFlags();
};
#endif
