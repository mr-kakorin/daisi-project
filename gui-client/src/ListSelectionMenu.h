#ifndef ListSelectionMenu_H
#define ListSelectionMenu_H
#include <QtWidgets/QWidget>
#include <functional>
class QPushButton;
class QListWidget;
class QGridLayout;
class QGroupBox;
class QPushButton;
class QListWidgetItem;
class QVBoxLayout;
class ListSelectionMenu : public QWidget
{
    Q_OBJECT
    QGridLayout*                  grid;
    QListWidget*                  leftList;
    QListWidget*                  rigthList;
    QGroupBox*                    out;
    QVBoxLayout*                  vbox;
    QPushButton*                  left;
    QPushButton*                  rigth;
    std::vector<QListWidgetItem*> items;

    std::function<void()>    listExec;
    std::function<void(int)> ClickExec;

  public:
    void GetListsData(std::vector<int>& list1, std::vector<int>& list2);
    void GetListsData(std::vector<std::string>& list1, std::vector<std::string>& list2);
    QGroupBox* GetPointer();
    ListSelectionMenu();
    ListSelectionMenu(const QString& name);
    ~ListSelectionMenu();
    void Create(const std::vector<int>&, const std::vector<int>&, const std::function<void()>& listExec,
                const std::function<void(int)>& ClickExec);
    void Create(const std::vector<std::string>&, const std::vector<std::string>&, const std::function<void()>& listExec,
                const std::function<void(int)>& ClickExec);
  public slots:
    void listClick(QListWidgetItem* item);
    void leftList2rigth();
    void rigthList2Left();
};
#endif
