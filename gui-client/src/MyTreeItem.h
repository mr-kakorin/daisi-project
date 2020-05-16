#ifndef MYTREEITEM_H
#define MYTREEITEM_H
#include <QTextEdit>
#include <QtWidgets>
class MyQListWidgetItem : public QListWidgetItem
{
  public:
    int flag;
    int flag1;
};
class MyTreeItem : public QTreeWidgetItem
{
  public:
    std::string flag;
    std::string flag3;
    std::string flagSearchBoundary;
    int         flag1;
    int         flag2;
    int         flag4;
    int         flag5;
    int         flag6;
    int         flagSearchBoundaryI;
    MyTreeItem()
    {
        flagSearchBoundaryI = 0;
    };
};
class MyQTextEdit : public QTextEdit
{
  public:
    double toDoubleMy(bool* ok);
};
#endif // MYTREEITEM_H
