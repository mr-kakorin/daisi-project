#include "ListSelectionMenu.h"
#include <QFileDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

ListSelectionMenu::~ListSelectionMenu()
{
    for (int i = 0; i < leftList->count(); i++)
        leftList->item(i);

    for (int i = 0; i < rigthList->count(); i++)
        rigthList->item(i);

    items.clear();

    delete grid;
    delete leftList;
    delete rigthList;
    //	delete out;
    //	delete vbox;
    delete left;
    delete rigth;
};

ListSelectionMenu::ListSelectionMenu(){

};

ListSelectionMenu::ListSelectionMenu(const QString& name)
{
    grid = new QGridLayout();

    out      = new QGroupBox(name);
    leftList = new QListWidget();
    grid->addWidget(leftList, 0, 0);

    vbox = new QVBoxLayout;

    left = new QPushButton(tr("L-R"));
    left->setMaximumSize(QSize(50, 30));
    connect(left, SIGNAL(clicked()), this, SLOT(leftList2rigth()));

    rigth = new QPushButton(tr("R-L"));
    rigth->setMaximumSize(QSize(50, 30));
    connect(rigth, SIGNAL(clicked()), this, SLOT(rigthList2Left()));

    vbox->addWidget(left);
    vbox->addWidget(rigth);

    grid->addLayout(vbox, 0, 1);

    rigthList = new QListWidget();

    grid->addWidget(rigthList, 0, 2);

    out->setLayout(grid);

    connect(rigthList, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listClick(QListWidgetItem*)));

    connect(leftList, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(listClick(QListWidgetItem*)));
};

void ListSelectionMenu::Create(const std::vector<int>& L1, const std::vector<int>& L2,
                               const std::function<void()>& listExecIn, const std::function<void(int)>& ClickExecIn)
{
    listExec  = listExecIn;
    ClickExec = ClickExecIn;

    out->setMaximumHeight(std::max(L1.size(), L2.size()) * 30 + 60);

    for (int i = 0; i < L1.size(); i++)
    {
        items.push_back(new QListWidgetItem());
        items.back()->setText(QString::number(L1[i]));
        leftList->addItem(items.back());
    };
    for (int i = 0; i < L2.size(); i++)
    {
        items.push_back(new QListWidgetItem());
        items.back()->setText(QString::number(L2[i]));
        rigthList->addItem(items.back());
    };
};

void ListSelectionMenu::Create(const std::vector<std::string>& L1, const std::vector<std::string>& L2,
                               const std::function<void()>& listExecIn, const std::function<void(int)>& ClickExecIn)
{

    listExec  = listExecIn;
    ClickExec = ClickExecIn;

    for (int i = 0; i < L1.size(); i++)
    {
        items.push_back(new QListWidgetItem());
        items.back()->setText(L1[i].c_str());
        leftList->addItem(items.back());
    };

    for (int i = 0; i < L2.size(); i++)
    {
        items.push_back(new QListWidgetItem());
        items.back()->setText(L2[i].c_str());
        rigthList->addItem(items.back());
    };
};

void ListSelectionMenu::listClick(QListWidgetItem* item)
{
    if (ClickExec)
        ClickExec(item->text().toInt());
};
void ListSelectionMenu::rigthList2Left()
{
    //	HighLigthBoundary(0);
    for (int i = 0; i < rigthList->count(); i++)
    {
        if (rigthList->item(i)->isSelected())
        {

            QListWidgetItem* newitem = new QListWidgetItem();
            newitem->setText(rigthList->item(i)->text());
            //	rigthList->addItem(leftList->item(i));
            leftList->addItem(newitem);
            delete rigthList->item(i);
            //			ShowBoundaries(List2Vector(), redcolor, mainVTK);

            break;
        };
    };
    if (listExec)
        listExec();
};
void ListSelectionMenu::leftList2rigth()
{
    for (int i = 0; i < leftList->count(); i++)
    {
        if (leftList->item(i)->isSelected())
        {

            QListWidgetItem* newitem = new QListWidgetItem();
            newitem->setText(leftList->item(i)->text());
            //	rigthList->addItem(leftList->item(i));
            rigthList->addItem(newitem);
            delete leftList->item(i);
            //		ShowBoundaries(List2Vector(), redcolor, mainVTK);
            break;
        };
    };
    //	ShowBoundaries(List2Vector(), redcolor, mainVTK);
    if (listExec)
        listExec();
};
void ListSelectionMenu::GetListsData(std::vector<int>& list1, std::vector<int>& list2)
{
    list1.clear();
    for (int i = 0; i < leftList->count(); i++)
        list1.push_back(leftList->item(i)->text().toInt());

    list2.clear();
    for (int i = 0; i < rigthList->count(); i++)
        list2.push_back(rigthList->item(i)->text().toInt());
};
void ListSelectionMenu::GetListsData(std::vector<std::string>& list1, std::vector<std::string>& list2)
{
    list1.clear();
    for (int i = 0; i < leftList->count(); i++)
        list1.push_back(leftList->item(i)->text().toStdString());

    list2.clear();
    for (int i = 0; i < rigthList->count(); i++)
        list2.push_back(rigthList->item(i)->text().toStdString());
};
QGroupBox* ListSelectionMenu::GetPointer()
{
    return out;
};
