#include <QtWidgets/QCheckBox>
#include <QtWidgets/QLabel>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

#ifdef NUCL
#include "namesNucl.h"
#else
#include "names.h"
#endif

class startWindow : public QWidget
{
    Q_OBJECT
    QWidget*     window;
    QCheckBox*   chb;
    QPushButton* button1;
    QLabel*      l;
    QVBoxLayout* vbl;

  public:
    void Create(const std::vector<double>& flags, QIcon icon)
    {
        if (flags[0])
            return;

        window = new QWidget();
        window->resize(160, 150);
        window->show();
        window->setWindowTitle("Information");
        window->setWindowIcon(icon);

        button1 = new QPushButton("ok");
        button1->setMaximumSize(70, 35);

        l = new QLabel(aboutStr0.c_str());

        vbl = new QVBoxLayout();
        window->setLayout(vbl);
        vbl->insertWidget(0, l, 0, Qt::AlignCenter);
        vbl->insertWidget(2, button1, 0, Qt::AlignCenter);

        chb = new QCheckBox("Do not show this message again");
        vbl->insertWidget(1, chb, 0, Qt::AlignLeft);
        chb->setChecked(flags[0]);

        connect(button1, SIGNAL(clicked()), this, SLOT(buttonClick()));
    };
  public slots:
    void buttonClick()
    {
        FILE* fp = fopen("daisi.ini", "w");

        int flag;
        if (chb->isChecked())
            flag = 1;
        else
            flag = 0;
        fprintf(fp, "%d", flag);

        fclose(fp);
        window->close();
    }
};
