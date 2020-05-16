#include "MyTreeItem.h"
double MyQTextEdit::toDoubleMy(bool* ok)
{
    QString str    = toPlainText();
    double  result = str.toDouble(ok);

    if (!(*ok))
    {
        blockSignals(true);
        QPalette palette = this->palette();
        palette.setColor(QPalette::Base, QColor(207, 206, 255));
        palette.setColor(QPalette::Text, QColor(Qt::red));
        setPalette(palette);
        blockSignals(false);
    }
    else
    {
        blockSignals(true);
        QPalette palette = this->palette();
        palette.setColor(QPalette::Base, QColor(Qt::white));
        palette.setColor(QPalette::Text, QColor(Qt::black));
        setPalette(palette);
        blockSignals(false);
    }

    return result;

    /*bool old = *ok;
    double result = 0;
    QString str = toPlainText();
    result = str.toDouble(ok);
    if (!(*ok))
    {
            this->clear();
            this->setText(str);
            this->setTextColor(Qt::red);
    };
    if (!old)
            (*ok) = false;
    return result;*/
};
