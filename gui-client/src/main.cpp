#include "DaisiClient.h"
#include "startWindow.h"
#include <QApplication>

#undef RELATIVE
#undef ABSOLUTE

char                line[256];
std::vector<double> flags;

void readini()
{

    FILE* fp = fopen("daisi.ini", "r");
    while (fgets(line, 256, fp))
    {
        int tmp;
        sscanf(line, "%d", &tmp);
        flags.push_back(double(tmp));
    };
    fclose(fp);
};

int main(int argc, char* argv[])
{
    QStringList paths = QCoreApplication::libraryPaths();
    paths.append(".");
    paths.append("imageformats");
    paths.append("platforms");
    paths.append("sqldrivers");
    QCoreApplication::setLibraryPaths(paths);
    qApp->addLibraryPath("plugins/");

    QApplication a(argc, argv);
    Daizy        w;

    // QPalette  palette = w.palette();
    // palette.setColor(QPalette::Base, QColor(128, 128, 128));
    // w.setPalette(palette);

    QIcon icon("ico.ico");
    w.setWindowIcon(icon);
    w.show();
// readini();
// startWindow strtWin;
// strtWin.Create(flags, icon);

#if DEBUG
//	w.move(-1550,160);
#endif
    return a.exec();
}
