#include "GeneralTools.h"
#include <QtWidgets/QLayout>
#include <QtWidgets/QLayoutItem>
#include <QtWidgets/QWidget>
// #include <windows.h>
#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>

#undef RELATIVE
#undef ABSOLUTE


float vectorMax(std::vector<float>& input)
{
    float max = input[0];
    for (int i = 1; i < input.size(); i++)
    {
        if (input[i] > max)
            max = input[i];
    };
    return max;
};

void vectorSubtraction(std::vector<int>& v1, std::vector<int> v2)
{
    for (int i = 0; i < v2.size(); i++)
        v1[v2[i]] = -1;

    int k = int(v1.size());
    for (int i = 0; i < k; i++)
    {
        if (v1[i] == -1)
        {
            v1.erase(v1.begin() + i);
            i--;
            k--;
        };
    };
};
void clearLayout(QLayout* layout, bool deleteWidgets)
{
    while (QLayoutItem* item = layout->takeAt(0))
    {
        if (deleteWidgets)
        {
            if (QWidget* widget = item->widget())
                delete widget;
        }
        if (QLayout* childLayout = item->layout())
            clearLayout(childLayout, deleteWidgets);
        delete item;
    }
}

void strsplit(char* string, char* split, std::vector<double>& result, std::string& error)
{
    // __try
    // {
    std::string tmpStr;
    result.clear();
    char* pch = strtok(string, split);
    while (pch != NULL)
    {
        tmpStr = pch;
        result.push_back(std::stod(tmpStr));
        pch = strtok(NULL, split);
    }
    // }
    // __except (EXCEPTION_EXECUTE_HANDLER)
    // {
    //     error = "Invalid input character \n please use the following format: t1;t2;t3";
    //     return;
    // }
};