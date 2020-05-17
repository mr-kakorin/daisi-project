#ifndef GeneralTools_H
#define GeneralTools_H
#include <string>
#include <vector>
void vectorSubtraction(std::vector<int>& v1, std::vector<int> v2);
class QLayout;
void clearLayout(QLayout* layout, bool deleteWidgets = true);
void strsplit(char* string, char* split, std::vector<double>& result, std::string& error);
float vectorMax(std::vector<float>& input);
#endif