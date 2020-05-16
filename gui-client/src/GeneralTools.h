#ifndef GeneralTools_H
#define GeneralTools_H
#include <string>
#include <vector>
std::string GetFileName(const std::string& fullPath);
void vectorSubtraction(std::vector<int>& v1, std::vector<int> v2);
class QLayout;
void clearLayout(QLayout* layout, bool deleteWidgets = true);
void strsplit(char* string, char* split, std::vector<double>& result, std::string& error);
float vectorMax(std::vector<float>& input);
bool getChart(std::vector<float>& data, int nCharts, 
	std::vector<float>& XChart, std::vector<float>& YChart, int flag, double wmax_in = -1);
bool convertChart(const std::vector<float>& XChart, const std::vector<float>& YChart,
	std::vector<float>& XChart_out, std::vector<float>& YChart_out);
std::vector<float> estimate_gauss(const std::vector<float>& XChart, const std::vector<float>& YChart, const double requiredAccuracy);
std::vector<float> gaussmf(const std::vector<float>& x, const float a, const float mu, const float sigma);
std::vector<float> linspace(float xmin, float xmax, size_t n_fm);

#endif