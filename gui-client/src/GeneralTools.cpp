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

#include <notk/controller.hpp>

std::vector<float> linspace(float xmin, float xmax, size_t n_fm)
{
    float d = (xmax - xmin) / (float(n_fm) - 1.0);

    std::vector<float> xfm(n_fm);

    for (size_t ii = 0; ii < n_fm; ii++)
    {
        xfm[ii] = xmin + ii * d;
    }
    return xfm;
}

std::string tmpStr;

std::vector<float> gaussmf(const std::vector<float>& x, const float mu, const float sigma)
{
    if (std::abs(sigma) < 1e-20)
    {
        throw std::runtime_error("gaussmf::too small sigma");
    }
    std::vector<float> result(x.size());

    float koef = 1.0 / (sigma * std::sqrt(2.0 * 3.14159265359));

    for (size_t i = 0; i < x.size(); i++)
    {
        result[i] = koef * std::exp(-(x[i] - mu) * (x[i] - mu) / (2.0 * sigma * sigma));
    }
    return result;
}

float fitness(const std::vector<float>& x, std::vector<float>& y, const std::vector<float>& param)
{
    std::vector<float> y_  = gaussmf(x, param[0], param[1]);
    float              fit = 0;
    for (size_t i = 0; i < y.size(); i++)
    {
        fit = fit + (y[i] - y_[i]) * (y[i] - y_[i]);
    }
    return fit;
}

std::vector<float> gauss_fit(const std::vector<float>& x, const std::vector<float>& y, const double requiredAccuracy,
                             const std::vector<float>& x1, const std::vector<float>& x2)
{
    notk::NOTKController<float, float> NOTKController;

    std::string config = "{\"version\": \"0.2\",\"log_each_iteration\": \"10\", \"allow_maximization\": true, \"initial_search_borders\": {\"type\": \"vector\"},\"methods\": "
                         "[{\"type\": \"Gauss search\", \"use_random\": true, \"maximal_iterations\": 100, "
                         "\"maximal_preprocess_iterations\": 10, \"required_accuracy\": 1e-3, \"accuracy_type\": "
                         "\"relative\", \"accuracy_source_type\": \"fitness\",\"searcher1d\": {\"type\": \"Uniform "
                         "1d search\", \"maximal_iterations\": 10, \"required_accuracy\": 1e-3, "
                         "\"n_divisions_first\": 20, \"n_divisions\": 10, \"accuracy_type\": \"relative\", "
                         "\"accuracy_source_type\": \"fitness\" }}]}";

    NOTKController.set_problem_config_str(config);
    NOTKController.set_borders_fitness(std::bind(fitness, x, y, std::placeholders::_1), x1, x2);
    bool flag_abort = true;

    auto result = NOTKController.process(flag_abort);

    /*const OptimizationAccuracy accuracy_mode(OptimizationAccuracy::AccuracyType::RELATIVE,
            OptimizationAccuracy::AccuracySourceType::FITNESS,
            requiredAccuracy);

    auto search1d = std::bind(uniformSearch1d, maxNumIterationsInner, nDivisions * 2, nDivisions,
            std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
            accuracy_mode);

    auto fit = std::bind(fitness, std::ref(x), std::ref(y), std::placeholders::_1);

    auto results = gaussMultiDimSearch(100, 0.8, x1, x2, fit, search1d, accuracy_mode);

    return results.best_fitness_argument;*/
    return result->get_last_argument();
}

std::vector<float> estimate_gauss(const std::vector<float>& x, const std::vector<float>& y,
                                  const double requiredAccuracy)
{
    std::vector<float> x1(2);
    std::vector<float> x2(2);
    std::vector<float> x_s = x;
    std::vector<float> y_s = y;

    for (size_t i = 0; i < x_s.size(); i++)
    {
        if (std::abs(y_s[i]) < 1e-17)
        {
            x_s.erase(x_s.begin() + i);
            y_s.erase(y_s.begin() + i);

            i--;
        }
    }

    x1[0] = *std::min_element(x_s.begin(), x_s.end());
    x2[0] = *std::max_element(x_s.begin(), x_s.end());

    float x_av = (x1[0] + x2[0]) / 2.0;

    x1[1] = 1e-9;
    x2[1] = 2 * (x2[0] - x_av);

    return gauss_fit(x, y, requiredAccuracy, x1, x2);
}

bool convertChart(const std::vector<float>& XChart, const std::vector<float>& YChart, std::vector<float>& XChart_out,
                  std::vector<float>& YChart_out)
{
    for (size_t i = 0; i < XChart.size() - 1; i++)
    {
        XChart_out.push_back(XChart[i]);
        XChart_out.push_back(XChart[i + 1]);
        YChart_out.push_back(YChart[i]);
        YChart_out.push_back(YChart[i]);
    }
    return true;
}

bool getChart(std::vector<float>& data, int nCharts, std::vector<float>& XChart, std::vector<float>& YChart, int flag,
              double wmax_in)
{
    if (data.size() == 0)
        return false;

    volatile float wmin = data[0];
    volatile float wmax = data[0];

    float sum = 0;

    for (int i = 0; i < data.size(); i++)
    {
        sum = sum + data[i];
        if (data[i] > wmax)
            wmax = data[i];
        if (data[i] < wmin)
            wmin = data[i];
    };

    if (std::abs(wmax - wmin) < 1e-7)
        return false;

    if (wmax_in != -1)
    {
        wmax = wmax_in;
    }

    if (wmin < 0)
        wmin = wmin * 1.04;
    else
        wmin = wmin * 0.96;

    if (wmax > 0)
        wmax = wmax * 1.04;
    else
        wmax = wmax * 0.96;

    float d = (wmax - wmin) / float(nCharts);
    XChart.resize(nCharts + 1);
    YChart.resize(nCharts + 1);

    XChart[0]     = wmin;
    XChart.back() = wmax;

    for (int i = 1; i < nCharts; i++)
        XChart[i] = wmin + i * d;

    for (int i = 1; i < nCharts; i++)
    {
        for (int j = 0; j < data.size(); j++)
        {
            if (flag == 0)
            {
                if (data[j] > XChart[i] && data[j] < XChart[i + 1])
                    YChart[i]++;
            }
            if (flag == 1)
            {
                if (data[j] < XChart[i])
                    YChart[i]++;
            }
        }
    }
    for (auto& y_ch : YChart)
    {
        y_ch = 100.0 * y_ch / data.size();
    }
    return true;
}

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
std::string GetFileName(const std::string& fullPath)
{
    std::string result = "";
    int         i      = fullPath.size() - 1;
    if (i == -1)
        return result;
    while (1)
    {
        if ((fullPath[i] == '/' && fullPath[i - 1] == '/') || (fullPath[i] == '\\'))
        {
            result = fullPath.substr(i + 1, fullPath.size() - 1);
            return result;
        };
        i--;
        if (i == 0)
            break;
    };
    return result;
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