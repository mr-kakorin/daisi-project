#ifndef DMATH_H
#define DMATH_H
#include "armadillo"
#include "boost/serialization/vector.hpp"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <Constants.h>

#include <functional>
namespace Dmath
{
template <class T>
long long vectorsize(T vec)
{
    return sizeof(T);
}

template <class T>
long long vectorsize(std::vector<T>& vec)
{
    volatile long long result = 0;

    if (sizeof(vec[0]) > 8)
    {
        for (int i = 0; i < vec.size(); i++)
        {
            result = result + vectorsize(vec[i]);
        }
        return result;
    }
    else
        return result + sizeof(T) * vec.capacity() + sizeof(vec);
}

template <class T>
long long vectorsize(const std::vector<T>& vec)
{
    volatile long long result = 0;

    if (sizeof(vec[0]) > 8)
    {
        for (int i = 0; i < vec.size(); i++)
        {
            result = result + vectorsize(vec[i]);
        }
    }
    else
        return result + sizeof(T) * vec.capacity() + sizeof(vec);
}

template <class PointType>
inline int sign(PointType val)
{
    if (val == 0)
        return 0;
    if (val < 0)
        return -1;
    return 1;
}

template <class PointType>
class CRSMatrix
{
  private:
    std::vector<PointType> value;
    std::vector<int>       col;
    std::vector<int>       rowIndex;

  public:
    PointType getElement(int i, int j);

    void setElement(int i, int j, PointType val);
    void initMatrix(int size);
    void returnColIndexes(int i, std::vector<int>& result);
    void destruct();

    int getStrSize(int i);
};

double Interpolate(const std::vector<double>& X, const std::vector<double>& Y, double x);

template <class _Fn, class... _Args>
double integral(double a, double b, _Fn&& _Fx, _Args&&... _Ax)
{
    std::function<double(double)> function;
    function =
        std::bind(std::forward<_Fn>(_Fx), std::placeholders::_1, std::forward<_Args>(_Ax)...);

    int    s      = int((b - a) * 200 / PI());
    double result = 0;
    double step   = (b - a) / s;
    // double f1 = f(a, param);
    double f1       = function(a);
    double currentP = a;
    double f2;
    for (int i = 0; i < s; i++)
    {
        currentP = currentP + step;
        //	f2 = f(currentP, param);
        f2     = function(currentP);
        result = result + (f2 + f1 + 4 * function(currentP - step / 2)) / 6;
        f1     = f2;
    };
    result = result * step;
    return result;
}

template <class _Fn, class... _Args>
double zerosSearchRFQ(double X0, double tolerance, double facSeacr, _Fn&& _Fx, _Args&&... _Ax)
{
    std::function<double(double)> function;
    function =
        std::bind(std::forward<_Fn>(_Fx), std::placeholders::_1, std::forward<_Args>(_Ax)...);

    double val = function(X0);
    int    i   = 0;
    double x1, x2;
    while (1)
    {
        double val1 = function(X0 + facSeacr * X0 * i);
        double val2 = function(X0 - facSeacr * X0 * i);
        if (val1 * val < 0)
        {
            x1 = X0;
            x2 = X0 + facSeacr * X0 * i;
            break;
        }
        if (val2 * val < 0)
        {
            x1 = X0 - facSeacr * X0 * i;
            x2 = X0;
            break;
        }
        i++;
        if (i > 10)
            return -1;
    }

    while (1)
    {
        double val1 = function(x1);
        double val2 = function(x2);
        double x3   = (x1 + x2) / 2;
        double val3 = function(x3);

        if (val1 * val3 < 0)
            x2 = x3;
        else
            x1 = x3;

        if (std::abs(x1 - x2) < tolerance)
            return x3;
    }
}

template <class PointType>
void Cartesian2Polar(PointType x, PointType y, PointType& r, PointType& phi);

template <class PointType>
void Polar2Cartesian(PointType* r, PointType* phi, PointType* x, PointType* y, int size);
template <class PointType>
void Polar2Cartesian(PointType r, PointType phi, PointType& x, PointType& y);

template <class PointType>
void SmoothMovingAverage(std::vector<PointType>& data, int N);

// template<class PointType>
// void erase(std::vector<PointType>& data, std::vector<unsigned int>& ind);

template <class PointType>
void erase(std::vector<PointType>& data, std::vector<unsigned int>& ind)
{

    if (data.size() == 0)
        return;

    for (int i = 0; i < ind.size(); i++)
    {
        data[ind[i]] = -1;
    };

    /*int j;
    for (int i = 0; i < data.size() - 1; i++)
    {
            if (data[i] == -1)
            {
                    j = i;
                    while (data[j])
            }
    }*/
    /*int flagBreak = 0;
    for (int i = 0; i < data.size()-1; i++)
    {
            if (data[i] == -1)
            {
                    int j = i + 1;
                    while (1)
                    {
                            if (data[j] != -1)
                            {
                                    std::swap(data[i], data[j]);
                                    break;
                            };
                            j++;
                            if (j == data.size())
                            {
                                    flagBreak = 1;
                                    break;
                            }
                    };
                    if (flagBreak == 1)
                            break;
            };
    }
    data.erase(data.begin() + data.size() - ind.size(), data.begin() + data.size());*/
}

class imat
{
    friend class boost::serialization::access;

  private:
    int nrow;
    int ncol;
    int nz;
    int size;

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& nrow;
        ar& ncol;
        ar& size;
        ar& nz;
        ar& data;
    };
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& nrow;
        ar& ncol;
        ar& size;
        ar& nz;
        ar& data;
    };

  public:
    std::vector<int> data;

    int getElemIndex(int i, int j, int k);
    int getElemIndex(int i, int j);

    std::vector<int> GetData();

    int getElem(int base, int d1, int d2, int d3);
    int getElem(int base, int d1, int d2);
    int GetNz();

    int GetNrow();

    int GetNcol();

    int ArSize() const;

    void clear();

    void ones(int nrowIn, int ncolIn);
    void ones(int nrowIn, int ncolIn, int nzIn);

    imat operator*(int right);
    imat operator+(int right);

    void Set(const imat& right, int startPos);

    int& operator()(int i, int j, int k);
    int& operator()(int i, int j);
    imat& operator=(const imat& right);
};

}
#endif
