#include "Dmath.h"

// template Dmath::CRSMatrix<double>;
// template Dmath::CRSMatrix<float>;

/*template void Dmath::erase<float>(std::vector<float>& data, std::vector<unsigned int>& ind);
template void Dmath::erase<double>(std::vector<double>& data, std::vector<unsigned int>& ind);
template void Dmath::erase<int>(std::vector<int>& data, std::vector<unsigned int>& ind);
template void Dmath::erase<char>(std::vector<char>& data, std::vector<unsigned int>& ind);
template void Dmath::erase<unsigned int>(std::vector<unsigned int>& data, std::vector<unsigned int>&
ind);

template<class PointType>
void erase(std::vector<PointType>& data, std::vector<unsigned int>& ind)
{
        for (int i = 0; i < ind.size(); i++)
        {
                data[ind[i]] = -1;
        };

        int flagBreak = 0;
        for (int i = 0; i < data.size(); i++)
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
        data.erase(data.begin() + data.size() - ind.size(), data.begin() + data.size());
};*/

double Dmath::Interpolate(const std::vector<double>& X, const std::vector<double>& Y, double x)
{
    double res;

    int I = -1;

    if (x <= X[0])
        I = 0;

    if (x >= X.back())
        I = X.size() - 2;

    if (-1 == I)
    {
        for (int i = 0; i < X.size() - 1; i++)
        {
            if (x >= X[i] && x <= X[i + 1])
            {
                I = i;
                break;
            }
        }
    }

    double w1 = (X[I + 1] - x) / (X[I + 1] - X[I]);
    double w2 = 1 - w1;
    res       = w1 * Y[I] + w2 * Y[I + 1];
    return res;
}

template <class PointType>
PointType Dmath::CRSMatrix<PointType>::getElement(int i, int j)
{
    int k1 = rowIndex[i];
    int k2;
    if (i >= rowIndex.size() - 1)
    {
        k2 = value.size();
    }
    else
    {
        k2 = rowIndex[i + 1];
    }

    if (k1 == k2)
        return 0.0;

    for (int k = k1; k < k2; k++)
    {
        if (j == col[k])
            return value[k];
    }
    return 0.0;
}

template <class PointType>
int Dmath::CRSMatrix<PointType>::getStrSize(int i)
{
    int begin = rowIndex[i];
    int end;
    if (i >= rowIndex.size() - 1)
    {
        end = value.size();
    }
    else
    {
        end = rowIndex[i + 1];
    }
    return end - begin - 1;
}

template <class PointType>
void Dmath::CRSMatrix<PointType>::setElement(int i, int j, PointType val)
{

    int k1 = rowIndex[i];
    int k2;
    if (i >= rowIndex.size() - 1)
    {
        k2 = value.size();
    }
    else
    {
        k2 = rowIndex[i + 1];
    }
    int k = k1;

    for (k = k1; (k < k2) && (col[k] <= j); k++)
    {
        if (col[k] == j)
        {
            value[k] = val;
            return;
        }
    }
    value.insert(value.begin() + k, val);
    col.insert(col.begin() + k, j);
    for (int m = i + 1; m < rowIndex.size(); ++m)
    {
        rowIndex[m]++;
    }
    return;
}

template <class PointType>
void Dmath::CRSMatrix<PointType>::initMatrix(int size)
{
    for (int i = 0; i < size; ++i)
    {
        rowIndex.push_back(i);
        col.push_back(0);
        value.push_back(0);
    }
}

template <class PointType>
void Dmath::CRSMatrix<PointType>::returnColIndexes(int i, std::vector<int>& result)
{

    int begin = rowIndex[i];
    int end;
    if (i >= rowIndex.size() - 1)
    {
        end = value.size();
    }
    else
    {
        end = rowIndex[i + 1];
    }
    result.resize(end - begin, 0);
    int m = 0;
    for (int j = begin; j < end; ++j)
    {
        result[m] = col[j];
        ++m;
    }
}

template <class PointType>
void Dmath::CRSMatrix<PointType>::destruct()
{
    value.clear();
    col.clear();
    rowIndex.clear();
}


template <class PointType>
void Dmath::Cartesian2Polar(PointType x, PointType y, PointType& r, PointType& phi)
{
    r = sqrt(double(x * x + y * y));
    if (x > 0 && y >= 0)
    {
        phi = atan(y / x);
        phi = 2 * PI() - phi;
        return;
    }
    if (x > 0 && y < 0)
    {
        phi = atan(y / x) + 2 * PI();
        phi = 2 * PI() - phi;
        return;
    }
    if (x < 0)
    {
        phi = atan(y / x) + PI();
        phi = 2 * PI() - phi;
        return;
    }
    if (std::abs(x) < 1e-7 && y > 0)
    {
        phi = PI() / 2;
        phi = 2 * PI() - phi;
        return;
    }
    if (std::abs(x) < 1e-7 && y < 0)
    {
        phi = 3 * PI() / 2;
        phi = 2 * PI() - phi;
        return;
    }
}

template <class PointType>
void Dmath::Polar2Cartesian(PointType* r, PointType* phi, PointType* x, PointType* y, int size)
{
    for (int i = 0; i < size; i++)
    {
        x[i] = r[i] * std::cos(double(2 * PI() - phi[i]));
        y[i] = r[i] * std::sin(double(2 * PI() - phi[i]));
    }
}

template <class PointType>
void Dmath::Polar2Cartesian(PointType r, PointType phi, PointType& x, PointType& y)
{
    x = r * std::cos(double(2 * PI() - phi));
    y = r * std::sin(double(2 * PI() - phi));
}

template <class PointType>
void Dmath::SmoothMovingAverage(std::vector<PointType>& data, int N)
{
    std::vector<PointType> dataTmp = data;
    for (int i = 0; i < data.size(); i++)
    {
        int    border = std::min(std::min(int(N), int(i)), int(data.size() - i - 1));
        double tmp    = 0;
        for (int ii = -border; ii <= border; ii++)
            tmp = tmp + dataTmp[i + border];

        data[i] = tmp / (2.0 * double(border) + 1.0);
    }
}

int Dmath::imat::getElemIndex(int i, int j)
{
    return i + j * nrow;
}

std::vector<int> Dmath::imat::GetData()
{
    return data;
}

int Dmath::imat::getElem(int base, int d1, int d2)
{
    if (base + d1 + d2 * nrow < 0)
        return -1;
    if (base + d1 + d2 * nrow >= size)
        return -1;

    return data[base + d1 + d2 * nrow];
}

int Dmath::imat::GetNrow()
{
    return nrow;
}

int Dmath::imat::GetNcol()
{
    return ncol;
}

int Dmath::imat::GetNz()
{
    return nz;
}

int Dmath::imat::ArSize() const
{
    return int(data.size());
}

void Dmath::imat::ones(int nrowIn, int ncolIn)
{
    nrow = nrowIn;
    ncol = ncolIn;
    size = nrowIn * ncolIn;
    data.clear();
    data = std::vector<int>(size);
    for (int i  = 0; i < size; i++)
        data[i] = 1;
}

Dmath::imat Dmath::imat::operator*(int right)
{
    imat result = *this;
    for (int i         = 0; i < size; i++)
        result.data[i] = result.data[i] * right;
    return result;
}

Dmath::imat Dmath::imat::operator+(int right)
{
    imat result = *this;
    for (int i = 0; i < size; i++)
    {
        if (result.data[i] != -1)
            result.data[i] = result.data[i] + right;
    }
    return result;
}

int& Dmath::imat::operator()(int i, int j)
{
    return data[i + j * nrow];
}

int Dmath::imat::getElemIndex(int i, int j, int k)
{
    return i + j * nrow + k * nrow * ncol;
}

int Dmath::imat::getElem(int base, int d1, int d2, int d3)
{
    if (base + d1 + d2 * nrow + d3 * nrow * ncol < 0)
        return -1;
    if (base + d1 + d2 * nrow + d3 * nrow * ncol >= size)
        return -1;

    return data[base + d1 + d2 * nrow + d3 * nrow * ncol];
}

void Dmath::imat::clear()
{
    nrow = 0;
    ncol = 0;
    nz   = 0;
    size = 0;
    data.clear();
}

void Dmath::imat::Set(const imat& right, int startPos)
{
    for (int i             = 0; i < right.ArSize(); i++)
        data[i + startPos] = right.data[i];
}

void Dmath::imat::ones(int nrowIn, int ncolIn, int nzIn)
{
    nrow = nrowIn;
    ncol = ncolIn;
    nz   = nzIn;
    size = nrowIn * ncolIn * nzIn;
    data.clear();
    data = std::vector<int>(size);
    for (int i  = 0; i < size; i++)
        data[i] = 1;
}

int& Dmath::imat::operator()(int i, int j, int k)
{
    return data[i + j * nrow + k * nrow * ncol];
}

Dmath::imat& Dmath::imat::operator=(const imat& right)
{
    if (this == &right)
    {
        return *this;
    }
    data = right.data;
    nrow = right.nrow;
    ncol = right.ncol;
    nz   = right.nz;
    size = right.size;
    return *this;
}

template void Dmath::Cartesian2Polar<float>(float x, float y, float& r, float& phi);
template void Dmath::Cartesian2Polar<double>(double x, double y, double& r, double& phi);

template void Dmath::Polar2Cartesian<float>(float* r, float* phi, float* x, float* y, int size);
template void Dmath::Polar2Cartesian<double>(double* r, double* phi, double* x, double* y,
                                             int size);
template void Dmath::Polar2Cartesian<float>(float r, float phi, float& x, float& y);
template void Dmath::Polar2Cartesian<double>(double r, double phi, double& x, double& y);

template void Dmath::SmoothMovingAverage<float>(std::vector<float>& data, int N);
template void Dmath::SmoothMovingAverage<double>(std::vector<double>& data, int N);
