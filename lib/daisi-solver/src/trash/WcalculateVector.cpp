#include "WcalculateVector.h"
#include "mkl_vml.h"
void test(double* x1, double* x2)
{
#pragma simd
    for (int i = 0; i < 100; i++)
    {
        x1[i] = x1[i] + x2[i];
    }
}

template void WcalculateVector<float>(int* cells, int* index, float* x1, float* x2, float (*W)[9], int size,
                                      const float* x1Array, const float* x2Array, unsigned int* levelHigh);
template void WcalculateVector<double>(int* cells, int* index, double* x1, double* x2, double (*W)[9], int size,
                                       const double* x1Array, const double* x2Array, unsigned int* levelHigh);

template void WcalculateVector1<float>(Dmath::imat& tempNumb, char* cellType, unsigned int* cells, float* x1, float* x2,
                                       float (*W)[9], int size, const float* x1Array, const float* x2Array,
                                       unsigned int* levelHigh);
template void WcalculateVector1<double>(Dmath::imat& tempNumb, char* cellType, unsigned int* cells, double* x1,
                                        double* x2, double (*W)[9], int size, const double* x1Array,
                                        const double* x2Array, unsigned int* levelHigh);

template <class PointType>
void WcalculateVector(int* cells, int* index, PointType* x1, PointType* x2, PointType (*W)[9], int size,
                      const PointType* x1Array, const PointType* x2Array, unsigned int* levelHigh)
{
#pragma simd
    for (int i = 0; i < size; i++)
    {
        PointType H2 = x2Array[levelHigh[cells[i]]] - x2Array[cells[i]];
        PointType H1 = x1Array[cells[i] + 1] - x1Array[cells[i]];

        PointType X1 = (x1Array[cells[i]] - x1[index[i]]) / H1;

        PointType wx1_1 = PointType(0.5 * (0.5 + X1) * (0.5 + X1));
        PointType wx11  = PointType(0.5 * (0.5 - X1) * (0.5 - X1));
        PointType wx10  = PointType(1 - wx1_1 - wx11);
        PointType wx10T = PointType(0.75 - X1 * X1);

        PointType X2 = (x2Array[cells[i]] - x2[index[i]]) / H2;

        PointType wx2_1 = PointType(0.5 * (0.5 + X2) * (0.5 + X2));
        PointType wx21  = PointType(0.5 * (0.5 - X2) * (0.5 - X2));
        PointType wx20  = PointType(1 - wx2_1 - wx21);
        PointType wx20T = PointType(0.75 - X2 * X2);

        W[index[i]][0] = wx2_1 * wx1_1;
        W[index[i]][1] = wx2_1 * wx10;
        W[index[i]][2] = wx2_1 * wx11;
        W[index[i]][3] = wx20 * wx1_1;
        W[index[i]][4] = wx20 * wx10;
        W[index[i]][5] = wx20 * wx11;
        W[index[i]][6] = wx21 * wx1_1;
        W[index[i]][7] = wx21 * wx10;
        W[index[i]][8] = wx21 * wx11;
    }
}
template <class PointType>
void WcalculateVector1(Dmath::imat& tempNumb, char* cellType, unsigned int* cells, PointType* x1, PointType* x2,
                       PointType (*W)[9], int size, const PointType* x1Array, const PointType* x2Array,
                       unsigned int* levelHigh)
{
    int* d = &tempNumb.data[0];
#pragma simd
    for (int i = 0; i < size; i++)
    {
        int index = d[cells[i]];

        if (cellType[index] != 0)
            continue;

        PointType H2 = x2Array[levelHigh[index]] - x2Array[index];
        PointType H1 = x1Array[index + 1] - x1Array[index];

        PointType X1 = (x1Array[index] - x1[i]) / H1;

        PointType wx1_1 = PointType(0.5 * (0.5 + X1) * (0.5 + X1));
        PointType wx11  = PointType(0.5 * (0.5 - X1) * (0.5 - X1));
        PointType wx10  = PointType(1 - wx1_1 - wx11);
        PointType wx10T = PointType(0.75 - X1 * X1);

        PointType X2 = (x2Array[index] - x2[i]) / H2;

        PointType wx2_1 = PointType(0.5 * (0.5 + X2) * (0.5 + X2));
        PointType wx21  = PointType(0.5 * (0.5 - X2) * (0.5 - X2));
        PointType wx20  = PointType(1 - wx2_1 - wx21);
        PointType wx20T = PointType(0.75 - X2 * X2);

        W[i][0] = wx2_1 * wx1_1;
        W[i][1] = wx2_1 * wx10;
        W[i][2] = wx2_1 * wx11;
        W[i][3] = wx20 * wx1_1;
        W[i][4] = wx20 * wx10;
        W[i][5] = wx20 * wx11;
        W[i][6] = wx21 * wx1_1;
        W[i][7] = wx21 * wx10;
        W[i][8] = wx21 * wx11;
    }
}

template void calc<float>(int* levelHigh, int i1, int i2, int emType, arma::mat& Icoef, float* kI, char* flag,
                          unsigned int* startCell, float* r1, float* z1, int* index1, float (*W1)[9], float* r2,
                          float* z2, int* index2, float (*W2)[9], float* dr, float* dz, float* I, float dt, float* rho,
                          float* cartesianX1, float* cartesianX2, int flagType, float* cartesianX12,
                          float* cartesianX22);
template void calc<double>(int* levelHigh, int i1, int i2, int emType, arma::mat& Icoef, double* kI, char* flag,
                           unsigned int* startCell, double* r1, double* z1, int* index1, double (*W1)[9], double* r2,
                           double* z2, int* index2, double (*W2)[9], double* dr, double* dz, double* I, double dt,
                           double* rho, double* cartesianX1, double* cartesianX2, int flagType, double* cartesianX12,
                           double* cartesianX22);

template <class PointType>
void calc(int* levelHigh, int i1, int i2, int emType, arma::mat& Icoef, PointType* kI, char* flag,
          unsigned int* startCell, PointType* r1, PointType* z1, int* index1, PointType (*W1)[9], PointType* r2,
          PointType* z2, int* index2, PointType (*W2)[9], PointType* dr, PointType* dz, PointType* I, PointType dt,
          PointType* rho, PointType* cartesianX1, PointType* cartesianX2, int flagType, PointType* cartesianX12,
          PointType* cartesianX22)
{
#pragma vector always
    for (int i = i1; i < i2; i++)
    {
        // if (index1[i - i1] == index2[i - i1])
        //{
        rho[index1[i - i1]]     = rho[index1[i - i1]] + (W1[i - i1][0] + W2[i - i1][0]) * I[i] * dt / 2;
        rho[index1[i - i1] + 1] = rho[index1[i - i1] + 1] + (W1[i - i1][1] + W2[i - i1][1]) * I[i] * dt / 2;
        rho[levelHigh[index1[i - i1]]] =
            rho[levelHigh[index1[i - i1]]] + (W1[i - i1][2] + W2[i - i1][2]) * I[i] * dt / 2;
        rho[levelHigh[index1[i - i1]] + 1] =
            rho[levelHigh[index1[i - i1]] + 1] + (W1[i - i1][3] + W2[i - i1][3]) * I[i] * dt / 2;
        //};
    };
};