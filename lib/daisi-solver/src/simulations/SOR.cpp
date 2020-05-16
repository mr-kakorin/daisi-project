#include "SOR.h"
#include <cmath>
template void SOR<float>(float* x, int* rights, float* c_rights, int* lefts, float* c_lefts, int* ups, float* c_ups,
                         int* downs, float* c_downs, int* middles, float* c_middles, float w, int size, float* vect,
                         double& diff_sum, float* tmp);
template void SOR<double>(double* x, int* rights, double* c_rights, int* lefts, double* c_lefts, int* ups,
                          double* c_ups, int* downs, double* c_downs, int* middles, double* c_middles, double w,
                          int size, double* vect, double& diff_sum, double* tmp);

template <class PointType>
void SOR(PointType* x, int* rights, PointType* c_rights, int* lefts, PointType* c_lefts, int* ups, PointType* c_ups,
         int* downs, PointType* c_downs, int* middles, PointType* c_middles, PointType w, int size, PointType* vect,
         double& diff_sum, PointType* tmp)
{
    PointType tmp_sum = 0;
    int       n;
    diff_sum = -1;
    PointType tol;
    for (int k = 0; k < size; ++k)
    {
        tmp_sum = 0;
        n       = middles[k];
        tmp_sum =
            -(x[rights[k]] * c_rights[k] + x[lefts[k]] * c_lefts[k] + x[ups[k]] * c_ups[k] + x[downs[k]] * c_downs[k]) *
            w / c_middles[k];
        tmp_sum = tmp_sum + (1 - w) * x[n] + w * vect[n] / c_middles[k];
        //	diff_sum = diff_sum + std::abs(x[n] - tmp_sum);
        tol = std::abs((x[n] - tmp_sum) / tmp_sum);
        if (tmp_sum < 0.00001)
            tol = 0;

        if (tol > diff_sum)
            diff_sum = tol;

        x[n] = tmp_sum;
    }
};
