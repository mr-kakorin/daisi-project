#include "updatePositionsVector.h"
#include <cmath>


template void updatePositionsVector<float>(float* r, float* z, float* phi, float* pr, float* pz, float* pphi, int size,
                                           float timeStep, float* tmp1, float* tmp2);
template void updatePositionsVector<double>(double* r, double* z, double* phi, double* pr, double* pz, double* pphi,
                                            int size, double timeStep, double* tmp1, double* tmp2);

template <class PointType>
void updatePositionsVector(PointType* r1, PointType* z1, PointType* phi1, PointType* pr1, PointType* pz1,
                           PointType* pphi1, int size, PointType timeStep, PointType* tmp11, PointType* tmp21)
{
/*double*r = (double*)r1;
double*z = (double*)z1;
double*phi = (double*)phi1;
double*pr = (double*)pr1;
double*pz = (double*)pz1;
double*pphi = (double*)pphi1;
double*tmp1 = (double*)tmp11;
double*tmp2 = (double*)tmp21;

vdDiv(size, pphi, r, tmp1);
vdMul(size, tmp1, tmp1, tmp1);

vdMul(size, pz, pz, tmp2);
vdAdd(size, tmp2, tmp1, tmp1);

vdMul(size, pr, pr, tmp2);
vdAdd(size, tmp2, tmp1, tmp1);

#pragma simd
for (int i = 0; i < size; i++)
{
        tmp1[i] = tmp1[i] + 1;
}

vdSqrt(size, tmp1, tmp1);

vdDiv(size, pr, tmp1, tmp2);


cblas_daxpy(size, timeStep, tmp2, 1, r, 1);

vdDiv(size, pz, tmp1, tmp2);

cblas_daxpy(size, timeStep, tmp2, 1, z, 1);

vdDiv(size, pphi, tmp1, tmp1);
vdDiv(size, tmp1, r, tmp1);
vdDiv(size, tmp1, r, tmp1);

cblas_daxpy(size, timeStep, tmp1, 1, phi, 1);*/

#pragma simd
    for (int i = 0; i < size; i++)
    {
        PointType gamma = sqrt(1 + pr1[i] * pr1[i] + pz1[i] * pz1[i] + (pphi1[i] / r1[i]) * (pphi1[i] / r1[i]));
        PointType rr    = r1[i] + 0.5 * timeStep * pr1[i] / gamma;
        r1[i]           = r1[i] + timeStep * pr1[i] / gamma;
        z1[i]           = z1[i] + timeStep * pr1[i] / gamma;
        phi1[i]         = phi1[i] + timeStep * pphi1[i] / (gamma * rr * rr);
    };
};

template void updateMomentumsVector<float>(float* r, float* z, float* phi, float* pr, float* pz, float* pphi, int size,
                                           float timeStep, float* tmp1, float* tmp2);
template void updateMomentumsVector<double>(double* r, double* z, double* phi, double* pr, double* pz, double* pphi,
                                            int size, double timeStep, double* tmp1, double* tmp2);
template <class PointType>
void updateMomentumsVector(PointType* r1, PointType* z1, PointType* phi1, PointType* pr1, PointType* pz1,
                           PointType* pphi1, int size, PointType timeStep, PointType* tmp11, PointType* tmp21){};
