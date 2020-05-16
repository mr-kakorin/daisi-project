#include "Dmath.h"
template <class PointType>
void WcalculateVector(int* cells, int* index, PointType* x1, PointType* x2, PointType (*W)[9], int size,
                      const PointType* x1Array, const PointType* x2Array, unsigned int* levelHigh);
template <class PointType>
void WcalculateVector1(Dmath::imat& tempNumb, char* cellType, unsigned int* cells, PointType* x1, PointType* x2,
                       PointType (*W)[9], int size, const PointType* x1Array, const PointType* x2Array,
                       unsigned int* levelHigh);
void test(double* x1, double* x2);
template <class PointType>
void calc(int* levelHigh, int i1, int i2, int emType, arma::mat& Icoef, PointType* kI, char* flag,
          unsigned int* startCell, PointType* r1, PointType* z1, int* index1, PointType (*W1)[9], PointType* r2,
          PointType* z2, int* index2, PointType (*W2)[9], PointType* dr, PointType* dz, PointType* I, PointType dt,
          PointType* rho, PointType* cartesianX1, PointType* cartesianX2, int flagType, PointType* cartesianX12,
          PointType* cartesianX22);