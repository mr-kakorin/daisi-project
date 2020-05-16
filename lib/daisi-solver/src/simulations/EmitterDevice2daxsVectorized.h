#include <cmath>
template <class PointType>
void EmitterDevice2daxsVectorized(PointType* a, PointType* b);
float IntegrateCurrent(float a, float b, float t0);
double IntegrateEnergy(double a, double b, double t0);
double IntegrateVelocity(double a, double b, double t0, double m, double charge, int thread);