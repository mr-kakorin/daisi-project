#ifndef GRIDDATAGPU_H
#define GRIDDATAGPU_H
#include <cuda_runtime.h>
template <class DataContainer> class GridData2d;
template <class DataContainer> class GridData2daxs;
template <class DataContainer> class GridData2dpolar;
class GridDataBaseGPU
{
  protected:
    int nPoints;
};
template <class PointType> class __align__(128) GridData1dGPU : GridDataBaseGPU
{
  public:
    PointType* x;
    PointType* Ex;
    PointType* V;
    PointType* rho;
};
template <class PointType> class __align__(128) GridData2dGPU : GridDataBaseGPU
{
  public:
    PointType* x;
    PointType* y;
    PointType* Ex;
    PointType* Ey;
    PointType* V;
    PointType* rho;
    GridData2dGPU(){

    };
    template <class DataContainer> GridData2dGPU(GridData2d<DataContainer> & obj)
    {
        nPoints = obj.x.size();
        int s   = nPoints * sizeof(PointType);

        cudaMalloc((void**)&x, s);
        cudaMemcpy(x, &obj.x[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&y, s);
        cudaMemcpy(y, &obj.y[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Ex, s);
        cudaMemcpy(Ex, &obj.Ex[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Ey, s);
        cudaMemcpy(Ey, &obj.Ey[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&V, s);
        cudaMemcpy(V, &obj.GetV()[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&rho, s);
        cudaMemcpy(rho, &obj.Getrho()[0], s, cudaMemcpyHostToDevice);
    };
};
template <class PointType> class __align__(128) GridData2daxsGPU : GridDataBaseGPU
{
  public:
    PointType* r;
    PointType* z;
    PointType* Er;
    PointType* Ez;
    PointType* Bphi;
    PointType* V;
    PointType* rho;
    GridData2daxsGPU(){

    };
    template <class DataContainer> GridData2daxsGPU(GridData2daxs<DataContainer> & obj)
    {
        nPoints = obj.r.size();
        int s   = nPoints * sizeof(PointType);

        cudaMalloc((void**)&r, s);
        cudaMemcpy(r, &obj.r[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&z, s);
        cudaMemcpy(z, &obj.z[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Er, s);
        cudaMemcpy(Er, &obj.Er[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Ez, s);
        cudaMemcpy(Ez, &obj.Ez[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Bphi, s);
        cudaMemcpy(Bphi, &obj.Bphi[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&V, s);
        cudaMemcpy(V, &obj.GetV()[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&rho, s);
        cudaMemcpy(rho, &obj.Getrho()[0], s, cudaMemcpyHostToDevice);
    };
};

template <class PointType> class __align__(128) GridData2dpolarGPU : GridDataBaseGPU
{
  public:
    PointType* r;
    PointType* phi;
    PointType* Er;
    PointType* Ephi;
    PointType* V;
    PointType* rho;
    GridData2dpolarGPU(){

    };
    template <class DataContainer> GridData2dpolarGPU(GridData2dpolar<DataContainer> & obj)
    {
        nPoints = obj.r.size();
        int s   = nPoints * sizeof(PointType);

        cudaMalloc((void**)&r, s);
        cudaMemcpy(r, &obj.r[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&phi, s);
        cudaMemcpy(phi, &obj.phi[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Er, s);
        cudaMemcpy(Er, &obj.Er[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&Ephi, s);
        cudaMemcpy(Ephi, &obj.Ephi[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&V, s);
        cudaMemcpy(V, &obj.GetV()[0], s, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&rho, s);
        cudaMemcpy(rho, &obj.Getrho()[0], s, cudaMemcpyHostToDevice);
    };
};

#endif