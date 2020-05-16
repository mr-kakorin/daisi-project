#ifndef PARTICLESGPU_H
#define PARTICLESGPU_H
#include "Particle.h"
#include <cuda_runtime.h>
class ParticlesBaseGPU
{
  protected:
    int    nParticles;
    double Time;
    int    AllocateSize;
};
template <class PointType> class Particles2dGPU : public ParticlesBaseGPU
{
  public:
    PointType* x;
    PointType* px;
    PointType* q;
    PointType* y;
    PointType* py;

    Particles2dGPU(){

    };
    template <class PointType> Particles2dGPU(Particles2d<PointType>& obj)
    {
        nParticles = obj.x.size();
        Time       = obj.Time;

        int nByte = 1e7;

        AllocateSize = nByte / sizeof(PointType);

        int sCopy = nParticles * sizeof(PointType);

        cudaMalloc((void**)&x, AllocateSize);
        cudaMemcpy(x, &obj.x[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&y, AllocateSize);
        cudaMemcpy(y, &obj.y[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&q, AllocateSize);
        cudaMemcpy(q, &obj.q[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&px, AllocateSize);
        cudaMemcpy(px, &obj.px[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&py, AllocateSize);
        cudaMemcpy(py, &obj.py[0], sCopy, cudaMemcpyHostToDevice);
    };
};
template <class PointType> class Particles3dcilGPU : public ParticlesBaseGPU
{
  public:
    PointType* r;
    PointType* pr;
    PointType* z;
    PointType* pz;
    PointType* phi;
    PointType* pphi;
    PointType* q;
    Particles3dcilGPU(){

    };
    template <class PointType> Particles3dcilGPU(Particles3dcil<PointType>& obj)
    {
        nParticles = obj.r.size();
        Time       = obj.Time;

        int nByte = 1e7;

        AllocateSize = nByte / sizeof(PointType);
        int sCopy    = nParticles * sizeof(PointType);

        cudaMalloc((void**)&r, AllocateSize);
        cudaMemcpy(r, &obj.r[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&z, AllocateSize);
        cudaMemcpy(z, &obj.z[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&pr, AllocateSize);
        cudaMemcpy(pr, &obj.pr[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&pz, AllocateSize);
        cudaMemcpy(pz, &obj.pz[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&phi, AllocateSize);
        cudaMemcpy(phi, &obj.phi[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&pphi, AllocateSize);
        cudaMemcpy(pphi, &obj.pphi[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&q, AllocateSize);
        cudaMemcpy(q, &obj.q[0], sCopy, cudaMemcpyHostToDevice);
    };
    Particles3dcilGPU(int size);
    Particles3dcilGPU* AllocateOnGPU();
};
template <class PointType> class Particles2dpolarGPU : public ParticlesBaseGPU
{
  public:
    PointType* r;
    PointType* pr;
    PointType* phi;
    PointType* pphi;
    PointType* q;
    Particles2dpolarGPU(){

    };
    template <class PointType> Particles2dpolarGPU(Particles2dpolar<PointType>& obj)
    {
        nParticles = obj.r.size();
        Time       = obj.Time;

        int nByte = 1e7;

        AllocateSize = nByte / sizeof(PointType);
        int sCopy    = nParticles * sizeof(PointType);

        cudaMalloc((void**)&r, AllocateSize);
        cudaMemcpy(r, &obj.r[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&pr, AllocateSize);
        cudaMemcpy(pr, &obj.pr[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&phi, AllocateSize);
        cudaMemcpy(phi, &obj.phi[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&pphi, AllocateSize);
        cudaMemcpy(pphi, &obj.pphi[0], sCopy, cudaMemcpyHostToDevice);

        cudaMalloc((void**)&q, AllocateSize);
        cudaMemcpy(q, &obj.q[0], sCopy, cudaMemcpyHostToDevice);
    };
    Particles2dpolarGPU(int size);
    Particles2dpolarGPU* AllocateOnGPU();
};
#endif