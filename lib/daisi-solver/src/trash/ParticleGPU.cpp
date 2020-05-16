#include "ParticleGPU.h"

/*Particles3dcilGPU::Particles3dcilGPU(int size)
{
        nParticles = size;
        cudaMalloc((void **)&positions, size*sizeof(float3));
        cudaMalloc((void **)&momentums, size*sizeof(float3));
}
Particles3dcilGPU* Particles3dcilGPU::AllocateOnGPU()
{
        Particles3dcilGPU *result;
        cudaMalloc((void**)&result, sizeof(Particles3dcilGPU));
        cudaMemcpy(result, this, sizeof(Particles3dcilGPU), cudaMemcpyHostToDevice);
        return result;
};*/