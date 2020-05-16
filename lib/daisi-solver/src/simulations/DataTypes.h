#ifndef DATATYPES_H
#define DATATYPES_H
#include "Device.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice2daxs.h"
#include "EmitterDevice3d.h"
#include "MeshContainer3d.h"
#include "ModelTemplate.h"

typedef DeviceStatus<Particles3dcil<float>, EmitterDevice2daxs<float>, GridData2daxs<float>, float,
                     BoundaryContainer2d<float>, MeshContainer2d<float>>
    device2daxsfloat;

typedef DeviceStatus<Particles3dcil<double>, EmitterDevice2daxs<double>, GridData2daxs<double>,
                     double, BoundaryContainer2d<double>, MeshContainer2d<double>>
    device2daxsdouble;

typedef DeviceStatus<Particles2d<float>, EmitterDevice2d<float>, GridData2d<float>, float,
                     BoundaryContainer2d<float>, MeshContainer2d<float>>
    device2dfloat;

typedef DeviceStatus<Particles2d<double>, EmitterDevice2d<double>, GridData2d<double>, double,
                     BoundaryContainer2d<double>, MeshContainer2d<double>>
    device2ddouble;

typedef DeviceStatus<Particles2dpolar<float>, EmitterDevice2d<float>, GridData2dpolar<float>, float,
                     BoundaryContainer2d<float>, MeshContainer2d<float>>
    device2dpolarfloat;
typedef DeviceStatus<Particles2dpolar<double>, EmitterDevice2d<double>, GridData2dpolar<double>,
                     double, BoundaryContainer2d<double>, MeshContainer2d<double>>
    device2dpolardouble;

typedef DeviceStatus<Particles3d<double>, EmitterDevice3d<double>, GridData3d<double>, double,
                     BoundaryContainer2d<double>, MeshContainer3d<double>>
    device3dExtrdouble;
typedef DeviceStatus<Particles3d<float>, EmitterDevice3d<float>, GridData3d<float>, float,
                     BoundaryContainer2d<float>, MeshContainer3d<float>>
    device3dExtrfloat;

typedef DeviceStatus<Particles3d<double>, EmitterDevice3d<double>, GridData3d<double>, double,
                     BoundaryContainer3d<double>, MeshContainer3d<double>>
    device3ddouble;

// DeviceStatus<Particles3d<float>, EmitterDevice3d<float>, const std::shared_ptr<GridData3d<float>,
// float, BoundaryContainer2d<float>, MeshContainer3d<float>> device3dfloat;

/*typedef DeviceStatusGPU<Particles3dcilGPU<float>, EmitterDevice2daxsGPU<float>,
GridData2daxsGPU<float>, float, BoundaryContainer2dGPU<float>> device2daxsfloatGPU; typedef
DeviceStatusGPU<Particles3dcilGPU<double>, EmitterDevice2daxsGPU<double>, GridData2daxsGPU<double>,
double, BoundaryContainer2dGPU<double>> device2daxsdoubleGPU;

typedef DeviceStatusGPU<Particles2dGPU<float>, EmitterDevice2dGPU<float>, GridData2dGPU<float>,
float, BoundaryContainer2dGPU<float>> device2dfloatGPU; typedef
DeviceStatusGPU<Particles2dGPU<double>, EmitterDevice2dGPU<double>, GridData2dGPU<double>, double,
BoundaryContainer2dGPU<double>> device2ddoubleGPU;

typedef DeviceStatusGPU<Particles2dpolarGPU<float>, EmitterDevice2dGPU<float>,
GridData2dpolarGPU<float>, float, BoundaryContainer2dGPU<float>> device2dpolarfloatGPU; typedef
DeviceStatusGPU<Particles2dpolarGPU<double>, EmitterDevice2dGPU<double>, GridData2dpolarGPU<double>,
double, BoundaryContainer2dGPU<double>> device2dpolardoubleGPU;*/

// typedef DeviceStatus<Particles3d<double>, EmitterDevice3d<double>,
// const std::shared_ptr<GridData3d<double>, double, BoundaryContainer3d<double>,
// MeshContainer3d<double>>
// RFQDesignModel;

/*typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> Solver2daxsfloat; typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> Solver2daxsfloat; typedef Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> Solver2daxsdouble; typedef Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
Solver2daxsdouble;


typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> Solver2dfloat; typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> Solver2dfloat; typedef Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> Solver2ddouble; typedef Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
Solver2ddouble;

typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> Solver2dpolarfloat; typedef Solver <ParticleGridInterface<float>,
EmissionCurrentSolver<float>, PoissonSolver<float>, float> Solver2dpolarfloat; typedef Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
Solver2dpolardouble; typedef Solver <ParticleGridInterface<double>, EmissionCurrentSolver<double>,
PoissonSolver<double>, double> Solver2dpolardouble;



typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> SolverPTI2dfloat; typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> SolverPIC2daxsfloat; typedef Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> SolverPTI2ddouble; typedef Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
SolverPIC2daxsdouble;

typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> SolverPIC2dfloat; typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> SolverPTI2daxsfloat; typedef Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> SolverPIC2ddouble; typedef Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
SolverPTI2daxsdouble;


typedef Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> SolverPIC2dpolarfloat; typedef Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> SolverPIC2dpolardouble; typedef Solver
<ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>, float>
SolverPTI2dpolarfloat; typedef Solver <ParticleGridInterface<double>, EmissionCurrentSolver<double>,
PoissonSolver<double>, double> SolverPTI2dpolardouble;*/

#endif