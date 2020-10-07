#ifndef DATATYPES_H
#define DATATYPES_H
#include "Device.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice2daxs.h"
#include "EmitterDevice3d.h"
#include "MeshContainer3d.h"
#include "ModelTemplate.h"

using device2daxsfloat = DeviceStatus<Particles3dcil<float>, EmitterDevice2daxs<float>, GridData2daxs<float>, float,
                     BoundaryContainer2d<float>, MeshContainer2d<float>>;

using device2daxsdouble = DeviceStatus<Particles3dcil<double>, EmitterDevice2daxs<double>, GridData2daxs<double>,
                     double, BoundaryContainer2d<double>, MeshContainer2d<double>>;

using device2dfloat = DeviceStatus<Particles2d<float>, EmitterDevice2d<float>, GridData2d<float>, float,
                     BoundaryContainer2d<float>, MeshContainer2d<float>>;

using device2ddouble = DeviceStatus<Particles2d<double>, EmitterDevice2d<double>, GridData2d<double>, double,
                     BoundaryContainer2d<double>, MeshContainer2d<double>>;

using device2dpolarfloat = DeviceStatus<Particles2dpolar<float>, EmitterDevice2d<float>, GridData2dpolar<float>, float,
                     BoundaryContainer2d<float>, MeshContainer2d<float>>;
using device2dpolardouble = DeviceStatus<Particles2dpolar<double>, EmitterDevice2d<double>, GridData2dpolar<double>,
                     double, BoundaryContainer2d<double>, MeshContainer2d<double>>;

using device3dExtrdouble = DeviceStatus<Particles3d<double>, EmitterDevice3d<double>, GridData3d<double>, double,
                     BoundaryContainer2d<double>, MeshContainer3d<double>>;
using device3dExtrfloat = DeviceStatus<Particles3d<float>, EmitterDevice3d<float>, GridData3d<float>, float,
                     BoundaryContainer2d<float>, MeshContainer3d<float>>;

using device3ddouble = DeviceStatus<Particles3d<double>, EmitterDevice3d<double>, GridData3d<double>, double,
                     BoundaryContainer3d<double>, MeshContainer3d<double>>;

// DeviceStatus<Particles3d<float>, EmitterDevice3d<float>, const std::shared_ptr<GridData3d<float>,
// float, BoundaryContainer2d<float>, MeshContainer3d<float>> device3dfloat;

/*using DeviceStatusGPU<Particles3dcilGPU<float>, EmitterDevice2daxsGPU<float>,
GridData2daxsGPU<float>, float, BoundaryContainer2dGPU<float>> device2daxsfloatGPU; using
DeviceStatusGPU<Particles3dcilGPU<double>, EmitterDevice2daxsGPU<double>, GridData2daxsGPU<double>,
double, BoundaryContainer2dGPU<double>> device2daxsdoubleGPU;

using DeviceStatusGPU<Particles2dGPU<float>, EmitterDevice2dGPU<float>, GridData2dGPU<float>,
float, BoundaryContainer2dGPU<float>> device2dfloatGPU; using
DeviceStatusGPU<Particles2dGPU<double>, EmitterDevice2dGPU<double>, GridData2dGPU<double>, double,
BoundaryContainer2dGPU<double>> device2ddoubleGPU;

using DeviceStatusGPU<Particles2dpolarGPU<float>, EmitterDevice2dGPU<float>,
GridData2dpolarGPU<float>, float, BoundaryContainer2dGPU<float>> device2dpolarfloatGPU; using
DeviceStatusGPU<Particles2dpolarGPU<double>, EmitterDevice2dGPU<double>, GridData2dpolarGPU<double>,
double, BoundaryContainer2dGPU<double>> device2dpolardoubleGPU;*/

// using DeviceStatus<Particles3d<double>, EmitterDevice3d<double>,
// const std::shared_ptr<GridData3d<double>, double, BoundaryContainer3d<double>,
// MeshContainer3d<double>>
// RFQDesignModel;

/*using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> Solver2daxsfloat; using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> Solver2daxsfloat; using Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> Solver2daxsdouble; using Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
Solver2daxsdouble;


using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> Solver2dfloat; using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> Solver2dfloat; using Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> Solver2ddouble; using Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
Solver2ddouble;

using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> Solver2dpolarfloat; using Solver <ParticleGridInterface<float>,
EmissionCurrentSolver<float>, PoissonSolver<float>, float> Solver2dpolarfloat; using Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
Solver2dpolardouble; using Solver <ParticleGridInterface<double>, EmissionCurrentSolver<double>,
PoissonSolver<double>, double> Solver2dpolardouble;



using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> SolverPTI2dfloat; using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> SolverPIC2daxsfloat; using Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> SolverPTI2ddouble; using Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
SolverPIC2daxsdouble;

using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> SolverPIC2dfloat; using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>,
PoissonSolver<float>, float> SolverPTI2daxsfloat; using Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> SolverPIC2ddouble; using Solver
<ParticleGridInterface<double>, EmissionCurrentSolver<double>, PoissonSolver<double>, double>
SolverPTI2daxsdouble;


using Solver <ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>,
float> SolverPIC2dpolarfloat; using Solver <ParticleGridInterface<double>,
EmissionCurrentSolver<double>, PoissonSolver<double>, double> SolverPIC2dpolardouble; using Solver
<ParticleGridInterface<float>, EmissionCurrentSolver<float>, PoissonSolver<float>, float>
SolverPTI2dpolarfloat; using Solver <ParticleGridInterface<double>, EmissionCurrentSolver<double>,
PoissonSolver<double>, double> SolverPTI2dpolardouble;*/

#endif