#ifndef PARTICLEGRIDINTERFACE2DGPU_H
#define PARTICLEGRIDINTERFACE2DGPU_H
#include "GPUTypes.h"
#include "ParticleGridInterface2d.h"
template <class PointType> class ParticleGridInterface2d;
class ParticleGridInterface2dGPU
{
    template <class PointType> friend ParticleGridInterface2d<PointType>;

  private:
    int            priorityParticleShapeType;
    GPUTypes::imat templNumb;
    int*           CICArray;

  public:
    template <class PointType> ParticleGridInterface2dGPU(ParticleGridInterface2d<PointType>& obj)
    {
        priorityParticleShapeType = obj.priorityParticleShapeType;
        templNumb                 = obj.templNumb;
        int s                     = obj.CICArray.size();
    };
};

#endif