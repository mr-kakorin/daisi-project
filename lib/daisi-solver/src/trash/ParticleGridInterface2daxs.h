#ifndef PARTICLEGRIDINTERFACE2DAXS_H
#define PARTICLEGRIDINTERFACE2DAXS_H
#include "MeshContainer2d.h"
#include "ParticleShape2d.h"
template <class PointType> class MeshContainer2d;
template <class PointType> class ParticleGridInterface2daxs
{
    friend class boost::serialization::access;
    friend MeshContainer2d<PointType>;

  private:
    std::vector<IParticleShape2d<PointType>> cellArray;
    void AddCell(unsigned short particleShapeType, unsigned int basePoint, unsigned int levelLow,
                 unsigned int levelHigh)
    {
        IParticleShape2d* newcell = new ParticleShapeCIC2d(basePoint, levelLow, levelHigh);
        switch (particleShapeType)
        {
        case 0:
            break;
        case 1:
            newcell = new ParticleShapeTSC2d(basePoint, levelLow, levelHigh);
            break;
        }
        cellArray.push_back(newcell);
    };

  public:
    void Particle2Grid();
};
#endif