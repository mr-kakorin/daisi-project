#ifndef PARTICLEGRIDINTERFACE1D_H
#define PARTICLEGRIDINTERFACE1D_H
#include "ParticleShape1d.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <vector>
class ParticleGridInterface1d
{
    friend class boost::serialization::access;
    std::vector<IParticleShape1d*> cellArray;

  public:
    void GenerateCells(int priorityType)
    {
        ParticleShapeCIC1d* cell1 = new ParticleShapeCIC1d();
        cellArray.push_back(cell1);
    };
    void Particle2Grid();
};
#endif