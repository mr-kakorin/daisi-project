#ifndef PARTICLESHAPE1D_H
#define PARTICLESHAPE1D_H
class IParticleShape1d
{
    friend class boost::serialization::access;

  protected:
    unsigned int basePoint;
    virtual bool InCell(double x) = 0;
};
class ParticleShapeCIC1d : public IParticleShape1d
{
    friend class boost::serialization::access;
    bool InCell(double x1)
    {
        return true;
    };
};
#endif