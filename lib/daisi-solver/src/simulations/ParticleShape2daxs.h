#ifndef PARTICLESHAPE2DAXS_H
#define PARTICLESHAPE2DAXS_H
class IParticleShape2daxs
{
  public:
    unsigned int basePoint;
    unsigned int levelLow;
    unsigned int levelHigh;
    IParticleShape2daxs(){};
    IParticleShape2daxs(unsigned int basePointIn, unsigned int levelLowIn, unsigned int levelHighIn)
    {
        basePoint = basePointIn;
        levelLow  = levelLowIn;
        levelHigh = levelHighIn;
    };
    virtual bool InCell(double x1, double x2) = 0;

  public:
    IParticleShape2daxs& operator=(const IParticleShape2daxs& right)
    {
        if (this == &right)
        {
            return *this;
        }
        basePoint = right.basePoint;
        levelLow  = right.levelLow;
        levelHigh = right.levelHigh;
        return *this;
    }
};
class ParticleShapeCIC2daxs : public IParticleShape2daxs
{
    bool InCell(double x1, double x2)
    {
        return true;
    };

  public:
    ParticleShapeCIC2daxs(){};
    ParticleShapeCIC2daxs(unsigned int basePointIn, unsigned int levelLowIn,
                          unsigned int levelHighIn);
};
ParticleShapeCIC2daxs::ParticleShapeCIC2daxs(unsigned int basePointIn, unsigned int levelLowIn,
                                             unsigned int levelHighIn)
    : IParticleShape2d(basePointIn, levelLowIn, levelHighIn){};
class ParticleShapeTSC2daxs : public IParticleShape2daxs
{
    bool InCell(double x1, double x2)
    {
        return true;
    };

  public:
    ParticleShapeTSC2daxs(unsigned int basePointIn, unsigned int levelLowIn,
                          unsigned int levelHighIn);
};
ParticleShapeTSC2daxs::ParticleShapeTSC2daxs(unsigned int basePointIn, unsigned int levelLowIn,
                                             unsigned int levelHighIn)
    : IParticleShape2d(basePointIn, levelLowIn, levelHighIn){};
#endif