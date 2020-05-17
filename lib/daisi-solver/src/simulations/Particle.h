#ifndef PARTICLES_H
#define PARTICLES_H

#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <Constants.h>
#include <cmath>

template <class PointType>
class particlesFields
{
    std::vector<std::vector<PointType>> fields;

  public:
    void Init(int nflows, int fieldSize, int blockSize)
    {
        fields.resize(fieldSize);
        for (int j = 0; j < fieldSize; j++)
            fields[j].resize(blockSize);
    };
    std::vector<PointType>& Get_Er()
    {
        return fields[0];
    };
    std::vector<PointType>& Get_Ephi()
    {
        return fields[1];
    };
    std::vector<PointType>& Get_Bphi()
    {
        return fields[3];
    };
    std::vector<PointType>& Get_Ex()
    {
        return fields[0];
    };
    std::vector<PointType>& Get_Ey()
    {
        return fields[1];
    };
    std::vector<PointType>& Get_Ez()
    {
        return fields[2];
    };

    std::vector<PointType>& Get_Bx()
    {
        return fields[3];
    };
    std::vector<PointType>& Get_By()
    {
        return fields[4];
    };
    std::vector<PointType>& Get_Bz()
    {
        return fields[5];
    };
};

template <class PointType>
class ParticlesBase
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& Time;
        ar& flagEmitted;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& Time;
        ar& flagEmitted;
    }

  public:
    int NParticles()
    {
        return q.size();
    };
    int  currentSolverType;
    void GetAveragePositions(std::vector<PointType>& data);
    void clear();
    void removeBadParticle();
    void CopyParticle(ParticlesBase<PointType>* object, int index);
    void InsertParticels(ParticlesBase* data);
    void InsertParticels(ParticlesBase* prev, const std::vector<unsigned int>& indexes);
    virtual PointType GetCartesianPX(int i)
    {
        return this->momentums[0][i];
    };
    virtual PointType GetCartesianPY(int i)
    {
        return this->momentums[1][i];
    };
    PointType* GetPointerToPosition1()
    {
        return &positions[0][0];
    };
    PointType* GetPointerToPosition2()
    {
        return &positions[1][0];
    };
    PointType* GetPointerToPosition3()
    {
        return &positions[1][0];
    };
    virtual PointType GetBeta(int i);
    virtual void GetBetaComponents(PointType& beta1, PointType& beta2, PointType& beta3, int i);
    virtual std::vector<int> CheckEmittanceCondition(double cond, ParticlesBase<PointType>* prev,
                                                     int i1, int i2);

    void resize(int totalParticles);
    void ReserveMemory(int size);
    void resizeLight(int totalParticles);
    void FastCopy(ParticlesBase* object, int i1, int i2);
    void removeParticle(std::vector<unsigned int>& indexes);
    std::vector<float> GetData(int i, float mass, float charge);
    std::vector<void*> GetData();
    int                SpaseSize();
    ParticlesBase& operator=(ParticlesBase* object);
    void GetParticlesCloud(int flag, std::vector<void*>& pointArray, int& sizeArray,
                           int& sizeElement);

    std::vector<std::vector<PointType>> positions;
    std::vector<std::vector<PointType>> momentums;
    std::vector<std::vector<PointType>> coloumbFields;
    std::vector<std::vector<PointType>> magneticFields;

    void searchBadParticle(std::vector<unsigned int>& Indexes);
    void GetSaveIndexes(std::vector<unsigned int>& saveIndexes, double outTime, int step,
                        int saveParam, double tracesSaveProbability, int& flag);
    double minStep;
    double Time;
    // std::vector<char> Get_isPeriodical();

    // std::vector<unsigned int> startCellNumbers;
    // std::vector<PointType> gamma;
    std::vector<PointType>     q;
    std::vector<unsigned int>  cellsNumbers;
    std::vector<unsigned char> flagEmitted;
    // std::vector<unsigned char> Get_isPeriodical();
    PointType avCharge;
    //	std::vector <std::vector< PointType>> Get_currentCoef();
    std::vector<std::vector<PointType>> additionalInfType; // Get_currentCoef() - 0
    //	std::vector <unsigned int> Get_startCellNumbers();
    std::vector<std::vector<unsigned int>>
        additionalInf; // Get_isPeriodical() - 1 Get_startCellNumbers()-0
  public:
    void InsertParticelsEmittances(ParticlesBase* prev, const std::vector<unsigned int>& indexes,
                                   double zEm);
    std::vector<unsigned int>& Get_startCellNumbers()
    {
        return additionalInf.back();
    };
    std::vector<unsigned int>& Get_isPeriodical()
    {
        return additionalInf[0];
    };
    /*	std::vector <unsigned int>& Get_cellsNumbersAdd()
            {
                    return additionalInf[1];
            };*/
    std::vector<PointType>& Get_currentCoef()
    {
        return additionalInfType[0];
    };
    virtual void GammaCalc(std::vector<PointType>& gamma);
    PointType                 GetTotalCurrent();
    std::vector<unsigned int> GetStartCellNumbersGrid();

    long long GetMemorySize();

    void markBadParticles(const std::vector<unsigned int>& indexes);

    void Reset();
    ParticlesBase();

    PointType* GetPointerToMomentum1();
    PointType* GetPointerToMomentum2();
    PointType* GetPointerToMomentum3();

    virtual PointType* GetPointerToCartesianX()
    {
        return &positions[0][0];
    };
    virtual PointType* GetPointerToCartesianY()
    {
        return &positions[1][0];
    };
    virtual void PeriodicalEvent(){

    };
    virtual void GetBeamMeasuriments(std::vector<std::vector<PointType>>& data, double mass,
                                     double Zav, double betaAv){

    };

    virtual void GetEmittanceData(std::vector<std::vector<float>>& data, int emFlag, double mass,
                                  double lambda);
    virtual PointType GetEnergy(int number, PointType mass);
    virtual PointType GetR(int i)
    {
        return 0;
    };
    void MultiplyMomentum(int i, double k);
    void SetCartesianMomentum(int i, const std::vector<PointType>& PX);
    void SetCartesianPosition(int i, const std::vector<PointType>& xin);
    void setNewParticles(const std::vector<unsigned int>& EmptyPlaces, ParticlesBase* newParticles);
    void RecombinateParticles(ParticlesBase* newParticles, int nParticles);
};

template <class PointType>
class Particles1d : public ParticlesBase<PointType>
{
    friend class boost::serialization::access;

  private:
    std::vector<PointType> x, px;

  public:
    void GetParticlesCloud(int flag, void* pointArray[3], int& sizeArray, int& sizeElement){};
    int GetParticles(int a)
    {
        return a;
    };
};

template <class PointType>
class Particles2d : public ParticlesBase<PointType>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
        // ar & x;
        // ar & px;
        // ar & y;
        // ar & py;
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
    }

  public:
    std::vector<PointType>& Get_ExCol()
    {
        return this->coloumbFields[0];
    };
    std::vector<PointType>& Get_EyCol()
    {
        return this->coloumbFields[1];
    };

    int getFieldSize()
    {
        return 2;
    };
    std::vector<PointType>& Get_x()
    {
        return this->positions[0];
    };
    std::vector<PointType>& Get_y()
    {
        return this->positions[1];
    };
    std::vector<PointType>& Get_z()
    {
        return this->positions[2];
    };
    std::vector<PointType>& Get_px()
    {
        return this->momentums[0];
    };
    std::vector<PointType>& Get_py()
    {
        return this->momentums[1];
    };
    std::vector<PointType>& Get_pz()
    {
        return this->momentums[2];
    };
    Particles2d(int currentSolverType);
    Particles2d(){

    };
};
template <class PointType>
class Particles3dcil : public ParticlesBase<PointType>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
    }

  public:
    Particles3dcil(int currentSolverType);
    Particles3dcil(){

    };
    int getFieldSize()
    {
        return 4;
    };

    std::vector<PointType>& Bphi()
    {
        return this->magneticFields[0];
    };
    std::vector<PointType>& Get_r()
    {
        return this->positions[0];
    };
    std::vector<PointType>& Get_z()
    {
        return this->positions[1];
    };
    std::vector<PointType>& Get_ErCol()
    {
        return this->coloumbFields[0];
    };
    std::vector<PointType>& Get_EzCol()
    {
        return this->coloumbFields[1];
    };

    std::vector<PointType>& Get_phi()
    {
        return this->positions[2];
    };
    std::vector<PointType>& Get_pr()
    {
        return this->momentums[0];
    };
    std::vector<PointType>& Get_pz()
    {
        return this->momentums[1];
    };
    std::vector<PointType>& Get_pphi()
    {
        return this->momentums[2];
    };
    void GammaCalc(std::vector<PointType>& gamma);
    void GetBeamMeasuriments(std::vector<std::vector<PointType>>& data, double mass, double Zav,
                             double betaAv);
    void GetEmittanceData(std::vector<std::vector<float>>& data, int emFlag, double mass,
                          double lambda);
    PointType GetBeta(int i);
    void GetBetaComponents(PointType& beta1, PointType& beta2, PointType& beta3, int i);
    PointType GetEnergy(int number, PointType mass);
    void SetCartesianMomentumPolar(int i, PointType PX, PointType PY)
    {
        this->momentums[0][i] = PX * std::cos(2 * PI() - this->positions[2][i]) +
                                PY * std::sin(2 * PI() - this->positions[2][i]);
        this->momentums[2][i] = PX * std::sin(2 * PI() - this->positions[2][i]) -
                                PY * std::cos(2 * PI() - this->positions[2][i]);
        this->momentums[2][i] = this->momentums[2][i] * this->positions[0][i];
    }
};

template <class PointType>
class Particles2dpolar : public ParticlesBase<PointType>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    std::vector<PointType>& Get_ErCol()
    {
        return this->coloumbFields[0];
    };
    std::vector<PointType>& Get_EphiCol()
    {
        return this->coloumbFields[1];
    };

    Particles2dpolar(int currentSolverType);
    Particles2dpolar(){

    };

    int getFieldSize()
    {
        return 3;
    };
    void GetBeamMeasuriments(std::vector<std::vector<PointType>>& data, double mass, double Zav,
                             double betaAv);
    void GetEmittanceData(std::vector<std::vector<float>>& data, int emFlag, double mass,
                          double lambda);
    void PeriodicalEvent();
    void GetBetaComponents(PointType& beta1, PointType& beta2, PointType& beta3, int i);
    PointType GetBeta(int i);
    PointType GetEnergy(int number, PointType mass);

    PointType GetCartesianPX(int i);
    PointType GetCartesianPY(int i);

    void SetCartesianMomentum(int i, const std::vector<PointType>& PX);
    void SetCartesianPosition(int i, const std::vector<PointType>& xin);
    PointType* GetPointerToCartesianX();
    PointType* GetPointerToCartesianY();
    void GetParticlesCloud(int flag, std::vector<void*>& pointArray, int& sizeArray,
                           int& sizeElement);
    void GammaCalc(std::vector<PointType>& gamma);

    //// 0 r 1 phi 2 x 3 y 4 phiReal 5 xReal 6 yReal,
    std::vector<PointType>& Get_phiReal()
    {
        return this->positions[4];
    };
    std::vector<PointType>& Get_r()
    {
        return this->positions[0];
    };
    std::vector<PointType>& Get_phi()
    {
        return this->positions[1];
    };
    std::vector<PointType>& Get_pr()
    {
        return this->momentums[0];
    };
    std::vector<PointType>& Get_pphi()
    {
        return this->momentums[1];
    };
    std::vector<PointType>& Get_x()
    {
        return this->positions[2];
    };
    std::vector<PointType>& Get_y()
    {
        return this->positions[3];
    };
    std::vector<PointType>& Get_xReal()
    {
        return this->positions[5];
    };
    std::vector<PointType>& Get_yReal()
    {
        return this->positions[6];
    };
};

template <class PointType>
class Particles3d : public ParticlesBase<PointType>
{
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const
    {
        ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
    }
    template <class Archive>
    void load(Archive& ar, const unsigned int)
    {
        ar& boost::serialization::base_object<ParticlesBase<PointType>>(*this);
    }

  public:
    int getFieldSize()
    {
        return 6;
    };
    std::vector<PointType>& Get_x()
    {
        return this->positions[0];
    };
    std::vector<PointType>& Get_y()
    {
        return this->positions[1];
    };
    std::vector<PointType>& Get_z()
    {
        return this->positions[2];
    };
    std::vector<PointType>& Get_px()
    {
        return this->momentums[0];
    };
    std::vector<PointType>& Get_py()
    {
        return this->momentums[1];
    };
    std::vector<PointType>& Get_pz()
    {
        return this->momentums[2];
    };
    Particles3d(int currentSolverType)
    {
        this->momentums.resize(3);
        this->positions.resize(3);

        this->coloumbFields.resize(3);
        this->magneticFields.resize(3);
        if (currentSolverType == 1)
        {
            this->additionalInfType.resize(1);
            this->additionalInf.resize(1);
        }
    };
    Particles3d(){};
};

#endif