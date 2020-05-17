#include "Particle.h"
#include "Dmath.h"

#include <Constants.h>

#ifdef USE_BLAS
#include <cblas.h>
#else
#include <mkl.h>
#endif

template class Particles2d<double>;
template class Particles2d<float>;
template class ParticlesBase<double>;
template class ParticlesBase<float>;

template <class PointType>
void ParticlesBase<PointType>::GetEmittanceData(std::vector<std::vector<float>>& data, int emFlag,
                                                double mass, double lambda)
{
    if (!this->positions.size())
        return;
    if (!this->positions[0].size())
        return;
    data.resize(2);
    data[0].resize(this->positions[0].size());
    data[1].resize(this->positions[0].size());
    PointType en = LIGHT_VELOCITY() * LIGHT_VELOCITY() * mass;

    double enAv  = 0;
    double zAv   = 0;
    double omega = 2 * PI() / lambda;

    float tmp;
    float beta2;

    switch (emFlag)
    {
    case 0:

        for (int i = 0; i < q.size(); i++)
        {
            double gamma = sqrt(1 + this->momentums[0][i] * this->momentums[0][i] +
                                this->momentums[1][i] * this->momentums[1][i] +
                                this->momentums[2][i] * this->momentums[2][i]);

            /*if (gamma < 1.002)
            {
                    double v = LIGHT_VELOCITY()*sqrt(this->momentums[0][i] *
            this->momentums[0][i] +
            this->momentums[1][i] *
            this->momentums[1][i] + this->momentums[2][i] * this->momentums[2][i]); data[1][i] =
            mass*v*v / 2;

            }
            else
                    data[1][i] = -(gamma - 1)*en;*/

            data[1][i] = this->momentums[2][i];

            beta2 = this->momentums[2][i] / gamma;

            // enAv = enAv + data[1][i];
            enAv       = enAv + this->momentums[2][i];
            data[0][i] = this->positions[2][i];
            /*data[0][i] = 180 * this->positions[2][i] * omega / PI();
            while (data[0][i] > 180)
                    data[0][i] = data[0][i] - 360;*/
            zAv = zAv + data[0][i];
        };

        if (q.size() != 0)
        {
            zAv  = zAv / q.size();
            enAv = enAv / q.size();
        }
        for (int i = 0; i < q.size(); i++)
        {
            data[1][i] = 100 * (data[1][i] - enAv) / enAv;
            data[0][i] = data[0][i];
        };
        break;
    case 1:
        for (int i = 0; i < q.size(); i++)
        {
            double gamma = sqrt(1 + this->momentums[0][i] * this->momentums[0][i] +
                                this->momentums[1][i] * this->momentums[1][i] +
                                this->momentums[2][i] * this->momentums[2][i]);
            beta2      = this->momentums[2][i] / gamma;
            data[0][i] = this->positions[0][i];
            data[1][i] = this->momentums[0][i] / (beta2);
        };
        break;
    case 2:
        for (int i = 0; i < q.size(); i++)
        {
            double gamma = sqrt(1 + this->momentums[0][i] * this->momentums[0][i] +
                                this->momentums[1][i] * this->momentums[1][i] +
                                this->momentums[2][i] * this->momentums[2][i]);
            beta2      = this->momentums[2][i] / gamma;
            data[0][i] = this->positions[1][i];
            data[1][i] = this->momentums[1][i] / (beta2);
        };
        break;
    case 3:
        for (int i = 0; i < q.size(); i++)
        {
            data[0][i] = this->positions[0][i];
            data[1][i] = this->positions[1][i];
        }
        break;
    };
};

template <class PointType>
void ParticlesBase<PointType>::RecombinateParticles(ParticlesBase<PointType>* newParticles,
                                                    int                       nParticles)
{
    std::vector<unsigned int> Indexes;
    for (int i = 0; i < nParticles; i++)
        Indexes.push_back(newParticles->q.size() - nParticles + i);

    InsertParticels(newParticles, Indexes);

    newParticles->removeParticle(Indexes);
};

template <class PointType>
PointType ParticlesBase<PointType>::GetTotalCurrent()
{
    double result = 0;
    for (int i = 0; i < q.size(); i++)
        result = result + q[i];
    return result;
};

template <class PointType>
long long ParticlesBase<PointType>::GetMemorySize()
{
    volatile long long result = 0;
    result = result + sizeof(cellsNumbers) + sizeof(cellsNumbers[0]) * cellsNumbers.capacity();
    result = result + sizeof(flagEmitted) + sizeof(flagEmitted[0]) * flagEmitted.capacity();
    result = result + sizeof(Get_currentCoef()) +
             sizeof(Get_currentCoef()[0]) * Get_currentCoef().capacity();

    return result;
};

template <class PointType>
void ParticlesBase<PointType>::markBadParticles(const std::vector<unsigned int>& indexes)
{
    for (int i                   = 0; i < indexes.size(); i++)
        cellsNumbers[indexes[i]] = -1;
};

template <class PointType>
void ParticlesBase<PointType>::Reset()
{
#pragma ivdep
    for (int i          = 0; i < q.size(); i++)
        cellsNumbers[i] = 0;
};
template <class PointType>
ParticlesBase<PointType>::ParticlesBase()
{
    Time     = 0;
    avCharge = 1e27;
    minStep  = -1;
};

template <class PointType>
PointType* ParticlesBase<PointType>::GetPointerToMomentum1()
{
    return &this->momentums[0][0];
};
template <class PointType>
PointType* ParticlesBase<PointType>::GetPointerToMomentum2()
{
    return &this->momentums[1][0];
};
template <class PointType>
PointType* ParticlesBase<PointType>::GetPointerToMomentum3()
{
    return &this->momentums[2][0];
};

template <class PointType>
void ParticlesBase<PointType>::searchBadParticle(std::vector<unsigned int>& Indexes)
{
    for (int i = 0; i < q.size(); i++)
    {
        if (-1 == cellsNumbers[i])
        {
            int flagAdd = 0;
            for (int j = 0; j < Indexes.size(); j++)
            {
                if (Indexes[j] == i)
                {
                    flagAdd = 1;
                    break;
                };
            }
            if (0 == flagAdd)
                Indexes.push_back(i);
        }
    };
};

template <class PointType>
void ParticlesBase<PointType>::GetSaveIndexes(std::vector<unsigned int>& saveIndexes,
                                              double outTime, int step, int saveParam,
                                              double tracesSaveProbability, int& flag)
{
    saveIndexes.clear();
    int t;
    if (outTime >= 0)
    {
        if (Time >= outTime && flag == 0)
        {
            for (int i = 0; i < q.size(); i++)
            {
                if (cellsNumbers[i] == 0)
                {
                    flag = 1;
                    t    = rand() % 100;
                    if (t <= tracesSaveProbability * 100)
                        saveIndexes.push_back(i);
                }
            };
        }
    }
    else
    {
        flag = 1;
        for (int i = 0; i < q.size(); i++)
        {
            if (cellsNumbers[i] == 0)
            {
                t = rand() % 100;
                if (t <= tracesSaveProbability * 100)
                    saveIndexes.push_back(i);
            }
        };
    };
};

template <class PointType>
std::vector<unsigned int> ParticlesBase<PointType>::GetStartCellNumbersGrid()
{
    std::vector<unsigned int> result;

    int i0 = 0;
    int i1 = -1;

    while (1)
    {
        do
        {
            i1++;
            if (i1 == q.size() - 1)
                break;
        } while (additionalInf.back()[i1] == additionalInf.back()[i1 + 1]);

        if (i1 - i0 > 1)
            result.push_back(additionalInf.back()[i1]);

        if (i1 == q.size() - 1)
            break;

        i0 = i1 + 1;
    };
    return result;
};

template <class PointType>
void ParticlesBase<PointType>::removeBadParticle()
{
    std::vector<unsigned int> Indexes;
    searchBadParticle(Indexes);
    removeParticle(Indexes);
};

template <class PointType>
PointType ParticlesBase<PointType>::GetBeta(int i)
{
    PointType pp = 0;
    for (int j = 0; j < this->momentums.size(); j++)
        pp = pp + this->momentums[j][i] * this->momentums[j][i];

    PointType gamma = sqrt(1 + pp);
    PointType beta  = sqrt(gamma * gamma - 1) / gamma;
    return beta;
}
template <class PointType>
void ParticlesBase<PointType>::GammaCalc(std::vector<PointType>& gamma)
{
    gamma.resize(this->momentums[0].size());

    for (int k = 0; k < this->momentums[0].size(); k++)
    {

        PointType pp = 0;
        for (int j = 0; j < this->momentums.size(); j++)
            pp = pp + this->momentums[j][k] * this->momentums[j][k];

        gamma[k] = sqrt(1 + pp);
    }
};

template <class PointType>
void ParticlesBase<PointType>::GetBetaComponents(PointType& beta1, PointType& beta2,
                                                 PointType& beta3, int i)
{
    PointType pp = 0;
    for (int j = 0; j < this->momentums.size(); j++)
        pp = pp + this->momentums[j][i] * this->momentums[j][i];

    PointType gamma = sqrt(1 + pp);
    beta1           = this->momentums[0][i] / gamma;
    beta2           = this->momentums[1][i] / gamma;
    if (this->momentums.size() == 3)
        beta3 = this->momentums[2][i] / gamma;
    else
        beta3 = 0;
};

template <class PointType>
void ParticlesBase<PointType>::resize(int totalParticles)
{
    for (int j = 0; j < this->momentums.size(); j++)
        this->momentums[j].resize(totalParticles);
    for (int j = 0; j < this->positions.size(); j++)
        this->positions[j].resize(totalParticles);

    for (int j = 0; j < additionalInf.size(); j++)
        additionalInf[j].resize(totalParticles);

    for (int j = 0; j < this->additionalInfType.size(); j++)
        this->additionalInfType[j].resize(totalParticles);

    for (int j = 0; j < this->coloumbFields.size(); j++)
        this->coloumbFields[j].resize(totalParticles);

    for (int j = 0; j < magneticFields.size(); j++)
        magneticFields[j].resize(totalParticles);

    q.resize(totalParticles);
    cellsNumbers.resize(totalParticles);
    flagEmitted.resize(totalParticles);
};

template <class PointType>
PointType ParticlesBase<PointType>::GetEnergy(int number, PointType mass)
{
    return 0;
};

template <class PointType>
void ParticlesBase<PointType>::ReserveMemory(int size)
{
    for (int j = 0; j < this->momentums.size(); j++)
        this->momentums[j].reserve(size);
    for (int j = 0; j < this->positions.size(); j++)
        this->positions[j].reserve(size);

    for (int j = 0; j < additionalInf.size(); j++)
        additionalInf[j].reserve(size);

    for (int j = 0; j < this->additionalInfType.size(); j++)
        this->additionalInfType[j].reserve(size);

    for (int j = 0; j < this->coloumbFields.size(); j++)
        this->coloumbFields[j].reserve(size);

    for (int j = 0; j < magneticFields.size(); j++)
        magneticFields[j].reserve(size);

    q.reserve(size);
    cellsNumbers.reserve(size);
    flagEmitted.reserve(size);
};

template <class PointType>
void ParticlesBase<PointType>::resizeLight(int totalParticles)
{
    for (int j = 0; j < this->positions.size(); j++)
        this->positions[j].resize(totalParticles);
};

template <class PointType>
void ParticlesBase<PointType>::FastCopy(ParticlesBase<PointType>* object, int i1, int i2)
{
    if (q.size() == 0)
        return;
    object->resizeLight(q.size());
    int size = sizeof(PointType);
    int ink  = 1;

    object->cellsNumbers.resize(i2 - i1);

    for (int i                       = i1; i < i2; i++)
        object->cellsNumbers[i - i1] = cellsNumbers[i];

    switch (size)
    {
    case 4:
        for (int k = 0; k < this->positions.size(); k++)
            cblas_scopy(i2 - i1, (float*)&this->positions[k][i1], ink,
                        (float*)&object->positions[k][0], ink);

    case 8:
        for (int k = 0; k < this->positions.size(); k++)
            cblas_dcopy(i2 - i1, (double*)&this->positions[k][i1], ink,
                        (double*)&object->positions[k][0], ink);
    };
};

template <class PointType>
void ParticlesBase<PointType>::removeParticle(std::vector<unsigned int>& indexes)
{
    //	searchBadParticle(indexes);

    int n = int(indexes.size());
    for (int i = 0; i < n; i++)
    {
        this->positions[0][indexes[i]] = -1;
    };

    int n1 = int(this->positions[0].size());
    int j;
    for (int i = 0; i < n1; i++)
    {
        if (this->positions[0][i] == -1)
        {
            j = i;
            while (this->positions[0][j] == -1)
            {
                j++;
                if (j >= n1)
                    break;
            }

            for (int k = 0; k < this->momentums.size(); k++)
                this->momentums[k].erase(this->momentums[k].begin() + i,
                                         this->momentums[k].begin() + j);

            for (int k = 0; k < this->positions.size(); k++)
                this->positions[k].erase(this->positions[k].begin() + i,
                                         this->positions[k].begin() + j);

            //	for (int k = 0; k < this->momentums.size(); k++)
            //		fields[k].erase(fields[k].begin() + i, fields[k].begin() + j);

            for (int k = 0; k < this->additionalInfType.size(); k++)
                this->additionalInfType[k].erase(this->additionalInfType[k].begin() + i,
                                                 this->additionalInfType[k].begin() + j);

            for (int k = 0; k < additionalInf.size(); k++)
                additionalInf[k].erase(additionalInf[k].begin() + i, additionalInf[k].begin() + j);

            for (int k = 0; k < this->coloumbFields.size(); k++)
                this->coloumbFields[k].erase(this->coloumbFields[k].begin() + i,
                                             this->coloumbFields[k].begin() + j);

            for (int k = 0; k < magneticFields.size(); k++)
                magneticFields[k].erase(magneticFields[k].begin() + i,
                                        magneticFields[k].begin() + j);

            q.erase(q.begin() + i, q.begin() + j);
            flagEmitted.erase(flagEmitted.begin() + i, flagEmitted.begin() + j);
            cellsNumbers.erase(cellsNumbers.begin() + i, cellsNumbers.begin() + j);

            n1 = n1 - (j - i);
            i  = i - 1;
        };
    }

    //	indexes.clear();
};

template <class PointType>
std::vector<float> ParticlesBase<PointType>::GetData(int i, float mass, float charge)
{
    std::vector<float> result(2 * this->positions.size() + 1);
    if (this->positions[0].size() == 0)
        return result;

    for (int k    = 0; k < this->positions.size(); k++)
        result[i] = this->positions[k][i];

    for (int k                             = 0; k < this->momentums.size(); k++)
        result[i + this->positions.size()] = this->momentums[k][i];

    result.back() = q[i];

    return result;
};
template <class PointType>
std::vector<void*> ParticlesBase<PointType>::GetData()
{
    std::vector<void*> result(2 * this->positions.size() + 1);
    if (this->positions.size() == 0)
        return result;

    if (this->positions[0].size() == 0)
        return result;

    for (int k    = 0; k < this->positions.size(); k++)
        result[k] = (void*)(&this->positions[k][0]);

    for (int k                             = 0; k < this->momentums.size(); k++)
        result[k + this->positions.size()] = (void*)(&this->momentums[k][0]);

    result.back() = (void*)(&q[0]);
    return result;
};
template <class PointType>
int ParticlesBase<PointType>::SpaseSize()
{
    return 2 * this->positions.size() + 1;
};
template <class PointType>
Particles2d<PointType>::Particles2d(int currentSolverType)
{
    this->momentums.resize(3);
    this->positions.resize(3);
    this->coloumbFields.resize(2);
    if (currentSolverType == 1)
    {
        this->additionalInfType.resize(1);
        this->additionalInf.resize(1);
    }
};

template <class PointType>
ParticlesBase<PointType>& ParticlesBase<PointType>::operator=(ParticlesBase<PointType>* object)
{
    if (this == object)
    {
        return *this;
    }
    this->positions = object->positions;
    this->momentums = object->momentums;
    return *this;
};
template <class PointType>
void ParticlesBase<PointType>::GetParticlesCloud(int flag, std::vector<void*>& pointArray,
                                                 int& sizeArray, int& sizeElement)
{
    if (q.size() <= 0)
    {
        sizeArray   = 0;
        sizeElement = 8;
        return;
    }
    pointArray.resize(this->positions.size());
    sizeArray = q.size();
    for (int k        = 0; k < this->positions.size(); k++)
        pointArray[k] = (void*)(&this->positions[k][0]);

    sizeElement = sizeof(this->positions[0][0]);
};
template <class PointType>
void ParticlesBase<PointType>::clear()
{
    for (int j = 0; j < this->momentums.size(); j++)
        this->momentums[j].clear();
    for (int j = 0; j < this->positions.size(); j++)
        this->positions[j].clear();

    for (int j = 0; j < this->coloumbFields.size(); j++)
        this->coloumbFields[j].clear();
    for (int j = 0; j < magneticFields.size(); j++)
        magneticFields[j].clear();

    q.clear();
    cellsNumbers.clear();
    flagEmitted.clear();
    Time    = 0;
    minStep = -1;
};

template <class PointType>
void ParticlesBase<PointType>::GetAveragePositions(std::vector<PointType>& data)
{
    data.resize(3);
    /*for (int i = 0; i < NParticles; i++)
    {
            data[0] = data[0] + z[i];
            data[1] = data[1] + pz[i];
    };
    data[0] = data[0] / NParticles;
    data[1] = data[1] / NParticles;*/
}

template <class PointType>
void ParticlesBase<PointType>::InsertParticels(ParticlesBase<PointType>* prev)
{
    int nOld = q.size();
    resize(q.size() + prev->NParticles());
    for (int i = 0; i < prev->NParticles(); i++)
    {

        for (int j                       = 0; j < this->momentums.size(); j++)
            this->momentums[j][i + nOld] = prev->momentums[j][i];
        for (int j                       = 0; j < this->positions.size(); j++)
            this->positions[j][i + nOld] = prev->positions[j][i];

        for (int j                     = 0; j < additionalInf.size(); j++)
            additionalInf[j][i + nOld] = prev->additionalInf[j][i];

        for (int j                               = 0; j < this->additionalInfType.size(); j++)
            this->additionalInfType[j][i + nOld] = prev->additionalInfType[j][i];

        for (int j                           = 0; j < this->coloumbFields.size(); j++)
            this->coloumbFields[j][i + nOld] = prev->coloumbFields[j][i];

        for (int j                      = 0; j < magneticFields.size(); j++)
            magneticFields[j][i + nOld] = prev->magneticFields[j][i];

        q[i + nOld]            = prev->q[i];
        cellsNumbers[i + nOld] = prev->cellsNumbers[i];
        flagEmitted[i + nOld]  = prev->flagEmitted[i];
    };
}
template <class PointType>
void ParticlesBase<PointType>::InsertParticels(ParticlesBase*                   prev,
                                               const std::vector<unsigned int>& indexes)
{
    int nOld = q.size();
    resize(q.size() + indexes.size());
    for (int i = 0; i < indexes.size(); i++)
    {

        for (int j                       = 0; j < this->momentums.size(); j++)
            this->momentums[j][i + nOld] = prev->momentums[j][indexes[i]];
        for (int j                       = 0; j < this->positions.size(); j++)
            this->positions[j][i + nOld] = prev->positions[j][indexes[i]];

        for (int j                     = 0; j < additionalInf.size(); j++)
            additionalInf[j][i + nOld] = prev->additionalInf[j][indexes[i]];

        for (int j                               = 0; j < this->additionalInfType.size(); j++)
            this->additionalInfType[j][i + nOld] = prev->additionalInfType[j][indexes[i]];

        for (int j                           = 0; j < this->coloumbFields.size(); j++)
            this->coloumbFields[j][i + nOld] = prev->coloumbFields[j][indexes[i]];

        for (int j                      = 0; j < magneticFields.size(); j++)
            magneticFields[j][i + nOld] = prev->magneticFields[j][indexes[i]];

        q[i + nOld]            = prev->q[indexes[i]];
        cellsNumbers[i + nOld] = prev->cellsNumbers[indexes[i]];
        flagEmitted[i + nOld]  = prev->flagEmitted[indexes[i]];
    };
};
template <class PointType>
void ParticlesBase<PointType>::InsertParticelsEmittances(ParticlesBase*                   prev,
                                                         const std::vector<unsigned int>& indexes,
                                                         double                           zEm)
{
    int nOld = q.size();
    resize(q.size() + indexes.size());
    for (int i = 0; i < indexes.size(); i++)
    {

        for (int j                       = 0; j < this->momentums.size(); j++)
            this->momentums[j][i + nOld] = prev->momentums[j][indexes[i]];

        for (int j                       = 0; j < this->positions.size() - 1; j++)
            this->positions[j][i + nOld] = prev->positions[j][indexes[i]];

        double gamma = sqrt(1 + this->momentums[0][i] * this->momentums[0][i] +
                            this->momentums[1][i] * this->momentums[1][i] +
                            this->momentums[2][i] * this->momentums[2][i]);
        double beta2 = this->momentums[2][i] / gamma;

        this->positions[2][i + nOld] = prev->Time + (zEm - prev->positions[2][indexes[i]]) / beta2;

        this->positions[2][i + nOld] =
            1e9 * this->positions[2][i + nOld] / LIGHT_VELOCITY();
    };
};

template <class PointType>
std::vector<int> ParticlesBase<PointType>::CheckEmittanceCondition(double                    cond,
                                                                   ParticlesBase<PointType>* prev,
                                                                   int i1, int i2)
{
    std::vector<int> result;
    for (int i = i1; i < i2; i++)
    {
        if ((this->positions[2][i] > cond && prev->positions[2][i - i1] < cond) ||
            (this->positions[2][i] < cond && prev->positions[2][i - i1] > cond))
            result.push_back(i);
    };
    return result;
};

template <class PointType>
void ParticlesBase<PointType>::CopyParticle(ParticlesBase<PointType>* object, int index)
{
    for (int i                      = 0; i < object->positions.size(); i++)
        object->positions[i][index] = this->positions[i][index];
};

template <class PointType>
void ParticlesBase<PointType>::MultiplyMomentum(int i, double kin)
{
    for (int k                = 0; k < this->momentums.size(); k++)
        this->momentums[k][i] = this->momentums[k][i] * kin;
};
template <class PointType>
void ParticlesBase<PointType>::SetCartesianMomentum(int i, const std::vector<PointType>& PX)
{
    for (int k                = 0; k < PX.size(); k++)
        this->momentums[k][i] = PX[k];
};

template <class PointType>
void ParticlesBase<PointType>::SetCartesianPosition(int i, const std::vector<PointType>& xin)
{
    for (int k                = 0; k < xin.size(); k++)
        this->positions[k][i] = xin[k];
};

template <class PointType>
void ParticlesBase<PointType>::setNewParticles(const std::vector<unsigned int>& EmptyPlaces,
                                               ParticlesBase*                   newParticles){
    /*int empty = int(EmptyPlaces.size());
    int totalParticles = newParticles->q.size();
    int nowParticles = NParticles;
    if (empty<totalParticles)
            resize(totalParticles + NParticles - empty);

    int index;
    for (int k = 0; k <totalParticles; k++)
    {
            if (k < empty)
                    index = EmptyPlaces[k];
            else
            {
                    index = nowParticles + k - empty;
            }


            r[index] = newParticles->GetPointerToPosition1()[k];
            z[index] = newParticles->GetPointerToPosition2()[k];
            pr[index] = newParticles->Get_pr()[k];
            pz[index] = newparticles->GetPointerToPosition2()[k];
            phi[index] = newParticles->Get_phi()[k];
            pphi[index] = newParticles->Get_pphi()[k];
            q[index] = newParticles->q[k];
            cellsNumbers[index] = 0;
    };

    std::vector<unsigned int> remove;
    for (int i = totalParticles; i < EmptyPlaces.size(); i++)
            remove.push_back(EmptyPlaces[i]);

    removeParticle(remove);*/
};