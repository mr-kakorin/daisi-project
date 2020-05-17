#include "ParticlesMover.h"
#include "Dmath.h"
#include "Particle.h"
#include "updatePositionsVector.h"

template class ParticlesMover<double>;
template class ParticlesMover<float>;

template void
ParticlesMover<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                             const unsigned int file_version);
template void ParticlesMover<double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void
ParticlesMover<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                              const unsigned int file_version);
template void
ParticlesMover<float>::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                             const unsigned int file_version) const;

template <class PointType>
template <class Archive>
void ParticlesMover<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& params;
}
template <class PointType>
template <class Archive>
void ParticlesMover<PointType>::load(Archive& ar, const unsigned int)
{
    ar& params;
    flagInit = 0;
    params.resize(2);
}

template <class PointType>
void ParticlesMover<PointType>::InitParallel(int numThreads)
{
    timeSteps.resize(numThreads);
    DT.resize(numThreads);
};

template <class PointType>
ParticlesMover<PointType>::ParticlesMover()
{
    params.resize(2);
    params[0] = 2; // timestep
    params[1] = 0; // type
};
template <class PointType>
PointType ParticlesMover<PointType>::GetTimeStep(int i, int thread)
{
    return DT[thread];
};
template <class PointType>
void ParticlesMover<PointType>::SetTimeSteps(std::vector<double> in)
{
    for (int i = 0; i < timeSteps.size(); i++)
    {
        timeSteps[i] = in;
    }
};
template <class PointType>
void ParticlesMover<PointType>::SetParameters(std::vector<double> in)
{
    params = in;
};
template <class PointType>
void ParticlesMover<PointType>::DoubleTimeSteps()
{
    //	for (int i = 0; i < timeSteps.size();i++)
    //		timeSteps[i] = timeSteps[i] * 2;
}
template <class PointType>
void ParticlesMover<PointType>::HalveTimeSteps()
{
    //	for (int i = 0; i < timeSteps.size();i++)
    //	timeSteps[i] = timeSteps[i] / 2;
}

template <class PointType>
std::vector<double> ParticlesMover<PointType>::GetParameters()
{
    return params;
};

template void ParticlesMover<float>::stepEstimate<Particles3dcil<float>>(
    const std::shared_ptr<Particles3dcil<float>>& particlesData, int nlow, int i1, int i2,
    int thread, int solverType);
template void ParticlesMover<float>::stepEstimate<Particles2d<float>>(
    const std::shared_ptr<Particles2d<float>>& particlesData, int nlow, int i1, int i2, int thread,
    int solverType);
template void ParticlesMover<float>::stepEstimate<Particles2dpolar<float>>(
    const std::shared_ptr<Particles2dpolar<float>>& particlesData, int nlow, int i1, int i2,
    int thread, int solverType);

template void ParticlesMover<double>::stepEstimate<Particles3dcil<double>>(
    const std::shared_ptr<Particles3dcil<double>>& particlesData, int nlow, int i1, int i2,
    int thread, int solverType);
template void ParticlesMover<double>::stepEstimate<Particles2d<double>>(
    const std::shared_ptr<Particles2d<double>>& particlesData, int nlow, int i1, int i2, int thread,
    int solverType);
template void ParticlesMover<double>::stepEstimate<Particles2dpolar<double>>(
    const std::shared_ptr<Particles2dpolar<double>>& particlesData, int nlow, int i1, int i2,
    int thread, int solverType);

template void ParticlesMover<double>::stepEstimate<Particles3d<double>>(
    const std::shared_ptr<Particles3d<double>>& particlesData, int nlow, int i1, int i2, int thread,
    int solverType);

template void ParticlesMover<float>::stepEstimate<Particles3d<float>>(
    const std::shared_ptr<Particles3d<float>>& particlesData, int nlow, int i1, int i2, int thread,
    int solverType);

template <class PointType>
template <class particlesdataType>
void ParticlesMover<PointType>::stepEstimate(
    const std::shared_ptr<particlesdataType>& particlesData, int nlow, int i1, int i2, int thread,
    int solverType)
{
    double    timeStep  = params[0];
    int       moverType = params[1];
    PointType dt;
    // if (flagInit == 0 || (flagInit == 1 && moverType == 1))
    if (flagInit == 0)
    {

        PointType timeStepLoc = particlesData->minStep / (timeStep);
        /*PointType cellSizeMin;
        for (int i = i1; i < i2; i++)
        {

                if (particlesData->Time< 2 * particlesData->cellSize[i] )
                        dt = particlesData->cellSize[i] / (timeStep * 4 );
                else
                        dt = particlesData->cellSize[i] / (timeStep);


                if (dt < timeStepLoc)
                {
                        timeStepLoc = dt;
                        cellSizeMin = particlesData->cellSize[i];
                }
        }*/
        /*	if (moverType == 0 || moverType == 2)
                        timeSteps[nlow] = timeStepLoc;

                if (moverType == 1)
                        timeSteps[nlow] = timeStepLoc / cellSizeMin;*/

        timeSteps[thread][nlow] = timeStepLoc;

        DT[thread] = timeStepLoc;
    }

    if (flagInit == 1)
    {

        if (moverType == 0 || moverType == 2)
        {
            DT[thread] = timeSteps[thread][nlow];

            /*if (solverType == 0)
            {
                    for (int i = 0; i < particlesData->NParticles(); i++)
                    {
                            if (particlesData->GetPointerToPosition1()[i]<0.11 &&
            particlesData->GetPointerToPosition1()[i]>0.09)
                            {
                                    DT[thread] = DT[thread] / 2;
                                    break;
                            }
                    }
            }*/
        }

        if (moverType == 1)
        {
            PointType timeStepLoc = particlesData->minStep / (timeStep);
            /*PointType cellSizeMin;
            for (int i = 0; i < particlesData->NParticles(); i++)
            {

                    if (particlesData->Time< 2 * particlesData->cellSize[i])
                            dt = particlesData->cellSize[i] / (timeStep * 4);
                    else
                            dt = particlesData->cellSize[i] / (timeStep);


                    if (dt < timeStepLoc)
                    {
                            timeStepLoc = dt;
                            cellSizeMin = particlesData->cellSize[i];
                    }
            }*/

            double rel = timeStepLoc / timeSteps[thread][nlow];

            std::vector<double> levels = {1, 10, 50, 100, 250, 500, 1e9};

            for (int i = 0; i < levels.size() - 1; i++)
            {

                if (levels[i] < rel && rel < levels[i + 1])
                    DT[thread] = timeSteps[thread][nlow] * levels[i];
            }
            /*PointType cellSizeMin = 1e9;

            for (int i = 0; i < particlesData->NParticles(); i++)
            {
                    if (cellSizeMin > particlesData->cellSize[i])
                    {
                            cellSizeMin = particlesData->cellSize[i];
                    }
            }
            DT = timeSteps[nlow] * cellSizeMin;*/
        };
    }
}

template <class PointType>
void ParticlesMover<PointType>::updatePositions(
    const std::shared_ptr<Particles3d<PointType>>& particlesData,
    const std::shared_ptr<Particles3d<PointType>>& particlesDataTmp, int nlow, int i1, int i2,
    int thread)
{
    PointType gamma;
    PointType pxn;
    PointType pyn;
    PointType pzn;
    PointType beta;

    for (int i = i1; i < i2; i++)
    {
        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            continue;

        pxn = particlesData->Get_px()[i];
        pyn = particlesData->Get_py()[i];
        pzn = particlesData->Get_pz()[i];

        gamma = sqrt(1 + pxn * pxn + pyn * pyn + pzn * pzn);

        particlesData->Get_x()[i] = particlesData->Get_x()[i] + DT[thread] * pxn / gamma;
        particlesData->Get_y()[i] = particlesData->Get_y()[i] + DT[thread] * pyn / gamma;
        particlesData->Get_z()[i] = particlesData->Get_z()[i] + DT[thread] * pzn / gamma;
    };
};
template <class PointType>
void ParticlesMover<PointType>::updateMomentums(
    const std::shared_ptr<Particles3d<PointType>>& particlesData,
    particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow, int i1, int i2,
    int thread)
{

    double K;

    if (params[2] == 0) // flagfriction
        K = 0;

    PointType B = params[3] / 10000.0; //������ � �����

    for (int i = i1; i < i2; i++)
    {

        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            DT[thread] = DT[thread] / 2;

        //		particlesDataTmp.Get_Ex()[i - i1] = 0;
        //		particlesDataTmp.Get_Ex()[i - i1] = 0;

        //	particlesData->gamma[i] = sqrt(1 + particlesData->Get_pz()[i] *
        // particlesData->Get_pz()[i] +
        // particlesData->Get_px()[i] * particlesData->Get_px()[i] + particlesData->Get_py()[i] *
        // particlesData->Get_py()[i]);

        //	particlesData->Get_px()[i] = particlesData->Get_px()[i] + DT[thread] *
        // alpha*particlesDataTmp.Get_Ex()[i
        //- i1]; 	particlesData->Get_py()[i] = particlesData->Get_py()[i] + DT[thread] *
        // alpha*particlesDataTmp.Get_Ey()[i - i1]; 	particlesData->Get_pz()[i] =
        // particlesData->Get_pz()[i] +
        // DT[thread] * alpha*particlesDataTmp.Get_Ez()[i - i1];

        PointType px12 = particlesData->Get_px()[i] +
                         0.5 * DT[thread] * alpha * particlesDataTmp.Get_Ex()[i - i1];
        PointType py12 = particlesData->Get_py()[i] +
                         0.5 * DT[thread] * alpha * particlesDataTmp.Get_Ey()[i - i1];
        PointType pz12 = particlesData->Get_pz()[i] +
                         0.5 * DT[thread] * alpha * particlesDataTmp.Get_Ez()[i - i1];

        PointType gamma12 = sqrt(1 + px12 * px12 + py12 * py12 + pz12 * pz12);

        PointType a = DT[thread] * alpha * particlesDataTmp.Get_Ex()[i - i1];

        PointType b = DT[thread] * alpha * particlesDataTmp.Get_Ey()[i - i1];

        PointType B1 = DT[thread] * B * alpha * LIGHT_VELOCITY() / (2 * gamma12);

        PointType B2 = -DT[thread] * B * alpha * LIGHT_VELOCITY() / (2 * gamma12);

        PointType pxold = particlesData->Get_px()[i];

        particlesData->Get_px()[i] =
            (particlesData->Get_px()[i] + a +
             (particlesData->Get_py()[i] + b + particlesData->Get_px()[i] * B2) * B1 +
             particlesData->Get_py()[i] * B1) /
            (1 - B1 * B2);

        particlesData->Get_py()[i] =
            particlesData->Get_py()[i] + b + (particlesData->Get_px()[i] + pxold) * B2;

        particlesData->Get_pz()[i] =
            particlesData->Get_pz()[i] + DT[thread] * alpha * particlesDataTmp.Get_Ez()[i - i1];

        if (particlesData->flagEmitted[i] < 1)
        {
            DT[thread] = DT[thread] * 2;
            continue;
        };
    };
};

template <class PointType>
void ParticlesMover<PointType>::updatePositions(
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesDataTmp, int nlow, int i1, int i2,
    int thread)
{
    PointType gamma;
    PointType prn;
    PointType pzn;
    PointType pphi;
    PointType r;
    PointType phiold;

    double dphi    = 0.058177641733144;
    double border1 = 1.5 * PI() - dphi / 2;
    double border2 = 1.5 * PI() + dphi / 2;

#pragma simd
    for (int i = i1; i < i2; i++)
    {
        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            continue;

        r = particlesData->Get_r()[i];

        prn                       = particlesData->Get_pr()[i];
        pphi                      = particlesData->Get_pphi()[i];
        gamma                     = sqrt(1 + prn * prn + (pphi / r) * (pphi / r));
        particlesData->Get_r()[i] = r + DT[thread] * prn / gamma;

        double dr = particlesData->Get_r()[i] - r;
        if (dr > 0)
        {
            int tt = 0;
        };

        PointType dPh                   = DT[thread] * pphi / (gamma * r * r);
        particlesData->Get_phiReal()[i] = particlesData->Get_phiReal()[i] + dPh;

        phiold                      = particlesData->Get_phi()[i];
        particlesData->Get_phi()[i] = particlesData->Get_phiReal()[i];

        particlesData->Get_isPeriodical()[i] = 0;

        if (particlesData->Get_phi()[i] < border1)
        {
            int c                       = (border1 - particlesData->Get_phi()[i]) / dphi + 1;
            particlesData->Get_phi()[i] = particlesData->Get_phi()[i] + dphi * c;
        }

        if (particlesData->Get_phi()[i] > border2)
        {
            int c                       = (particlesData->Get_phi()[i] - border2) / dphi + 1;
            particlesData->Get_phi()[i] = particlesData->Get_phi()[i] - dphi * c;
        }

        if (particlesData->Get_phi()[i] < 1.5 * PI())
            particlesData->Get_phi()[i] =
                1.5 * PI() + (1.5 * PI() - particlesData->Get_phi()[i]);

        PointType dd = std::abs(phiold - particlesData->Get_phi()[i]);

        if ((std::abs(dPh) - dd) > std::abs(dPh) * 0.1)
        {
            if (dPh < 0)
                particlesData->Get_isPeriodical()[i] = 1;
            else
                particlesData->Get_isPeriodical()[i] = 2;
        };
    };
    // particlesData->PeriodicalEvent(dPh);
    Dmath::Polar2Cartesian(&particlesData->Get_r()[0], &particlesData->Get_phi()[0],
                           &particlesData->Get_x()[0], &particlesData->Get_y()[0],
                           particlesData->NParticles());
    Dmath::Polar2Cartesian(&particlesData->Get_r()[0], &particlesData->Get_phiReal()[0],
                           &particlesData->Get_xReal()[0], &particlesData->Get_yReal()[0],
                           particlesData->NParticles());
};

template <class PointType>
void ParticlesMover<PointType>::updateMomentums(
    const std::shared_ptr<Particles2dpolar<PointType>>& particlesData,
    particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow, int i1, int i2,
    int thread)
{
    PointType gamma;
    PointType prn;
    PointType pphin;
    PointType pphi0;
    PointType pr0;
    PointType pphi;
    PointType r;
    PointType r12;

    double dphi    = 0.058177641733144;
    double border1 = 1.5 * PI() - dphi / 2;
    double border2 = 1.5 * PI() + dphi / 2;

#pragma simd
    for (int i = i1; i < i2; i++)
    {

        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            DT[thread] = DT[thread] / 2;

        r     = particlesData->Get_r()[i];
        prn   = particlesData->Get_pr()[i];
        pphin = particlesData->Get_pphi()[i];

        pr0 = prn +
              DT[thread] * std::abs(alpha) * Dmath::sign(particlesData->q[i]) *
                  particlesDataTmp.Get_Er()[i - i1] / 2;
        pphi0 = pphin +
                r * DT[thread] * std::abs(alpha) * Dmath::sign(particlesData->q[i]) *
                    particlesDataTmp.Get_Ephi()[i - i1] / 2;
        // particlesData->gamma[i] = sqrt(1 + pr0*pr0 + (pphi0 / r)*(pphi0 / r));

        PointType gamma = sqrt(1 + pr0 * pr0 + (pphi0 / r) * (pphi0 / r));
        // particlesData->Get_Er()[i] = 0;

        //	particlesData->Ephi[i] = 0;

        particlesData->Get_pr()[i] = prn +
                                     DT[thread] * (std::abs(alpha) * Dmath::sign(particlesData->q[i]) *
                                                       particlesDataTmp.Get_Er()[i - i1] +
                                                   (pphi0 * pphi0) / (r * r * r * gamma));

        if (prn < 0 && particlesData->Get_pr()[i] > 0)
        {
            int tt = 0;
        }

        PointType ph = particlesData->Get_phiReal()[i];

        if (ph < border1)
        {
            int c = (border1 - ph) / dphi + 1;
            ph    = ph + dphi * c;
        }

        if (ph > border2)
        {
            int c = (ph - border2) / dphi + 1;
            ph    = ph - dphi * c;
        }

        if (ph < 1.5 * PI())
            particlesDataTmp.Get_Ephi()[i - i1] = -particlesDataTmp.Get_Ephi()[i - i1];

        particlesData->Get_pphi()[i] = pphin +
                                       DT[thread] * std::abs(alpha) * Dmath::sign(particlesData->q[i]) *
                                           particlesDataTmp.Get_Ephi()[i - i1];

        if (particlesData->flagEmitted[i] < 1)
        {
            DT[thread] = DT[thread] * 2;
            continue;
        };
    };
};

template <class PointType>
void ParticlesMover<PointType>::init(std::vector<PointType> alpha, int solverType)
{
    double timeStep = params[0];
    for (int thread = 0; thread < timeSteps.size(); thread++)
        timeSteps[thread].resize(alpha.size());
    if (solverType == 0)
    {
        for (int i = 0; i < alpha.size(); i++)
            for (int thread          = 0; thread < timeSteps.size(); thread++)
                timeSteps[thread][i] = timeStep;
        return;
    }

    if (solverType == 1)
    {
        PointType amin = std::abs(alpha[0]);
        int       imin = 0;
        for (int i = 1; i < alpha.size(); i++)
        {
            if (alpha[i] < amin)
            {
                imin = i;
                amin = std::abs(alpha[i]);
            };
        };
        for (int i = 0; i < alpha.size(); i++)
        {
            for (int thread          = 0; thread < timeSteps.size(); thread++)
                timeSteps[thread][i] = 0.5 * timeStep * pow(amin / std::abs(alpha[i]), 0.4);
        }
        return;
    }
};

template <class PointType>
void ParticlesMover<PointType>::updatePositions(
    const std::shared_ptr<Particles2d<PointType>>& particlesData,
    const std::shared_ptr<Particles2d<PointType>>& particlesDataTmp, int nlow, int i1, int i2,
    int thread)
{
    PointType gamma;
    PointType pxn;
    PointType pyn;
    PointType pzn;

    PointType beta;

#pragma simd
    for (int i = i1; i < i2; i++)
    {
        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            continue;

        pxn = particlesData->Get_px()[i];
        pyn = particlesData->Get_py()[i];
        pzn = particlesData->Get_pz()[i];

        gamma = sqrt(1 + pxn * pxn + pyn * pyn + pzn * pzn);

        particlesData->Get_x()[i] = particlesData->Get_x()[i] + DT[thread] * pxn / gamma;
        particlesData->Get_y()[i] = particlesData->Get_y()[i] + DT[thread] * pyn / gamma;
        particlesData->Get_z()[i] = particlesData->Get_z()[i] + DT[thread] * pzn / gamma;
    };
};
template <class PointType>
void ParticlesMover<PointType>::updateMomentums(
    const std::shared_ptr<Particles2d<PointType>>& particlesData,
    particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow, int i1, int i2,
    int thread)
{

    for (int i = i1; i < i2; i++)
    {
        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            DT[thread] = DT[thread] / 2;

        //	particlesData->gamma[i] = sqrt(1 + particlesData->Get_px()[i] *
        // particlesData->Get_px()[i] +
        // particlesData->Get_py()[i] * particlesData->Get_py()[i]);

        particlesData->Get_px()[i] =
            particlesData->Get_px()[i] + DT[thread] * alpha * particlesDataTmp.Get_Ex()[i - i1];
        particlesData->Get_py()[i] =
            particlesData->Get_py()[i] + DT[thread] * alpha * particlesDataTmp.Get_Ey()[i - i1];

        if (particlesData->flagEmitted[i] < 1)
        {
            DT[thread] = DT[thread] * 2;
            continue;
        };
    };
};

template <class PointType>
void ParticlesMover<PointType>::updatePositions(
    const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
    const std::shared_ptr<Particles3dcil<PointType>>& particlesDataTmp, int nlow, int i1, int i2,
    int thread)
{
    PointType gamma;
    PointType prn;
    PointType pzn;
    PointType pphi;
    PointType r;
    PointType r12;

// updatePositionsVector(&particlesData->Get_r()[0], &particlesData->Get_z()[0],
// &particlesData->Get_phi()[0],
// &particlesData->Get_pr()[0], &particlesData->Get_pz()[0], &particlesData->Get_pphi()[0],
// particlesData->NParticles(),
// DT[thread], &tmp1[0], &tmp2[0]);

#pragma simd
    for (int i = i1; i < i2; i++)
    {
        if (particlesData->cellsNumbers[i] == -1)
            continue;

        if (particlesData->flagEmitted[i] < 1)
            continue;

        if (i == 1545)
        {
            int tt = 0;
        }

        if (std::abs(particlesData->Get_pr()[i]) > 1e-3)
        {
            int tt = 0;
        }

        r    = particlesData->Get_r()[i];
        prn  = particlesData->Get_pr()[i];
        pzn  = particlesData->Get_pz()[i];
        pphi = particlesData->Get_pphi()[i];

        gamma = sqrt(1 + prn * prn + pzn * pzn + (pphi / r) * (pphi / r));

        particlesData->Get_r()[i] = r + DT[thread] * prn / gamma;

        r12   = (particlesData->Get_r()[i] + r) / 2;
        gamma = sqrt(1 + prn * prn + pzn * pzn + (pphi / r12) * (pphi / r12));

        particlesData->Get_z()[i] = particlesData->Get_z()[i] + DT[thread] * pzn / gamma;

        particlesData->Get_phi()[i] =
            particlesData->Get_phi()[i] + DT[thread] * pphi / (gamma * r12 * r12);

        if (particlesData->Get_r()[i] < 0)
        {
            particlesData->Get_r()[i]    = -particlesData->Get_r()[i];
            particlesData->Get_pr()[i]   = -particlesData->Get_pr()[i];
            particlesData->Get_pphi()[i] = 0;
        }
    };
};
template <class PointType>
void ParticlesMover<PointType>::updateMomentums(
    const std::shared_ptr<Particles3dcil<PointType>>& particlesData,
    particlesFields<PointType>& particlesDataTmp, PointType alpha, int nlow, int i1, int i2,
    int thread)
{
    // return;
    PointType gamma;
    PointType gamma0;

    PointType pr0;
    PointType pz0;
    PointType B;
    PointType A_r;
    PointType A_z;
    PointType r;
    PointType prn;
    PointType pzn;

    double dphi    = 0.058177641733144;
    double border1 = 1.5 * PI() - dphi / 2;
    double border2 = 1.5 * PI() + dphi / 2;

    for (int i = i1; i < i2; i++)
    {

        //	particlesDataTmp.Get_Ez()[i] = 0;
        //	particlesDataTmp.Get_Er()[i] = 0;
        //	particlesDataTmp.Get_Bphi()[i] = 0;

        if (particlesData->cellsNumbers[i] == -1)
            continue;

        r   = particlesData->Get_r()[i];
        prn = particlesData->Get_pr()[i];
        pzn = particlesData->Get_pz()[i];

        if (particlesData->flagEmitted[i] < 1)
            DT[thread] = DT[thread] / 2;

        // particlesData->Get_pr()[i] = (prn + DT[thread] * alpha*(particlesDataTmp.Get_Er()[i - i1]
        // -
        // particlesDataTmp.Get_Bphi()[i - i1]*LIGHT_VELOCITY()*pzn/ gamma));
        // particlesData->Get_pz()[i] = (pzn +
        // DT[thread] * alpha*(particlesDataTmp.Get_Ez()[i - i1] + particlesDataTmp.Get_Bphi()[i -
        // i1]*
        // LIGHT_VELOCITY()*prn / gamma));

        PointType gamma =
            sqrt(1 + prn * prn + pzn * pzn +
                 (particlesData->Get_pphi()[i] / r) * (particlesData->Get_pphi()[i] / r));

        pr0 = prn +
              (DT[thread] / 2) * (alpha * particlesDataTmp.Get_Er()[i - i1] +
                                  (particlesData->Get_pphi()[i] * particlesData->Get_pphi()[i]) /
                                      (r * r * r * gamma));
        pz0   = pzn + DT[thread] * alpha * particlesDataTmp.Get_Ez()[i - i1] / 2;
        gamma = sqrt(1 + pr0 * pr0 + pz0 * pz0 +
                     (particlesData->Get_pphi()[i] / r) * (particlesData->Get_pphi()[i] / r));
        B = alpha * LIGHT_VELOCITY() * particlesDataTmp.Get_Bphi()[i - i1] / (2 * gamma);
        A_r = alpha * particlesDataTmp.Get_Er()[i - i1] +
              particlesData->Get_pphi()[i] * particlesData->Get_pphi()[i] / (r * r * r * gamma);
        A_z = alpha * particlesDataTmp.Get_Ez()[i - i1];

        particlesData->Get_pr()[i] =
            (prn + DT[thread] * (A_r - (2 * pzn + DT[thread] * (A_z + prn * B)) * B)) /
            (1 + DT[thread] * DT[thread] * B * B);
        particlesData->Get_pz()[i] =
            (pzn + DT[thread] * (A_z + (2 * prn + DT[thread] * (A_r - pzn * B)) * B)) /
            (1 + DT[thread] * DT[thread] * B * B);
        particlesData->Get_pphi()[i] = particlesData->Get_pphi()[i];

        /*PointType ph = particlesData->Get_phi()[i];

        if (ph < border1)
        {
                int c = (border1 - ph) / dphi + 1;
                ph = ph + dphi * c;
        }


        if (ph > border2)
        {
                int c = (ph - border2) / dphi + 1;
                ph = ph - dphi * c;
        }


        if (ph < 1.5 * PI())
                particlesData->Ephi[i] = -particlesData->Ephi[i];


        particlesData->Get_pphi()[i] = particlesData->Get_pphi()[i] + DT[thread] *
        alpha*particlesData->Ephi[i];*/

        if (particlesData->flagEmitted[i] < 1)
        {
            DT[thread] = DT[thread] * 2;
            continue;
        };
    };
};