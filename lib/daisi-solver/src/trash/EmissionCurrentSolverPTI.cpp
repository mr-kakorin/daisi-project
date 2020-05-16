#include "EmissionCurrentSolver.h"

template class EmissionCurrentSolverPTI<float>;
template class EmissionCurrentSolverPTI<double>;

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::SetEmissionCurrent(
    const std::shared_ptr<EmitterDevice3d<PointType>>& emitter, PointType current)
{
}

template <class PointType>
std::vector<double> EmissionCurrentSolverPTI<PointType>::CalcCathodeFieldsPTI(
    ParticleSource2d<PointType>* source, const std::shared_ptr<GridData2d<PointType>>& gridData, int nE,
    int flowNumber)
{
    std::vector<double> result;
    return result;
}
template <class PointType>
std::vector<double> EmissionCurrentSolverPTI<PointType>::CalcCathodeFieldsPTI(
    ParticleSource2d<PointType>* source, const std::shared_ptr<GridData2dpolar<PointType>>& gridData,
    int nE, int flowNumber)
{
    std::vector<double> result;
    return result;
}

template <class PointType>
std::vector<double> EmissionCurrentSolverPTI<PointType>::CalcCathodeFieldsPTI(
    ParticleSource2d<PointType>* source, const std::shared_ptr<GridData2daxs<PointType>>& gridData,
    int nE, int flowNumber)
{

    CalculateCathodeFields(source, gridData, flowNumber);

    double              r0;
    double              z0;
    int                 k = 0;
    std::vector<double> result(nE);
    std::vector<int>    pointsPer(nE);

    int np = points1[flowNumber].size() / nE;

    for (int i = 0; i < points1[flowNumber].size(); i++)
    {

        r0 = points1[flowNumber][i].x;
        z0 = points1[flowNumber][i].y;

        int ind = i / np;
        if (ind > nE - 1)
            ind = nE - 1;

        result[ind] = result[ind] + CathodeFields[flowNumber][i];

        pointsPer[ind]++;

        k++;
    };
    for (int i = 0; i < result.size(); i++)
        result[i] = result[i] / pointsPer[i];

    return result;
}

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::UpdateEmissionCurrent(
    int& flagE, const std::shared_ptr<EmitterDevice2d<PointType>>& emitter,
    ParticleGridInterface<PointType>*            particleGridInterface,
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge)
{
    if (emitter->GetEmissionType() == 1)
        return;

    PointType x1;
    PointType x0;
    PointType y1;
    PointType y0;

    std::vector<double> jv(points1[flowNumber].size());

    for (int i = 0; i < points1[flowNumber].size(); i++)
    {

        x1 = points2[flowNumber][i].x;
        x0 = points1[flowNumber][i].x;

        y1 = points2[flowNumber][i].y;
        y0 = points1[flowNumber][i].y;

        double d = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));

        double V1 = gridData->interpolatePoint(x1, y1, 0, std::string("V, V"), 0);
        double V0 = gridData->interpolatePoint(x0, y0, 0, std::string("V, V"), 0);

        double dV = 0.85 * std::abs(V1 - V0);

        // if (CathodeFields[flowNumber][i] * charge < 0)
        //	dV = dV / 2;

        jv[i] = -Dmath::sign(charge) * Dmath::sign(V1 - V0) * K[flowNumber][i] * 4.0 *
                commtools::PI() * commtools::VACUUM_PERMITTIVITY() * (2.0 / 9.0) *
                sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) / (d * d);

        if (jv[i] < 0)
            jv[i] = 0;

        double decl = CathodeFields[flowNumber][i] / E0[flowNumber][i];

        if (decl > 0.2 && K[flowNumber][i] < 1.0)
            K[flowNumber][i] = K[flowNumber][i] + 0.05;
    };

    SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);
}

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::UpdateEmissionCurrent(
    int& flagE, const std::shared_ptr<EmitterDevice3d<PointType>>& emitter,
    ParticleGridInterface<PointType>*       particleGridInterface,
    const std::shared_ptr<GridData3d<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge)
{
}

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::UpdateEmissionCurrent(
    int& flagE, arma::mat Icoef, std::vector<std::shared_ptr<EmitterDevice3d<PointType>>&> emitter,
    ParticleGridInterface<PointType>*       particleGridInterface,
    const std::shared_ptr<GridData3d<PointType>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign)
{
}

template void
EmissionCurrentSolverPTI<double>::UpdateEmissionCurrent<GridData2d<double>,
                                                        std::vector<EmitterDevice2d<double>*>>(
    int& flagE, arma::mat Icoef, std::vector<EmitterDevice2d<double>*> emitter,
    ParticleGridInterface<double>*       particleGridInterface,
    const std::shared_ptr<GridData2d<double>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);

template void
EmissionCurrentSolverPTI<double>::UpdateEmissionCurrent<GridData2daxs<double>,
                                                        std::vector<EmitterDevice2daxs<double>*>>(
    int& flagE, arma::mat Icoef, std::vector<EmitterDevice2daxs<double>*> emitter,
    ParticleGridInterface<double>*          particleGridInterface,
    const std::shared_ptr<GridData2daxs<double>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);

template void
EmissionCurrentSolverPTI<double>::UpdateEmissionCurrent<GridData2dpolar<double>,
                                                        std::vector<EmitterDevice2d<double>*>>(
    int& flagE, arma::mat Icoef, std::vector<EmitterDevice2d<double>*> emitter,
    ParticleGridInterface<double>*            particleGridInterface,
    const std::shared_ptr<GridData2dpolar<double>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);

template void
EmissionCurrentSolverPTI<float>::UpdateEmissionCurrent<GridData2d<float>,
                                                       std::vector<EmitterDevice2d<float>*>>(
    int& flagE, arma::mat Icoef, std::vector<EmitterDevice2d<float>*> emitter,
    ParticleGridInterface<float>*       particleGridInterface,
    const std::shared_ptr<GridData2d<float>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);

template void
EmissionCurrentSolverPTI<float>::UpdateEmissionCurrent<GridData2daxs<float>,
                                                       std::vector<EmitterDevice2daxs<float>*>>(
    int& flagE, arma::mat Icoef, std::vector<EmitterDevice2daxs<float>*> emitter,
    ParticleGridInterface<float>*          particleGridInterface,
    const std::shared_ptr<GridData2daxs<float>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);

template void
EmissionCurrentSolverPTI<float>::UpdateEmissionCurrent<GridData2dpolar<float>,
                                                       std::vector<EmitterDevice2d<float>*>>(
    int& flagE, arma::mat Icoef, std::vector<EmitterDevice2d<float>*> emitter,
    ParticleGridInterface<float>*            particleGridInterface,
    const std::shared_ptr<GridData2dpolar<float>>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign);

template <class PointType>
template <class gridDataType, class emittersType>
void EmissionCurrentSolverPTI<PointType>::UpdateEmissionCurrent(
    int& flagE, arma::mat Icoef, emittersType emitter,
    ParticleGridInterface<PointType>* particleGridInterface,
    const std::shared_ptr<gridDataType>& gridData, int flagStopUpdateCurrent,
    std::vector<double> ChargeSign)
{
    PointType J;

    arma::vec B;

    int size = Icoef.n_cols;

    B.zeros(Icoef.n_cols);

    int k = 0;

    std::vector<PointType> v1(size);
    std::vector<PointType> v2(size);
    std::vector<PointType> v3(size);

    double Ex;
    double Ey;

    std::vector<std::vector<double>> nearCathField(points1.size());

    int N = 0;

    for (int flowNumber = 0; flowNumber < points1.size(); flowNumber++)
    {
        nearCathField[flowNumber].resize(points1[flowNumber].size());
        for (int i = 0; i < points1[flowNumber].size(); i++)
        {
            nearCathField[flowNumber][i] = 0;
            for (int ii = 1; ii < 4; ii++)
            {
                gridData->interpolatePoint(nearCathodeVolumes[flowNumber][i].fieldPointsX[ii],
                                           nearCathodeVolumes[flowNumber][i].fieldPointsY[ii], 0,
                                           Ex, Ey);
                nearCathField[flowNumber][i] = nearCathField[flowNumber][i] +
                                               (Ex * nearCathodeVolumes[flowNumber][i].normalX[ii] +
                                                Ey * nearCathodeVolumes[flowNumber][i].normalY[ii]);
            }
        }
        Dmath::SmoothMovingAverage(nearCathField[flowNumber], N);
    }

    for (int flowNumber = 0; flowNumber < points1.size(); flowNumber++)
    {

        for (int i = 0; i < points1[flowNumber].size(); i++)
        {

            PointType CellFlow = 0;
            CellFlow           = CellFlow +
                       nearCathField[flowNumber][i] * nearCathodeVolumes[flowNumber][i].Areas[1];

            //	if (CathodeFields[flowNumber][i] * ChargeSign[flowNumber] < 0)
            //			CellFlow = - CellFlow * 0.8;

            if (emitter[flowNumber]->GetEmissionType() != 1 && flagStopUpdateCurrent != 1)
            {
                B(k) = K[flowNumber][i] * CellFlow * commtools::VACUUM_PERMITTIVITY();
                k++;
            }

            double decl = CathodeFields[flowNumber][i] / E0[flowNumber][i];

            if (K[flowNumber][i] < 1.0 && decl > 0.15)
                K[flowNumber][i] = K[flowNumber][i] + 0.1;

            if (K[flowNumber][i] > 1.0)
                K[flowNumber][i] = 1.0;
            // if (decl<0)
            //	K[flowNumber][i] = K[flowNumber][i] - 0.2;

            /*	if (K[flowNumber][i] <= 0)
                            K[flowNumber][i] = 0.1;

                    if (K[flowNumber][i] > 0.9)
                            K[flowNumber][i] = 1.0;*/

            // K[flowNumber][i] = K[flowNumber][i] + 0.75*(decl - 0.15);

            // K[flowNumber][i] = 1.0;

            //	if (K[flowNumber][i] < 1.0 && decl>0.1)
            //		K[flowNumber][i] = K[flowNumber][i] + 0.05;

            //	if (decl<0.05)
            //		K[flowNumber][i] = K[flowNumber][i] - 0.05;

            //	if (K[flowNumber][i] > 1.0)
            //		K[flowNumber][i] = 1.0;

            // K[flowNumber][i] = 1.0;
        }
    };

    std::vector<int> zeros;
    for (int i = 0; i < size; i++)
    {
        v1[i] = Icoef(i, i);
        if (std::abs(v1[i]) < 1e-15)
        {
            zeros.push_back(i);
            Icoef(i, i) = 1;
        };
        v2[i] = B(i);
    }

    arma::vec sol = inv(Icoef) * B;

    for (int i = 0; i < size; i++)
    {
        v3[i] = sol(i);
    }

    for (int i = 0; i < zeros.size(); i++)
    {
        sol(zeros[i]) = 1e-10 * Dmath::sign(B(zeros[i]));
    }

    k = 0;
    for (int flowNumber = 0; flowNumber < points1.size(); flowNumber++)
    {
        std::vector<double> jv(points1[flowNumber].size());

        if (emitter[flowNumber]->GetEmissionType() == 1 || flagStopUpdateCurrent == 1)
            continue;

        for (int i = 0; i < points1[flowNumber].size(); i++)
        {
            jv[i] = sol(k) / nearCathodeVolumes[flowNumber][i].Areas[0];
            // if (jv[i] <= 0)
            //	jv[i] =
            //0.5*emitter[flowNumber]->GetParticleSource()->sourceSurface[int((istarts[flowNumber][i
            //+ 1] +
            // istarts[flowNumber][i]) / 2)].currentDensity;

            // if (nearCathField[flowNumber][i] *ChargeSign[flowNumber] < 0)
            //	jv[i] = 1e-10;
            k++;
        };

        Dmath::SmoothMovingAverage(jv, N);

        SetValueOnSource(emitter[flowNumber]->GetParticleSource(), jv, flowNumber, 1);
    };
};

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::ChangePolinom(
    const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter, std::vector<double> dJ, int p1, int p2)
{
    for (int i = p1; i < p2; i++)
        emitter->GetParticleSource()->polinom[i] =
            emitter->GetParticleSource()->polinom[i] - dJ[i - p1];

    for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++)
    {
        emitter->GetParticleSource()->sourceSurface[i].currentDensity =
            polinomCalc(emitter->GetParticleSource()->polinom,
                        emitter->GetParticleSource()->sourceSurface[i].curveLength);
    };

    double I = emitter->GetParticleSource()->GetEmissionCurrent(1);
};

template <class PointType>
PointType EmissionCurrentSolverPTI<PointType>::polinomCalc(const std::vector<PointType>& polinom,
                                                           PointType                     x)
{
    double result = 0;
    /*double dl = 0.015 / (polinom.size() - 1);
    for (int i = 0; i < polinom.size()-1; i++)
    {
    if (x > i*dl && x <= (i+1)*dl)
    {
    double w2 = (x - i*dl) / dl;
    double w1 = 1 - w2;
    return w1*polinom[i] + w2*polinom[i + 1];
    };
    }*/
    for (int i = 0; i < polinom.size(); i++)
    {
        result = result + polinom[i] * pow(std::abs(x), double(i));
    };
    return result;
};
template <class PointType>
void EmissionCurrentSolverPTI<PointType>::SetEmissionCurrent(
    const std::shared_ptr<EmitterDevice2d<PointType>>& emitter, PointType current)
{
    PointType dcurrentDensity = 0;

    dcurrentDensity = (current - emitter->GetParticleSource()->GetEmissionCurrent(0)) /
                      emitter->GetParticleSource()->length();

    emitter->GetParticleSource()->polinom[0] =
        emitter->GetParticleSource()->polinom[0] + dcurrentDensity;

    for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++)
    {
        emitter->GetParticleSource()->sourceSurface[i].currentDensity =
            polinomCalc(emitter->GetParticleSource()->polinom,
                        emitter->GetParticleSource()->sourceSurface[i].curveLength);
    };
};
template <class PointType>
void EmissionCurrentSolverPTI<PointType>::SetEmissionCurrent(
    const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter, PointType current)
{
    PointType dcurrentDensity = 0;
    dcurrentDensity           = (current - emitter->GetParticleSource()->GetEmissionCurrent(1)) /
                      (2 * commtools::PI() *
                       emitter->GetParticleSource()->sourceSurface[0].extractingEdge.Middle().x *
                       emitter->GetParticleSource()->length());

    emitter->GetParticleSource()->polinom[0] =
        emitter->GetParticleSource()->polinom[0] + dcurrentDensity;

    // for (int i = 0; i < emitter->GetParticleSource()->polinom.size(); i++)
    //	emitter->GetParticleSource()->polinom[i] = emitter->GetParticleSource()->polinom[i] +
    //dcurrentDensity;

    for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++)
    {
        emitter->GetParticleSource()->sourceSurface[i].currentDensity =
            polinomCalc(emitter->GetParticleSource()->polinom,
                        emitter->GetParticleSource()->sourceSurface[i].curveLength);
    };

    zz
};

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::UpdateEmissionCurrent(
    int& flagE, const std::shared_ptr<EmitterDevice2d<PointType>>& emitter,
    ParticleGridInterface<PointType>*       particleGridInterface,
    const std::shared_ptr<GridData2d<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge)
{
    if (emitter->GetEmissionType() == 1)
        return;

    PointType x1;
    PointType x0;
    PointType y1;
    PointType y0;

    std::vector<double> jv(points1[flowNumber].size());

    for (int i = 0; i < points1[flowNumber].size(); i++)
    {

        x1 = points2[flowNumber][i].x;
        x0 = points1[flowNumber][i].x;

        y1 = points2[flowNumber][i].y;
        y0 = points1[flowNumber][i].y;

        double d = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));

        double V1 = gridData->interpolatePoint(x1, y1, 0, std::string("V, V"), 0);
        double V0 = gridData->interpolatePoint(x0, y0, 0, std::string("V, V"), 0);

        double dV = std::abs(V1 - V0);

        // if (CathodeFields[flowNumber][i] * charge < 0)
        //	dV = dV / 2;

        jv[i] = K[flowNumber][i] * 4.0 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() *
                (2.0 / 9.0) * sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) / (d * d);

        double decl = CathodeFields[flowNumber][i] / E0[flowNumber][i];

        if (decl > 0.2 && K[flowNumber][i] < 1.0)
            K[flowNumber][i] = K[flowNumber][i] + 0.05;
    };

    SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);
};
template <class PointType>
void EmissionCurrentSolverPTI<PointType>::UpdateEmissionCurrent(
    int& flagE, const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    ParticleGridInterface<PointType>*          particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge)
{
    if (emitter->GetEmissionType() == 1)
        return;

    PointType x1;
    PointType x0;
    PointType y1;
    PointType y0;

    std::vector<double> jv(points1[flowNumber].size());

    for (int i = 0; i < points1[flowNumber].size(); i++)
    {

        x1 = points2[flowNumber][i].x;
        x0 = points1[flowNumber][i].x;

        y1 = points2[flowNumber][i].y;
        y0 = points1[flowNumber][i].y;

        double angle = atan((y1 - y0) / (x1 - x0));

        double V1 = gridData->interpolatePoint(x1, y1, 0, std::string("V, V"), 0);
        double V0 = gridData->interpolatePoint(x0, y0, 0, std::string("V, V"), 0);

        double dV = std::abs(V1 - V0);

        if (std::abs(angle) < 0.1)
        {
            double g     = x0 / x1;
            double gamma = log(g);

            double beta = gamma - 0.4 * pow(gamma, 2.0) + 0.091667 * pow(gamma, 3.0) -
                          0.014242 * pow(gamma, 4.0) + 0.001679 * pow(gamma, 5.0);

            // if (CathodeFields[flowNumber][i] * charge < 0)
            //	dV = dV / 2;

            jv[i] = 4.0 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * (2.0 / 9.0) *
                    sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) /
                    (x1 * beta * beta * 2.0 * commtools::PI() * x0);
        }
        else
        {

            double r0 = sqrt(x0 * x0 + y0 * y0);
            double r1 = sqrt(x1 * x1 + y1 * y1);

            double gs     = r0 / r1;
            double gammas = log(gs);

            double alpha = gammas - 0.3 * pow(gammas, 2.0) + 0.075 * pow(gammas, 3.0) -
                           0.014318 * pow(gammas, 4.0) + 0.0022161 * pow(gammas, 5.0);

            jv[i] = K[flowNumber][i] * 8.0 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() *
                    (2.0 / 9.0) * sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) /
                    (alpha * alpha * 4 * commtools::PI() * r0 * r0);
        }

        double decl = CathodeFields[flowNumber][i] / E0[flowNumber][i];

        if (decl > 0.2 && K[flowNumber][i] < 1.0)
            K[flowNumber][i] = K[flowNumber][i] + 0.05;
    };

    SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);
};

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::ChargeConserving(
    const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    ParticleGridInterface<PointType>*               particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge){

};

template <class PointType>
void EmissionCurrentSolverPTI<PointType>::VirtualDiode(
    int& flagE, const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter,
    ParticleGridInterface<PointType>*          particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge)
{
    PointType r;
    PointType z;

    PointType dr;
    PointType dz;
    PointType Er;
    PointType J;
    PointType Er1;
    PointType jm[2]   = {4000.0, 7000.0};
    int       k       = 0;
    int       baseOld = -1;
    PointType Jold;

    PointType ErAverage = 0;

    double Kmin = 0.99;
    double Kmax = 1.002;

    // if (flowNumber == 0)
    //	return;

    PointType r1;
    PointType z1;

    PointType r0;
    PointType z0;

    for (int i = 0; i < points1[flowNumber].size(); i++)
    {

        r1 = points2[flowNumber][i].x;
        r0 = points1[flowNumber][i].x;

        z1 = points2[flowNumber][i].y;
        z0 = points1[flowNumber][i].y;

        double g     = r0 / r1;
        double gamma = log(g);

        double beta = gamma - 0.4 * pow(gamma, 2.0) + 0.091667 * pow(gamma, 3.0) -
                      0.014242 * pow(gamma, 4.0) + 0.001679 * pow(gamma, 5.0);

        double V1 = gridData->interpolatePoint(r1, z1, 0, std::string("V, V"), 0);
        double V0 = gridData->interpolatePoint(r0, z0, 0, std::string("V, V"), 0);

        double Er;
        double Ez;

        gridData->interpolatePoint(nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
                                   nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 0, Er, Ez);
        double Er0 = (Er * nearCathodeVolumes[flowNumber][i].normalX[0] +
                      Ez * nearCathodeVolumes[flowNumber][i].normalY[0]);

        ErAverage = ErAverage + std::abs(Er0);

        double dV = std::abs(V1 - V0);

        if (Er0 < 0)
        {
            dV = dV / 2;
        }

        double jnew = 4.0 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * (2.0 / 9.0) *
                      sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) /
                      (r1 * beta * beta * 2.0 * commtools::PI() * r0);

        double omega = 0.3;

        J = jnew;

        //	J = emitter->GetParticleSource()->sourceSurface[i].currentDensity*(1 - omega) +
        //omega*jnew;

        J = std::min(1.0e6, double(J));

        double Jold;
        double k0;
        double k1;

        if (i > 0)
        {
            k0   = (istarts[flowNumber][i] + istarts[flowNumber][i - 1]) / 2;
            k1   = (istarts[flowNumber][i + 1] + istarts[flowNumber][i]) / 2;
            Jold = emitter->GetParticleSource()->sourceSurface[k0 - 1].currentDensity;

            int kend;

            if (i == points1[flowNumber].size() - 1)
                kend = istarts[flowNumber][i + 1];
            else
                kend = k1;

            for (int kk = k0; kk < kend; kk++)
            {
                double jn = (J - Jold) * (kk - k0) / (k1 - k0) + Jold;
                emitter->GetParticleSource()->sourceSurface[kk].currentDensity = jn;
            }
        }
        if (i == 0)
        {
            for (int kk = 0; kk < istarts[flowNumber][i + 1] / 2; kk++)
                emitter->GetParticleSource()->sourceSurface[kk].currentDensity = J;
        }

        k++;
    };

    ErAverage = ErAverage / k;
    emitter->GetParticleSource()->setErAverage(ErAverage);
}

template <class PointType>
double EmissionCurrentSolverPTI<PointType>::ChangeCurrent(
    const std::shared_ptr<EmitterDevice2daxs<PointType>>& emitter, int i, double per)
{
    double Itot = emitter->GetParticleSource()->GetEmissionCurrent(0);
    double dI   = Itot * per;
    double dp   = dI * (i + 1) / (pow(emitter->GetParticleSource()->length(), double(i + 1)));

    emitter->GetParticleSource()->polinom[i] = emitter->GetParticleSource()->polinom[i] + dp;

    // emitter->GetParticleSource()->polinom[i] = emitter->GetParticleSource()->polinom[i] +
    // emitter->GetParticleSource()->polinom[i] * per;

    for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++)
    {
        emitter->GetParticleSource()->sourceSurface[i].currentDensity =
            polinomCalc(emitter->GetParticleSource()->polinom,
                        emitter->GetParticleSource()->sourceSurface[i].curveLength);
    };
    //	return emitter->GetParticleSource()->polinom[i] * per;
    return dp;
};