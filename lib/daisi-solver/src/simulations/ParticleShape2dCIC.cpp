#include "ParticleShape2dCIC.h"
#include "BoundaryContainer2d.h"
#include "GridData.h"
#include "WcalculateVector.h"

#ifdef USE_BLAS
#include <cblas.h>
#else
#include <mkl.h>
#endif

template <class PointType>
int ParticleShape2dCIC<PointType>::InitEmCells(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, PointType r1,
    PointType z1, int index, int flagInit)
{

    int       n1  = 0;
    int       n   = 0;
    PointType rr1 = r1;
    PointType zz1 = z1;

    for (int i = 0; i < nearCathodeVolumes.size(); i++)
    {
        int flagBreak = 0;

        for (int j = 0; j < nearCathodeVolumes[i].size(); j++)
        {
            if (0 == nearCathodeVolumes[i][j].flagCurrentLimited)
            {
                if (nearCathodeVolumes[i][j].InCell(rr1, zz1))
                {
                    int flagAdd = 1;
                    for (int k = 0; k < emCells[index].size(); k++)
                    {
                        if (emCells[index][k] == n1)
                            flagAdd = 0;
                    }
                    if (flagAdd && flagInit)
                        emCells[index].push_back(n1);
                    return n;
                }
                n++;
            }
            n1++;
        }
    }
    return -1;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::ChargeCalculate(
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    arma::mat& Icoef, PointType kI, unsigned char flag,
    const std::vector<unsigned int>& emissionCells, unsigned int startCell, PointType r1,
    PointType z1, int index1, const PointType* W1, PointType r2, PointType z2, int index2,
    const PointType* W2, PointType I, PointType dt, PointType* rho, PointType cartesianX1,
    PointType cartesianX2, int flagType, PointType cartesianX12, PointType cartesianX22)
{
    PointType rInt;
    PointType zInt;
    int       flagSearch = 1;

    if (emType == 0)
    {
        for (int kk = 0; kk < emCells[index1].size(); kk++)
        {
            int n  = emCells[index1][kk];
            int nn = n;
            int i  = 0;
            while (1)
            {
                nn = nn - nearCathodeVolumes[i].size();
                if (nn < 0)
                    break;
                i++;
            }

            int j = nn + nearCathodeVolumes[i].size();

            for (int s = 0; s < i; s++)
            {
                if (nearCathodeVolumes[s][0].flagCurrentLimited)
                    n = n - nearCathodeVolumes[s].size();
            }
            if (nearCathodeVolumes[i][j].InCell(cartesianX1, cartesianX2) &&
                0 == nearCathodeVolumes[i][j].flagCurrentLimited)
            {
                if (!nearCathodeVolumes[i][j].InCell(cartesianX12, cartesianX22))
                {
                    int tt = 0;
                }

                Icoef(n, startCell) = Icoef(n, startCell) + kI * dt;
                break;
            }
        }
    }

    if (index1 == index2)
    {
        rho[index1]                = rho[index1] + (W1[0] + W2[0]) * I * dt / 2;
        rho[index1 + 1]            = rho[index1 + 1] + (W1[1] + W2[1]) * I * dt / 2;
        rho[levelHigh[index1]]     = rho[levelHigh[index1]] + (W1[2] + W2[2]) * I * dt / 2;
        rho[levelHigh[index1] + 1] = rho[levelHigh[index1] + 1] + (W1[3] + W2[3]) * I * dt / 2;

        return;
    }
    int flagNeigh = 0;

    if (index1 + 1 == index2 && levelHigh[index1] + 1 == levelHigh[index2])
    {
        flagNeigh = 1;
        rInt      = x1Array[index1 + 1];
        if (std::abs(r1 - r2) < 1e-12)
            zInt = z2;
        else
            zInt = z1 + (z2 - z1) * (rInt - r1) / (r2 - r1);
    }

    if (index1 == index2 + 1 && levelHigh[index1] == levelHigh[index2] + 1)
    {
        flagNeigh = 1;
        rInt      = x1Array[index1];
        if (std::abs(r1 - r2) < 1e-12)
            zInt = z2;
        else
            zInt = z1 + (z2 - z1) * (rInt - r1) / (r2 - r1);
    }

    if (levelHigh[index1] == index2 && levelHigh[index1] + 1 == index2 + 1)
    {
        flagNeigh = 1;
        zInt      = x2Array[levelHigh[index1]];
        if (std::abs(z1 - z2) < 1e-12)
            rInt = r2;
        else
            rInt = r1 + (r2 - r1) * (zInt - z1) / (z2 - z1);
    }

    if (index1 == levelHigh[index2] && index1 + 1 == levelHigh[index2] + 1)
    {
        flagNeigh = 1;
        zInt      = x2Array[index1];
        if (std::abs(z1 - z2) < 1e-12)
            rInt = r2;
        else
            rInt = r1 + (r2 - r1) * (zInt - z1) / (z2 - z1);
    }

    if (levelHigh[index1] + 1 == index2)
    {
        flagNeigh = 1;
        zInt      = x2Array[index2];
        rInt      = x1Array[index2];
    }

    if (levelHigh[index1] == index2 + 1)
    {
        flagNeigh = 1;
        zInt      = x2Array[index2 + 1];
        rInt      = x1Array[index2 + 1];
    }

    if (index1 == levelHigh[index2] + 1)
    {
        flagNeigh = 1;
        zInt      = x2Array[index1];
        rInt      = x1Array[index1];
    }

    if (index1 + 1 == levelHigh[index2])
    {
        flagNeigh = 1;
        zInt      = x2Array[index1 + 1];
        rInt      = x1Array[index1 + 1];
    }

    if (flagNeigh == 0)
    {
        int i1, i2;

        if (index1 < index2)
        {
            i1 = index1;
            i2 = index2;
        }
        else
        {
            i2 = index1;
            i1 = index2;
        }

        int s;

        PointType rnew = (r1 + r2) / 2;
        PointType znew = (z1 + z2) / 2;

        for (s = i1 - 1; s < i2 + 1; s++)
        {
            if (InCell(s, rnew, znew))
                break;
        }

        if (s == i2 + 1)
            return;

        PointType W1tmp[4];

        if (flagType == 3)
        {
            WcalculatePolar(s, rnew, znew, W1tmp);
        }
        else
        {
            Wcalculate(s, rnew, znew, W1tmp);
        }

        ChargeCalculate(nearCathodeVolumes, emType, Icoef, kI, flag, emissionCells, startCell, r1,
                        z1, index1, W1, rnew, znew, s, W1tmp, I, dt / 2, rho, cartesianX1,
                        cartesianX2, flagType, (cartesianX1 + cartesianX12) / 2,
                        (cartesianX2 + cartesianX22) / 2);

        ChargeCalculate(nearCathodeVolumes, emType, Icoef, kI, flag, emissionCells, startCell, rnew,
                        znew, s, W1tmp, r2, z2, index2, W2, I, dt / 2, rho,
                        (cartesianX1 + cartesianX12) / 2, (cartesianX2 + cartesianX22) / 2, flag,
                        cartesianX12, cartesianX22);

        return;
    }

    PointType W1tmp[4];
    PointType W2tmp[4];

    PointType dt1;
    PointType dt2;

    if (flagType == 3)
    {
        WcalculatePolar(index1, rInt, zInt, W1tmp);
        WcalculatePolar(index2, rInt, zInt, W2tmp);
        Dmath::Polar2Cartesian(r1, z1, r1, z1);
        Dmath::Polar2Cartesian(r2, z2, r2, z2);
        Dmath::Polar2Cartesian(rInt, zInt, rInt, zInt);

        PointType L  = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
        PointType l1 = sqrt((r1 - rInt) * (r1 - rInt) + (z1 - zInt) * (z1 - zInt));
        dt1          = dt * l1 / L;
        dt2          = dt - dt1;
    }
    else
    {
        Wcalculate(index1, rInt, zInt, W1tmp);
        Wcalculate(index2, rInt, zInt, W2tmp);
        PointType L  = sqrt((r2 - r1) * (r2 - r1) + (z2 - z1) * (z2 - z1));
        PointType l1 = sqrt((r1 - rInt) * (r1 - rInt) + (z1 - zInt) * (z1 - zInt));
        dt1          = dt * l1 / L;
        dt2          = dt - dt1;
    }
    rho[index1]                = rho[index1] + (W1[0] + W1tmp[0]) * I * dt1 / 2;
    rho[index1 + 1]            = rho[index1 + 1] + (W1[1] + W1tmp[1]) * I * dt1 / 2;
    rho[levelHigh[index1]]     = rho[levelHigh[index1]] + (W1[2] + W1tmp[2]) * I * dt1 / 2;
    rho[levelHigh[index1] + 1] = rho[levelHigh[index1] + 1] + (W1[3] + W1tmp[3]) * I * dt1 / 2;

    rho[index2]                = rho[index2] + (W2[0] + W2tmp[0]) * I * dt2 / 2;
    rho[index2 + 1]            = rho[index2 + 1] + (W2[1] + W2tmp[1]) * I * dt2 / 2;
    rho[levelHigh[index2]]     = rho[levelHigh[index2]] + (W2[2] + W2tmp[2]) * I * dt2 / 2;
    rho[levelHigh[index2] + 1] = rho[levelHigh[index2] + 1] + (W2[3] + W2tmp[3]) * I * dt2 / 2;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::ChargeCalculate(
    unsigned int* isPeriodical, int i1, int i2,
    const std::vector<std::vector<NearCathodeVolume<PointType>>>& nearCathodeVolumes, int emType,
    arma::mat& Icoef, PointType* kI, unsigned char* flag,
    const std::vector<unsigned int>& emissionCells, unsigned int* startCell, PointType* r1,
    PointType* z1, unsigned int* index1, PointType (*W1)[9], PointType* r2, PointType* z2,
    unsigned int* index2, PointType (*W2)[9], PointType* I, PointType dt, PointType* rho,
    PointType* cartesianX1, PointType* cartesianX2, int flagType, PointType* cartesianX12,
    PointType* cartesianX22)
{
    PointType x, y;

    if (emType == 0)
    {
        for (int ii = i1; ii < i2; ii++)
        {

            if (index1[ii - i1] == -1)
                continue;

            {
                for (int kk = 0; kk < emCells[index1[ii - i1]].size(); kk++)
                {
                    int n  = emCells[index1[ii - i1]][kk];
                    int nn = n;
                    int i  = 0;
                    while (1)
                    {
                        nn = nn - nearCathodeVolumes[i].size();
                        if (nn < 0)
                            break;
                        i++;
                    }

                    int j = nn + nearCathodeVolumes[i].size();

                    for (int s = 0; s < i; s++)
                    {
                        if (nearCathodeVolumes[s][0].flagCurrentLimited)
                            n = n - nearCathodeVolumes[s].size();
                    };

                    if (flagType == 3)
                    {
                        Dmath::Polar2Cartesian(r1[ii], z1[ii], x, y);
                    }
                    else
                    {
                        x = cartesianX1[ii];
                        y = cartesianX2[ii];
                    }
                    if (nearCathodeVolumes[i][j].InCell(x, y) &&
                        0 == nearCathodeVolumes[i][j].flagCurrentLimited)
                    {
                        Icoef(n, startCell[ii]) = Icoef(n, startCell[ii]) + kI[ii] * dt;
                        break;
                    }
                };
            }
        }
    }

    //#pragma vector always
    for (int i = i1; i < i2; i++)
    {

        if (index1[i - i1] == -1 || index2[i - i1] == -1)
            continue;
        /*	rho[index1[i - i1]] = rho[index1[i - i1]] + (W1[i - i1][0]) *I[i] * dt / 2;
                rho[index1[i - i1] + 1] = rho[index1[i - i1] + 1] + (W1[i - i1][1]) * I[i] * dt / 2;
                rho[levelHigh[index1[i - i1]]] = rho[levelHigh[index1[i - i1]]] + (W1[i - i1][2]) *
           I[i] * dt / 2;
                rho[levelHigh[index1[i - i1]] + 1] = rho[levelHigh[index1[i - i1]] + 1] + (W1[i -
           i1][3] ) * I[i] * dt /
           2;

                rho[index2[i - i1]] = rho[index2[i - i1]] + (W2[i - i1][0]) *I[i] * dt / 2;
                rho[index2[i - i1] + 1] = rho[index2[i - i1] + 1] + (W2[i - i1][1]) * I[i] * dt / 2;
                rho[levelHigh[index2[i - i1]]] = rho[levelHigh[index2[i - i1]]] + (W2[i - i1][2]) *
           I[i] * dt / 2;
                rho[levelHigh[index2[i - i1]] + 1] = rho[levelHigh[index2[i - i1]] + 1] + (W2[i -
           i1][3]) * I[i] * dt /
           2;

                continue;*/

        if (index1[i - i1] == index2[i - i1])
        {
            rho[index1[i - i1]] =
                rho[index1[i - i1]] + (W1[i - i1][0] + W2[i - i1][0]) * I[i] * dt / 2;
            rho[index1[i - i1] + 1] =
                rho[index1[i - i1] + 1] + (W1[i - i1][1] + W2[i - i1][1]) * I[i] * dt / 2;
            rho[levelHigh[index1[i - i1]]] =
                rho[levelHigh[index1[i - i1]]] + (W1[i - i1][2] + W2[i - i1][2]) * I[i] * dt / 2;
            rho[levelHigh[index1[i - i1]] + 1] = rho[levelHigh[index1[i - i1]] + 1] +
                                                 (W1[i - i1][3] + W2[i - i1][3]) * I[i] * dt / 2;
        }
        else
        {

            if (flagType == 3 && isPeriodical[i] != 0)
            {
                /*rho[index1[i - i1]] = rho[index1[i - i1]] + (W1[i - i1][0]) *I[i] * dt / 2;
                rho[index1[i - i1] + 1] = rho[index1[i - i1] + 1] + (W1[i - i1][1]) * I[i] * dt / 2;
                rho[levelHigh[index1[i - i1]]] = rho[levelHigh[index1[i - i1]]] + (W1[i - i1][2]) *
                I[i] * dt / 2;
                rho[levelHigh[index1[i - i1]] + 1] = rho[levelHigh[index1[i - i1]] + 1] + (W1[i -
                i1][3]) * I[i] * dt /
                2;

                rho[index2[i - i1]] = rho[index2[i - i1]] + (W2[i - i1][0]) *I[i] * dt / 2;
                rho[index2[i - i1] + 1] = rho[index2[i - i1] + 1] + (W2[i - i1][1]) * I[i] * dt / 2;
                rho[levelHigh[index2[i - i1]]] = rho[levelHigh[index2[i - i1]]] + (W2[i - i1][2]) *
                I[i] * dt / 2;
                rho[levelHigh[index2[i - i1]] + 1] = rho[levelHigh[index2[i - i1]] + 1] + (W2[i -
                i1][3]) * I[i] * dt /
                2;

                continue;*/

                int ii1, ii2;

                if (index1[i - i1] < index2[i - i1])
                {
                    ii1 = index1[i - i1];
                    ii2 = index2[i - i1];
                }
                else
                {
                    ii2 = index1[i - i1];
                    ii1 = index2[i - i1];
                }

                ii1 = 1;
                ii2 = levelHigh.size() - 1;

                int s1;

                double dphi = 0.058177641733144;

                PointType rnew1 = (r1[i] + r2[i]) / 2;
                PointType znew1 = 0;

                if (isPeriodical[i] == 1)
                    znew1 = 1.5 * PI();
                else
                    znew1 = 1.5 * PI() + dphi / 2;

                PointType x1;
                PointType y1;
                Dmath::Polar2Cartesian(rnew1, znew1, x1, y1);
                PointType W1tmp[4];

                for (s1 = ii1 - 1; s1 < ii2 + 1; s1++)
                {
                    if (InCell(s1, rnew1, znew1))
                        break;
                }

                if (s1 == ii2 + 1)
                {
                    W1tmp[0] = 0;
                    W1tmp[1] = 0;
                    W1tmp[2] = 0;
                    W1tmp[3] = 0;
                }
                else
                    WcalculatePolar(s1, rnew1, znew1, W1tmp);

                Dmath::Polar2Cartesian(rnew1, znew1, x1, y1);

                double k1 = std::abs(znew1 - z1[i]) / (std::abs(znew1 - z1[i]) + std::abs(znew1 - z2[i]));

                ChargeCalculate(nearCathodeVolumes, 1, Icoef, kI[i], flag[i], emissionCells,
                                startCell[i], r1[i], z1[i], index1[i - i1], &W1[i - i1][0], rnew1,
                                znew1, s1, W1tmp, I[i], dt * k1, rho, cartesianX1[i],
                                cartesianX2[i], flagType, x1, y1);

                ChargeCalculate(nearCathodeVolumes, 1, Icoef, kI[i], flag[i], emissionCells,
                                startCell[i], rnew1, znew1, s1, W1tmp, r2[i], z2[i], index2[i - i1],
                                &W2[i - i1][0], I[i], dt * (1 - k1), rho, x1, y1, flagType,
                                cartesianX12[i], cartesianX22[i]);

                continue;
            }

            ChargeCalculate(nearCathodeVolumes, 1, Icoef, kI[i], flag[i], emissionCells,
                            startCell[i], r1[i], z1[i], index1[i - i1], &W1[i - i1][0], r2[i],
                            z2[i], index2[i - i1], &W2[i - i1][0], I[i], dt, rho, cartesianX1[i],
                            cartesianX2[i], flagType, cartesianX12[i], cartesianX22[i]);
        }
    }
    /*

                    int flagNeigh = 0;
                    if (index1[i - i1] + 1 == index2[i - i1] && levelHigh[index1[i - i1]] + 1 ==
    levelHigh[index2[i -
    i1]])
                    {
                            flagNeigh = 1;
                            rInt = x1Array[index1[i - i1] + 1];
                            if (std::abs(r1[i] - r2[i])<1e-12)
                                    zInt = z2[i];
                            else
                                    zInt = z1[i] + (z2[i] - z1[i])*(rInt - r1[i]) / (r2[i] - r1[i]);
                    }

                    if (index1[i - i1] == index2[i - i1] + 1 && levelHigh[index1[i - i1]] ==
    levelHigh[index2[i - i1]] +
    1)
                    {
                            flagNeigh = 1;
                            rInt = x1Array[index1[i - i1]];
                            if (std::abs(r1[i] - r2[i])<1e-12)
                                    zInt = z2[i];
                            else
                                    zInt = z1[i] + (z2[i] - z1[i])*(rInt - r1[i]) / (r2[i] - r1[i]);
                    }

                    if (levelHigh[index1[i - i1]] == index2[i - i1] && levelHigh[index1[i - i1]] + 1
    == index2[i - i1] +
    1)
                    {
                            flagNeigh = 1;
                            zInt = x2Array[levelHigh[index1[i - i1]]];
                            if (std::abs(z1[i] - z2[i])<1e-12)
                                    rInt = r2[i];
                            else
                                    rInt = r1[i] + (r2[i] - r1[i])*(zInt - z1[i]) / (z2[i] - z1[i]);
                    }

                    if (index1[i - i1] == levelHigh[index2[i - i1]] && index1[i - i1] + 1 ==
    levelHigh[index2[i - i1]] +
    1)
                    {
                            flagNeigh = 1;
                            zInt = x2Array[index1[i - i1]];
                            if (std::abs(z1[i] - z2[i])<1e-12)
                                    rInt = r2[i];
                            else
                                    rInt = r1[i] + (r2[i] - r1[i])*(zInt - z1[i]) / (z2[i] - z1[i]);
                    }

                    if (levelHigh[index1[i - i1]] + 1 == index2[i - i1])
                    {
                            flagNeigh = 1;
                            zInt = x2Array[index2[i - i1]];
                            rInt = x1Array[index2[i - i1]];
                    }

                    if (levelHigh[index1[i - i1]] == index2[i - i1] + 1)
                    {
                            flagNeigh = 1;
                            zInt = x2Array[index2[i - i1] + 1];
                            rInt = x1Array[index2[i - i1] + 1];
                    }

                    if (index1[i - i1] == levelHigh[index2[i - i1]] + 1)
                    {
                            flagNeigh = 1;
                            zInt = x2Array[index1[i - i1]];
                            rInt = x1Array[index1[i - i1]];
                    }

                    if (index1[i - i1] + 1 == levelHigh[index2[i - i1]])
                    {
                            flagNeigh = 1;
                            zInt = x2Array[index1[i - i1] + 1];
                            rInt = x1Array[index1[i - i1] + 1];
                    }

                    if (flagNeigh == 0)
                    {
                            int ii1, ii2;

                            if (index1[i - i1] < index2[i - i1])
                            {
                                    ii1 = index1[i - i1];
                                    ii2 = index2[i - i1];
                            }
                            else
                            {
                                    ii2 = index1[i - i1];
                                    ii1 = index2[i - i1];
                            }

                            int s;

                            PointType rnew = (r1[i] + r2[i]) / 2;
                            PointType znew = (z1[i] + z2[i]) / 2;

                            for (s = ii1 - 1; s < ii2 + 1; s++)
                            {
                                    if (InCell(s, rnew, znew))
                                            break;
                            };

                            if (s == ii2 + 1)
                                    return;

                            PointType W1tmp[4];

                            if (flagType == 3)
                            {
                                    WcalculatePolar(s, rnew, znew, W1tmp);
                            }
                            else
                            {
                                    Wcalculate(s, rnew, znew, W1tmp);
                            };

                            ChargeCalculate(nearCathodeVolumes, emType, Icoef, kI[i], flag[i],
    emissionCells,
                                    startCell[i], r1[i], z1[i], index1[i - i1], &W1[i - i1][0],
    rnew, znew, s, W1tmp,
    I[i], dt / 2, rho, cartesianX1[i], cartesianX2[i], flagType, (cartesianX1[i] + cartesianX12[i])
    / 2, (cartesianX2[i]
    + cartesianX22[i]) / 2);

                            ChargeCalculate(nearCathodeVolumes, emType, Icoef, kI[i], flag[i],
    emissionCells,
                                    startCell[i], rnew, znew, s, W1tmp, r2[i], z2[i],
                                    index2[i - i1], &W2[i - i1][0], I[i], dt / 2, rho,
    (cartesianX1[i] +
    cartesianX12[i]) / 2, (cartesianX2[i] + cartesianX22[i]) / 2, flagType, cartesianX12[i],
    cartesianX22[i]);

                            return;
                    };

                    PointType W1tmp[4];
                    PointType W2tmp[4];

                    PointType dt1;
                    PointType dt2;

                    PointType rtmp1;
                    PointType ztmp1;

                    PointType rtmp2;
                    PointType ztmp2;

                    if (flagType == 3)
                    {
                            WcalculatePolar(index1[i - i1], rInt, zInt, W1tmp);
                            WcalculatePolar(index2[i - i1], rInt, zInt, W2tmp);
                            Dmath::Polar2Cartesian(r1[i], z1[i], rtmp1, ztmp1);
                            Dmath::Polar2Cartesian(r2[i], z2[i], rtmp2, ztmp2);
                            Dmath::Polar2Cartesian(rInt, zInt, rInt, zInt);

                            PointType L = sqrt((rtmp2 - rtmp1)*(rtmp2 - rtmp1) + (ztmp2 -
    ztmp1)*(ztmp2 - ztmp1));
                            PointType l1 = sqrt((rtmp1 - rInt)*(rtmp1 - rInt) + (ztmp1 -
    zInt)*(ztmp1 - zInt));
                            dt1 = dt*l1 / L;
                            dt2 = dt - dt1;
                    }
                    else
                    {
                            Wcalculate(index1[i - i1], rInt, zInt, W1tmp);
                            Wcalculate(index2[i - i1], rInt, zInt, W2tmp);
                            PointType L = sqrt((rtmp2 - rtmp1)*(rtmp2 - rtmp1) + (ztmp2 -
    ztmp1)*(ztmp2 - ztmp1));
                            PointType l1 = sqrt((rtmp1 - rInt)*(rtmp1 - rInt) + (ztmp1 -
    zInt)*(ztmp1 - zInt));
                            dt1 = dt*l1 / L;
                            dt2 = dt - dt1;
                    }

                    for (int ll = 0; ll < 4; ll++)
                    {
                            if (std::isnan(W1tmp[ll]) || std::isinf(W1tmp[ll]) || std::isnan(W2tmp[ll]) ||
    std::isinf(W2tmp[ll]) ||
    std::isnan(dt1) || std::isinf(dt1) || std::isnan(dt2) || std::isinf(dt2))
                            {
                                    int tt = 0;
                            }
                    };

                    for (int ll = 0; ll < 4; ll++)
                    {
                            if ((std::abs(W1tmp[ll])>1.01) || (std::abs(W2tmp[ll])>1.01) || dt1>0.1 ||
    dt2>0.1)
                            {
                                    int tt = 0;

                            };
                    }

                    rho[index1[i - i1]] = rho[index1[i - i1]] + (W1[i - i1][0] + W1tmp[0]) *I[i] *
    dt1 / 2;
                    rho[index1[i - i1] + 1] = rho[index1[i - i1] + 1] + (W1[i - i1][1] + W1tmp[1]) *
    I[i] * dt1 / 2;
                    rho[levelHigh[index1[i - i1]]] = rho[levelHigh[index1[i - i1]]] + (W1[i - i1][2]
    + W1tmp[2]) * I[i]
    * dt1 / 2; rho[levelHigh[index1[i - i1]] + 1] = rho[levelHigh[index1[i - i1]] + 1] + (W1[i -
    i1][3] + W1tmp[3]) *
    I[i] * dt1 / 2;

                    rho[index2[i - i1]] = rho[index2[i - i1]] + (W2[i - i1][0] + W2tmp[0]) *I[i] *
    dt2 / 2;
                    rho[index2[i - i1] + 1] = rho[index2[i - i1] + 1] + (W2[i - i1][1] + W2tmp[1]) *
    I[i] * dt2 / 2;
                    rho[levelHigh[index2[i - i1]]] = rho[levelHigh[index2[i - i1]]] + (W2[i - i1][2]
    + W2tmp[2]) * I[i]
    * dt2 / 2; rho[levelHigh[index2[i - i1]] + 1] = rho[levelHigh[index2[i - i1]] + 1] + (W2[i -
    i1][3] + W2tmp[3]) *
    I[i] * dt2 / 2;

            };


    };

    //calc(&levelHigh[0], i1, i2,  emType, Icoef, kI, flag,  startCell, r1, z1, index1, W1, r2, z2,
    index2, W2, dr, dz,
    I, dt, rho, cartesianX1, cartesianX2, flagType, cartesianX12, cartesianX22);


    /*
    int flagNeigh = 0;;
    if (index1 + 1 == index2 && levelHigh[index1] + 1 == levelHigh[index2])
    {
            flagNeigh = 1;
            rInt = x1Array[index1 + 1];
            if (std::abs(r1 - r2)<1e-12)
                    zInt = z2;
            else
                    zInt = z1 + (z2 - z1)*(rInt - r1) / (r2 - r1);
    }

    if (index1 == index2 + 1 && levelHigh[index1] == levelHigh[index2] + 1)
    {
            flagNeigh = 1;
            rInt = x1Array[index1];
            if (std::abs(r1 - r2)<1e-12)
                    zInt = z2;
            else
                    zInt = z1 + (z2 - z1)*(rInt - r1) / (r2 - r1);
    }

    if (levelHigh[index1] == index2 && levelHigh[index1] + 1 == index2 + 1)
    {
            flagNeigh = 1;
            zInt = x2Array[levelHigh[index1]];
            if (std::abs(z1-z2)<1e-12)
                    rInt = r2;
            else
                    rInt = r1 + (r2 - r1)*(zInt - z1) / (z2 - z1);
    }

    if (index1 == levelHigh[index2] && index1 + 1 == levelHigh[index2] + 1)
    {
            flagNeigh = 1;
            zInt = x2Array[index1];
            if (std::abs(z1 - z2)<1e-12)
                    rInt = r2;
            else
                    rInt = r1 + (r2 - r1)*(zInt - z1) / (z2 - z1);
    }

    if (levelHigh[index1] + 1 == index2)
    {
            flagNeigh = 1;
            zInt = x2Array[index2];
            rInt = x1Array[index2];
    }

    if (levelHigh[index1] == index2 + 1)
    {
            flagNeigh = 1;
            zInt = x2Array[index2 + 1];
            rInt = x1Array[index2 + 1];
    }

    if (index1 == levelHigh[index2] + 1)
    {
            flagNeigh = 1;
            zInt = x2Array[index1];
            rInt = x1Array[index1];
    }

    if (index1 + 1 == levelHigh[index2])
    {
            flagNeigh = 1;
            zInt = x2Array[index1 + 1];
            rInt = x1Array[index1 + 1];
    }

    if (flagNeigh == 0)
    {
            int i1,i2;

            if (index1 < index2)
            {
                    i1 = index1;
                    i2 = index2;
            }
            else
            {
                    i2 = index1;
                    i1 = index2;
            }

            int s;

            PointType rnew = (r1 + r2) / 2;
            PointType znew = (z1 + z2) / 2;

            for (s = i1-1; s < i2+1; s++)
            {
                    if (InCell(s, rnew, znew))
                            break;
            };

            if (s == i2 + 1)
                    return;

            PointType W1tmp[4];

            if (flagType == 3)
            {
                    WcalculatePolar(s, rnew, znew, W1tmp);
            }
            else
            {
                    Wcalculate(s, rnew, znew, W1tmp);
            };

            ChargeCalculate(nearCathodeVolumes, emType, Icoef, kI, flag, emissionCells,
                     startCell, r1, z1, index1, W1, rnew, znew, s, W1tmp, dr, dz, I, dt / 2, rho,
    cartesianX1,
    cartesianX2, flagType, (cartesianX1 + cartesianX12) / 2, (cartesianX2 + cartesianX22) / 2);

            ChargeCalculate(nearCathodeVolumes, emType, Icoef, kI, flag, emissionCells,
                    startCell, rnew, znew, s, W1tmp, r2, z2,
                    index2, W2, dr, dz, I, dt / 2, rho, (cartesianX1 + cartesianX12) / 2,
    (cartesianX2 + cartesianX22) /
    2, flag, cartesianX12, cartesianX22);

            return;
    };

    PointType W1tmp[4];
    PointType W2tmp[4];

    PointType dt1;
    PointType dt2;

    if (flagType == 3)
    {
            WcalculatePolar(index1, rInt, zInt, W1tmp);
            WcalculatePolar(index2, rInt, zInt, W2tmp);
            Dmath::Polar2Cartesian(r1,z1,r1,z1);
            Dmath::Polar2Cartesian(r2, z2, r2, z2);
            Dmath::Polar2Cartesian(rInt, zInt, rInt, zInt);

            PointType L = sqrt((r2 - r1)*(r2 - r1) + (z2 - z1)*(z2 - z1));
            PointType l1 = sqrt((r1 - rInt)*(r1 - rInt) + (z1 - zInt)*(z1 - zInt));
            dt1 = dt*l1 / L;
            dt2 = dt - dt1;
    }
    else
    {
            Wcalculate(index1, rInt, zInt, W1tmp);
            Wcalculate(index2, rInt, zInt, W2tmp);
            PointType L = sqrt((r2 - r1)*(r2 - r1) + (z2 - z1)*(z2 - z1));
            PointType l1 = sqrt((r1 - rInt)*(r1 - rInt) + (z1 - zInt)*(z1 - zInt));
            dt1 = dt*l1 / L;
            dt2 = dt - dt1;
    }
    rho[index1] = rho[index1] + (W1[0] + W1tmp[0]) *I*dt1 / 2;
    rho[index1 + 1] = rho[index1 + 1] + (W1[1] + W1tmp[1]) * I*dt1 / 2;
    rho[levelHigh[index1]] = rho[levelHigh[index1]] + (W1[2] + W1tmp[2]) * I*dt1 / 2;
    rho[levelHigh[index1] + 1] = rho[levelHigh[index1] + 1] + (W1[3] + W1tmp[3]) * I*dt1 / 2;

    if (int(isBoundary[index1]) == 1 && emType == 0 && int(flag) == 1)
    {
            for (int k = 0; k < emissionCells.size(); k++)
            {
                    if (emissionCells[k] == index1)
                    {
                            for (int k1 = 0; k1 < emissionCells.size(); k1++)
                            {
                                    if (emissionCells[k1] == startCell)
                                    {
                                    //	Icoef(k, k1) = Icoef(k, k1) + kI *dt1;
                                            break;
                                    };
                            }
                    };
            };
    };


    rho[index2] = rho[index2] + (W2[0] + W2tmp[0]) *I*dt2 / 2;
    rho[index2 + 1] = rho[index2 + 1] + (W2[1] + W2tmp[1]) * I*dt2/ 2;
    rho[levelHigh[index2]] = rho[levelHigh[index2]] + (W2[2] + W2tmp[2]) * I*dt2 / 2;
    rho[levelHigh[index2] + 1] = rho[levelHigh[index2] + 1] + (W2[3] + W2tmp[3]) * I*dt2 / 2;

    if (int(isBoundary[index2]) == 1 && emType == 0 && int(flag) == 1)
    {
            for (int k = 0; k < emissionCells.size(); k++)
            {
                    if (emissionCells[k] == index2)
                    {
                            for (int k1 = 0; k1 < emissionCells.size(); k1++)
                            {
                                    if (emissionCells[k1] == startCell)
                                    {
                                    //	Icoef(k, k1) = Icoef(k, k1) + kI *dt2;
                                            break;
                                    };
                            }
                    };
            };
    };*/
}

template <class PointType>
void ParticleShape2dCIC<PointType>::init(int size, const PointType* x1ArrayIn,
                                         const PointType* x2ArrayIn)
{
    x1Array = x1ArrayIn;
    x2Array = x2ArrayIn;
    levelHigh.resize(size);
    isBoundary.resize(size);
    cellVolume.resize(size);
    emCells.resize(size);
    cellNumber = 0;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::init(int size, const PointType* x1ArrayIn,
                                         const PointType* x2ArrayIn, const PointType* x3ArrayIn)
{
    x1Array = x1ArrayIn;
    x2Array = x2ArrayIn;
    x3Array = x3ArrayIn;

    levelHigh.resize(size);
    levelZ.resize(size);
    isBoundary.resize(size);
    cellVolume.resize(size);
    emCells.resize(size);
    cellNumber = 0;
}

template <class PointType>
PointType ParticleShape2dCIC<PointType>::GetH1(int index)
{
    return x1Array[index + 1] - x1Array[index];
}

template <class PointType>
PointType ParticleShape2dCIC<PointType>::GetH2(int index)
{
    return x2Array[levelHigh[index]] - x2Array[index];
}

template <class PointType>
int ParticleShape2dCIC<PointType>::InCell(PointType x1, PointType x2)
{
    for (int i = 0; i < cellNumber; i++)
    {
        if (InCell(i, x1, x2))
            return i;
    }
    return -1;
}

template <class PointType>
int ParticleShape2dCIC<PointType>::InCellWithEps(PointType x1, PointType x2)
{
    for (int i = 0; i < cellNumber; i++)
    {
        if (InCellWithEps(i, x1, x2))
            return i;
    }
    return -1;
}

template <class PointType>
int ParticleShape2dCIC<PointType>::InCellWithEps(PointType x1, PointType x2, PointType x3)
{
    for (int i = 0; i < cellNumber; i++)
    {
        if (InCellWithEps(i, x1, x2, x3))
            return i;
    }
    return -1;
}

template <class PointType>
bool ParticleShape2dCIC<PointType>::InCellWithEps(int basePoint, PointType x1, PointType x2)
{
    if (levelHigh[basePoint] == -1)
        return false;

    PointType epsx1 = std::abs(1e-4 * (x1Array[basePoint + 1] - x1Array[basePoint]));
    PointType epsx2 = std::abs(1e-4 * (x2Array[levelHigh[basePoint]] - x2Array[basePoint]));

    return (x1 >= x1Array[basePoint] || std::abs(x1 - x1Array[basePoint]) < epsx1) &&
           (x1 <= x1Array[basePoint + 1] || std::abs(x1 - x1Array[basePoint + 1]) < epsx1) &&
           (x2 >= x2Array[basePoint] || std::abs(x2 - x2Array[basePoint]) < epsx2) &&
           (x2 <= x2Array[levelHigh[basePoint]] || std::abs(x2 - x2Array[levelHigh[basePoint]]) < epsx2);
}

template <class PointType>
bool ParticleShape2dCIC<PointType>::InCellWithEps(int basePoint, PointType x1, PointType x2,
                                                  PointType x3)
{
    if (levelHigh[basePoint] == -1)
        return false;

    if (levelZ[basePoint] == -1)
        return false;

    PointType epsx1 = std::abs(1e-4 * (x1Array[basePoint + 1] - x1Array[basePoint]));
    PointType epsx2 = std::abs(1e-4 * (x2Array[levelHigh[basePoint]] - x2Array[basePoint]));
    PointType epsx3 = std::abs(1e-4 * (x3Array[levelZ[basePoint]] - x3Array[basePoint]));

    return (x1 >= x1Array[basePoint] || std::abs(x1 - x1Array[basePoint]) < epsx1) &&
           (x1 <= x1Array[basePoint + 1] || std::abs(x1 - x1Array[basePoint + 1]) < epsx1) &&
           (x2 >= x2Array[basePoint] || std::abs(x2 - x2Array[basePoint]) < epsx2) &&
           (x2 <= x2Array[levelHigh[basePoint]] || std::abs(x2 - x2Array[levelHigh[basePoint]]) < epsx2) &&
           (x3 >= x3Array[basePoint] || std::abs(x3 - x3Array[basePoint]) < epsx3) &&
           (x3 <= x3Array[levelZ[basePoint]] || std::abs(x3 - x3Array[levelZ[basePoint]]) < epsx3);
}

template <class PointType>
void ParticleShape2dCIC<PointType>::AddCell(
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, unsigned int levelHighIn,
    const std::vector<int>& flagOut, std::vector<int> pointsIn,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType)
{
    levelHigh[cellNumber] = levelHighIn;
	auto const p = emCells.size();
    isBoundary[cellNumber] = 0;
    PointType hx1          = 0;
    PointType hx2          = 0;

    if (pointsIn[5] != -1 && pointsIn[7] != -1 && pointsIn[8] != -1)
    {
        if (flagOut[levelHighIn] != 2 || flagOut[levelHighIn + 1] != 2 ||
            flagOut[cellNumber + 1] != 2 || flagOut[cellNumber] != 2)
            isBoundary[cellNumber] = 1;

        for (int k = 0; k < boundarypoints.size(); k++)
        {
            if (cellNumber == boundarypoints[k] || cellNumber + 1 == boundarypoints[k] ||
                levelHighIn == boundarypoints[k] || levelHighIn + 1 == boundarypoints[k])
                isBoundary[cellNumber] = 1;
        }
    }
    else
    {
        levelHigh[cellNumber] = -1;
    }

    int flagBound = 0;
    for (int k = 0; k < 9; k++)
    {
        if (pointsIn[k] == -1)
            flagBound = 1;
    }

    if (flagBound == 0)
    {
        hx1 = x1Array[cellNumber + 1] - x1Array[cellNumber];
        hx2 = x2Array[levelHighIn] - x2Array[cellNumber];
        cellVolume[cellNumber] =
            PI() * hx2 *
            ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
             (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
    }
    else
    {

        if (pointsIn[1] != -1 && pointsIn[7] != -1)
            hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]);

        if (pointsIn[1] == -1)
            hx2 = (x2Array[pointsIn[7]] - x2Array[pointsIn[4]]) / 2;

        if (pointsIn[7] == -1)
            hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]) / 2;

        if (pointsIn[5] == -1)
        {
            hx1 = x1Array[pointsIn[4]] - x1Array[pointsIn[3]];
            cellVolume[cellNumber] =
                PI() * hx2 *
                ((x1Array[cellNumber]) * (x1Array[cellNumber]) -
                 (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
        }

        if (pointsIn[3] == -1)
        {
            hx1 = x1Array[pointsIn[5]] - x1Array[pointsIn[4]];
            cellVolume[cellNumber] =
                PI() * hx2 *
                ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                 (x1Array[cellNumber]) * (x1Array[cellNumber]));
        }

        if (pointsIn[3] != -1 && pointsIn[5] != -1)
        {
            hx1 = x1Array[pointsIn[5]] - x1Array[pointsIn[4]];
            cellVolume[cellNumber] =
                PI() * hx2 *
                ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                 (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
        }
    }

    cellNumber++;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::AddCell(
    const std::shared_ptr<GridData2d<PointType>>& gridData, unsigned int levelHighIn,
    const std::vector<int>& flagOut, std::vector<int> pointsIn,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType)
{
    levelHigh[cellNumber] = levelHighIn;

    isBoundary[cellNumber] = 0;
    PointType hx1          = 0;
    PointType hx2          = 0;

    if (pointsIn[5] != -1 && pointsIn[7] != -1 && pointsIn[8] != -1)
    {
        if (flagOut[levelHighIn] != 2 || flagOut[levelHighIn + 1] != 2 ||
            flagOut[cellNumber + 1] != 2 || flagOut[cellNumber] != 2)
            isBoundary[cellNumber] = 1;

        for (int k = 0; k < boundarypoints.size(); k++)
        {
            if (cellNumber == boundarypoints[k] || cellNumber + 1 == boundarypoints[k] ||
                levelHighIn == boundarypoints[k] || levelHighIn + 1 == boundarypoints[k])
                isBoundary[cellNumber] = 1;
        }
    }
    else
    {
        levelHigh[cellNumber] = -1;
    }

    int flagBound = 0;
    for (int k = 0; k < 9; k++)
    {
        if (pointsIn[k] == -1)
            flagBound = 1;
    }

    if (flagBound == 0)
    {
        hx1                    = x1Array[cellNumber + 1] - x1Array[cellNumber];
        hx2                    = x2Array[levelHighIn] - x2Array[cellNumber];
        cellVolume[cellNumber] = hx1 * hx2;
    }
    else
    {

        if (pointsIn[1] != -1 && pointsIn[7] != -1)
            hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]);

        if (pointsIn[1] == -1)
            hx2 = (x2Array[pointsIn[7]] - x2Array[pointsIn[4]]) / 2;

        if (pointsIn[7] == -1)
            hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]) / 2;

        if (pointsIn[5] == -1)
            hx1 = (x1Array[pointsIn[4]] - x1Array[pointsIn[3]]) / 2;

        if (pointsIn[3] == -1)
            hx1 = (x1Array[pointsIn[5]] - x1Array[pointsIn[4]]) / 2;

        /*if (pointsIn[1] == -1)
                hx2 = (x2Array[pointsIn[7]] - x2Array[pointsIn[4]]);

        if (pointsIn[7] == -1)
                hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]);

        if (pointsIn[5] == -1)
                hx1 = (x1Array[pointsIn[4]] - x1Array[pointsIn[3]]);

        if (pointsIn[3] == -1)
                hx1 = (x1Array[pointsIn[5]] - x1Array[pointsIn[4]]);

        if (pointsIn[3] != -1 && pointsIn[5] != -1)
                hx1 = (x1Array[pointsIn[5]] - x1Array[pointsIn[4]]);*/

        cellVolume[cellNumber] = hx2 * hx1;
    }

    cellNumber++;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::AddCell(
    const std::shared_ptr<GridData3d<PointType>>& gridData, unsigned int levelHighIn,
    unsigned int levelZIn, const std::vector<int>& flagOut, std::vector<int> pointsIn,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType)
{
    levelHigh[cellNumber] = levelHighIn;
    levelZ[cellNumber]    = levelZIn;

    isBoundary[cellNumber] = 0;
    PointType hx1          = 0;
    PointType hx2          = 0;
    PointType hx3          = 0;

    for (int i = 0; i < 8; i++)
    {
        if (-1 == pointsIn[i])
        {
            levelHigh[cellNumber] = -1;
            //	levelZ[cellNumber] = -1;
            break;
        }
    }

    if (-1 != levelHigh[cellNumber])
    {

        for (int i = 0; i < 8; i++)
        {
            if (-1 != pointsIn[i] && flagOut[pointsIn[i]] != 2)
            {
                isBoundary[cellNumber] = 1;
                break;
            }
        }

        if (isBoundary[cellNumber] == 0)
        {
            for (int k = 0; k < boundarypoints.size(); k++)
            {
                for (int i = 0; i < 8; i++)
                {
                    if (pointsIn[i] == boundarypoints[k])
                    {
                        isBoundary[cellNumber] = 1;
                        break;
                    }
                }
                if (1 == isBoundary[cellNumber])
                    break;
            }
        }
    }

    if (levelZIn != -1 && levelHighIn != -1)
    {
        hx1 = x1Array[cellNumber + 1] - x1Array[cellNumber];
        hx2 = x2Array[levelHighIn] - x2Array[cellNumber];
        hx3 = x3Array[levelZIn] - x3Array[cellNumber];

        cellVolume[cellNumber] = hx1 * hx2 * hx3;
    }
    cellNumber++;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::AddCell(
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, unsigned int levelHighIn,
    const std::vector<int>& flagOut, std::vector<int> pointsIn,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType)
{
    levelHigh[cellNumber] = levelHighIn;

    isBoundary[cellNumber] = 0;
    PointType hx1          = 0;
    PointType hx2          = 0;

    if (pointsIn[5] != -1 && pointsIn[7] != -1 && pointsIn[8] != -1)
    {
        if (flagOut[levelHighIn] != 2 || flagOut[levelHighIn + 1] != 2 ||
            flagOut[cellNumber + 1] != 2 || flagOut[cellNumber] != 2)
            isBoundary[cellNumber] = 1;

        for (int k = 0; k < boundarypoints.size(); k++)
        {
            if (cellNumber == boundarypoints[k] || cellNumber + 1 == boundarypoints[k] ||
                levelHighIn == boundarypoints[k] || levelHighIn + 1 == boundarypoints[k])
                isBoundary[cellNumber] = 1;
        }
    }
    else
    {
        levelHigh[cellNumber] = -1;
    }

    int flagBound = 0;
    for (int k = 0; k < 9; k++)
    {
        if (pointsIn[k] == -1)
            flagBound = 1;
    }

    if (flagBound == 0)
    {
        hx1 = x1Array[cellNumber + 1] - x1Array[cellNumber];
        hx2 = x2Array[levelHighIn] - x2Array[cellNumber];
        cellVolume[cellNumber] =
            0.5 * hx2 * ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                         (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
    }
    else
    {

        if (pointsIn[1] != -1 && pointsIn[7] != -1)
            hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]);

        if (pointsIn[1] == -1)
            hx2 = (x2Array[pointsIn[7]] - x2Array[pointsIn[4]]) / 2;

        if (pointsIn[7] == -1)
            hx2 = (x2Array[pointsIn[4]] - x2Array[pointsIn[1]]) / 2;

        if (pointsIn[5] == -1)
        {
            hx1 = x1Array[pointsIn[4]] - x1Array[pointsIn[3]];
            cellVolume[cellNumber] =
                0.5 * hx2 * ((x1Array[cellNumber]) * (x1Array[cellNumber]) -
                             (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
        }

        if (pointsIn[3] == -1)
        {
            hx1 = x1Array[pointsIn[5]] - x1Array[pointsIn[4]];
            cellVolume[cellNumber] =
                0.5 * hx2 * ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                             (x1Array[cellNumber]) * (x1Array[cellNumber]));
        }

        if (pointsIn[3] != -1 && pointsIn[5] != -1)
        {
            hx1 = x1Array[pointsIn[5]] - x1Array[pointsIn[4]];
            cellVolume[cellNumber] =
                0.5 * hx2 * ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                             (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
        }
    }
    cellNumber++;
}

template <class PointType>
bool ParticleShape2dCIC<PointType>::InCell(int basePoint, PointType x1, PointType x2)
{
    if (levelHigh[basePoint] == -1)
        return false;
    return x1 >= x1Array[basePoint] && x1 <= x1Array[basePoint + 1] && x2 >= x2Array[basePoint] &&
           x2 <= x2Array[levelHigh[basePoint]];
}

template <class PointType>
bool ParticleShape2dCIC<PointType>::InCell(int basePoint, PointType x1, PointType x2, PointType x3)
{
    if (levelHigh[basePoint] == -1)
        return false;
    if (levelZ[basePoint] == -1)
        return false;

    return x1 >= x1Array[basePoint] && x1 <= x1Array[basePoint + 1] && x2 >= x2Array[basePoint] &&
           x2 <= x2Array[levelHigh[basePoint]] && x3 >= x3Array[basePoint] &&
           x3 <= x3Array[levelZ[basePoint]];
}

template <class PointType>
void ParticleShape2dCIC<PointType>::Wcalculate(int basePoint, PointType x1, PointType x2,
                                               PointType* W)
{
    PointType H1   = x1Array[basePoint + 1] - x1Array[basePoint];
    PointType H2   = x2Array[levelHigh[basePoint]] - x2Array[basePoint];
    PointType wx11 = (x1Array[basePoint + 1] - x1) / H1;
    PointType wx12 = 1 - wx11;
    PointType wx21 = (x2Array[levelHigh[basePoint]] - x2) / H2;
    PointType wx22 = 1 - wx21;
    W[0]           = wx11 * wx21;
    W[1]           = wx12 * wx21;
    W[2]           = wx11 * wx22;
    W[3]           = wx12 * wx22;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::Wcalculate(int basePoint, PointType x1, PointType x2,
                                               PointType x3, PointType* W)
{
    PointType H1 = x1Array[basePoint + 1] - x1Array[basePoint];
    PointType H2 = x2Array[levelHigh[basePoint]] - x2Array[basePoint];
    PointType H3 = x3Array[levelZ[basePoint]] - x3Array[basePoint];

    PointType wx11 = (x1Array[basePoint + 1] - x1) / H1;
    PointType wx12 = 1 - wx11;

    PointType wx21 = (x2Array[levelHigh[basePoint]] - x2) / H2;
    PointType wx22 = 1 - wx21;

    PointType wx31 = (x3Array[levelZ[basePoint]] - x3) / H3;
    PointType wx32 = 1 - wx31;

    W[0] = wx11 * wx21 * wx31;
    W[1] = wx12 * wx21 * wx31;
    W[2] = wx11 * wx22 * wx31;
    W[3] = wx12 * wx22 * wx31;

    W[4] = wx11 * wx21 * wx32;
    W[5] = wx12 * wx21 * wx32;
    W[6] = wx11 * wx22 * wx32;
    W[7] = wx12 * wx22 * wx32;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::Wcalculate(unsigned int* basePoint, PointType* x1,
                                               PointType* x2, PointType* x3, PointType (*W)[9],
                                               int i1, int i2, PointType* p1, PointType* p2,
                                               PointType* p3, PointType* tmp1, PointType* tmp2,
                                               PointType* gamma, double& size)
{
    PointType sizeM;
    size = 1e9;
    for (int i = i1; i < i2; i++)
    {
        PointType H1 = x1Array[basePoint[i - i1] + 1] - x1Array[basePoint[i - i1]];
        PointType H2 = x2Array[levelHigh[basePoint[i - i1]]] - x2Array[basePoint[i - i1]];
        PointType H3 = x3Array[levelZ[basePoint[i - i1]]] - x3Array[basePoint[i - i1]];
        sizeM        = std::min(std::abs(H1 / p1[i]) * gamma[i],
                         std::min(std::abs(H2 / p2[i]) * gamma[i], std::abs(H3 / p3[i]) * gamma[i]));

        if (sizeM < size)
            size = sizeM;

        PointType wx11 = (x1Array[basePoint[i - i1] + 1] - x1[i]) / H1;
        PointType wx12 = 1 - wx11;

        PointType wx21 = (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / H2;
        PointType wx22 = 1 - wx21;

        PointType wx31 = (x3Array[levelZ[basePoint[i - i1]]] - x3[i]) / H3;
        PointType wx32 = 1 - wx31;

        W[i - i1][0] = wx11 * wx21 * wx31;
        W[i - i1][1] = wx12 * wx21 * wx31;
        W[i - i1][2] = wx11 * wx22 * wx31;
        W[i - i1][3] = wx12 * wx22 * wx31;

        W[i - i1][4] = wx11 * wx21 * wx32;
        W[i - i1][5] = wx12 * wx21 * wx32;
        W[i - i1][6] = wx11 * wx22 * wx32;
        W[i - i1][7] = wx12 * wx22 * wx32;
    }
}

template <class PointType>
void ParticleShape2dCIC<PointType>::Wcalculate(unsigned int* basePoint, PointType* x1,
                                               PointType* x2, PointType* x3, PointType (*W)[9],
                                               int i1, int i2)
{
    for (int i = i1; i < i2; i++)
    {
        PointType H1 = x1Array[basePoint[i - i1] + 1] - x1Array[basePoint[i - i1]];
        PointType H2 = x2Array[levelHigh[basePoint[i - i1]]] - x2Array[basePoint[i - i1]];
        PointType H3 = x3Array[levelZ[basePoint[i - i1]]] - x3Array[basePoint[i - i1]];

        PointType wx11 = (x1Array[basePoint[i - i1] + 1] - x1[i]) / H1;
        PointType wx12 = 1 - wx11;

        PointType wx21 = (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / H2;
        PointType wx22 = 1 - wx21;

        PointType wx31 = (x3Array[levelZ[basePoint[i - i1]]] - x3[i]) / H3;
        PointType wx32 = 1 - wx31;

        W[i - i1][0] = wx11 * wx21 * wx31;
        W[i - i1][1] = wx12 * wx21 * wx31;
        W[i - i1][2] = wx11 * wx22 * wx31;
        W[i - i1][3] = wx12 * wx22 * wx31;

        W[i - i1][4] = wx11 * wx21 * wx32;
        W[i - i1][5] = wx12 * wx21 * wx32;
        W[i - i1][6] = wx11 * wx22 * wx32;
        W[i - i1][7] = wx12 * wx22 * wx32;
    }
}

template <class PointType>
void ParticleShape2dCIC<PointType>::Wcalculate(unsigned int* basePoint, PointType* x1,
                                               PointType* x2, PointType (*W)[9], int i1, int i2,
                                               PointType* p1, PointType* p2, PointType* tmp1,
                                               PointType* tmp2, PointType* gamma, double& size)
{
    for (int i = i1; i < i2; i++)
    {
        PointType H1 = x1Array[basePoint[i - i1] + 1] - x1Array[basePoint[i - i1]];
        PointType H2 = x2Array[levelHigh[basePoint[i - i1]]] - x2Array[basePoint[i - i1]];
        if (std::abs(p1[i]) < 1e-12)
            tmp1[i - i1] = 1e15;
        else
            tmp1[i - i1] = std::abs(H1 / p1[i]) * gamma[i];

        if (std::abs(p2[i]) < 1e-12)
            tmp2[i - i1] = 1e15;
        else
            tmp2[i - i1] = std::abs(H2 / p2[i]) * gamma[i];

        PointType wx11 = (x1Array[basePoint[i - i1] + 1] - x1[i]) / H1;
        PointType wx12 = 1 - wx11;
        PointType wx21 = (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / H2;
        PointType wx22 = 1 - wx21;
        W[i - i1][0]   = wx11 * wx21;
        W[i - i1][1]   = wx12 * wx21;
        W[i - i1][2]   = wx11 * wx22;
        W[i - i1][3]   = wx12 * wx22;
    }
    int im1;
    int im2;
    int s  = i2 - i1;
    int ss = 1;

    im1 = std::distance(tmp1, std::min_element(tmp1, tmp1 + s));
    im2 = std::distance(tmp2, std::min_element(tmp2, tmp2 + s));

    // if (sizeof(x1[0]) == 8)
    // {
    //     im1 = idamin(&s, (double*)tmp1, &ss);
    //     im2 = idamin(&s, (double*)tmp2, &ss);
    // }
    // else
    // {
    //     im1 = isamin(&s, (float*)tmp1, &ss);
    //     im2 = isamin(&s, (float*)tmp2, &ss);
    // }

    size = std::min(std::abs(tmp1[im1]), std::abs(tmp2[im2]));
    if (size < 1e-15)
    {
        int tt = 0;
    }
}

template <class PointType>
void ParticleShape2dCIC<PointType>::Wcalculate(unsigned int* basePoint, PointType* x1,
                                               PointType* x2, PointType (*W)[9], int i1, int i2)
{
    for (int i = i1; i < i2; i++)
    {
        PointType H1   = x1Array[basePoint[i - i1] + 1] - x1Array[basePoint[i - i1]];
        PointType H2   = x2Array[levelHigh[basePoint[i - i1]]] - x2Array[basePoint[i - i1]];
        PointType wx11 = (x1Array[basePoint[i - i1] + 1] - x1[i]) / H1;
        PointType wx12 = 1 - wx11;
        PointType wx21 = (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / H2;
        PointType wx22 = 1 - wx21;
        W[i - i1][0]   = wx11 * wx21;
        W[i - i1][1]   = wx12 * wx21;
        W[i - i1][2]   = wx11 * wx22;
        W[i - i1][3]   = wx12 * wx22;
    }
}

template <class PointType>
void ParticleShape2dCIC<PointType>::WcalculatePolar(unsigned int* basePoint, PointType* x1,
                                                    PointType* x2, PointType (*W)[9], int i1,
                                                    int i2, PointType* p1, PointType* p2,
                                                    PointType* tmp1, PointType* tmp2,
                                                    PointType* gamma, double& size)
{
    for (int i = i1; i < i2; i++)
    {
        PointType H1 = x1Array[basePoint[i - i1] + 1] - x1Array[basePoint[i - i1]];
        PointType H2 = (x2Array[levelHigh[basePoint[i - i1]]] - x2Array[basePoint[i - i1]]) * x1[i];

        tmp1[i - i1] = std::abs(H1 / p1[i]) * gamma[i];
        tmp2[i - i1] = std::abs(H2 / p2[i]) * gamma[i] * x1[i];

        PointType vol1 = (x1[i] * x1[i] - x1Array[basePoint[i - i1]] * x1Array[basePoint[i - i1]]) *
                         (x2[i] - x2Array[basePoint[i - i1]]) / 2;

        PointType vol2 = (x1[i] * x1[i] - x1Array[basePoint[i - i1]] * x1Array[basePoint[i - i1]]) *
                         (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / 2;

        PointType vol3 =
            (x1Array[basePoint[i - i1] + 1] * x1Array[basePoint[i - i1] + 1] - x1[i] * x1[i]) *
            (x2[i] - x2Array[basePoint[i - i1]]) / 2;

        PointType vol4 =
            (x1Array[basePoint[i - i1] + 1] * x1Array[basePoint[i - i1] + 1] - x1[i] * x1[i]) *
            (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / 2;

        PointType V = vol1 + vol2 + vol3 + vol4;

        W[i - i1][0] = vol4 / V;
        W[i - i1][1] = vol2 / V;
        W[i - i1][2] = vol3 / V;
        W[i - i1][3] = vol1 / V;
    }
    int im1;
    int im2;
    int s  = i2 - i1;
    int ss = 1;

    im1 = std::distance(tmp1, std::min_element(tmp1, tmp1 + s));
    im2 = std::distance(tmp2, std::min_element(tmp2, tmp2 + s));

    // if (sizeof(x1[0]) == 8)
    // {
    //     im1 = idamin(&s, (double*)tmp1, &ss);
    //     im2 = idamin(&s, (double*)tmp2, &ss);
    // }
    // else
    // {
    //     im1 = isamin(&s, (float*)tmp1, &ss);
    //     im2 = isamin(&s, (float*)tmp2, &ss);
    // }

    size = std::min(std::abs(tmp1[im1 - 1]), std::abs(tmp1[im2 - 1]));
}

template <class PointType>
void ParticleShape2dCIC<PointType>::WcalculatePolar(unsigned int* basePoint, PointType* x1,
                                                    PointType* x2, PointType (*W)[9], int i1,
                                                    int i2)
{
    for (int i = i1; i < i2; i++)
    {
        PointType H1 = x1Array[basePoint[i - i1] + 1] - x1Array[basePoint[i - i1]];
        PointType H2 = (x2Array[levelHigh[basePoint[i - i1]]] - x2Array[basePoint[i - i1]]) * x1[i];
        PointType vol1 = (x1[i] * x1[i] - x1Array[basePoint[i - i1]] * x1Array[basePoint[i - i1]]) *
                         (x2[i] - x2Array[basePoint[i - i1]]) / 2;

        PointType vol2 = (x1[i] * x1[i] - x1Array[basePoint[i - i1]] * x1Array[basePoint[i - i1]]) *
                         (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / 2;

        PointType vol3 =
            (x1Array[basePoint[i - i1] + 1] * x1Array[basePoint[i - i1] + 1] - x1[i] * x1[i]) *
            (x2[i] - x2Array[basePoint[i - i1]]) / 2;

        PointType vol4 =
            (x1Array[basePoint[i - i1] + 1] * x1Array[basePoint[i - i1] + 1] - x1[i] * x1[i]) *
            (x2Array[levelHigh[basePoint[i - i1]]] - x2[i]) / 2;

        PointType V = vol1 + vol2 + vol3 + vol4;

        W[i - i1][0] = vol4 / V;
        W[i - i1][1] = vol2 / V;
        W[i - i1][2] = vol3 / V;
        W[i - i1][3] = vol1 / V;
    }
}

template <class PointType>
void ParticleShape2dCIC<PointType>::WcalculatePolar(int basePoint, PointType x1, PointType x2,
                                                    PointType* W)
{

    PointType vol1 =
        (x1 * x1 - x1Array[basePoint] * x1Array[basePoint]) * (x2 - x2Array[basePoint]) / 2;

    PointType vol2 = (x1 * x1 - x1Array[basePoint] * x1Array[basePoint]) *
                     (x2Array[levelHigh[basePoint]] - x2) / 2;

    PointType vol3 =
        (x1Array[basePoint + 1] * x1Array[basePoint + 1] - x1 * x1) * (x2 - x2Array[basePoint]) / 2;

    PointType vol4 = (x1Array[basePoint + 1] * x1Array[basePoint + 1] - x1 * x1) *
                     (x2Array[levelHigh[basePoint]] - x2) / 2;

    PointType V = vol1 + vol2 + vol3 + vol4;

    W[0] = vol4 / V;
    W[1] = vol2 / V;
    W[2] = vol3 / V;
    W[3] = vol1 / V;
}

template <class PointType>
void ParticleShape2dCIC<PointType>::ValueInterpolate(int basePoint, const PointType* W,
                                                     const std::vector<PointType>& ValueArray,
                                                     PointType&                    result)
{
    result = W[0] * ValueArray[basePoint] + W[1] * ValueArray[basePoint + 1] +
             W[2] * ValueArray[levelHigh[basePoint]] + W[3] * ValueArray[levelHigh[basePoint] + 1];
}

template <class PointType>
void ParticleShape2dCIC<PointType>::ValueInterpolate3d(int basePoint, const PointType* W,
                                                       const std::vector<PointType>& ValueArray,
                                                       PointType&                    result)
{
    result = W[0] * ValueArray[basePoint] + W[1] * ValueArray[basePoint + 1] +
             W[2] * ValueArray[levelHigh[basePoint]] + W[3] * ValueArray[levelHigh[basePoint] + 1] +
             W[4] * ValueArray[levelZ[basePoint]] + W[5] * ValueArray[levelZ[basePoint] + 1] +
             W[6] * ValueArray[levelZ[levelHigh[basePoint]]] +
             W[7] * ValueArray[levelZ[levelHigh[basePoint]] + 1];
}

template <class PointType>
void ParticleShape2dCIC<PointType>::ChargeCalculate(int basePoint, const PointType* W,
                                                    PointType particleCharge, PointType* rho)
{
    rho[basePoint]                = rho[basePoint] + W[0] * particleCharge;
    rho[basePoint + 1]            = rho[basePoint + 1] + W[1] * particleCharge;
    rho[levelHigh[basePoint]]     = rho[levelHigh[basePoint]] + W[2] * particleCharge;
    rho[levelHigh[basePoint] + 1] = rho[levelHigh[basePoint] + 1] + W[3] * particleCharge;
    /*	if (std::isnan(rho[basePoint]) || std::isnan(rho[basePoint + 1]) || std::isnan(rho[levelHigh[basePoint]])
       ||
       std::isnan(rho[levelHigh[basePoint] + 1]))
            {
                    __debugbreak();
                    int tt = 0;
            }
            if (std::isinf(rho[basePoint]) || std::isinf(rho[basePoint + 1]) ||
       std::isinf(rho[levelHigh[basePoint]]) ||
       std::isinf(rho[levelHigh[basePoint] + 1]))
            {
                    __debugbreak();
                    int tt = 0;
            }*/
}

template <class PointType>
void ParticleShape2dCIC<PointType>::ChargeCalculate3d(int basePoint, const PointType* W,
                                                      PointType particleCharge, PointType* rho)
{
    rho[basePoint]                = rho[basePoint] + W[0] * particleCharge;
    rho[basePoint + 1]            = rho[basePoint + 1] + W[1] * particleCharge;
    rho[levelHigh[basePoint]]     = rho[levelHigh[basePoint]] + W[2] * particleCharge;
    rho[levelHigh[basePoint] + 1] = rho[levelHigh[basePoint] + 1] + W[3] * particleCharge;

    rho[levelZ[basePoint]]            = rho[levelZ[basePoint]] + W[4] * particleCharge;
    rho[levelZ[basePoint] + 1]        = rho[levelZ[basePoint] + 1] + W[5] * particleCharge;
    rho[levelZ[levelHigh[basePoint]]] = rho[levelZ[levelHigh[basePoint]]] + W[6] * particleCharge;
    rho[levelZ[levelHigh[basePoint] + 1]] =
        rho[levelZ[levelHigh[basePoint] + 1]] + W[7] * particleCharge;
};

template <class PointType>
PointType ParticleShape2dCIC<PointType>::CellVolume(int index)
{
    return cellVolume[index];
}

template class ParticleShape2dCIC<float>;
template class ParticleShape2dCIC<double>;
