#include "MagneticFieldSolver.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "Geom.h"
#include "GridData.h"

template class MagneticFieldSolver<float>;
template class MagneticFieldSolver<double>;

template <class PointType>
void MagneticFieldSolver<PointType>::FieldSimulate(
    std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType dt)
{
    PointType r0;
    PointType r;
    PointType kphN = 20;
    PointType phi0 = 0;
    PointType phi  = 0;
    PointType dphi = 2 * commtools::PI() / kphN;
    PointType tmp1;
    PointType tmp2;

    PointType cosPh;
    PointType sinPh;

    PointType xtmp;
    PointType K = commtools::VACUUM_PERMEABILITY() / (4 * commtools::PI());

    PointType dlz;

    PointType tB;

    int j      = 0;
    int jstart = 0;

    if (electrodes.size() != 0)
    {
        for (int k                 = 0; k < gridData->Getr().size(); k++)
            gridData->GetBphi()[k] = 0;
    }

    for (int i = 0; i < electrodes.size(); i++)
    {
        if (electrodes[i]->averageCollectedCurrentDensity.back() == 0)
            continue;
        for (int k = 0; k < gridData->Getr().size(); k++)
        {
            if (k == 15000)
            {
                int tt = 0;
            };
            jstart = 0;
            if (k > 1 && gridData->Getz()[k] == gridData->Getz()[k - 1])
                jstart = j;
            for (j = jstart; j < electrodes[i]->ElectrodeEdges.size() - 1; j++)
            {

                if ((gridData->Getz()[k] >= electrodes[i]->ElectrodeEdges[j].point1.y &&
                     gridData->Getz()[k] <= electrodes[i]->ElectrodeEdges[j + 1].point1.y) ||
                    (gridData->Getz()[k] <= electrodes[i]->ElectrodeEdges[j].point1.y &&
                     gridData->Getz()[k] >= electrodes[i]->ElectrodeEdges[j + 1].point1.y))
                {
                    if (gridData->Getr()[k] < electrodes[i]->ElectrodeEdges[j].point1.x)
                        break;

                    dlz = (electrodes[i]->ElectrodeEdges[j + 1].point1.y -
                           electrodes[i]->ElectrodeEdges[j].point1.y);

                    //		gridData->GetBphi()[k] = gridData->GetBphi()[k] + (dlz /
                    // electrodes[i]->ElectrodeEdges[j
                    //+  1].length())*commtools::VACUUM_PERMEABILITY()*tmp[j] / (2 *
                    // commtools::PI()*gridData->Getr()[k]);
                    volatile PointType I = electrodes[i]->averageCollectedCurrentDensitySim[j] * 2 *
                                           commtools::PI() *
                                           electrodes[i]->ElectrodeEdges[j].point1.x;

                    gridData->GetBphi()[k] = gridData->GetBphi()[k] +
                                             commtools::VACUUM_PERMEABILITY() * I /
                                                 (2 * commtools::PI() * gridData->Getr()[k]);
                    break;
                }
            }
        }
    };

    /*	for (int i = 0; i < electrodes.size(); i++)
            {
                    for (int j = 1; j < electrodes[i]->CurrentAtPoints.size(); j++)
                            electrodes[i]->CurrentAtPoints[j] = electrodes[i]->CurrentAtPoints[j] /
       dt;

                    for (int j = 1; j < electrodes[i]->CurrentAtPoints.size(); j++)
                            electrodes[i]->CurrentAtPoints[j] = electrodes[i]->CurrentAtPoints[j] +
       electrodes[i]->CurrentAtPoints[j-1];

            //	if (electrodes[i]->CurrentAtPoints.back() == 0)
            //		continue;

                    for (int j = 0; j < electrodes[i]->CurrentAtPoints.size()-1; j++)
                    {
                            for (int kphi = 0; kphi < kphN; kphi++)
                            {
                                    phi = phi0 + kphi*dphi + dphi/2;
                                    cosPh=std::cos(phi);
                                    sinPh=std::sin(phi);
                                    tmp2 = (electrodes[i]->ElectrodeEdges[j].point1.x*sinPh);
                                    tmp2 = tmp2*tmp2;
                                    for (int k = 0; k < gridData->Getr().size(); k++)
                                    {
                                            tmp1 = (electrodes[i]->ElectrodeEdges[j].point1.x*cosPh
       - gridData->Getr()[k]); tmp1 = tmp1 + 0.001*Dmath::sign(tmp1);

                                            r0 = sqrt(tmp1*tmp1 + tmp2);

                                    //	if (r0 > 0.05)
                                    //		continue;

                                            r = sqrt(tmp1*tmp1 + tmp2 +
       (electrodes[i]->ElectrodeEdges[j].point1.y -
       gridData->Getz()[k])*(electrodes[i]->ElectrodeEdges[j].point1.y - gridData->Getz()[k]));




                                            dlz = (electrodes[i]->ElectrodeEdges[j + 1].Middle().y -
       electrodes[i]->ElectrodeEdges[j].Middle().y);

                                            electrodes[i]->CurrentAtPoints[j] = 1;
                                            tB = -Dmath::sign(tmp1)*dlz*
       electrodes[i]->CurrentAtPoints[j] * r0 / (kphN*r*r*r);

                                            gridData->GetBphi()[k] = gridData->GetBphi()[k] + tB;
                                    };
                            };
                    };
            };
            for (int k = 0; k < gridData->Getr().size(); k++)
            {
                    gridData->GetBphi()[k] = gridData->GetBphi()[k] * K;
            };*/
};

template <class PointType>
void MagneticFieldSolver<PointType>::FieldSimulate(
    std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
    const std::shared_ptr<GridData2d<PointType>>& gridData, PointType dt){

};

template <class PointType>
void MagneticFieldSolver<PointType>::FieldSimulate(
    std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, PointType dt){
    /*for (int i = 0; i < electrodes.size(); i++)
    {
            std::vector<PointType> tmp = electrodes[i]->CurrentAtPoints[0];

            for (int j = 1; j < electrodes[i]->CurrentAtPoints[0].size(); j++)
                    tmp[j] = tmp[j] + tmp[j - 1];


            electrodes[i]->SetCurrent(std::abs(tmp.back()));
    }*/
};