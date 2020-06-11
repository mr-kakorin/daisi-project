#include "ParticleShape2dTSC.h"
#include "BoundaryContainer2d.h"
#include "GridData.h"
#include <Constants.h>


template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate(PointType r1, PointType z1, int index1,
                                                    const PointType* W1, PointType r2, PointType z2,
                                                    int index2, const PointType* W2, PointType dr,
                                                    PointType dz, PointType I, PointType dt,
                                                    PointType* rho)
{
    /*PointType rInt;
    PointType zInt;

    if (index1 == index2)
    {
            tauIncell = dt;
            rInt = r2;
            zInt = z2;
            rho
            return
    }
    else
    {
            if (index1 + 1 == index2 && levelHigh[index1] + 1 == levelHigh[index2])
            {
                    rInt = x1Array(index1 + 1);
                    zInt = z1 + (z2 - z1)*(rInt - r1) / (r2 - r1);
            }

            if (index1 == index2+1 && levelHigh[index1] + 1 == levelHigh[index2])
            {
                    rInt = x1Array(index1 + 1);
                    zInt = z1 + (z2 - z1)*(rInt - r1) / (r2 - r1);
            }

            if (levelHigh[index1] == index2 && levelHigh[index1] + 1 == index2 + 1)
            {
                    zInt = x2Array(levelHigh[index1]);
                    rInt = r1 + (r2 - r1)*(zInt - z1) / (z2 - z1);
            }

            if (index1 == levelHigh[index2] && index1 + 1 == levelHigh[index2] + 1)
            {
                    zInt = x2Array(index1);
                    rInt = r1 + (r2 - r1)*(zInt - z1) / (z2 - z1);
            }

            if (levelHigh[index1] + 1 == index2)
            {
                    zInt = x2Array[index2];
                    rInt = x1Array[index2];
            }

            if (levelHigh[index1] == index2+1)
            {
                    zInt = x2Array[index2 + 1];
                    rInt = x1Array[index2 + 1];
            }

            if (index1 == levelHigh[index2] + 1)
            {
                    zInt = x2Array[index1];
                    rInt = x1Array[index1];
            }

            if (index1+1 == levelHigh[index2])
            {
                    zInt = x2Array[index1 + 1];
                    rInt = x1Array[index1 + 1];
            }
    };*/
}

template <class PointType>
void ParticleShape2dTSC<PointType>::init(int size, const PointType* x1ArrayIn,
                                         const PointType* x2ArrayIn)
{
    x1Array = x1ArrayIn;
    x2Array = x2ArrayIn;
    levelHigh.resize(size);
    levelLow.resize(size);
    cellType.resize(size);
    points.resize(size);
    isBoundary.resize(size);
    cellVolume.resize(size);
    cellNumber = 0;
}

template <class PointType>
PointType ParticleShape2dTSC<PointType>::GetH1(int index)
{
    if (cellType[index] == 0)
        return x1Array[index + 1] - x1Array[index];

    if (cellType[index] == 4)
        return x1Array[points[index][5]] - x1Array[index];

    return x1Array[points[index][1]] - x1Array[points[index][0]];
}

template <class PointType>
PointType ParticleShape2dTSC<PointType>::GetH2(int index)
{
    if (cellType[index] == 0)
        return x2Array[levelHigh[index]] - x2Array[index];

    if (cellType[index] == 4)
        return x2Array[points[index][7]] - x2Array[index];

    if (cellType[index] == 2)
        return x2Array[points[index][3]] - x2Array[points[index][0]];

    return x2Array[points[index][2]] - x2Array[points[index][0]];
}

template <class PointType>
void ParticleShape2dTSC<PointType>::AddCell(
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, unsigned int levelLowIn,
    unsigned int levelHighIn, char cellTypeIn, std::vector<int> pointsIn,
    const std::vector<int>&                               flagOut,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType)
{
    levelHigh[cellNumber] = levelHighIn;
    levelLow[cellNumber]  = levelLowIn;
    cellType[cellNumber]  = cellTypeIn;
    points[cellNumber]    = pointsIn;

    isBoundary[cellNumber] = 0;
    PointType hx1          = 0;
    PointType hx2          = 0;
    if (int(cellTypeIn) == 0)
    {
        hx1 = x1Array[cellNumber] - x1Array[cellNumber - 1];
        hx2 = x2Array[cellNumber] - x2Array[levelLowIn];

        if (flagOut[levelLowIn] != 2 || flagOut[levelLowIn - 1] != 2 ||
            flagOut[levelLowIn + 1] != 2)
            isBoundary[cellNumber] = 1;
        if (flagOut[cellNumber] != 2 || flagOut[cellNumber - 1] != 2 ||
            flagOut[cellNumber + 1] != 2)
            isBoundary[cellNumber] = 1;
        if (flagOut[levelHighIn] != 2 || flagOut[levelHighIn - 1] != 2 ||
            flagOut[levelHighIn + 1] != 2)
            isBoundary[cellNumber] = 1;
        for (int k = 0; k < boundarypoints.size(); k++)
        {
            if (levelLowIn == boundarypoints[k] || levelLowIn - 1 == boundarypoints[k] ||
                levelLowIn + 1 == boundarypoints[k])
                isBoundary[cellNumber] = 1;
            if (cellNumber == boundarypoints[k] || cellNumber - 1 == boundarypoints[k] ||
                cellNumber + 1 == boundarypoints[k])
                isBoundary[cellNumber] = 1;
            if (levelHighIn == boundarypoints[k] || levelHighIn - 1 == boundarypoints[k] ||
                levelHighIn + 1 == boundarypoints[k])
                isBoundary[cellNumber] = 1;
        }
        switch (problemType)
        {

        case 0:
            cellVolume[cellNumber] = hx1 * hx2;
            break;
        case 1:
            cellVolume[cellNumber] =
                PI() * hx2 *
                ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                 (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
            break;
        }
    }
    else
    {
        isBoundary[cellNumber] = 1;
        switch (cellTypeIn)
        {
        case 1:
            hx1 = x1Array[pointsIn[1]] - x1Array[pointsIn[0]];
            hx2 = x2Array[pointsIn[2]] - x2Array[pointsIn[0]];
            switch (problemType)
            {
            case 0:
                cellVolume[cellNumber] = hx1 * hx2 / 4;
                break;
            case 1:
                if (cellNumber == pointsIn[1] || cellNumber == pointsIn[3])
                    cellVolume[cellNumber] =
                        0.5 * PI() * hx2 *
                        (x1Array[cellNumber] * x1Array[cellNumber] -
                         (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
                if (cellNumber == pointsIn[0] || cellNumber == pointsIn[2])
                    cellVolume[cellNumber] =
                        0.5 * PI() * hx2 *
                        ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                         x1Array[cellNumber] * x1Array[cellNumber]);
                break;
            }
            break;
        case 2:
            hx1 = x1Array[pointsIn[1]] - x1Array[pointsIn[0]];
            hx2 = x2Array[pointsIn[3]] - x2Array[pointsIn[0]];
            switch (problemType)
            {
            case 0:
                cellVolume[cellNumber] = hx1 * hx2 / 2;
                break;
            case 1:
                cellVolume[cellNumber] =
                    0.5 * PI() * hx2 *
                    ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                     (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
                break;
            }
            break;
        case 3:
            hx1 = x1Array[pointsIn[1]] - x1Array[pointsIn[0]];
            hx2 = x2Array[pointsIn[2]] - x2Array[pointsIn[0]];
            switch (problemType)
            {
            case 0:
                cellVolume[cellNumber] = hx1 * hx2 / 4;
                break;
            case 1:
                if (cellNumber == pointsIn[3])
                    cellVolume[cellNumber] =
                        PI() * hx2 *
                        (x1Array[cellNumber] * x1Array[cellNumber] -
                         (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
                if (cellNumber == pointsIn[2])
                    cellVolume[cellNumber] =
                        PI() * hx2 *
                        ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                         x1Array[cellNumber] * x1Array[cellNumber]);
                break;
            };
            break;
        case 4:
            hx1 = x1Array[cellNumber] - x1Array[pointsIn[3]];
            hx2 = x2Array[pointsIn[7]] - x2Array[cellNumber];
            switch (problemType)
            {
            case 0:
                cellVolume[cellNumber] = 3 * hx1 * hx2 / 4;
                break;
            case 1:
                if (pointsIn[0] == -1 || pointsIn[6] == -1)
                    cellVolume[cellNumber] =
                        PI() * hx2 *
                            (x1Array[cellNumber] * x1Array[cellNumber] -
                             (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2)) +
                        0.5 * PI() * hx2 *
                            ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                             x1Array[cellNumber] * x1Array[cellNumber]);
                if (pointsIn[2] == -1 || pointsIn[8] == -1)
                    cellVolume[cellNumber] =
                        PI() * hx2 *
                            ((x1Array[cellNumber] + hx1 / 2) * (x1Array[cellNumber] + hx1 / 2) -
                             x1Array[cellNumber] * x1Array[cellNumber]) +
                        0.5 * PI() * hx2 *
                            (x1Array[cellNumber] * x1Array[cellNumber] -
                             (x1Array[cellNumber] - hx1 / 2) * (x1Array[cellNumber] - hx1 / 2));
                break;
            }
            break;
        }
    }

    /*if (flagBoundaryCell)
    {
            if (int(cellTypeIn) == 3)
            {
                    std::vector<int> outPoints;

                    if (flagOut[pointsIn[0]] != 2)
                            outPoints.push_back(pointsIn[0]);

                    if (flagOut[pointsIn[1]] != 2)
                            outPoints.push_back(pointsIn[1]);

                    if (flagOut[pointsIn[3]] != 2)
                            outPoints.push_back(pointsIn[3]);

                    if (flagOut[pointsIn[5]] != 2)
                            outPoints.push_back(pointsIn[5]);

                    std::vector<DGeo::Edge<PointType>> cellEdges (4);

                    cellEdges[0].point1.x = x1Array[pointsIn[0]];
                    cellEdges[0].point1.y = x2Array[pointsIn[0]];
                    cellEdges[0].point1.z = 0;
                    cellEdges[0].point2.x = x1Array[pointsIn[1]];
                    cellEdges[0].point2.y = x2Array[pointsIn[1]];
                    cellEdges[0].point2.z = 0;

                    cellEdges[1].point1.x = x1Array[pointsIn[1]];
                    cellEdges[1].point1.y = x2Array[pointsIn[1]];
                    cellEdges[1].point1.z = 0;
                    cellEdges[1].point2.x = x1Array[pointsIn[5]];
                    cellEdges[1].point2.y = x2Array[pointsIn[5]];
                    cellEdges[1].point2.z = 0;

                    cellEdges[2].point1.x = x1Array[pointsIn[5]];
                    cellEdges[2].point1.y = x2Array[pointsIn[5]];
                    cellEdges[2].point1.z = 0;
                    cellEdges[2].point2.x = x1Array[pointsIn[4]];
                    cellEdges[2].point2.y = x2Array[pointsIn[4]];
                    cellEdges[2].point2.z = 0;

                    cellEdges[3].point1.x = x1Array[pointsIn[4]];
                    cellEdges[3].point1.y = x2Array[pointsIn[4]];
                    cellEdges[3].point1.z = 0;
                    cellEdges[3].point2.x = x1Array[pointsIn[0]];
                    cellEdges[3].point2.y = x2Array[pointsIn[0]];
                    cellEdges[3].point2.z = 0;

                    std::vector<DGeo::Point<PointType>> intersectionPoints;
                    for (int i = 0; i < domainBoundary->ContainerSize; i++)
                    {
                                    for (int j = 0; j < 4; j++)
                                    {
                                            DGeo::Point<PointType> tmp;
                                            if
    (domainBoundary->EdgesData[i].IsIntersectionEdge(cellEdges[j],
    cellEdges[j].length()*1e-7, &tmp)) intersectionPoints.push_back(tmp);
                                    }
                    }
                    switch (outPoints.size())
                    {
                            case 2:
                                    PointType r1 = intersectionPoints[0].x;
                                    PointType r2 = intersectionPoints[1].x;
                                    PointType vol1 = (1.0 /
    3.0)*PI()*std::abs(intersectionPoints[0].y -
    intersectionPoints[1].y)*(r1*r1 + r1*r2 + r2*r2); PointType vol2 =
    PI()*x1Array[outPoints[0]] *
    x1Array[outPoints[0]] * std::abs(x2Array[outPoints[0]] - x2Array[outPoints[1]]); PointType dVol =
    vol2 - vol1;

                                    PointType oldVol = PI()*std::abs(x2Array[outPoints[0]] -
    x2Array[outPoints[1]])*std::abs(x1Array[pointsIn[0]] * x1Array[pointsIn[0]] - x1Array[pointsIn[1]] *
    x1Array[pointsIn[1]]); PointType newVol = oldVol - dVol;
                    };
            }



    };*/

    cellNumber++;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::AddCell(
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, unsigned int levelLowIn,
    unsigned int levelHighIn, char cellTypeIn, std::vector<int> pointsIn,
    const std::vector<int>&                               flagOut,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType){}

template <class PointType>
void ParticleShape2dTSC<PointType>::AddCell(
    const std::shared_ptr<GridData2d<PointType>>& gridData, unsigned int levelLowIn,
    unsigned int levelHighIn, char cellTypeIn, std::vector<int> pointsIn,
    const std::vector<int>&                               flagOut,
    const std::shared_ptr<BoundaryContainer2d<PointType>> domainBoundary,
    const std::vector<int>& boundarypoints, int problemType){

}

template <class PointType>
int ParticleShape2dTSC<PointType>::InCell(PointType x1, PointType x2)
{
    for (int i = 0; i < cellNumber; i++)
    {
        if (InCell(i, x1, x2))
            return i;
    }
    return -1;
}

template <class PointType>
bool ParticleShape2dTSC<PointType>::InCell(int index, PointType x1, PointType x2)
{
    switch (cellType[index])
    {
    case 0:
        return InCell0(index, x1, x2);
    case 1:
        return InCell1(index, x1, x2);
    case 2:
        return InCell2(index, x1, x2);
    case 3:
        return InCell3(index, x1, x2);
    case 4:
        return InCell4(index, x1, x2);
    }
    return false;
}

template <class PointType>
bool ParticleShape2dTSC<PointType>::InCell0(int index, PointType x1, PointType x2)
{
    PointType H2 = x2Array[levelHigh[index]] - x2Array[index];
    PointType H1 = x1Array[index + 1] - x1Array[index];
    if (std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
        return true;
    return false;
}

template <class PointType>
bool ParticleShape2dTSC<PointType>::InCell1(int index, PointType x1, PointType x2)
{
    PointType H2 = x2Array[points[index][2]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    if (x1 >= x1Array[points[index][0]] && x1 <= x1Array[points[index][1]] &&
        x2 >= x2Array[points[index][0]] && x2 <= x2Array[points[index][2]] &&
        std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
        return true;
    return false;
}

template <class PointType>
bool ParticleShape2dTSC<PointType>::InCell2(int index, PointType x1, PointType x2)
{
    PointType H2 = x2Array[points[index][3]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    if (std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2 &&
        x2 <= x2Array[points[index][3]] && x2 >= x2Array[points[index][0]])
        return true;
    return false;
}

template <class PointType>
bool ParticleShape2dTSC<PointType>::InCell3(int index, PointType x1, PointType x2)
{
    PointType H2 = x2Array[points[index][2]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    if (std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2 &&
        x1 <= x1Array[points[index][1]] && x1 >= x1Array[points[index][0]])
        return true;
    return false;
}

template <class PointType>
bool ParticleShape2dTSC<PointType>::InCell4(int index, PointType x1, PointType x2)
{
    PointType H2 = x2Array[points[index][7]] - x2Array[index];
    PointType H1 = x1Array[points[index][5]] - x1Array[index];

    if (points[index][0] != -1)
    {
        if (x1 >= x1Array[points[index][0]] && x1 <= x1Array[points[index][1]] &&
            x2 >= x2Array[points[index][0]] && x2 <= x2Array[points[index][3]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
            return true;
    }

    if (points[index][2] != -2)
    {
        if (x1 >= x1Array[points[index][1]] && x1 <= x1Array[points[index][2]] &&
            x2 >= x2Array[points[index][0]] && x2 <= x2Array[points[index][4]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
            return true;
    }

    if (points[index][6] != -1)
    {
        if (x1 >= x1Array[points[index][0]] && x1 <= x1Array[points[index][1]] &&
            x2 >= x2Array[points[index][3]] && x2 <= x2Array[points[index][6]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
            return true;
    }

    if (points[index][8] != -1)
    {
        if (x1 >= x1Array[points[index][1]] && x1 <= x1Array[points[index][2]] &&
            x2 >= x2Array[points[index][3]] && x2 <= x2Array[points[index][6]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
            return true;
    }

    return false;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::WcalculatePolar(int index, PointType x1, PointType x2,
                                                    PointType* W)
{
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate(int index, PointType x1, PointType x2, PointType* W)
{
    switch (cellType[index])
    {
    case 0:
        return Wcalculate0(index, x1, x2, W);
    case 1:
        return Wcalculate1(index, x1, x2, W);
    case 2:
        return Wcalculate2(index, x1, x2, W);
    case 3:
        return Wcalculate3(index, x1, x2, W);
    case 4:
        return Wcalculate4(index, x1, x2, W);
    }
    return;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate0(int index, PointType x1, PointType x2, PointType* W)
{
    PointType H2 = x2Array[levelHigh[index]] - x2Array[index];
    PointType H1 = x1Array[index + 1] - x1Array[index];

    PointType X1 = (x1Array[index] - x1) / H1;

    PointType wx1_1 = PointType(0.5 * (0.5 + X1) * (0.5 + X1));
    PointType wx11  = PointType(0.5 * (0.5 - X1) * (0.5 - X1));
    PointType wx10  = PointType(1 - wx1_1 - wx11);
    PointType wx10T = PointType(0.75 - X1 * X1);

    PointType X2 = (x2Array[index] - x2) / H2;

    PointType wx2_1 = PointType(0.5 * (0.5 + X2) * (0.5 + X2));
    PointType wx21  = PointType(0.5 * (0.5 - X2) * (0.5 - X2));
    PointType wx20  = PointType(1 - wx2_1 - wx21);
    PointType wx20T = PointType(0.75 - X2 * X2);

    W[0] = wx2_1 * wx1_1;
    W[1] = wx2_1 * wx10;
    W[2] = wx2_1 * wx11;
    W[3] = wx20 * wx1_1;
    W[4] = wx20 * wx10;
    W[5] = wx20 * wx11;
    W[6] = wx21 * wx1_1;
    W[7] = wx21 * wx10;
    W[8] = wx21 * wx11;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate0Vector(int* index, PointType* x1, PointType* x2,
                                                      PointType (*W)[9], int size)
{
    for (int i = 0; i < size; i++)
    {
        PointType H2 = x2Array[levelHigh[index[i]]] - x2Array[index[i]];
        PointType H1 = x1Array[index[i] + 1] - x1Array[index[i]];

        PointType X1 = (x1Array[index[i]] - x1[i]) / H1;

        PointType wx1_1 = PointType(0.5 * (0.5 + X1) * (0.5 + X1));
        PointType wx11  = PointType(0.5 * (0.5 - X1) * (0.5 - X1));
        PointType wx10  = PointType(1 - wx1_1 - wx11);
        PointType wx10T = PointType(0.75 - X1 * X1);

        PointType X2 = (x2Array[index[i]] - x2[i]) / H2;

        PointType wx2_1 = PointType(0.5 * (0.5 + X2) * (0.5 + X2));
        PointType wx21  = PointType(0.5 * (0.5 - X2) * (0.5 - X2));
        PointType wx20  = PointType(1 - wx2_1 - wx21);
        PointType wx20T = PointType(0.75 - X2 * X2);

        W[i][0] = wx2_1 * wx1_1;
        W[i][1] = wx2_1 * wx10;
        W[i][2] = wx2_1 * wx11;
        W[i][3] = wx20 * wx1_1;
        W[i][4] = wx20 * wx10;
        W[i][5] = wx20 * wx11;
        W[i][6] = wx21 * wx1_1;
        W[i][7] = wx21 * wx10;
        W[i][8] = wx21 * wx11;
    }
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate1(int index, PointType x1, PointType x2, PointType* W)
{
    PointType H2 = x2Array[points[index][2]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    PointType XX = (x2Array[index] - x2) / H2;

    PointType wx21;
    PointType wx22;

    PointType wx11;
    PointType wx12;

    PointType xx1 = std::abs(x1 - x1Array[index]) / H1;
    PointType xx2 = std::abs(x2 - x2Array[index]) / H2;

    if (index == points[index][0])
    {
        wx12 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
        wx11 = 1 - wx12;

        wx22 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
        wx22 =
            PointType(0.5 * (0.5 + xx2) * (0.5 + xx2)) + PointType(0.5 * (0.5 - xx2) * (0.5 - xx2));
        wx21 = 1 - wx22;
    }

    if (index == points[index][1])
    {
        wx11 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
        wx12 = 1 - wx11;

        wx22 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
        wx22 =
            PointType(0.5 * (0.5 + xx2) * (0.5 + xx2)) + PointType(0.5 * (0.5 - xx2) * (0.5 - xx2));
        wx21 = 1 - wx22;
    }

    if (index == points[index][2])
    {
        wx12 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
        wx11 = 1 - wx12;

        wx21 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
        wx22 = 1 - wx21;
    }

    if (index == points[index][3])
    {
        wx11 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
        wx12 = 1 - wx11;

        wx21 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
        wx22 = 1 - wx21;
    }

    W[0] = wx11 * wx21;
    W[1] = wx12 * wx21;
    W[2] = wx11 * wx22;
    W[3] = wx12 * wx22;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate2(int index, PointType x1, PointType x2, PointType* W)
{
    PointType H2 = x2Array[points[index][3]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    PointType xx2 = std::abs(x2Array[index] - x2) / H2;

    PointType wx21;
    PointType wx22;

    if (index == points[index][1])
    {
        //	wx22 = 8.0 * xx2*xx2*xx2 - 8.0*xx2*xx2 + 3 * xx2;
        wx22 =
            PointType(0.5 * (0.5 + xx2) * (0.5 + xx2)) + PointType(0.5 * (0.5 - xx2) * (0.5 - xx2));
        wx21 = 1 - wx22;
    }

    if (index == points[index][4])
    {
        wx21 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
        wx22 = 1 - wx21;
    }

    PointType X = (x1Array[points[index][1]] - x1) / H1;

    PointType wx1_1 = PointType(0.5 * (0.5 + X) * (0.5 + X));
    PointType wx11  = PointType(0.5 * (0.5 - X) * (0.5 - X));
    PointType wx10  = PointType(1 - wx1_1 - wx11);

    W[0] = wx21 * wx1_1;
    W[1] = wx21 * wx10;
    W[2] = wx21 * wx11;
    W[3] = wx22 * wx1_1;
    W[4] = wx22 * wx10;
    W[5] = wx22 * wx11;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate3(int index, PointType x1, PointType x2, PointType* W)
{
    PointType H2 = x2Array[points[index][2]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    PointType xx1 = std::abs(x1 - x1Array[index]) / H1;

    PointType wx11;
    PointType wx12;

    if (index == points[index][3])
    {
        wx11 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
        wx12 = 1 - wx11;
    }

    if (index == points[index][2])
    {
        wx12 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
        wx11 = 1 - wx12;
    }

    PointType X = (x2Array[points[index][2]] - x2) / H2;

    PointType wx2_1 = PointType(0.5 * (0.5 + X) * (0.5 + X));
    PointType wx21  = PointType(0.5 * (0.5 - X) * (0.5 - X));
    PointType wx20  = PointType(1 - wx2_1 - wx21);
    PointType wx20T = PointType(0.75 - X * X);

    W[0] = wx11 * wx2_1;
    W[1] = wx12 * wx2_1;
    W[2] = wx11 * wx20;
    W[3] = wx12 * wx20;
    W[4] = wx11 * wx21;
    W[5] = wx12 * wx21;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::Wcalculate4(int index, PointType x1, PointType x2, PointType* W)
{
    PointType H2 = x2Array[points[index][2]] - x2Array[points[index][0]];
    PointType H1 = x1Array[points[index][1]] - x1Array[points[index][0]];

    PointType xx1 = std::abs(x1 - x1Array[index]) / H1;
    PointType xx2 = std::abs(x2 - x2Array[index]) / H2;

    PointType wx21;
    PointType wx22;

    PointType wx11;
    PointType wx12;

    if (points[index][0] != -1)
    {
        if (x1 >= x1Array[points[index][0]] && x1 <= x1Array[points[index][1]] &&
            x2 >= x2Array[points[index][0]] && x2 <= x2Array[points[index][3]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
        {
            wx11 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
            wx12 = 1 - wx11;

            wx21 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
            wx22 = 1 - wx21;

            W[0] = wx11 * wx21;
            W[1] = wx12 * wx21;
            W[3] = wx11 * wx22;
            W[4] = wx12 * wx22;

            W[2] = 0;
            W[5] = 0;
            W[6] = 0;
            W[7] = 0;
            W[8] = 0;
        };
    }

    if (points[index][2] != -2)
    {
        if (x1 >= x1Array[points[index][1]] && x1 <= x1Array[points[index][2]] &&
            x2 >= x2Array[points[index][0]] && x2 <= x2Array[points[index][4]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
        {
            wx12 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
            wx11 = 1 - wx12;

            wx21 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
            wx22 = 1 - wx21;

            W[1] = wx11 * wx21;
            W[2] = wx12 * wx21;
            W[4] = wx11 * wx22;
            W[5] = wx12 * wx22;

            W[0] = 0;
            W[3] = 0;
            W[6] = 0;
            W[7] = 0;
            W[8] = 0;
        }
    }

    if (points[index][6] != -1)
    {
        if (x1 >= x1Array[points[index][0]] && x1 <= x1Array[points[index][1]] &&
            x2 >= x2Array[points[index][3]] && x2 <= x2Array[points[index][6]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
        {
            wx11 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
            wx12 = 1 - wx11;

            wx22 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
            wx21 = 1 - wx22;

            W[3] = wx11 * wx21;
            W[4] = wx12 * wx21;
            W[6] = wx11 * wx22;
            W[7] = wx12 * wx22;

            W[0] = 0;
            W[1] = 0;
            W[2] = 0;
            W[5] = 0;
            W[8] = 0;
        }
    }

    if (points[index][8] != -1)
    {
        if (x1 >= x1Array[points[index][1]] && x1 <= x1Array[points[index][2]] &&
            x2 >= x2Array[points[index][3]] && x2 <= x2Array[points[index][6]] &&
            std::abs(x1 - x1Array[index]) < H1 / 2 && std::abs(x2 - x2Array[index]) < H2 / 2)
        {
            wx12 = 8.0 * xx1 * xx1 * xx1 - 8.0 * xx1 * xx1 + 3 * xx1;
            wx11 = 1 - wx12;

            wx22 = 8.0 * xx2 * xx2 * xx2 - 8.0 * xx2 * xx2 + 3 * xx2;
            wx21 = 1 - wx22;

            W[4] = wx11 * wx21;
            W[5] = wx12 * wx21;
            W[6] = wx11 * wx22;
            W[8] = wx12 * wx22;

            W[0] = 0;
            W[1] = 0;
            W[2] = 0;
            W[3] = 0;
            W[7] = 0;
        }
    }
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ValueInterpolate(int index, const PointType* W,
                                                     const std::vector<PointType>& ValueArray,
                                                     PointType&                    result)
{
    switch (cellType[index])
    {
    case 0:
        return ValueInterpolate0(index, W, ValueArray, result);
    case 1:
        return ValueInterpolate1(index, W, ValueArray, result);
    case 2:
        return ValueInterpolate2(index, W, ValueArray, result);
    case 3:
        return ValueInterpolate3(index, W, ValueArray, result);
    case 4:
        return ValueInterpolate4(index, W, ValueArray, result);
    }
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ValueInterpolate0(int index, const PointType* W,
                                                      const std::vector<PointType>& ValueArray,
                                                      PointType&                    result)
{
    result = W[0] * ValueArray[levelLow[index] - 1] + W[1] * ValueArray[levelLow[index]] +
             W[2] * ValueArray[levelLow[index] + 1] + W[3] * ValueArray[index - 1] +
             W[4] * ValueArray[index] + W[5] * ValueArray[index + 1] +
             W[6] * ValueArray[levelHigh[index] - 1] + W[7] * ValueArray[levelHigh[index]] +
             W[8] * ValueArray[levelHigh[index] + 1];
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ValueInterpolate1(int index, const PointType* W,
                                                      const std::vector<PointType>& ValueArray,
                                                      PointType&                    result)
{
    result = 0;
    for (int i = 0; i < 4; i++)
        result = result + W[i] * ValueArray[points[index][i]];
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ValueInterpolate2(int index, const PointType* W,
                                                      const std::vector<PointType>& ValueArray,
                                                      PointType&                    result)
{
    result = 0;
    for (int i = 0; i < 6; i++)
        result = result + W[i] * ValueArray[points[index][i]];
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ValueInterpolate3(int index, const PointType* W,
                                                      const std::vector<PointType>& ValueArray,
                                                      PointType&                    result)
{
    result = 0;
    for (int i = 0; i < 6; i++)
        result = result + W[i] * ValueArray[points[index][i]];
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ValueInterpolate4(int index, const PointType* W,
                                                      const std::vector<PointType>& ValueArray,
                                                      PointType&                    result)
{
    result = 0;
    for (int i = 0; i < 9; i++)
    {
        if (points[index][i] != -1)
            result = result + W[i] * ValueArray[points[index][i]];
    }
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate(int index, const PointType* W,
                                                    PointType particleCharge, PointType* rho)
{
    switch (cellType[index])
    {
    case 0:
        return ChargeCalculate0(index, W, particleCharge, rho);
    case 1:
        return ChargeCalculate1(index, W, particleCharge, rho);
    case 2:
        return ChargeCalculate2(index, W, particleCharge, rho);
    case 3:
        return ChargeCalculate3(index, W, particleCharge, rho);
    case 4:
        return ChargeCalculate4(index, W, particleCharge, rho);
    };
    return;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate0(int index, const PointType* W,
                                                     PointType particleCharge, PointType* rho)
{
    rho[levelLow[index] - 1]  = rho[levelLow[index] - 1] + W[0] * particleCharge;
    rho[levelLow[index]]      = rho[levelLow[index]] + W[1] * particleCharge;
    rho[levelLow[index] + 1]  = rho[levelLow[index] + 1] + W[2] * particleCharge;
    rho[index - 1]            = rho[index - 1] + W[3] * particleCharge;
    rho[index]                = rho[index] + W[4] * particleCharge;
    rho[index + 1]            = rho[index + 1] + W[5] * particleCharge;
    rho[levelHigh[index] - 1] = rho[levelHigh[index] - 1] + W[6] * particleCharge;
    rho[levelHigh[index]]     = rho[levelHigh[index]] + W[7] * particleCharge;
    rho[levelHigh[index] + 1] = rho[levelHigh[index] + 1] + W[8] * particleCharge;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate1(int index, const PointType* W,
                                                     PointType particleCharge, PointType* rho)
{
    for (int i                = 0; i < 4; i++)
        rho[points[index][i]] = rho[points[index][i]] + W[i] * particleCharge;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate2(int index, const PointType* W,
                                                     PointType particleCharge, PointType* rho)
{
    for (int i                = 0; i < 6; i++)
        rho[points[index][i]] = rho[points[index][i]] + W[i] * particleCharge;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate3(int index, const PointType* W,
                                                     PointType particleCharge, PointType* rho)
{
    for (int i                = 0; i < 6; i++)
        rho[points[index][i]] = rho[points[index][i]] + W[i] * particleCharge;
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeCalculate4(int index, const PointType* W,
                                                     PointType particleCharge, PointType* rho)
{
}

template <class PointType>
void ParticleShape2dTSC<PointType>::ChargeZeroing(){

}

template <class PointType>
PointType ParticleShape2dTSC<PointType>::CellArea()
{
    return 0;
}

template <class PointType>
PointType ParticleShape2dTSC<PointType>::CellVolume(int index)
{
    return cellVolume[index];
}

template class ParticleShape2dTSC<float>;
template class ParticleShape2dTSC<double>;