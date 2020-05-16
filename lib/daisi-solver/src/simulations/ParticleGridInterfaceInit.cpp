#include "ParticleGridInterface.h"
#include "BoundaryConditions.h"
#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "GridData.h"
#include "ParticleShape2dCIC.h"
#include "ParticleShape2dTSC.h"
#include "WcalculateVector.h"

template class ParticleGridInterface<float>;
template class ParticleGridInterface<double>;

bool isTemplateEqual(std::vector<int> templ1, std::vector<int> templ2)
{
    for (int i = 0; i < 8; i++)
    {
        if ((templ1[i] >= 0 && templ2[i] < 0) || (templ1[i] < 0 && templ2[i] >= 0))
            return false;
    }
    return true;
};

template void ParticleGridInterface<double>::init<GridData2daxs<double>>(
    int nflow, const std::shared_ptr<GridData2daxs<double>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<double>>& domainBoundary, int size, int numThreads,
    int blockSize);

template void ParticleGridInterface<float>::init<GridData2d<float>>(
    int nflow, const std::shared_ptr<GridData2d<float>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<float>>& domainBoundary, int size, int numThreads,
    int blockSize);

template void ParticleGridInterface<float>::init<GridData2daxs<float>>(
    int nflow, const std::shared_ptr<GridData2daxs<float>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<float>>& domainBoundary, int size, int numThreads,
    int blockSize);

template void ParticleGridInterface<float>::init<GridData2dpolar<float>>(
    int nflow, const std::shared_ptr<GridData2dpolar<float>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<float>>& domainBoundary, int size, int numThreads,
    int blockSize);

template void ParticleGridInterface<double>::init<GridData2d<double>>(
    int nflow, const std::shared_ptr<GridData2d<double>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<double>>& domainBoundary, int size, int numThreads,
    int blockSize);

template void ParticleGridInterface<double>::init<GridData2dpolar<double>>(
    int nflow, const std::shared_ptr<GridData2dpolar<double>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<double>>& domainBoundary, int size, int numThreads,
    int blockSize);

template <class PointType>
template <class gridDatatype>
void ParticleGridInterface<PointType>::init(
    int nflow, const std::shared_ptr<gridDatatype>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<PointType>>& domainBoundary, int size, int numThreads,
    int blockSize)
{

    gamma.resize(numThreads);
    startCellNumbers.resize(nflow);
    searchIndexes.resize(numThreads);
    searchIndexesLarge.resize(numThreads);
    searchIndexesLarge1.resize(numThreads);
    searchIndexesSpecial.resize(numThreads);

    tmpData1.resize(numThreads);
    tmpData2.resize(numThreads);
    // list.resize(numThreads);
    templNumb  = std::shared_ptr<Dmath::imat>(new Dmath::imat());
    *templNumb = templNumbIn;
    CICcellArray =
        std::shared_ptr<ParticleShape2dCIC<PointType>>(new ParticleShape2dCIC<PointType>());

    int nrow = templNumbIn.GetNrow();

    for (int thread = 0; thread < numThreads; thread++)
    {
        searchIndexesSpecial[thread].resize(templNumbIn.ArSize());

        for (int k                          = 0; k < templNumbIn.ArSize(); k++)
            searchIndexesSpecial[thread][k] = k;
    }

    for (int i = 0; i < numThreads; i++)
    {
        searchIndexes[i].resize(9);
        searchIndexesLarge[i].reserve(blockSize);
        searchIndexesLarge[i].resize(blockSize);

        searchIndexesLarge1[i].reserve(blockSize);
        tmpData1[i].reserve(blockSize);
        tmpData2[i].reserve(blockSize);
    }
    templNumb  = std::shared_ptr<Dmath::imat>(new Dmath::imat());
    *templNumb = templNumbIn;
    std::vector<std::vector<int>> templCICBoundary;
    std::vector<int>              vtmp;

    vtmp = {-1, -1, -1, -1, 0, 0, -1, 0, 0};
    templCICBoundary.push_back(vtmp);
    vtmp = {-1, -1, -1, 0, 0, -1, 0, 0, -1};
    templCICBoundary.push_back(vtmp);
    vtmp = {0, 0, -1, 0, 0, -1, -1, -1, -1};
    templCICBoundary.push_back(vtmp);
    vtmp = {-1, 0, 0, -1, 0, 0, -1, -1, -1};
    templCICBoundary.push_back(vtmp);

    std::vector<int> templTSC = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<std::vector<int>> templCICx2TSCx1Boundary;
    vtmp = {0, 0, 0, 0, 0, 0, -1, -1, -1};
    templCICx2TSCx1Boundary.push_back(vtmp);
    vtmp = {-1, -1, -1, 0, 0, 0, 0, 0, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {0, -1, -1, 0, 0, 0, 0, 0, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {-1, -1, 0, 0, 0, 0, 0, 0, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, 0, 0, 0, 0, -1, -1};
    templCICx2TSCx1Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, 0, 0, 0, -1, -1, 0};
    templCICx2TSCx1Boundary.push_back(vtmp);

    std::vector<std::vector<int>> templCICx1TSCx2Boundary;
    vtmp = {0, 0, -1, 0, 0, -1, 0, 0, -1};
    templCICx1TSCx2Boundary.push_back(vtmp);
    vtmp = {-1, 0, 0, -1, 0, 0, -1, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {0, 0, -1, 0, 0, -1, 0, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, 0, 0, -1, 0, 0, -1};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {-1, 0, 0, -1, 0, 0, 0, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    vtmp = {0, 0, 0, -1, 0, 0, -1, 0, 0};
    templCICx1TSCx2Boundary.push_back(vtmp);

    std::vector<std::vector<int>> templ3CellBoundary;
    vtmp = {0, 0, 0, 0, 0, 0, 0, 0, -1};
    templ3CellBoundary.push_back(vtmp);
    vtmp = {0, 0, 0, 0, 0, 0, -1, 0, 0};
    templ3CellBoundary.push_back(vtmp);
    vtmp = {-1, 0, 0, 0, 0, 0, 0, 0, 0};
    templ3CellBoundary.push_back(vtmp);
    vtmp = {0, 0, -1, 0, 0, 0, 0, 0, 0};
    templ3CellBoundary.push_back(vtmp);

    int n         = 0;
    int levelLow  = 0;
    int levelHigh = 0;

    cellArray.clear();
    CICArray.clear();

    int flagBreak;

    std::vector<int> flagOut;
    for (int j = 1; j < templNumb->GetNcol() - 1; j++)
    {

        for (int i = 1; i < templNumb->GetNrow() - 1; i++)
        {
            if ((*templNumb)(i, j) != -1)
                flagOut.push_back(flagMatrix(i, j));
        };
    };

    switch (priorityParticleShapeType)
    {
    case 0:
        CICcellArray->init(gridData->Getrho().size(), gridData->GetX1Pointer(),
                           gridData->GetX2Pointer());
        for (int j = 1; j < templNumb->GetNcol() - 1; j++)
        {

            for (int i = 1; i < templNumb->GetNrow() - 1; i++)
            {
                if ((*templNumb)(i, j) != -1)
                {
                    std::vector<int> v(9);
                    v[0] = (*templNumb)(i - 1, j - 1);
                    v[1] = (*templNumb)(i, j - 1);
                    v[2] = (*templNumb)(i + 1, j - 1);
                    v[3] = (*templNumb)(i - 1, j);
                    v[4] = (*templNumb)(i, j);
                    v[5] = (*templNumb)(i + 1, j);
                    v[6] = (*templNumb)(i - 1, j + 1);
                    v[7] = (*templNumb)(i, j + 1);
                    v[8] = (*templNumb)(i + 1, j + 1);

                    levelHigh = (*templNumb)(i, j + 1);
                    CICcellArray->AddCell(gridData, levelHigh, flagOut, v, domainBoundary,
                                          boundarypoints, 1);
                };
                n++;
            };
        };
        break;
    case 1:
        TSCcellArray->init(gridData->Getrho().size(), gridData->GetX1Pointer(),
                           gridData->GetX2Pointer());
        for (int j = 1; j < templNumb->GetNcol() - 1; j++)
        {

            for (int i = 1; i < templNumb->GetNrow() - 1; i++)
            {

                if ((*templNumb)(i, j) != -1)
                {
                    flagBreak = 0;

                    std::vector<int> v(9);
                    v[0]                  = (*templNumb)(i - 1, j - 1);
                    v[1]                  = (*templNumb)(i, j - 1);
                    v[2]                  = (*templNumb)(i + 1, j - 1);
                    v[3]                  = (*templNumb)(i - 1, j);
                    v[4]                  = (*templNumb)(i, j);
                    v[5]                  = (*templNumb)(i + 1, j);
                    v[6]                  = (*templNumb)(i - 1, j + 1);
                    v[7]                  = (*templNumb)(i, j + 1);
                    v[8]                  = (*templNumb)(i + 1, j + 1);
                    std::vector<int> vtmp = v;
                    int              flag = 1;

                    levelHigh = (*templNumb)(i, j + 1);
                    levelLow  = (*templNumb)(i, j - 1);

                    while (1)
                    {

                        if (isTemplateEqual(v, templTSC))
                        {
                            std::vector<int> t;
                            TSCcellArray->AddCell(gridData, levelLow, levelHigh, 0, t, flagOut,
                                                  domainBoundary, boundarypoints, 1);
                            // AddCell(gridData, 1, (*templNumb)(i, j), levelLow, levelHigh);
                            flag = 0;
                            break;
                        };

                        if (isTemplateEqual(v, templCICBoundary[0]) ||
                            isTemplateEqual(v, templCICBoundary[1]) ||
                            isTemplateEqual(v, templCICBoundary[2]) ||
                            isTemplateEqual(v, templCICBoundary[3]))
                        {
                            size_t n = v.size();
                            for (int k = 0; k < n; k++)
                            {
                                if (v[k] < 0)
                                {
                                    v.erase(v.begin() + k);
                                    k--;
                                    n = v.size();
                                };
                            }
                            TSCcellArray->AddCell(gridData, levelLow, levelHigh, 1, v, flagOut,
                                                  domainBoundary, boundarypoints, 1);
                            // AddCell(gridData, (*templNumb)(i, j), v[0], v[1], v[2], v[3]);
                            flag = 0;
                            break;
                        };

                        if (isTemplateEqual(v, templCICx2TSCx1Boundary[0]) ||
                            isTemplateEqual(v, templCICx2TSCx1Boundary[1]))
                        {
                            size_t n = v.size();
                            for (int k = 0; k < n; k++)
                            {
                                if (v[k] < 0)
                                {
                                    v.erase(v.begin() + k);
                                    k--;
                                    n = v.size();
                                };
                            }
                            TSCcellArray->AddCell(gridData, levelLow, levelHigh, 2, v, flagOut,
                                                  domainBoundary, boundarypoints, 1);
                            // AddCell(gridData, 1, (*templNumb)(i, j), v[0], v[1], v[2], v[3],
                            // v[4], v[5]);
                            flag = 0;
                            break;
                        };

                        for (int k = 2; k < 6; k++)
                        {
                            if (isTemplateEqual(v, templCICx2TSCx1Boundary[k]))
                            {
                                switch (k)
                                {
                                case 2:
                                    v.erase(v.begin() + 0);
                                    break;
                                case 3:
                                    v.erase(v.begin() + 2);
                                    break;
                                case 4:
                                    v.erase(v.begin() + 6);
                                    break;
                                case 5:
                                    v.erase(v.begin() + 8);
                                    break;
                                }
                                size_t n = v.size();
                                for (int k = 0; k < n; k++)
                                {
                                    if (v[k] < 0)
                                    {
                                        v.erase(v.begin() + k);
                                        k--;
                                        n = v.size();
                                    };
                                }
                                TSCcellArray->AddCell(gridData, levelLow, levelHigh, 2, v, flagOut,
                                                      domainBoundary, boundarypoints, 1);
                                //	AddCell(gridData, 1, (*templNumb)(i, j), v[0], v[1], v[2],
                                // v[3], v[4], v[5]);
                                flag      = 0;
                                flagBreak = 1;
                                break;
                            };
                        }
                        if (flagBreak)
                            break;

                        if (isTemplateEqual(v, templCICx1TSCx2Boundary[0]) ||
                            isTemplateEqual(v, templCICx1TSCx2Boundary[1]))
                        {
                            size_t n = v.size();
                            for (int k = 0; k < n; k++)
                            {
                                if (v[k] < 0)
                                {
                                    v.erase(v.begin() + k);
                                    k--;
                                    n = v.size();
                                };
                            }
                            TSCcellArray->AddCell(gridData, levelLow, levelHigh, 3, v, flagOut,
                                                  domainBoundary, boundarypoints, 1);
                            // AddCell(gridData, 0, (*templNumb)(i, j), v[0], v[1], v[2], v[3],
                            // v[4], v[5]);
                            flag = 0;
                            break;
                        };

                        for (int k = 2; k < 6; k++)
                        {
                            if (isTemplateEqual(v, templCICx1TSCx2Boundary[k]))
                            {
                                switch (k)
                                {
                                case 2:
                                    v.erase(v.begin() + 8);
                                    break;
                                case 3:
                                    v.erase(v.begin() + 2);
                                    break;
                                case 4:
                                    v.erase(v.begin() + 6);
                                    break;
                                case 5:
                                    v.erase(v.begin() + 0);
                                    break;
                                }
                                size_t n = v.size();
                                for (int k = 0; k < n; k++)
                                {
                                    if (v[k] < 0)
                                    {
                                        v.erase(v.begin() + k);
                                        k--;
                                        n = v.size();
                                    };
                                }
                                TSCcellArray->AddCell(gridData, levelLow, levelHigh, 3, v, flagOut,
                                                      domainBoundary, boundarypoints, 1);
                                // AddCell(gridData, 0, (*templNumb)(i, j), v[0], v[1], v[2], v[3],
                                // v[4], v[5]);
                                flag      = 0;
                                flagBreak = 1;
                                break;
                            };
                        }
                        if (flagBreak)
                            break;

                        if (isTemplateEqual(v, templ3CellBoundary[0]) ||
                            isTemplateEqual(v, templ3CellBoundary[1]) ||
                            isTemplateEqual(v, templ3CellBoundary[2]) ||
                            isTemplateEqual(v, templ3CellBoundary[3]))
                        {
                            TSCcellArray->AddCell(gridData, levelLow, levelHigh, 4, v, flagOut,
                                                  domainBoundary, boundarypoints, 1);
                            // AddCell(gridData, (*templNumb)(i, j), v);
                            flag = 0;
                            break;
                        };
                    }
                    if (flag == 0)
                    {
                        int flag1 = 0;
                        for (int s = 0; s < boundarypoints.size(); s++)
                        {
                            for (int s1 = 0; s1 < 9; s1++)
                            {
                                if (boundarypoints[s] == vtmp[s1])
                                    flag1 = 1;
                            }
                        }
                        if (flag1 == 1)
                            isBoundary.push_back(true);
                        else
                            isBoundary.push_back(false);
                    };
                };
                n++;
            };
        };
        break;
    }
    // particleGridInterface->AddCell(0, 0, 0, 0);
    // particleGridInterface->AddCell(1, 0, 0, 0);
    //	ParticleGridInterface->cellArray[0]->InCell(0, 0);
    //	ParticleGridInterface->cellArray[1]->InCell(0, 0);
};

template <class PointType>
void ParticleGridInterface<PointType>::init(
    int nflow, const std::shared_ptr<GridData3d<PointType>>& gridData, Dmath::imat templNumbIn,
    Dmath::imat flagMatrix, const std::vector<int>& boundarypoints,
    const std::shared_ptr<BoundaryContainer2d<PointType>>& domainBoundary, int size, int numThreads,
    int blockSize)
{
    startCellNumbers.resize(nflow);
    searchIndexes.resize(numThreads);
    searchIndexesLarge.resize(numThreads);
    searchIndexesLarge1.resize(numThreads);

    tmpData1.resize(numThreads);
    tmpData2.resize(numThreads);

    for (int i = 0; i < numThreads; i++)
    {
        searchIndexes[i].resize(9);
        searchIndexesLarge[i].reserve(blockSize);
        searchIndexesLarge1[i].reserve(blockSize);
        tmpData1.reserve(blockSize);
        tmpData2.reserve(blockSize);
    }

    templNumb  = std::shared_ptr<Dmath::imat>(new Dmath::imat());
    *templNumb = templNumbIn;
    CICcellArray =
        std::shared_ptr<ParticleShape2dCIC<PointType>>(new ParticleShape2dCIC<PointType>());
    std::vector<std::vector<int>> templCICBoundary;
    std::vector<int>              vtmp;

    int n         = 0;
    int levelLow  = 0;
    int levelHigh = 0;
    int levelZ    = 0;

    cellArray.clear();
    CICArray.clear();

    int flagBreak;

    std::vector<int> flagOut;
    for (int k = 1; k < templNumb->GetNz() - 1; k++)
    {
        for (int j = 1; j < templNumb->GetNcol() - 1; j++)
        {

            for (int i = 1; i < templNumb->GetNrow() - 1; i++)
            {

                if ((*templNumb)(i, j, k) != -1)
                    flagOut.push_back(flagMatrix(i, j, k));
            };
        };
    };
    CICcellArray->init(gridData->Getrho().size(), gridData->GetX1Pointer(),
                       gridData->GetX2Pointer(), gridData->GetX3Pointer());

    std::vector<int> v(8);

    for (int k = 1; k < templNumb->GetNz() - 1; k++)
    {
        for (int j = 1; j < templNumb->GetNcol() - 1; j++)
        {

            for (int i = 1; i < templNumb->GetNrow() - 1; i++)
            {

                if ((*templNumb)(i, j, k) != -1)
                {
                    v[0] = (*templNumb)(i, j, k);
                    v[1] = (*templNumb)(i, j + 1, k);
                    v[2] = (*templNumb)(i + 1, j, k);
                    v[3] = (*templNumb)(i, j, k + 1);
                    v[4] = (*templNumb)(i + 1, j + 1, k);
                    v[5] = (*templNumb)(i + 1, j, k + 1);
                    v[6] = (*templNumb)(i, j + 1, k + 1);
                    v[7] = (*templNumb)(i + 1, j + 1, k + 1);

                    levelHigh = (*templNumb)(i, j + 1, k);
                    levelZ    = (*templNumb)(i, j, k + 1);

                    if (levelZ == -1)
                    {
                        int tt = 0;
                    }

                    CICcellArray->AddCell(gridData, levelHigh, levelZ, flagOut, v, domainBoundary,
                                          boundarypoints, 1);
                }
            };
            n++;
        };
    };
};

template <class PointType>
void ParticleGridInterface<PointType>::AddCell(unsigned short particleShapeType,
                                               unsigned int basePoint, unsigned int levelLow,
                                               unsigned int levelHigh){
    /*IParticleShape2d<PointType> * newcell = new ParticleShapeCIC2d<PointType>(basePoint, levelLow,
    levelHigh);
    switch (particleShapeType)
    {
    case 0:
            break;
    case 1:
            newcell = new ParticleShapeTSC2d<PointType>(basePoint, levelLow, levelHigh); break;
    }
    cellArray.push_back(newcell);*/
};
template <class PointType>
void ParticleGridInterface<PointType>::AddCell(unsigned int basePoint, unsigned int p1,
                                               unsigned int p2, unsigned int p3, unsigned int p4){
    /*IParticleShape2d<PointType> * newcell = new ParticleShapeBoundaryCIC2d<PointType>(basePoint,
p1, p2, p3, p4);
    cellArray.push_back(newcell);
};
template <class PointType>
void ParticleGridInterface <PointType>::AddCell(unsigned short particleShapeType, unsigned int
basePoint, unsigned int
p1, unsigned int p2, unsigned int p3, unsigned int p4, unsigned int p5, unsigned int p6)
{
    /*IParticleShape2d<PointType> * newcell = new
ParticleShapeBoundaryCICx1TSCx2_2d<PointType>(basePoint, p1, p2, p3,
p4, p5, p6);

    switch (particleShapeType)
    {
    case 0:
            break;
    case 1:
            newcell = new ParticleShapeBoundaryCICx2TSCx1_2d<PointType>(basePoint, p1, p2, p3, p4,
p5, p6); break;
    }

    cellArray.push_back(newcell);*/
};
template <class PointType>
void ParticleGridInterface<PointType>::AddCell(unsigned int basePoint, const std::vector<int>& v){
    /*ParticleShapeBoundary3CellCIC2d<PointType> * newcell = new
    ParticleShapeBoundary3CellCIC2d<PointType>();

    if (v[0] == -1)
    {
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell1 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[1], v[2], v[4], v[5]); newcell->AddCell(*newSubcell1); ParticleShapeBoundaryCIC2d<PointType> *
    newSubcell2 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint, v[3], v[4], v[6], v[7]);
    newcell->AddCell(*newSubcell2);
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell3 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[4], v[5], v[7], v[8]); newcell->AddCell(*newSubcell3);
    }
    if (v[2] == -1)
    {
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell1 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[0], v[1], v[3], v[4]); newcell->AddCell(*newSubcell1); ParticleShapeBoundaryCIC2d<PointType> *
    newSubcell2 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint, v[3], v[4], v[6], v[7]);
    newcell->AddCell(*newSubcell2);
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell3 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[4], v[5], v[7], v[8]); newcell->AddCell(*newSubcell3);
    }
    if (v[6] == -1)
    {
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell1 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[0], v[1], v[3], v[4]); newcell->AddCell(*newSubcell1); ParticleShapeBoundaryCIC2d<PointType> *
    newSubcell2 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint, v[1], v[2], v[4], v[5]);
    newcell->AddCell(*newSubcell2);
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell3 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[4], v[5], v[7], v[8]); newcell->AddCell(*newSubcell3);
    }
    if (v[8] == -1)
    {
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell1 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[0], v[1], v[3], v[4]); newcell->AddCell(*newSubcell1); ParticleShapeBoundaryCIC2d<PointType> *
    newSubcell2 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint, v[1], v[2], v[4], v[5]);
    newcell->AddCell(*newSubcell2);
            ParticleShapeBoundaryCIC2d<PointType> * newSubcell3 = new
    ParticleShapeBoundaryCIC2d<PointType>(basePoint,
    v[3], v[4], v[6], v[7]); newcell->AddCell(*newSubcell3);
    }
    cellArray.push_back((IParticleShape2d<PointType>*)(newcell));*/
};
