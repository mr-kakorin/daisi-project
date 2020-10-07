#include "MeshContainer3d.h"
#include "Dmath.h"
#include "FlagStringsSolver.h"
#include "Geom.h"
#include "GridData.h"
#include "ParticleGridInterface.h"


template <class PointType>
vtkUnstructuredGrid* MeshContainer3d<PointType>::GetVTKGrid(
    int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
    const std::shared_ptr<GridData3d<PointType>>& gridData, std::string name)
{
    if (resultPlane != NULL)
        resultPlane->Delete();

    resultPlane = vtkUnstructuredGrid::New();
    int       nearestIndex;
    PointType Mindist = 1e9;
    ;

    int       k1 = 0;
    int       k2 = templNumb.GetNz();
    int       i1 = 0;
    int       i2 = templNumb.GetNrow();
    int       j1 = 0;
    int       j2 = templNumb.GetNcol();
    PointType dist;

    Dmath::imat templNumbLoc;
    templNumbLoc.ones(i2, j2, k2);
    templNumbLoc = templNumbLoc * (-1);

    for (int k = k1; k < k2; k++)
    {
        for (int i = i1; i < i2; i++)
        {
            for (int j = j1; j < j2; j++)
            {
                if (templNumb(i, j, k) != -1)
                {
                    switch (flag)
                    {
                    case 0:
                        dist = std::abs(serialMeshData[templNumb(i, j, k)].x - param);
                        if (dist < Mindist)
                        {
                            Mindist      = dist;
                            nearestIndex = i;
                        }
                        break;
                    case 1:
                        dist = std::abs(serialMeshData[templNumb(i, j, k)].y - param);
                        if (dist < Mindist)
                        {
                            Mindist      = dist;
                            nearestIndex = j;
                        }
                        break;
                    case 2:
                        dist = std::abs(serialMeshData[templNumb(i, j, k)].z - param);
                        if (dist < Mindist)
                        {
                            Mindist      = dist;
                            nearestIndex = k;
                        }
                        break;
                    }
                }
            }
        }
    }

    switch (flag)
    {
    case 0:
        i1 = nearestIndex;
        i2 = nearestIndex + 1;
        break;
    case 1:
        j1 = nearestIndex;
        j2 = nearestIndex + 1;
        break;
    case 2:
        k1 = nearestIndex;
        k2 = nearestIndex + 1;
        break;
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    int n = 0;

    int index;

    for (int i = 0; i < flagStringsSolver::PlotFlags3d.size(); i++)
    {
        if (flagStringsSolver::PlotFlags3d[i] == name)
        {
            index = i;
            break;
        };
    }

    float valTmp;
    for (int k = k1; k < k2; k++)
    {
        for (int i = i1; i < i2; i++)
        {
            for (int j = j1; j < j2; j++)
            {
                if (templNumb(i, j, k) != -1)
                {
                    points->InsertNextPoint(serialMeshData[templNumb(i, j, k)].x,
                                            serialMeshData[templNumb(i, j, k)].y,
                                            serialMeshData[templNumb(i, j, k)].z);
                    switch (index)
                    {
                    case 0:
                        valTmp = gridData->Get_Ex()[templNumb(i, j, k)];
                        break;

                    case 1:
                        valTmp = gridData->Get_Ey()[templNumb(i, j, k)];
                        break;

                    case 2:
                        valTmp = gridData->Get_Ez()[templNumb(i, j, k)];
                        break;

                    case 3:
                        valTmp = sqrt(gridData->Get_Ez()[templNumb(i, j, k)] >
                                      gridData->Get_Ez()[templNumb(i, j, k)] +
                                          gridData->Get_Ex()[templNumb(i, j, k)] >
                                      gridData->Get_Ex()[templNumb(i, j, k)] +
                                          gridData->Get_Ey()[templNumb(i, j, k)] >
                                      gridData->Get_Ey()[templNumb(i, j, k)]);
                        break;

                    case 4:
                        valTmp = gridData->GetBx()[templNumb(i, j, k)];
                        break;

                    case 5:
                        valTmp = gridData->GetBy()[templNumb(i, j, k)];
                        break;

                    case 6:
                        valTmp = gridData->GetBz()[templNumb(i, j, k)];
                        break;

                    case 7:
                        valTmp = gridData->GetV()[templNumb(i, j, k)];
                        break;

                    case 8:
                        valTmp = gridData->Getrho()[templNumb(i, j, k)];
                        break;
                    };
                    vtkData->InsertNextValue(valTmp);
                    templNumbLoc(i, j, k) = n;
                    n++;
                }
            }
        }
    }

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    switch (flag)
    {
    case 0:
        j2 = j2 - 1;
        k2 = k2 - 1;
        break;
    case 1:
        i2 = i2 - 1;
        k2 = k2 - 1;
        break;
    case 2:
        j2 = j2 - 1;
        i2 = i2 - 1;
        break;
    }

    for (int k = k1; k < k2; k++)
    {
        for (int i = i1; i < i2; i++)
        {
            for (int j = j1; j < j2; j++)
            {
                if (templNumb(i, j, k) != -1)
                {

                    switch (flag)
                    {

                    case 0:
                        if (templNumb(i, j + 1, k) != -1 && templNumb(i, j + 1, k + 1) != -1 &&
                            templNumb(i, j, k + 1) != -1)
                        {
                            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                            quad->GetPointIds()->SetId(0, templNumbLoc(i, j, k));
                            quad->GetPointIds()->SetId(1, templNumbLoc(i, j + 1, k));
                            quad->GetPointIds()->SetId(3, templNumbLoc(i, j, k + 1));
                            quad->GetPointIds()->SetId(2, templNumbLoc(i, j + 1, k + 1));
                            cellArray->InsertNextCell(quad);
                        }
                        break;

                    case 1:
                        if (templNumb(i + 1, j, k) != -1 && templNumb(i + 1, j, k + 1) != -1 &&
                            templNumb(i, j, k + 1) != -1)
                        {
                            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                            quad->GetPointIds()->SetId(0, templNumbLoc(i, j, k));
                            quad->GetPointIds()->SetId(1, templNumbLoc(i + 1, j, k));
                            quad->GetPointIds()->SetId(3, templNumbLoc(i, j, k + 1));
                            quad->GetPointIds()->SetId(2, templNumbLoc(i + 1, j, k + 1));
                            cellArray->InsertNextCell(quad);
                        }
                        break;

                    case 2:
                        if (templNumb(i + 1, j, k) != -1 && templNumb(i + 1, j + 1, k) != -1 &&
                            templNumb(i, j + 1, k) != -1)
                        {
                            vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
                            quad->GetPointIds()->SetId(0, templNumbLoc(i, j, k));
                            quad->GetPointIds()->SetId(1, templNumbLoc(i + 1, j, k));
                            quad->GetPointIds()->SetId(3, templNumbLoc(i, j + 1, k));
                            quad->GetPointIds()->SetId(2, templNumbLoc(i + 1, j + 1, k));
                            cellArray->InsertNextCell(quad);
                        }
                        break;
                    }
                }
            }
        }
    }

    resultPlane->SetCells(VTK_QUAD, cellArray);
    resultPlane->SetPoints(points);

    return resultPlane;
}

template <class PointType>
void MeshContainer3d<PointType>::Convert2GridData(
    const std::shared_ptr<GridData3d<PointType>>& gridData){
    /*	if (mesh.size() == 0)
                    return;

            gridData->Getx().resize(serialMeshData.size());
            gridData->Gety().resize(serialMeshData.size());
            gridData->Getz().resize(serialMeshData.size());

            gridData->Get_Ex().resize(serialMeshData.size());;
            gridData->Get_Ey().resize(serialMeshData.size());
            gridData->Get_Ez().resize(serialMeshData.size());;

            gridData->GetBx().resize(serialMeshData.size());;
            gridData->GetBy().resize(serialMeshData.size());
            gridData->GetBz().resize(serialMeshData.size());;

            gridData->GetV().resize(serialMeshData.size());;
            gridData->Getrho().resize(serialMeshData.size());;

            gridData->CICArray.clear();
            gridData->flagOut.clear();

            int levelHigh = 0;
            int levelZ = 0;






            for (int k = 0; k < mesh.size() - 1; k++)
            {
                    for (int j = 0; j < mesh[0].getnVertexY() - 1; j++)
                    {
                            for (int i = 0; i < mesh[0].getnVertexX() - 1; i++)
                            {

                                    if (templNumb(i, j, k) != -1 && templNumb(i, j, k + 1) != -1)
                                    {
                                            if (templNumb(i + 1, j, k) != -1 && templNumb(i + 1, j +
       1, k) != -1 && templNumb(i, j + 1, k) != -1 && templNumb(i + 1, j, k + 1) != -1 &&
       templNumb(i + 1, j + 1, k + 1) != -1 && templNumb(i, j + 1, k + 1) != -1)
                                            {

                                                    int  p = templNumb(i, j, k);
                                                    levelHigh = templNumb(i, j + 1, k);
                                                    levelZ = templNumb(i, j, k + 1);

                                                    DGeo::Point<PointType> pp1 = serialMeshData[p];
                                                    DGeo::Point<PointType> pp2 =
       serialMeshData[levelHigh]; DGeo::Point<PointType> pp3 = serialMeshData[levelZ];


                                                    gridData->CICArray.push_back(levelHigh);
                                                    gridData->CICArrayZ.push_back(levelZ);

                                            }
                                            else
                                            {
                                                    gridData->CICArray.push_back(-1);
                                                    gridData->CICArrayZ.push_back(-1);

                                            };
                                    };
                            };
                    };
            }


            for (int i = 0; i < serialMeshData.size(); i++)
            {
                    gridData->Getx()[i] = serialMeshData[i].x;
                    gridData->Gety()[i] = serialMeshData[i].y;
                    gridData->Getz()[i] = serialMeshData[i].z;
            };*/
}

template <class PointType>
vtkUnstructuredGrid* MeshContainer3d<PointType>::GetVTKGrid()
{
    return VTKgrid;
}

template <class PointType>
vtkSmartPointer<vtkPolyData> MeshContainer3d<PointType>::GetVTKBoundaryPoints()
{
    vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
    return pointCloud;
}

template <class PointType>
void MeshContainer3d<PointType>::CreateCells()
{
    templNumb.ones(mesh[0].getnVertexX(), mesh[0].getnVertexY(), mesh.size() + 2);
    templNumb    = templNumb * (-1);
    int startPos = mesh[0].templNumb.ArSize();
    for (int k = 0; k < mesh.size(); k++)
    {
        templNumb.Set(mesh[k].templNumb + serialMeshData.size(), startPos);
        for (int i = 0; i < mesh[k].meshData.size(); i++)
        {
            for (int j = 0; j < mesh[k].meshData[i].size(); j++)
                serialMeshData.push_back(mesh[k].meshData[i][j]);
        }
        startPos = startPos + mesh[k].templNumb.ArSize();
    }

    nPoints = serialMeshData.size();
    std::vector<int> celltmp(8);
    std::vector<int> neighbour_tmp(8);

    int n = 0;

    for (int k = 0; k < templNumb.GetNz() - 1; k++)
    {
        for (int j = 0; j < templNumb.GetNcol() - 1; j++)
        {
            for (int i = 0; i < templNumb.GetNrow() - 1; i++)
            {

                if (templNumb(i, j, k) != -1 && templNumb(i, j, k + 1) != -1)
                {
                    n++;
                    if (templNumb(i + 1, j, k) != -1 && templNumb(i + 1, j + 1, k) != -1 &&
                        templNumb(i, j + 1, k) != -1 && templNumb(i + 1, j, k + 1) != -1 &&
                        templNumb(i + 1, j + 1, k + 1) != -1 && templNumb(i, j + 1, k + 1) != -1)
                    {
                        celltmp[0] = templNumb(i, j, k);
                        celltmp[1] = templNumb(i + 1, j, k);
                        celltmp[2] = templNumb(i + 1, j + 1, k);
                        celltmp[3] = templNumb(i, j + 1, k);

                        celltmp[4] = templNumb(i, j, k + 1);
                        celltmp[5] = templNumb(i + 1, j, k + 1);
                        celltmp[6] = templNumb(i + 1, j + 1, k + 1);
                        celltmp[7] = templNumb(i, j + 1, k + 1);
                        cellContainer.push_back(celltmp);
                    }
                }
            }
        }
    }

    nCell = cellContainer.size();
}

template <class PointType>
void MeshContainer3d<PointType>::ConvertMesh2VTKUnstructuredGrid()
{
    VTKgrid                                 = vtkUnstructuredGrid::New();
    vtkSmartPointer<vtkPoints>    points    = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    int ar[] = {0, 1, 2, 3, 4, 5, 6, 7};

    for (int i = 0; i < serialMeshData.size(); i++)
    {
        points->InsertNextPoint(serialMeshData[i].x, serialMeshData[i].y, serialMeshData[i].z);
    }

    for (int i = 0; i < nCell; i++)
    {
        vtkSmartPointer<vtkHexahedron> quad = vtkSmartPointer<vtkHexahedron>::New();
        for (int j = 0; j < 8; j++)
            quad->GetPointIds()->SetId(j, cellContainer[i][ar[j]]);
        cellArray->InsertNextCell(quad);
    };

    VTKgrid->SetPoints(points);
    VTKgrid->SetCells(VTK_HEXAHEDRON, cellArray);
}

template <class PointType>
DGeo::Point<int> MeshContainer3d<PointType>::doStepY(DGeo::Point<int> mpoint, int step)
{

    int i, j, k;
    i = mpoint.x;
    j = mpoint.y;
    k = mpoint.z;

    DGeo::Point<int> result;

    if ((mpoint.y == 0) && (step == -1))
    {
        return mpoint;
    }

    if ((mpoint.y == templNumb.GetNcol() - 1) && (step == 1))
    {
        return mpoint;
    }

    if (templNumb(i, j + step, k) != -1)
    {
        result        = mpoint;
        result.y      = mpoint.y + step;
        result.Number = templNumb(i, j + step, k);
        return result;
    }
    else
    {
        return mpoint;
    }
}

template <class PointType>
DGeo::Point<int> MeshContainer3d<PointType>::doStepX(DGeo::Point<int> mpoint, int step)
{

    DGeo::Point<int> result;

    int i, j, k;
    i = mpoint.x;
    j = mpoint.y;
    k = mpoint.z;

    if ((mpoint.x == 0) && (step == -1))
    {
        return mpoint;
    }

    if ((mpoint.x == templNumb.GetNrow() - 1) && (step == 1))
    {
        return mpoint;
    }

    if (templNumb(i + step, j, k) != -1)
    {
        result        = mpoint;
        result.x      = mpoint.x + step;
        result.Number = templNumb(i + step, j, k);
        return result;
    }
    else
    {
        return mpoint;
    }
}

template <class PointType>
void MeshContainer3d<PointType>::clear()
{
    for (int i = 0; i < mesh.size(); i++)
    {
        mesh[i].NullVtkGrid();
        mesh[i].clear();
    }
    neighbourContainer.clear();
    mesh.clear();
    serialMeshData.clear();
    cellContainer.clear();
    templNumb.clear();
    if (VTKgrid != NULL)
        VTKgrid->Delete();
    VTKgrid     = NULL;
    resultPlane = NULL;
}

template <class PointType>
MeshContainer3d<PointType>::MeshContainer3d()
{
    VTKgrid       = NULL;
    flagTemplNumb = 0;
    nPoints       = 0;
    nCell         = 0;
}

template <class PointType>
DGeo::Point<int> MeshContainer3d<PointType>::doStepZ(DGeo::Point<int> mpoint, int step)
{

    int i, j, k;
    i = mpoint.x;
    j = mpoint.y;
    k = mpoint.z;

    DGeo::Point<int> result;

    if ((mpoint.z == 0) && (step == -1))
    {
        return mpoint;
    }

    if ((mpoint.z == templNumb.GetNz() - 1) && (step == 1))
    {
        return mpoint;
    }

    if (templNumb(i, j, k + step) != -1)
    {
        result        = mpoint;
        result.z      = mpoint.z + step;
        result.Number = templNumb(i, j, k + step);
        return result;
    }
    else
    {
        return mpoint;
    }
}

template <class PointType>
template <class Archive>
void MeshContainer3d<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& mesh;
    ar& flagTemplNumb;
    ar& templNumb;
    ar& xmin;
    ar& xmax;
    ar& ymin;
    ar& ymax;
    ar& nPoints;
    ar& nCell;
    ar& h1;
    ar& h2;
    ar& flagMatrix;
}

template <class PointType>
template <class Archive>
void MeshContainer3d<PointType>::load(Archive& ar, const unsigned int)
{
    ar& mesh;
    ar& flagTemplNumb;
    ar& templNumb;
    ar& xmin;
    ar& xmax;
    ar& ymin;
    ar& ymax;
    ar& nPoints;
    ar& nCell;
    ar& h1;
    ar& h2;
    ar& flagMatrix;
    resultPlane = NULL;

    //	Dconst::epsFrac = h1;

    if (flagTemplNumb != 0)
    {
        CreateCells();
        ConvertMesh2VTKUnstructuredGrid();
    }
}

template class MeshContainer3d<float>;
template class MeshContainer3d<double>;

template void MeshContainer3d<float>::serialize<boost::archive::binary_iarchive>(
        boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void MeshContainer3d<double>::serialize<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version);

template void MeshContainer3d<double>::serialize<boost::archive::binary_iarchive>(
        boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void MeshContainer3d<float>::serialize<boost::archive::binary_oarchive>(
        boost::archive::binary_oarchive& ar, const unsigned int file_version);