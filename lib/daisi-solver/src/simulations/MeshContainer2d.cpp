#include "MeshContainer2d.h"
#include "Dmath.h"
#include "Geom.h"
#include "GridData.h"
#include "ParticleGridInterface.h"

template class MeshContainer2d<float>;
template class MeshContainer2d<double>;

template <class PointType>
DGeo::Point<PointType> MeshContainer2d<PointType>::getNextY(int& x, int& y, bool& isEnd)
{
    double x_ = meshData[y][x].x;
    while (y + 1 < meshData.size())
    {
        ++y;
        for (int i = 0; i < meshData[y].size(); ++i)
        {
            if (meshData[y][i].x == x_)
            {
                x = i;
                return meshData[y][x];
            }
        }
    }
    ++y;
    isEnd = true;
};

template <class PointType>
DGeo::Point<PointType> MeshContainer2d<PointType>::getNextX(int& x, int& y, bool& isEnd)
{

    if (x + 1 < meshData[y].size())
    {
        ++x;
        return meshData[y][x];
    }
    else
    {
        ++x;
        isEnd = true;
    }
};

template <class PointType>
vtkSmartPointer<vtkPolyData> MeshContainer2d<PointType>::GetVTKBoundaryPoints()
{
    vtkSmartPointer<vtkPoints>    points   = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType                     pid;

    for (int i = 0; i < meshData.size(); i++)
    {
        for (int j = 0; j < meshData[i].size(); j++)
        {
            if (meshData[i][j].isOut)
            {
                pid = points->InsertNextPoint(meshData[i][j].x, meshData[i][j].y, 0.0);
                vertices->InsertNextCell(1, &pid);
            }
        }
    }

    vtkSmartPointer<vtkPolyData> pointCloud = vtkSmartPointer<vtkPolyData>::New();
    pointCloud->SetPoints(points);
    pointCloud->SetVerts(vertices);
    return pointCloud;
}
template <class PointType>
void MeshContainer2d<PointType>::ConvertMesh2VTKUnstructuredGrid()
{
    VTKgrid                           = vtkUnstructuredGrid::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < nPoints; i++)
        points->InsertNextPoint(serialMeshData[i].x, serialMeshData[i].y, serialMeshData[i].z);

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    for (int i = 0; i < nCell; i++)
    {
        vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();
        quad->GetPointIds()->SetId(0, cellContainer[i][0]);
        quad->GetPointIds()->SetId(1, cellContainer[i][1]);
        quad->GetPointIds()->SetId(2, cellContainer[i][2]);
        quad->GetPointIds()->SetId(3, cellContainer[i][3]);
        cellArray->InsertNextCell(quad);
    };

    VTKgrid->SetPoints(points);
    VTKgrid->SetCells(VTK_QUAD, cellArray);
    //	VTKgrid->SetFieldData();
};

template <class PointType>
void MeshContainer2d<PointType>::Convert2GridData(
    const std::shared_ptr<GridData2d<PointType>>& gridData)
{
    gridData->Init(serialMeshData, 1);
    gridData->SetCells(flagMatrix, templNumb, nVertexX, nVertexY, 1);
};
template <class PointType>
void MeshContainer2d<PointType>::Convert2GridData(
    const std::shared_ptr<GridData2daxs<PointType>>& gridData)
{
    gridData->Init(serialMeshData, 2);
    gridData->SetCells(flagMatrix, templNumb, nVertexX, nVertexY, 2);
};

template <class PointType>
void MeshContainer2d<PointType>::Convert2GridData(
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData)
{
    gridData->Init(serialMeshData, 3);
    gridData->SetCells(flagMatrix, templNumb, nVertexX, nVertexY, 3);
};

template <class PointType>
void MeshContainer2d<PointType>::Extrude(double dz)
{
    for (int i = 0; i < meshData.size(); i++)
        for (int j           = 0; j < meshData[i].size(); j++)
            meshData[i][j].z = meshData[i][j].z + dz;

    for (int i              = 0; i < serialMeshData.size(); i++)
        serialMeshData[i].z = serialMeshData[i].z + dz;
}

template <class PointType>
void MeshContainer2d<PointType>::VtkWrite(std::string OutputFile)
{

    FILE* fpOut = fopen(OutputFile.c_str(), "w");
    fprintf(fpOut,
            "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d "
            "double\n\n",
            nPoints);
    for (int i = 0; i < nPoints; i++)
    {
        fprintf(fpOut, "%lf %lf %lf\n", serialMeshData[i].x, serialMeshData[i].y,
                serialMeshData[i].z);
    };
    int nc = int(nCell * 5);
    fprintf(fpOut, "\nCELLS %d %d\n\n", nCell, nc);
    for (int i = 0; i < nCell; i++)
    {
        fprintf(fpOut, "4 %d %d %d %d\n", cellContainer[i][0], cellContainer[i][1],
                cellContainer[i][2], cellContainer[i][3]);
    };
    fprintf(fpOut, "\nCELL_TYPES %d \n\n", nCell);
    for (int i = 0; i < nCell; i++)
    {
        fprintf(fpOut, "9 \n");
    }
    fclose(fpOut);
}
template <class PointType>
void MeshContainer2d<PointType>::CreateCells()
{

    for (int i = 0; i < meshData.size(); i++)
    {
        for (int j = 0; j < meshData[i].size(); j++)
            serialMeshData.push_back(meshData[i][j]);
    }
    nPoints = serialMeshData.size();
    std::vector<int> celltmp(4);
    std::vector<int> neighbour_tmp(4);

    int n = 0;
    for (int i = 0; i < nVertexX - 1; i++)
    {
        for (int j = 0; j < nVertexY - 1; j++)
        {
            if (templNumb(i, j) != -1)
            {
                n++;
                if (templNumb(i + 1, j) != -1 && templNumb(i + 1, j + 1) != -1 &&
                    templNumb(i, j + 1) != -1)
                {
                    celltmp[0] = templNumb(i, j);
                    celltmp[1] = templNumb(i, j + 1);
                    celltmp[2] = templNumb(i + 1, j + 1);
                    celltmp[3] = templNumb(i + 1, j);
                    cellContainer.push_back(celltmp);
                };
            };
        };
    };

    for (int i = 0; i < nVertexY; i++)
    {
        for (int j = 0; j < nVertexX; j++)
        {
            if (templNumb(j, i) != -1)
            {
                neighbour_tmp[0] = templNumb(j - 1, i);
                neighbour_tmp[1] = templNumb(j, i + 1);
                neighbour_tmp[2] = templNumb(j + 1, i);
                neighbour_tmp[3] = templNumb(j, i - 1);
                neighbourContainer.push_back(neighbour_tmp);
            }
        }
    }

    nCell = cellContainer.size();
};

template <class PointType>
DGeo::Point<int> MeshContainer2d<PointType>::doStepY(DGeo::Point<int> mpoint, int step)
{
    DGeo::Point<int> result;
    double           eps   = 1e-16;
    int              flagx = -1;
    if (step == 0)
        return mpoint;
    int x_ = mpoint.x;

    int node_number;
    if (step < 0)
    {
        node_number = neighbourContainer[meshData[mpoint.y][mpoint.x].Number][3];
    }
    else
    {
        node_number = neighbourContainer[meshData[mpoint.y][mpoint.x].Number][1];
    }

    if ((meshData.size() > (mpoint.y + step)) && ((mpoint.y + step) >= 0))
    {
        for (int i = 0; i < meshData[mpoint.y + step].size(); ++i)
        {
            if (std::abs(meshData[mpoint.y][x_].x - meshData[mpoint.y + step][i].x) < eps)
            {
                mpoint.x = i;
                flagx    = 1;
            }
        }
        if (flagx == 1)
        {
            if (meshData[mpoint.y + step][mpoint.x].Number == node_number)
            {
                result.x      = mpoint.x;
                result.y      = mpoint.y + step;
                result.Number = meshData[result.y][result.x].Number;
                return result;
            }
        }
    }

    return mpoint;
};

template <class PointType>
DGeo::Point<int> MeshContainer2d<PointType>::doStepX(DGeo::Point<int> mpoint, int step)
{
    DGeo::Point<int> result;

    if (step == 0)
        return mpoint;

    int node_number;
    if (step > 0)
    {
        node_number = neighbourContainer[meshData[mpoint.y][mpoint.x].Number][2];
    }
    else
    {
        node_number = neighbourContainer[meshData[mpoint.y][mpoint.x].Number][0];
    }

    if ((meshData[mpoint.y].size() > (mpoint.x + step)) && ((mpoint.x + step) >= 0))
    {
        if (meshData[mpoint.y][mpoint.x + step].Number == node_number)
        {
            result.x      = mpoint.x + step;
            result.y      = mpoint.y;
            result.Number = meshData[result.y][result.x].Number;
            return result;
        }
    }
    return mpoint;
};

template <class PointType>
vtkUnstructuredGrid* MeshContainer2d<PointType>::GetVTKGrid(
    int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
    const std::shared_ptr<GridData2d<PointType>>& gridData, std::string name)
{
    return VTKgrid;
};
template <class PointType>
vtkUnstructuredGrid* MeshContainer2d<PointType>::GetVTKGrid(
    int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
    const std::shared_ptr<GridData2daxs<PointType>>& gridData, std::string name)
{
    return VTKgrid;
};
template <class PointType>
vtkUnstructuredGrid* MeshContainer2d<PointType>::GetVTKGrid(
    int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
    const std::shared_ptr<GridData2dpolar<PointType>>& gridData, std::string name)
{
    return VTKgrid;
};
template <class PointType>
void MeshContainer2d<PointType>::NullVtkGrid()
{
    VTKgrid = NULL;
};
template <class PointType>
void MeshContainer2d<PointType>::clear()
{
    meshData.clear();
    serialMeshData.clear();
    cellContainer.clear();
    neighbourContainer.clear();
    templNumb.clear();
    flagMatrix.clear();
    linearToTemplNumb.clear();
    if (VTKgrid != NULL)
        VTKgrid->Delete();
    VTKgrid = NULL;
};
template <class PointType>
MeshContainer2d<PointType>::MeshContainer2d()
{
    VTKgrid       = NULL;
    flagTemplNumb = 0;
    nVertexX      = 0;
    nVertexY      = 0;
    nPoints       = 0;
    nCell         = 0;
    //	templNumb = std::shared_ptr<Dmath::imat> ( new Dmath::imat());
    //	flagMatrix = std::shared_ptr<Dmath::imat>(new Dmath::imat());
    // templNumb.zeros(5, 5);
};
template <class PointType>
vtkUnstructuredGrid* MeshContainer2d<PointType>::GetVTKGrid()
{
    return VTKgrid;
};
template <class PointType>
int MeshContainer2d<PointType>::getnVertexX() const
{
    return nVertexX;
}
template <class PointType>
int MeshContainer2d<PointType>::getnVertexY() const
{
    return nVertexY;
}

template void MeshContainer2d<float>::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void MeshContainer2d<double>::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version);

template void MeshContainer2d<double>::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void MeshContainer2d<float>::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version);

template <class PointType>
template <class Archive>
void MeshContainer2d<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& meshData;
    ar& nVertexX;
    ar& nVertexY;
    ar& flagTemplNumb;
    ar& templNumb;
    ar& h1;
    ar& h2;
    ar& linearToTemplNumb;
    ar& flagMatrix;

    ar& xmin;
    ar& xmax;
    ar& ymin;
    ar& ymax;
};
template <class PointType>
template <class Archive>
void MeshContainer2d<PointType>::load(Archive& ar, const unsigned int)
{
    ar& meshData;
    ar& nVertexX;
    ar& nVertexY;
    ar& flagTemplNumb;
    ar& templNumb;
    ar& h1;
    ar& h2;
    ar& linearToTemplNumb;
    ar& flagMatrix;

    //	Dconst::epsFrac = h1;

    ar& xmin;
    ar& xmax;
    ar& ymin;
    ar& ymax;

    if (flagTemplNumb != 0)
    {
        CreateCells();
        ConvertMesh2VTKUnstructuredGrid();
    }
};
