#ifndef MESHCONTAINER3D_H
#define MESHCONTAINER3D_H
#include "Dmath.h"
#include "Geom.h"
#include "GridData.h"
#include "MeshContainer2d.h"
#include "ParticleGridInterface.h"
#include "VTKIncludeSolver.h"
#include <iostream>
#include <vector>
template <class PointType>
class MeshGenerator;
template <class PointType>
class MeshContainer3d
{
    friend MeshGenerator<PointType>;
    friend class boost::serialization::access;
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);
    void                 CreateCells();
    vtkUnstructuredGrid* VTKgrid;
    int                  flagTemplNumb;
    vtkUnstructuredGrid* resultPlane;

  public:
    vtkUnstructuredGrid* GetVTKGrid(int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
                                    const std::shared_ptr<GridData3d<PointType>>& gridData,
                                    std::string                                   name);
    std::vector<MeshContainer2d<PointType>> mesh;
    DGeo::Point<int> doStepX(DGeo::Point<int> mpoint, int step);

    DGeo::Point<int> doStepY(DGeo::Point<int> mpoint, int step);

    DGeo::Point<int> doStepZ(DGeo::Point<int> mpoint, int step);

    PointType xmin;
    PointType xmax;
    PointType ymin;
    PointType ymax;

    PointType h1;
    PointType h2;

    std::vector<std::vector<int>>       cellContainer;
    std::vector<std::vector<int>>       neighbourContainer;
    Dmath::imat                         templNumb;
    Dmath::imat                         flagMatrix;
    size_t                              nPoints;
    size_t                              nCell;
    std::vector<DGeo::Point<PointType>> serialMeshData;

    void Convert2GridData(const std::shared_ptr<GridData3d<PointType>>&);

    vtkUnstructuredGrid* GetVTKGrid();

    vtkSmartPointer<vtkPolyData> GetVTKBoundaryPoints();

    void clear();
    MeshContainer3d();
    void ConvertMesh2VTKUnstructuredGrid();
};

#endif