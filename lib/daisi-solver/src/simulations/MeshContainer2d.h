#ifndef MESHCONTAINER2D_H
#define MESHCONTAINER2D_H
#include "Dmath.h"
#include "VTKIncludeSolver.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <iostream>
#include <vector>
namespace DGeo
{
template <class PointType>
class Edge;
template <class PointType>
class Point;
};
namespace Dmath
{
class imat;
};
template <class PointType>
class GridData2daxs;
template <class PointType>
class GridData2d;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;
template <class PointType>
class GridData3dpolar;

template <class PointType>
class MeshGenerator;
template <class PointType>
class MeshContainer2d
{
    friend class boost::serialization::access;
    friend MeshGenerator<PointType>;

  private:
    std::vector<std::vector<int>> cellContainer;

    int flagTemplNumb;

    vtkUnstructuredGrid* VTKgrid;

    int    nVertexX;
    int    nVertexY;
    size_t nPoints;
    size_t nCell;
    void   CreateCells();
    // vtkSmartPointer<vtkUnstructuredGrid> VTKgrid;

    BOOST_SERIALIZATION_SPLIT_MEMBER()
    template <class Archive>
    void save(Archive& ar, const unsigned int) const;
    template <class Archive>
    void load(Archive& ar, const unsigned int);

  public:
    vtkUnstructuredGrid* GetVTKGrid(int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
                                    const std::shared_ptr<GridData2d<PointType>>& gridData,
                                    std::string                                   name);
    vtkUnstructuredGrid* GetVTKGrid(int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
                                    const std::shared_ptr<GridData2daxs<PointType>>& gridData,
                                    std::string                                      name);
    vtkUnstructuredGrid* GetVTKGrid(int flag, double param, vtkSmartPointer<vtkFloatArray>& vtkData,
                                    const std::shared_ptr<GridData2dpolar<PointType>>& gridData,
                                    std::string                                        name);
    void NullVtkGrid();
    void clear();
    MeshContainer2d();
    vtkUnstructuredGrid* GetVTKGrid();
    int                  getnVertexX() const;
    int                  getnVertexY() const;
    PointType            xmin;
    PointType            xmax;
    PointType            ymin;
    PointType            ymax;

    Dmath::imat templNumb;
    Dmath::imat flagMatrix;

    std::vector<int>              linearToTemplNumb;
    std::vector<std::vector<int>> neighbourContainer;

    PointType h1;
    PointType h2;

    std::vector<DGeo::Point<PointType>> serialMeshData;

    std::vector<std::vector<DGeo::Point<PointType>>> meshData;

    DGeo::Point<int> doStepX(DGeo::Point<int> mpoint, int step);
    DGeo::Point<int> doStepY(DGeo::Point<int> mpoint, int step);
    DGeo::Point<PointType> getNextX(int& x, int& y, bool& isEnd);
    DGeo::Point<PointType> getNextY(int& x, int& y, bool& isEnd);

    void Extrude(double dz);

    void VtkWrite(std::string);
    void Convert2GridData(const std::shared_ptr<GridData2d<PointType>>&);
    void Convert2GridData(const std::shared_ptr<GridData2daxs<PointType>>&);
    void Convert2GridData(const std::shared_ptr<GridData2dpolar<PointType>>&);

    void                         ConvertMesh2VTKUnstructuredGrid();
    vtkSmartPointer<vtkPolyData> GetVTKBoundaryPoints();
};

#endif