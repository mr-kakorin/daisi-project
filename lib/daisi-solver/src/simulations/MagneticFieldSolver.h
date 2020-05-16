#ifndef MAGNETICFIELDSOLVER_H
#define MAGNETICFIELDSOLVER_H
#include <memory>
#include <vector>

template <class PointType>
class GridData2d;
template <class PointType>
class GridData2dpolar;
template <class PointType>
class GridData3d;
template <class PointType>
class GridData2daxs;
template <class PointType>
class ElectrodeCurrent;
template <class PointType>
class MagneticFieldSolver
{
  public:
    void FieldSimulate(std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
                       const std::shared_ptr<GridData2daxs<PointType>>& gridData, PointType dt);
    void FieldSimulate(std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
                       const std::shared_ptr<GridData2d<PointType>>& gridData, PointType dt);
    void FieldSimulate(std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
                       const std::shared_ptr<GridData2dpolar<PointType>>& gridData, PointType dt);
    void FieldSimulate(std::vector<std::shared_ptr<ElectrodeCurrent<PointType>>>& electrodes,
                       const std::shared_ptr<GridData3d<PointType>>& gridData, PointType dt){

    };
};
#endif