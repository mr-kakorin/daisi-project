#ifndef geomtool_h
#define geomtool_h

#include <memory>
#include <string>
#include <vector>

template <class PointType>
class BoundaryContainer2d;

namespace DGeo
{
template <class PointType>
class Edge;
}

template <class PointType>
void mergeSortResize(int nElements, std::vector<int> list,
                     std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                     std::vector<DGeo::Edge<PointType>>& result, std::string& errorMsg);

#endif