#include "geomTools.h"
#include "BoundaryContainer2d.h"
#include "Geom.h"

template <class PointType>
void mergeSortResize(int nElements, std::vector<int> list,
                     std::vector<std::shared_ptr<BoundaryContainer2d<PointType>>>& boundaries,
                     std::vector<DGeo::Edge<PointType>>& result, std::string& errorMsg)
{

    std::vector<DGeo::Edge<PointType>> edges = boundaries[list[0]]->EdgesData;

    for (int i = 1; i < list.size(); i++)
    {
        for (int j = 0; j < boundaries[list[i]]->EdgesData.size(); j++)
            edges.push_back(boundaries[list[i]]->EdgesData[j]);
    };

    double length = 0;

    for (int i = 0; i < edges.size(); i++)
        length = length + edges[i].length();

    double avLength = length / nElements;

    DGeo::Edge<PointType>  lastEdge;
    int                    flag1;
    int                    flag2;
    DGeo::Point<PointType> nextPoint      = edges[0].point2;
    int                    lastEdgeNumber = 0;
    int                    nonConnected   = 0;
    lastEdge                              = edges[0];
    int flagisLast                        = 0;
    for (int i = 0; i < edges.size(); i++)
    {
        flag1      = 0;
        flag2      = 0;
        double tol = edges[i].length() * 1e-3;

        for (int j = 0; j < edges.size(); j++)
        {
            if (i != j)
            {
                if (edges[i].point1.isEqual(edges[j].point1, tol) ||
                    edges[i].point1.isEqual(edges[j].point2, tol))
                    flag1++;
                if (edges[i].point2.isEqual(edges[j].point1, tol) ||
                    edges[i].point2.isEqual(edges[j].point2, tol))
                    flag2++;
            }
        }
        if (flag1 == 0 && flag2 != 0)
        {
            if (flagisLast == 0)
            {
                nextPoint      = edges[i].point2;
                lastEdge       = edges[i];
                lastEdgeNumber = i;
                flagisLast     = 1;
            }
            nonConnected++;
        }
        if (flag2 == 0 && flag1 != 0)
        {
            if (flagisLast == 0)
            {
                nextPoint       = edges[i].point1;
                lastEdge.point2 = edges[i].point1;
                lastEdge.point1 = edges[i].point2;
                flagisLast      = 1;
                lastEdgeNumber  = i;
            }
            nonConnected++;
        }
    }
    if (nonConnected > 2)
    {
        errorMsg = "Boundary is not continuous!";
        return;
    };

    std::vector<DGeo::Edge<PointType>> resArr;
    int                                resize;

    for (int k = 0; k < edges.size(); k++)
    {
        resize     = ceil(lastEdge.length() / avLength);
        resArr     = lastEdge.resize(resize);
        double tol = lastEdge.length() * 1e-3;
        for (int i = 0; i < resArr.size(); i++)
        {
            if (result.size() && (result.back().point2.Dist2Point(resArr[i].point1) >
                                  result.back().point2.Dist2Point(resArr[i].point2)))
            {
                resArr[i].SwapPoints();
            };
            result.push_back(resArr[i]);
        }

        flag1 = 0;
        flag2 = 0;
        for (int i = 0; i < edges.size(); i++)
        {
            if (i != lastEdgeNumber)
            {
                if (nextPoint.isEqual(edges[i].point1, tol))
                    flag2++;
                if (nextPoint.isEqual(edges[i].point2, tol))
                    flag1++;
            };
            if (flag1 == 0 && flag2 != 0)
            {
                nextPoint      = edges[i].point2;
                lastEdge       = edges[i];
                lastEdgeNumber = i;
                break;
            }
            if (flag2 == 0 && flag1 != 0)
            {
                nextPoint       = edges[i].point1;
                lastEdge.point2 = edges[i].point1;
                lastEdge.point1 = edges[i].point2;
                lastEdgeNumber  = i;
                break;
            }
        }
        if (flag2 == 0 && flag1 == 0)
            break;
    }

    for (int i = 0; i < result.size() - 1; i++)
    {
        if (result[i].IsEqual(result[i + 1], 1e-12))
            result.erase(result.begin() + i);
    }
}

template void
mergeSortResize<float>(int nElements, std::vector<int> list,
                       std::vector<std::shared_ptr<BoundaryContainer2d<float>>>& boundaries,
                       std::vector<DGeo::Edge<float>>& result, std::string& errorMsg);

template void
mergeSortResize<double>(int nElements, std::vector<int> list,
                        std::vector<std::shared_ptr<BoundaryContainer2d<double>>>& boundaries,
                        std::vector<DGeo::Edge<double>>& result, std::string& errorMsg);