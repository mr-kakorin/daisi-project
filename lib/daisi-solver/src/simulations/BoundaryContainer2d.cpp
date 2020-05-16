#include <exception>

#include <common_tools/constants.h>

#include "BoundaryContainer2d.h"
#include "Dmath.h"
#include "Geom.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <iostream>
DGeo::Edge<double>  EdgeTmp;
DGeo::Point<double> PointTmp;
void swap(float* a, float* b)
{
    float tmp = *a;
    *a        = *b;
    *b        = tmp;
};
void qs(float* s_arr, float* s1_arr, int first, int last)
{
    int   i = first, j = last;
    float x = s_arr[(first + last) / 2];

    do
    {
        while (s_arr[i] < x)
            i++;
        while (s_arr[j] > x)
            j--;

        if (i <= j)
        {
            if (s_arr[i] > s_arr[j])
            {
                swap(&s_arr[i], &s_arr[j]);
                swap(&s1_arr[i], &s1_arr[j]);
            }
            i++;
            j--;
        }
    } while (i <= j);

    if (i < last)
        qs(s_arr, s1_arr, i, last);
    if (first < j)
        qs(s_arr, s1_arr, first, j);
}

bool min(double a, double b)
{
    if (a <= b)
        return true;
    return false;
}

template class BoundaryContainer2d<float>;
template class BoundaryContainer2d<double>;

template <class PointType>
void BoundaryContainer2d<PointType>::RemoveFineEdges(float k)
{
    float l;
    for (int i = 0; i < ContainerSize - 1; i++)
    {
        l = EdgesData[i].length();
        for (int j = 0; j < ContainerSize; j++)
        {
            if (i > ContainerSize - 1)
                break;
            if (i != j)
            {
                if (EdgesData[i].point1.isEqual(EdgesData[j].point1, l * 1e-5))
                {
                    if (EdgesData[i].IsParallelFixEps(EdgesData[j], l * k))
                    {
                        EdgesData[i].point1 = EdgesData[j].point2;
                        EdgesData.erase(EdgesData.begin() + j);
                        ContainerSize--;
                        j--;
                        continue;
                    };
                }

                if (EdgesData[i].point1.isEqual(EdgesData[j].point2, l * 1e-5))
                {
                    if (EdgesData[i].IsParallelFixEps(EdgesData[j], l * k))
                    {
                        EdgesData[i].point1 = EdgesData[j].point1;
                        EdgesData.erase(EdgesData.begin() + j);
                        ContainerSize--;
                        j--;
                        continue;
                    };
                }

                if (EdgesData[i].point2.isEqual(EdgesData[j].point2, l * 1e-5))
                {
                    if (EdgesData[i].IsParallelFixEps(EdgesData[j], l * k))
                    {
                        EdgesData[i].point2 = EdgesData[j].point1;
                        EdgesData.erase(EdgesData.begin() + j);
                        ContainerSize--;
                        j--;
                        continue;
                    };
                }

                if (EdgesData[i].point2.isEqual(EdgesData[j].point1, l * 1e-5))
                {
                    if (EdgesData[i].IsParallelFixEps(EdgesData[j], l * k))
                    {
                        EdgesData[i].point2 = EdgesData[j].point2;
                        EdgesData.erase(EdgesData.begin() + j);
                        ContainerSize--;
                        j--;
                        continue;
                    };
                }
            };
        }
    };

    ConvertBoundary2VTKUnstructuredGrid();
};

template <class PointType>
void BoundaryContainer2d<PointType>::Merge(
    const std::shared_ptr<BoundaryContainer2d<PointType>> obj)
{

    for (int i = 0; i < obj->VertexData.size(); i++)
        VertexData.push_back(obj->VertexData[i]);

    for (int i = 0; i < obj->EdgesData.size(); i++)
        EdgesData.push_back(obj->EdgesData[i]);

    ContainerSize = ContainerSize + obj->ContainerSize;

    xmin = std::min(xmin, obj->xmin);
    ymin = std::min(ymin, obj->ymin);
    xmax = std::max(xmax, obj->xmax);
    ymax = std::max(ymax, obj->ymax);

    ConvertBoundary2VTKUnstructuredGrid();
};

template <class PointType>
void BoundaryContainer2d<PointType>::Reverse()
{
    BoundaryContainer2d<PointType> tmp = *this;
    int                            s   = EdgesData.size();
    for (int i = 0; i < s; i++)
    {
        EdgesData[i].point1 = tmp.EdgesData[s - 1 - i].point2;
        EdgesData[i].point2 = tmp.EdgesData[s - 1 - i].point1;
    };
}

/*template <class PointType>
int BoundaryContainer2d<PointType>::MergeWithSort(std::shared_ptr <BoundaryContainer2d<PointType>>
obj)
{

        ContainerSize = ContainerSize + obj->ContainerSize;
        double eps;
        DGeo::Edge<PointType> EdgeTmp;

        int flag = 0;

        while (obj->EdgesData.size()>0)
        {
                for (int i = 0; i < obj->EdgesData.size(); i++)
                {
                        if (obj->EdgesData.size() == 0)
                                break;

                        eps = obj->EdgesData[i].length()*1e-3;

                        if (EdgesData.back().point1.isEqual(obj->EdgesData[i].point1, eps) ||
EdgesData.back().point2.isEqual(obj->EdgesData[i].point2, eps))
                        {
                                EdgeTmp.point1 = obj->EdgesData[i].point2;
                                EdgeTmp.point2 = obj->EdgesData[i].point1;
                                EdgesData.push_back(EdgeTmp);
                                if (obj->EdgesData.size() == 1)
                                {
                                        obj->EdgesData.clear();
                                        break;
                                }
                                else
                                        obj->EdgesData.erase(obj->EdgesData.begin() + i);
                                flag = 1;
                                continue;
                        }

                        if (EdgesData.back().point1.isEqual(obj->EdgesData[i].point2, eps) ||
EdgesData.back().point2.isEqual(obj->EdgesData[i].point1, eps))
                        {
                                EdgesData.push_back(obj->EdgesData[i]);
                                if (obj->EdgesData.size() == 1)
                                {
                                        obj->EdgesData.clear();
                                        break;
                                }
                                else
                                        obj->EdgesData.erase(obj->EdgesData.begin() + i);
                                flag = 1;
                                continue;
                        }
                }
                if (!flag)
                        return flag;

                if (obj->EdgesData.size() == 0)
                        break;
        }

        return flag;

        xmin = std::min(xmin, obj->xmin);
        ymin = std::min(ymin, obj->ymin);
        xmax = std::max(xmax, obj->xmax);
        ymax = std::max(ymax, obj->ymax);

        ConvertBoundary2VTKUnstructuredGrid();

};*/

template <class PointType>
void BoundaryContainer2d<PointType>::sort()
{
    double eps;
    int    j;
    int    i;
    int    flag1 = 0;
    int    flag2 = 0;

    for (i = 0; i < EdgesData.size(); i++)
    {
        eps   = EdgesData[i].length() * 1e-3;
        flag1 = 0;
        flag2 = 0;
        for (j = 0; j < EdgesData.size(); j++)
        {
            if (i == j)
                continue;

            if (EdgesData[i].point1.isEqual(EdgesData[j].point1, eps) ||
                EdgesData[i].point1.isEqual(EdgesData[j].point2, eps))
            {
                flag1 = 1;
            }
            if (EdgesData[i].point2.isEqual(EdgesData[j].point1, eps) ||
                EdgesData[i].point2.isEqual(EdgesData[j].point2, eps))
            {
                flag2 = 1;
            }
        }
        if (0 == flag1 || 0 == flag2)
            break;
    };
    if (i == EdgesData.size())
        i == 0;

    BoundaryContainer2d<PointType> tmp = *this;

    EdgesData.clear();
    EdgesData.push_back(tmp.EdgesData[i]);
    tmp.EdgesData.erase(tmp.EdgesData.begin() + i);
    if (flag2 == 0)
    {
        DGeo::Point<PointType> t = EdgesData[0].point1;
        EdgesData[0].point1      = EdgesData[0].point2;
        EdgesData[0].point2      = t;
    };

    for (i = 0; i < ContainerSize - 1; i++)
    {
        for (j = 0; j < tmp.EdgesData.size(); j++)
        {

            if (EdgesData.back().point2.isEqual(tmp.EdgesData[j].point1, eps))
            {
                EdgesData.push_back(tmp.EdgesData[j]);
                tmp.EdgesData.erase(tmp.EdgesData.begin() + j);
                j--;
                continue;
            }
            if (EdgesData.back().point2.isEqual(tmp.EdgesData[j].point2, eps))
            {
                EdgesData.push_back(tmp.EdgesData[j]);

                DGeo::Point<PointType> t = tmp.EdgesData[j].point1;
                EdgesData.back().point1  = EdgesData.back().point2;
                EdgesData.back().point2  = t;

                tmp.EdgesData.erase(tmp.EdgesData.begin() + j);
                j--;
                continue;
            }
        }
    };

    tmp.ContainerSize = tmp.EdgesData.size();

    int ContainerSizeOld = ContainerSize;

    ContainerSize = EdgesData.size();

    while (ContainerSize != ContainerSizeOld)
    {
        tmp.sort();
        Merge(std::shared_ptr<BoundaryContainer2d<PointType>>(&tmp));
    };
};

template <class PointType>
void BoundaryContainer2d<PointType>::GetPlotXY(std::vector<float>& Xarr, std::vector<float>& Yarr)
{
    BoundaryContainer2d<PointType> tmp = *this;
    //	tmp.sort();

    /*for (int i = 0; i < ContainerSize; i++)
    {
    Xarr.push_back(EdgesData[i].point1.x);
    //Xarr.push_back(EdgesData[i].point2.x);

    Yarr.push_back(EdgesData[i].point1.y);
    //Yarr.push_back(EdgesData[i].point2.y);
    };
    Xarr.push_back(EdgesData.back().point2.x);
    Yarr.push_back(EdgesData.back().point2.y);*/

    for (int i = 0; i < ContainerSize; i++)
    {
        Xarr.push_back(tmp.EdgesData[i].point1.x);
        Yarr.push_back(tmp.EdgesData[i].point1.y);
    };
    Xarr.push_back(tmp.EdgesData.back().point2.x);
    Yarr.push_back(tmp.EdgesData.back().point2.y);

    //	qs(&Xarr[0], &Yarr[0], 0, Xarr.size()-1);
};
template <class PointType>
void BoundaryContainer2d<PointType>::GetPlotRotate(std::vector<std::vector<float>>& Xarr,
                                                   std::vector<std::vector<float>>& Yarr,
                                                   int                              angles)
{
    float                  dangle = 2 * commtools::PI() / angles;
    DGeo::Point<PointType> rotatePoint;
    rotatePoint.x = 0;
    rotatePoint.z = 0;
    rotatePoint.y = 0;
    DGeo::Point<PointType> point;

    std::vector<float> xtmp;
    std::vector<float> ytmp;
    for (int i = 0; i < ContainerSize; i++)
    {
        xtmp.clear();
        ytmp.clear();
        point   = EdgesData[i].point1;
        point.y = 0;

        for (int j = 0; j <= angles; j++)
        {
            xtmp.push_back(point.x);
            ytmp.push_back(point.y);
            point = point.rotate(dangle, rotatePoint);
        }
        Xarr.push_back(xtmp);
        Yarr.push_back(ytmp);

        xtmp.clear();
        ytmp.clear();
        point   = EdgesData[i].point2;
        point.y = 0;

        for (int j = 0; j <= angles; j++)
        {
            xtmp.push_back(point.x);
            ytmp.push_back(point.y);
            point = point.rotate(dangle, rotatePoint);
        }
        Xarr.push_back(xtmp);
        Yarr.push_back(ytmp);
    };
};

template <class PointType>
void BoundaryContainer2d<PointType>::GetPlotRotate1(std::vector<std::vector<float>>& Xarr,
                                                    std::vector<std::vector<float>>& Yarr,
                                                    int                              angles)
{
    float                  dangle = 2 * commtools::PI() / angles;
    DGeo::Point<PointType> rotatePoint;
    rotatePoint.x = 0;
    rotatePoint.z = 0;
    rotatePoint.y = 0;
    DGeo::Point<PointType> point1;
    DGeo::Point<PointType> point2;

    std::vector<float> xtmp;
    std::vector<float> ytmp;

    for (int j = 0; j <= angles; j++)
    {
        xtmp.clear();
        ytmp.clear();
        for (int i = 0; i < ContainerSize; i++)
        {

            point1 = EdgesData[i].point1;

            point2 = EdgesData[i].point2;

            point1 = point1.rotate(dangle * j, rotatePoint);
            xtmp.push_back(point1.x);
            ytmp.push_back(point1.y);

            point2 = point2.rotate(dangle * j, rotatePoint);
            xtmp.push_back(point2.x);
            ytmp.push_back(point2.y);
        };

        Xarr.push_back(xtmp);
        Yarr.push_back(ytmp);
    }

    for (int j = 0; j <= angles; j++)
    {
        xtmp.clear();
        ytmp.clear();
        for (int i = 0; i < ContainerSize; i++)
        {

            point1   = EdgesData[i].point1;
            point1.x = -point1.x;

            point2   = EdgesData[i].point2;
            point2.x = -point2.x;

            point1 = point1.rotate(dangle * j, rotatePoint);
            xtmp.push_back(point1.x);
            ytmp.push_back(point1.y);

            point2 = point2.rotate(dangle * j, rotatePoint);
            xtmp.push_back(point2.x);
            ytmp.push_back(point2.y);
        };

        Xarr.push_back(xtmp);
        Yarr.push_back(ytmp);
    }
};

template <class PointType>
BoundaryContainer2d<PointType>& BoundaryContainer2d<PointType>::
operator=(const BoundaryContainer2d<PointType>& right)
{
    if (this == &right)
    {
        return *this;
    }
    xmin          = right.xmin;
    xmax          = right.xmax;
    ymin          = right.ymin;
    ymax          = right.ymax;
    EdgesData     = right.EdgesData;
    VertexData    = right.VertexData;
    ContainerSize = right.ContainerSize;
    flagInit      = right.flagInit;
    if (flagInit != 0)
        ConvertBoundary2VTKUnstructuredGrid();
    return *this;
};
template <class PointType>
BoundaryContainer2d<PointType>::BoundaryContainer2d()
{
    VTKgrid       = NULL;
    flagInit      = 0;
    xmin          = 0;
    xmax          = 0;
    ymin          = 0;
    ymax          = 0;
    ContainerSize = 0;
};
template <class PointType>
vtkSmartPointer<vtkUnstructuredGrid>
BoundaryContainer2d<PointType>::GetBoundaryVTKUnstructuredGrid()
{
    return VTKgrid;
};
template <class PointType>
void BoundaryContainer2d<PointType>::ConvertBoundary2VTKUnstructuredGrid()
{
    VTKgrid                           = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < ContainerSize; i++)
    {
        points->InsertNextPoint(EdgesData[i].point1.x, EdgesData[i].point1.y,
                                EdgesData[i].point1.z);
        points->InsertNextPoint(EdgesData[i].point2.x, EdgesData[i].point2.y,
                                EdgesData[i].point2.z);
    };

    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();

    for (int i = 0; i < ContainerSize; i++)
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, 2 * i);
        line->GetPointIds()->SetId(1, 2 * i + 1);
        cellArray->InsertNextCell(line);
    };

    VTKgrid->SetPoints(points);
    VTKgrid->SetCells(VTK_LINE, cellArray);
};
template <class PointType>
void BoundaryContainer2d<PointType>::WriteBoundary2VTK(std::string InputFileName)
{
    FILE* fpOut = fopen(InputFileName.c_str(), "w");
    fprintf(fpOut,
            "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d "
            "double\n\n",
            2 * ContainerSize);
    for (int i = 0; i < ContainerSize; i++)
    {
        fprintf(fpOut, "%lf %lf %lf\n", EdgesData[i].point1.x, EdgesData[i].point1.y,
                EdgesData[i].point1.z);
        fprintf(fpOut, "%lf %lf %lf\n", EdgesData[i].point2.x, EdgesData[i].point2.y,
                EdgesData[i].point2.z);
    };
    int nc = ContainerSize * 3;
    fprintf(fpOut, "\nCELLS %d %d\n\n", ContainerSize, nc);
    for (int i = 0; i < ContainerSize; i++)
    {
        fprintf(fpOut, "2 %d %d\n", 2 * i, 2 * i + 1);
    };
    fprintf(fpOut, "\nCELL_TYPES %d \n\n", ContainerSize);
    for (int i = 0; i < ContainerSize; i++)
    {
        fprintf(fpOut, "3 \n");
    }
    fclose(fpOut);
};

template <class PointType>
int BoundaryContainer2d<PointType>::OnBoundary(double epsilon, DGeo::Point<PointType> Point)
{
    int val;
    for (int i = 0; i < ContainerSize; i++)
    {
        val = EdgesData[i].IsOnEdge(Point, epsilon);
        //	val = EdgesData[i].IsOnEdge(Point, std::min(double(epsilon), double(1e-3 *
        // EdgesData[i].length()))); 	val = EdgesData[i].IsOnEdge(Point, epsilon);

        if (val < 2)
            return val;
    };
    return 2; // �� �� �������
};
template <class PointType>
int BoundaryContainer2d<PointType>::InBoundary(double epsilon, DGeo::Point<PointType> StartPoint,
                                               PointType h)
{
    int numberIntersections = 0;

    std::vector<int> results;

    DGeo::Point<PointType> CurrentPoint = StartPoint;

    double xmaxloc = StartPoint.x;

    int val = OnBoundary(epsilon, StartPoint);

    if (val != 2)
        return val; //�� �������

    int flag = 1;
    while (xmaxloc <= xmax + epsilon)
    {
        // flag++;

        int val = OnBoundary(epsilon, StartPoint);
        if (val != 2)
        {
            results.push_back(0); //����� �� ������� - ���������
        }

        DGeo::Point<PointType> CurrentPoint = DGeo::Point<PointType>(StartPoint, h, 0, 0);

        DGeo::Edge<PointType> GridEdge;
        GridEdge.point1 = StartPoint;
        GridEdge.point2 = CurrentPoint;

        //	if (flag > 1)
        //	{
        flag                = 0;
        numberIntersections = numberIntersections + NumberIntersections(epsilon, GridEdge);
        //	}

        xmaxloc    = CurrentPoint.x;
        StartPoint = CurrentPoint;
    };
    if (numberIntersections % 2)
        return 2; //������
    return 3;     //�������
};

template <class PointType>
int BoundaryContainer2d<PointType>::NumberIntersections(double                       epsilon,
                                                        const DGeo::Edge<PointType>& edge)
{

    // double epsilonLoc = std::min(double(epsilon), double(edge.length()*1e-3));
    // double epsilonLoc = double(edge.length()*Dconst::epsFrac);
    int lastNumber          = -1;
    int numberIntersections = 0;
    for (int i = 0; i < ContainerSize; i++)
    {
        //	if (EdgesData[i].point1.y > edge.point1.y && EdgesData[i].point2.y > edge.point1.y)
        //		continue;

        //	if (EdgesData[i].point1.y < edge.point1.y && EdgesData[i].point2.y < edge.point1.y)
        //		continue;

        if (edge.IsIntersectionEdgeMesh(EdgesData[i], epsilon))
        {

            if (lastNumber != -1)
            {
                DGeo::Point<PointType> p1;
                DGeo::Point<PointType> p2;
                EdgesData[i].IsIntersectionEdge(edge, epsilon, &p1);
                EdgesData[lastNumber].IsIntersectionEdge(edge, epsilon, &p2);

                // if (p1.isEqual(p2, epsilonLoc))
                //{

                DGeo::Point<PointType> pp1 = EdgesData[lastNumber].Middle();
                DGeo::Point<PointType> pp2 = EdgesData[i].Middle();

                if (edge.IsOnOneSide(pp1, pp2, epsilon))
                    return 2;
                else
                    return 1;
                //	}
                //	else
                return 2;
            }
            lastNumber = i;
            numberIntersections++;
        }
    };
    return numberIntersections;
};
template <class PointType>
bool BoundaryContainer2d<PointType>::IsIntersection(double epsilon, DGeo::Edge<PointType> edge,
                                                    DGeo::Point<PointType>& intersectionPoint,
                                                    int& intersectionEdgeNumber) const
{
    DGeo::Point<PointType> point;

    for (int i = 0; i < ContainerSize; i++)
    {

        //	if (edge.IsIntersectionEdge(EdgesData[i], 1e-6, &point))
        //	if (edge.IsIntersectionEdgeLight(EdgesData[i], std::min(float(0.1*epsilon),
        // float(edge.length()*1e-6)), &point))//5e-3
        if (edge.IsIntersectionEdge(EdgesData[i], epsilon, &point))
        {
            intersectionEdgeNumber = i;
            intersectionPoint      = point;
            return true;
        }
    };
    return false;
};
template <class PointType>
bool BoundaryContainer2d<PointType>::IsIntersectionLight(double epsilon, DGeo::Edge<PointType> edge,
                                                         DGeo::Point<PointType>& intersectionPoint,
                                                         int& intersectionEdgeNumber) const
{
    DGeo::Point<PointType> point;
    // double eps = double(edge.length()*Dconst::epsFrac);

    for (int i = 0; i < ContainerSize; i++)
    {

        if (edge.IsIntersectionEdgeLight(EdgesData[i], epsilon, &point))
        {
            intersectionEdgeNumber = i;
            intersectionPoint      = point;
            return true;
        }
    };
    return false;
};
template <class PointType>
void BoundaryContainer2d<PointType>::FindAllIntersections(
    double epsilon, DGeo::Edge<PointType> edge,
    std::vector<DGeo::Point<PointType>>& intersectionPoint,
    std::vector<DGeo::Edge<PointType>>&  intersectionEdgesTmp) const
{
    intersectionPoint.clear();
    intersectionEdgesTmp.clear();
    DGeo::Point<PointType> point;
    for (int i = 0; i < ContainerSize; i++)
    {

        if (EdgesData[i].IsIntersectionEdgePoints(edge, epsilon, &point))
        {
            point.z = 0;
            intersectionPoint.push_back(point);
            intersectionEdgesTmp.push_back(EdgesData[i]);
        }
    };
};

template <class PointType>
void BoundaryContainer2d<PointType>::MaxVertex()
{
    int number = 0;
    for (int i = 1; i < ContainerSize; i++)
    {
        if (PointCmp(VertexData[i], VertexData[number], 0) < 0)
        {
            number = i;
        }
    };

    number = 0;
    for (int i = 1; i < ContainerSize; i++)
    {
        if (PointCmp(VertexData[i], VertexData[number], 0) > 0)
        {
            number = i;
        }
    };
};

template <class PointType>
void BoundaryContainer2d<PointType>::ttf(std::string& InputFile, std::string& errorMsg)
{
    FILE* fp = fopen(InputFile.c_str(), "r");
    char  ss[250];

    if (!fp)
    {
        errorMsg = "Invalid file path.";
        return;
        // this = NULL;
    };
    errorMsg = "";

    int nNodes;
    int nCells;
    fgets(ss, 250, fp);

    //std::ofstream myfile;
    //myfile.open ("/home/mrkakorin/diploma/Projects/razorblade_triode_unix_uni/mesh.txt", std::ios::app);
    //myfile << InputFile << std::endl;
    //myfile << ss <<std::endl;

    sscanf(ss, "%d %d", &nNodes, &nCells);
   //std::istringstream inp_str3(ss);
    //char c;
    //inp_str3 >> nNodes >> c >> nCells; // Eat the space and comma.
    //myfile << nNodes << ' ' << nCells <<std::endl;
    ContainerSize = nCells;
    double tmpx, tmpy, tmpz;
    int    num;
    fgets(ss, 250, fp);
    sscanf(ss, "%d %lf %lf %lf", &num, &tmpx, &tmpy, &tmpz);
   // std::istringstream inp_str(ss);
    //inp_str >> num >> c >> tmpx >> c >> tmpy >> c >> tmpz; // Eat the space and comma.
    // read amount of cells
    //myfile << num << ' '<< tmpx << ' '<< tmpy << ' '<< tmpz <<std::endl;

    // read index of point

    PointTmp.Number = num;
    // read x coordinate of point
    PointTmp.x = tmpx;

    xmax = tmpx;
    xmin = tmpx;
    // read y coordinate of point
    PointTmp.y = tmpy;

    ymax = tmpy;
    ymin = tmpy;

    // read z coordinate of point
    PointTmp.z = tmpz;

    VertexData.push_back(PointTmp);

    for (int i = 1; i < nNodes; i++)
    {
        fgets(ss, 250, fp);
        sscanf(ss, "%d %lf %lf %lf", &num, &tmpx, &tmpy, &tmpz);
        //std::istringstream inp_str1(ss);
        //inp_str1 >> num >> c >> tmpx >> c >> tmpy >> c >> tmpz;
        //myfile << num << ' '<< tmpx << ' '<< tmpy<< ' '<<tmpz<<std::endl;

        // read x coordinate of point
        PointTmp.Number = num;

        PointTmp.x = tmpx;

        if (tmpx > xmax)
            xmax = tmpx;
        if (tmpx < xmin)
            xmin = tmpx;
        // read y coordinate of point
        PointTmp.y = tmpy;

        if (tmpy > ymax)
            ymax = tmpy;
        if (tmpy < ymin)
            ymin = tmpy;
        // read z coordinate of point
        PointTmp.z = tmpz;

        VertexData.push_back(PointTmp);
    };

    int tmpI;
    for (int i = 0; i < nCells; i++)
    {
        // read index of edge
        // read number for reference
        int ind, nedge1, nedge2;
        fgets(ss, 250, fp);
        sscanf(ss, "%d %d %d %d", &num, &ind, &nedge1, &nedge2);
        //std::istringstream inp_str2(ss);
        //inp_str2 >> num >> c >> ind >> c >> nedge1 >> c >> nedge2;
        //myfile << num << ' '<< ind<< ' ' << nedge1 << ' '<< nedge2 << std::endl;

        EdgeTmp.point1 = findPoint(nedge1);
        PointTmp       = findPoint(nedge1);

        PointTmp       = findPoint(nedge2);
        EdgeTmp.point2 = findPoint(nedge2);

        EdgesData.push_back(EdgeTmp);
    };

    //	sort();
    MaxVertex();
    ConvertBoundary2VTKUnstructuredGrid();
    flagInit = 1;
    //myfile.close();
    fclose(fp);
}

template <class PointType>
BoundaryContainer2d<PointType>::BoundaryContainer2d(std::string& InputFile, std::string& errorMsg)
{
    ttf(InputFile, errorMsg);
};

template <class PointType>
BoundaryContainer2d<PointType>::BoundaryContainer2d(std::vector<std::vector<double>> vertex,
                                                    bool isTube, double roundingRadius)
{
    // ContainerSize = vertex[0].size() - 1;
    DGeo::Point<PointType>              PointTmp1;
    DGeo::Point<PointType>              PointTmp;
    DGeo::Edge<PointType>               EdgeTmp;
    std::vector<DGeo::Point<PointType>> VertexDataTmp1(3);
    std::vector<DGeo::Point<PointType>> VertexDataTmp2(3);
    int                                 num = 1;
    PointTmp1.x                             = vertex[0][0];
    xmax                                    = vertex[0][0];
    xmin                                    = vertex[0][0];
    PointTmp1.y                             = vertex[1][0];
    ymax                                    = vertex[1][0];
    ymin                                    = vertex[1][0];
    PointTmp1.z                             = 0;
    PointTmp1.Number                        = num;
    VertexData.push_back(PointTmp1);
    for (int i = 1; i < vertex[0].size(); i++)
    {
        num++;
        PointTmp.Number = num;
        PointTmp.x      = vertex[0][i];
        if (vertex[0][i] > xmax)
            xmax = vertex[0][i];
        if (vertex[0][i] < xmin)
            xmin   = vertex[0][i];
        PointTmp.y = vertex[1][i];
        if (vertex[1][i] > ymax)
            ymax = vertex[1][i];
        if (vertex[1][i] < ymin)
            ymin   = vertex[1][i];
        PointTmp.z = 0;
        VertexData.push_back(PointTmp);
    };
    if (isTube == true)
    {
        num = 2;
        VertexDataTmp1.assign(VertexData.begin(), VertexData.end() - 1);
        AngleRounding(VertexDataTmp1, roundingRadius, num);
        VertexDataTmp2.assign(VertexData.begin() + 1, VertexData.end());
        AngleRounding(VertexDataTmp2, roundingRadius, num);
        VertexData.erase(VertexData.begin() + 1, VertexData.end() - 1);
        VertexData.insert(VertexData.begin() + 1, VertexDataTmp1.begin(), VertexDataTmp1.end());
        VertexData.insert(VertexData.begin() + VertexDataTmp1.size() + 1, VertexDataTmp2.begin(),
                          VertexDataTmp2.end());
        VertexData[VertexDataTmp1.size() + VertexDataTmp1.size() + 1].Number = num;
    }
    EdgeTmp.point1 = VertexData[0];
    // num = 0;
    for (int i = 1; i < VertexData.size(); i++)
    {
        EdgeTmp.point2 = VertexData[i];
        EdgesData.push_back(EdgeTmp);
        EdgeTmp.point1 = VertexData[i];
    }
    // EdgeTmp.point2 = VertexData[0]; ��������� ������ � ��������� �����
    // EdgesData.push_back(EdgeTmp);
    ContainerSize = EdgesData.size();
    MaxVertex();
    ConvertBoundary2VTKUnstructuredGrid();
    flagInit = 1;
}
template <class PointType>
BoundaryContainer2d<PointType>::BoundaryContainer2d(std::vector<std::vector<double>> vertex,
                                                    bool IsFirst, bool IsLast,
                                                    double roundingRadius)
{
    int                                 TubeNumber = vertex[0].size() / 4;
    DGeo::Point<PointType>              PointTmp1;
    DGeo::Point<PointType>              PointTmp2;
    DGeo::Point<PointType>              PointTmp;
    DGeo::Edge<PointType>               EdgeTmp;
    std::vector<DGeo::Point<PointType>> VertexDataTmp1(3);
    std::vector<DGeo::Point<PointType>> VertexDataTmp2(3);
    std::vector<DGeo::Point<PointType>> VertexDataTmp3(3);
    std::vector<DGeo::Point<PointType>> VertexDataTmp4(3);
    int                                 currSize;
    int                                 num;
    int                                 num1 = 1;
    VertexData.resize(vertex[0].size());
    for (int TubeCount = 0; TubeCount < TubeNumber; TubeCount++)
    {
        num                 = TubeCount * 4 + 1;
        PointTmp1.x         = vertex[0][TubeCount * 4];
        xmax                = vertex[0][TubeCount * 4];
        xmin                = vertex[0][TubeCount * 4];
        PointTmp1.y         = vertex[1][TubeCount * 4];
        ymax                = vertex[1][TubeCount * 4];
        ymin                = vertex[1][TubeCount * 4];
        PointTmp1.z         = 0;
        PointTmp1.Number    = num;
        VertexData[num - 1] = PointTmp1;
        for (int i = 1; i < 4; i++)
        {
            num++;
            PointTmp.Number = num;
            PointTmp.x      = vertex[0][TubeCount * 4 + i];
            if (vertex[0][TubeCount * 4 + i] > xmax)
                xmax = vertex[0][TubeCount * 4 + i];
            if (vertex[0][TubeCount * 4 + i] < xmin)
                xmin   = vertex[0][TubeCount * 4 + i];
            PointTmp.y = vertex[1][TubeCount * 4 + i];
            if (vertex[1][TubeCount * 4 + i] > ymax)
                ymax = vertex[1][TubeCount * 4 + i];
            if (vertex[1][TubeCount * 4 + i] < ymin)
                ymin            = vertex[1][TubeCount * 4 + i];
            PointTmp.z          = 0;
            VertexData[num - 1] = PointTmp;
        };
        if ((TubeCount == 0 && IsFirst == 1) || (TubeCount == TubeNumber - 1 && IsLast == 1))
        {
            VertexData[TubeCount * 4].Number = num1;
            VertexDataTmp1 = {VertexData[TubeCount * 4 + 0], VertexData[TubeCount * 4 + 1],
                              VertexData[TubeCount * 4 + 2]};
            AngleRounding(VertexDataTmp1, roundingRadius, num1);
            VertexDataTmp2 = {VertexData[TubeCount * 4 + 1], VertexData[TubeCount * 4 + 2],
                              VertexData[TubeCount * 4 + 3]};
            AngleRounding(VertexDataTmp2, roundingRadius, num1);
            VertexData[TubeCount * 4 + 3].Number = num1;
            VertexData.push_back(VertexData[TubeCount * 4]);
            VertexData.insert(VertexData.end(), VertexDataTmp1.begin(), VertexDataTmp1.end());
            VertexData.insert(VertexData.end(), VertexDataTmp2.begin(), VertexDataTmp2.end());
            VertexData.push_back(VertexData[TubeCount * 4 + 3]);
            EdgeTmp.point1 = VertexData[TubeCount * 4];
            if (TubeCount == 0)
            {
                for (int i = TubeNumber * 4 + 1; i < VertexData.size(); i++)
                {
                    EdgeTmp.point2 = VertexData[i];
                    EdgesData.push_back(EdgeTmp);
                    EdgeTmp.point1 = VertexData[i];
                }
            }
            else
            {
                for (int i = TubeCount * 4 * VertexDataTmp1.size() + TubeNumber * 4 -
                             (2 * VertexDataTmp1.size() - 2) * IsFirst + 1;
                     i < VertexData.size(); i++)
                {
                    EdgeTmp.point2 = VertexData[i];
                    EdgesData.push_back(EdgeTmp);
                    EdgeTmp.point1 = VertexData[i];
                }
            }
        }
        else
        {
            VertexDataTmp1 = {VertexData[TubeCount * 4 + 0], VertexData[TubeCount * 4 + 1],
                              VertexData[TubeCount * 4 + 2]};
            AngleRounding(VertexDataTmp1, roundingRadius, num1);
            VertexDataTmp2 = {VertexData[TubeCount * 4 + 1], VertexData[TubeCount * 4 + 2],
                              VertexData[TubeCount * 4 + 3]};
            AngleRounding(VertexDataTmp2, roundingRadius, num1);
            VertexDataTmp3 = {VertexData[TubeCount * 4 + 2], VertexData[TubeCount * 4 + 3],
                              VertexData[TubeCount * 4 + 0]};
            AngleRounding(VertexDataTmp3, roundingRadius, num1);
            VertexDataTmp4 = {VertexData[TubeCount * 4 + 3], VertexData[TubeCount * 4 + 0],
                              VertexData[TubeCount * 4 + 1]};
            AngleRounding(VertexDataTmp4, roundingRadius, num1);

            // currSize = VertexData.size();
            // VertexData.resize(currSize + VertexDataTmp1.size() * 4);
            VertexData.insert(VertexData.end(), VertexDataTmp1.begin(), VertexDataTmp1.end());
            VertexData.insert(VertexData.end(), VertexDataTmp2.begin(), VertexDataTmp2.end());
            VertexData.insert(VertexData.end(), VertexDataTmp3.begin(), VertexDataTmp3.end());
            VertexData.insert(VertexData.end(), VertexDataTmp4.begin(), VertexDataTmp4.end());

            EdgeTmp.point1 = VertexData[TubeNumber * 4 + TubeCount * 4 * VertexDataTmp1.size() -
                                        (2 * VertexDataTmp1.size() - 2) * IsFirst];
            for (int i = TubeCount * 4 * VertexDataTmp1.size() + TubeNumber * 4 -
                         (2 * VertexDataTmp1.size() - 2) * IsFirst + 1;
                 i < VertexData.size(); i++)
            {
                EdgeTmp.point2 = VertexData[i];
                EdgesData.push_back(EdgeTmp);
                EdgeTmp.point1 = VertexData[i];
            }
            /*if ((TubeCount == 0 && IsFirst == 1) || (TubeCount == TubeNumber-1 && IsLast == 1))
            {}
            else*/
            //{
            EdgeTmp.point2 = VertexData[TubeNumber * 4 + TubeCount * 4 * VertexDataTmp1.size() -
                                        (2 * VertexDataTmp1.size() - 2) *
                                            IsFirst]; //��������� ������ � ��������� �����
            EdgesData.push_back(EdgeTmp);
        }
    }
    VertexData.erase(VertexData.begin(), VertexData.begin() + 4 * TubeNumber);
    ContainerSize = EdgesData.size();
    MaxVertex();
    ConvertBoundary2VTKUnstructuredGrid();
    flagInit = 1;
}

template <class PointType>
DGeo::Point<PointType> BoundaryContainer2d<PointType>::findPoint(int n)
{
    for (int i = 0; i < VertexData.size(); ++i)
    {
        if (VertexData[i].Number == n)
            return VertexData[i];
    }
};
template <class PointType>
void BoundaryContainer2d<PointType>::AddEdge(DGeo::Edge<PointType> edg)
{
    EdgesData.push_back(edg);
};
template <class PointType>
void BoundaryContainer2d<PointType>::AddVertex(DGeo::Point<PointType> ver)
{
    VertexData.push_back(ver);
};
template <class PointType>
void BoundaryContainer2d<PointType>::AngleRounding(
    std::vector<DGeo::Point<PointType>>& vertex, double roundingRadius,
    int& num) //������ ������������ 3 �������(2� �������) � ������ �����������
{
    bool Check = false;
    if (vertex.size() == 3)
        Check                            = true;
    DGeo::Point<PointType> PointTmp1     = vertex[0];
    DGeo::Point<PointType> PointTmpAngle = vertex[1];
    DGeo::Point<PointType> PointTmp3     = vertex[2];
    // DGeo::Point<PointType> PointTmp;
    DGeo::Point<PointType>              Center;
    DGeo::Point<PointType>              PointTmp1new;
    DGeo::Point<PointType>              PointTmp3new;
    int                                 numberOfPointsOnCurve = 10;
    std::vector<DGeo::Point<PointType>> VertexRound(numberOfPointsOnCurve + 1);
    bool                                Up;
    bool                                Left;
    if (PointTmpAngle.x == PointTmp1.x && PointTmpAngle.y > PointTmp1.y) // left up
    {
        PointTmp1new.z = 0;
        PointTmp3new.z = 0;
        Center.z       = 0;
        PointTmp1new.x = PointTmpAngle.x;
        PointTmp1new.y = PointTmpAngle.y - roundingRadius;
        PointTmp3new.x = PointTmpAngle.x + roundingRadius;
        PointTmp3new.y = PointTmpAngle.y;
        Center.x       = PointTmpAngle.x + roundingRadius;
        Center.y       = PointTmpAngle.y - roundingRadius;
        for (int i = 0; i < numberOfPointsOnCurve + 1; i++)
        {
            VertexRound[i].x =
                Center.x -
                roundingRadius * std::cos(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].y =
                Center.y +
                roundingRadius * std::sin(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].z      = 0;
            VertexRound[i].Number = num;
            num++;
        }
    }
    if (PointTmpAngle.x == PointTmp1.x && PointTmpAngle.y < PointTmp1.y) // right down
    {
        PointTmp1new.z = 0;
        PointTmp3new.z = 0;
        Center.z       = 0;
        PointTmp1new.x = PointTmpAngle.x;
        PointTmp1new.y = PointTmpAngle.y + roundingRadius;
        PointTmp3new.x = PointTmpAngle.x - roundingRadius;
        PointTmp3new.y = PointTmpAngle.y;
        Center.x       = PointTmpAngle.x - roundingRadius;
        Center.y       = PointTmpAngle.y + roundingRadius;
        for (int i = 0; i < numberOfPointsOnCurve + 1; i++)
        {
            VertexRound[i].x =
                Center.x +
                roundingRadius * std::cos(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].y =
                Center.y -
                roundingRadius * std::sin(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].z      = 0;
            VertexRound[i].Number = num;
            num++;
        }
    }
    if (PointTmpAngle.x > PointTmp1.x && PointTmpAngle.y == PointTmp1.y) // right up
    {
        PointTmp1new.z = 0;
        PointTmp3new.z = 0;
        Center.z       = 0;
        PointTmp1new.x = PointTmpAngle.x - roundingRadius;
        PointTmp1new.y = PointTmpAngle.y;
        PointTmp3new.x = PointTmpAngle.x;
        PointTmp3new.y = PointTmpAngle.y - roundingRadius;
        Center.x       = PointTmpAngle.x - roundingRadius;
        Center.y       = PointTmpAngle.y - roundingRadius;
        for (int i = 0; i < numberOfPointsOnCurve + 1; i++)
        {
            VertexRound[i].x =
                Center.x +
                roundingRadius * std::sin(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].y =
                Center.y +
                roundingRadius * std::cos(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].z      = 0;
            VertexRound[i].Number = num;
            num++;
        }
    }
    if (PointTmpAngle.x < PointTmp1.x && PointTmpAngle.y == PointTmp1.y) // left down
    {
        PointTmp1new.z = 0;
        PointTmp3new.z = 0;
        Center.z       = 0;
        PointTmp1new.x = PointTmpAngle.x + roundingRadius;
        PointTmp1new.y = PointTmpAngle.y;
        PointTmp3new.x = PointTmpAngle.x;
        PointTmp3new.y = PointTmpAngle.y + roundingRadius;
        Center.x       = PointTmpAngle.x + roundingRadius;
        Center.y       = PointTmpAngle.y + roundingRadius;
        for (int i = 0; i < numberOfPointsOnCurve + 1; i++)
        {
            VertexRound[i].x =
                Center.x -
                roundingRadius * std::sin(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].y =
                Center.y -
                roundingRadius * std::cos(i * commtools::PI() / 180 * (90 / numberOfPointsOnCurve));
            VertexRound[i].z      = 0;
            VertexRound[i].Number = num;
            num++;
        }
    }
    // VertexRound[0] = PointTmp1new;
    // VertexRound[numberOfPointsOnCurve + 1] = PointTmp3new;
    vertex.clear();
    vertex.resize(VertexRound.size());
    vertex = VertexRound;
}

template void BoundaryContainer2d<float>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void BoundaryContainer2d<double>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template void BoundaryContainer2d<double>::load<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive& ar, const unsigned int file_version);
template void BoundaryContainer2d<float>::save<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive& ar, const unsigned int file_version) const;

template <class PointType>
template <class Archive>
void BoundaryContainer2d<PointType>::save(Archive& ar, const unsigned int) const
{
    ar& EdgesData;
    ar& xmin;
    ar& xmax;
    ar& ymin;
    ar& ymax;
    ar& ContainerSize;
    ar& flagInit;
    ar& VertexData;
};
template <class PointType>
template <class Archive>
void BoundaryContainer2d<PointType>::load(Archive& ar, const unsigned int)
{
    ar& EdgesData;
    ar& xmin;
    ar& xmax;
    ar& ymin;
    ar& ymax;
    ar& ContainerSize;
    ar& flagInit;
    ar& VertexData;
    if (flagInit != 0)
        ConvertBoundary2VTKUnstructuredGrid();
};
