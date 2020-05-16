#include "BoundaryContainer2dpolar.h"
bool min(double a, double b)
{
    if (a <= b)
        return true;
    return false;
};
template BoundaryContainer2dpolar<float>;
template BoundaryContainer2dpolar<double>;

template <class PointType> void BoundaryContainer2dpolar<PointType>::RemoveFineEdges()
{
    float l;
    float k = 0.00005;
    for (int i = 0; i < ContainerSize; i++)
    {
        l = EdgesData[i].length();
        for (int j = 0; j < ContainerSize; j++)
        {
            if (i != j)
            {
                if (EdgesData[i].point1.isEqual(EdgesData[j].point1, l * 1e-5))
                {
                    if (EdgesData[i].IsParallel(EdgesData[j], l * k))
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
                    if (EdgesData[i].IsParallel(EdgesData[j], l * k))
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
                    if (EdgesData[i].IsParallel(EdgesData[j], l * k))
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
                    if (EdgesData[i].IsParallel(EdgesData[j], l * k))
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

template <class PointType> void BoundaryContainer2dpolar<PointType>::Merge(BoundaryContainer2d<PointType> obj)
{

    for (int i = 0; i < obj.VertexData.size(); i++)
        VertexData.push_back(obj.VertexData[i]);

    for (int i = 0; i < obj.EdgesData.size(); i++)
        EdgesData.push_back(obj.EdgesData[i]);

    ContainerSize = ContainerSize + obj.ContainerSize;

    xmin    = std::min(xmin, obj.xmin);
    ymin    = std::min(ymin, obj.ymin);
    xmax    = std::max(xmax, obj.xmax);
    ymax    = std::max(ymax, obj.ymax);
    epsilon = std::min(epsilon, obj.epsilon);

    ConvertBoundary2VTKUnstructuredGrid();
};

template <class PointType>
void BoundaryContainer2dpolar<PointType>::GetPlotXY(std::vector<float>& Xarr, std::vector<float>& Yarr)
{
    for (int i = 0; i < ContainerSize; i++)
    {
        Xarr.push_back(EdgesData[i].point1.x);
        Xarr.push_back(EdgesData[i].point2.x);

        Yarr.push_back(EdgesData[i].point1.y);
        Yarr.push_back(EdgesData[i].point2.y);
    };
};
template <class PointType>
void BoundaryContainer2dpolar<PointType>::GetPlotRotate(std::vector<std::vector<float>>& Xarr,
                                                        std::vector<std::vector<float>>& Yarr)
{
    float                  angles = 100;
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
BoundaryContainer2dpolar<PointType>& BoundaryContainer2d<PointType>::
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
    epsilon       = right.epsilon;
    flagInit      = right.flagInit;
    if (flagInit != 0)
        ConvertBoundary2VTKUnstructuredGrid();
    return *this;
};
template <class PointType> BoundaryContainer2dpolar<PointType>::BoundaryContainer2d()
{
    VTKgrid       = NULL;
    flagInit      = 0;
    xmin          = 0;
    xmax          = 0;
    ymin          = 0;
    ymax          = 0;
    epsilon       = 0;
    ContainerSize = 0;
};
template <class PointType>
vtkSmartPointer<vtkUnstructuredGrid> BoundaryContainer2dpolar<PointType>::GetBoundaryVTKUnstructuredGrid()
{
    return VTKgrid;
};
template <class PointType> void BoundaryContainer2dpolar<PointType>::ConvertBoundary2VTKUnstructuredGrid()
{
    VTKgrid                           = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    for (int i = 0; i < ContainerSize; i++)
    {
        points->InsertNextPoint(EdgesData[i].point1.x, EdgesData[i].point1.y, EdgesData[i].point1.z);
        points->InsertNextPoint(EdgesData[i].point2.x, EdgesData[i].point2.y, EdgesData[i].point2.z);
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
template <class PointType> void BoundaryContainer2dpolar<PointType>::WriteBoundary2VTK(std::string InputFileName)
{
    FILE* fpOut;
    fopen_s(&fpOut, InputFileName.c_str(), "w");
    fprintf_s(fpOut, "# vtk DataFile Version 3.0\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS %d double\n\n",
              2 * ContainerSize);
    for (int i = 0; i < ContainerSize; i++)
    {
        fprintf_s(fpOut, "%lf %lf %lf\n", EdgesData[i].point1.x, EdgesData[i].point1.y, EdgesData[i].point1.z);
        fprintf_s(fpOut, "%lf %lf %lf\n", EdgesData[i].point2.x, EdgesData[i].point2.y, EdgesData[i].point2.z);
    };
    int nc = ContainerSize * 3;
    fprintf_s(fpOut, "\nCELLS %d %d\n\n", ContainerSize, nc);
    for (int i = 0; i < ContainerSize; i++)
    {
        fprintf_s(fpOut, "2 %d %d\n", 2 * i, 2 * i + 1);
    };
    fprintf_s(fpOut, "\nCELL_TYPES %d \n\n", ContainerSize);
    for (int i = 0; i < ContainerSize; i++)
    {
        fprintf_s(fpOut, "3 \n");
    }
    fclose(fpOut);
};

template <class PointType> int BoundaryContainer2dpolar<PointType>::OnBoundary(DGeo::Point<PointType> Point)
{
    int val;
    for (int i = 0; i < ContainerSize; i++)
    {
        val = EdgesData[i].IsOnEdge(Point, epsilon);
        if (EdgesData[i].IsOnEdge(Point, epsilon) < 2)
            return val;
    };
    return 2; // �� �� �������
};
template <class PointType>
int BoundaryContainer2dpolar<PointType>::InBoundary(DGeo::Point<PointType> StartPoint, PointType h)
{
    int numberIntersections = 0;

    std::vector<int> results;

    DGeo::Point<PointType> CurrentPoint = StartPoint;

    double xmaxloc = StartPoint.x;

    int val = OnBoundary(StartPoint);

    if (val != 2)
        return val; //�� �������

    int flag = 1;
    while (xmaxloc <= xmax + epsilon)
    {
        // flag++;

        int val = OnBoundary(StartPoint);
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
        numberIntersections = numberIntersections + NumberIntersections(GridEdge);
        //	}

        xmaxloc    = CurrentPoint.x;
        StartPoint = CurrentPoint;
    };
    if (numberIntersections % 2)
        return 2; //������
    return 3;     //�������
};

template <class PointType> int BoundaryContainer2dpolar<PointType>::NumberIntersections(DGeo::Edge<PointType> edge)
{
    edge.point1.x -= 0 * epsilon * Dmath::sign(edge.point1.x);
    edge.point1.y += 0.0 * epsilon * Dmath::sign(edge.point1.y);
    edge.point2.x += 0.0000 * epsilon * Dmath::sign(edge.point2.x);
    edge.point2.y += 0.0 * epsilon * Dmath::sign(edge.point2.y);
    int numberIntersections = 0;
    int lastNumber          = -1;
    for (int i = 0; i < ContainerSize; i++)
    {

        if (edge.IsIntersectionEdge(EdgesData[i], epsilon))
        {

            if (lastNumber != -1)
            {
                DGeo::Point<PointType> p1;
                DGeo::Point<PointType> p2;
                edge.IsIntersectionEdge(EdgesData[i], epsilon, &p1);
                edge.IsIntersectionEdge(EdgesData[lastNumber], epsilon, &p2);

                // if (EdgesData[i].IsIntersectionEdge(EdgesData[lastNumber], epsilon))
                if (p1.isEqual(p2, epsilon))
                {

                    edge.point1.x -= edge.length() * 0.05 * Dmath::sign(edge.point1.x);
                    edge.point1.y -= edge.length() * 0.05 * Dmath::sign(edge.point1.y);
                    edge.point2.x += edge.length() * 0.05 * Dmath::sign(edge.point2.x);
                    edge.point2.y += edge.length() * 0.05 * Dmath::sign(edge.point2.y);
                    numberIntersections = NumberIntersections(edge);
                    return numberIntersections;
                }
                //���� 2 edge ���������� 2 ������������ edge - �� ��� 1 �����������
                //����������� 2� ������������� - 2 �����������
            }
            lastNumber = i;
            numberIntersections++;
        }
    };
    return numberIntersections;
};
template <class PointType>
bool BoundaryContainer2dpolar<PointType>::IsIntersection(DGeo::Edge<PointType>   edge,
                                                         DGeo::Point<PointType>& intersectionPoint,
                                                         int&                    intersectionEdgeNumber) const
{
    DGeo::Point<PointType> point;
    for (int i = 0; i < ContainerSize; i++)
    {

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
bool BoundaryContainer2dpolar<PointType>::IsIntersectionLight(DGeo::Edge<PointType>   edge,
                                                              DGeo::Point<PointType>& intersectionPoint,
                                                              int&                    intersectionEdgeNumber) const
{
    DGeo::Point<PointType> point;
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
void BoundaryContainer2dpolar<PointType>::FindAllIntersections(
    DGeo::Edge<PointType> edge, std::vector<DGeo::Point<PointType>>& intersectionPoint,
    std::vector<DGeo::Edge<PointType>>& intersectionEdgesTmp) const
{
    intersectionPoint.clear();
    intersectionEdgesTmp.clear();
    DGeo::Point<PointType> point;
    for (int i = 0; i < ContainerSize; i++)
    {

        if (edge.IsIntersectionEdgePoints(EdgesData[i], epsilon, &point))
        {
            point.z = 0;
            intersectionPoint.push_back(point);
            intersectionEdgesTmp.push_back(EdgesData[i]);
        }
    };
};

template <class PointType> void BoundaryContainer2dpolar<PointType>::MaxVertex()
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

template <class PointType> BoundaryContainer2dpolar<PointType>::BoundaryContainer2d(std::string InputFile)
{

    FILE* fp;
    char  ss[250];
    fopen_s(&fp, InputFile.c_str(), "r");
    int nNodes;
    int nCells;
    fgets(ss, 250, fp);
    sscanf_s(ss, "%d %d", &nNodes, &nCells);
    ContainerSize = nCells;
    double                 tmpx, tmpy, tmpz;
    int                    num;
    DGeo::Point<PointType> PointTmp;
    fgets(ss, 250, fp);
    sscanf_s(ss, "%d %lf %lf %lf", &num, &tmpx, &tmpy, &tmpz);
    // read amount of cells

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
        sscanf_s(ss, "%d %lf %lf %lf", &num, &tmpx, &tmpy, &tmpz);
        DGeo::Point<PointType> PointTmp;
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
        sscanf_s(ss, "%d %d %d %d", &num, &ind, &nedge1, &nedge2);

        DGeo::Edge<PointType>  EdgeTmp;
        DGeo::Point<PointType> PointTmp;

        EdgeTmp.point1 = findPoint(nedge1);
        PointTmp       = findPoint(nedge1);

        PointTmp       = findPoint(nedge2);
        EdgeTmp.point2 = findPoint(nedge2);

        EdgesData.push_back(EdgeTmp);
    };

    MaxVertex();
    ConvertBoundary2VTKUnstructuredGrid();
    flagInit = 1;
};

template <class PointType> DGeo::Point<PointType> BoundaryContainer2dpolar<PointType>::findPoint(int n)
{
    for (int i = 0; i < VertexData.size(); ++i)
    {
        if (VertexData[i].Number == n)
            return VertexData[i];
    }
};
