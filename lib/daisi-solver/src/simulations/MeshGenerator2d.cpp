#include "BoundaryContainer2d.h"
#include "BoundaryContainer3d.h"
#include "Geom.h"
#include "MeshContainer3d.h"
#include "MeshGenerator.h"
#include "ParticleGridInterface.h"
#include <omp.h>

MeshParams::MeshParams(std::string filename, double& eps, std::string& errorMsg)
{
  char  ss[250];
  FILE* fp = fopen(filename.c_str(), "r");

  if (!fp)
  {
    errorMsg = "Unable to open mesh.dat file";
    return;
    // this = NULL;
  }

  int i = -1;
  int j = 0;

  fgets(ss, 250, fp);
  int er = sscanf(ss, "%lf", &eps);
  if (!er)
  {
    errorMsg = "Tolerance is not defined";
    fclose(fp);
    return;
  }
  while (fgets(ss, 250, fp))
  {
    if (ss[0] == 'x')
      i++;
    else
    {
      std::vector<double> tmp(4);
      for (int ti = 0; ti < 4; ti++)
        tmp[ti] = NAN;
      int tmpi;
      int er = sscanf(ss, "%lf %lf %lf %lf", &tmp[0], &tmp[1], &tmp[2], &tmp[3]);
      for (int ti = 0; ti < 4; ti++)
      {
        if (std::isnan(tmp[ti]))
        {
          fclose(fp);
          errorMsg = "Incorrect input mesh.dat file";
          return;
        }
      }
      if (!er)
      {
        fclose(fp);
        errorMsg = "Incorrect input character in mesh.dat file";
        return;
      };
      leftBorder[i].push_back(tmp[0]);
      rightBorder[i].push_back(tmp[1]);
      params[i].push_back(std::vector<double>{tmp[2], tmp[3]});
      j++;
    }
  }
  if (i < 1)
  {
    errorMsg = "Incorrect input mesh description";
    fclose(fp);
    return;
  }
  fclose(fp);
}

double MeshParams::GetStep(double point, int flagX) const
{
  for (int i = 0; i < leftBorder[flagX].size(); i++)
  {
    if (point >= leftBorder[flagX][i] && point < rightBorder[flagX][i])
      return params[flagX][i][0] + (params[flagX][i][1] - params[flagX][i][0]) *
                                   (point - leftBorder[flagX][i]) /
                                   (rightBorder[flagX][i] - leftBorder[flagX][i]);
  }
}

int MeshParams::GetNumberOfSteps(int flagX, double min, double max, double& hMin) const
{
  double yCurLoc = min;
  hMin           = 1e18;
  int res        = 0;
  while (max > yCurLoc)
  {
    double h2 = GetStep(yCurLoc, flagX);
    if (h2 < hMin)
      hMin = h2;
    yCurLoc = yCurLoc + h2;
    res++;
  }
  return res;
}

template <class PointType>
void MeshGenerator<PointType>::MeshGenerate(
        std::string meshParamFile, double& progress,
        const std::shared_ptr<BoundaryContainer3d<PointType>> boundary,
        const std::shared_ptr<MeshContainer3d<PointType>>& mesh, int flagVTK, std::string& errorMsg){};

template <class PointType>
void MeshGenerator<PointType>::MeshGenerate(
        std::string meshParamFile, double& progress,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        const std::shared_ptr<MeshContainer3d<PointType>>& mesh, int flagVTK, std::string& errorMsg)
{
  mesh->clear();
  int    k = 0;
  double epsilon;
  MeshParams meshParams(meshParamFile, epsilon, errorMsg);

  if (errorMsg.size())
    return;

  mesh->xmin = boundary->xmin;
  mesh->xmax = boundary->xmax;
  mesh->ymin = boundary->ymin;
  mesh->ymax = boundary->ymax;

  mesh->flagTemplNumb = 1;

  mesh->mesh.resize(1);
  double tmp;
  mesh->mesh[0].nVertexX =
          meshParams.GetNumberOfSteps(0, boundary->xmin, boundary->xmax, tmp) + 4;
  mesh->h1 = tmp;
  mesh->h1 = meshParams.GetStep(boundary->xmin, 0);
  mesh->mesh[0].nVertexY =
          meshParams.GetNumberOfSteps(1, boundary->ymin, boundary->ymax, tmp) + 4;
  mesh->h2 = tmp;
  mesh->h2 = meshParams.GetStep(boundary->ymin, 1);

  const std::shared_ptr<MeshContainer2d<PointType>> tmp_mesh(&mesh->mesh[0]);

  MeshGenerate(meshParamFile, progress, boundary, tmp_mesh, flagVTK, errorMsg);

  double zCur = meshParams.leftBorder[2][0];
  double zMax = meshParams.rightBorder[2].back();

  int nZ = meshParams.GetNumberOfSteps(2, zCur, zMax, tmp) + 2;

  mesh->flagMatrix.ones(mesh->mesh[0].getnVertexX(), mesh->mesh[0].getnVertexY(), nZ);
  mesh->flagMatrix = mesh->flagMatrix * (-1);

  for (int i = 0; i < mesh->mesh[0].templNumb.GetNrow() - 1; i++)
  {
    for (int j = 0; j < mesh->mesh[0].templNumb.GetNcol() - 1; j++)
    {
      mesh->flagMatrix(i, j, 1) = mesh->mesh[0].flagMatrix(i, j);
    }
  }

  progress = 0;
  while (zCur < zMax)
  {
    progress = zCur / zMax;
    k++;

    for (int i = 0; i < mesh->mesh[0].templNumb.GetNrow() - 1; i++)
    {
      for (int j = 0; j < mesh->mesh[0].templNumb.GetNcol() - 1; j++)
      {
        mesh->flagMatrix(i, j, k + 1) = mesh->mesh[0].flagMatrix(i, j);
      }
    }

    mesh->mesh.resize(mesh->mesh.size() + 1);
    mesh->mesh[k] = mesh->mesh[k - 1];
    double hz     = meshParams.GetStep(zCur, 2);
    mesh->mesh[k].Extrude(hz);
    zCur = zCur + hz;
  }

  // mesh->CreateSerialMesh();
  mesh->CreateCells();
  mesh->ConvertMesh2VTKUnstructuredGrid();
  progress = 1.0;
}

template <class PointType>
void MeshGenerator<PointType>::MeshGeneratePolar(
        std::string meshParamFile, double& progress,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        const std::shared_ptr<MeshContainer2d<PointType>>& mesh, int flagVTK, std::string& errorMsg)
{
  mesh->clear();
  double     epsilon;
  MeshParams meshParams(meshParamFile, epsilon, errorMsg);
  if (errorMsg.size())
    return;
  mesh->xmin = boundary->xmin;
  mesh->xmax = boundary->xmax;
  mesh->ymin = boundary->ymin;
  mesh->ymax = boundary->ymax;
  PointType r[4];
  PointType phi[4];

  Dmath::Cartesian2Polar(boundary->xmin, boundary->ymin, r[0], phi[0]);
  Dmath::Cartesian2Polar(boundary->xmin, boundary->ymax, r[1], phi[1]);
  Dmath::Cartesian2Polar(boundary->xmax, boundary->ymin, r[2], phi[2]);
  Dmath::Cartesian2Polar(boundary->xmax, boundary->ymax, r[3], phi[3]);
  double phiMin = phi[0];
  double phiMax = phi[0];

  double rMin = r[0];
  double rMax = r[0];

  for (int i = 1; i < 4; i++)
  {
    if (phi[i] < phiMin)
      phiMin = phi[i];
    if (phi[i] > phiMax)
      phiMax = phi[i];

    if (r[i] < rMin)
      rMin = r[i];
    if (r[i] > rMax)
      rMax = r[i];
  }

  double yCurLoc = boundary->ymin;
  double xCurLoc = boundary->xmin;

  double tmp;
  rMax           = rMax * 1.2;
  mesh->nVertexX = meshParams.GetNumberOfSteps(0, rMin, rMax, tmp) + 4;
  mesh->h1       = tmp;
  // mesh->h1 = meshParams.GetStep(rMin, 0);

  phiMax         = phiMin + PI() / 60;
  mesh->nVertexY = meshParams.GetNumberOfSteps(1, phiMin, phiMax, tmp) + 4;
  mesh->h2       = tmp;
  // mesh->h2 = meshParams.GetStep(phiMin, 1);

  mesh->xmin = rMin;
  mesh->xmax = rMax;
  mesh->ymin = phiMin;
  mesh->ymax = phiMax;

  MeshGeneratePolar(epsilon, meshParams, progress, boundary, mesh, errorMsg);

  mesh->CreateCells();
  mesh->ConvertMesh2VTKUnstructuredGrid();
  progress = 1.0;
}

template <class PointType>
void MeshGenerator<PointType>::MeshGeneratePolar(
        double epsilon, const MeshParams& meshParams, double& progress,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        const std::shared_ptr<MeshContainer2d<PointType>>& mesh, std::string& errorMsg)
{

  mesh->flagMatrix.ones(mesh->nVertexX, mesh->nVertexY);
  mesh->templNumb.ones(mesh->nVertexX, mesh->nVertexY);
  mesh->templNumb     = mesh->templNumb * (-1);
  mesh->flagMatrix    = mesh->flagMatrix * (-1);
  mesh->flagTemplNumb = 1;

  flagMeshInit                = true;
  int                    yNum = 0;
  DGeo::Point<PointType> CurrentPoint;

  PointType yCur   = boundary->ymin;
  int       number = 0;

  CurrentPoint.x = 0;
  CurrentPoint.y = 0;
  CurrentPoint.z = 0;

  double phiStart = mesh->ymin;
  double phiEnd   = mesh->ymax;

  double phiCurrent = phiStart;

  double h2;

  while (phiCurrent < phiEnd)
  {
    progress = std::abs(phiCurrent - phiStart) / std::abs(phiEnd - phiStart);

    if (yNum == 20)
    {
      int tt = 0;
    };

    RayTracePolar(epsilon, meshParams, CurrentPoint, mesh->flagMatrix,
                  std::cos(double(2 * PI() - phiCurrent)),
                  std::sin(double(2 * PI() - phiCurrent)), yNum, boundary);
    yNum++;
    h2         = meshParams.GetStep(phiCurrent, 1);
    phiCurrent = phiCurrent + h2;
  }

  yNum   = 0;
  yCur   = boundary->ymin;
  number = 0;

  phiCurrent = phiStart;
  while (phiCurrent < phiEnd)
  {

    std::vector<DGeo::Point<PointType>>* tmpL1 = new std::vector<DGeo::Point<PointType>>();

    MeshAssemblyPolar(meshParams, CurrentPoint, tmpL1, mesh->templNumb, mesh->flagMatrix,
                      std::cos(2. * PI() - phiCurrent),
                      std::sin(2. * PI() - phiCurrent), yNum, number, boundary,
                      mesh->linearToTemplNumb);
    if (tmpL1->size() != 0)
      mesh->meshData.push_back(*tmpL1);
    yNum++;
    h2         = meshParams.GetStep(phiCurrent, 1);
    phiCurrent = phiCurrent + h2;
  }
}

template <class PointType>
bool MeshGenerator<PointType>::IsInitBoundary()
{
  if (flagBoundaryInit)
    return true;
  return false;
}

template <class PointType>
bool MeshGenerator<PointType>::IsMeshGen()
{
  if (flagMeshInit)
    return true;
  return false;
}

template <class PointType>
std::vector<double> MeshGenerator<PointType>::GetMeshParam()
{
  return meshParam;
}

template <class PointType>
void MeshGenerator<PointType>::SetMeshParam(std::vector<double> in)
{
  meshParam = in;
}

template <class PointType>
MeshGenerator<PointType>::MeshGenerator()
{
  flagBoundaryInit = false;
  flagMeshInit     = false;
  meshParam.push_back(0);
  meshParam.push_back(0);
  t1 = 0;
  t2 = 0;
}

template <class PointType>
void MeshGenerator<PointType>::MeshGenerate(
        std::string meshParamFile, double& progress,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        const std::shared_ptr<MeshContainer2d<PointType>>& mesh, int flagVTK, std::string& errorMsg)
{
  mesh->clear();
  double     epsilon;
  MeshParams meshParams(meshParamFile, epsilon, errorMsg);
  if (errorMsg.size())
    return;
  mesh->xmin = boundary->xmin;
  mesh->xmax = boundary->xmax;
  mesh->ymin = boundary->ymin;
  mesh->ymax = boundary->ymax;

  double yCurLoc = boundary->ymin;
  double xCurLoc = boundary->xmin;

  double tmp;
  mesh->nVertexX = meshParams.GetNumberOfSteps(0, boundary->xmin, boundary->xmax, tmp) + 4;
  mesh->h1       = tmp;
  mesh->nVertexY = meshParams.GetNumberOfSteps(1, boundary->ymin, boundary->ymax, tmp) + 4;
  mesh->h2       = tmp;

  MeshGenerate(epsilon, meshParams, progress, boundary, mesh, errorMsg);

  mesh->CreateCells();
  if (flagVTK == 0)
    mesh->ConvertMesh2VTKUnstructuredGrid();
  progress = 1.0;
}

template <class PointType>
void MeshGenerator<PointType>::MeshGenerate(
        double epsilon, const MeshParams& meshParams, double& progress,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        const std::shared_ptr<MeshContainer2d<PointType>>& mesh, std::string& errorMsg)
{
  mesh->flagMatrix.ones(mesh->nVertexX, mesh->nVertexY);
  mesh->templNumb.ones(mesh->nVertexX, mesh->nVertexY);
  mesh->templNumb     = mesh->templNumb * (-1);
  mesh->flagMatrix    = mesh->flagMatrix * (-1);
  mesh->flagTemplNumb = 1;

  //	mesh = MeshContainer2d <PointType>();
  // mesh = MeshContainer2d <PointType>();
  flagMeshInit                = true;
  int                    yNum = 0;
  DGeo::Point<PointType> CurrentPoint;

  PointType yCur   = boundary->ymin;
  int       number = 0;
  double    h2     = 0;
  double    h1     = 0;

  std::vector<PointType> yC;
  yC.push_back(boundary->ymin);

  while (yC.back() <= boundary->ymax + h2)
  {
    h2 = meshParams.GetStep(yC.back(), 1);
    yC.push_back(yC.back() + h2);
  }

  int              numThreads = 1;
  std::vector<int> threadIndexes(numThreads + 1);
  int              perThread = yC.size() / numThreads;
  threadIndexes[0]           = 0;

  for (int i = 0; i < numThreads - 1; i++)
    threadIndexes[i + 1] = threadIndexes[i] + perThread;

  threadIndexes[numThreads] = yC.size();

#pragma omp parallel num_threads(numThreads)
  {
    volatile int thread = omp_get_thread_num();

    for (int yNum = threadIndexes[thread]; yNum < threadIndexes[thread + 1]; yNum++)
    {
      if (thread == numThreads - 1)
        progress = double(std::abs(yNum - threadIndexes[thread])) /
                   double(std::abs(threadIndexes[thread + 1] - threadIndexes[thread]));

      CurrentPoint.x = boundary->xmin;
      CurrentPoint.y = yC[yNum];
      CurrentPoint.z = 0;

      if (yNum == yC.size() - 1)
        h2 = yC[yNum] - yC[yNum - 1];
      else
        h2 = yC[yNum + 1] - yC[yNum];

      RayTrace(epsilon, meshParams, CurrentPoint, mesh->flagMatrix, h1, h2, yNum, boundary);

    }

  }

  /*while (yCur <= boundary->ymax + h2)
  {
          progress = std::abs(yCur - boundary->ymin) / std::abs(boundary->ymax - boundary->ymin);

          CurrentPoint.x = boundary->xmin;
          CurrentPoint.y = yCur;
          CurrentPoint.z = 0;

          if (yNum == 24)
          {
                  int a = 0;
          };
          RayTrace(epsilon, meshParams, CurrentPoint, mesh->flagMatrix, h1, 0, yNum, boundary);
          yNum++;
          h2 = meshParams.GetStep(yCur, 1);
          yCur = yCur + h2;
  }*/

  yNum   = 0;
  yCur   = boundary->ymin;
  number = 0;

  while (yCur < boundary->ymax + h2)
  {
    progress = yCur / (boundary->ymax + h2);

    CurrentPoint.x = boundary->xmin;
    CurrentPoint.y = yCur;
    CurrentPoint.z = 0;

    std::vector<DGeo::Point<PointType>>* tmpL1 = new std::vector<DGeo::Point<PointType>>();

    if (yNum == 79)
    {
      int a = 0;
    }

    MeshAssembly(meshParams, CurrentPoint, tmpL1, mesh->templNumb, mesh->flagMatrix, h1, 0,
                 yNum, number, boundary, mesh->linearToTemplNumb);

    if (tmpL1->size() != 0)
      mesh->meshData.push_back(*tmpL1);
    yNum++;
    h2   = meshParams.GetStep(yCur, 1);
    yCur = yCur + h2;
  }
}

template <class PointType>
void MeshGenerator<PointType>::MeshAssembly(
        const MeshParams& meshParams, DGeo::Point<PointType> StartPoint,
        std::vector<DGeo::Point<PointType>>* PointContainer, Dmath::imat& templNumb,
        Dmath::imat& flagMartix, PointType h1, PointType h2, int yNum, int& number,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        std::vector<int>&                                     linearToTemplNumb)
{
  DGeo::Point<PointType> CurrentPoint = StartPoint;
  int                    npoints;
  for (int j = 0; j < templNumb.GetNrow() - 2; j++)
  {
    int npoints = 0;
    if (flagMartix(j + 1 - 1, yNum + 1 - 1) == 2)
      npoints++;

    if (flagMartix(j + 1, yNum + 1 - 1) == 2)
      npoints++;

    if (flagMartix(j + 1 - 1, yNum + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 + 1, yNum + 1 + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 + 1, yNum + 1) == 2)
      npoints++;

    if (flagMartix(j + 1, yNum + 1 + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 - 1, yNum + 1 + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 + 1, yNum + 1 - 1) == 2)
      npoints++;

    if (flagMartix(j + 1, yNum + 1) == 0 || flagMartix(j + 1, yNum + 1) == 2 ||
        flagMartix(j + 1, yNum + 1) == 1 || npoints > 0)
    {

      //	h2 = meshParams.GetStep(CurrentPoint.y, 1);

      CurrentPoint.isOut = false;
      if (flagMartix(j + 1, yNum + 1) != 2)
        CurrentPoint.isOut = true;
      PointContainer->push_back(CurrentPoint);
      templNumb(j + 1, yNum + 1) = number;
      number++;
      linearToTemplNumb.push_back(templNumb.getElemIndex(j + 1, yNum + 1));
    }
    h1           = meshParams.GetStep(CurrentPoint.x, 0);
    CurrentPoint = DGeo::Point<PointType>(CurrentPoint, h1, h2, 0);
  }
}

template <class PointType>
void MeshGenerator<PointType>::MeshAssemblyPolar(
        const MeshParams& meshParams, DGeo::Point<PointType> StartPoint,
        std::vector<DGeo::Point<PointType>>* PointContainer, Dmath::imat& templNumb,
        Dmath::imat& flagMartix, PointType h1, PointType h2, int yNum, int& number,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary,
        std::vector<int>&                                     linearToTemplNumb)
{
  DGeo::Point<PointType> CurrentPoint = StartPoint;
  int                    npoints;
  for (int j = 0; j < templNumb.GetNrow() - 2; j++)
  {
    int npoints = 0;
    if (flagMartix(j + 1 - 1, yNum + 1 - 1) == 2)
      npoints++;

    if (flagMartix(j + 1, yNum + 1 - 1) == 2)
      npoints++;

    if (flagMartix(j + 1 - 1, yNum + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 + 1, yNum + 1 + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 + 1, yNum + 1) == 2)
      npoints++;

    if (flagMartix(j + 1, yNum + 1 + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 - 1, yNum + 1 + 1) == 2)
      npoints++;

    if (flagMartix(j + 1 + 1, yNum + 1 - 1) == 2)
      npoints++;

    if (flagMartix(j + 1, yNum + 1) == 0 || flagMartix(j + 1, yNum + 1) == 2 ||
        flagMartix(j + 1, yNum + 1) == 1 || npoints > 0)
    {

      //	h2 = meshParams.GetStep(CurrentPoint.y, 1);

      CurrentPoint.isOut = false;
      if (flagMartix(j + 1, yNum + 1) != 2)
        CurrentPoint.isOut = true;
      PointContainer->push_back(CurrentPoint);
      templNumb(j + 1, yNum + 1) = number;
      number++;
      linearToTemplNumb.push_back(templNumb.getElemIndex(j + 1, yNum + 1));
    }

    double r     = sqrt(CurrentPoint.x * CurrentPoint.x + CurrentPoint.y * CurrentPoint.y);
    double hr    = meshParams.GetStep(r, 0);
    CurrentPoint = DGeo::Point<PointType>(CurrentPoint, hr * h1, hr * h2, 0);
  }
}

template <class PointType>
void MeshGenerator<PointType>::RayTrace(
        double epsilon, const MeshParams& meshParams, DGeo::Point<PointType> StartPoint,
        Dmath::imat& flagMartix, PointType h1, PointType h2, int yNum,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary)
{
  int j = 1;

  PointType xmax = StartPoint.x;

  DGeo::Point<PointType> CurrentPoint = StartPoint;
  xmax                                = CurrentPoint.x;

  DGeo::Point<PointType> PointTmp;

  int flag2;
  int intersections = -1;

  int val = boundary->OnBoundary(epsilon, StartPoint);

  flagMartix(1, yNum + 1)     = -1;
  int                    f1   = 1;
  DGeo::Point<PointType> ptmp = StartPoint;

  if (val != 2)
  {
    intersections           = -1; //�� �������
    flagMartix(1, yNum + 1) = 1;
    // f1 = 2;
    h1     = meshParams.GetStep(StartPoint.x, 0);
    ptmp.x = ptmp.x - h1;
    ptmp.y = ptmp.y + 0.01 * h2;
  }

  while (xmax < boundary->xmax + h1)
  {
    h1         = meshParams.GetStep(StartPoint.x, 0);
    StartPoint = DGeo::Point<PointType>(StartPoint, h1, 0, 0);

    xmax = StartPoint.x;
    val  = boundary->OnBoundary(epsilon, StartPoint);

    if (val != 2)
    {
      //�� �������
      // ptmp = StartPoint;
      //	f1 = 2;
      flagMartix(j + 1, yNum + 1) = 0;
      j++;
      continue;
    }

    if (f1 == 3 || f1 == 1)
      intersections = -1;

    if (f1 == 2)
      intersections = 0;

    DGeo::Edge<PointType> GridEdge;
    GridEdge.point1 = ptmp;

    GridEdge.point2   = StartPoint;
    GridEdge.point2.y = GridEdge.point2.y + epsilon * h2;
    // GridEdge.point1.y = GridEdge.point1.y - 0.01*h2;

    if (std::abs(StartPoint.x - 0.205) < 0.00001 && std::abs(StartPoint.y - 0.5) < 0.00001)
    {
      int tt = 0;
    }

    intersections = intersections + boundary->NumberIntersections(epsilon, GridEdge);

    flagMartix(j + 1, yNum + 1) = 2; // ������ �������
    if (intersections % 2)
      flagMartix(j + 1, yNum + 1) = 3; // ������� �������
    else
      int aa = 0;

    ptmp = StartPoint;
    f1   = flagMartix(j + 1, yNum + 1);
    j++;
  }
}

template <class PointType>
void MeshGenerator<PointType>::RayTracePolar(
        double epsilon, const MeshParams& meshParams, DGeo::Point<PointType> StartPoint,
        Dmath::imat& flagMartix, PointType h1, PointType h2, int yNum,
        const std::shared_ptr<BoundaryContainer2d<PointType>> boundary)
{

  int j = 1;

  PointType xmax = StartPoint.x;

  DGeo::Point<PointType> CurrentPoint = StartPoint;
  xmax                                = CurrentPoint.x;

  DGeo::Point<PointType> PointTmp;

  int flag2;
  int intersections = -1;

  int val = boundary->OnBoundary(epsilon, StartPoint);

  flagMartix(1, yNum + 1) = -1;

  if (val != 2)
  {
    intersections           = -1; //�� �������
    flagMartix(1, yNum + 1) = 1;
  }
  double Rmax = sqrt(boundary->xmax * boundary->xmax + boundary->ymax * boundary->ymax);
  double r    = 0;
  double rTot = 0;
  DGeo::Point<PointType> ptmp = StartPoint;
  int                    f1   = 1;
  while (r < Rmax + h1)
  {

    if (j == 148)
    {
      int tt = 0;
    }
    r = sqrt(double(StartPoint.x * StartPoint.x + StartPoint.y * StartPoint.y));

    if (sqrt(double(StartPoint.x * StartPoint.x + StartPoint.y * StartPoint.y)) > 0.15)
    {
      int tt = 0;
    }

    double hr  = meshParams.GetStep(r, 0);
    StartPoint = DGeo::Point<PointType>(StartPoint, hr * h1, hr * h2, 0);

    if (sqrt(double(StartPoint.x * StartPoint.x + StartPoint.y * StartPoint.y)) > 0.15)
    {
      int tt = 0;
    }

    if (std::abs(StartPoint.x - 0.053000000000000005) < 0.00001 &&
        std::abs(StartPoint.y - 0.12500000000000006) < 0.00001)
    {
      int kk = 0;
    }
    val = boundary->OnBoundary(epsilon, StartPoint);

    if (val != 2)
    {
      //�� �������
      flagMartix(j + 1, yNum + 1) = 0;
      j++;
      continue;
    }

    if (f1 == 3 || f1 == 1)
      intersections = -1;

    if (f1 == 2)
      intersections = 0;

    DGeo::Edge<PointType> GridEdge;
    GridEdge.point1 = ptmp;
    GridEdge.point2 = StartPoint;

    intersections = intersections + boundary->NumberIntersections(epsilon, GridEdge);

    flagMartix(j + 1, yNum + 1) = 2; // ������ �������
    if (intersections % 2)
      flagMartix(j + 1, yNum + 1) = 3; // ������� �������
    else
      int aa = 0;

    ptmp = StartPoint;
    f1   = flagMartix(j + 1, yNum + 1);
    j++;
  }
}

template <class PointType>
template <class Archive>
void MeshGenerator<PointType>::save(Archive& ar, const unsigned int) const
{
  // ar & mesh;
  ar& t1;
  ar& t2;
  ar& boundaryList;
  ar& meshParam;
}
template <class PointType>
template <class Archive>
void MeshGenerator<PointType>::load(Archive& ar, const unsigned int)
{
  //
  ar& t1;
  ar& t2;
  ar& boundaryList;
  ar& meshParam;
}

template class MeshGenerator<float>;
template class MeshGenerator<double>;

template void
MeshGenerator<float>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                            const unsigned int file_version);
template void
MeshGenerator<double>::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                             const unsigned int file_version) const;

template void
MeshGenerator<double>::load<boost::archive::binary_iarchive>(boost::archive::binary_iarchive& ar,
                                                             const unsigned int file_version);
template void
MeshGenerator<float>::save<boost::archive::binary_oarchive>(boost::archive::binary_oarchive& ar,
                                                            const unsigned int file_version) const;