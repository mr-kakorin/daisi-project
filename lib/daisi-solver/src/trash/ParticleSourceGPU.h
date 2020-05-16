#ifndef PARTICLESOURCEGPU_H
#define PARTICLESOURCEGPU_H
#include "BoundaryContainer2d.h"
#include "GPUTypes.h"
#include "Geom.h"
#include "GeomGPU.h"
#include "ParticleSource.h"
template <class PointType> class ParticleSourceEdgeGPU
{
  public:
    DGeoGPU::Edge<PointType> extractingEdge;
    double                   curveLength;
    double                   currentDensity;

    ParticleSourceEdgeGPU(){

    };
    ParticleSourceEdgeGPU(const ParticleSourceEdge<PointType>& obj)
    {
        curveLength    = obj.curveLength;
        currentDensity = obj.currentDensity;
        extractingEdge = obj.extractingEdge;
    };
    /*	ParticleSourceEdge(DGeo::Edge <PointType> Edge, double curveLengthPrev)
            {
                    extractingEdge = Edge;
                    currentDensity = 0;
                    curveLength = curveLengthPrev + Edge.length() / 2;
            };
    */
};
template <class PointType> class ParticleSource2dGPU
{
    /*static bool EdgeCmp(DGeo::Edge <PointType> & Edge1, DGeo::Edge <PointType> & Edge2)
    {
            double r1 = Edge1.Middle().Radius();
            double r2 = Edge2.Middle().Radius();
            return  r1< r2;
    };*/
  public:
    GPUTypes::vector<ParticleSourceEdgeGPU<PointType>> sourceSurface;

    ParticleSource2dGPU(){

    };
    ParticleSource2dGPU(ParticleSource2d<PointType>& obj)
    {
        sourceSurface = obj.sourceSurface;
    };

    /*
    std::vector<PointType>  GetParticle(PointType L1, PointType L2, int flag);

    double length()
    {
            return sourceSurface.back().curveLength + sourceSurface.back().extractingEdge.length() / 2; // for 2daxs
    case
    };*/
};

#endif