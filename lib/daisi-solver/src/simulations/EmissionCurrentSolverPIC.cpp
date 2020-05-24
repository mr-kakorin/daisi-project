#include "EmissionCurrentSolver.h"
#include "Dmath.h"
#include "ElectrodeCurrent.h"
#include "EmitterDevice2d.h"
#include "EmitterDevice2daxs.h"
#include "EmitterDevice3d.h"
#include "Geom.h"
#include "GridData.h"
#include "ParticleGridInterface.h"
#include "ParticleSource.h"
#include <common_tools/constants.h>

template
class EmissionCurrentSolverPIC<float>;

template
class EmissionCurrentSolverPIC<double>;

template void EmissionCurrentSolverPIC<double>::UpdateEmissionCurrent<GridData2d<double>,
    EmitterDevice2d<double>>(
    const std::shared_ptr<EmitterDevice2d<double>> &emitter,
    const std::shared_ptr<ParticleGridInterface<double>> &particleGridInterface,
    const std::shared_ptr<GridData2d<double>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void EmissionCurrentSolverPIC<double>::UpdateEmissionCurrent<GridData2daxs<double>,
    EmitterDevice2daxs<double>>(
    const std::shared_ptr<EmitterDevice2daxs<double>> &emitter,
    const std::shared_ptr<ParticleGridInterface<double>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<double>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void
EmissionCurrentSolverPIC<float>::UpdateEmissionCurrent<GridData2d<float>, EmitterDevice2d<float>>(
    const std::shared_ptr<EmitterDevice2d<float>> &emitter,
    const std::shared_ptr<ParticleGridInterface<float>> &particleGridInterface,
    const std::shared_ptr<GridData2d<float>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void EmissionCurrentSolverPIC<float>::UpdateEmissionCurrent<GridData2daxs<float>,
    EmitterDevice2daxs<float>>(
    const std::shared_ptr<EmitterDevice2daxs<float>> &emitter,
    const std::shared_ptr<ParticleGridInterface<float>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<float>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void EmissionCurrentSolverPIC<float>::UpdateEmissionCurrent<GridData2dpolar<float>,
    EmitterDevice2d<float>>(
    const std::shared_ptr<EmitterDevice2d<float>> &emitter,
    const std::shared_ptr<ParticleGridInterface<float>> &particleGridInterface,
    const std::shared_ptr<GridData2dpolar<float>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void EmissionCurrentSolverPIC<double>::UpdateEmissionCurrent<GridData2dpolar<double>,
    EmitterDevice2d<double>>(
    const std::shared_ptr<EmitterDevice2d<double>> &emitter,
    const std::shared_ptr<ParticleGridInterface<double>> &particleGridInterface,
    const std::shared_ptr<GridData2dpolar<double>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void
EmissionCurrentSolverPIC<float>::UpdateEmissionCurrent<GridData3d<float>, EmitterDevice3d<float>>(
    const std::shared_ptr<EmitterDevice3d<float>> &emitter,
    const std::shared_ptr<ParticleGridInterface<float>> &particleGridInterface,
    const std::shared_ptr<GridData3d<float>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template void EmissionCurrentSolverPIC<double>::UpdateEmissionCurrent<GridData3d<double>,
    EmitterDevice3d<double>>(
    const std::shared_ptr<EmitterDevice3d<double>> &emitter,
    const std::shared_ptr<ParticleGridInterface<double>> &particleGridInterface,
    const std::shared_ptr<GridData3d<double>> &gridData, double timeStep, int flowNumber,
    int stepNumber, double mass, double charge, int emissionType);

template<class PointType>
template<class gridDataType, class emittersType>
void EmissionCurrentSolverPIC<PointType>::UpdateEmissionCurrent(
    const std::shared_ptr<emittersType> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<gridDataType> &gridData, double timeStep, int flowNumber, int stepNumber,
    double mass, double charge, int emissionType) {
  if (emissionType == 0) {
    // double T    = emitter->GetDistribParams()[0];
    // double n    = emitter->GetDistribParams()[1] * 1e6;
    //double jmax = 0.43 * n * std::abs(charge) *
    //              sqrt(2 * std::abs(charge) * T / mass); //������� ���� ��� ������� �� ������
    double jmax = 1.1 * 1e15;
    for (int i = 0; i < emitter->GetParticleSources().size(); i++)
      this->SetValueOnSource(emitter->GetParticleSources()[i], {jmax}, flowNumber, 2);

    UpdateEmissionCurrent(emitter, particleGridInterface, gridData, PointType(timeStep),
                          flowNumber, stepNumber, mass, charge);
  };
  if (emissionType == 1) {
    std::vector<std::vector<double>> cd =
        emitter->GetAssignedElectrode()->GetElectrodeValueD(1);
    std::vector<double> currentDensAtPoints(this->istarts[flowNumber].size() - 1);
    std::vector<double> jmax(this->istarts[flowNumber].size() - 1);
    double n = 1e10 * 1e6;

    for (int i = 0; i < this->istarts[flowNumber].size() - 1; i++) {
      int k = (this->istarts[flowNumber][i] + this->istarts[flowNumber][i + 1]) / 2;
      double ss = emitter->GetParticleSources()[0]->sourceSurface[k].curveLength;
      currentDensAtPoints[i] = Dmath::Interpolate(cd[0], cd[1], ss);
      if (std::abs(currentDensAtPoints[i]) > 1)
        jmax[i] = 1e22;
        // jmax[i] = 0.43*n*std::abs(charge)*sqrt(2 * std::abs(charge) * 100 / mass);
      else
        jmax[i] = 0;
    }
    for (int i = 0; i < emitter->GetParticleSources().size(); i++)
      this->SetValueOnSource(emitter->GetParticleSources()[i], jmax, flowNumber, 2);

    UpdateEmissionCurrent(emitter, particleGridInterface, gridData, PointType(timeStep),
                          flowNumber, stepNumber, mass, charge);
  }
  if (emissionType != 5)
    for (int i = 0; i < emitter->GetParticleSources().size(); i++)
      this->CalculateCathodeFields(emitter->GetParticleSources()[i], gridData, flowNumber);
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::SetEmissionCurrent(
    const std::shared_ptr<EmitterDevice3d<PointType>> &emitter, PointType current) {
  PointType dcurrentDensity = 0;

  for (int k = 0; k < emitter->GetParticleSource().size(); k++) {
    dcurrentDensity = (current / emitter->GetParticleSource().size() -
                       emitter->GetParticleSource()[k]->GetEmissionCurrent(0)) /
                      emitter->GetParticleSource()[k]->length();
    for (int i = 0; i < emitter->GetParticleSource()[k]->sourceSurface.size(); i++) {
      emitter->GetParticleSource()[k]->sourceSurface[i].currentDensity =
          emitter->GetParticleSource()[k]->sourceSurface[i].currentDensity + dcurrentDensity;
    };
  }
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::UpdateEmissionCurrent(
    const std::shared_ptr<EmitterDevice2d<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2dpolar<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {

};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::UpdateEmissionCurrent(
    const std::shared_ptr<EmitterDevice3d<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData3d<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {

};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::SetEmissionCurrent(
    const std::shared_ptr<EmitterDevice2d<PointType>> &emitter, PointType current) {
  PointType dcurrentDensity = 0;
  dcurrentDensity = (current - emitter->GetParticleSource()->GetEmissionCurrent(0)) /
                    emitter->GetParticleSource()->length();
  for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++) {
    emitter->GetParticleSource()->sourceSurface[i].currentDensity =
        emitter->GetParticleSource()->sourceSurface[i].currentDensity + dcurrentDensity;
  };
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::SetEmissionCurrent(
    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter, PointType current) {
  if (!emitter->GetParticleSource())
    return;
  PointType dcurrentDensity = 0;
  dcurrentDensity = (current - emitter->GetParticleSource()->GetEmissionCurrent(1)) /
                    (2 * commtools::PI() *
                     emitter->GetParticleSource()->sourceSurface[0].extractingEdge->Middle().x *
                     emitter->GetParticleSource()->length());
  for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++) {
    emitter->GetParticleSource()->sourceSurface[i].currentDensity =
        emitter->GetParticleSource()->sourceSurface[i].currentDensity + dcurrentDensity;
  };
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::UpdateEmissionCurrent(
    const std::shared_ptr<EmitterDevice2d<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2d<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {

  PointType z;

  PointType dr;
  PointType dz;
  PointType Er;
  PointType J;
  PointType Er1;
  PointType jm[2] = {4000.0, 7000.0};
  int k = 0;
  int baseOld = -1;
  PointType Jold;

  PointType ErAverage = 0;

  double Kmin = 0.99;
  double Kmax = 1.002;

  // if (flowNumber == 0)
  //	return;
  double Ex;
  double Ey;

  PointType x1;
  PointType y1;

  PointType x0;
  PointType y0;

  std::vector<double> jv;
  std::vector<double> dvs;
  for (int i = 0; i < this->points1[flowNumber].size(); i++) {

    //	r1 = this->points2[flowNumber][i].x;
    // r0 = this->points1[flowNumber][i].x;
    // r1 = 0.0005;
    // r0 = 0.0055;
    //	z1 = this->points2[flowNumber][i].y;
    //	z0 = this->points1[flowNumber][i].y;

    x0 = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0];
    x1 = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[1];

    y0 = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0];
    y1 = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[1];

    double d = sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));

    // double V1 = gridData->interpolatePoint(r1, z1, 0, std::string("V, V"), 0)+
    // gridData->interpolatePoint(r1, z1, 0, std::string("V, V"), 1);
    //	double V0 = gridData->interpolatePoint(r0, z1, 0, std::string("V, V"), 0)+
    // gridData->interpolatePoint(r0, z1, 0, std::string("V, V"), 1);

    double V1 = gridData->interpolatePoint(x1, y1, 0, std::string("V, V"), 2);
    double V0 = gridData->interpolatePoint(x0, y0, 0, std::string("V, V"), 2);

    gridData->interpolatePoint(this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
                               this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 0, Ex,
                               Ey);

    double cathField = (Ex * this->nearCathodeVolumes[flowNumber][i].normalX[0] +
                        Ey * this->nearCathodeVolumes[flowNumber][i].normalY[0]);

    double dV = std::abs(V1 - V0);
    dvs.push_back(dV);
    double jnew = 4.0 * commtools::VACUUM_PERMITTIVITY() * (1.0 / 9.0) *
                  sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) / (d * d);

    if (cathField * charge > 0)
      jnew = 0;

    jv.push_back(jnew);
  };

  this->SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);

  /*double Ex;
  double Ey;
  double Ex1;
  double Ey1;
  double Ex11;
  double Ey11;
  double Ex111;
  double ExCol;
  double EyCol;
  std::vector<double > nearCathField;

  int N = 0;


  nearCathField.resize(this->points1[flowNumber].size());
  for (int i = 0; i < this->points1[flowNumber].size(); i++)
  {
          nearCathField[i] = 0;
          for (int ii = 1; ii < 4; ii++)
          {
                  gridData->interpolatePoint(this->nearCathodeVolumes[flowNumber][i].fieldPointsX[ii],
  this->nearCathodeVolumes[flowNumber][i].fieldPointsY[ii], 0, Ex, Ey);

                  nearCathField[i] = nearCathField[i] +
  (Ex*this->nearCathodeVolumes[flowNumber][i].normalX[ii] +
  Ey*this->nearCathodeVolumes[flowNumber][i].normalY[ii]);
          }
  }

  std::vector<double> jvC(this->points1[flowNumber].size());



  std::vector<double> eerrT;
  std::vector<double> eerrT1;

  for (int it = 0; it < 50; it++)
  {
          double x1 = 0.0005 + it*0.0001;
          double x0 = x1 + 0.0001;

          double y1 = 0.0002;
          double y0 = 0.0003;

          gridData->interpolatePoint(x0, (y0 + y1) / 2, 0, Ex, Ey);
          gridData->interpolatePoint(x1, (y0 + y1) / 2, 0, Ex1, Ey1);

          double CellFlow1 = (Ex - Ex1)*(y0 - y1)*commtools::VACUUM_PERMITTIVITY();



          PointType rho2 = gridData->interpolatePoint((x0+ x1)/2, (y1 + y0) / 2, 0,
  std::string("Charge density, cl/m^3"), 0); PointType rho11 = gridData->interpolatePoint(x1 -
  0.0001, (y1 + y0) / 2, 0, std::string("Charge density, cl/m^3"), 0);

          PointType rho12 = gridData->interpolatePoint(x1, (y1 + y0) / 2, 0, std::string("Charge
  density, cl/m^3"), 0);

          PointType chargeL = rho2*(x0 - x1)*(y1 - y0);
          PointType chargeL1 = (2* rho12- rho11)*(x0 - x1)*(y1 - y0);

          eerrT.push_back((std::abs(CellFlow1) - std::abs(chargeL)) / std::abs(chargeL));
          eerrT1.push_back((std::abs(CellFlow1) - std::abs(chargeL1)) / std::abs(chargeL1));

  }



  for (int i = 0; i < this->points1[flowNumber].size(); i++)
  {

          PointType CellFlow = 0;
          CellFlow = CellFlow + nearCathField[i] *
  this->nearCathodeVolumes[flowNumber][i].Areas[1];

          gridData->interpolatePoint(this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
  this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 0, Ex, Ey);

          double cathField = (Ex*this->nearCathodeVolumes[flowNumber][i].normalX[0] +
  Ey*this->nearCathodeVolumes[flowNumber][i].normalY[0]);

          double cathFlow = 4*cathField*this->nearCathodeVolumes[flowNumber][i].Areas[1] *
  commtools::VACUUM_PERMITTIVITY();

          CellFlow = CellFlow*commtools::VACUUM_PERMITTIVITY();

          double x1 = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0];
          double x0 = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[1];

  //x1 = 0.0054;
          //x0 = 0.0055;

          double y1 = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0];
          double y0 = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[1];

          PointType rho2 = gridData->interpolatePoint((x1 + x0) / 2, (y0 + y1) / 2, 0,
  std::string("Charge density, cl/m^3"), 0);

          PointType chargeL = rho2*this->nearCathodeVolumes[flowNumber][i].volume;
          PointType chargeNew = std::abs(CellFlow - chargeL);


          double decl = CathodeFields[flowNumber][i] / E0[flowNumber][i];

          jvC[i] = (K[flowNumber][i]+0.1) *chargeNew /
  (this->nearCathodeVolumes[flowNumber][i].Areas[0]
  * timeStep);

          //if(decl>0.25)
  //		K[flowNumber][i] = K[flowNumber][i] * 1.1;

  //	if (decl<0.1)
  //		K[flowNumber][i] = K[flowNumber][i] * 0.9;

  //	y1 = 0.0002;
          /*y0 = 0.0003;
          double dy = 0.0001;
          y1 = 0.0002;

          double y00 = 0.00025;
          double y01 = 0.00035;

          double x00 = 0.0054;
          double x01 = 0.0055;

          gridData->interpolatePoint(x0, y0, 0, Ex, Ey);
          gridData->interpolatePoint(x1, y0, 0, Ex1, Ey1);

          gridData->interpolatePoint(x00, (y00 + y01) / 2, 0, Ex11, Ey);
          gridData->interpolatePoint(x01, (y00 + y01) / 2, 0, Ex111, Ey1);

          double CellFlow1 = (Ex - Ex1)* dy *commtools::VACUUM_PERMITTIVITY();

          double CellFlow11 = (Ex11 - Ex111)* (y00 - y01) *commtools::VACUUM_PERMITTIVITY();


          PointType rho2 = gridData->interpolatePoint((x1+ x0)/2, y0, 0, std::string("Charge
  density, cl/m^3"), 0); PointType rho21 = gridData->interpolatePoint((x00 + x01) / 2, (y00 +
  y01)/2, 0, std::string("Charge density, cl/m^3"), 0);

          PointType chargeL = rho2*(x0 - x1)*dy;
          PointType chargeL1 = rho21*(x00 - x01)*(y00 - y01);

  //	PointType chargeNew = std::abs(CellFlow - chargeL);
          PointType chargeNew = std::abs(cathFlow - chargeL);*/
  //	jvC[i] = chargeNew/( this->nearCathodeVolumes[flowNumber][i].Areas[0]* timeStep);

  /*if (cathField*charge > 0)
  {
          jvC[i] = 0;
          K[flowNumber][i] = K[flowNumber][i] * 0.8;
  }
}*/

  // for (int i = 0; i<emitter->GetParticleSources().size(); i++)
  //	this->SetValueOnSource(emitter->GetParticleSources()[i], jvC, flowNumber, 1);

  // double decl = CathodeFields[flowNumber][i] / E0[flowNumber][i];
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::UpdateEmissionCurrent(
    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {
  PointType r1;
  PointType r0;
  PointType z1;
  PointType z0;

  PointType Er1;

  PointType ErAverage = 0;
  int k = 0;
  for (int i = 0; i < this->points1[flowNumber].size(); i++) {

    r1 = this->points2[flowNumber][i].x;
    r0 = this->points1[flowNumber][i].x;

    z1 = this->points2[flowNumber][i].y;
    z0 = this->points1[flowNumber][i].y;

    Er1 = gridData->interpolatePoint(r0, z0, 0, std::string("Er, V/m"), 0);

    ErAverage = ErAverage + Er1;
    k++;
  }

  ErAverage = ErAverage / k;
  emitter->GetParticleSource()->setErAverage(ErAverage);

  if (emitter->GetEmissionType() == 1)
    return;

  if (stepNumber % emitter->GetEmitPeriod())
    return;

  switch (this->algorithm) {
    case 0:
      VirtualDiode(emitter, particleGridInterface, gridData, timeStep, flowNumber, stepNumber,
                   mass, charge);
      break;
    case 1:
      ChargeConserving(emitter, particleGridInterface, gridData, timeStep, flowNumber, stepNumber,
                       mass, charge);
      break;
    case 2:
      Poisson(emitter, particleGridInterface, gridData, timeStep, flowNumber, stepNumber, mass,
              charge);
      break;
  }
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::ChargeConserving(
    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {
  PointType r;
  PointType z;

  PointType dr;
  PointType dz;
  PointType Er;
  PointType r1;
  PointType r2;
  PointType J;
  PointType Er1;
  PointType jm[2] = {10000.0, 7000.0};
  int k = 0;
  int baseOld = -1;
  PointType Jold;
  PointType base;

  PointType ErAverage = 0;

  double Kmin = 0.99;
  double Kmax = 1.002;

  // if (flowNumber == 0)
  //	return;
  for (int i = 0; i < this->points1[flowNumber].size(); i++) {

    PointType r1 = this->points2[flowNumber][i].x;
    PointType r0 = this->points1[flowNumber][i].x;

    PointType z1 = this->points2[flowNumber][i].y;
    PointType z0 = this->points1[flowNumber][i].y;

    base = this->cellNumbers[flowNumber][i];

    if (base == baseOld) {
      emitter->GetParticleSource()->sourceSurface[i].currentDensity = std::abs(Jold);
      continue;
    }

    dr = particleGridInterface->GetH1(base);
    dz = particleGridInterface->GetH2(base);

    if (flowNumber == 0)
      r1 = gridData->Getr()[base] - dr / 2;

    if (flowNumber == 1)
      r1 = gridData->Getr()[base] + dr;

    r2 = gridData->Getr()[base] + dr / 2;

    Er = gridData->interpolatePoint(r1, z1, 0, std::string("Er, V/m"), 0) +
         gridData->interpolatePoint(r1, z1, 0, std::string("Er, V/m"), 1);

    PointType v1 = commtools::VACUUM_PERMITTIVITY() * Er * 2 * r1;
    PointType v2 = 1 > gridData->Getrho()[base] * (r2 * r2 - r1 * r1);

    PointType qold = 1 > gridData->Getrho()[base] * (r2 * r2 - r1 * r1) * commtools::PI() * dz;
    PointType qnew =
        -commtools::VACUUM_PERMITTIVITY() * Er * 2 * commtools::PI() * r1 * dz - qold;

    PointType rho1 = gridData->interpolatePoint(r1, gridData->Getz()[base], 0,
                                                std::string("Charge density, cl/m^3"), 0);
    PointType rho2 = gridData->interpolatePoint(r2, gridData->Getz()[base], 0,
                                                std::string("Charge density, cl/m^3"), 0);

    //	PointType rho = (rho1 + rho2)/2;

    PointType rho =
        gridData->interpolatePoint(gridData->Getr()[base] + dr / 2, gridData->Getz()[base], 0,
                                   std::string("Charge density, cl/m^3"), 0);

    //	J = (-commtools::VACUUM_PERMITTIVITY()*Er * 2 * r1 - 1 > gridData->rho()[base] *
    //(r2*r2 - r1*r1)) / (timeStep
    //* 2 * r2);

    if (flowNumber == 0)
      J = (-commtools::VACUUM_PERMITTIVITY() * Er * 2 * r1 -
           1 * rho * (gridData->Getr()[base] > gridData->Getr()[base] - r1 * r1)) /
          (emitter->GetEmitPeriod() * timeStep * 2 > gridData->Getr()[base]);

    if (flowNumber == 1)
      J = (commtools::VACUUM_PERMITTIVITY() * Er * 2 * (r1) - 1 * (rho) *
                                                              ((r1) * (r1) - gridData->Getr()[base] >
                                                               gridData->Getr()[base])) /
          (emitter->GetEmitPeriod() * timeStep * 2 > gridData->Getr()[base]);

    // phoOld[k] = rho;
    //	J = (commtools::VACUUM_PERMITTIVITY()*Er * 2 * r1 - 1 * rho * (r1*r1 -
    // gridData->Getr()[base] >
    // gridData->Getr()[base])) / (emitter->GetEmitPeriod()*timeStep * 2 >
    // gridData->Getr()[base]);

    if (std::abs(J) < jm[flowNumber])
      J = std::abs(J);
    else
      J = jm[flowNumber];

    //	sr0 = this->points1[flowNumber][i].x;

    z0 = this->points1[flowNumber][i].y;

    Er1 = gridData->interpolatePoint(r0, z0, 0, std::string("Er, V/m"), 0);

    if (Er < 0 || Er1 < 0)
      J = 0;
    k++;

    double K = emitter->GetParticleSource()->sourceSurface[i].currentDensity / std::abs(J);

    if (K > Kmax)
      K = Kmax;

    if (K < Kmin)
      K = Kmin;

    //	double jnew = emitter->GetParticleSource()->sourceSurface[i].currentDensity / K;
    double omega = 0.3;
    //	emitter->GetParticleSource()->sourceSurface[i].currentDensity =
    // emitter->GetParticleSource()->sourceSurface[i].currentDensity*(1 - omega) + omega*jnew;

    emitter->GetParticleSource()->sourceSurface[i].currentDensity = J;

    baseOld = base;
    Jold = J;
  };
};

#include <fstream>

std::vector<double> gradient(std::vector<double> &x) {
  std::vector<double> result(x.size());
  //result.reserve(x.size());
  double h = 1e-8;
  result[0] = (x[1] - x[0]) / h;
  result[result.size() - 1] = (x[result.size() - 1] - x[result.size() - 2]) / h;
  for (int i = 1; i < result.size() - 1; ++i) {
    result[i] = (x[i + 1] - x[i - 1]) / (2 * h);
  }
  return result;
}


template<class PointType>
void EmissionCurrentSolverPIC<PointType>::VirtualDiode(
    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {
  PointType ErAverage = 0;
  std::vector<double> jv;
  double jnew;
  PointType r0, r1, z0, z1;

  double wf = 3.5;
  double A = 1.54141e-6;
  double B = 6.8309e9;
  //double t2 = 1.11619;
  double yc = 0.3795e-4;
  double max = 0;
  double pol[5] = {-7.535610479528366e-39, 1.094842614232595e-29, 7.580173203931757e-20, -2.183803187849293e-10,
                   2.848487403712323 + 50};
  //double powerlaw[2] = {0.601106264666250  -0.018556094759490};
  auto polylaw = [pol](double E) {
    double K = 0.;
    for (int j = 0; j < 5; ++j) {
      K += pol[j] * pow(E, 4 - j);
    }
    return K;
  };
  double powers[2] = {9.592000846639715, -1.032789333272806};
//  auto powerlaw = [powers](double E) {
//    return pow(10, powers[0]) * pow(E, powers[1]);
//  };
  double powers2[2] = {6.23830494009342,
                       -0.661449170433808};//{3.11703131996135, -0.305484211026773};//{3.77621334008232,	-0.117650789685004};
//  auto powerlaw = [powers](double E) {
//    return pow(10, powers[0]) * pow(E, powers[1]);
//  };
  auto powerlaw = [powers2](double E) {
    double lin[2] = {-0.000000013559322,
                     5.067796610169492};//{-0.000000010169492,   5.050847457627119};//{-0.000000015384615,   5.076923076923077};//{-0.000000020000000,   5.448608000000000};//{-0.000000050000000,   5.621520000000000};//{-0.000000080000000,   6.794432000000001};//{-0.000000100000000,   9.243040000000001};
    auto c = lin[0] * E + lin[1];
    return c * (pow(10, powers2[0]) * pow(E, powers2[1]));// + (powers2[1]) * log10(E)));
    //return pow(10, powers[0]) * pow(E, powers[1]);
  };
//  std::vector<double> x;
//  std::vector<double> y;
//
//  for (int i = 0; i < this->points1[flowNumber].size(); i++) {
//    x.push_back(this->points1[flowNumber][i].x);
//    y.push_back(this->points1[flowNumber][i].y);
//  }
//  auto dx = gradient(x);
//  auto dy = gradient(y);
  //std::ofstream out;
  //out.open("nearCathodeVolumes_polylaw");
  for (int i = 0; i < this->points1[flowNumber].size(); i++) {
    double Er;
    double Ez;
    //double r = this->points1[flowNumber][i].x;
    //double z = this->points1[flowNumber][i].y;
    double r = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0];
    double z = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0];
    //gridData->interpolatePoint(this->points1[flowNumber][i].x, this->points1[flowNumber][i].y, 0, Er,
    //                           Ez);
    gridData->interpolatePoint(r,
                               z, 0, Er,
                               Ez);
    double cathField = (Er * this->nearCathodeVolumes[flowNumber][i].normalX[0] +
                        Ez * this->nearCathodeVolumes[flowNumber][i].normalY[0]);
//    double cathField = (Er * (-dy[i]) +
//                        Ez * (dx[i]));
    double E_napryazennost = sqrt(Er * Er + Ez * Ez);
    double K_rough = powerlaw(E_napryazennost);

    if (cathField * charge > 0) {
      jnew = 0;
//      out << r << " "
//          << z << " " << E_napryazennost << " " << K_rough
//          << " " << jnew << std::endl;
      ErAverage = ErAverage + std::abs(cathField);
      jv.push_back(jnew);
      continue;
    }


    double KE = E_napryazennost * K_rough;
    //double v  = 0.9387 - yc * yc * KE / (wf * wf);
    //double t2 = 0.9387 + yc * yc * KE / (wf * wf) / 3;
    double y_coeff = yc * sqrt(KE) / wf;
	double v  = 1 - y_coeff* y_coeff / 3 * (3 - log(y_coeff));
	double t2 = 1 + y_coeff* y_coeff / 9 * (1 - log(y_coeff));

    jnew = A * KE * KE * exp(-B * (sqrt(wf) * wf) * v / KE) / (wf * t2);
//    out << r << " "
//        << z << " " << E_napryazennost << " " << K_rough
//        << " " << jnew << std::endl;
    if (std::abs(jnew) > max)
      max = std::abs(jnew);
    ErAverage = ErAverage + std::abs(cathField);

    jv.push_back(jnew);
  };
//  out << "END" << std::endl;
//  out.close();
  this->SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);
};

//template<class PointType>
//void EmissionCurrentSolverPIC<PointType>::VirtualDiode(
//    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter,
//    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
//    const std::shared_ptr<GridData2daxs<PointType>> &gridData, PointType timeStep, int flowNumber,
//    int stepNumber, double mass, double charge) {
//  PointType ErAverage = 0;
//  std::vector<double> jv;
//  double jnew;
//  PointType r0, r1, z0, z1;
//
//  double wf = 3.5;
//  double A = 1.54141e-6;
//  double B = 6.8309e9;
//  double t2 = 1.11619;
//  double yc = 0.3795e-4;
//  double max = 0;
//  double pol[5] = {-7.535610479528366e-39, 1.094842614232595e-29, 7.580173203931757e-20, -2.183803187849293e-10,
//                   2.848487403712323 + 50};
//  //double powerlaw[2] = {0.601106264666250  -0.018556094759490};
//  auto polylaw = [pol](double E) {
//    double K = 0.;
//    for (int j = 0; j < 5; ++j) {
//      K += pol[j] * pow(E, 4 - j);
//    }
//    return K;
//  };
//  double powers[2] = {9.592000846639715, -1.032789333272806};
//  double powers2[2] = {6.23830494009342,
//                       -0.661449170433808};//{3.11703131996135, -0.305484211026773};//{3.77621334008232,	-0.117650789685004};
//  //double hypdata[3] = {-0.000000013856303 * 1e9, 5.108321354040379 * 1e9, 0.012775077318344 * 1e9};
//  //double hypdata[3] = {-0.000000000253663 * 1e9,
//  //                     3.360986090872413 * 1e9,
//  //                    0.012857880887026 * 1e9};
//  double hypdata[3] = {-0.000000000941878 * 1e9,
//                       3.257770424128304 * 1e9,
//                       0.006056767808293 * 1e9};
//  auto hypb = [hypdata](double E) {
//    return hypdata[0] + hypdata[1] / (E + hypdata[2]);
//  };
//  //  auto powerlaw = [powers](double E) {
////    return pow(10, powers[0]) * pow(E, powers[1]);
////  };
//  auto powerlaw = [powers2](double E) {
//    double lin[2] = {-0.000000013559322,
//                     5.067796610169492};//{-0.000000010169492,   5.050847457627119};//{-0.000000015384615,   5.076923076923077};//{-0.000000020000000,   5.448608000000000};//{-0.000000050000000,   5.621520000000000};//{-0.000000080000000,   6.794432000000001};//{-0.000000100000000,   9.243040000000001};
//    //double lin[2] = {-0.000000004512563, 4.522562814070352};
//    auto c = lin[0] * E + lin[1];
//    return c * (pow(10, powers2[0]) * pow(E, powers2[1]));// + (powers2[1]) * log10(E)));
//    //return pow(10, powers[0]) * pow(E, powers[1]);
//  };
//  std::ofstream out;
//  out.open("nearCathodeVolumes_polylaw");
//  std::vector<double> x;
//  std::vector<double> y;
//
//  for (int i = 0; i < this->points1[flowNumber].size(); i++) {
//    x.push_back(this->points1[flowNumber][i].x);
//    y.push_back(this->points1[flowNumber][i].y);
//  }
//  auto dx = gradient(x);
//  auto dy = gradient(y);
//
//  for (int i = 0; i < this->points1[flowNumber].size(); i++) {
//    double Er;
//    double Ez;
//    double r = this->points1[flowNumber][i].x;
//    double z = this->points1[flowNumber][i].y;
//    //double r = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0];
//    //double z = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0];
//    //gridData->interpolatePoint(this->points1[flowNumber][i].x, this->points1[flowNumber][i].y, 0, Er,
//    //                           Ez);
//    gridData->interpolatePoint(r,
//                               z, 0, Er,
//                               Ez);
//    //double cathField = (Er * this->nearCathodeVolumes[flowNumber][i].normalX[0] +
//    //                    Ez * this->nearCathodeVolumes[flowNumber][i].normalY[0]);
//    double cathField = (Er * (-dy[i]) +
//                        Ez * (dx[i]));
//
//    double E_napryazennost = sqrt(Er * Er + Ez * Ez);
//    double K_rough = powerlaw(E_napryazennost);//powerlaw(E_napryazennost, this->points1[flowNumber].size() / 2, i);
//
//    if (cathField * charge > 0) {
//      jnew = 0;
//      out << r << " "
//          << z << " " << E_napryazennost << " " << K_rough
//          << " " << jnew << std::endl;
//      ErAverage = ErAverage + std::abs(cathField);
//      jv.push_back(jnew);
//      continue;
//    }
//
//
//    double KE = E_napryazennost * K_rough;
//    double v = 0.95 - yc * yc * KE / (wf * wf);
//    jnew = A * KE * KE * exp(-B * (sqrt(wf) * wf) * v / KE) / (wf * t2);
//    out << r << " "
//        << z << " " << E_napryazennost << " " << K_rough
//        << " " << jnew << std::endl;
//    if (std::abs(jnew) > max)
//      max = std::abs(jnew);
//    ErAverage = ErAverage + std::abs(cathField);
//
//    jv.push_back(jnew);
//  };
//  out << "END" << std::endl;
//  out.close();
//  this->SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);
//};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::VirtualDiode1(
    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {

  // return;
  PointType r;
  PointType z;

  PointType dr;
  PointType dz;
  PointType Er;
  PointType J;
  PointType Er1;
  PointType jm[2] = {4000.0, 7000.0};
  int k = 0;
  int baseOld = -1;
  PointType Jold;

  PointType ErAverage = 0;

  double Kmin = 0.99;
  double Kmax = 1.002;

  // if (flowNumber == 0)
  //	return;

  PointType r1;
  PointType z1;

  PointType r0;
  PointType z0;

  std::vector<double> jv;
  double jnew;

  for (int i = 0; i < this->points1[flowNumber].size(); i++) {

    r1 = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[1];
    r0 = this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0];

    z1 = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[1];
    z0 = this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0];

    double zCenter = -r0 * (z1 - z0) / (r1 - r0) + z0;
    /*if (r1<r0)
    r1 = r1 + 0.8*(r0 - r1);
    else
    r1 = r1 - 0.8*(r1 - r0);*/

    /*double g = r0 / r1;
    double	gamma = log(g);
    double angle = atan((z1 - z0) / (r1 - r0));


    double beta = gamma - 0.4*pow(gamma, 2.0) + 0.091667*pow(gamma, 3.0) -
    0.014242*pow(gamma, 4.0) + 0.001679*pow(gamma, 5.0);

    */

    double V1 = gridData->interpolatePoint(r1, z1, 0, std::string("V, V"), 0) +
                gridData->interpolatePoint(r1, z1, 0, std::string("V, V"), 1);
    double V0 = gridData->interpolatePoint(r0, z0, 0, std::string("V, V"), 0) +
                gridData->interpolatePoint(r0, z0, 0, std::string("V, V"), 1);

    double dV = std::abs(V1 - V0);
    double angle = atan((z1 - z0) / (r1 - r0));

    if (std::abs(r1 - r0) < 1e-10) {
      zCenter = 0;
      angle = 3.14;
    }
    if (std::abs(angle) < 0.1) {
      double g = r1 / r0;
      double gamma = log(g);

      double beta = gamma - 0.4 * pow(gamma, 2.0) + 0.091667 * pow(gamma, 3.0) -
                    0.014242 * pow(gamma, 4.0) + 0.001679 * pow(gamma, 5.0);

      /*if (CathodeFields[flowNumber][i] * charge < 0)
       dV = dV / 2;*/

      jnew = 4.0 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * (2.0 / 9.0) *
             sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) /
             (r1 * beta * beta * 2.0 * commtools::PI() * r0);
    } else {

      double rr0 = sqrt(r0 * r0 + (z0 - zCenter) * (z0 - zCenter));
      double rr1 = sqrt(r1 * r1 + (z1 - zCenter) * (z1 - zCenter));

      double gs = rr1 / rr0;
      double gammas = log(gs);

      double alpha = gammas - 0.3 * pow(gammas, 2.0) + 0.075 * pow(gammas, 3.0) -
                     0.014318 * pow(gammas, 4.0) + 0.0022161 * pow(gammas, 5.0);

      jnew = 8.0 * commtools::PI() * commtools::VACUUM_PERMITTIVITY() * (2.0 / 9.0) *
             sqrt(2.0 * std::abs(charge) / mass) * (sqrt(dV) * dV) /
             (alpha * alpha * 4 * commtools::PI() * rr0 * rr0);
    }
    /* rr0 = r0;
     rr1 = r1;

    double gs = rr0 / rr1;
    double	gammas = log(gs);

    double alpha = gammas - 0.3*pow(gammas, 2.0) + 0.075*pow(gammas, 3.0) -
    0.014318*pow(gammas, 4.0) + 0.0022161*pow(gammas, 5.0);


    double jnew1 = 8.0 * commtools::PI()*commtools::VACUUM_PERMITTIVITY()*(2.0 / 9.0)*sqrt(2.0 *
    std::abs(charge) / mass)*(sqrt(dV)*dV) / (alpha * alpha * 4 * commtools::PI() * rr0 * rr0);

    */

    double Er;
    double Ez;

    gridData->interpolatePoint(this->nearCathodeVolumes[flowNumber][i].fieldPointsX[0],
                               this->nearCathodeVolumes[flowNumber][i].fieldPointsY[0], 0, Er,
                               Ez);

    double cathField = (Er * this->nearCathodeVolumes[flowNumber][i].normalX[0] +
                        Ez * this->nearCathodeVolumes[flowNumber][i].normalY[0]);

    if (cathField * charge > 0)
      jnew = 0;

    ErAverage = ErAverage + std::abs(cathField);

    jv.push_back(jnew);
  };

  this->SetValueOnSource(emitter->GetParticleSource(), jv, flowNumber, 1);
};

template<class PointType>
void EmissionCurrentSolverPIC<PointType>::Poisson(
    const std::shared_ptr<EmitterDevice2daxs<PointType>> &emitter,
    const std::shared_ptr<ParticleGridInterface<PointType>> &particleGridInterface,
    const std::shared_ptr<GridData2daxs<PointType>> &gridData, PointType timeStep, int flowNumber,
    int stepNumber, double mass, double charge) {

  // return;
  PointType r;
  PointType z;

  PointType dr;
  PointType dz;
  PointType Er;
  PointType r1;
  PointType r2;
  PointType J;
  PointType Er1;
  PointType jm[2] = {4000.0, 7000.0};
  int k = 0;
  int baseOld = -1;
  PointType Jold;
  PointType base;

  PointType ErAverage = 0;

  double Kmin = 0.99;
  double Kmax = 1.002;

  // if (flowNumber == 0)
  //	return;

  /*
  for (int i = 0; i < emitter->GetParticleSource()->sourceSurface.size(); i++)
  {
          base = this->cellNumbers[flowNumber][i];

          if (base == baseOld)
          {
                  emitter->GetParticleSource()->sourceSurface[i].currentDensity = std::abs(Jold);
                  continue;
          }

          dr = particleGridInterface->GetH1(base);
          dz = particleGridInterface->GetH2(base);

          if (flowNumber == 0)
                  r1 = gridData->Getr()[base] - dr / 2;

          if (flowNumber == 1)
                  r1 = gridData->Getr()[base] + 2 * dr;

          double g = r1 / gridData->Getr()[base];

          double b;
          for (int i = 0; i < rmass.size() - 1; i++)
          {
                  if (rmass[i] < g && rmass[i + 1] >= g)
                  {
                          double h = rmass[i + 1] - rmass[i];
                          double v2 = (g - rmass[i]) / h;
                          double v1 = 1 - v2;
                          b = v1*bmass1[i] + v2*bmass1[i + 1];
                          break;
                  };
          };

          double V1 = gridData->interpolatePoint(r1, gridData->Getz()[base], 0, std::string("V,
  V")); double V0 = gridData->interpolatePoint(gridData->Getr()[base], gridData->Getz()[base], 0,
  std::string("V, V"));

          double dV = std::abs(V1 - V0);

          double A = 4.0 * commtools::VACUUM_PERMITTIVITY()*sqrt(2.0 * charge / mass) / (9 *
  commtools::PI());

          double omega = 0.3;

          double jold = emitter->GetParticleSource()->sourceSurface[i].currentDensity;

          double jnew = sqrt(jold*jold + A*A*V1*(b*0.3*Er1 / (2 * r1))*(0.3*Er1 / (2 * r1)));

          if (Er1<0)
                  jnew = sqrt(jold*jold - 10 * A*A*V1*(b*0.3*Er1 / (2 * r1))*(0.3*Er1 / (2 *
  r1)));

          J = emitter->GetParticleSource()->sourceSurface[i].currentDensity*(1 - omega) +
  omega*jnew;

          J = std::min(1.0e6,double(J));

          if (Er1 < 0)
                  J = 0;

          emitter->GetParticleSource()->sourceSurface[i].currentDensity = J;


          k++;

          baseOld = base;
          Jold = J;
  };*/
};