/*









void ttt(const double* t)
{

};
freq::freq(std::string parFile)
{
double a = 2;
ttt(&a);
FILE* fp;
char ss[250];
fopen_s(&fp, parFile.c_str(), "r");
double tmp1, tmp2, tmp3, tmp4, tmp5;

while (fgets(ss, 250, fp))
{
int tmpi;
sscanf_s(ss, "%lf %lf %lf %lf %lf", &tmp1, &tmp2, &tmp3, &tmp4, &tmp5);
Lcell.push_back(tmp1);
kList.push_back(tmp2);
Rleft.push_back(tmp3);
Rright.push_back(tmp4);
frequency.push_back(tmp5);
};
fclose(fp);

LDifarray.push_back(Lcell[0]);

for (int i = 1; i < Lcell.size(); i++)
{
if (Lcell[i] != Lcell[i - 1])
{
int flag = 0;
for (int j = 0; j < LDifarray.size(); j++)
{
if (Lcell[i] == LDifarray[j])
{
flag = 1;
break;
}
};
if (flag == 0)
LDifarray.push_back(Lcell[i]);
};
};
std::sort(LDifarray.begin(), LDifarray.end());

};
void freq::InterpFreq(double Linput, double Rleftinput, double Rrightinput, double frequencyInpt,
double LgapEstimate,
std::vector<double>& out)
{
int  flag = 0;
int k = 0;
{
while (1)
{
if (k >= LDifarray.size())
break;

if (LDifarray[k] < Linput && LDifarray[k + 1] > Linput)
{
flag = 1;
break;
}
k++;
}
}

if (Linput < LDifarray[0])
{
k = 0;
flag = 1;
};

if (flag == 0)
{
out[0] = LgapEstimate / Linput;
out[1] = Rrightinput;
return;
}

double L1 = LDifarray[k];
double L2 = LDifarray[k + 1];



double w1 = (L2 - Linput) / (L2 - L1);
double w2 = 1 - w1;

int j = 0;
std::vector <double> kInterpList1;
std::vector <double> kInterpList2;

std::vector <double> freqList1;
std::vector <double> freqList2;

for (int ii = 0; ii < Lcell.size(); ii++)
{
if (std::abs(Rleftinput - Rleft[ii])<1e-8 && std::abs(Rrightinput - Rright[ii])<1e-8 && std::abs(L1 -
Lcell[ii])<1e-8)
{
kInterpList1.push_back(kList[ii]);
freqList1.push_back(frequency[ii]);
}
}

for (int ii = 0; ii < Lcell.size(); ii++)
{
if (std::abs(Rleftinput - Rleft[ii])<1e-8 && std::abs(Rrightinput - Rright[ii])<1e-8 &&  std::abs(L2 -
Lcell[ii])<1e-8)
{
kInterpList2.push_back(kList[ii]);
freqList2.push_back(frequency[ii]);
}
}

std::vector <double> kInterpList;
std::vector <double> freqList;


for (int i = 0; i < kInterpList1.size(); i++)
{
for (int j = 0; j < kInterpList2.size(); j++)
{
if (kInterpList1[i] == kInterpList2[j])
{
kInterpList.push_back(kInterpList1[i]);
freqList.push_back(freqList1[i] * w1 + freqList2[j] * w2);
};
}
}
if (freqList.size() == 0)
return;
//	freqList.push_back(frequency[k1 + ii] * w1 + frequency[k + ii + 1] * w2);

double  Res = Dmath::Interpolate(freqList, kInterpList, frequencyInpt);


out[0] = Res;
out[1] = Rrightinput;


};
void freq::GetCellParameters(double Linput, double Rleftinput, double frequencyInpt, double
LgapEstimate,
std::vector<double>& out)
{

InterpFreq(Linput, Rleftinput, Rleftinput, frequencyInpt, LgapEstimate, out);

if (out[0] * Linput < LgapEstimate)
InterpFreq(Linput, Rleftinput, Rleftinput + 0.001, frequencyInpt, LgapEstimate, out);

if (out[0]>0.7 && Rleftinput>0.014)
InterpFreq(Linput, Rleftinput, Rleftinput - 0.001, frequencyInpt, LgapEstimate, out);

if (out[0] > 0.75)
out[0] = 0.75;

out[0] = out[0] * Linput;
};















std::vector<double> GetAcceleratorParams()
{
        return LinacParams;
}
void SetBoundaries(double particleMomentum, double& dt, double& L, int period)
{
        double frequency = LinacParams[0];
        double initEnergy = LinacParams[1];
        double factorL1 = LinacParams[2];
        double factorL2 = LinacParams[12];

        double dtInnerRadius1 = LinacParams[3];
        double dtOuterRadius1 = LinacParams[4];
        double dtInnerRadius2 = LinacParams[5];
        double dtOuterRadius2 = LinacParams[6];
        double accInnerRadius = LinacParams[7];
        int numberOfPeriods = LinacParams[11];

        double dtInnerRadius = dtInnerRadius1 + (double(period) /
double(numberOfPeriods))*(dtInnerRadius2 -
dtInnerRadius1); double dtOuterRadius = dtOuterRadius1 + (double(period) /
double(numberOfPeriods))*(dtOuterRadius2 -
dtOuterRadius1);

        double dtInnerRadius11 = dtInnerRadius1 + (double(period + 1) /
double(numberOfPeriods))*(dtInnerRadius2 -
dtInnerRadius1); double dtOuterRadius11 = dtOuterRadius1 + (double(period + 1) /
double(numberOfPeriods))*(dtOuterRadius2 - dtOuterRadius1);

        double factorL = factorL1 + (double(period) / double(numberOfPeriods))*(factorL2 -
factorL1);


        boundaries.clear();
        PointType restEnergy;
        PointType gamma;
        PointType beta;
        double lambda = commtools::LIGHT_VELOCITY() / frequency;
        dt = lambda / 1000.0;
        beta = particleMomentum / sqrt(1 + particleMomentum*particleMomentum);

        if (L == -1)
        {
                if (LinacParams[13] == 0)
                        L = beta*lambda;
                if (LinacParams[13] == 1)
                        L = beta*lambda / 2;
        }

        double down;

        if (LinacParams[14] == 0)
                down = (1 - factorL) / 2 * L;

        if (LinacParams[14] == 1)
                down = (L - factorL) / 2;

        //down = (1 - factorL) / 2 * L;
        double up = L - down;
        std::vector<std::vector<double>> dtDown = { { dtInnerRadius, dtInnerRadius, dtOuterRadius,
dtOuterRadius },{ 0,
down, down, 0 } }; std::vector<std::vector<double>> dtUp = { { dtOuterRadius11, dtOuterRadius11,
dtInnerRadius11,
dtInnerRadius11 },{ L, up, up, L } }; std::vector<std::vector<double>> accLeft = { { dtInnerRadius,
0, 0,
dtInnerRadius11 },{ 0, 0, L, L } }; std::vector<std::vector<double>> accRigth = { { dtOuterRadius,
accInnerRadius,
accInnerRadius, dtOuterRadius11 },{ 0, 0, L, L, } };

        //std::vector<std::vector<double>> vertex = { { 0, dtInnerRadius, dtInnerRadius,
dtOuterRadius, dtOuterRadius,
accInnerRadius, accInnerRadius, dtOuterRadius, dtOuterRadius, dtInnerRadius, dtInnerRadius, 0 }, {
0, 0, down, down, 0,
0, L, L, up, up, L, L } }; boundaries.push_back(BoundaryContainer2d<PointType>(accLeft));
        boundaries.push_back(BoundaryContainer2d<PointType>(accRigth));
        boundaries.push_back(BoundaryContainer2d<PointType>(dtDown));
        boundaries.push_back(BoundaryContainer2d<PointType>(dtUp));

        //boundariesForFlows.push_back(bound);
        //boundariesForFlows.back().RemoveFineEdges(0.0001);
}


void SetBoundaries(int nCells, const std::vector<double>& L, const std::vector<double>& Lgaps, const
std::vector<double>& TubesInner, const std::vector<double>& TubesOuter)
{

        boundaries.clear();

        double dtInnerRadius;
        double dtOuterRadius;
        double accInnerRadius = LinacParams[6];

        double roundingRadius = LinacParams[5];


        double Lsum = 0;
        std::vector<std::vector<double>> dtOdd(2);//�������� ������ ������ ���� � 1�
        std::vector<std::vector<double>> dtEven(2);//������
        boundariesForFlows.clear();

        double down = (L[0] - Lgaps[0]) / 2;

        dtOdd[0].push_back(TubesInner[0]);
        dtOdd[0].push_back(TubesInner[0]);
        dtOdd[0].push_back(TubesOuter[0]);
        dtOdd[0].push_back(TubesOuter[0]);
        dtOdd[1].push_back(0);
        dtOdd[1].push_back(down);
        dtOdd[1].push_back(down);
        dtOdd[1].push_back(0);

        double up = down;
        Lsum += L[0];

        for (int period = 1; period < nCells; period++)
        {
                down = (L[period] - Lgaps[period]) / 2;

                if (period & 1)
                {
                        dtEven[0].push_back(TubesInner[period]);
                        dtEven[0].push_back(TubesInner[period]);
                        dtEven[0].push_back(TubesOuter[period]);
                        dtEven[0].push_back(TubesOuter[period]);
                        dtEven[1].push_back(Lsum - up);
                        dtEven[1].push_back(Lsum + down);
                        dtEven[1].push_back(Lsum + down);
                        dtEven[1].push_back(Lsum - up);
                }
                else
                {
                        dtOdd[0].push_back(TubesInner[period]);
                        dtOdd[0].push_back(TubesInner[period]);
                        dtOdd[0].push_back(TubesOuter[period]);
                        dtOdd[0].push_back(TubesOuter[period]);
                        dtOdd[1].push_back(Lsum - up);
                        dtOdd[1].push_back(Lsum + down);
                        dtOdd[1].push_back(Lsum + down);
                        dtOdd[1].push_back(Lsum - up);
                }
                up = down;
                Lsum += L[period];
        }

        //Lsum += L[p2];
        std::vector<std::vector<double>> accLeft = { { TubesInner[0], 0, 0, TubesInner[nCells] },{
0, 0, Lsum, Lsum } };
        std::vector<std::vector<double>> accRigth = { { TubesOuter[nCells], accInnerRadius,
accInnerRadius,
TubesOuter[0] },{ Lsum, Lsum, 0, 0 } }; boundaries.push_back(BoundaryContainerType(accLeft, 0,
roundingRadius));
        boundaries.push_back(BoundaryContainerType(accRigth, 0, roundingRadius));
        boundariesForFlows.push_back(BoundaryContainerType(accLeft, 0, roundingRadius));
        boundariesForFlows.push_back(BoundaryContainerType(accRigth, 0, roundingRadius));


        if (nCells & 1)//��������� ��������� ��������� ����� - �������� �� ����������. ����������
��������� ������ �
����� ��������� ������ � ����� ����� �� ��������� - ������ ��� ��������
        {
                dtEven[0].push_back(TubesOuter[nCells]);
                dtEven[0].push_back(TubesOuter[nCells]);
                dtEven[0].push_back(TubesInner[nCells]);
                dtEven[0].push_back(TubesInner[nCells]);
                dtEven[1].push_back(Lsum);
                dtEven[1].push_back(Lsum - up);
                dtEven[1].push_back(Lsum - up);
                dtEven[1].push_back(Lsum);
                boundaries.push_back(BoundaryContainerType(dtOdd, 1, 0, roundingRadius));
                boundaries.push_back(BoundaryContainerType(dtEven, 0, 1, roundingRadius));
                boundariesForFlows.push_back(BoundaryContainerType(dtOdd, 1, 0, roundingRadius));
                boundariesForFlows.push_back(BoundaryContainerType(dtEven, 1, 0, roundingRadius));

        }
        else
        {
                dtOdd[0].push_back(TubesOuter[nCells]);
                dtOdd[0].push_back(TubesOuter[nCells]);
                dtOdd[0].push_back(TubesInner[nCells]);
                dtOdd[0].push_back(TubesInner[nCells]);
                dtOdd[1].push_back(Lsum);
                dtOdd[1].push_back(Lsum - up);
                dtOdd[1].push_back(Lsum - up);
                dtOdd[1].push_back(Lsum);
                boundaries.push_back(BoundaryContainerType(dtOdd, 1, 1, roundingRadius));
                boundaries.push_back(BoundaryContainerType(dtEven, 0, 0, roundingRadius));
                boundariesForFlows.push_back(BoundaryContainerType(dtOdd, 1, 1, roundingRadius));
                boundariesForFlows.push_back(BoundaryContainerType(dtEven, 0, 0, roundingRadius));
        }

        //boundariesForFlows.push_back(bound);
        //boundariesForFlows.back().RemoveFineEdges(0.0001);
}


void SetBoundariesMany(double particleMomentum, double& dt, std::vector<double>& L)
{
        double frequency = LinacParams[0];
        double initEnergy = LinacParams[1];
        double factorL1 = LinacParams[2];
        double factorL2 = LinacParams[12];
        double dtInnerRadius1 = LinacParams[3];
        double dtOuterRadius1 = LinacParams[4];
        double dtInnerRadius2 = LinacParams[5];
        double dtOuterRadius2 = LinacParams[6];
        double accInnerRadius = LinacParams[7];
        int numberOfPeriods = LinacParams[11];
        double dtInnerRadius;
        double dtOuterRadius;
        double factorL;
        double roundingRadius = LinacParams[13];

        boundaries.clear();
        PointType restEnergy;
        PointType gamma;
        PointType beta;
        double lambda = commtools::LIGHT_VELOCITY() / frequency;
        dt = lambda / 1000.0;
        beta = particleMomentum / sqrt(1 + particleMomentum*particleMomentum);
        double Lcurr;
        double down;
        double up = 0;
        double Lsum = 0;
        //std::vector<std::vector<double>> dtOdd(2);//�������� ������ ������ ���� � 1�
        //std::vector<std::vector<double>> dtEven(2);//������
        std::vector<std::vector<double>> tubes(2);
        tubes[0].resize(4);
        tubes[1].resize(4);

        boundariesForFlows.clear();
        for (int period = 0; period < numberOfPeriods; period++)
        {
                factorL = factorL1 + (double(period) / double(numberOfPeriods))*(factorL2 -
factorL1);
                Lcurr = L[period];


                if (LinacParams[15] == 0)
                        down = (1 - factorL) / 2 * Lcurr;

                if (LinacParams[15] == 1)
                        down = (Lcurr - factorL) / 2;

                dtInnerRadius = dtInnerRadius1 + (double(period) /
double(numberOfPeriods))*(dtInnerRadius2 -
dtInnerRadius1); dtOuterRadius = dtOuterRadius1 + (double(period) /
double(numberOfPeriods))*(dtOuterRadius2 -
dtOuterRadius1); tubes[0] = { dtInnerRadius, dtInnerRadius, dtOuterRadius, dtOuterRadius }; tubes[1]
= { Lsum - up, Lsum
+ down, Lsum + down, Lsum - up }; up = down; Lsum += L[period]; if (period == 0)
                {
                        boundaries.push_back(BoundaryContainerType(tubes, 1, 0, roundingRadius));
                        boundariesForFlows.push_back(BoundaryContainerType(tubes, 1, 0,
roundingRadius));

                }
                else
                {
                        boundaries.push_back(BoundaryContainerType(tubes, 0, 0, roundingRadius));
                        boundariesForFlows.push_back(BoundaryContainerType(tubes, 0, 0,
roundingRadius));
                }
        }
        tubes[0] = { dtInnerRadius2, dtInnerRadius2, dtOuterRadius2, dtOuterRadius2 };
        tubes[1] = { Lsum, Lsum - up, Lsum - up, Lsum };
        boundaries.push_back(BoundaryContainerType(tubes, 0, 1));
        boundariesForFlows.push_back(BoundaryContainerType(tubes, 0, 1));

        std::vector<std::vector<double>> accLeft = { { dtInnerRadius1, 0, 0, dtInnerRadius2 },{ 0,
0, Lsum, Lsum } };
        std::vector<std::vector<double>> accRigth = { { dtOuterRadius1, accInnerRadius,
accInnerRadius, dtOuterRadius2
},{ 0, 0, Lsum, Lsum } }; boundaries.push_back(BoundaryContainerType(accLeft, 1, 0,
roundingRadius));
        boundaries.push_back(BoundaryContainerType(accRigth, 1, 0, roundingRadius));

        boundariesForFlows.push_back(BoundaryContainerType(accLeft, 1, 0, roundingRadius));
        boundariesForFlows.push_back(BoundaryContainerType(accRigth, 1, 0, roundingRadius));


        //boundariesForFlows.back().RemoveFineEdges(0.0001);
}

void SetAcceleratorParams(std::vector<double> param, std::string filename)
{
        LinacParams = param;

        if (filename.size() == 0)
                return;

        FILE* fp;

        std::string filenameTmp;
        char ss[1000];

        double tmp;
        fopen_s(&fp, filename.c_str(), "r");
        SyncPhases.clear();


        while (fscanf(fp, "%lf", &tmp) > 0) {
                SyncPhases.push_back(tmp);
        }

        fclose(fp);
};

template  void  Solver<double>::GenerateGeometry  <device3ddouble>
(double& progress, bool& flagAbort, std::shared_ptr<device3ddouble>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<double>::GenerateGeometry  <device2daxsdouble>
(double& progress, bool& flagAbort, std::shared_ptr<device2daxsdouble>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void   Solver<float>::GenerateGeometry  <device2daxsfloat>
(double& progress, bool& flagAbort, std::shared_ptr<device2daxsfloat>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<double>::GenerateGeometry  <device2ddouble>
(double& progress, bool& flagAbort, std::shared_ptr<device2ddouble>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<float>::GenerateGeometry  <device2dfloat>
(double& progress, bool& flagAbort, std::shared_ptr<device2dfloat>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<double>::GenerateGeometry  <device2dpolardouble>
(double& progress, bool& flagAbort, std::shared_ptr<device2dpolardouble>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<float>::GenerateGeometry  <device2dpolarfloat>
(double& progress, bool& flagAbort, std::shared_ptr<device2dpolarfloat>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<double>::GenerateGeometry  <device3dExtrdouble>
(double& progress, bool& flagAbort, std::shared_ptr<device3dExtrdouble>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);

template  void  Solver<float>::GenerateGeometry  <device3dExtrfloat>
(double& progress, bool& flagAbort, std::shared_ptr<device3dExtrfloat>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder);



template <class PointType>
template <class deviceType>
void Solver<PointType>::GenerateGeometry(double& progress, bool& flagAbort,
std::shared_ptr<deviceType>&
deviceStatus,std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>&  outputData,
std::mutex& plotMutex,
std::shared_ptr<SimulationData>& simulationData, std::string meshParamFile, std::string
projectFolder)
{

std::string parFile = "..//frequencies.txt";

freq freParams(parFile);


if (deviceStatus->GetNumberParticlesFlows() == 0)
{
deviceStatus->resetFlows();
deviceStatus->AddFlow(1, 1, deviceStatus->LinacParams[8], deviceStatus->LinacParams[9]);
simulationData.dataFlags.resize(1);
std::vector<std::string> names = flagStringsSolver::simulationDataNamesFlowPIC;
for (int i = 0; i < names.size(); i++)
names[i] = names[i] + "0";
//	if (simulationData.dataFlags.size()<3)
simulationData.AddFlags(names);

std::vector<double> v(20);
v[1] = deviceStatus->LinacParams[1];
deviceStatus->GetFlow(0)->GetEmitterDevice()->SetDistributionParameters(v);
}



simulationData.reset();


std::vector<unsigned int> t;
W.resize(1);
for (int i = 0; i < 1; i++)
W[i] = new PointType[2][9];

double dt;
particleMover->InitParallel(1);
int numThreads = 1;

double frequency = deviceStatus->LinacParams[0];
double lambda = commtools::LIGHT_VELOCITY() / frequency;

dt = lambda / 1000.0;

std::vector<double> L(deviceStatus->LinacParams[10]);
std::vector<double> Lgaps(deviceStatus->LinacParams[10]);

std::vector<double> ROuter(deviceStatus->LinacParams[10]+1);
std::vector<double> RInner(deviceStatus->LinacParams[10] + 1);

std::vector<double> Ltmp(3);
std::vector<double> Lgapstmp(3);

std::vector<double> RInnerTmp = { deviceStatus->LinacParams[4], deviceStatus->LinacParams[4],
deviceStatus->LinacParams[4], deviceStatus->LinacParams[4] };

std::vector<double>  ROuterTmp(4);
ROuterTmp[0] = deviceStatus->LinacParams[3];
ROuter[0] = deviceStatus->LinacParams[3];
RInner[0] = deviceStatus->LinacParams[4];

std::vector<double> bS;


deviceStatus->GetFlow(0)->GenerateSyncParticle();

double periodMomentum = deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToMomentum2()[0];

double pr;
PointType beta1, beta2, beta3;
std::vector<std::vector<double>> cond;
double Lend;
std::vector<double> out(2);

for (int period = 0; period < deviceStatus->LinacParams[10]; period++)
{
L[period] = -1;
for (int iter = 0; iter < 1; iter++)
{
deviceStatus->GetFlow(0)->GetDynamicsData()->Time = 0;
deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToMomentum2()[0] = periodMomentum;


/*deviceStatus->SetBoundaries(deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToMomentum2()[0],
dt, L[period],
period); deviceStatus->GetboundaryConditions()->clear();

std::vector<std::vector<double>> cond;
if (deviceStatus->LinacParams[13] == 0)
cond = { { 0, 0, 0, 0 }, { deviceStatus->LinacParams[8], 0, 0, 0 } };
if (deviceStatus->LinacParams[13] == 1)
{
cond = { { 0, 0, 0, 0 }, { -deviceStatus->LinacParams[8], 0, 0, 0 } };
}

deviceStatus->GetboundaryConditions()->SetDefaultConditionsList({ 0, 1 });
for (int i = 0; i < 2; i++)
{
deviceStatus->GetboundaryConditions()->AddPropertyCondition(std::string("potential"), 0);
deviceStatus->GetboundaryConditions()->SetPropertyConditionsBoundariesList(i, { i + 2 });
deviceStatus->GetboundaryConditions()->SetConditionProperties(i, cond[i]);
}
deviceStatus->SetDomainBoundaryList({ 0, 1, 2, 3 });*/

/*
double phi1 = deviceStatus->SyncPhases[period];
double phi2 = deviceStatus->SyncPhases[period + 1];

double T;

if (deviceStatus->LinacParams[11] == 0)
T = commtools::LIGHT_VELOCITY()*(phi2 - phi1 + 2 * commtools::PI()) / (2 *
commtools::PI()*deviceStatus->LinacParams[0]);

if (deviceStatus->LinacParams[11] == 1)
T = commtools::LIGHT_VELOCITY()*(phi2 - phi1 + commtools::PI()) / (2 *
commtools::PI()*deviceStatus->LinacParams[0]);

double Tnorm = commtools::LIGHT_VELOCITY()*commtools::PI() / (2 *
commtools::PI()*deviceStatus->LinacParams[0]);


if (L[period] == -1)
{
        double lambda = commtools::LIGHT_VELOCITY() / frequency;
        double beta = periodMomentum / sqrt(1 + periodMomentum*periodMomentum);
        if (deviceStatus->LinacParams[11] == 0)
                L[period] = beta*lambda;
        if (deviceStatus->LinacParams[11] == 1)
                L[period] = (T / Tnorm)*beta*lambda / 2;
}

freParams.GetCellParameters(L[period], ROuter[period], frequency, deviceStatus->LinacParams[2],
out);

Lgaps[period] = out[0];
ROuter[period + 1] = out[1];
RInner[period + 1] = RInnerTmp[0];

if (period == 0)
{
        Ltmp[0] = L[period];
        Ltmp[1] = L[period];

        ROuterTmp[0] = ROuter[period];
        ROuterTmp[1] = ROuter[period + 1];
        ROuterTmp[2] = ROuter[period + 1];

        Lgapstmp[0] = Lgaps[period];
        Lgapstmp[1] = Lgaps[period];

        deviceStatus->SetBoundaries(2, Ltmp, Lgapstmp, RInnerTmp, ROuterTmp);
        deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToPosition2()[0] = 0;
        Lend = Ltmp[0];
        cond = { { 0, 0, 0, 0 },{ -deviceStatus->LinacParams[7], 0, 0, 0 } };
}
else
{
        Ltmp[0] = L[period - 1];
        Ltmp[1] = L[period];
        Ltmp[2] = L[period];

        ROuterTmp[0] = ROuter[period - 1];
        ROuterTmp[1] = ROuter[period];
        ROuterTmp[2] = ROuter[period + 1];
        ROuterTmp[3] = ROuter[period + 1];

        Lgapstmp[0] = Lgaps[period - 1];
        Lgapstmp[1] = Lgaps[period];
        Lgapstmp[2] = Lgaps[period];

        deviceStatus->SetBoundaries(3, Ltmp, Lgapstmp, RInnerTmp, ROuterTmp);
        deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToPosition2()[0] = Ltmp[0];
        Lend = Ltmp[0] + Ltmp[1];
        cond = { { -deviceStatus->LinacParams[7], 0, 0, 0 },{ 0, 0, 0, 0 } };
};





deviceStatus->GetboundaryConditions()->clear();
cond.resize(2);

deviceStatus->GetboundaryConditions()->SetDefaultConditionsList({ 0, 1 });
for (int i = 0; i < 2; i++)
{
        deviceStatus->GetboundaryConditions()->AddPropertyCondition(std::string("potential"), 0);
        deviceStatus->GetboundaryConditions()->SetPropertyConditionsBoundariesList(i, { i + 2 });
        deviceStatus->GetboundaryConditions()->SetConditionProperties(i, cond[i]);
}
deviceStatus->SetDomainBoundaryList({ 0, 1, 2, 3 });



meshGenerator->MeshGenerate(meshParamFile, pr, deviceStatus->GetDomainBoundary(),
deviceStatus->Getmesh(), 1);

//return;

deviceStatus->Getmesh()->Convert2GridData(deviceStatus->GetGridData());
deviceStatus->GetGridData()->densityReset();
InitFieldSolver(deviceStatus);

particleGridInterface->init(deviceStatus->GetGridData(), deviceStatus->Getmesh()->templNumb,
deviceStatus->Getmesh()->flagMatrix, {}, deviceStatus->GetDomainBoundary(), 4, 1, 4);
particleGridInterface->SearchStartCells(deviceStatus->GetFlow(0)->GetDynamicsData());

particleMover->flagInit = 1;
particleMover->SetTimeSteps({ dt });

int step = 0;
fieldSolver->FieldSimulate(deviceStatus->GetGridData(), deviceStatus->Getmesh(),
deviceStatus->GetboundaryConditions(),
0);



std::vector<double> betaS;

//	if (period==1)
//		break;

while (1)
{
        deviceStatus->GetGridData()->ApplyTimeDepending(deviceStatus->LinacParams[0], phi1 -
commtools::PI() / 2,
deviceStatus->GetFlow(0)->GetDynamicsData()->Time / (commtools::LIGHT_VELOCITY()));

        particleGridInterface->InCell(deviceStatus->GetFlow(0)->GetDynamicsData(), t, 0, 1, 0);

        particleGridInterface->Wcalculate(deviceStatus->GetFlow(0)->GetDynamicsData(), W[0], 0, 1,
0, 0);

        particleGridInterface->Grid2Particles(deviceStatus->GetFlow(0)->GetDynamicsData(),
deviceStatus->GetFlow(0)->GetDynamicsData(), deviceStatus->GetGridData(), W[0], 0, 1);

        particleMover->stepEstimate(deviceStatus->GetFlow(0)->GetDynamicsData(), 0, 0, 1, 0, 1);
        particleMover->updateMomentums(deviceStatus->GetFlow(0)->GetDynamicsData(),
deviceStatus->GetFlow(0)->GetDynamicsData(), deviceStatus->GetFlow(0)->GetAlpha(), 0, 0, 1, 0);
        deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToMomentum1()[0] = 0;
        particleMover->updatePositions(deviceStatus->GetFlow(0)->GetDynamicsData(),
deviceStatus->GetFlow(0)->GetDynamicsData(), 0, 0, 1, 0);
        deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToPosition1()[0] = 0.0001;


        deviceStatus->GetFlow(0)->GetDynamicsData()->Time =
deviceStatus->GetFlow(0)->GetDynamicsData()->Time +
particleMover->GetTimeStep(0, 0); deviceStatus->GetFlow(0)->GetDynamicsData()->flagEmitted[0] = 2;

        //deviceStatus->GetFlow(0)->GetDynamicsData()->GetBetaComponents(beta1, beta2, beta3, 0);
        //betaS.push_back(beta2);

        step++;
        if (deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToPosition2()[0] > Lend ||
deviceStatus->GetFlow(0)->GetDynamicsData()->cellsNumbers[0] == -1) break;
}
double tt = deviceStatus->GetFlow(0)->GetDynamicsData()->Time;
deviceStatus->GetFlow(0)->GetDynamicsData()->GetBetaComponents(beta1, beta2, beta3, 0);

double dL = beta2*(T - tt);
L[period] = L[period] + dL;
                }
                //	break;
                simulationData.addData(period, double(beta2));
                bS.push_back(beta2);
                progress = double(period) / double(deviceStatus->LinacParams[10]);
                periodMomentum =
deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToMomentum2()[0];
        }

        //	return;

        cond.clear();
        if (deviceStatus->LinacParams[11] == 0)
        {
                deviceStatus->SetBoundariesMany(deviceStatus->GetFlow(0)->GetDynamicsData()->GetPointerToMomentum2()[0],
dt, L); deviceStatus->GetboundaryConditions()->clear(); std::vector<int>
DomainBoundaryList(deviceStatus->LinacParams[10] + 3); int k = deviceStatus->LinacParams[10]; for
(int period = 0;
period < deviceStatus->LinacParams[10] + 1; period++)
                {
                        cond.push_back({ deviceStatus->LinacParams[7] * period, 0, 0, 0 });
                        deviceStatus->GetboundaryConditions()->SetDefaultConditionsList({ k + 1, k +
2 });
                        deviceStatus->GetboundaryConditions()->AddPropertyCondition(std::string("potential"),
0);
                        deviceStatus->GetboundaryConditions()->SetPropertyConditionsBoundariesList(period,
{ period });
                        deviceStatus->GetboundaryConditions()->SetConditionProperties(period,
cond[period]);
                        DomainBoundaryList[period] = period;
                }
                DomainBoundaryList[k + 1] = k + 1;
                DomainBoundaryList[k + 2] = k + 2;
                deviceStatus->SetDomainBoundaryList(DomainBoundaryList);
        }

        if (deviceStatus->LinacParams[11] == 1)
        {
                deviceStatus->SetBoundaries(deviceStatus->LinacParams[10], L, Lgaps, RInner,
ROuter);
                deviceStatus->GetboundaryConditions()->clear();
                cond.resize(2);
                cond = { { 0, 0, 0, 0 },{ -deviceStatus->LinacParams[7], 0, 0, 0 } };

                deviceStatus->GetboundaryConditions()->SetDefaultConditionsList({ 0, 1 });
                for (int i = 0; i < 2; i++)
                {
                        deviceStatus->GetboundaryConditions()->AddPropertyCondition(std::string("potential"),
0);
                        deviceStatus->GetboundaryConditions()->SetPropertyConditionsBoundariesList(i,
{ i + 2 });
                        deviceStatus->GetboundaryConditions()->SetConditionProperties(i, cond[i]);
                }
                deviceStatus->SetDomainBoundaryList({ 0, 1, 2, 3 });
        }

        std::string outFile = projectFolder + "//Linac.txt";

        FILE* fp1;

        fopen_s(&fp1, outFile.c_str(), "w");

        L.push_back(0);
        Lgaps.push_back(0);

        for (int i = 0; i < L.size(); i++)
        {
                fprintf(fp1, "%.8lf \t %.8lf \t %.8lf \t %.8lf \n", L[i], Lgaps[i], RInner[i],
ROuter[i]);
        };
        fclose(fp1);
        progress = 1;
}





*/
