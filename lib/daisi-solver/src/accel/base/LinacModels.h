#ifndef LinacModels_H
#define LinacModels_H

//#include "DataTypes.h"
#include "LinacFlow.h"

/*class RFQDesignModel : public ModelTemplate <device3ddouble>
{

        friend class boost::serialization::access;
public:
        std::vector<double> RFQParameters;
        std::vector<double> CellsLengths;
        std::vector<double> MinumalRadiiRegular;
        std::vector<double> MatcherRadii;
        std::vector<double> AvEnergies;
        std::vector<double> AccEff;

        std::vector<double> MinumalRadii;

        std::vector<double> Modulations;
        std::vector<double> SyncPhases;
        std::vector<double> OutRadii;

        std::vector<std::vector<double>> RFQApproxParameters;

        std::vector<LinacFlow> LinacFlows;
        std::vector<DynamicsData*> outputDataLinac;
        LinacTools::LinacDynSimulator DynSimulator;

        BOOST_SERIALIZATION_SPLIT_MEMBER()
        template<class Archive>
        void save(Archive & ar, const unsigned int) const
        {
                ar & boost::serialization::base_object<ModelTemplate>(*this);
                ar & RFQParameters;
                ar & LinacFlows;
                ar & CellsLengths;
                ar & MinumalRadii;
                ar & DynSimulator;
                ar & MatcherRadii;
                ar & MinumalRadiiRegular;
                ar & AvEnergies;
                ar & AccEff;
                ar & Modulations;
                ar & SyncPhases;
                ar & RFQApproxParameters;
                ar & OutRadii;

        }
        template<class Archive>
        void load(Archive & ar, const unsigned int)
        {
                ar & boost::serialization::base_object<ModelTemplate>(*this);
                ar & RFQParameters;
                ar & LinacFlows;
                ar & CellsLengths;
                ar & MinumalRadii;
                ar & DynSimulator;
                ar & MatcherRadii;
                ar & MinumalRadiiRegular;
                ar & AvEnergies;
                ar & AccEff;
                ar & Modulations;
                ar & SyncPhases;
                ar & RFQApproxParameters;
                ar & OutRadii;

                RFQParameters.resize(2);
                RFQApproxParameters.resize(6);
                RFQApproxParameters[5].resize(2);
        }
        void RFQEmOpt(bool& flagAbort, double& progress, std::vector<double> params);
        void RFQAccOpt(bool& flagAbort, double& progress, std::vector<double> params);
        void RFQMatcherOpt(bool& flagAbort, double& progress, std::vector<double> params);
        void SetLinacApproximationParameters(std::vector<std::vector<double>> in)
        {
                RFQApproxParameters = in;
        };
        std::vector<std::vector<double>> GetLinacApproximationParameters()
        {
                return RFQApproxParameters;
        };
        void ExportToLidos(std::string filename, int currentFlow);

        void SearchStartParameters(bool& flagAbort, double& progress, std::vector<double> params);

        std::vector<DynamicsData*>  GetDynamicsDataLinac()
        {
                return outputDataLinac;
        };
        void CalculateRFQAcc(int currentFlow);

        void SetLinacDynSimParameters(std::vector<double> in)
        {
                DynSimulator.SetParameters(in);
        };


        std::vector<double> GetLinacDynSimParameters()
        {
                return DynSimulator.GetParameters();
        };


        void SimulateLinacFlowDynamics(int flow);

        std::vector<double>&  GetRFQCavityParameters(int flag)
        {
                if (flag==0)
                        return CellsLengths;
                if (flag == 1)
                        return MinumalRadiiRegular;
                if (flag == 2)
                        return MatcherRadii;
                if (flag == 3)
                        return MinumalRadii;
                if (flag == 4)
                        return AvEnergies;
                if (flag == 5)
                        return AccEff;
        };
        void GenerateRFQForFlow(int& succes, int flow);


        std::vector<double>  GetLinacFlowParameters(int flow)
        {
                return LinacFlows[flow].GetParametersRFQ();
        };
        void  SetLinacFlowParameters(int flow, std::vector<double> in)
        {
                LinacFlows[flow].SetParametersRFQ(in);
        };


        RFQDesignModel();


        std::vector<double>& GetControl(int flag);


        void SetControlFunction(std::string file, int flag);


        void SetLinacParameters(std::vector<double> in, std::vector<std::vector<double>> inParams,
std::string folder);




        std::vector<double> GetLinacParameters();


        void AddLinacFlow();

        int GetNumberOfLinacFlows();


        void GenerateMesh(std::string meshParam, double& progress);


        void AddFlow(int ParticleType, int DistributionStyle, double massNumber, double
chargeNumber);


        void AddFlow(int ParticleType, int DistributionStyle);

};*/

#endif