//
// Created by artoria on 3/15/20.
//

#ifndef DAISI_CLI_DAISIRESULT_H
#define DAISI_CLI_DAISIRESULT_H

#include <vector>
#include <string>
#include <daisi-solver/project.h>
#include <daisi-solver/Results.h>
#include <memory>

using dynamicsDaisiData = std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>;
using simulationDaisiData = std::vector<std::shared_ptr<SimulationData>>;

class DaisiEnv {
private:
	static std::string get_result_path() {
		time_t rawtime;
		struct tm *timeinfo;
		char buffer[80];

		time( &rawtime );
		timeinfo = localtime( &rawtime );

		strftime( buffer, sizeof( buffer ), "%d-%m-%Y %H:%M:%S", timeinfo );
		std::string time_of_start( buffer );
		return "./" + time_of_start + "/";
	}
public :
	inline const static std::string version = "1.2.0";
	inline const static std::string result_path = get_result_path();
};

class DaisiPlot {
public:
	int nSimData;
	int plotType;
	int nFlow;
	int nPlot;
	std::string yName;

	DaisiPlot( int _nSimData, int _plotType, int _nFlow, int _nPlot, std::string _yName ) : nSimData( _nSimData ),
	                                                                                        plotType( _plotType ),
	                                                                                        nFlow( _nFlow ),
	                                                                                        nPlot( _nPlot ),
	                                                                                        yName( std::move(
			                                                                                        _yName )) {}

	std::string toString() {
		return "nSimData = " + std::to_string( this->nSimData ) + ", plotType = " +
		       std::to_string( this->plotType ) + ", nFlow = " + std::to_string( this->nFlow ) + ", nPlot = " +
		       std::to_string( this->nPlot ) + ", yName = " + this->yName + ";";
	}
};

class DaisiResult {
	std::vector<DaisiPlot> _flowPlots;
	std::vector<DaisiPlot> _electrodePlots;

public:
	void loadFlowData( Dproject::project *dProject ) {
		dynamicsDaisiData dynamicData = dProject->currentModel->GetDynamicsData();
		simulationDaisiData simData = dProject->currentModel->GetSimulationData();

		for (int nd = 0; nd < dynamicData.size(); nd++) {
			//std::cout << "Time plots" << std::endl;
			for (int k1 = 0; k1 < simData[nd]->YDataFlow.size(); k1++) {
				//std::cout << "Flow number " << k1 << std::endl;
				for (int k2 = 0; k2 < simData[nd]->YDataFlow[k1].size(); k2++) {
					_flowPlots.emplace_back( nd, 0, k1, k2, simData[nd]->dataFlags[0][k2] );
				}
			}
		}
	}

	void loadElectrodeData( Dproject::project *dProject ) {
		dynamicsDaisiData dynamicData = dProject->currentModel->GetDynamicsData();
		simulationDaisiData simData = dProject->currentModel->GetSimulationData();

		for (int nd = 0; nd < dynamicData.size(); nd++) {
			//std::cout << "Time plots" << std::endl;
			for (int k1 = 0; k1 < simData[nd]->YDataElectrode.size(); k1++) {
				//std::cout << "Flow number " << k1 << std::endl;
				for (int k2 = 0; k2 < simData[nd]->YDataElectrode[k1].size(); k2++) {
					_electrodePlots.emplace_back( nd, 1, k1, k2, simData[nd]->dataFlags[1][k2] );
				}
			}
		}
	}

	void saveCoordinatesData( Dproject::project *dProject, const std::string &p ) {
		void* ArrayX[1];
		void* ArrayY[1];
		int   size;
		int   elemSize;
		dProject->currentModel->GetGridData(ArrayX, size, elemSize, 0, 3);
		dProject->currentModel->GetGridData(ArrayY, size, elemSize, 1, 3);
		std::ofstream outfile( DaisiEnv::result_path + p + "Coordinates" + std::string( ".csv" ));

		for (int i = 0; i < size-1; ++i) {
				outfile << double(((double*)ArrayX[0])[i]) << ',';
		}
		outfile << double(((double*)ArrayX[0])[size-1]) << std::endl;
		for (int i = 0; i < size-1; ++i) {
			outfile << double(((double*)ArrayY[0])[i]) << ',';
		}
		outfile << double(((double*)ArrayY[0])[size-1]) << std::endl;
		outfile.close();
	}


	void saveFieldNormData( Dproject::project *dProject, const std::string &p ) {
		void* ArrayX[1];
		void* ArrayY[1];
		void* ArrayE[1];
		int   size;
		int   elemSize;
		dProject->currentModel->GetGridData(ArrayX, size, elemSize, 0, 0);
		dProject->currentModel->GetGridData(ArrayY, size, elemSize, 1, 0);
		ArrayE[0] = malloc(size*elemSize);
		std::ofstream outfile( DaisiEnv::result_path + p + "FieldNorm" + std::string( ".csv" ));
		for (int i = 0; i < size; i++)
		{
			double valtmp1                               = *((double*)((char*)ArrayX[0] + i * elemSize));
			double valtmp2                               = *((double*)((char*)ArrayY[0] + i * elemSize));
			*((double*)((char*)ArrayE[0] + i * elemSize)) = sqrt(valtmp1 * valtmp1 + valtmp2 * valtmp2);
		}
		for (int i = 0; i < size-1; ++i) {
			outfile << double(((double*)ArrayE[0])[i]) << ',';
		}
		outfile << double(((double*)ArrayE[0])[size-1]) << std::endl;
		outfile.close();
	}

	void saveEmitterField( Dproject::project *dProject, const std::string &p ) {
		auto result = dProject->currentModel->GetEmitterField( 0 );
		std::ofstream outfile( DaisiEnv::result_path + p + "EmitterField" + std::string( ".csv" ));
		auto& curve = result[0];
		auto& E = result[1];
		auto size = result[0].size();
		for (auto i = 0; i < size; ++i){
			outfile << curve[i] << ',' << E[i] << std::endl;
		}
		outfile.close();
	}


	void saveElectrodeData( Dproject::project *dProject, const std::string &p ) {
		for (int i = 0; i < _electrodePlots.size(); ++i) {
			std::ofstream outfile( DaisiEnv::result_path + p + _electrodePlots.at( i ).yName + std::string( ".csv" ));

			int nsd = _electrodePlots.at( i ).nSimData;
			int nf = _electrodePlots.at( i ).nFlow;
			int np = _electrodePlots.at( i ).nPlot;
			auto datax = dProject->currentModel->GetSimulationData()[nsd]->XData;
			auto datay = dProject->currentModel->GetSimulationData()[nsd]->YDataElectrode[np][0];
			int resLength = std::min( datax.size(), datay.size());
			outfile << _electrodePlots.at( i ).yName << std::endl;
			for (int j = 0; j < resLength; ++j) {
				outfile << datax[j] << ',' << datay[j] << ',' << std::endl;
			}
			outfile.close();
		}
	}

	void saveFlowData( Dproject::project *dProject, const std::string &p ) {
		for (int i = 0; i < _flowPlots.size(); ++i) {
			std::ofstream outfile( DaisiEnv::result_path + p + _flowPlots.at( i ).yName + std::string( ".csv" ));
			int nsd = _flowPlots.at( i ).nSimData;
			int nf = _flowPlots.at( i ).nFlow;
			int np = _flowPlots.at( i ).nPlot;
			auto datax = dProject->currentModel->GetSimulationData()[nsd]->XData;
			auto datay = dProject->currentModel->GetSimulationData()[nsd]->YDataFlow[nf][np];
			int resLength = std::min( datax.size(), datay.size());
			outfile << _flowPlots.at( i ).yName << std::endl;
			for (int j = 0; j < resLength; ++j) {
				outfile << datax[j] << ',' << datay[j] << ',' << std::endl;
			}
			outfile.close();
		}
	}
};
#endif //DAISI_CLI_DAISIRESULT_H
