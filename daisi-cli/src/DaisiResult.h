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
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include "vtkBrush.h"
#include "vtkContext2D.h"
#include "vtkContextActor.h"
#include "vtkContextItem.h"
#include <vtkLookupTable.h>
#include <vtkPolyDataMapper.h>
#include <vtkGL2PSExporter.h>
#include "vtkContextScene.h"
#include "vtkObjectFactory.h"
#include "vtkOpenGLContextDevice2D.h"
#include "vtkPen.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkStdString.h"
#include "vtkTextProperty.h"
using dynamicsDaisiData = std::vector<std::vector<std::vector<std::shared_ptr<DynamicsData>>>>;
using simulationDaisiData = std::vector<std::shared_ptr<SimulationData>>;

namespace colors
{
	const float                     whitecolor[3] = {255, 255, 255};
	const float                     redcolor[3]   = {1.0, 0.0, 0.0};
	const float                     greencolor[3] = {0.0, 1.0, 0.0};
	const float                     blackcolor[3] = {0.0, 0.0, 0.0};
	const float                     bluecolor[3]  = {0.0, 0.0, 1.0};
	const std::vector<const float*> colors        = {redcolor, greencolor, bluecolor};
	const std::vector<const float*> colorsPoints  = {blackcolor, greencolor, bluecolor};
}

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
		return "/home/artoria/results/" + time_of_start + "/";
	}
public :
	inline const static std::string version = "1.2.0";
	inline const static std::string result_path = get_result_path();
};

struct DaisiPlot {
	int nSimData;
	int plotType;
	int nFlow;
	int nPlot;
	std::string yName;
};

class DaisiResult {
public:
	static inline std::mutex vtk_mutex{};

    DaisiResult( Dproject::project& dProject )
        : proj( dProject )
    {
	    double background[]         = {255, 255, 255};
	    renderer = vtkSmartPointer<vtkRenderer>::New() ;
	    renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	    renderWindow->AddRenderer(renderer);
	    renderer->SetBackground(background);
        dynamicsDaisiData dynamicData = proj.currentModel->GetDynamicsData();
        simulationDaisiData simData = proj.currentModel->GetSimulationData();
        for (int nd = 0; nd < dynamicData.size(); nd++) {
            for (int k1 = 0; k1 < simData[nd]->YDataElectrode.size(); k1++) {
                for (int k2 = 0; k2 < simData[nd]->YDataElectrode[k1].size(); k2++) {
                    _electrodePlots.emplace_back( DaisiPlot{ nd, 1, k1, k2, simData[nd]->dataFlags[1][k2] } );
                }
            }
            for (int k1 = 0; k1 < simData[nd]->YDataFlow.size(); k1++) {
                //std::cout << "Flow number " << k1 << std::endl;
                for (int k2 = 0; k2 < simData[nd]->YDataFlow[k1].size(); k2++) {
                    _flowPlots.emplace_back( DaisiPlot{ nd, 0, k1, k2, simData[nd]->dataFlags[0][k2] } );
                }
            }
        }
    }

    void save( std::string const& fname_prefix )
    {
        saveCoordinatesData( fname_prefix );
        saveFieldNormData( fname_prefix );
        saveElectrodeData( fname_prefix );
        saveFlowData( fname_prefix );
        saveEmitterField( fname_prefix );
        saveParticlesPos( fname_prefix );
        saveParticlesMom( fname_prefix );
        saveParticlesQ( fname_prefix );
	    saveVisualizationDataMesh( fname_prefix );
	    saveCloud( fname_prefix );
	    saveMesh ( fname_prefix );
    }

private:
	void saveCoordinatesData( const std::string &p ) {
		void* ArrayX[1];
		void* ArrayY[1];
		int   size;
		int   elemSize;
		proj.currentModel->GetGridData(ArrayX, size, elemSize, 0, 3);
		proj.currentModel->GetGridData(ArrayY, size, elemSize, 1, 3);
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

	void saveFieldNormData( const std::string &p ) {
		void* ArrayX[1];
		void* ArrayY[1];
		void* ArrayE[1];
		int   size;
		int   elemSize;
		proj.currentModel->GetGridData(ArrayX, size, elemSize, 0, 0);
		proj.currentModel->GetGridData(ArrayY, size, elemSize, 1, 0);
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

	void saveEmitterField( const std::string &p ) {
		auto result = proj.currentModel->GetEmitterField( 0 );
		std::ofstream outfile( DaisiEnv::result_path + p + "EmitterField" + std::string( ".csv" ));
		auto& curve = result[0];
		auto& E = result[1];
		auto size = result[0].size();
		for (auto i = 0; i < size; ++i){
			outfile << curve[i] << ',' << E[i] << std::endl;
		}
		outfile.close();
	}

	void saveElectrodeData( const std::string &p ) {
		for (int i = 0; i < _electrodePlots.size(); ++i) {
			std::ofstream outfile( DaisiEnv::result_path + p + _electrodePlots.at( i ).yName + std::string( ".csv" ));
			int nsd = _electrodePlots.at( i ).nSimData;
			int nf = _electrodePlots.at( i ).nFlow;
			int np = _electrodePlots.at( i ).nPlot;
			auto datax = proj.currentModel->GetSimulationData()[nsd]->XData;
			auto datay = proj.currentModel->GetSimulationData()[nsd]->YDataElectrode[np][0];
			int resLength = std::min( datax.size(), datay.size());
			outfile << _electrodePlots.at( i ).yName << std::endl;
			for (int j = 0; j < resLength; ++j) {
				outfile << datax[j] << ',' << datay[j] << ',' << std::endl;
			}
			outfile.close();
		}
	}

	void saveFlowData( const std::string &p ) {
		for (int i = 0; i < _flowPlots.size(); ++i) {
			std::ofstream outfile( DaisiEnv::result_path + p + _flowPlots.at( i ).yName + std::string( ".csv" ));
			int nsd = _flowPlots.at( i ).nSimData;
			int nf = _flowPlots.at( i ).nFlow;
			int np = _flowPlots.at( i ).nPlot;
			auto datax = proj.currentModel->GetSimulationData()[nsd]->XData;
			auto datay = proj.currentModel->GetSimulationData()[nsd]->YDataFlow[nf][np];
			int resLength = std::min( datax.size(), datay.size());
			outfile << _flowPlots.at( i ).yName << std::endl;
			for (int j = 0; j < resLength; ++j) {
				outfile << datax[j] << ',' << datay[j] << ',' << std::endl;
			}
			outfile.close();
		}
	}

	void saveParticlesPos( const std::string &p)
	{
		auto pos = proj.currentModel->get_particles_positions();
		std::ofstream outfile( DaisiEnv::result_path + p + "Positions" + std::string( ".csv" ));
		for (int i = 0; i < pos.size(); ++i){
			for (int j = 0; j< pos[i].size(); ++j)
			outfile << pos[i][j] << ',' << std::endl;
		}
		outfile.close();
	}

	void saveParticlesMom( const std::string &p)
	{
		auto mom = proj.currentModel->get_particles_moments();
		std::ofstream outfile( DaisiEnv::result_path + p + "Momentums" + std::string( ".csv" ));
		for (int i = 0; i < mom.size(); ++i){
			for (int j = 0; j< mom[i].size(); ++j)
				outfile << mom[i][j] << ',' << std::endl;
		}
		outfile.close();
	}

	void saveParticlesQ( const std::string &p)
	{
		auto q = proj.currentModel->get_particles_charge();
		std::ofstream outfile( DaisiEnv::result_path + p + "Charge" + std::string( ".csv" ));
		for (int i = 0; i < q.size(); ++i) {
			outfile << q[i] << ',' << std::endl;
		}
		outfile.close();
	}

	void saveCloud( const std::string &p) {
		std::vector<int> list = proj.currentModel->GetBoundariesList();
		std::ofstream boundaries_x( DaisiEnv::result_path + p + "Boundaries_x" + std::string( ".csv" ));
		std::ofstream boundaries_y( DaisiEnv::result_path + p + "Boundaries_y" + std::string( ".csv" ));
		for (int i = 0; i < list.size(); i++)
		{
			std::vector<float> X;
			std::vector<float> Y;
			proj.currentModel->GetPlotXYBoundary( list[i], X, Y );
			for (int i=0;i<X.size(); ++i)
				boundaries_x << X[i] << ',';
			for (int i=0;i<Y.size(); ++i)
				boundaries_y << Y[i] << ',';
			boundaries_x<<std::endl;
			boundaries_y<<std::endl;
		}

		std::vector<int>                arraySize;
		int                             elemSize;
		std::vector<std::vector<void*>> pointArray;
		proj.currentModel->GetMomentumsParticlesCloud( 0, 0, pointArray, arraySize, elemSize );
		std::ofstream cloud( DaisiEnv::result_path + p + "ParticlesCloud" + std::string( ".csv" ));

		for ( int thread = 0; thread < pointArray.size(); thread++ )
		{
			for ( int i = 0; i < arraySize[thread]; i++ )
			{
				for ( int k = 0; k < pointArray[thread].size(); k++ )
					cloud << ( *((double*)((char*)pointArray[thread][k] + i * elemSize)) ) << ',';
				cloud << std::endl;
			}
		}
	}

	void saveMesh( std::string const& p )
	{
		std::vector<std::vector<std::pair<double,double>>> meshData = proj.currentModel->get_mesh();
		std::ofstream mesh_x( DaisiEnv::result_path + p + "Mesh_x" + std::string( ".csv" ));
		std::ofstream mesh_y( DaisiEnv::result_path + p + "Mesh_y" + std::string( ".csv" ));
		for (int i=0;i<meshData.size();++i) {
			for (int j = 0; j < meshData[i].size(); ++j)
				mesh_x << meshData[i][j].first << ',';
			mesh_x << std::endl;
		}
		for (int i=0;i<meshData.size();++i) {
			for (int j = 0; j < meshData[i].size(); ++j)
				mesh_y << meshData[i][j].second << ',';
			mesh_y << std::endl;
		}
	}

	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkRenderWindow> renderWindow;
	void saveVisualizationDataMesh( std::string const& fname_prefix )
	{
		std::lock_guard<std::mutex> lock( DaisiResult::vtk_mutex );
		void* Array[1];
		int   size;
		int   elemSize;
		proj.currentModel->GetGridData(Array, size, elemSize, "Charge density, cl/m^3", 0);
		vtkSmartPointer<vtkFloatArray> vtkData = vtkSmartPointer<vtkFloatArray>::New();

		vtkData->SetName("Charge density, cl/m^3");

		vtkUnstructuredGrid*                 VTKgrid    = proj.currentModel->GetVTKGrid();
		vtkSmartPointer<vtkUnstructuredGrid> newVTKgrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
		newVTKgrid->ShallowCopy(VTKgrid);

		double minVal = *((double*)((char*)Array[0]));
		double maxVal = *((double*)((char*)Array[0]));
		double valtmp;
		for (int i = 1; i < size; i++)
		{
			valtmp = *((double*)((char*)Array[0] + i * elemSize));
			if (valtmp < minVal && !std::isnan(valtmp) && !std::isinf(valtmp))
				minVal = valtmp;
			if (valtmp > maxVal && !std::isnan(valtmp) && !std::isinf(valtmp))
				maxVal = valtmp;
		}
		double newmaxVal = maxVal - minVal;
		for (int i = 0; i < size; i++)
		{
			valtmp = (*((double*)((char*)Array[0] + i * elemSize)) - minVal) / newmaxVal;
			vtkData->InsertNextValue(valtmp);
		}

		newVTKgrid->GetPointData()->SetScalars(vtkData);

		vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
		mapper->SetInputData(newVTKgrid);

		actor->SetMapper(mapper);
		actor->GetProperty()->SetColor(0.0, 1.0, 0.0);

		renderer->AddActor(actor);

		vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
		poly->DeepCopy(newVTKgrid);
		poly->GetPointData()->SetScalars(vtkData);

		vtkSmartPointer<vtkPolyDataMapper> poly_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

		poly_mapper->SetInputData(poly);
		poly_mapper->ScalarVisibilityOn();
		poly_mapper->SetScalarModeToUsePointData();
		poly_mapper->SetColorModeToMapScalars();

		vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
		scalarBar->SetLookupTable(mapper->GetLookupTable());
		scalarBar->SetTitle("Charge density, cl/m^3");
		scalarBar->GetTitleTextProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1], colors::blackcolor[2]);
		scalarBar->SetNumberOfLabels(10);
		scalarBar->GetLabelTextProperty()->SetColor(colors::blackcolor[0], colors::blackcolor[1], colors::blackcolor[2]);

		vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
		hueLut->SetTableRange(minVal, maxVal);
		hueLut->Build();

		poly_mapper->SetLookupTable(hueLut);
		scalarBar->SetLookupTable(hueLut);

		renderer->AddActor2D(scalarBar);

		vtkSmartPointer<vtkGL2PSExporter> vtext;
		vtext = vtkSmartPointer<vtkGL2PSExporter>::New();
		vtext->SetFileFormatToEPS();
		vtext->SetLineWidthFactor(1.5);
		vtext->SetRenderWindow(renderWindow);
		vtext->SetFilePrefix( ( DaisiEnv::result_path + fname_prefix + "_ChargeDensity").c_str());
		vtext->Write();
	};


	Dproject::project& proj;
    std::vector<DaisiPlot> _flowPlots;
    std::vector<DaisiPlot> _electrodePlots;
};
#endif //DAISI_CLI_DAISIRESULT_H
