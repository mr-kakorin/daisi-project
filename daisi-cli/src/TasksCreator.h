//
// Created by artoria on 3/15/20.
//

#ifndef DAISI_CLI_TASKSCREATOR_H
#define DAISI_CLI_TASKSCREATOR_H

#include <iostream>
#include <daisi-solver/project.h>
#include <daisi-solver/Results.h>
#include <algorithm>
#include "DaisiResult.h"
#include "EnergyDistribution.h"

using CalculationTask = std::function<void()>;

class ITasksCreator {
public:
	ITasksCreator( unsigned int num_cores )
		: _num_cores( num_cores )
	{}
	virtual ~ITasksCreator() = default;
	virtual std::vector<CalculationTask> operator()() = 0;
	[[nodiscard]] unsigned int get_num_cores() const { return _num_cores; }
protected:
	unsigned int _num_cores;
};

//Diode is Base
class BaseTasksCreator
		: public ITasksCreator
{
public:
	BaseTasksCreator( std::vector<double> &&anode_V, std::vector<double> &&emitter_V, std::string &&project_path,
	                  double timeout_ns, unsigned int num_cores )
		: ITasksCreator( num_cores ),
		  _anode_V( std::move( anode_V )), _emitter_V( std::move( emitter_V )),
		  _project_path( std::move( project_path )),
		  _timeout_ns( timeout_ns )
	{}
	virtual ~BaseTasksCreator() override = default;
protected:
	std::vector<double> _anode_V;
	std::vector<double> _emitter_V;
	std::string _project_path;
	double _timeout_ns;
};

class DiodeTasksCreator
		: public BaseTasksCreator {
public:
	DiodeTasksCreator( std::vector<double> &&anode_V, std::vector<double> &&emitter_V, std::string &&project_path,
	                   double timeout_ns, unsigned int num_cores  )
			: BaseTasksCreator( std::move(anode_V), std::move(emitter_V), std::move(project_path), timeout_ns, num_cores )
	{}
	~DiodeTasksCreator() override = default;

	std::vector<CalculationTask> operator()() override {
		std::vector<CalculationTask> result;
		for (auto &e_v: _emitter_V) {
			for (auto &a_v: _anode_V) {
				result.emplace_back(
						[ e_v, a_v, project_path = this->_project_path, timeout_ns = this->_timeout_ns ]() {
							std::cout << "Starting task: " << project_path << ", timeout_ns = " << timeout_ns
							          << ", e_v = " << e_v << ", a_v = " << a_v << std::endl;
							Dproject::project currentProject;
							std::string errorMsgOpen;
							currentProject.LoadProject( project_path, DaisiEnv::version, errorMsgOpen );
							double progress = 0, progressLoc = 0;
							std::vector<std::string> status;
							bool flagAbort = true;
							std::string errorMsg;
							currentProject.currentModel->SetConditionProperties( "poisson", 0, 0, {{e_v, 0., 0., 0.}} );
							currentProject.currentModel->SetConditionProperties( "poisson", 0, 1, {{a_v, 0., 0., 0.}} );
							auto numbers = currentProject.currentModel->GetFlowMCNumbers( 0 );
							numbers[2] = timeout_ns;
							currentProject.currentModel->SetFlowMCNumbers( 0, numbers );
							currentProject.currentModel->setEnergyDistribution( &get_initial_energy_bimodal );
							currentProject.currentModel->SimulateCPUPIC( progress, progressLoc, status, flagAbort, 1,
							                                             errorMsg );
							std::cout << "Ending task: " << project_path << ", timeout_ns = " << timeout_ns
							          << ", e_v = " << e_v << ", a_v = " << a_v << std::endl;
							DaisiResult daisiResult;
							daisiResult.loadFlowData( &currentProject );
							daisiResult.loadElectrodeData( &currentProject );
							daisiResult.saveEmitterField(  &currentProject,
							                               std::to_string( e_v ).substr( 0, 4 ) + "_" +
							                               std::to_string( a_v ).substr( 0, 4 ) );
							daisiResult.saveCoordinatesData( &currentProject,
							                          std::to_string( e_v ).substr( 0, 4 ) + "_" +
							                          std::to_string( a_v ).substr( 0, 4 ));
							daisiResult.saveFieldNormData( &currentProject,
							                          std::to_string( e_v ).substr( 0, 4 ) + "_" +
							                          std::to_string( a_v ).substr( 0, 4 ));
							daisiResult.saveFlowData( &currentProject,
							                          std::to_string( e_v ).substr( 0, 4 ) + "_" +
							                          std::to_string( a_v ).substr( 0, 4 ));
							daisiResult.saveElectrodeData( &currentProject,
							                               std::to_string( e_v ).substr( 0, 4 ) + "_" +
							                               std::to_string( a_v ).substr( 0, 4 ));
						} );
			}
		}
		return result;
	}
};

class TriodeTasksCreator
		: public BaseTasksCreator {
public:
	TriodeTasksCreator( std::vector<double> &&anode_V, std::vector<double> &&emitter_V, std::vector<double> &&gateway_V,
	                    std::string &&project_path,
	                    double timeout_ns, unsigned int num_cores  )
			: BaseTasksCreator( std::move(anode_V), std::move(emitter_V), std::move(project_path), timeout_ns, num_cores ),
			  _gateway_V( std::move( gateway_V ) )
			  {}
	~TriodeTasksCreator() override = default;

	std::vector<CalculationTask> operator()() override {
		std::vector<CalculationTask> result;
		for (auto &e_v: _emitter_V) {
			for (auto &a_v: _anode_V) {
				for (auto &g_v: _gateway_V) {
					result.emplace_back(
							[ e_v, a_v, g_v, project_path = this->_project_path, timeout_ns = this->_timeout_ns ]() {
								std::cout << "Starting task: " << project_path << ", timeout_ns = " << timeout_ns
								          << ", e_v = " << e_v << ", a_v = " << a_v << ", g_v = "<< g_v << std::endl;
								Dproject::project currentProject;
								std::string errorMsgOpen;
								currentProject.LoadProject( project_path, DaisiEnv::version, errorMsgOpen );
								double progress = 0, progressLoc = 0;
								std::vector<std::string> status;
								bool flagAbort = true;
								std::string errorMsg;
								currentProject.currentModel->SetConditionProperties( "poisson", 0, 0,
								                                                     {{e_v, 0., 0., 0.}} );
								currentProject.currentModel->SetConditionProperties( "poisson", 0, 1,
								                                                     {{g_v, 0., 0., 0.}} );
								currentProject.currentModel->SetConditionProperties( "poisson", 0, 2,
								                                                     {{a_v, 0., 0., 0.}} );
								auto numbers = currentProject.currentModel->GetFlowMCNumbers( 0 );
								numbers[2] = timeout_ns;
								currentProject.currentModel->SetFlowMCNumbers( 0, numbers );
								currentProject.currentModel->setEnergyDistribution( &get_initial_energy_bimodal );
								currentProject.currentModel->SimulateCPUPIC( progress, progressLoc, status, flagAbort,
								                                             1,
								                                             errorMsg );
								std::cout << "Ending task: " << project_path << ", timeout_ns = " << timeout_ns
								          << ", e_v = " << e_v << ", a_v = " << a_v << std::endl;
								DaisiResult daisiResult;
								daisiResult.loadFlowData( &currentProject );
								daisiResult.loadElectrodeData( &currentProject );
								daisiResult.saveCoordinatesData(&currentProject,
								                                std::to_string( e_v ).substr( 0, 4 ) + "_" +
								                                std::to_string( g_v ).substr( 0, 4 ) + "_" +
								                                std::to_string( a_v ).substr( 0, 4 ));
								daisiResult.saveEmitterField(&currentProject,
								                            std::to_string( e_v ).substr( 0, 4 ) + "_" +
								                            std::to_string( g_v ).substr( 0, 4 ) + "_" +
								                            std::to_string( a_v ).substr( 0, 4 ));
								daisiResult.saveFieldNormData(&currentProject,
								                                std::to_string( e_v ).substr( 0, 4 ) + "_" +
								                                std::to_string( g_v ).substr( 0, 4 ) + "_" +
								                                std::to_string( a_v ).substr( 0, 4 ));
								daisiResult.saveFlowData( &currentProject,
								                          std::to_string( e_v ).substr( 0, 4 ) + "_" +
								                          std::to_string( g_v ).substr( 0, 4 ) + "_" +
								                          std::to_string( a_v ).substr( 0, 4 ));
								daisiResult.saveElectrodeData( &currentProject,
								                               std::to_string( e_v ).substr( 0, 4 ) + "_" +
								                               std::to_string( g_v ).substr( 0, 4 ) + "_" +
								                               std::to_string( a_v ).substr( 0, 4 ));
							} );
				}
			}
		}
		return result;
	}

private:
	std::vector<double> _gateway_V;
};


class TasksFlowFactory {
public:
	static std::vector<CalculationTask> createTasksFlow( ITasksCreator *tasks_creator ) {
		return tasks_creator->operator()();
	}
};

#endif //DAISI_CLI_TASKSCREATOR_H
