#include <iostream>
#include <sys/stat.h>
#include "taskflow/taskflow.hpp"
#include "ITasksCreator.h"
#include "config.h"
#include "validate.h"

int solver( int argc, char *argv[] ) {
	std::string config_path( "./config.json" );
	if (argc > 1) {
		config_path = argv[1];
	}

	mkdir( DaisiEnv::result_path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	auto tasks_creator = read_config( config_path.c_str());
	if (!tasks_creator)
		return 1;
	auto num_cores = tasks_creator->get_num_cores();
	auto tasks = TasksFlowFactory::createTasksFlow( tasks_creator );
	tf::Taskflow taskflow;
	for (auto &task: tasks) {
		taskflow.emplace( task );
	}
	tf::Executor executor( num_cores );
	executor.run( taskflow );
	executor.wait_for_all();
	return 0;
}

int main( int argc, char *argv[] ) {
	solver(argc, argv);
	return 0;
}
