//
// Created by artoria on 3/15/20.
//

#ifndef DAISI_CLI_ITASKSCREATOR_H
#define DAISI_CLI_ITASKSCREATOR_H

#include <iostream>
#include <daisi-solver/project.h>
#include <daisi-solver/Results.h>
#include <algorithm>
#include <sstream>
#include "DaisiResult.h"

using CalculationTask = std::function<void()>;

struct boundary_values_t
{
    double emitter;
    double gateway;
    double anode;
};

class ITasksCreator {
public:
	ITasksCreator( std::string &&project_path, std::string &&energy_distribution_type,
                   double timeout_ns, unsigned int num_cores )
		: _num_cores( num_cores )
		, _project_path( std::move( project_path ) )
		, _energy_distribution_type( std::move( energy_distribution_type ) )
		, _timeout_ns( timeout_ns )
	{}
	virtual ~ITasksCreator() = default;

    unsigned int                 get_num_cores() const { return _num_cores; }
    std::vector<CalculationTask> operator()();
protected:
    CalculationTask     boundary_to_task( boundary_values_t const& boundary_values );
    virtual void        log_status( std::string const& code, boundary_values_t const& boundary_values, std::string const& path, double timeout ) = 0;
    virtual std::string fname_prefix( boundary_values_t const& boundary_values ) = 0;
    virtual void        set_boundaries( Dproject::project& proj, boundary_values_t const& boundary_values ) = 0;
    void                set_timeout( Dproject::project& proj, double timeout );

    virtual std::vector<boundary_values_t> get_boundaries() = 0;

    unsigned int _num_cores;
    std::string _project_path;
    std::string _energy_distribution_type;
    double      _timeout_ns;
};


class TasksFlowFactory {
public:
	static std::vector<CalculationTask> createTasksFlow( ITasksCreator *tasks_creator ) {
		return tasks_creator->operator()();
	}
};

#endif //DAISI_CLI_ITASKSCREATOR_H
