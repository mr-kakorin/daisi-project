#include "ITasksCreator.h"

CalculationTask ITasksCreator::boundary_to_task( boundary_values_t const& boundary_values )
{
    return [ this, project_path = this->_project_path, timeout = this->_timeout_ns, boundary_values ]()
    {
        log_status("Starting", boundary_values, project_path, timeout );
        Dproject::project currentProject;
        std::string errorMsg;
        currentProject.LoadProject( project_path, DaisiEnv::version, errorMsg );
        set_boundaries( currentProject, boundary_values );
        set_timeout( currentProject, timeout );
        double progress = 0, progressLoc = 0; std::vector<std::string> status; bool flagAbort = true;

        currentProject.currentModel->setEnergyDistribution( &get_initial_energy_bimodal );
        currentProject.currentModel->SimulateCPUPIC( progress, progressLoc, status, flagAbort, 1,
                                                     errorMsg );

        log_status("Ending", boundary_values, project_path, timeout );

        DaisiResult daisiResult( currentProject );
        daisiResult.save( fname_prefix( boundary_values ) );
    };
}

std::vector<CalculationTask> ITasksCreator::operator ()()  {
    std::vector<CalculationTask> result;
    auto boundaries = get_boundaries();
    std::transform( boundaries.begin(), boundaries.end(), std::back_inserter( result ),
                    std::bind( &ITasksCreator::boundary_to_task, this, std::placeholders::_1 ) );
    return result;
}

void ITasksCreator::set_timeout( Dproject::project& proj, double timeout )
{
    auto numbers = proj.currentModel->GetFlowMCNumbers( 0 );
    numbers[2] = timeout;
    proj.currentModel->SetFlowMCNumbers( 0, numbers );
}