#include "DiodeTasksCreator.h"

/*virtual*/ std::vector<boundary_values_t> DiodeTasksCreator::get_boundaries()
{
    std::vector<boundary_values_t> boundaries;
    boundaries.reserve( _emitter_V.size() * _anode_V.size() );
    for (auto const&e_v: _emitter_V) {
        for (auto const&a_v: _anode_V) {
            boundaries.emplace_back( boundary_values_t{ e_v, 0.0, a_v } );
        }
    }
    return boundaries;
}

/*virtual*/ void DiodeTasksCreator::set_boundaries( Dproject::project& proj, boundary_values_t const& boundary_values )
{
    proj.currentModel->SetConditionProperties( "poisson", 0, 0, {{boundary_values.emitter, 0., 0., 0.}} );
    proj.currentModel->SetConditionProperties( "poisson", 0, 1, {{boundary_values.anode, 0., 0., 0.}} );
}

/*virtual*/ void DiodeTasksCreator::log_status( std::string const& code, boundary_values_t const& boundary_values, std::string const& path, double timeout )
{
    std::stringstream ss;
    ss << code << " task: " << path << ", timeout_ns = " << timeout
    << ", e_v = " << boundary_values.emitter << ", a_v = " << boundary_values.anode << std::endl;
    std::cout << ss.str();
}

/*virtual*/ std::string DiodeTasksCreator::fname_prefix( boundary_values_t const& boundary_values )
{
    return std::to_string( boundary_values.emitter ).substr( 0, 4 ) + "_" + std::to_string( boundary_values.anode ).substr( 0, 4 );
}