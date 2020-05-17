#include "TriodeTasksCreator.h"

/*virtual*/ std::vector<boundary_values_t> TriodeTasksCreator::get_boundaries()
{
    std::vector<boundary_values_t> boundaries;
    boundaries.reserve( _emitter_V.size() * _anode_V.size() );
    for (auto &e_v: _emitter_V) {
        for (auto &a_v: _anode_V) {
            for (auto &g_v: _gateway_V) {
                boundaries.emplace_back(boundary_values_t{e_v, g_v, a_v});
            }
        }
    }
    return boundaries;
}

/*virtual*/ void TriodeTasksCreator::set_boundaries( Dproject::project& proj, boundary_values_t const& boundary_values )
{
    proj.currentModel->SetConditionProperties( "poisson", 0, 0,
    {{boundary_values.emitter, 0., 0., 0.}} );
    proj.currentModel->SetConditionProperties( "poisson", 0, 1,
    {{boundary_values.gateway, 0., 0., 0.}} );
    proj.currentModel->SetConditionProperties( "poisson", 0, 2,
    {{boundary_values.anode, 0., 0., 0.}} );
}

/*virtual*/ void TriodeTasksCreator::log_status( std::string const& code, boundary_values_t const& boundary_values, std::string const& path, double timeout )
{
    std::stringstream ss;
    ss << code << " task: " << path << ", timeout_ns = " << timeout
    << ", e_v = " << boundary_values.emitter << ", a_v = " << boundary_values.anode << ", g_v = "<< boundary_values.gateway << std::endl;
    std::cout << ss.str();
}

/*virtual*/ std::string TriodeTasksCreator::fname_prefix( boundary_values_t const& boundary_values )
{
    int len = 5;
    return std::to_string( boundary_values.emitter ).substr( 0, len ) + "_" + std::to_string( boundary_values.gateway ).substr( 0, len ) +
           "_" + std::to_string( boundary_values.anode ).substr( 0, len );
}