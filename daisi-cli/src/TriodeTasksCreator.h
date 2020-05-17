//
// Created by nikita on 5/18/20.
//

#ifndef DAISI_CLI_TRIODETASKSCREATOR_H
#define DAISI_CLI_TRIODETASKSCREATOR_H

#include "ITasksCreator.h"

class TriodeTasksCreator
        : public ITasksCreator {
public:
    TriodeTasksCreator( std::vector<double> &&anode_V, std::vector<double> &&emitter_V, std::vector<double> &&gateway_V,
                        std::string &&project_path,
                        double timeout_ns, unsigned int num_cores  )
            : ITasksCreator( std::move(project_path), timeout_ns, num_cores )
            , _anode_V( std::move( anode_V ) )
            , _emitter_V( std::move( emitter_V ) )
            , _gateway_V( std::move( gateway_V ) )
    {}
    ~TriodeTasksCreator() override = default;
protected:
    /*virtual*/ std::vector<boundary_values_t> get_boundaries() override;
    /*virtual*/ void                           set_boundaries( Dproject::project& proj, boundary_values_t const& boundary_values ) override;
    /*virtual*/ void                           log_status( std::string const& code, boundary_values_t const& boundary_values,
                                                           std::string const& path, double timeout ) override;
    /*virtual*/ std::string                    fname_prefix( boundary_values_t const& boundary_values ) override;
private:
    std::vector<double> _anode_V;
    std::vector<double> _emitter_V;
    std::vector<double> _gateway_V;
};

#endif //DAISI_CLI_TRIODETASKSCREATOR_H
