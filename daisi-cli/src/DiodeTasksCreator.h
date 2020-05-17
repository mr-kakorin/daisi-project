//
// Created by nikita on 5/18/20.
//

#ifndef DAISI_CLI_DIODETASKSCREATOR_H
#define DAISI_CLI_DIODETASKSCREATOR_H

#include "ITasksCreator.h"

class DiodeTasksCreator
        : public ITasksCreator {
public:
    DiodeTasksCreator( std::vector<double> &&anode_V, std::vector<double> &&emitter_V, std::string &&project_path,
                       double timeout_ns, unsigned int num_cores  )
            : ITasksCreator( std::move(project_path), timeout_ns, num_cores )
            , _anode_V( std::move(anode_V) )
            , _emitter_V( std::move(emitter_V) )
    {}
    ~DiodeTasksCreator() override = default;
protected:
    /*virtual*/ std::vector<boundary_values_t> get_boundaries() override;
    /*virtual*/ void                           set_boundaries( Dproject::project& proj, boundary_values_t const& boundary_values ) override;
    /*virtual*/ void                           log_status( std::string const& code, boundary_values_t const& boundary_values,
                                                           std::string const& path, double timeout ) override;
    /*virtual*/ std::string                    fname_prefix( boundary_values_t const& boundary_values ) override;

    std::vector<double> _anode_V;
    std::vector<double> _emitter_V;
};

#endif //DAISI_CLI_DIODETASKSCREATOR_H
