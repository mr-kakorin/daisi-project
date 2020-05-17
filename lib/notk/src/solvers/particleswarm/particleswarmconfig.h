
#ifndef NOTK_PARTICLE_SWARM_CONFIG
#define NOTK_PARTICLE_SWARM_CONFIG

#include "../../base/iconfig.h"

namespace notk
{
class ParticleSwarmConfig final : public BaseOptConfig
{
  public:
    size_t get_n_agents();
    bool   get_is_gauss();

    size_t n_agents;
    bool   is_gauss;
};
}

SERIALIZIBLE_STRUCT(notk::ParticleSwarmConfig, srfl::CheckModes::FATAL,
                    SER_BASE()(size_t, n_agents, srfl::nan, 1.0, srfl::inf)(bool, is_gauss, DEF_D()))

#endif
