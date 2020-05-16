
#ifndef NOTK_GAUSS_CONFIG_H
#define NOTK_GAUSS_CONFIG_H

#include "../../base/iconfig.h"

namespace notk
{
class GaussSearchConfig : public BaseOptConfig
{
  public:
    bool   get_use_random();
    size_t get_max_preprocess_iterations();

    size_t maximal_preprocess_iterations;
    bool   use_random;
};
}

SERIALIZIBLE_STRUCT(notk::GaussSearchConfig, srfl::CheckModes::FATAL,
                    SER_BASE()(size_t, maximal_preprocess_iterations, srfl::nan, 1,
                               srfl::inf)(bool, use_random, DEF_D()))

#endif
