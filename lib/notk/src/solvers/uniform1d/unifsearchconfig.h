#include "../../base/iconfig.h"

namespace notk
{
class UniformSearch1dConfig : public BaseOptConfig
{
  public:
    size_t get_n_divisions_first() const;

    size_t get_n_divisions() const;

    size_t n_divisions_first;
    size_t n_divisions;
};
}

SERIALIZIBLE_STRUCT(notk::UniformSearch1dConfig, srfl::CheckModes::FATAL,
                    SER_BASE()(size_t, n_divisions_first, srfl::nan, 2,
                               srfl::inf)(size_t, n_divisions, srfl::nan, 2, srfl::inf))
