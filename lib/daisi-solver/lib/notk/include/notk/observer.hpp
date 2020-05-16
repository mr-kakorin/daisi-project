#ifndef NOTK_OBSERVER_HPP
#define NOTK_OBSERVER_HPP

#include <string>

namespace notk
{
class NOTKObserver
{
  public:
    virtual ~NOTKObserver() = default;

    virtual void handle_event(const std::string& status_message) const = 0;
};
}

#endif
