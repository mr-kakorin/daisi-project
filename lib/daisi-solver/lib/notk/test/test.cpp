#include <gtest/gtest.h>

#include <boostlog0/boostlog.h>
#include <boostlog0/boostlogconf.h>

namespace logging = boost::log;

int main(int argc, char** argv)
{
    initboostlog::init(sev_lvl::info, "test");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
