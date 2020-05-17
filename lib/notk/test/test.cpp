#include <gtest/gtest.h>

#include <common_tools/boostlog.hpp>

namespace logging = boost::log;

int main(int argc, char** argv)
{
    commtools::BoostLog::get(commtools::sev_lvl::debug, "test");
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
