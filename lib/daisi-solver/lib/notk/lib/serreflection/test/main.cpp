#include <gtest/gtest.h>

#include <boostlog0/boostlogconf.h>

int main(int argc, char** argv)
{
    initboostlog::init(sev_lvl::trace, "test");

    ::testing::InitGoogleTest(&argc, argv);

    return RUN_ALL_TESTS();
}
