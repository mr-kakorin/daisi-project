#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include <armadillo>
#include <gtest/gtest.h>

#include <notk/notkbuilder.hpp>
#include <notk/testfunctions.hpp>

#define STR1(x) #x
#define STR(x) STR1(x)

class Rastrigin : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
        path = STR(TEST_PATH);
    }
    void test(const std::string& config, const double val)
    {
        bool flag_abort = true;

        auto NOTKController =
            notk::NOTKBuilder<double, double>::build_notk_file_config(config, notk::Rastrigin);

        ASSERT_TRUE(NOTKController);

        auto result_1 = NOTKController->process(flag_abort);

        ASSERT_TRUE(result_1);

        EXPECT_TRUE(result_1->get_last_it_res().second <= val);
    }
    std::string path;
};

TEST_F(Rastrigin, gauss_multidim)
{
    test(path + "/test_func/gauss_test.json", 9.1532115220616106e-11);
}

TEST_F(Rastrigin, heavy_ball)
{
    test(path + "/test_func/heavyball_test.json", 8.5835267782385927e-09);
}

TEST_F(Rastrigin, particle_swarm)
{
    test(path + "/test_func/swarm_test.json", 1.05 * 0.00028709706293739146);
}

TEST_F(Rastrigin, swarm_and_ball)
{
    test(path + "/test_func/swarm&ball_test.json", 1.05 * 8.206999481785715e-06);
}
