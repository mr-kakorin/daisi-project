#include <gtest/gtest.h>

#include <notk/controller.hpp>
#include <notk/testfunctions.hpp>

#define STR1(x) #x
#define STR(x) STR1(x)

class Statemachine : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
        path = STR(TEST_PATH);
    }
    std::string path;
    bool        flag_abort = true;
};

TEST_F(Statemachine, StateMachineTest1)
{
    notk::NOTKController<double, double> NOTKController;

    ASSERT_FALSE(NOTKController.process(flag_abort));
    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2}));
    ASSERT_FALSE(NOTKController.set_fitness(notk::Rastrigin));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/test_func/gauss_test.json"));
    ASSERT_TRUE(NOTKController.set_problem_config(path + "/test_func/gauss_test.json"));

    ASSERT_FALSE(NOTKController.process(flag_abort));
    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2}));

    ASSERT_TRUE(NOTKController.set_fitness(notk::Rastrigin));

    ASSERT_TRUE(NOTKController.process(flag_abort));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2}));

    ASSERT_FALSE(NOTKController.process(flag_abort));

    ASSERT_TRUE(NOTKController.set_fitness(notk::Rastrigin));

    ASSERT_TRUE(NOTKController.process(flag_abort));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/test_func/gauss_test.json"));

    ASSERT_FALSE(NOTKController.process(flag_abort));
    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2}));
}

TEST_F(Statemachine, StateMachineTest2)
{
    notk::NOTKController<double, double> NOTKController;

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(NOTKController.process(flag_abort));
    ASSERT_TRUE(NOTKController.set_fitness(notk::Rastrigin));
    ASSERT_FALSE(NOTKController.process(flag_abort));

    ASSERT_TRUE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2}));
    ASSERT_TRUE(NOTKController.process(flag_abort));
    ASSERT_FALSE(NOTKController.process(flag_abort));
    ASSERT_TRUE(NOTKController.set_fitness(notk::Rastrigin));
    ASSERT_FALSE(NOTKController.process(flag_abort));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1}, {2, 2}));
    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_TRUE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 2}, 2));

    ASSERT_TRUE(NOTKController.process(flag_abort));
    ASSERT_TRUE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 2}, 2));
    ASSERT_TRUE(NOTKController.process(flag_abort));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {2, 2, 3}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1}, {1, 1}, {2, 2, 3}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1, 1}, {1, 1}, {2, 2, 3}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_TRUE(
        NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1, 1}, {3, 3, 3}, {2, 2, 2}));

    ASSERT_TRUE(NOTKController.process(flag_abort));
}

TEST_F(Statemachine, StateMachineTest3)
{
    notk::NOTKController<double, double> NOTKController;

    ASSERT_FALSE(NOTKController.set_problem_config(path + ""));

    ASSERT_FALSE(NOTKController.set_problem_config_str("{}"));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(
        NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1, 1}, {3, 3, 3}, {4, 2, 2}));

    ASSERT_FALSE(
        NOTKController.set_borders_fitness(notk::Rastrigin, {1, 1, 1}, {3, 3, 3}, {2, 2, 2}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(
        NOTKController.set_borders_fitness(notk::Rastrigin, {5, 1, 1}, {3, 3, 3}, {2, 2, 2}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(
        NOTKController.set_borders_fitness(notk::Rastrigin, {1, 6, 1}, {3, 3, 3}, {4, 4, 4}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_TRUE(NOTKController.set_borders_fitness(notk::Rastrigin, {3, 3, 3}, {4, 4, 4}));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {3, 3, 5}, {4, 4, 4}));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {3, 3, 3}, {4, 4, 4}));

    ASSERT_TRUE(NOTKController.set_problem_config(path + "/gauss_fit/gauss_test.json"));

    ASSERT_FALSE(NOTKController.set_borders_fitness(notk::Rastrigin, {3, 3, 5}, {4, 4, 4}));
}
