#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

#include <armadillo>
#include <gtest/gtest.h>

#include <notk/notkbuilder.hpp>

#define STR1(x) #x
#define STR(x) STR1(x)

void get_gauss_borders(const arma::vec& x, const arma::vec& y, std::vector<double>& x1,
                       std::vector<double>& x2)
{
    x1.resize(3);
    x2.resize(3);

    x1[0] = arma::min(x);
    x2[0] = arma::max(x);

    x1[1] = 1e-9 * (arma::max(x) - arma::min(x));
    x2[1] = 1.1 * (arma::max(x) - arma::min(x));

    x1[2] = 0;
    x2[2] = 3 * arma::max(y);
}

template <class T>
T gaussmf(const T& x, const double mu, const double sigma, const double a)
{
    if (std::abs(sigma) < 1e-20)
    {
        throw std::runtime_error("gaussmf::too small sigma");
    }
    T result(x.n_rows);
    for (size_t i = 0; i < x.n_rows; i++)
    {
        result(i) = a * std::exp(-(x(i) - mu) * (x(i) - mu) / (2.0 * sigma * sigma));
    }
    return result;
}

double fitness(const arma::vec& x, const arma::vec& y, const std::vector<double>& param)
{
    arma::vec y_ = gaussmf(x, param[0], param[1], param[2]);
    return arma::sum((y_ - y) % (y_ - y));
}

class TestGaussFit : public ::testing::Test
{
  protected:
    virtual void SetUp()
    {
        path = STR(TEST_PATH);

        x_1.load(path + "/gauss_fit/gauss_x_1", arma::raw_binary);
        y_1.load(path + "/gauss_fit/gauss_y_1", arma::raw_binary);

        x_2.load(path + "/gauss_fit/gauss_x_2", arma::raw_binary);
        y_2.load(path + "/gauss_fit/gauss_y_2", arma::raw_binary);

        get_gauss_borders(x_1, y_1, border_x1_1, border_x2_1);
        get_gauss_borders(x_2, y_2, border_x1_2, border_x2_2);

        fitness_1 = std::bind(fitness, x_1, y_1, std::placeholders::_1);
        fitness_2 = std::bind(fitness, x_2, y_2, std::placeholders::_1);
    }
    std::string path;

    arma::vec x_1;
    arma::vec y_1;
    arma::vec x_2;
    arma::vec y_2;

    std::vector<double> border_x1_1;
    std::vector<double> border_x2_1;
    std::vector<double> border_x1_2;
    std::vector<double> border_x2_2;

    std::function<double(const std::vector<double>&)> fitness_1;
    std::function<double(const std::vector<double>&)> fitness_2;

    bool flag_abort = true;

    void test(const std::string& config, const std::array<double, 4>& expected, const double tol)
    {
        auto NOTKController = notk::NOTKBuilder<double, double>::build_notk_file_config(
            config, fitness_1, border_x1_1, border_x2_1);

        ASSERT_TRUE(NOTKController);

        auto result_1 = NOTKController->process(flag_abort);

        ASSERT_TRUE(result_1);

        std::cout << " tol = " << tol << std::endl;

        for (size_t i = 0; i < 3; i++)
        {
            EXPECT_NEAR(expected[i], result_1->get_last_argument()[i], expected[i] * tol);
        }

        NOTKController->set_borders_fitness(fitness_2, border_x1_2, border_x2_2);

        auto result_2 = NOTKController->process(flag_abort);

        ASSERT_TRUE(result_2);

        EXPECT_NEAR(expected[3], result_2->get_last_argument()[0], expected[3] * tol);
    }
};

TEST_F(TestGaussFit, gauss_multidim)
{
    test(path + "/gauss_fit/gauss_test.json",
         std::array<double, 4>{21.847052993004247, 0.4917746155735323, 262.27086413032117,
                               18.022617945262549},
         0);
}

TEST_F(TestGaussFit, particle_swarm)
{

    test(path + "/gauss_fit/swarm_test.json",
         std::array<double, 4>{21.838229576388134, 0.48405155140051387, 260.83346944643898,
                               18.025435901975236},
         5e-2);
}

TEST_F(TestGaussFit, heavy_ball)
{
    test(path + "/gauss_fit/heavyball_test.json",
         std::array<double, 4>{21.674007564354305, 0.36300851549203667, 368.93273202212077,
                               18.645157712406714},
         0);
}

TEST_F(TestGaussFit, swarm_and_ball)
{
    test(path + "/gauss_fit/swarm&ball_test.json",
         std::array<double, 4>{21.804476702281104, 0.466282315234399, 270.95529147101951,
                               18.022042385320702},
         5e-2);
}
