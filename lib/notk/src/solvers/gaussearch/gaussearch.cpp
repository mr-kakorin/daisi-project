#include <numeric>

#include "../src/tools/tools.h"
#include "gaussconfig.h"
#include "gaussearch.h"

#include <serreflection/read_json.hpp>

#include <notk/notkresult.hpp>

namespace notk
{

template class GaussSearch<double, double>;
template class GaussSearch<float, float>;

REGISTER_CHILD(GaussSearch<double COMMA double>, IOptimizationStep<double COMMA double>,
               "Gauss search")
REGISTER_CHILD(GaussSearch<float COMMA float>, IOptimizationStep<float COMMA float>, "Gauss search")

template <class Targ, class Tfit>
Tfit GaussSearch<Targ, Tfit>::fitness1dWrapper(const double x, const std::vector<Targ>& all,
                                               const size_t dim) const
{
    std::vector<Targ> current = all;

    current[dim] = x;

    return this->fitness(current);
}

template <class Targ, class Tfit>
Tfit GaussSearch<Targ, Tfit>::single_iteration(const std::vector<size_t>& indexes,
                                               std::vector<Targ>&     x_current,
                                               const borders_t<Targ>& current_borders) const
{
    size_t problem_dim = indexes.size();
    NOTKResults<Targ, Tfit> result_tmp;
    for (size_t dim_i = 0; dim_i < problem_dim; dim_i++)
    {
        const size_t dim = indexes[dim_i];

        m_searcher1d->set_fitness(
            std::bind(&GaussSearch::fitness1dWrapper, this, std::placeholders::_1, x_current, dim));

        borders_t<Targ> borders_1d = std::make_pair(std::vector<Targ>{current_borders.first[dim]},
                                                    std::vector<Targ>{current_borders.second[dim]});

        std::vector<Targ> x_current_tmp = {x_current[dim]};
        bool              flag_abort    = true;

        m_searcher1d->opt_routine(result_tmp, x_current_tmp, borders_1d, 1, true, flag_abort,
                                  false);

        x_current[dim] = x_current_tmp[0];
    }
    return this->fitness(x_current);
}

template <class Targ, class Tfit>
it_res_t<Targ, Tfit> GaussSearch<Targ, Tfit>::do_preprocess(const it_res_t<Targ, Tfit>& iter_result,
                                                            const borders_t<Targ>& current_borders)
{
    const size_t problem_dim = current_borders.first.size();

    std::vector<std::vector<size_t>> initial_indexes =
        tools::get_array_of_rand_indexes(problem_dim, config()->get_max_preprocess_iterations());

    std::vector<Tfit>              start_fits(initial_indexes.size() + 1);
    std::vector<std::vector<Targ>> start_x(initial_indexes.size() + 1);

    start_fits[0] = iter_result.second;
    start_x[0]    = iter_result.first;

    for (size_t start_iter = 0; start_iter < initial_indexes.size(); start_iter++)
    {
        start_x[start_iter + 1] = iter_result.first;
        start_fits[start_iter + 1] =
            single_iteration(initial_indexes[start_iter], start_x[start_iter + 1], current_borders);
    }

    size_t n_best =
        std::distance(start_fits.begin(), std::min_element(start_fits.begin(), start_fits.end()));

    it_res_t<Targ, Tfit> result = std::make_pair(start_x[n_best], this->fitness(start_x[n_best]));

    return result;
}

template <class Targ, class Tfit>
it_res_t<Targ, Tfit> GaussSearch<Targ, Tfit>::do_iteration(const size_t iter_counter,
                                                           const it_res_t<Targ, Tfit>& iter_result,
                                                           const borders_t<Targ>& current_borders)
{
    const size_t problem_dim = current_borders.first.size();

    // use random
    std::vector<size_t> new_indexes(problem_dim);
    std::iota(new_indexes.begin(), new_indexes.end(), 0);

    size_t old_last_index = 0;

    // use non-random
    if (config()->get_use_random())
    {
        new_indexes = tools::get_rand_indexes(problem_dim);

        // не стоит проводить оптимизацию дважды по одному и тому же измерению
        if (old_last_index == new_indexes[0])
        {
            std::swap(new_indexes[0], new_indexes.back());
        }
        old_last_index = new_indexes.back();
    }
    std::vector<Targ> x_best = iter_result.first;
    single_iteration(new_indexes, x_best, current_borders);
    it_res_t<Targ, Tfit> result = std::make_pair(x_best, this->fitness(x_best));
    return result;
}

template <class Targ, class Tfit>
borders_t<Targ> GaussSearch<Targ, Tfit>::squeez_borders(const size_t                iter_counter,
                                                        const it_res_t<Targ, Tfit>& iter_result,
                                                        const borders_t<Targ>& current_borders)
{
    return current_borders;
}
template <class Targ, class Tfit>
bool GaussSearch<Targ, Tfit>::read_config(const boost::property_tree::ptree& config)
{
    this->m_config = srfl::read_json<GaussSearchConfig>(config);

    try
    {
        auto pt_searcher1d = config.get_child("searcher1d");

        m_searcher1d = IOptimizationStep<Targ, Tfit>::create_step(pt_searcher1d);

        m_searcher1d->read_config(pt_searcher1d);

        return this->m_config && m_searcher1d->config() ? true : false;
    }
    catch (const std::exception& ex)
    {
        BL_ERROR() << "Error read searcher 1d config: " << ex.what();
        return false;
    }
}
template <class Targ, class Tfit>
std::shared_ptr<GaussSearchConfig> GaussSearch<Targ, Tfit>::config() const
{
    return std::static_pointer_cast<GaussSearchConfig>(this->m_config);
}
}
