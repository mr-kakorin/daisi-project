#include <notk/notkresult.hpp>

namespace notk
{
template class NOTKResults<double, double>;
template class NOTKResults<float, float>;

template <class Targ, class Tfit>
std::vector<Targ> NOTKResults<Targ, Tfit>::get_last_argument() const noexcept
{
    if (m_result.empty())
    {
        return std::vector<Targ>();
    }
    return m_result.back().first;
}
template <class Targ, class Tfit>
Tfit NOTKResults<Targ, Tfit>::get_fitness() const noexcept
{
    if (m_result.empty())
    {
        return 0;
    }
    return m_result.back().second;
}
template <class Targ, class Tfit>
it_res_t<Targ, Tfit> NOTKResults<Targ, Tfit>::get_last_it_res() const noexcept
{
    return m_result.back();
}
template <class Targ, class Tfit>
it_res_t<Targ, Tfit> NOTKResults<Targ, Tfit>::get_res(const size_t number) const noexcept
{
    return m_result[number];
}

template std::vector<float>  NOTKResults<double, double>::get_fit_array() const noexcept;
template std::vector<float>  NOTKResults<float, float>::get_fit_array() const noexcept;
template std::vector<double> NOTKResults<double, double>::get_fit_array() const noexcept;
template std::vector<double> NOTKResults<float, float>::get_fit_array() const noexcept;

template <class Targ, class Tfit>
template <class Tres>
std::vector<Tres> NOTKResults<Targ, Tfit>::get_fit_array() const noexcept
{
    std::vector<Tres> result;
    for (auto& val : m_result)
    {
        result.push_back(val.second);
    }
    return result;
}

template <class Targ, class Tfit>
std::vector<std::vector<Targ>> NOTKResults<Targ, Tfit>::get_arg_arr() const noexcept
{
    std::vector<std::vector<Targ>> result(m_result[0].first.size());

    for (auto& val : m_result)
    {
        for (size_t i = 0; i < val.first.size(); i++)
        {
            result[i].push_back(val.first[i]);
        }
    }
    return result;
}

template <class Targ, class Tfit>
const std::vector<it_res_t<Targ, Tfit>>& NOTKResults<Targ, Tfit>::result() const noexcept
{
    return m_result;
}
}