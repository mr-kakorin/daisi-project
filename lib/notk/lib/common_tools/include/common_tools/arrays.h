
#ifndef COMMON_TOOLS_JSON_ARRAYS_H
#define COMMON_TOOLS_JSON_ARRAYS_H

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace commtools
{
enum class Operations
{
    MULT,
    PLUS,
    DIV,
    MINUS
};

template <class C>
constexpr auto size_(const C& c) -> decltype(c.size())
{
    return c.size();
}

#define applier(apply_inner, oper, type)                                                           \
    switch (oper)                                                                                  \
    {                                                                                              \
    case Operations::MULT:                                                                         \
        apply_inner(std::multiplies<type>());                                                      \
        break;                                                                                     \
    case Operations::PLUS:                                                                         \
        apply_inner(std::plus<type>());                                                            \
        break;                                                                                     \
    case Operations::DIV:                                                                          \
        apply_inner(std::divides<type>());                                                         \
        break;                                                                                     \
    case Operations::MINUS:                                                                        \
        apply_inner(std::minus<type>());                                                           \
        break;                                                                                     \
    default:                                                                                       \
        throw std::runtime_error("apply_op::unexpected operator");                                 \
    }

#define apply_inner_rval1(op)                                                                      \
    std::transform(std::begin(result), std::end(result), std::begin(v2), std::begin(result), op)

#define apply_inner_rval2(op)                                                                      \
    std::transform(std::begin(result), std::end(result), std::begin(v1), std::begin(result), op)

#define apply_inner_lval(op)                                                                       \
    std::transform(std::begin(v1), std::end(v1), std::begin(v2), std::begin(result), op)

template <class T1, class T2>
auto apply_op(T1&& v1, T2&& v2, const Operations oper) ->
    typename std::remove_cv<typename std::remove_reference<T1>::type>::type
{
    if (size_(v1) != size_(v2))
    {
        throw std::runtime_error("apply_op::nonequal input sizes");
    }

    if (std::is_rvalue_reference<decltype(std::forward<T1>(v1))>::value)
    {
        // if first argument is rvalue-reference, we use move it to result to improve performance
        auto result = std::forward<T1>(v1);
        applier(apply_inner_rval1, oper, typename decltype(result)::value_type);
        return result;
    }
    else if (std::is_rvalue_reference<decltype(std::forward<T2>(v2))>::value &&
             (oper == Operations::PLUS || oper == Operations::MULT))
    {
        // if second argument is rvalue-reference for commutativus operations,
        // we use move it to result to improve performance
        auto result = std::forward<T2>(v2);
        applier(apply_inner_rval2, oper, typename decltype(result)::value_type);
        return result;
    }

    // in other cases, the best way is preallocate memory and the write
    // result

    typename std::remove_cv<typename std::remove_reference<T1>::type>::type result(size_(v1));

    applier(apply_inner_lval, oper, typename decltype(result)::value_type);

    return result;
}
}

#endif // COMMON_TOOLS_JSON_ARRAYS_H