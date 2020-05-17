#ifndef COMM_TOOLS_PROPAGATE_CONST_HPP_
#define COMM_TOOLS_PROPAGATE_CONST_HPP_


#include <memory>
#include <type_traits>

namespace std
{

namespace experimental
{
inline namespace fundamentals_v2
{
/**
 * @defgroup propagate_const Const-propagating wrapper
 * @ingroup experimental
 *
 * A const-propagating wrapper that propagates const to pointer-like members,
 * as described in n4388 "A Proposal to Add a Const-Propagating Wrapper
 * to the Standard Library".
 *
 * @{
 */

/// Const-propagating wrapper.
template <typename _Tp>
class propagate_const
{
  public:
    typedef remove_reference_t<decltype(*std::declval<_Tp&>())> element_type;

  private:
    template <typename _Up>
    struct __is_propagate_const : false_type
    {
    };

    template <typename _Up>
    struct __is_propagate_const<propagate_const<_Up>> : true_type
    {
    };

    template <typename _Up>
    friend const _Up& get_underlying(const propagate_const<_Up>& __pt) noexcept;
    template <typename _Up>
    friend _Up& get_underlying(propagate_const<_Up>& __pt) noexcept;

    template <typename _Up>
    static element_type* __to_raw_pointer(_Up* __u)
    {
        return __u;
    }

    template <typename _Up>
    static element_type* __to_raw_pointer(_Up& __u)
    {
        return __u.get();
    }

    template <typename _Up>
    static const element_type* __to_raw_pointer(const _Up* __u)
    {
        return __u;
    }

    template <typename _Up>
    static const element_type* __to_raw_pointer(const _Up& __u)
    {
        return __u.get();
    }

  public:
    static_assert(__and_<is_object<typename remove_pointer<_Tp>::type>, __not_<is_array<_Tp>>,
                         __or_<is_class<_Tp>, is_pointer<_Tp>>>::value,
                  "propagate_const requires a class or a pointer to an"
                  " object type");

    // [propagate_const.ctor], constructors
    propagate_const()                           = default;
    propagate_const(const propagate_const& __p) = delete;
    propagate_const(propagate_const&& __p)      = default;
    template <
        typename _Up,
        typename enable_if<__and_<is_constructible<_Tp, _Up&&>, is_convertible<_Up&&, _Tp>>::value,
                           bool>::type = true>
    propagate_const(propagate_const<_Up>&& __pu) : _M_t(std::move(get_underlying(__pu)))
    {
    }
    template <typename _Up,
              typename enable_if<
                  __and_<is_constructible<_Tp, _Up&&>, __not_<is_convertible<_Up&&, _Tp>>>::value,
                  bool>::type = false>
    explicit propagate_const(propagate_const<_Up>&& __pu) : _M_t(std::move(get_underlying(__pu)))
    {
    }
    template <
        typename _Up,
        typename enable_if<__and_<is_constructible<_Tp, _Up&&>, is_convertible<_Up&&, _Tp>,
                                  __not_<__is_propagate_const<typename decay<_Up>::type>>>::value,
                           bool>::type = true>
    propagate_const(_Up&& __u) : _M_t(std::forward<_Up>(__u))
    {
    }
    template <
        typename _Up,
        typename enable_if<__and_<is_constructible<_Tp, _Up&&>, __not_<is_convertible<_Up&&, _Tp>>,
                                  __not_<__is_propagate_const<typename decay<_Up>::type>>>::value,
                           bool>::type = false>
    explicit propagate_const(_Up&& __u) : _M_t(std::forward<_Up>(__u))
    {
    }

    // [propagate_const.assignment], assignment
    propagate_const& operator=(const propagate_const& __p) = delete;
    propagate_const& operator=(propagate_const&& __p) = default;

    template <typename _Up, typename = typename enable_if<is_convertible<_Up&&, _Tp>::value>::type>
    propagate_const& operator=(propagate_const<_Up>&& __pu)
    {
        _M_t = std::move(get_underlying(__pu));
        return *this;
    }

    template <typename _Up,
              typename = typename enable_if<
                  __and_<is_convertible<_Up&&, _Tp>,
                         __not_<__is_propagate_const<typename decay<_Up>::type>>>::value>::type>
    propagate_const& operator=(_Up&& __u)
    {
        _M_t = std::forward<_Up>(__u);
        return *this;
    }

    // [propagate_const.const_observers], const observers
    explicit operator bool() const
    {
        return bool(_M_t);
    }

    const element_type* operator->() const
    {
        return get();
    }

    template <
        typename _Up = _Tp,
        typename enable_if<__or_<is_pointer<_Up>, is_convertible<_Up, const element_type*>>::value,
                           bool>::type = true>
    operator const element_type*() const
    {
        return get();
    }

    const element_type& operator*() const
    {
        return *get();
    }

    const element_type* get() const
    {
        return __to_raw_pointer(_M_t);
    }

    // [propagate_const.non_const_observers], non-const observers
    element_type* operator->()
    {
        return get();
    }

    template <
        typename _Up = _Tp,
        typename enable_if<__or_<is_pointer<_Up>, is_convertible<_Up, const element_type*>>::value,
                           bool>::type = true>
    operator element_type*()
    {
        return get();
    }

    element_type& operator*()
    {
        return *get();
    }

    element_type* get()
    {
        return __to_raw_pointer(_M_t);
    }

  private:
    _Tp _M_t;
};

// [propagate_const.relational], relational operators
template <typename _Tp>
bool operator==(const propagate_const<_Tp>& __pt, nullptr_t)
{
    return get_underlying(__pt) == nullptr;
}

template <typename _Tp>
bool operator==(nullptr_t, const propagate_const<_Tp>& __pu)
{
    return nullptr == get_underlying(__pu);
}

template <typename _Tp>
bool operator!=(const propagate_const<_Tp>& __pt, nullptr_t)
{
    return get_underlying(__pt) != nullptr;
}

template <typename _Tp>
bool operator!=(nullptr_t, const propagate_const<_Tp>& __pu)
{
    return nullptr != get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator==(const propagate_const<_Tp>& __pt, const propagate_const<_Up>& __pu)
{
    return get_underlying(__pt) == get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator!=(const propagate_const<_Tp>& __pt, const propagate_const<_Up>& __pu)
{
    return get_underlying(__pt) != get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator<(const propagate_const<_Tp>& __pt, const propagate_const<_Up>& __pu)
{
    return get_underlying(__pt) < get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator>(const propagate_const<_Tp>& __pt, const propagate_const<_Up>& __pu)
{
    return get_underlying(__pt) > get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator<=(const propagate_const<_Tp>& __pt, const propagate_const<_Up>& __pu)
{
    return get_underlying(__pt) <= get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator>=(const propagate_const<_Tp>& __pt, const propagate_const<_Up>& __pu)
{
    return get_underlying(__pt) >= get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator==(const propagate_const<_Tp>& __pt, const _Up& __u)
{
    return get_underlying(__pt) == __u;
}

template <typename _Tp, typename _Up>
bool operator!=(const propagate_const<_Tp>& __pt, const _Up& __u)
{
    return get_underlying(__pt) != __u;
}

template <typename _Tp, typename _Up>
bool operator<(const propagate_const<_Tp>& __pt, const _Up& __u)
{
    return get_underlying(__pt) < __u;
}

template <typename _Tp, typename _Up>
bool operator>(const propagate_const<_Tp>& __pt, const _Up& __u)
{
    return get_underlying(__pt) > __u;
}

template <typename _Tp, typename _Up>
bool operator<=(const propagate_const<_Tp>& __pt, const _Up& __u)
{
    return get_underlying(__pt) <= __u;
}

template <typename _Tp, typename _Up>
bool operator>=(const propagate_const<_Tp>& __pt, const _Up& __u)
{
    return get_underlying(__pt) >= __u;
}

template <typename _Tp, typename _Up>
bool operator==(const _Tp& __t, const propagate_const<_Up>& __pu)
{
    return __t == get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator!=(const _Tp& __t, const propagate_const<_Up>& __pu)
{
    return __t != get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator<(const _Tp& __t, const propagate_const<_Up>& __pu)
{
    return __t < get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator>(const _Tp& __t, const propagate_const<_Up>& __pu)
{
    return __t > get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator<=(const _Tp& __t, const propagate_const<_Up>& __pu)
{
    return __t <= get_underlying(__pu);
}

template <typename _Tp, typename _Up>
bool operator>=(const _Tp& __t, const propagate_const<_Up>& __pu)
{
    return __t >= get_underlying(__pu);
}

// [propagate_const.underlying], underlying pointer access
template <typename _Tp>
const _Tp& get_underlying(const propagate_const<_Tp>& __pt) noexcept
{
    return __pt._M_t;
}

template <typename _Tp>
_Tp& get_underlying(propagate_const<_Tp>& __pt) noexcept
{
    return __pt._M_t;
}

// @} group propagate_const
} // namespace fundamentals_v2
} // namespace experimental

// [propagate_const.hash], hash support
template <typename _Tp>
struct hash<experimental::propagate_const<_Tp>>
{
    using result_type   = size_t;
    using argument_type = experimental::propagate_const<_Tp>;

    size_t operator()(const experimental::propagate_const<_Tp>& __t) const
        noexcept(noexcept(hash<_Tp>{}(get_underlying(__t))))
    {
        return hash<_Tp>{}(get_underlying(__t));
    }
};

// [propagate_const.comparison_function_objects], comparison function objects
template <typename _Tp>
struct equal_to<experimental::propagate_const<_Tp>>
{
    bool operator()(const experimental::propagate_const<_Tp>& __x,
                    const experimental::propagate_const<_Tp>& __y) const
    {
        return equal_to<_Tp>{}(get_underlying(__x), get_underlying(__y));
    }

    typedef experimental::propagate_const<_Tp> first_argument_type;
    typedef experimental::propagate_const<_Tp> second_argument_type;
    typedef bool                               result_type;
};

template <typename _Tp>
struct not_equal_to<experimental::propagate_const<_Tp>>
{
    bool operator()(const experimental::propagate_const<_Tp>& __x,
                    const experimental::propagate_const<_Tp>& __y) const
    {
        return not_equal_to<_Tp>{}(get_underlying(__x), get_underlying(__y));
    }

    typedef experimental::propagate_const<_Tp> first_argument_type;
    typedef experimental::propagate_const<_Tp> second_argument_type;
    typedef bool                               result_type;
};

template <typename _Tp>
struct less<experimental::propagate_const<_Tp>>
{
    bool operator()(const experimental::propagate_const<_Tp>& __x,
                    const experimental::propagate_const<_Tp>& __y) const
    {
        return less<_Tp>{}(get_underlying(__x), get_underlying(__y));
    }

    typedef experimental::propagate_const<_Tp> first_argument_type;
    typedef experimental::propagate_const<_Tp> second_argument_type;
    typedef bool                               result_type;
};

template <typename _Tp>
struct greater<experimental::propagate_const<_Tp>>
{
    bool operator()(const experimental::propagate_const<_Tp>& __x,
                    const experimental::propagate_const<_Tp>& __y) const
    {
        return greater<_Tp>{}(get_underlying(__x), get_underlying(__y));
    }

    typedef experimental::propagate_const<_Tp> first_argument_type;
    typedef experimental::propagate_const<_Tp> second_argument_type;
    typedef bool                               result_type;
};

template <typename _Tp>
struct less_equal<experimental::propagate_const<_Tp>>
{
    bool operator()(const experimental::propagate_const<_Tp>& __x,
                    const experimental::propagate_const<_Tp>& __y) const
    {
        return less_equal<_Tp>{}(get_underlying(__x), get_underlying(__y));
    }

    typedef experimental::propagate_const<_Tp> first_argument_type;
    typedef experimental::propagate_const<_Tp> second_argument_type;
    typedef bool                               result_type;
};

template <typename _Tp>
struct greater_equal<experimental::propagate_const<_Tp>>
{
    bool operator()(const experimental::propagate_const<_Tp>& __x,
                    const experimental::propagate_const<_Tp>& __y) const
    {
        return greater_equal<_Tp>{}(get_underlying(__x), get_underlying(__y));
    }

    typedef experimental::propagate_const<_Tp> first_argument_type;
    typedef experimental::propagate_const<_Tp> second_argument_type;
    typedef bool                               result_type;
};

} // namespace std
#endif
