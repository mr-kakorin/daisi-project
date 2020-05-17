#ifndef COMM_TOOLS_CHILD_FACTORY_HPP
#define COMM_TOOLS_CHILD_FACTORY_HPP

#include <map>
#include <memory>

// How to use this factory:
// struct Base
// {
//     BUILD_CHILD_FACTORY(std::string, Base)
// };

// struct Der : Base
// {
//     BUILD_CHILD(Der, Base, "Der")
// };

// REGISTER_CHILD_NON_TEMPLATE(Der, Base)

// auto der = Base::make_child("Der");

#define COMMA ,

#define BUILD_CHILD_FACTORY(TKey, TParent)                                                         \
    struct ChildFactory                                                                            \
    {                                                                                              \
        static ChildFactory& instance()                                                            \
        {                                                                                          \
            static ChildFactory instance;                                                          \
            return instance;                                                                       \
        }                                                                                          \
        bool add(const TKey& key, std::unique_ptr<TParent>&& child_ptr)                            \
        {                                                                                          \
            map[key] = std::move(child_ptr);                                                       \
            return true;                                                                           \
        }                                                                                          \
        std::shared_ptr<TParent> make_child(const TKey& key)                                       \
        {                                                                                          \
            auto it = map.find(key);                                                               \
            if (it == map.end())                                                                   \
            {                                                                                      \
                return nullptr;                                                                    \
            }                                                                                      \
            return it->second->clone();                                                            \
        }                                                                                          \
        std::unique_ptr<TParent> make_unique_child(const TKey& key)                                \
        {                                                                                          \
            auto it = map.find(key);                                                               \
            if (it == map.end())                                                                   \
            {                                                                                      \
                return nullptr;                                                                    \
            }                                                                                      \
            return it->second->clone_unique();                                                     \
        }                                                                                          \
                                                                                                   \
        std::map<TKey, std::unique_ptr<TParent>> map;                                              \
    };                                                                                             \
    static std::shared_ptr<TParent> make_child(const TKey& key)                                    \
    {                                                                                              \
        return ChildFactory::instance().make_child(key);                                           \
    }                                                                                              \
    static std::unique_ptr<TParent> make_unique_child(const TKey& key)                             \
    {                                                                                              \
        return ChildFactory::instance().make_unique_child(key);                                    \
    }                                                                                              \
    virtual TKey                     get_key()      = 0;                                           \
    virtual std::shared_ptr<TParent> clone()        = 0;                                           \
    virtual std::unique_ptr<TParent> clone_unique() = 0;                                           \
    const static TKey                key;

#define BUILD_CHILD(T, TParent, key_value)                                                         \
    virtual std::shared_ptr<TParent> clone() override final                                        \
    {                                                                                              \
        return std::make_shared<T>(*this);                                                         \
    }                                                                                              \
    virtual std::unique_ptr<TParent> clone_unique() override final                                 \
    {                                                                                              \
        return std::make_unique<T>(*this);                                                         \
    }                                                                                              \
                                                                                                   \
    virtual typename std::remove_const<decltype(TParent::key)>::type get_key() override final      \
    {                                                                                              \
        return key_value;                                                                          \
    }                                                                                              \
    const static decltype(TParent::key) key;                                                       \
    const static bool                   registered;

#define REGISTER_CHILD(T, TParent)                                                                 \
    template <>                                                                                    \
    const bool T::registered = TParent::ChildFactory::instance().add(                              \
        std::make_shared<T>()->get_key(), std::make_unique<T>());                                  \
    template <>                                                                                    \
    const decltype(TParent::key) T::key = std::make_shared<T>()->get_key();

#define REGISTER_CHILD_NON_TEMPLATE(T, TParent)                                                    \
    const bool T::registered = TParent::ChildFactory::instance().add(                              \
        std::make_shared<T>()->get_key(), std::make_unique<T>());                                  \
    const decltype(TParent::key) T::key = std::make_shared<T>()->get_key();

#endif
