#ifndef SERREFLECTION_BASE_READER_H
#define SERREFLECTION_BASE_READER_H

#include <cmath>
#include <memory>

#include <boost/fusion/algorithm.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/size.hpp>
#include <boost/fusion/tuple.hpp>

#include <boost/mpl/range_c.hpp>

#include <serreflection/private/defines_private.hpp>

namespace srfl
{
    template <typename T, typename Tdata, typename Base>
    class BaseReader
    {
      public:
        template <typename TTdata>
        BaseReader(TTdata &&data) : m_value(make_unique<T>()), m_data(std::forward<TTdata>(data))
        {
        }

        template <typename Targ>
        void operator()(const Targ &x) const
        {
            constexpr size_t index = Targ::value;

            read(boost::fusion::at_c<index>(*m_value),
                 boost::fusion::extension::struct_member_name<T, index>::call());
        }

      protected:
        void read(bool &element, const std::string &name) const
        {
            static_cast<const Base *>(this)->read_fundamental(element, name);
        }

        template <typename Targ>
        void read(Targ &element, const std::string &name,
                  typename std::enable_if<std::is_fundamental<Targ>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->read_fundamental(element, name);
            srfl::StructChecker<T>::check_value(element, name);
        }

        template <typename Targ>
        void read(Targ &element, const std::string &name,
                  typename std::enable_if<std::is_enum<Targ>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->read_enum(element, name);
        }

        template <typename Targ>
        void read(Targ &element, const std::string &name,
                  typename std::enable_if<std::is_class<Targ>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->read_class(element, name);
        }

        template <typename Tvector>
        void read(std::vector<Tvector> &element, const std::string &name,
                  typename std::enable_if<std::is_fundamental<Tvector>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->read_fundamental_vector(element, name);
            srfl::StructChecker<T>::check_value(element, name);
        }

        template <typename Tvector>
        void read(std::vector<Tvector> &element, const std::string &name,
                  typename std::enable_if<std::is_class<Tvector>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->read_objects_vector(element, name);
        }

        void read(std::string &element, const std::string &name) const
        {
            static_cast<const Base *>(this)->read_string(element, name);
        }
        std::unique_ptr<T> m_value;  // pointer to serializible object

        Tdata m_data;  // source of serialization data

      public:
        std::unique_ptr<T> &get()
        {
            return m_value;
        }

      private:
    };
}

#endif
