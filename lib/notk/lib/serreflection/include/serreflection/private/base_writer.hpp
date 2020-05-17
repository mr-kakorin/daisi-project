#ifndef SERREFLECTION_BASE_WRITER_H
#define SERREFLECTION_BASE_WRITER_H

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
    class BaseWriter
    {
      public:
        BaseWriter(const T &object) : m_value(object), m_data(make_unique<Tdata>())
        {
        }

        template <typename Targ>
        void operator()(const Targ &x) const
        {
            constexpr size_t index = Targ::value;

            write(boost::fusion::at_c<index>(m_value),
                  boost::fusion::extension::struct_member_name<T, index>::call());
        }

      protected:
        void write(const bool &element, const std::string &name) const
        {
            static_cast<const Base *>(this)->write_fundamental(element, name);
        }

        template <typename Targ>
        void write(const Targ &element, const std::string &name,
                   typename std::enable_if<std::is_fundamental<Targ>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->write_fundamental(element, name);
        }

        template <typename Targ>
        void write(const Targ &element, const std::string &name,
                   typename std::enable_if<std::is_enum<Targ>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->write_enum(element, name);
        }

        template <typename Targ>
        void write(const Targ &element, const std::string &name,
                   typename std::enable_if<std::is_class<Targ>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->write_class(element, name);
        }

        template <typename Tvector>
        void write(const std::vector<Tvector> &element, const std::string &name,
                   typename std::enable_if<std::is_fundamental<Tvector>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->write_fundamental_vector(element, name);
        }

        template <typename Tvector>
        void write(const std::vector<Tvector> &element, const std::string &name,
                   typename std::enable_if<std::is_class<Tvector>::value>::type * = 0) const
        {
            static_cast<const Base *>(this)->write_objects_vector(element, name);
        }

        void write(const std::string &element, const std::string &name) const
        {
            static_cast<const Base *>(this)->write_string(element, name);
        }

        const T &              m_value;  // serializable object
        std::unique_ptr<Tdata> m_data;   // source of deserialization data

      public:
        std::unique_ptr<Tdata> &get()
        {
            return m_data;
        }

      private:
    };
}

#endif
