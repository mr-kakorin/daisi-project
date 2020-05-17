#ifndef SERREFLECTION_JSON_WRITER_H
#define SERREFLECTION_JSON_WRITER_H

#include <serreflection/private/base_writer.hpp>

namespace srfl
{
    template <typename T>
    std::unique_ptr<boost::property_tree::ptree> write_json_p(const T &object);

    template <typename T>
    class JSONWriter : public BaseWriter<T, boost::property_tree::ptree, JSONWriter<T>>
    {
      public:
        JSONWriter(const T &object)
            : BaseWriter<T, boost::property_tree::ptree, JSONWriter<T>>(object)
        {
        }

        template <typename Tvector>
        void write_fundamental_vector(const std::vector<Tvector> &element,
                                      const std::string &         name) const
        {
            boost::property_tree::ptree children;
            for (const auto &val : element)
            {
                boost::property_tree::ptree tmp;
                tmp.put("", val);
                children.push_back(std::make_pair("", std::move(tmp)));
            }
            this->m_data->add_child(name, children);
        }

        template <typename Tvector>
        void write_objects_vector(const std::vector<Tvector> &element,
                                  const std::string &         name) const
        {
            boost::property_tree::ptree children;
            for (const auto &val : element)
            {
                auto tmp = write_json_p(val);
                if (!tmp)
                {
                    throw std::runtime_error("Error write object " + name);
                }
                children.push_back(std::make_pair("", std::move(*tmp)));
            }
            this->m_data->add_child(name, children);
        }

        void write_string(const std::string &element, const std::string &name) const
        {
            this->m_data->put(name, element);
        }

        template <typename Targ>
        void write_fundamental(const Targ &element, const std::string &name) const
        {
            this->m_data->put(name, element);
        }

        template <typename Targ>
        void write_enum(const Targ &element, const std::string &name) const
        {
            auto it = srfl::EnumParser<Targ>::inv_parse_map.find(element);

            if (srfl::EnumParser<Targ>::inv_parse_map.end() == it)
            {
                auto msg =
                    "Failed to write variable: " + srfl::StructChecker<T>::class_name + "::" + name;

                throw std::runtime_error(msg);
            }
            this->m_data->put(name, it->second);
        }

        template <typename Targ>
        void write_class(const Targ &element, const std::string &name) const
        {
            auto tmp = write_json_p(element);
            if (!tmp)
            {
                throw std::runtime_error("Error write object " + name);
            }
            this->m_data->add_child(name, *tmp);
        }
    };
}

#endif
