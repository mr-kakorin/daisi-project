#ifndef SERREFLECTION_JSON_READER_H
#define SERREFLECTION_JSON_READER_H

#include <serreflection/private/base_reader.hpp>

namespace srfl
{
    template <typename T, typename Tpt>
    std::unique_ptr<T> read_json_p(Tpt &&pt);

    template <typename T>
    class JSONReader : public BaseReader<T, boost::property_tree::ptree, JSONReader<T>>
    {
      public:
        template <typename Tpt>
        JSONReader(Tpt &&pt)
            : BaseReader<T, boost::property_tree::ptree, JSONReader<T>>(std::forward<Tpt>(pt))
        {
        }

        template <typename Targ, typename TargDef>
        bool read_optional(Targ &element, TargDef &default_set, const std::string &name,
                           const bool is_fundamental = true) const
        {
            const boost::property_tree::ptree &pt = this->m_data;

            try
            {
                if (is_fundamental)
                {
                    element = pt.get<Targ>(name);
                }
                else
                {
                    pt.get_child(name);
                }
                return true;
            }
            catch (const std::exception &ex)
            {
                auto msg = "Failed to read variable: " + srfl::StructChecker<T>::class_name + "::" +
                           name + " : " + ex.what();

                auto default_val = srfl::StructChecker<T>::get_default_value(name);

                if (std::isnan(default_val))
                {
                    throw std::runtime_error(msg);
                }

                LOG(sev_lvl::warning) << msg << ". Use default value: " << default_val;

                default_set = static_cast<TargDef>(static_cast<int>(default_val));
                return false;
            }
        }

        template <typename Tvector>
        void read_fundamental_vector(std::vector<Tvector> &element, const std::string &name) const
        {
            const boost::property_tree::ptree &pt = this->m_data;

            Tvector default_val;

            if (!read_optional(default_val, default_val, name, false))
            {
                element.resize(1);
                element[0] = default_val;
                return;
            }

            for (auto &item : pt.get_child(name))
            {
                boost::property_tree::ptree tmp = item.second;
                element.push_back(tmp.get<Tvector>(""));
            }
        }

        template <typename Tvector>
        void read_objects_vector(std::vector<Tvector> &element, const std::string &name) const
        {
            const boost::property_tree::ptree &pt = this->m_data;

            for (auto &item : pt.get_child(name))
            {
                boost::property_tree::ptree tmp = item.second;
                element.push_back(*read_json_p<Tvector>(std::move(tmp)));
            }
        }

        void read_string(std::string &element, const std::string &name) const
        {
            const boost::property_tree::ptree &pt = this->m_data;

            element = pt.get<std::string>(name);
        }

        template <typename Targ>
        void read_fundamental(Targ &element, const std::string &name) const
        {
            read_optional(element, element, name);
        }

        template <typename Targ>
        void read_enum(Targ &element, const std::string &name) const
        {
            std::string value;

            if (!read_optional(value, element, name))
            {
                return;
            }

            auto default_val = srfl::StructChecker<T>::get_default_value(name);

            auto it = srfl::EnumParser<Targ>::parse_map.find(value);

            if (srfl::EnumParser<Targ>::parse_map.end() == it)
            {
                auto msg = "Failed to read variable: " + srfl::StructChecker<T>::class_name + "::" +
                           name + " : unexpected value \"" + value + "\"";

                if (std::isnan(default_val))
                {
                    throw std::runtime_error(msg);
                }

                LOG(sev_lvl::warning) << msg << ". Use default value: " << default_val;

                element = static_cast<Targ>(static_cast<int>(default_val));
            }
            else
            {
                element = it->second;
            }
        }

        template <typename Targ>
        void read_class(Targ &element, const std::string &name) const
        {
            const boost::property_tree::ptree &pt = this->m_data;
            element = *read_json_p<Targ>(std::move(pt.get_child(name)));
        }
    };
}

#endif
