//
// Created by artoria on 3/15/20.
//

#include "config.h"


std::map<std::string, TypeId> g_types_map = {{"diode",  TypeId::Diode},
                                             {"triode", TypeId::Triode}};


ITasksCreator *read_diode( const json &json_object ) {
	auto emitter_v = json_object["potentials"]["emitter"].get<std::vector<double>>();
	auto anode_v = json_object["potentials"]["anode"].get<std::vector<double>>();
	auto file_name = json_object["project_path"].get<std::string>();
	auto timeout_ns = json_object["timeout"].get<double>();
	auto num_cores = json_object["num_cores"].get<unsigned int>();
	auto energy_type = json_object["energy_distribution_type"].get<std::string>();
	return new DiodeTasksCreator( std::move( anode_v ), std::move( emitter_v ), std::move( file_name ), std::move( energy_type ), timeout_ns, num_cores );
}

ITasksCreator *read_triode( const json &json_object ) {
	auto emitter_v = json_object["potentials"]["emitter"].get<std::vector<double>>();
	auto anode_v = json_object["potentials"]["anode"].get<std::vector<double>>();
	auto gateway_v = json_object["potentials"]["gateway"].get<std::vector<double>>();
	auto file_name = json_object["project_path"].get<std::string>();
	auto timeout_ns = json_object["timeout"].get<double>();
	auto num_cores = json_object["num_cores"].get<unsigned int>();
	auto energy_type = json_object["energy_distribution_type"].get<std::string>();
	return new TriodeTasksCreator( std::move( anode_v ), std::move( emitter_v ), std::move( gateway_v ),
	                               std::move( file_name ), std::move( energy_type ), timeout_ns, num_cores );
}


ITasksCreator *read_config( char const *config_path ) {
	std::ifstream input_file_stream;
	input_file_stream.open( config_path );
	std::string json_string((std::istreambuf_iterator<char>( input_file_stream )),
	                        (std::istreambuf_iterator<char>()));
	auto json_object = json::parse( json_string );

	if (json_object.find( "type" ) == json_object.end())
		return nullptr;
	auto type = json_object["type"].get<std::string>();
	if (g_types_map.find( type ) == g_types_map.end())
		return nullptr;
	TypeId type_id = g_types_map[type];

	switch (type_id) {
		case TypeId::Diode: {
			return read_diode( json_object );
		}
		case TypeId::Triode: {
			return read_triode( json_object );
		}
		default:
			return nullptr;
	}
}