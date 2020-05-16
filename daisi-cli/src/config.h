//
// Created by artoria on 3/15/20.
//

#ifndef DAISI_CLI_CONFIG_H
#define DAISI_CLI_CONFIG_H

#include <vector>
#include <string>
#include "json/json.hpp"
#include "TasksCreator.h"

using json = nlohmann::json;


enum class TypeId {
	Diode, Triode
};

ITasksCreator *read_config( char const *config_path );

#endif //DAISI_CLI_CONFIG_H
