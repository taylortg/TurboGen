#ifndef CONFIG_H
#define CONFIG_H

#include <fstream>
#include <map>
#include <string>

std::map<std::string, std::string> ReadUserInput(std::istream& file);

#endif  // CONFIG_H