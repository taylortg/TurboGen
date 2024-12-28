#ifndef TGPARSER_H
#define TGPARSER_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>


namespace tgparser {
std::map<std::string, std::string> readInputFile(const std::string& f);
}

#endif // TGPARSER_H
