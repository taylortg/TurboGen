#include "config.h"

#include "..\..\..\externals\fmt\include\fmt\core.h"

std::map<std::string, std::string> ReadUserInput(std::istream& file) {
    std::map<std::string, std::string> userInput;
    std::string line;

    while (std::getline(file, line)) {
        if (line[0] != '#') {
            if (line.find('=') != std::string::npos) {
                int index = line.find(' = ');
                if (index != std::string::npos && index > 0 && index < line.size() - 1) {
                    std::string key = line.substr(0, line.find('='));
                    std::string value = line.substr(line.find('=') + 1);
                    // get rid of whitespace around the string
                    key.erase(std::remove_if(key.begin(), key.end(), std::isspace), key.end());
                    value.erase(std::remove_if(value.begin(), value.end(), std::isspace), value.end());

                    userInput[key] = value;
                }
            }
        }
    }

    return userInput;
}