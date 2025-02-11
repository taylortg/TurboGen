#include "../include/tgParser.h"

#include <fstream>
#include <iostream>
#include <map>
#include <string>

// Trim function to remove leading and trailing whitespaces
std::string trim(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r\f\v");
    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    return (start == std::string::npos || end == std::string::npos) ? "" : str.substr(start, end - start + 1);
}

namespace tgparser {
std::map<std::string, std::string> readInputFile(const std::string& f) {
    std::string line;
    std::map<std::string, std::string> fileContent;
    std::ifstream file(f);

    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        exit(1);
    }

    while (getline(file, line)) {
        if (line.find(':') != std::string::npos) {
            size_t pos = line.find(':');
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 2);  // +2 for 1 whitespace

            // Trim the value to remove leading and trailing whitespaces
            value = trim(value);

            fileContent[key] = value;
        }
    }

    std::cout << f << " has been successfully parsed.\n";
    return fileContent;
}
}  // namespace tgparser
