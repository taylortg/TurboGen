#include "makeBlade.h"

#include <filesystem>
#include <fstream>

#include "..\..\externals\fmt\include\fmt\core.h"
#include "..\common\config.h"

namespace fs = std::filesystem;
int main(int argc, char* argv[]) {
    std::ifstream file;
    try {
        fs::path currentPath = fs::current_path();
        file.open(currentPath.string() + "\\" + argv[1], std::ios::in);
    } catch (fs::filesystem_error& err) {
        fmt::print("Error: {}\n", err.what());
        fmt::print("Add a config file after the executable file, i.e. makeBlade input.cfg\n");
        exit(EXIT_FAILURE);
    }

    if (file.is_open()) {
        fmt::print("Opened file {}\n", argv[1]);
    }
    userInput = ReadUserInput(file);
    file.close();

    for (auto in : userInput) {
        fmt::println("Key: {}\nValue: {}", in.first, in.second);
    }
}