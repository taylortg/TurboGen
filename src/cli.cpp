#include "../include/cli.h"

#include <iostream>

CLITool::CLITool(Flags* flags) : fileName(""), flags(flags) {}

void CLITool::run() {
    std::string input;

    for (auto prompt : prompts) {
        switch (prompt) {
        case RUN_CORRELATION:
            getOptionsFromUser("> Do you want to run the correlation based preliminary sizing tool (y/n)? ",
                               flags->preliminarySizingFlag);
            break;
        case INLET_OPTIMIZATION:
            getOptionsFromUser("> Do you want to run the inlet inducer optimization loop (y/n)? ",
                               flags->inletInducerOptFlag);
            break;
        default:
            std::cout << "How did we end up at the default case?\n";
            break;
        }
    }

    // while (true) {
    //     std::cout << "> Do you want to provide an input file (y/n)? ";
    //     std::getline(std::cin, input);
    //     if (input == "y" || input == "Y") {
    //         std::cout << "> File name: ";
    //         std::getline(std::cin, input);
    //         fileName = input;
    //         break;
    //     } else if (input == "n" || input == "N") {
    //         break;
    //     } else {
    //         std::cout << "Invalid input. Please enter 'y' or 'n'.\n";
    //     }
    // }
}

bool CLITool::getOptionsFromUser(std::string message, bool& flag) {
    std::string input;

    while (true) {
        std::cout << message;
        std::getline(std::cin, input);
        if (input == "y" || input == "Y") {
            flag = true;
            break;
        } else if (input == "n" || input == "N") {
            flag = false;
            break;
        } else {
            std::cout << "Invalid input. Please enter 'y' or 'n'.\n";
        }
    }
    return flag;
}

bool CLITool::fileIsEmpty() { return fileName.empty(); }

const std::string CLITool::getFileName() { return fileName; }
