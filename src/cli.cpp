#include "../include/cli.h"

#include <iostream>

CLITool::CLITool() : fileName(""), preliminarySizingFlag(false) {}

void CLITool::run() {
    std::string input;
    std::cout << "TurboGen CLI Tool\n";

    while (true) {
        std::cout << "> Do you want to run the correlation based preliminary sizing tool (y/n)? ";
        std::getline(std::cin, input);
        if (input == "y" || input == "Y") {
            preliminarySizingFlag = true;
            break;
        } else if (input == "n" || input == "N") {
            preliminarySizingFlag = false;
            break;
        } else {
            std::cout << "Invalid input. Please enter 'y' or 'n'.\n";
        }
    }

    while (true) {
        std::cout << "> Do you want to provide an input file (y/n)? ";
        std::getline(std::cin, input);
        if (input == "y" || input == "Y") {
            std::cout << "> File name: ";
            std::getline(std::cin, input);
            fileName = input;
            break;
        } else if (input == "n" || input == "N") {
            break;
        } else {
            std::cout << "Invalid input. Please enter 'y' or 'n'.\n";
        }
    }
}

bool CLITool::fileIsEmpty() { return fileName.empty(); }

const std::string CLITool::getFileName() { return fileName; }

const bool CLITool::getSizeFlag() { return preliminarySizingFlag; }