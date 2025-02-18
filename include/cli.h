#ifndef CLI_H
#define CLI_H

#include <array>
#include <string>

#include "common.h"

enum Prompts { RUN_CORRELATION, INPUT_FILE, INLET_OPTIMIZATION };

class CLITool {
   private:
    std::string fileName;

   public:
    CLITool(Flags* flags);
    void run();
    bool getOptionsFromUser(std::string message, bool& flag);

    const std::string getFileName();
    bool fileIsEmpty();

    Flags* flags;
    std::array<Prompts, 2> prompts = {RUN_CORRELATION, INLET_OPTIMIZATION};
};

#endif  // CLI_H