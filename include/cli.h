#ifndef CLI_H
#define CLI_H

#include <string>

#include "common.h"

enum Prompts {
    // RUN_CORRELATION,
    // INPUT_FILE,
    INLET_OPTIMIZATION
};

class CLITool {
   private:
    std::string fileName;

   public:
    CLITool(Flags* flags);
    void run();

    const std::string getFileName();
    bool fileIsEmpty();

    Flags* flags;
};

#endif  // CLI_H