#ifndef CLI_H
#define CLI_H

#include <string>

class CLITool {
   private:
    std::string fileName;
    bool preliminarySizingFlag;

   public:
    CLITool();
    void run();

    const std::string getFileName();
    bool fileIsEmpty();

    const bool getSizeFlag();
};

#endif  // CLI_H