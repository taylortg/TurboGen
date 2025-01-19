#include <thread>

// #include "../externals/CoolProp/include/CoolProp.h"
#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"
#include "../include/common.h"
#include "../include/impeller.h"
// #include "../include/plotter.h"
#include "../include/tgParser.h"
#include "../include/thermo.h"

int main(int argc, char** argv) {
#ifdef DEBUG
    std::cout << "DEBUG mode is active.\n";
#endif

    //const std::map<std::string, std::string> inputData = tgparser::readInputFile("../input.in");
    const std::map<std::string, std::string> inputData = tgparser::readInputFile("C:/Users/mastodon/Documents/TurboGen/input.in");

    OperatingCondition op{};
    Geometry geom{};
    ThermoProps thermo;

    std::thread opThread([&]() {
        op = filterFileContent_op(inputData);
#ifdef DEBUG
        std::cout << "Filtering op struct.\n";
#endif
    });
    std::thread goemThread([&]() {
        geom = filterFileContent_geom(inputData);
#ifdef DEBUG
        std::cout << "Filtering geom struct.\n";
#endif
    });
    std::thread thermoThread([&]() {
        thermo = filterFileContent_thermo(inputData);
#ifdef DEBUG
        std::cout << "Filtering ThermoProps class.\n";
#endif
    });

    opThread.join();
    goemThread.join();
    thermoThread.join();

    fmt::print(fg(fmt::color::green), "All filtering completed.\n");

    Impeller impeller(thermo, geom, op);
    impeller.calculateInletCondition("Aungier");
    // impeller.calculateOutletCondition("Japikse", "Wiesner");
    // impeller.estimateAxialLength();
    // plotVelocityTriangle(impeller, true);

    return 0;
}