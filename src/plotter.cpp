#include "../include/plotter.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

void plotVelocityTriangle(const Impeller& impeller) {
    namespace fs = std::filesystem;
    fs::path dirPath = "./../tmp";

    try {
        if (!fs::exists(dirPath)) {
            fs::create_directory(dirPath);
        }
        std::ofstream file(dirPath.string() + "/velocities.csv", std::ios::out);
        if (file.is_open()) {
            file << fmt::format("Category,X,Y\n");
            file << fmt::format("{},{},{}\n", "Shroud", 0, 0);
            file << fmt::format("{},{},{}\n", "Shroud", -impeller.inlet.tip.C_theta, impeller.inlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Shroud", -impeller.inlet.tip.C_theta, impeller.inlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Shroud", impeller.inlet.tip.U - impeller.inlet.tip.C_theta,
                                impeller.inlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Shroud", 0, 0);
            file << fmt::format("{},{},{}\n", "Shroud", impeller.inlet.tip.U - impeller.inlet.tip.C_theta,
                                impeller.inlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Hub", 0, 0);
            file << fmt::format("{},{},{}\n", "Hub", -impeller.inlet.hub.C_theta, impeller.inlet.hub.C_m);
            file << fmt::format("{},{},{}\n", "Hub", -impeller.inlet.hub.C_theta, impeller.inlet.hub.C_m);
            file << fmt::format("{},{},{}\n", "Hub", impeller.inlet.hub.U - impeller.inlet.hub.C_theta,
                                impeller.inlet.hub.C_m);
            file << fmt::format("{},{},{}\n", "Hub", 0, 0);
            file << fmt::format("{},{},{}\n", "Hub", impeller.inlet.hub.U - impeller.inlet.hub.C_theta,
                                impeller.inlet.hub.C_m);
        }
        file.close();
    } catch (std::exception& e) {
        fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
    }

    auto result = std::system("python ./../scripts/velTrianglePlot.py");
    if (result == 0) {
        fmt::print(fg(fmt::color::green), "Python script for velocity triangles ran successfully\n");
    } else {
        fmt::print(fg(fmt::color::red), "Error running velocity triangle script\n");
    }

    try {
        fs::remove_all(dirPath);
    } catch (const fs::filesystem_error& err) {
        fmt::print("Error: {}\n", err.what());
    }
}