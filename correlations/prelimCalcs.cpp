#include "prelimCalcs.h"

#include <filesystem>
#include <fstream>
#include <sstream>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

/*
 *  General functions
 */
DesignCorrelations::DesignCorrelations(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op)
    : thermo(thermo), geom(geom), op(op) {
    U2 = op.omega * (geom.r2 * MM_M);
    phi = op.mfr / (thermo.props.D * U2 * std::pow((2 * geom.r2 * MM_M), 2));
    fmt::print("phi: {:.4f}\n", phi);
}

/*
 *  Casey and Robinson specific
 */
CaseyRobinsonCorrelations::CaseyRobinsonCorrelations(const ThermoProps& thermo, const Geometry& geom,
                                                     const OperatingCondition& op)
    : DesignCorrelations(thermo, geom, op) {
    Mu2 = U2 / thermo.props.A;
    fmt::print("Mu2: {:.4f}\n", Mu2);
    if (Mu2 >= 0.8) {
        P = phi * (Mu2 - 0.8);
        deta_p = -k4 * P - k5 * P * P;
    } else {
        deta_p = 0.0;
    }

    if (phi >= 0.08) {
        eta_p = eta_max * (1 - k3 * std::pow((phi - phi_max), 2)) + deta_p;
    } else {
        eta_p = eta_max * (1 - k1 * std::pow((phi_max - phi), 2) - k2 * std::pow((phi_max - phi), 4)) + deta_p;
    }
    fmt::print("eta_p: {:.4f}\n\n\n", eta_p * 100);
    double D2 = std::sqrt(op.mfr / (thermo.props.D * U2 * phi));
    fmt::print("eta_p: {:.4f}\n\n\n", eta_p * 100);

    namespace fs = std::filesystem;
    fs::path dirPath = "./../tmp";

    try {
        if (!fs::exists(dirPath)) {
            fs::create_directory(dirPath);
        }

        std::ofstream file(dirPath.string() + "/correlation.csv", std::ios::out);
        if (file.is_open()) {
            file << fmt::format("Category,Data\n");
            file << fmt::format("{},{}\n", "Mu2", Mu2);
            file << fmt::format("{},{}\n", "gamma", thermo.props.Y);
            file << fmt::format("{},{}\n", "eta_p", eta_p);
        }
        file.close();

        auto result = std::system("python ./../scripts/correlation.py");
        if (result == 0) {
            fmt::print(fg(fmt::color::green), "Python script for plotting correlation data ran successfully\n");
        } else {
            fmt::print(fg(fmt::color::red), "Error running correlation script\n");
        }

        // remove tmp folder that was created after the script executes
        try {
            fs::remove_all(dirPath);
        } catch (const fs::filesystem_error& err) {
            fmt::print("Error: {}\n", err.what());
        }
    } catch (std::exception& e) {
        fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
    }
}