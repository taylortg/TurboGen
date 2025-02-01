#include "../include/plotter.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

std::string formatHeaderRow(const std::string& variable) {
    return fmt::format(
        "\t\t\t<tr>\n"
        "\t\t\t\t<th colspan=\"3\">{}</th>\n"
        "\t\t\t</tr>\n",
        variable);
}

std::string formatSubHeaderRow(const std::string& variable) {
    return fmt::format(
        "\t\t\t<tr>\n"
        "\t\t\t\t<td colspan=\"3\">{}</td>\n"
        "\t\t\t</tr>\n",
        variable);
}

std::string formatVariableRow() {
    return fmt::format(
        "\t\t\t<tr>\n"
        "\t\t\t\t<th>Variable</th>\n"
        "\t\t\t\t<th>Value</th>\n"
        "\t\t\t\t<th>Units</th>\n"
        "\t\t\t</tr>\n");
}

std::string formatTableRow(const std::string& variable, const std::string& units, float value, int precision) {
    std::stringstream valueStream;
    valueStream << std::fixed << std::setprecision(precision) << value;

    return fmt::format(
        "\t\t\t<tr>\n"
        "\t\t\t\t<td>{:<15}</td>\n"
        "\t\t\t\t<td>{:<15}</td>\n"
        "\t\t\t\t<td>{:<15}</td>\n"
        "\t\t\t</tr>\n",
        variable, valueStream.str(), units);
}

void plotVelocityTriangle(const Impeller& impeller, bool htmlFlag) {
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
            file << fmt::format("{},{},{}\n", "Outlet", 0, 0);
            file << fmt::format("{},{},{}\n", "Outlet", impeller.outlet.tip.W_theta, impeller.outlet.tip.W_m);
            file << fmt::format("{},{},{}\n", "Outlet", 0, 0);
            file << fmt::format("{},{},{}\n", "Outlet", 0, impeller.outlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Outlet", impeller.outlet.tip.W_theta, impeller.outlet.tip.W_m);
            file << fmt::format("{},{},{}\n", "Outlet", impeller.outlet.tip.W_theta, 0);
            file << fmt::format("{},{},{}\n", "Outlet", 0, 0);
            file << fmt::format("{},{},{}\n", "Outlet", impeller.outlet.tip.C_theta, impeller.outlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Outlet", 0, impeller.outlet.tip.C_m);
            file << fmt::format("{},{},{}\n", "Outlet", impeller.outlet.tip.C_theta, 0);
        }
        file.close();

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

        if (htmlFlag) {
            dirPath = "./../out";
            std::ofstream file(dirPath.string() + "/index.html", std::ios::out);
            if (file.is_open()) {
                file << fmt::format("<!DOCTYPE html>\n");
                file << fmt::format("<html lang=\"en\">\n");
                file << fmt::format("\t<head>\n");
                file << fmt::format("\t\t<meta charset=\"UTF-8\">\n");
                file << fmt::format("\t\t<meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">\n");
                file << fmt::format("\t\t<title>TurboGen Report</title>\n");
                file << fmt::format("\t\t<link rel=\"stylesheet\" href=\"style.css\">\n");
                file << fmt::format("\t</head>\n");
                file << fmt::format("\t<body>\n");

                file << fmt::format("\t\t<table>\n");
                file << formatHeaderRow("Geometrical Data");
                file << formatVariableRow();
                file << formatSubHeaderRow("Impeller leading edge");
                file << formatTableRow("Hub radius", "mm", impeller.geom.r1h, 2);
                file << formatTableRow("RMS radius", "mm", impeller.geom.r1rms, 2);
                file << formatTableRow("Tip radius", "mm", impeller.geom.r1t, 2);
                file << formatTableRow("Hub metal angle", "°", impeller.geom.beta1h, 2);
                file << formatTableRow("RMS metal angle", "°", impeller.geom.beta1rms, 2);
                file << formatTableRow("Tip metal angle", "°", impeller.geom.beta1t, 2);
                file << formatSubHeaderRow("Impeller exit");
                file << formatTableRow("Radius", "mm", impeller.geom.r2, 2);
                file << formatTableRow("Blade tip width", "mm", impeller.geom.b2, 2);
                file << formatTableRow("Blade metal angle", "°", impeller.geom.beta2, 2);
                file << formatSubHeaderRow("General");
                file << formatTableRow("Inlet area", "mm^2", impeller.geom.area1 * 1000000.0, 0);
                file << formatTableRow("Inlet throat area", "mm^2", impeller.geom.throatArea * 1000000.0, 0);
                file << formatTableRow("Outlet area", "mm^2", impeller.geom.area2 * 1000000.0, 0);
                file << formatTableRow("Axial length - OD ratio", "", impeller.geom.axialLengthRatio, 4);
                // file << formatTableRow("Axial length", "mm", impeller.geom.axialLength, 4);
                file << formatTableRow("Inlet blockage factor", "", impeller.geom.blockage1, 4);
                file << formatTableRow("Outlet blockage factor", "", impeller.geom.blockage2, 4);
                file << formatTableRow("Number of full blades", "", impeller.geom.ZFull, 0);
                file << formatTableRow("Number of splitter blades", "", impeller.geom.ZSplit, 0);
                file << fmt::format("\t\t</table>\n");

                file << fmt::format("\t\t<table>\n");
                file << formatHeaderRow("Thermodynamic Data");
                file << formatVariableRow();
                file << formatSubHeaderRow("Impeller leading edge");
                file << formatTableRow("Total pressure", "Pa", impeller.inlet.total.props.P, 0);
                file << formatTableRow("Total temperature", "K", impeller.inlet.total.props.T, 2);
                file << formatTableRow("Total enthalpy", "J/kg", impeller.inlet.total.props.H, 0);
                file << formatTableRow("Total entropy", "J/kg-K", impeller.inlet.total.props.S, 2);
                file << formatTableRow("Static pressure", "Pa", impeller.inlet.static_.props.P, 0);
                file << formatTableRow("Static temperature", "K", impeller.inlet.static_.props.T, 2);
                file << formatTableRow("Static enthalpy", "J/kg", impeller.inlet.static_.props.H, 0);
                file << formatTableRow("Static entropy", "J/kg-K", impeller.inlet.static_.props.S, 2);
                file << formatSubHeaderRow("Impeller exit");
                file << formatTableRow("Total pressure", "Pa", impeller.outlet.total.props.P, 0);
                file << formatTableRow("Total temperature", "K", impeller.outlet.total.props.T, 2);
                file << formatTableRow("Total enthalpy", "J/kg", impeller.outlet.total.props.H, 0);
                file << formatTableRow("Total entropy", "J/kg-K", impeller.outlet.total.props.S, 2);
                file << formatTableRow("Static pressure", "Pa", impeller.outlet.static_.props.P, 0);
                file << formatTableRow("Static temperature", "K", impeller.outlet.static_.props.T, 2);
                file << formatTableRow("Static enthalpy", "J/kg", impeller.outlet.static_.props.H, 0);
                file << formatTableRow("Static entropy", "J/kg-K", impeller.outlet.static_.props.S, 2);
                file << fmt::format("\t\t</table>\n");

                file << fmt::format("\t\t<table>\n");
                file << formatHeaderRow("Impeller Velocities");
                file << formatVariableRow();
                file << formatSubHeaderRow("Impeller leading edge");
                file << formatSubHeaderRow("Hub");
                file << formatTableRow("C", "m/s", impeller.inlet.hub.C, 3);
                file << formatTableRow("Cm", "m/s", impeller.inlet.hub.C_m, 3);
                file << formatTableRow("Cθ", "m/s", impeller.inlet.hub.C_theta, 3);
                file << formatTableRow("W", "m/s", impeller.inlet.hub.W, 3);
                file << formatTableRow("Wm", "m/s", impeller.inlet.hub.W_m, 3);
                file << formatTableRow("Wθ", "m/s", impeller.inlet.hub.W_theta, 3);
                file << formatTableRow("U", "m/s", impeller.inlet.hub.U, 3);
                file << formatTableRow("Mach number", "m/s", impeller.inlet.hub.M, 4);
                file << formatTableRow("Relative Mach number", "m/s", impeller.inlet.hub.M_rel, 4);
                file << formatSubHeaderRow("RMS");
                file << formatTableRow("C", "m/s", impeller.inlet.rms.C, 3);
                file << formatTableRow("Cm", "m/s", impeller.inlet.rms.C_m, 3);
                file << formatTableRow("Cθ", "m/s", impeller.inlet.rms.C_theta, 3);
                file << formatTableRow("W", "m/s", impeller.inlet.rms.W, 3);
                file << formatTableRow("Wm", "m/s", impeller.inlet.rms.W_m, 3);
                file << formatTableRow("Wθ", "m/s", impeller.inlet.rms.W_theta, 3);
                file << formatTableRow("U", "m/s", impeller.inlet.rms.U, 3);
                file << formatTableRow("Mach number", "m/s", impeller.inlet.rms.M, 4);
                file << formatTableRow("Relative Mach number", "m/s", impeller.inlet.rms.M_rel, 4);
                file << formatSubHeaderRow("Tip");
                file << formatTableRow("C", "m/s", impeller.inlet.tip.C, 3);
                file << formatTableRow("Cm", "m/s", impeller.inlet.tip.C_m, 3);
                file << formatTableRow("Cθ", "m/s", impeller.inlet.tip.C_theta, 3);
                file << formatTableRow("W", "m/s", impeller.inlet.tip.W, 3);
                file << formatTableRow("Wm", "m/s", impeller.inlet.tip.W_m, 3);
                file << formatTableRow("Wθ", "m/s", impeller.inlet.tip.W_theta, 3);
                file << formatTableRow("U", "m/s", impeller.inlet.tip.U, 3);
                file << formatTableRow("Mach number", "m/s", impeller.inlet.tip.M, 4);
                file << formatTableRow("Relative Mach number", "m/s", impeller.inlet.tip.M_rel, 4);

                file << formatSubHeaderRow("Impeller exit");
                file << formatTableRow("C", "m/s", impeller.outlet.tip.C, 3);
                file << formatTableRow("Cm", "m/s", impeller.outlet.tip.C_m, 3);
                file << formatTableRow("Cθ", "m/s", impeller.outlet.tip.C_theta, 3);
                file << formatTableRow("Slip velocity", "m/s", impeller.outlet.tip.C_slip, 3);
                file << formatTableRow("W", "m/s", impeller.outlet.tip.W, 3);
                file << formatTableRow("Wm", "m/s", impeller.outlet.tip.W_m, 3);
                file << formatTableRow("Wθ", "m/s", impeller.outlet.tip.W_theta, 3);
                file << formatTableRow("U", "m/s", impeller.outlet.tip.U, 3);
                file << formatTableRow("Absolute flow angle, α", "°", impeller.outlet.tip.alpha * (1.0 / DEG_RAD), 3);
                file << formatTableRow("Relative flow angle, β", "°", impeller.outlet.tip.beta * (1.0 / DEG_RAD), 3);
                file << formatTableRow("Mach number", "m/s", impeller.outlet.tip.M, 4);
                file << formatTableRow("Relative Mach number", "m/s", impeller.outlet.tip.M_rel, 4);
                file << fmt::format("\t\t</table>\n");

                file << fmt::format("\t\t<h3>Inlet shroud velocity triangle</h3>\n");
                file << fmt::format("\t\t<img src=\"vel_inlet_shroud.png\" alt=\"Inlet shroud velocity triangle.\">\n");
                file << fmt::format("\t\t<h3>Inlet hub velocity triangle</h3>\n");
                file << fmt::format("\t\t<img src=\"vel_inlet_hub.png\" alt=\"Inlet hub velocity triangle.\">\n");
                file << fmt::format("\t\t<h3>Outlet velocity triangle</h3>\n");
                file << fmt::format("\t\t<img src=\"vel_outlet.png\" alt=\"Outlet velocity triangle.\">\n");
                file << fmt::format("\t</body>\n");
                file << fmt::format("</html>\n");
            }
            file.close();
        }
    } catch (std::exception& e) {
        fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
    }
}