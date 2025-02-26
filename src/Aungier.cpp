#include "../include/Aungier.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          Aungier class
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Aungier::Aungier(const ImpellerState& inlet, const ImpellerState& outlet, const OperatingCondition& op,
                 const Geometry& geom, const bool optimizationFlag)
    : inlet(inlet),
      outlet(outlet),
      throat(inlet),
      op(op),
      geom(geom),
      pr_tt(0.0),
      isenEff(0.0),
      flowCoeff(0.0),
      workCoeff(0.0),
      Re_b2(0.0),
      Re_r2(0.0),
      dH0(0.0),
      optimizationFlag(optimizationFlag) {}

Aungier::Aungier(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op,
                 const bool optimizationFlag)
    : inlet(thermo),
      outlet(thermo),
      throat(thermo),
      op(op),
      geom(geom),
      pr_tt(0.0),
      isenEff(0.0),
      flowCoeff(0.0),
      workCoeff(0.0),
      Re_b2(0.0),
      Re_r2(0.0),
      dH0(0.0),
      optimizationFlag(optimizationFlag) {}

Aungier::Aungier(const Impeller& impeller, const bool optimizationFlag)
    : inlet(impeller.inlet),
      outlet(impeller.outlet),
      throat(impeller.inlet),
      op(impeller.op),
      geom(impeller.geom),
      pr_tt(impeller.pr_tt),
      isenEff(impeller.isenEff),
      flowCoeff(impeller.flowCoeff),
      workCoeff(impeller.workCoeff),
      Re_b2(impeller.Re_b2),
      Re_r2(impeller.Re_r2),
      dH0(impeller.dH0),
      optimizationFlag(optimizationFlag) {}

void Aungier::runCalculations() {
    inletCalcs();
    inletInducerOptimization();
    throatCalcs();
}

void Aungier::inletCalcs() {
    int maxIterations = 100;
    int iteration = 0;
    double error = 1.0;
    double r1h = geom.r1h * MM_M;
    double r1t = geom.r1t * MM_M;
    double b1 = r1t - r1h;
    double r1 = getRMS(r1t, r1h);
    double rho1_guess = inlet.total.props.D;
    double tol = 10e-6;
    inlet.static_ = inlet.total;

    // Set blade thickness to 2mm by default
    double tb1 = 2.0 * MM_M;
    std::string temp;
    fmt::print("Enter the inlet blade thickness in mm: ");
    std::getline(std::cin, temp);
    try {
        tb1 = std::stod(temp);
    } catch (const std::invalid_argument&) {
        fmt::println("No conversion could be performed");
    } catch (const std::out_of_range&) {
        fmt::println("Could not convert string to double, value is out of range");
    }
    // Use the updated blade thickness from user to calculate the inlet area
    double A1 = PI * (std::pow(r1t, 2) - std::pow(r1h, 2)) - b1 * geom.ZFull * (tb1 * MM_M);

    Impeller::printBorder("Aungier", tol, maxIterations, INLET);
    do {
        iteration++;
        rho1_guess = inlet.static_.props.D;
        inlet.rms.C_m = op.mfr / (rho1_guess * A1 * (1.0 - geom.blockage1));
        if (op.alpha == 0) {
            inlet.rms.C = inlet.rms.C_m;
            // fmt::println("inlet.rms.C: {}", inlet.rms.C);
        } else {
            inlet.rms.C = inlet.rms.C_m / std::tan(op.alpha * DEG_RAD);
        }
        inlet.static_.props.H = inlet.total.props.H - (std::pow(inlet.rms.C, 2) / 2.0);
        inlet.static_.set_props("HS", inlet.static_.props.H, inlet.static_.props.S);
        error = std::abs((inlet.static_.props.D / rho1_guess) - 1.0);
        Impeller::printIteration(iteration, inlet.static_.props.P, inlet.static_.props.T, inlet.static_.props.D, error,
                                 "density");
    } while (error > tol && iteration <= maxIterations);

    inlet.hub.C = inlet.rms.C;
    inlet.tip.C = inlet.rms.C;
    inlet.hub.C_m = inlet.rms.C_m;
    inlet.tip.C_m = inlet.rms.C_m;
    inlet.hub.M = inlet.rms.M;
    inlet.tip.M = inlet.rms.M;

    if (error <= tol) {
        fmt::print(fg(fmt::color::green), "Solution converged!\n");
    } else {
        fmt::print(fg(fmt::color::red), "Solution did not converge, results could be inaccurate!\n");
    }
}

void Aungier::inletInducerOptimization() {
    if (!optimizationFlag) {
        return;  // if flag was set to false, just do nothing with this function
    }

    std::array<AungierOptData, 11> optData;
    for (auto& data : optData) {
        data.thermo = inlet.total;
        data.vel = inlet.rms;
    }

    // Cm Range should be 50-250 m/s
    std::array<double, 11> CmRange;
    for (int i = 0; i < CmRange.size(); i++) {
        CmRange[i] = i * 20 + 50;
    }

    for (int i = 0; i < CmRange.size(); i++) {
        optData[i].vel.C_m = CmRange[i];
        optData[i].vel.C_theta = CmRange[i] * std::tan(op.alpha * DEG_RAD);
        optData[i].vel.C = std::sqrt(CmRange[i] * CmRange[i] + optData[i].vel.C_theta * optData[i].vel.C_theta);
        optData[i].thermo.props.T =
            inlet.total.props.T - (optData[i].vel.C * optData[i].vel.C) / (2.0 * optData[i].thermo.props.CP);
        double k = optData[i].thermo.props.Y / (optData[i].thermo.props.Y - 1.0);
        optData[i].thermo.props.P = inlet.total.props.P * std::pow(optData[i].thermo.props.T / inlet.total.props.T, k);
        optData[i].thermo.set_props("PT", optData[i].thermo.props.P, optData[i].thermo.props.T);
        optData[i].vel.M = optData[i].vel.C / optData[i].thermo.props.A;
        optData[i].area = op.mfr / optData[i].thermo.props.D / optData[i].vel.C_m / (1 - geom.blockage1);
        optData[i].r1s = std::sqrt(optData[i].area / PI + ((geom.r1h * MM_M) * (geom.r1h * MM_M)));
        optData[i].vel.U = op.omega * optData[i].r1s;
        optData[i].vel.W =
            std::sqrt(optData[i].vel.C_m * optData[i].vel.C_m +
                      ((optData[i].vel.U * optData[i].vel.U) - (optData[i].vel.C_theta * optData[i].vel.C_theta)));
        fmt::println("r1s: {:.4f}, Cm: {:.3f}, W1s: {:.3f}", optData[i].r1s, optData[i].vel.C_m, optData[i].vel.W);
    }

    namespace fs = std::filesystem;
    fs::path dirPath = "./../tmp";
    try {
        if (!fs::exists(dirPath)) {
            fs::create_directory(dirPath);
        }

        std::ofstream file(dirPath.string() + "/optVelocities.csv", std::ios::out);
        if (file.is_open()) {
            file << fmt::format("r1s,Cm,W1s\n");
            for (auto data : optData) {
                file << fmt::format("{},{},{}\n", data.r1s, data.vel.C_m, data.vel.W);
            }
        }
        file.close();
        auto result = std::system("python ./../scripts/inletOptimizationAungier.py");
        if (result == 0) {
            fmt::print(fg(fmt::color::green), "Python script for velocity triangles ran successfully\n");
        } else {
            fmt::print(fg(fmt::color::red), "Error running velocity triangle script\n");
        }

        fs::remove_all(dirPath);
    } catch (std::exception& e) {
        fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
    }
}

void Aungier::calculateInletVelocities() {
    std::array<Velocities, 3> velArr = {inlet.hub, inlet.rms, inlet.tip};
    std::array<double, 3> radius = {geom.r1h * MM_M, geom.r1rms * MM_M, geom.r1t * MM_M};
    std::array<double, 3> angles = {geom.beta1h, geom.beta1rms, geom.beta1t};

    for (int i = 0; i < velArr.size(); i++) {
        velArr[i].alpha = op.alpha;
        velArr[i].U = (radius[i]) * op.omega;

        velArr[i].C_theta = velArr[i].C * std::sin(velArr[i].alpha);
        velArr[i].W_m = velArr[i].C_m;
        velArr[i].W_theta = velArr[i].U - velArr[i].C_theta;
        velArr[i].W = sqrt(pow(velArr[i].W_m, 2.) + pow(velArr[i].W_theta, 2.));
        velArr[i].beta = -std::atan(velArr[i].W_theta / velArr[i].W_m);

        velArr[i].M_rel = velArr[i].W / inlet.static_.props.A;
        velArr[i].inc = angles[i] - velArr[i].beta * (1.0 / DEG_RAD);
    }

    inlet.hub = velArr[0];
    inlet.rms = velArr[1];
    inlet.tip = velArr[2];
}

void Aungier::throatCalcs() {
    double x = ((geom.r2 - geom.r1h) * (1.0 / std::cos(PI / 18.0) - std::tan(PI / 18.0)) - geom.b2) * MM_M;
    double y = (geom.r2 - geom.r1t) * MM_M;
    fmt::println("x: {}, y: {}", x, y);

    // mj is the value of m at the junction of arcs
    double mj_m4 = (x + (1.0 - sqrt(2.0)) * y) / ((2.0 - sqrt(2.0)) * (x + y));
    fmt::println("mj/mtip: {}", mj_m4);

    double m4_tip = (PI / 4.0) * (x + y);
    double m4_hub = ((4 * PI / 9.0) * (1.0 / std::cos(PI / 18.0)) * (geom.r2 - geom.r1h)) * MM_M;
    fmt::println("m4_tip: {}, m4_hub: {}", m4_tip, m4_hub);

    // set up arrays for meridional distance from 0.0-1.0 and respective calculations
    std::array<double, 101> m_m4;
    std::array<double, 101> phi_tip;
    std::array<double, 101> phi_hub;
    std::array<double, 101> r_tip;
    std::array<double, 101> r_hub;

    for (int i = 0; i < m_m4.size(); i++) {
        m_m4[i] = i / 100.0;
        phi_hub[i] = (PI / 18.0) + (4 * PI / 9.0) * m_m4[i];
        r_hub[i] = (geom.r1h + (geom.r2 - geom.r1h) * (1.0 - (std::cos(phi_hub[i]) / std::cos(PI / 18.0)))) * MM_M;
        if (m_m4[i] < mj_m4) {
            phi_tip[i] = (PI / 4.0) * (m_m4[i] / mj_m4);
            r_tip[i] =
                (geom.r1t * MM_M) + (1.0 - std::cos(phi_tip[i])) * ((x + (1.0 - sqrt(2.0)) * y) / (2.0 - sqrt(2.0)));
        } else {
            phi_tip[i] = (PI / 4.0) * (1.0 + ((m_m4[i] - mj_m4) / (1.0 - mj_m4)));
            r_tip[i] = (geom.r1t * MM_M) +
                       (1.0 - std::cos(PI / 4.0)) * ((x + (1.0 - sqrt(2.0)) * y) / (2.0 - sqrt(2.0))) +
                       (std::cos(PI / 4.0) - std::cos(phi_tip[i])) * ((y + (1.0 - sqrt(2.0)) * x) / (2.0 - sqrt(2.0)));
        }
    }
}