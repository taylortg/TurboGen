#include "../include/impeller.h"

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <optional>
#include <sstream>
#include <string>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          ImpellerState class
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ImpellerState::ImpellerState(const std::string& fluid_name)
    : total(fluid_name), static_(fluid_name), isentropic(fluid_name), hub(), rms(), tip(), geom() {}

ImpellerState::ImpellerState(const ThermoProps& thermo, const Geometry& geom)
    : total(thermo), static_(thermo), isentropic(thermo), hub(), rms(), tip(), geom(geom) {}

// Copy constructor
ImpellerState::ImpellerState(const ImpellerState& other)
    : total(other.total),
      static_(other.static_),
      isentropic(other.isentropic),
      hub(other.hub),
      rms(other.rms),
      tip(other.tip),
      geom(other.geom) {}

// Assignment operator
ImpellerState& ImpellerState::operator=(const ImpellerState& other) {
    if (this != &other) {
        total = other.total;
        static_ = other.static_;
        isentropic = other.isentropic;
        hub = other.hub;
        rms = other.rms;
        tip = other.tip;
        geom = other.geom;
    }
    return *this;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          Impeller class
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Impeller::Impeller(ImpellerState& inlet, ImpellerState& outlet, const OperatingCondition& op)
    : inlet(inlet), outlet(outlet), op(op) {}

Impeller::Impeller(ThermoProps thermo, Geometry geom, const OperatingCondition& op)
    : inlet(thermo, geom), outlet(thermo, geom), op(op) {}

void Impeller::calculateInletCondition(std::string solverType) {
    // To-do: catch if user wants to continue when solution does not converge
    // use if (inletSolverResult) {}
    if (solverType == "Japikse") {
        auto inletSolverResult = inletJapikseSolver(0.3, 1000, 1e-8);
    } else {
        fmt::print(fg(fmt::color::red), "Solver type not implemented. Solver type passed: {}\n", solverType);
        exit(EXIT_FAILURE);
    }
}

std::optional<double> Impeller::inletJapikseSolver(double Mguess = 0.3, int maxIterations = 1000,
                                                   double tolerance = 1e-6) {
    int iteration = 1;
    double error = 1.0;

    std::string border = fmt::format("{:=<90}\n", "");
    fmt::print(
        "{}Entered impeller inlet solver loop!\nSolver: Japikse\nConvergence tolerance: {:.2e}\nMaximum iterations: "
        "{}\n{}",
        border, tolerance, maxIterations, border);

// This is only here to indicate if there are issues with the solver. The
// execution time of the while loop should be incredibly small. If you see it is
// taking ms to run, there might be a problem.
#ifdef DEBUG
    auto start = std::chrono::high_resolution_clock::now();
#endif

    do {
        inlet.static_.props.T = inlet.total.props.T / (1. + (inlet.static_.props.Y - 1.) / 2. * pow(Mguess, 2));
        double k = inlet.static_.props.Y / (inlet.static_.props.Y - 1.);
        inlet.static_.props.P = inlet.total.props.P / (pow((inlet.total.props.T / inlet.static_.props.T), k));
        inlet.static_.set_props("PT", inlet.static_.props.P, inlet.static_.props.T);
        double Rgas = inlet.static_.props.CP - inlet.static_.props.CV;
        inlet.tip.C_m = (op.mfr * Rgas * inlet.static_.props.T) / (inlet.static_.props.P * inlet.geom.area1);
        inlet.tip.M = inlet.tip.C_m / inlet.static_.props.A;
        error = std::abs((Mguess - inlet.tip.M) / Mguess);
        Mguess = inlet.tip.M;

        fmt::print("Iteration {}\n", iteration);
        fmt::print("Pressure: {:<10.0f}Temperature: {:<10.2f}Mach: {:<10.6f}Error:{:<10.4e}\n", inlet.static_.props.P,
                   inlet.static_.props.T, inlet.tip.M, error);

        iteration++;
    } while (iteration <= maxIterations && error > tolerance);

#ifdef DEBUG
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Solver execution time: " << std::setw(10) << std::fixed << std::setprecision(6)
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
#endif

    fmt::print("{}", border);
    inlet.rms.C_m = inlet.tip.C_m;
    inlet.hub.C_m = inlet.tip.C_m;
    inlet.rms.M = inlet.tip.M;
    inlet.hub.M = inlet.tip.M;

    calculateInletVelocities();

    if (error <= tolerance) {
        fmt::print(fg(fmt::color::green), "Solution converged!\n");
        return inlet.static_.props.P;
    } else {
        fmt::print(fg(fmt::color::red), "Solution did not converge, results could be inaccurate!\n");
        return std::nullopt;
    }
}

void Impeller::calculateInletVelocities() {
    std::array<Velocities, 3> velArr = {this->inlet.hub, this->inlet.rms, this->inlet.tip};
    std::array<double, 3> radius = {this->inlet.geom.r1h * MM_M, this->inlet.geom.r1rms * MM_M,
                                    this->inlet.geom.r1t * MM_M};
    std::array<double, 3> angles = {this->inlet.geom.beta1h, this->inlet.geom.beta1rms, this->inlet.geom.beta1t};

    for (int i = 0; i < velArr.size(); i++) {
        velArr[i].alpha = this->op.alpha;
        velArr[i].U = (radius[i]) * this->op.omega;
        if (velArr[i].alpha == 0) {
            velArr[i].C = velArr[i].C_m;
            velArr[i].C_theta = 0.0;
            velArr[i].W = sqrt(pow(velArr[i].C, 2.) + pow(velArr[i].U, 2.));
            velArr[i].W_m = velArr[i].C_m;
            velArr[i].W_theta = velArr[i].U;
            velArr[i].beta = -std::atan(velArr[i].U / velArr[i].C);
        } else {
            velArr[i].C = velArr[i].C_m / std::cos(velArr[i].alpha);
            velArr[i].C_theta = velArr[i].C * std::sin(velArr[i].alpha);
            velArr[i].W_m = velArr[i].C_m;
            velArr[i].W_theta = velArr[i].U = velArr[i].C_theta;
            velArr[i].W = sqrt(pow(velArr[i].W_m, 2.) + pow(velArr[i].W_theta, 2.));
            velArr[i].beta = -std::atan(velArr[i].W_theta / velArr[i].W_m);
        }
        velArr[i].M_rel = velArr[i].W / this->inlet.static_.props.A;
        velArr[i].inc = angles[i] - velArr[i].beta * (1.0 / DEG_RAD);
    }

    this->inlet.hub = velArr[0];
    this->inlet.rms = velArr[1];
    this->inlet.tip = velArr[2];
}

void Impeller::calculateOutletCondition(std::string solverType, std::string slipModel) {
    calculateSlip(slipModel);
    if (solverType == "Japikse") {
        auto outletSolverResult = outletJapikseSolver(1000, 1e-8, 0.92);
    } else {
        fmt::print(fg(fmt::color::red), "Solver type not implemented. Solver type passed: {}\n", solverType);
        exit(EXIT_FAILURE);
    }
    printOutputFile();
}

std::optional<double> Impeller::outletJapikseSolver(int maxIterations = 1000, double tolerance = 1e-6,
                                                    double rotorEfficiency = 0.92) {
    int iteration = 1;
    double error = 1.0;
    std::array<Velocities, 3> velArr = {this->outlet.hub, this->outlet.rms, this->outlet.tip};
    std::array<double, 3> dh0 = {};
    double Rgas = 0.0;
    this->outlet.total = this->inlet.total;
    this->outlet.static_ = this->inlet.static_;

    std::string border = fmt::format("{:=<90}\n", "");
    fmt::print(
        "{}Entered impeller outlet solver loop!\nSolver: Japikse\nConvergence tolerance: {:.2e}\nMaximum iterations: "
        "{}\n{}",
        border, tolerance, maxIterations, border);

    for (int i = 0; i < velArr.size(); ++i) {
        velArr[i].U = this->op.omega * (this->outlet.geom.r2 * MM_M);
        // Initial guess of M and static props will not account for backsweep
        velArr[i].C_theta = velArr[i].sigma * velArr[i].U;
        dh0[i] = velArr[i].U * velArr[i].C_theta;
        this->outlet.total.props.T = this->inlet.total.props.T + dh0[i] / this->outlet.total.props.CP;
        this->outlet.total.props.P =
            this->inlet.total.props.P *
            pow(1.0 + ((rotorEfficiency * dh0[i]) / (this->outlet.total.props.CP * this->inlet.total.props.T)),
                this->outlet.total.props.Y / (this->outlet.total.props.Y - 1.0));
        this->outlet.total.set_props("PT", this->outlet.total.props.P, this->outlet.total.props.T);
        this->outlet.static_.set_props(
            "PT", this->outlet.total.props.P,
            this->outlet.total.props.T);  // do this only for the initial guess at Mach number
    }

    const double initialMachGuess = 0.6;
    double M2guess = initialMachGuess;
    do {
        this->outlet.static_.props.T =
            this->outlet.total.props.T / (1. + (this->outlet.static_.props.Y - 1.) / 2. * pow(M2guess, 2));
        double k = this->outlet.static_.props.Y / (this->outlet.static_.props.Y - 1.);
        this->outlet.static_.props.P =
            this->outlet.total.props.P / (pow((this->outlet.total.props.T / this->outlet.static_.props.T), k));
        try {
            this->outlet.static_.set_props("PT", this->outlet.static_.props.P, this->outlet.static_.props.T);
        } catch (const std::runtime_error& e) {
            fmt::print(fg(fmt::color::red), "Runtime error: {}\n", e.what());
            return std::nullopt;
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
            return std::nullopt;
        }
        Rgas = this->outlet.static_.props.CP - this->outlet.static_.props.CV;
        for (int i = 0; i < velArr.size(); i++) {
            if (this->outlet.static_.props.P == 0 || this->outlet.geom.area2 == 0) {
                fmt::print(fg(fmt::color::red),
                           "Outlet static P or outlet area are set to 0\noutlet.static_.props.P: "
                           "{:<10.0f}\ngeom.area2: {:<10.4f}",
                           this->outlet.static_.props.P, this->outlet.geom.area2);
                return std::nullopt;
            }
            try {
                velArr[i].C_m = (this->op.mfr * Rgas * this->outlet.static_.props.T) /
                                (this->outlet.static_.props.P * this->outlet.geom.area2);
            } catch (const std::exception& e) {
                fmt::print(fg(fmt::color::red), "Error in C_m calculation: {}\n", e.what());
                return std::nullopt;
            }
            velArr[i].C = sqrt(pow(velArr[i].C_m, 2.) + pow(velArr[i].C_theta, 2.));
            velArr[i].M = velArr[i].C / this->outlet.static_.props.A;
            // Update C_theta since C_m should be more accurate and account for backsweep
            if (this->outlet.geom.beta2 == 0) {
                velArr[i].C_theta = velArr[i].sigma * velArr[i].U;
            } else {
                velArr[i].C_theta =
                    velArr[i].sigma * velArr[i].U + velArr[i].C_m * tan(this->outlet.geom.beta2 * DEG_RAD);
            }
            dh0[i] = velArr[i].U * velArr[i].C_theta;
            this->outlet.total.props.T = this->inlet.total.props.T + dh0[i] / this->outlet.total.props.CP;
            this->outlet.total.props.P =
                this->inlet.total.props.P *
                pow(1.0 + ((rotorEfficiency * dh0[i]) / (this->outlet.total.props.CP * this->inlet.total.props.T)),
                    this->outlet.total.props.Y / (this->outlet.total.props.Y - 1.0));
            this->outlet.total.set_props("PT", this->outlet.total.props.P, this->outlet.total.props.T);
        }
        error = fabs(M2guess - velArr[2].M) / M2guess;
        M2guess = velArr[2].M;

        fmt::print("Iteration {}\n", iteration);
        fmt::print("Pressure: {:<10.0f}Temperature: {:<10.2f}Mach: {:<10.6f}Error:{:<10.4e}\n",
                   this->outlet.static_.props.P, this->outlet.static_.props.T, M2guess, error);

        iteration++;
    } while (iteration <= maxIterations && error > tolerance);
    fmt::print("{}", border);

    this->outlet.hub = velArr[0];
    this->outlet.rms = velArr[1];
    this->outlet.tip = velArr[2];

    calculateOutletVelocities();
    this->pr_tt = this->outlet.total.props.P / this->inlet.total.props.P;
    this->dH0 = this->outlet.total.props.H - this->inlet.total.props.H;
    this->workCoefficient = this->dH0 / pow(this->outlet.tip.U, 2);
    this->Re_b2 = this->outlet.total.props.D * this->outlet.tip.C_m * this->outlet.geom.b2 / this->outlet.total.props.V;
    this->Re_r2 = this->outlet.total.props.D * this->outlet.tip.C_m * this->outlet.geom.r2 / this->outlet.total.props.V;

    if (error <= tolerance) {
        fmt::print(fg(fmt::color::green), "Solution converged!\n");
        return this->inlet.static_.props.P;
    } else {
        fmt::print(fg(fmt::color::red), "Solution did not converge, results could be inaccurate!\n");
        return std::nullopt;
    }
}

void Impeller::calculateOutletVelocities() {
    std::array<Velocities, 3> vel = {this->outlet.hub, this->outlet.rms, this->outlet.tip};

    for (int i = 0; i < vel.size(); i++) {
        vel[i].W_theta = vel[i].U - vel[i].C_theta;
        vel[i].W_m = vel[i].C_m;
        vel[i].W = sqrt(pow(vel[i].W_theta, 2.) + pow(vel[i].W_m, 2.));
        vel[i].beta = -atan(vel[i].W_theta / vel[i].W_m);
        vel[i].M_rel = vel[i].W / this->outlet.static_.props.A;
        vel[i].C_slip = vel[i].U - vel[i].C_theta;
        vel[i].alpha = atan(vel[i].C_theta / vel[i].C_m);
    }

    this->outlet.hub = vel[0];
    this->outlet.rms = vel[1];
    this->outlet.tip = vel[2];
}

void Impeller::calculateSlip(std::string slipModel) {
    if (slipModel == "Wiesner") {
        double sigma = 1.0 - sqrt(cos(this->outlet.geom.beta2 * DEG_RAD)) / pow(this->outlet.geom.Z, 0.7);
        this->outlet.hub.sigma = sigma;
        this->outlet.rms.sigma = sigma;
        this->outlet.tip.sigma = sigma;
    } else {
        fmt::print(fg(fmt::color::red), "Slip model not implemented. Slip model passed: {}\n", slipModel);
        exit(EXIT_FAILURE);
    }
}

void Impeller::printOutputFile() {
    namespace fs = std::filesystem;
    fs::path dirPath = "./../out";

    try {
        if (!fs::exists(dirPath)) {
            fs::create_directory(dirPath);
        }
        std::ofstream file(dirPath.string() + "/out.txt", std::ios::out);
        if (file.is_open()) {
            file << fmt::format("Impeller leading edge:\n");
            file << fmt::format("{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n", "P(Pa)", "T(K)", "H0(kJ/kg)", "s0(kJ/kg-K)",
                                "A(mm^2)", "Athroat(mm^2)");
            file << fmt::format("{:<15.0f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}\n\n",
                                this->inlet.static_.props.P, this->inlet.static_.props.T, this->inlet.total.props.H,
                                this->inlet.total.props.S, (this->inlet.geom.area1 * 1000000.0),
                                (this->inlet.geom.throatArea * 1000000.0));
            file << fmt::format("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<12}{:<12}{:<12}\n", "Station", "Dia(mm)",
                                "C(m/s)", "W(m/s)", "U(m/s)", "Mrel", "Beta(deg)", "Beta'(deg)", "Inc(deg)");

            const std::array<Velocities, 3> vel = {this->inlet.hub, this->inlet.rms, this->inlet.tip};
            const std::array<double, 3> diameters = {this->inlet.geom.r1h, this->inlet.geom.r1rms,
                                                     this->inlet.geom.r1t};
            const std::array<double, 3> angles = {this->inlet.geom.beta1h, this->inlet.geom.beta1rms,
                                                  this->inlet.geom.beta1t};
            const std::array<std::string, 3> names = {"Hub:", "RMS:", "Shroud:"};
            for (size_t i = 0; i < vel.size(); i++) {
                file << fmt::format("{:<10}{:<10.1f}{:<10.2f}{:<10.2f}{:<10.2f}{:<10.4f}{:<12.2f}{:<12.2f}{:<12.2f}\n",
                                    names[i], diameters[i], vel[i].C, vel[i].W, vel[i].U, vel[i].M_rel, angles[i],
                                    (vel[i].beta * (1 / DEG_RAD)), vel[i].inc);
            }
            file << fmt::format("\nImpeller exit\n");
            file << fmt::format("{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n", "Dia(mm)", "Tip width(mm)",
                                "P0(Pa)", "P(Pa)", "T0(K)", "H0(kJ/kg)", "s0(kJ/kg-K)", "U(m/s)", "M_U");
            file << fmt::format("{:<15.2f}{:<15.3f}{:<15.0f}{:<15.0f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.4f}\n",
                                this->outlet.geom.r2, this->outlet.geom.b2, this->outlet.total.props.P,
                                this->outlet.static_.props.P, this->outlet.total.props.T, this->outlet.total.props.H,
                                this->outlet.total.props.S, this->outlet.tip.U,
                                (this->outlet.tip.U / this->outlet.total.props.A));
            file << fmt::format("\n{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n", "M_rms", "W_rms(m/s)", "C_rms(m/s)",
                                "Alpha(deg)", "Beta(deg)", "W ratio");
            file << fmt::format("{:<15.3f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.3f}\n", this->outlet.rms.M,
                                this->outlet.rms.W, this->outlet.rms.C, this->outlet.rms.alpha * (1 / DEG_RAD),
                                this->outlet.rms.beta * (1 / DEG_RAD), (this->outlet.rms.W / this->inlet.rms.W));

            file << fmt::format("\nOverall Performance\n");
            file << fmt::format("{:<15}{:<15}{:<15}{:<15}{:<15}\n", "pr_tt", "dH/U^2", "dH0(kJ/kg)", "Re_b2", "Re_r2");
            file << fmt::format("{:<15.3f}{:<15.4f}{:<15.2f}{:<15.3e}{:<15.3e}\n", this->pr_tt, this->workCoefficient,
                                this->dH0, this->Re_b2, this->Re_r2);
            file.close();
        } else {
            fmt::print(fg(fmt::color::red), "Error creating output file.\n");
        }
    } catch (const fs::filesystem_error& e) {
        fmt::print(fg(fmt::color::red), "Error create output file: {}\n", e.what());
    }
}
