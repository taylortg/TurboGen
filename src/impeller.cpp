#include "../include/impeller.h"

#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <utility>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          ImpellerState class
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ImpellerState::ImpellerState(const std::string& fluid_name)
    : total(fluid_name), static_(fluid_name), isentropic(fluid_name), hub(), rms(), tip() {}

ImpellerState::ImpellerState(const ThermoProps& thermo)
    : total(thermo), static_(thermo), isentropic(thermo), hub(), rms(), tip() {}

// Copy constructor
ImpellerState::ImpellerState(const ImpellerState& other)
    : total(other.total),
      static_(other.static_),
      isentropic(other.isentropic),
      hub(other.hub),
      rms(other.rms),
      tip(other.tip) {}

// Assignment operator
ImpellerState& ImpellerState::operator=(const ImpellerState& other) {
    if (this != &other) {
        total = other.total;
        static_ = other.static_;
        isentropic = other.isentropic;
        hub = other.hub;
        rms = other.rms;
        tip = other.tip;
    }
    return *this;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          Impeller class
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Impeller::Impeller(const ImpellerState& inlet, const ImpellerState& outlet, const OperatingCondition& op,
                   const Geometry& geom)
    : inlet(inlet),
      outlet(outlet),
      op(op),
      geom(geom),
      pr_tt(0.0),
      isenEff(0.0),
      flowCoeff(0.0),
      workCoeff(0.0),
      Re_b2(0.0),
      Re_r2(0.0),
      dH0(0.0) {}

Impeller::Impeller(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op)
    : inlet(thermo),
      outlet(thermo),
      op(op),
      geom(geom),
      pr_tt(0.0),
      isenEff(0.0),
      flowCoeff(0.0),
      workCoeff(0.0),
      Re_b2(0.0),
      Re_r2(0.0),
      dH0(0.0) {}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          Driver function definitions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void Impeller::calculateInletCondition(std::string solverType) {
    // To-do: catch if user wants to continue when solution does not converge
    // use if (inletSolverResult) {}
    if (solverType == "Japikse") {
        inletJapikseSolver(0.3, 1000, 1e-8);
    } else if (solverType == "Aungier") {
        inletAungierSolver(1000, 1e-8);
    } else {
        fmt::print(fg(fmt::color::red), "Solver type not implemented. Solver type passed: {}\n", solverType);
        exit(EXIT_FAILURE);
    }
}

void Impeller::calculateOutletCondition(std::string solverType, std::string slipModel) {
    calculateSlip(std::move(slipModel));
    if (solverType == "Japikse") {
        outletJapikseSolver(1000, 1e-8, 0.92);
    } else if (solverType == "Aungier") {
        outletAungierSolver(1000, 1e-8, 0.9);
    } else {
        fmt::print(fg(fmt::color::red), "Solver type not implemented. Solver type passed: {}\n", solverType);
        exit(EXIT_FAILURE);
    }
    printOutputFile();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          Japikse function definitions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::optional<double> Impeller::inletJapikseSolver(double Mguess = 0.3, int maxIterations = 1000,
                                                   double tolerance = 1e-6) {
    int iteration = 1;
    double error;

    printBorder("Japikse", tolerance, maxIterations, INLET);

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
        try {
            inlet.static_.set_props("PT", inlet.static_.props.P, inlet.static_.props.T);
        } catch (const std::runtime_error& e) {
            fmt::print(fg(fmt::color::red), "Runtime error: {}\n", e.what());
            return std::nullopt;
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
            return std::nullopt;
        }
        double Rgas = inlet.static_.props.CP - inlet.static_.props.CV;
        inlet.tip.C = (op.mfr * Rgas * inlet.static_.props.T) / (inlet.static_.props.P * geom.area1);
        inlet.tip.M = inlet.tip.C / inlet.static_.props.A;
        error = std::abs((Mguess - inlet.tip.M) / Mguess);
        Mguess = inlet.tip.M;

        printIteration(iteration, inlet.static_.props.P, inlet.static_.props.T, inlet.tip.M, error, "Mach");

        iteration++;
    } while (iteration <= maxIterations && error > tolerance);

#ifdef DEBUG
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Solver execution time: " << std::setw(10) << std::fixed << std::setprecision(6)
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
#endif

    fmt::print("{:=<90}\n", "");
    inlet.rms.C = inlet.tip.C;
    inlet.hub.C = inlet.tip.C;
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

std::optional<double> Impeller::outletJapikseSolver(int maxIterations = 1000, double tolerance = 1e-6,
                                                    double rotorEfficiency = 0.92) {
    int iteration = 1;
    double error;
    std::array<Velocities, 3> velArr = {outlet.hub, outlet.rms, outlet.tip};
    std::array<double, 3> dh0 = {};
    double Rgas = 0.0;
    outlet.total = inlet.total;
    outlet.static_ = inlet.static_;

    printBorder("Japikse", tolerance, maxIterations, OUTLET);

    for (int i = 0; i < velArr.size(); ++i) {
        velArr[i].U = op.omega * (geom.r2 * MM_M);
        velArr[i].C_theta = velArr[i].sigma * velArr[i].U;
        dh0[i] = velArr[i].U * velArr[i].C_theta;
        outlet.total.props.T = inlet.total.props.T + dh0[i] / outlet.total.props.CP;
        outlet.total.props.P = inlet.total.props.P *
                               pow(1.0 + ((rotorEfficiency * dh0[i]) / (outlet.total.props.CP * inlet.total.props.T)),
                                   outlet.total.props.Y / (outlet.total.props.Y - 1.0));
    }

    outlet.total.set_props("PT", outlet.total.props.P, outlet.total.props.T);
    outlet.static_.set_props("PT", outlet.total.props.P,
                             outlet.total.props.T);  // do this only for the initial guess at Mach number

    constexpr double initialMachGuess = 0.6;
    double M2guess = initialMachGuess;
    do {
        outlet.static_.props.T = outlet.total.props.T / (1. + (outlet.static_.props.Y - 1.) / 2. * pow(M2guess, 2));
        double k = outlet.static_.props.Y / (outlet.static_.props.Y - 1.);
        outlet.static_.props.P = outlet.total.props.P / (pow((outlet.total.props.T / outlet.static_.props.T), k));
        try {
            outlet.static_.set_props("PT", outlet.static_.props.P, outlet.static_.props.T);
        } catch (const std::runtime_error& e) {
            fmt::print(fg(fmt::color::red), "Runtime error: {}\n", e.what());
            return std::nullopt;
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
            return std::nullopt;
        }
        Rgas = outlet.static_.props.CP - outlet.static_.props.CV;
        for (int i = 0; i < velArr.size(); i++) {
            if (outlet.static_.props.P == 0 || geom.area2 == 0) {
                fmt::print(fg(fmt::color::red),
                           "Outlet static P or outlet area are set to 0\noutlet.static_.props.P: "
                           "{:<10.0f}\ngeom.area2: {:<10.4f}",
                           outlet.static_.props.P, geom.area2);
                return std::nullopt;
            }

            velArr[i].C_m = (op.mfr * Rgas * outlet.static_.props.T) / (outlet.static_.props.P * geom.area2);
            velArr[i].C = sqrt(pow(velArr[i].C_m, 2.) + pow(velArr[i].C_theta, 2.));
            velArr[i].M = velArr[i].C / outlet.static_.props.A;
            if (geom.beta2 == 0) {
                velArr[i].C_theta = velArr[i].sigma * velArr[i].U;
            } else {
                velArr[i].C_theta = velArr[i].sigma * velArr[i].U + velArr[i].C_m * tan(geom.beta2 * DEG_RAD);
            }
            dh0[i] = velArr[i].U * velArr[i].C_theta;
            outlet.total.props.T = inlet.total.props.T + dh0[i] / outlet.total.props.CP;
            outlet.total.props.P =
                inlet.total.props.P *
                pow(1.0 + ((rotorEfficiency * dh0[i]) / (outlet.total.props.CP * inlet.total.props.T)),
                    outlet.total.props.Y / (outlet.total.props.Y - 1.0));
            outlet.total.set_props("PT", outlet.total.props.P, outlet.total.props.T);
        }
        error = fabs(M2guess - velArr[2].M) / M2guess;
        M2guess = velArr[2].M;

        printIteration(iteration, outlet.static_.props.P, outlet.static_.props.T, M2guess, error, "Mach");

        iteration++;
    } while (iteration <= maxIterations && error > tolerance);
    fmt::print("{:=<90}\n", "");

    outlet.hub = velArr[0];
    outlet.rms = velArr[1];
    outlet.tip = velArr[2];

    calculateOutletVelocities();
    pr_tt = outlet.total.props.P / inlet.total.props.P;
    dH0 = outlet.total.props.H - inlet.total.props.H;
    workCoeff = dH0 / pow(outlet.tip.U, 2);
    Re_b2 = outlet.total.props.D * outlet.tip.C_m * geom.b2 / outlet.total.props.V;
    Re_r2 = outlet.total.props.D * outlet.tip.C_m * geom.r2 / outlet.total.props.V;

    if (error <= tolerance) {
        fmt::print(fg(fmt::color::green), "Solution converged!\n");
        return inlet.static_.props.P;
    } else {
        fmt::print(fg(fmt::color::red), "Solution did not converge, results could be inaccurate!\n");
        return std::nullopt;
    }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          CCMD functions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
std::optional<double> Impeller::inletAungierSolver(int maxIterations, double tolerance) {
    int iteration = 1;
    double error;

    printBorder("Aungier", tolerance, maxIterations, INLET);

#ifdef DEBUG
    auto start = std::chrono::high_resolution_clock::now();
#endif

    inlet.static_ = inlet.total;
    double area = geom.area1;
    double rhoOld = inlet.static_.props.D;
    do {
        inlet.tip.C = op.mfr / (inlet.static_.props.D * area * std::cos(op.alpha * DEG_RAD));
        inlet.tip.M = inlet.tip.C / inlet.static_.props.A;
        inlet.static_.props.T = inlet.total.props.T / (1. + (inlet.static_.props.Y - 1.) / 2. * pow(inlet.tip.M, 2));
        inlet.static_.props.H = inlet.total.props.H - std::pow(inlet.tip.C, 2.0) / 2.0;
        try {
            inlet.static_.set_props("HS", inlet.static_.props.H, inlet.static_.props.S);
        } catch (const std::runtime_error& e) {
            fmt::print(fg(fmt::color::red), "Runtime error: {}\n", e.what());
            return std::nullopt;
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
            return std::nullopt;
        }
        error = std::abs((inlet.static_.props.D - rhoOld) / rhoOld);
        rhoOld = inlet.static_.props.D;

        printIteration(iteration, inlet.static_.props.P, inlet.static_.props.T, inlet.static_.props.D, error,
                       "Density");

        iteration++;
    } while (iteration <= maxIterations && error > tolerance);

#ifdef DEBUG
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Solver execution time: " << std::setw(10) << std::fixed << std::setprecision(6)
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
#endif

    fmt::print("{:=<90}\n", "");

    inlet.rms.C = inlet.tip.C;
    inlet.hub.C = inlet.tip.C;
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

std::optional<double> Impeller::outletAungierSolver(int maxIterations, double tolerance, double rotorEfficiency) {
    int pressIteration = 1;
    int effIteration = 1;
    int velIteration = 1;
    double pressError;
    double effError;
    double velError;

    printBorder("Aungier", tolerance, maxIterations, OUTLET);

#ifdef DEBUG
    auto start = std::chrono::high_resolution_clock::now();
#endif
    estimateAxialLength();
    constexpr double PRguess = 2.0;
    double area = geom.area2;
    outlet.hub.U = op.omega * (geom.r2 * MM_M);
    outlet.rms.U = op.omega * (geom.r2 * MM_M);
    outlet.tip.U = op.omega * (geom.r2 * MM_M);
    do {
        // reset efficiency loop
        effIteration = 1;
        effError = 1.0;
        while (effIteration <= maxIterations && effError > tolerance) {
            double effOld = rotorEfficiency;
            outlet.total.props.P = PRguess * inlet.total.props.P;
            try {
                outlet.isentropic.set_props("PS", outlet.total.props.P, inlet.total.props.S);
            } catch (const std::runtime_error& e) {
                fmt::print(fg(fmt::color::red), "Runtime error: {}\n", e.what());
                return std::nullopt;
            } catch (const std::exception& e) {
                fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
                return std::nullopt;
            }
            outlet.total.props.H =
                inlet.total.props.H + (outlet.isentropic.props.H - inlet.total.props.H) / rotorEfficiency;

            try {
                outlet.total.set_props("HP", outlet.total.props.H, outlet.total.props.P);
            } catch (const std::runtime_error& e) {
                fmt::print(fg(fmt::color::red), "Runtime error: {}\n", e.what());
                return std::nullopt;
            } catch (const std::exception& e) {
                fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
                return std::nullopt;
            }

            dH0 = outlet.total.props.H - inlet.total.props.H;
            outlet.rms.C_theta = (1.0 / outlet.rms.U) * (dH0 + inlet.rms.U * inlet.rms.C_theta);

            // reset velocity loop
            outlet.rms.C_m = op.mfr / (outlet.total.props.D * area);
            if (!aungierVelocityLoop(maxIterations, tolerance)) {
                return std::nullopt;
            }

            // Get outlet velocities
            outlet.hub.C_m = outlet.rms.C_m;
            outlet.tip.C_m = outlet.rms.C_m;
            outlet.hub.C_theta = outlet.rms.C_theta;
            outlet.tip.C_theta = outlet.rms.C_theta;
            outlet.hub.C = outlet.rms.C;
            outlet.tip.C = outlet.rms.C;
            outlet.hub.M = outlet.rms.M;
            outlet.tip.M = outlet.rms.M;
            calculateOutletVelocities();

            ImpellerLosses losses = internalLosses();
            double H02real = outlet.isentropic.props.H + 5000.0;
            // fmt::println("outlet.isentropic.props.H: {}", outlet.isentropic.props.H);
            // fmt::println("H02real: {}", H02real);
            rotorEfficiency = (outlet.isentropic.props.H - inlet.total.props.H) / (H02real - inlet.total.props.H);
            effError = std::abs((rotorEfficiency - effOld) / effOld);

            printIteration(effIteration, outlet.static_.props.P, outlet.static_.props.T, rotorEfficiency, effError,
                           "Efficiency");

            effIteration++;
        }

        pressError = 1.0;
        pressIteration++;
    } while (pressIteration <= 1 && pressError > tolerance);

#ifdef DEBUG
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Solver execution time: " << std::setw(10) << std::fixed << std::setprecision(6)
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << " ms\n";
#endif

    return std::nullopt;
}

// This is the velocity while loop to converge on static density
bool Impeller::aungierVelocityLoop(int maxIterations, double tolerance) {
    int velIteration = 1;
    double velError = 1.0;
    double rhoOld = outlet.total.props.D;

    while (velIteration <= maxIterations && velError > tolerance) {
        outlet.rms.C = std::sqrt(std::pow(outlet.rms.C_m, 2) + std::pow(outlet.rms.C_theta, 2));
        outlet.static_.props.H = outlet.total.props.H - std::pow(outlet.rms.C, 2) / 2.0;

        try {
            outlet.static_.set_props("HS", outlet.static_.props.H, outlet.total.props.S);
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "Error: {}\n", e.what());
            return false;
        }

        outlet.rms.C_m = op.mfr / (outlet.static_.props.D * geom.area2);
        velError = std::abs((outlet.static_.props.D - rhoOld) / rhoOld);

        rhoOld = outlet.static_.props.D;
        outlet.rms.M = outlet.rms.C / outlet.static_.props.A;

        velIteration++;
    }

    if (velIteration < maxIterations) {
        fmt::print(fg(fmt::color::green), "Inner velocity loop converged with error {:.4e}\n", velError);
    }

    return true;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          General functions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Impeller::calculateInletVelocities() {
    std::array<Velocities, 3> velArr = {inlet.hub, inlet.rms, inlet.tip};
    std::array<double, 3> radius = {geom.r1h * MM_M, geom.r1rms * MM_M, geom.r1t * MM_M};
    std::array<double, 3> angles = {geom.beta1h, geom.beta1rms, geom.beta1t};

    for (int i = 0; i < velArr.size(); i++) {
        velArr[i].alpha = op.alpha;
        velArr[i].U = (radius[i]) * op.omega;
        if (velArr[i].alpha == 0) {
            velArr[i].C_m = velArr[i].C;
            velArr[i].C_theta = 0.0;
            velArr[i].W = sqrt(pow(velArr[i].C, 2.) + pow(velArr[i].U, 2.));
            velArr[i].W_m = velArr[i].C_m;
            velArr[i].W_theta = velArr[i].U;
            velArr[i].beta = -std::atan(velArr[i].U / velArr[i].C);
        } else {
            velArr[i].C_m = velArr[i].C * std::cos(velArr[i].alpha);
            velArr[i].C_theta = velArr[i].C * std::sin(velArr[i].alpha);
            velArr[i].W_m = velArr[i].C_m;
            velArr[i].W_theta = velArr[i].U = velArr[i].C_theta;
            velArr[i].W = sqrt(pow(velArr[i].W_m, 2.) + pow(velArr[i].W_theta, 2.));
            velArr[i].beta = -std::atan(velArr[i].W_theta / velArr[i].W_m);
        }
        velArr[i].M_rel = velArr[i].W / inlet.static_.props.A;
        velArr[i].inc = angles[i] - velArr[i].beta * (1.0 / DEG_RAD);
    }

    inlet.hub = velArr[0];
    inlet.rms = velArr[1];
    inlet.tip = velArr[2];
}

void Impeller::calculateOutletVelocities() {
    std::array<Velocities, 3> vel = {outlet.hub, outlet.rms, outlet.tip};

    for (auto& i : vel) {
        i.W_theta = i.U - i.C_theta;
        i.W_m = i.C_m;
        i.W = sqrt(pow(i.W_theta, 2.) + pow(i.W_m, 2.));
        i.beta = -atan(i.W_theta / i.W_m);
        i.M_rel = i.W / outlet.static_.props.A;
        i.C_slip = i.U - i.C_theta;
        i.alpha = atan(i.C_theta / i.C_m);
    }

    outlet.hub = vel[0];
    outlet.rms = vel[1];
    outlet.tip = vel[2];
}

void Impeller::calculateSlip(std::string slipModel) {
    if (slipModel == "Wiesner") {
        double sigma = 1.0 - sqrt(cos(geom.beta2 * DEG_RAD)) / pow(geom.ZFull, 0.7);
        outlet.hub.sigma = sigma;
        outlet.rms.sigma = sigma;
        outlet.tip.sigma = sigma;
    } else {
        fmt::print(fg(fmt::color::red), "Slip model not implemented. Slip model passed: {}\n", slipModel);
        exit(EXIT_FAILURE);
    }
}

void Impeller::estimateAxialLength() {
    // Aungier 2000 pg.113
    flowCoeff = op.mfr / (outlet.total.props.D * geom.area2 * outlet.tip.U);
    geom.deltaZ = (2 * geom.r2 * MM_M) * (0.014 + 0.023 * (geom.r2 / geom.r1h) + 1.58 * flowCoeff);
    fmt::println("(flowCoeff): {:.6f}", flowCoeff);
    fmt::println("(2 * geom.r2 * MM_M): {:.6f}", (2 * geom.r2 * MM_M));
    fmt::println("deltaZ: {:.6f}", geom.deltaZ);
}

ImpellerLosses Impeller::internalLosses() {
    ImpellerLosses losses = {};
    double cf = 0.0;
    double lambda = 0.0;
    double bb2 = 0.0;

    // Skin friction loss
    double passageA1 = geom.area1 * (std::abs(1.0 + std::cos(geom.beta1rms * DEG_RAD)) / 2.0) / float(geom.ZFull);
    double tb1 = 2.11;
    double tb = 2.11;
    double throatArea = passageA1 - (tb1 * MM_M * (geom.r1t * MM_M - geom.r1h * MM_M));
    double tipArea = 0.0;
    double tipPerimeter = 0.0;

    if (geom.ZSplit != 0) {
        tipArea = std::abs(geom.b2 * (2 * PI * (geom.r2 * MM_M) / geom.ZFull - (geom.ZFull + geom.ZSplit) * tb) *
                           std::cos(outlet.rms.beta * DEG_RAD));
        tipPerimeter = 4 * (geom.b2 + 2 * PI * (geom.r2 * MM_M) / (geom.ZFull + geom.ZSplit));
    } else {
        tipArea = std::abs(geom.b2 * (2 * PI * (geom.r2 * MM_M) / geom.ZFull - (geom.ZFull) * tb) *
                           std::cos(outlet.rms.beta * DEG_RAD));
        tipPerimeter = 2 * (geom.b2 + 2 * PI * (geom.r2 * MM_M) / (geom.ZFull));
    }

    double throatPerimeter =
        2 * ((geom.r1t * MM_M) - (geom.r1h * MM_M) + PI * (geom.r1t * MM_M) + (geom.r1h * MM_M) / geom.ZFull);
    losses.impeller = 2 * (throatArea / throatPerimeter + tipArea / tipPerimeter);

    double WBar = (inlet.rms.C * ((geom.r1t * MM_M) / (geom.r1rms * MM_M)) + outlet.rms.C + inlet.tip.W +
                   2.0 * inlet.hub.W + 3.0 * outlet.rms.W) /
                  8.0;
    double visc = (inlet.static_.props.V + outlet.static_.props.V) / 2.0;
    double ReImpeller = (inlet.static_.props.D + outlet.static_.props.V) / 2.0 * WBar * losses.impeller / visc;

    cf = skinFrictionCoefficient(ReImpeller, losses.impeller, 0.0000050, 1e-8);

    estimateAxialLength();
    double d1 = std::sqrt(std::pow((2 * geom.r1t * MM_M), 2) + std::pow((2 * geom.r1h * MM_M), 2));
    geom.lb =
        (geom.deltaZ - (geom.b2 * MM_M) / 2.0) + ((2 * geom.r2 * MM_M) - d1) / (2.0 * std::cos(geom.beta2 * DEG_RAD));
    fmt::print("((2 * geom.r2 * MM_M) - d1): {:.6f}\n", ((2 * geom.r2 * MM_M) - d1));
    // Lh = (r4 * (1 - r2rms * 2 / 0.3048) / (math.cos(beta4 / 180 * math.pi)));

    losses.skinFriction = (2.0 * cf * geom.lb / losses.impeller * std::pow(WBar, 2)) * 1000;
    fmt::print("losses.skinFriction: {:.4f}\n", losses.skinFriction);

    // Blade work input loss
    double AR = tipArea / throatArea;
    double phi2 = op.mfr / (outlet.static_.props.D * geom.area2 * outlet.rms.U);

    if (std::abs((outlet.total.props.P / outlet.static_.props.P) / outlet.total.props.P) > 1e-3) {
        // Tip blockage (Aungier, 2000)
        bb2 = losses.skinFriction * (inlet.total.props.P - inlet.static_.props.P) /
                  (outlet.total.props.P - outlet.static_.props.P) *
                  std::sqrt(std::abs(inlet.rms.W * losses.impeller / (outlet.rms.W * (geom.b2 * MM_M)))) +
              (0.3 + std::pow((geom.b2 * MM_M) / geom.lb, 2)) * std::pow(AR, 2) * outlet.static_.props.D *
                  (geom.b2 * MM_M) / (inlet.static_.props.D * geom.lb) +
              0.000372 / (2 * (geom.b2 * MM_M));
        fmt::print("bb2: {:.4f}\n", bb2);
        lambda = (1.0 / (1.0 - bb2));
        fmt::print("lambda: {:.4f}\n", lambda);
        losses.bladeWork = outlet.rms.sigma * (1.0 + lambda * phi2 * std::tan(geom.beta2 * DEG_RAD)) -
                           inlet.rms.U * inlet.rms.C_theta / std::pow(outlet.rms.U, 2);
        fmt::print("losses.bladeWork: {:.4f}\n", losses.bladeWork);
    } else {
        losses.bladeWork =
            (outlet.rms.U * outlet.rms.C_theta - inlet.rms.U * inlet.rms.C_theta) / std::pow(outlet.rms.U, 2);
        fmt::print("losses.bladeWork: {:.4f}\n", losses.bladeWork);
    }

    double Z = 0.0;
    if (geom.ZSplit > 0) {
        Z = geom.ZFull + (geom.lbSplitter / geom.lb) * geom.ZSplit;
    } else {
        Z = geom.ZFull;
    }
    double Df =
        1.0 - outlet.rms.W / inlet.tip.W +
        (0.75 * (outlet.rms.U * outlet.rms.C_theta - inlet.rms.U * inlet.rms.C_theta) / std::pow(outlet.rms.U, 2)) /
            ((inlet.rms.W / outlet.rms.W) * ((Z / PI) * (1 - geom.r1t / geom.r2) + (2 * geom.r1t / geom.r2)));

    losses.bladeLoading = 0.05 * std::pow(Df, 2) * std::pow(outlet.rms.U, 2);
    fmt::print("losses.bladeLoading: {:.4f}\n", losses.bladeLoading);

    double temp = skinFrictionLosses(2.e-5, 1e-8);
    fmt::print("temp: {}\n", temp);
    return losses;
}

double Impeller::skinFrictionCoefficient(double Re, double dH, double E, double tol) {
    double cf, cfl, cft, error;
    int iteration = 1;
    // Laminar flow if (actually <2000 but incase transitional)
    if (Re < 4000) {
        cf = 16.0 / Re;
    } else {
        double X = 1.0;
        double last = 0.0;
        // Turbulent flow over a smooth surfuace given by (Aungier, 2000)
        do {
            last = X;
            X = -2.0 * std::log10(2.51 * X / Re);
            error = std::abs((last - X) / X);
            iteration++;
        } while (error > tol && iteration <= 100);
        double cfts = 0.25 * std::pow(X, -2);
        // Turbulent flow over a fully rough surface given by (Aungier, 2000)
        X = -2.0 * std::log10(E / (3.71 * dH));
        double cftr = 0.25 * std::pow(X, -2);
        // Surface roughness becomes significant when Re_e < 60
        double Re_e = (Re - 2000.0) * E / dH;
        if (Re_e < 60) {
            cft = cfts;
        } else {
            cft = cfts + (cftr - cfts) * (1.0 - 60 / Re_e);
        }
    }

    // If transitional from laminar to turbulent, set cf to a weighted average of the two given by (Aungier, 2000)
    if (Re < 4000 && Re > 2000) {
        cfl = 16.0 / Re;
        cf = cfl + (cft - cfl) * (Re / 2000.0 - 1.0);
    } else {
        cf = cft;
    }

    return cf;
}

double colebrook(double x, double r, double Re) {
    return -2.0 * std::log10(r / 3.72 + 2.51 / Re / std::sqrt(x)) - 1.0 / std::sqrt(x);
}

double Impeller::skinFrictionLosses(double E, double tol) {
    double Cf, fx, dfx;
    double la = geom.r1h / geom.r1t;
    geom.Lh = (geom.r2 * MM_M) * (1.0 - (geom.r1rms * MM_M) * 2.0 / 0.3048) / (std::cos(geom.beta2 * DEG_RAD));
    geom.Dh = (2 * (geom.r2 * MM_M) *
               (1.0 / (geom.ZFull / PI / std::cos(geom.beta2 * DEG_RAD) + 2.0 * (geom.r2 * MM_M) / (geom.b2 * MM_M)) +
                (geom.r1t * MM_M) / (geom.r2 * MM_M) /
                    (2.0 / (1.0 - la) +
                     2.0 * (geom.ZFull) / PI / (1 + la) *
                         (std::sqrt(1 + (1 + std::pow(la, 2) / 2.0) * std::pow(std::tan(geom.beta1t * DEG_RAD), 2))))));

    double Re = (geom.Dh * (inlet.rms.W + outlet.rms.W) / 2 * (inlet.static_.props.D + outlet.static_.props.D) / 2 /
                 ((inlet.static_.props.V + outlet.static_.props.V) / 2));
    fmt::print("Re: {}\n", Re);
    double r = E / geom.Dh;
    int iteration = 1;
    double error = 1.0;
    double x = 0.01;
    double xnew = 1.0;
    double step = 1e-4;
    do {
        fx = colebrook(x, r, Re);
        dfx = colebrook(x + step, r, Re) - colebrook(x - step, r, Re) / (2 * step);
        if (std::abs(dfx) < 1e-12) {  // Prevent division by near-zero
            throw std::runtime_error("Derivative is too small, convergence not possible.");
        }
        xnew = x - fx / dfx;
        // fmt::print("iteration: {}\tfx: {}\tdfx: {}\txnew: {}\nerror: {}\n", iteration, fx, dfx, xnew, error);
        error = std::abs(xnew - x);
        x = xnew;
        iteration++;
    } while (error > 1e-4 && iteration <= 1000);
    Cf = x;
    fmt::print("Cf: {}\n", Cf);
    return 4 * Cf * geom.Lh * std::pow(((inlet.rms.W + outlet.rms.W) / 2.0), 2) / (2 * geom.Dh);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          Output functions
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Impeller::printBorder(std::string solver, double& tolerance, int& maxIterations, ImpellerStation station) {
    std::string border = fmt::format("{:=<90}\n", "");
    switch (station) {
    case INLET:
        fmt::print(
            "{}Entered impeller inlet solver loop!\nSolver: {}\nConvergence tolerance: {:.2e}\nMaximum iterations: "
            "{}\n{}",
            border, solver, tolerance, maxIterations, border);
        break;

    case OUTLET:
        fmt::print(
            "{}Entered impeller outlet solver loop!\nSolver: {}\nConvergence tolerance: {:.2e}\nMaximum iterations: "
            "{}\n{}",
            border, solver, tolerance, maxIterations, border);
        break;
    }
}

void Impeller::printIteration(int& iteration, double& P, double& T, double& M, double& error, std::string varName) {
    fmt::print("Iteration {}\n", iteration);
    fmt::print("Pressure: {:<10.0f} Temperature: {:<10.2f} {}: {:<10.6f} Error: {:<10.4e}\n", P, T, varName, M, error);
}

void Impeller::printOutputFile() {
    namespace fs = std::filesystem;
    const fs::path dirPath = "./../out";

    try {
        if (!fs::exists(dirPath)) {
            fs::create_directory(dirPath);
        }
        std::ofstream file(dirPath.string() + "/out.txt", std::ios::out);
        if (file.is_open()) {
            file << fmt::format("Impeller leading edge:\n");
            file << fmt::format("{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n", "P(Pa)", "T(K)", "H0(kJ/kg)", "s0(kJ/kg-K)",
                                "A(mm^2)", "Athroat(mm^2)");
            file << fmt::format("{:<15.0f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}\n\n", inlet.static_.props.P,
                                inlet.static_.props.T, inlet.total.props.H, inlet.total.props.S,
                                (geom.area1 * 1000000.0), (geom.throatArea * 1000000.0));
            file << fmt::format("{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<12}{:<12}{:<12}\n", "Station", "Dia(mm)",
                                "C(m/s)", "W(m/s)", "U(m/s)", "Mrel", "Beta(deg)", "Beta'(deg)", "Inc(deg)");

            const std::array<Velocities, 3> vel = {inlet.hub, inlet.rms, inlet.tip};
            const std::array<double, 3> diameters = {geom.r1h, geom.r1rms, geom.r1t};
            const std::array<double, 3> angles = {geom.beta1h, geom.beta1rms, geom.beta1t};
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
                                geom.r2, geom.b2, outlet.total.props.P, outlet.static_.props.P, outlet.total.props.T,
                                outlet.total.props.H, outlet.total.props.S, outlet.tip.U,
                                (outlet.tip.U / outlet.total.props.A));
            file << fmt::format("\n{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n", "M_rms", "W_rms(m/s)", "C_rms(m/s)",
                                "Alpha(deg)", "Beta(deg)", "W ratio");
            file << fmt::format("{:<15.3f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.2f}{:<15.3f}\n", outlet.rms.M, outlet.rms.W,
                                outlet.rms.C, outlet.rms.alpha * (1 / DEG_RAD), outlet.rms.beta * (1 / DEG_RAD),
                                (outlet.rms.W / inlet.rms.W));

            file << fmt::format("\nOverall Performance\n");
            file << fmt::format("{:<15}{:<15}{:<15}{:<15}{:<15}\n", "pr_tt", "dH/U^2", "dH0(kJ/kg)", "Re_b2", "Re_r2");
            file << fmt::format("{:<15.3f}{:<15.4f}{:<15.2f}{:<15.3e}{:<15.3e}\n", pr_tt, workCoeff, dH0, Re_b2, Re_r2);
            file.close();
        } else {
            fmt::print(fg(fmt::color::red), "Error creating output file.\n");
        }
    } catch (const fs::filesystem_error& e) {
        fmt::print(fg(fmt::color::red), "Error create output file: {}\n", e.what());
    }
}
