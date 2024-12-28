#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <filesystem>

#include "../include/impeller.h"

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          ImpellerState class
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ImpellerState::ImpellerState(const std::string& fluid_name)
    : total(fluid_name),
      static_(fluid_name),
      isentropic(fluid_name),
      hub(),
      rms(),
      tip(),
      geom() {}


ImpellerState::ImpellerState(ThermoProps thermo, const Geometry& geom)
    : total(thermo),
      static_(thermo),
      isentropic(thermo),
      hub(),
      rms(),
      tip(),
      geom(geom) {}


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
    : inlet(inlet),
      outlet(outlet),
      op(op) {}


Impeller::Impeller(ThermoProps thermo, Geometry geom, const OperatingCondition& op)
    : inlet(thermo, geom),
      outlet(thermo, geom),
      op(op) {}


void Impeller::calculateInletCondition(std::string solverType) {
    // To-do: catch if user wants to continue when solution does not converge
    // use if (inletSolverResult) {}
    if (solverType == "Japikse") {
        auto inletSolverResult = inletJapikseSolver(0.3, 1000, 1e-8);
    } else {
        std::cerr << "Solver type not implemented. Solver type pass:\n" << solverType << "\n";
        exit(EXIT_FAILURE);
    }
}


std::optional<double> Impeller::inletJapikseSolver(double Mguess = 0.3, int maxIterations = 1000, double tolerance = 1e-6) {
    int iteration = 1;
    double error = 1.0;
    std::ostringstream oss;
    oss << "===============================================================================================\n"
        << "Entered impeller inlet solver loop!\n"
        << "Solver: Japikse\n"
        << "Convergence tolerance: " << tolerance << "\n"
        << "Maximum iterations: " << maxIterations << "\n"
        << "===============================================================================================\n";
    std::cout << oss.str();

    // This is only here to indicate if there are issues with the solver. The execution time of the while loop
    // should be incredibly small. If you see it is taking ms to run, there might be a problem.
    #ifdef DEBUG
    auto start = std::chrono::high_resolution_clock::now();
	#endif

    do {
        inlet.static_.props.T = inlet.total.props.T / (1. + (inlet.static_.props.Y - 1.)/2. * pow(Mguess, 2));
        double k = inlet.static_.props.Y / (inlet.static_.props.Y - 1.);
        inlet.static_.props.P = inlet.total.props.P / (pow((inlet.total.props.T/inlet.static_.props.T), k));
        inlet.static_.set_props("PT", inlet.static_.props.P, inlet.static_.props.T);
        double Rgas = inlet.static_.props.CP - inlet.static_.props.CV;
        inlet.tip.C_m = (op.mfr * Rgas * inlet.static_.props.T) / (inlet.static_.props.P * inlet.geom.area1);
        inlet.tip.M = inlet.tip.C_m / inlet.static_.props.A;
        error = std::abs((Mguess - inlet.tip.M) / Mguess);
        Mguess = inlet.tip.M;

        std::cout << "Iteration: "<< iteration << "\n";
        std::cout << std::setw(10) << std::fixed << std::setprecision(0) << "Pressure: " <<  inlet.static_.props.P
                  << std::setw(10) << std::fixed << std::setprecision(2) << "  Temperature: " <<  inlet.static_.props.T
                  << std::setw(10) << std::fixed << std::setprecision(6) << "  Mach: " <<  inlet.tip.M
                  << std::setw(10) << std::fixed << std::setprecision(8) << "  Error: " <<  error << "\n";

        iteration++;
    } while (iteration <= maxIterations && error > tolerance);

    #ifdef DEBUG
    auto end = std::chrono::high_resolution_clock::now();
	std::cout << "Solver execution time: " << std::setw(10) << std::fixed << std::setprecision(6)
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms\n";
	#endif

    std::cout << "===============================================================================================\n";
    inlet.rms.C_m = inlet.tip.C_m;
    inlet.hub.C_m = inlet.tip.C_m;
    inlet.rms.M = inlet.tip.M;
    inlet.hub.M = inlet.tip.M;

    calculateInletVelocities(&inlet, op);

    if (error <= tolerance) {
        std::cout << "Solution converged!\n";
        return inlet.static_.props.P;
    } else {
        std::cout << "Solution did not converge, results could be inaccurate!\n";
        return std::nullopt;
    }
}


void Impeller::calculateInletVelocities(ImpellerState* inlet, const OperatingCondition& op) {
    std::array<Velocities, 3> velArr = {inlet->hub, inlet->rms, inlet->tip};
    std::array<double, 3> radius = {inlet->geom.r1h*MM_M, inlet->geom.r1rms*MM_M, inlet->geom.r1t*MM_M};
    std::array<double, 3> angles = {inlet->geom.beta1h, inlet->geom.beta1rms, inlet->geom.beta1t};

    for(int i = 0; i < velArr.size(); i++) {
        velArr[i].alpha = op.alpha;
        velArr[i].U = (radius[i]) * op.omega;
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
        velArr[i].M_rel = velArr[i].W / inlet->static_.props.A;
        velArr[i].inc = angles[i] - velArr[i].beta * (1.0/DEG_RAD);
    }

    inlet->hub = velArr[0];
    inlet->rms = velArr[1];
    inlet->tip = velArr[2];
    printInletVelocities();
}


void Impeller::printInletVelocities() {
    namespace fs = std::filesystem;
    fs::path dirPath = "./../out";

    try {
        if (!fs::exists(dirPath)) {
            fs::create_directory(dirPath);
        }
        std::ofstream file(dirPath.string() + "/out.log", std::ios::out);
        if (file.is_open()) {
            std::ostringstream oss;
            oss << "Impeller leading edge:\n";
            oss << std::left << std::setw(15) << "Pressure(Pa)"
                << std::left << std::setw(15) << "Temperature(K)"
                << std::left << std::setw(15) << "H0(kJ/kg)"
                << std::left << std::setw(15) << "s0(kJ/kg-K)"
                << std::left << std::setw(15) << "A(mm^2)"
                << std::left << std::setw(15) << "Athroat(mm^2)"
                << "\n";
            oss << std::left << std::setw(15) << std::fixed << std::setprecision(0) << inlet.static_.props.P
                << std::left << std::setw(15) << std::fixed << std::setprecision(2) << inlet.static_.props.T
                << std::left << std::setw(15) << std::fixed << std::setprecision(2) << inlet.total.props.H
                << std::left << std::setw(15) << std::fixed << std::setprecision(2) << inlet.total.props.S
                << std::left << std::setw(15) << std::fixed << std::setprecision(2) << inlet.geom.area1 * 1000000.0
                << std::left << std::setw(15) << std::fixed << std::setprecision(2) << inlet.geom.throatArea * 1000000.0
                << "\n\n";
            oss << std::left << std::setw(10) << "Station"
                << std::left << std::setw(10) << "Dia(mm)"
                << std::left << std::setw(10) << "C(m/s)"
                << std::left << std::setw(10) << "W(m/s)"
                << std::left << std::setw(10) << "U(m/s)"
                << std::left << std::setw(10) << "Mrel"
                << std::left << std::setw(12) << "Beta(deg)"
                << std::left << std::setw(12) << "Beta'(deg)"
                << std::left << std::setw(12) << "Inc(deg)" << "\n";
            const std::array<Velocities, 3> vel = {inlet.hub, inlet.rms, inlet.tip};
            const std::array<double, 3> diameters = {inlet.geom.r1h, inlet.geom.r1rms, inlet.geom.r1t};
            const std::array<double, 3> angles = {inlet.geom.beta1h, inlet.geom.beta1rms, inlet.geom.beta1t};
            const std::array<std::string, 3> names = {"Hub:", "RMS:", "Shroud:"};
            for (size_t i = 0; i < vel.size(); i++) {
                oss << std::left << std::setw(10) << names[i]
                << std::left << std::setw(10) << std::fixed << std::setprecision(1) << diameters[i]
                << std::left << std::setw(10) << std::fixed << std::setprecision(2) << vel[i].C
                << std::left << std::setw(10) << std::fixed << std::setprecision(2) << vel[i].W
                << std::left << std::setw(10) << std::fixed << std::setprecision(2) << vel[i].U
                << std::left << std::setw(10) << std::fixed << std::setprecision(4) << vel[i].M_rel
                << std::left << std::setw(12) << std::fixed << std::setprecision(2) << angles[i]
                << std::left << std::setw(12) << std::fixed << std::setprecision(2) << vel[i].beta * (1/DEG_RAD)
                << std::left << std::setw(12) << std::fixed << std::setprecision(2) << vel[i].inc
                << "\n";
            }
            file << oss.str();
            file.close();
        } else {
            std::cerr << "Error creating log file.\n";
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Error create log file: " << e.what() << "\n";
    }
}
