#ifndef IMPELLER_H
#define IMPELLER_H

#include <optional>

#include "common.h"
#include "thermo.h"

struct ImpellerLosses {
    double skinFriction = 0.0;
    double bladeLoading = 0.0;
    double clearance = 0.0;
    double incidence = 0.0;
    double discFriction = 0.0;
    double recirculation = 0.0;
};

class ImpellerState {
   public:
    ThermoProps total;
    ThermoProps static_;
    ThermoProps isentropic;
    Velocities hub;
    Velocities rms;
    Velocities tip;

    ImpellerState(const std::string& fluid_name = "AIR");
    ImpellerState(const ThermoProps& thermo);
    ImpellerState(const ImpellerState& other);
    ImpellerState& operator=(const ImpellerState& other);
};

class Impeller {
   public:
    ImpellerState inlet;
    ImpellerState outlet;
    OperatingCondition op;
    Geometry geom;
    double pr_tt, isenEff, flowCoeff, workCoeff, Re_b2, Re_r2, dH0;

    Impeller(ImpellerState& inlet, ImpellerState& outlet, const OperatingCondition& op, const Geometry& geom);
    Impeller(ThermoProps thermo, const Geometry& geom, const OperatingCondition& op);

    // Driver functions
    void calculateInletCondition(std::string solverType);
    void calculateOutletCondition(std::string solverType, std::string slipModel);

    // Japikse function
    std::optional<double> inletJapikseSolver(double Mguess, int maxIterations, double tolerance);
    std::optional<double> outletJapikseSolver(int maxIterations, double tolerance, double rotorEfficiency);
    void calculateOutletVelocities();

    // CCMD functions
    std::optional<double> inletAungierSolver(int maxIterations, double tolerance);

    // General/Loss model functions
    void calculateInletVelocities();
    void calculateSlip(std::string slipModel);
    void estimateAxialLength();

    // Output functions
    void printBorder(std::string solver, double& tolerance, int& maxIterations);
    void printIteration(int& iteration, double& P, double& T, double& M, double& error, std::string varName);
    void printOutputFile();
};

#endif  // IMPELLER_H
