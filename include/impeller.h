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
    Geometry geom;

    ImpellerState(const std::string& fluid_name = "AIR");
    ImpellerState(const ThermoProps& thermo, const Geometry& geom);
    ImpellerState(const ImpellerState& other);
    ImpellerState& operator=(const ImpellerState& other);
};

class Impeller {
   public:
    ImpellerState inlet;
    ImpellerState outlet;
    OperatingCondition op;

    Impeller(ImpellerState& inlet, ImpellerState& outlet, const OperatingCondition& op);
    Impeller(ThermoProps thermo, Geometry geom, const OperatingCondition& op);

    void calculateInletCondition(std::string solverType);
    std::optional<double> inletJapikseSolver(double Mguess, int maxIterations, double tolerance);
    void calculateInletVelocities(ImpellerState* inlet, const OperatingCondition& op);
    void printInletVelocities();
};

#endif  // IMPELLER_H
