#include "../include/Aungier.h"

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

void Aungier::runCalculations() { inletCalcs(); }

void Aungier::inletCalcs() {
    int maxIterations = 100;
    int iteration = 0;
    double error = 1.0;
    double r1h = geom.r1h * MM_M;
    double r1t = geom.r1t * MM_M;
    double b1 = r1t - r1h;
    double A1 = PI * (std::pow(r1t, 2) - std::pow(r1h, 2)) - b1 * geom.ZFull * (2.11 * MM_M);
    fmt::println("A1: {}", A1);
    double r1 = getRMS(r1t, r1h);
    double rho1_guess = inlet.total.props.D;
    double tol = 10e-6;
    inlet.static_ = inlet.total;

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