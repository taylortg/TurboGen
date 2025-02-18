#ifndef AUNGIER_H
#define AUNGIER_H

#include "common.h"
#include "impeller.h"

class Aungier {
   public:
    ImpellerState inlet;
    ImpellerState outlet;
    ImpellerState throat;
    OperatingCondition op;
    Geometry geom;
    double pr_tt, isenEff, flowCoeff, workCoeff, Re_b2, Re_r2, dH0;
    bool optimizationFlag;

    Aungier(const ImpellerState& inlet, const ImpellerState& outlet, const OperatingCondition& op, const Geometry& geom,
            const bool optimizationFlag);
    Aungier(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op, const bool optimizationFlag);
    Aungier(const Impeller& impeller, const bool optimizationFlag);

    void runCalculations();
    void inletCalcs();
    void inletInducerOptimization();
    void calculateInletVelocities();
};
#endif  // AUNGIER_H