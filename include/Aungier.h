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

    Aungier(const ImpellerState& inlet, const ImpellerState& outlet, const OperatingCondition& op,
            const Geometry& geom);
    Aungier(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op);
    Aungier(const Impeller& impeller);

    void runCalculations();
    void inletCalcs();
    void calculateInletVelocities();
};
#endif  // AUNGIER_H