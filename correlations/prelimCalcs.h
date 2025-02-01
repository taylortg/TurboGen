#ifndef PRELIM_CALCS_H
#define PRELIM_CALCS_H

#include "../include/common.h"
#include "../include/thermo.h"

class DesignCorrelations {
   public:
    DesignCorrelations(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op);

   protected:
    ThermoProps thermo;
    Geometry geom;
    OperatingCondition op;
    double phi;
    double U2;
};

/*
 *  Casey and Robinson specific
 */
class CaseyRobinsonCorrelations : public DesignCorrelations {
   public:
    CaseyRobinsonCorrelations(const ThermoProps& thermo, const Geometry& geom, const OperatingCondition& op);

   private:
    double k1 = 27.0;
    double k2 = 5000.0;
    double k3 = 10.0;
    double k4 = 0.05;
    double k5 = 3.0;
    double eta_max = 0.86;
    double phi_max = 0.08;

    double Mu2;
    double eta_p;
    double deta_p;
    double P;
};

#endif  // PRELIM_CALCS_h