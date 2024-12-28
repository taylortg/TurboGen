#ifndef THERMO_H
#define THERMO_H

#include <limits>
#include <ostream>

#include "../externals/CoolProp/include/CoolProp.h"
#include "../externals/CoolProp/include/AbstractState.h"
#include "../externals/CoolProp/include/crossplatform_shared_ptr.h"


struct Props {
    double P = std::numeric_limits<double>::quiet_NaN();    // Pressure
    double T = std::numeric_limits<double>::quiet_NaN();    // Temperature
    double D = std::numeric_limits<double>::quiet_NaN();    // Density
    double H = std::numeric_limits<double>::quiet_NaN();    // Enthalpy
    double S = std::numeric_limits<double>::quiet_NaN();    // Entropy
    double V = std::numeric_limits<double>::quiet_NaN();    // Viscosity
    double A = std::numeric_limits<double>::quiet_NaN();    // Speed of sound
    double CP = std::numeric_limits<double>::quiet_NaN();    // Specific heat ratio
    double CV = std::numeric_limits<double>::quiet_NaN();    // Specific heat ratio
    double Y = std::numeric_limits<double>::quiet_NaN();    // Specific heat ratio
};


class ThermoProps {
public:
    Props props;
    std::shared_ptr<CoolProp::AbstractState> state;

    ThermoProps(const std::string& fluid_name = "AIR")
        : state(CoolProp::AbstractState::factory("HEOS", fluid_name)) {}

    ThermoProps(const ThermoProps& other)
        : props(other.props),
        state(other.state ?
            CoolProp::AbstractState::factory(
                "HEOS", other.state->fluid_names()) : nullptr) {}

    void set_props(std::string const& input_pair_, double const& val1, double const& val2);
};

#endif //THERMO_H
