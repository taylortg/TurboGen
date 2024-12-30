#include "../include/thermo.h"

#include <iomanip>
#include <string>


void ThermoProps::set_props(std::string const& input_pair_, double const& val1, double const& val2) {
    if (input_pair_ == "PT") {
        props.P = val1;
        props.T = val2;
        try {
            if (!state) {
                std::cerr << "CoolProp state is not initialized properly." << std::endl;
                return;
            }
            state->update(CoolProp::PT_INPUTS, val1, val2);
        } catch (const std::exception& e) {
            std::cerr << "CoolProp error: " << e.what() << std::endl;
            return;
        }
        props.D = state->keyed_output(CoolProp::iDmass);
        props.H = state->keyed_output(CoolProp::iHmass);
        props.S = state->keyed_output(CoolProp::iSmass);
        props.V = state->keyed_output(CoolProp::iviscosity);
        props.A = state->keyed_output(CoolProp::ispeed_sound);
        props.CP = state->keyed_output(CoolProp::iCpmass);
        props.CV = state->keyed_output(CoolProp::iCvmass);
        props.Y = props.CP / props.CV;
    }
}
