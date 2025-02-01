#include "../include/thermo.h"

#include <iomanip>
#include <string>

#include "../externals/fmt/include/fmt/color.h"
#include "../externals/fmt/include/fmt/core.h"

void ThermoProps::set_props(std::string const& input_pair_, double const& val1, double const& val2) {
    if (input_pair_ == "PT") {
        props.P = val1;
        props.T = val2;
        try {
            if (!state) {
                fmt::print(fg(fmt::color::red), "CoolProp state is not initialized properly.\n");
                return;
            }
            state->update(CoolProp::PT_INPUTS, val1, val2);
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "CoolProp error: {}\n", e.what());
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
    } else if (input_pair_ == "HS") {
        props.H = val1;
        props.S = val2;
        try {
            if (!state) {
                fmt::print(fg(fmt::color::red), "CoolProp state is not initialized properly.\n");
                return;
            }
            state->update(CoolProp::HmassSmass_INPUTS, val1, val2);
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "CoolProp error: {}\n", e.what());
            return;
        }
        props.D = state->keyed_output(CoolProp::iDmass);
        props.P = state->keyed_output(CoolProp::iP);
        props.T = state->keyed_output(CoolProp::iT);
        props.V = state->keyed_output(CoolProp::iviscosity);
        props.A = state->keyed_output(CoolProp::ispeed_sound);
        props.CP = state->keyed_output(CoolProp::iCpmass);
        props.CV = state->keyed_output(CoolProp::iCvmass);
        props.Y = props.CP / props.CV;
    } else if (input_pair_ == "TH") {
        props.T = val2;
        props.H = val1;
        try {
            if (!state) {
                fmt::print(fg(fmt::color::red), "CoolProp state is not initialized properly.\n");
                return;
            }
            state->update(CoolProp::HmassT_INPUTS, val1, val2);
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "CoolProp error: {}\n", e.what());
            return;
        }
        props.D = state->keyed_output(CoolProp::iDmass);
        props.P = state->keyed_output(CoolProp::iP);
        props.S = state->keyed_output(CoolProp::iSmass);
        props.V = state->keyed_output(CoolProp::iviscosity);
        props.A = state->keyed_output(CoolProp::ispeed_sound);
        props.CP = state->keyed_output(CoolProp::iCpmass);
        props.CV = state->keyed_output(CoolProp::iCvmass);
        props.Y = props.CP / props.CV;
    } else if (input_pair_ == "PS") {
        props.P = val1;
        props.S = val2;
        try {
            if (!state) {
                fmt::print(fg(fmt::color::red), "CoolProp state is not initialized properly.\n");
                return;
            }
            state->update(CoolProp::PSmass_INPUTS, val1, val2);
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "CoolProp error: {}\n", e.what());
            return;
        }
        props.D = state->keyed_output(CoolProp::iDmass);
        props.T = state->keyed_output(CoolProp::iT);
        props.H = state->keyed_output(CoolProp::iHmass);
        props.V = state->keyed_output(CoolProp::iviscosity);
        props.A = state->keyed_output(CoolProp::ispeed_sound);
        props.CP = state->keyed_output(CoolProp::iCpmass);
        props.CV = state->keyed_output(CoolProp::iCvmass);
        props.Y = props.CP / props.CV;
    } else if (input_pair_ == "HP") {
        props.P = val2;
        props.H = val1;
        try {
            if (!state) {
                fmt::print(fg(fmt::color::red), "CoolProp state is not initialized properly.\n");
                return;
            }
            state->update(CoolProp::HmassP_INPUTS, val1, val2);
        } catch (const std::exception& e) {
            fmt::print(fg(fmt::color::red), "CoolProp error: {}\n", e.what());
            return;
        }
        props.D = state->keyed_output(CoolProp::iDmass);
        props.T = state->keyed_output(CoolProp::iT);
        props.S = state->keyed_output(CoolProp::iSmass);
        props.V = state->keyed_output(CoolProp::iviscosity);
        props.A = state->keyed_output(CoolProp::ispeed_sound);
        props.CP = state->keyed_output(CoolProp::iCpmass);
        props.CV = state->keyed_output(CoolProp::iCvmass);
        props.Y = props.CP / props.CV;
    }
}
