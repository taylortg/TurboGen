#ifndef COMMON_H
#define COMMON_H

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <optional>
#include <string>

#include "../include/thermo.h"

constexpr double PI = 3.14159265358979323846;
constexpr double DEG_RAD = PI / 180.0;
constexpr double MM_M = 1. / 1000.;

enum ThermoType { TOTAL, STATIC, ISENTROPIC, NUM_THERMO_TYPES };

enum VelType { HUB, RMS, TIP, NUM_VEL_TYPES };

enum class FileKey {
    SPEED,
    MASS_FLOW,
    INLET_FLOW_ANGLE,  // OperatingConditions
    R1H,
    R1T,
    R2,
    B2,
    BETA1H,
    BETA1T,
    BETA2,
    BLOCKAGE1,
    BLOCKAGE2,
    NUM_BLADES,
    AXIAL_LENGTH_RATIO,  // Geometry
    FLUID,
    PRESSURE,
    TEMPERATURE
};

constexpr std::array<std::pair<FileKey, const char*>, 17> fileKeyStrings = {
    {{FileKey::SPEED, "SPEED"},
     {FileKey::MASS_FLOW, "MASS_FLOW"},
     {FileKey::INLET_FLOW_ANGLE, "INLET_FLOW_ANGLE"},
     {FileKey::R1H, "R1H"},
     {FileKey::R1T, "R1T"},
     {FileKey::R2, "R2"},
     {FileKey::B2, "B2"},
     {FileKey::BETA1H, "BETA1H"},
     {FileKey::BETA1T, "BETA1T"},
     {FileKey::BETA2, "BETA2"},
     {FileKey::BLOCKAGE1, "BLOCKAGE1"},
     {FileKey::BLOCKAGE2, "BLOCKAGE2"},
     {FileKey::NUM_BLADES, "NUM_BLADES"},
     {FileKey::AXIAL_LENGTH_RATIO, "AXIAL_LENGTH_RATIO"},
     {FileKey::FLUID, "FLUID"},
     {FileKey::PRESSURE, "PRESSURE"},
     {FileKey::TEMPERATURE, "TEMPERATURE"}}};

struct Velocities {
    // absolute
    double C_theta;
    double C_m;
    double C;
    double C_slip;

    // relative
    double W_theta;
    double W_m;
    double W;

    // tip
    double U;

    // Mach numbers
    double M;
    double M_rel;

    // angles
    double alpha;  // absolute
    double beta;   // relative
    double sigma;  // slip
    double inc;    // incidence
};

struct OperatingCondition {
    double speed;  // rev/min
    double omega;  // rad/s
    double mfr;    // mass flow rate
    double alpha;  // intlet flow angle
};

struct Geometry {
    double r1h;
    double r1rms;
    double r1t;
    double r2;
    double b2;
    double beta1h;
    double beta1rms;
    double beta1t;
    double beta2;
    double blockage1;
    double blockage2;
    double Z;
    double axialLengthRatio;
    double axialLength;
    double area1;
    double area2;
    double throatArea;
};

std::optional<double> tryConvertingStrToDouble(const std::string& var);
OperatingCondition filterFileContent_op(const std::map<std::string, std::string>& fileContent);
Geometry filterFileContent_geom(const std::map<std::string, std::string>& fileContent);
ThermoProps filterFileContent_thermo(const std::map<std::string, std::string>& fileContent);
double getRMS(const double& val1, const double& val2);
std::string keyToString(FileKey key);

template <typename T>
bool setPropertyFromFileContent(const std::map<std::string, std::string>& fileContent, const std::string& key,
                                T& property) {
    auto it = fileContent.find(key);
    if (it != fileContent.end()) {
        auto valueOpt = tryConvertingStrToDouble(it->second);
        if (valueOpt.has_value()) {
            property = valueOpt.value();
            return true;
        }
    }
    std::cerr << "Error: Key '" << key << "' is missing or has an invalid value.\n";
    return false;
}

#endif  // COMMON_H
