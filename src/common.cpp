#include "../include/common.h"

#include <cmath>
#include <iostream>
#include <map>
#include <optional>
#include <string>

std::optional<double> tryConvertingStrToDouble(const std::string& var) {
    try {
        return std::stod(var);
    } catch (const std::exception& e) {
        std::cerr << "Error converting \"" << var << "\" to double: " << e.what() << " (Input length: " << var.size()
                  << ")\n";
        return std::nullopt;
    }
}

std::string keyToString(FileKey key) {
    for (const auto& pair : fileKeyStrings) {
        if (pair.first == key) {
            return pair.second;
        }
    }
    throw std::invalid_argument("Unknown FileKey");
}

OperatingCondition filterFileContent_op(const std::map<std::string, std::string>& fileContent) {
    OperatingCondition op{};
    auto speedKey = keyToString(FileKey::SPEED);
    auto mfrKey = keyToString(FileKey::MASS_FLOW);
    auto alphaKey = keyToString(FileKey::INLET_FLOW_ANGLE);

    setPropertyFromFileContent(fileContent, speedKey, op.speed);
    setPropertyFromFileContent(fileContent, mfrKey, op.mfr);
    setPropertyFromFileContent(fileContent, alphaKey, op.alpha);

    if (op.speed > 0) {
        op.omega = op.speed * (PI / 30.0);
    }

    return op;
}

Geometry filterFileContent_geom(const std::map<std::string, std::string>& fileContent) {
    Geometry geom{};
    auto key = keyToString(FileKey::R1H);
    setPropertyFromFileContent(fileContent, key, geom.r1h);
    key = keyToString(FileKey::R1T);
    setPropertyFromFileContent(fileContent, key, geom.r1t);
    key = keyToString(FileKey::R2);
    setPropertyFromFileContent(fileContent, key, geom.r2);
    key = keyToString(FileKey::B2);
    setPropertyFromFileContent(fileContent, key, geom.b2);

    if (geom.r1h > 0 && geom.r1t > 0) {
        geom.r1rms = getRMS(geom.r1h, geom.r1t);
    }

    key = keyToString(FileKey::BETA1H);
    setPropertyFromFileContent(fileContent, key, geom.beta1h);
    key = keyToString(FileKey::BETA1T);
    setPropertyFromFileContent(fileContent, key, geom.beta1t);
    geom.beta1rms = -getRMS(geom.beta1h, geom.beta1t);

    key = keyToString(FileKey::BETA2);
    setPropertyFromFileContent(fileContent, key, geom.beta2);
    key = keyToString(FileKey::BLOCKAGE1);
    setPropertyFromFileContent(fileContent, key, geom.blockage1);
    key = keyToString(FileKey::BLOCKAGE2);
    setPropertyFromFileContent(fileContent, key, geom.blockage2);
    key = keyToString(FileKey::NUM_FULL_BLADES);
    setPropertyFromFileContent(fileContent, key, geom.ZFull);
    key = keyToString(FileKey::NUM_SPLITTER_BLADES);
    setPropertyFromFileContent(fileContent, key, geom.ZSplit);
    key = keyToString(FileKey::AXIAL_LENGTH_RATIO);
    setPropertyFromFileContent(fileContent, key, geom.axialLengthRatio);

    geom.area1 = (1. - geom.blockage1) * PI * (pow((geom.r1t * MM_M), 2) - pow((geom.r1h * MM_M), 2));
    geom.area2 = (2 * PI * (geom.r2 * MM_M) * (geom.b2 * MM_M)) * (1.0 - geom.blockage2);
    geom.throatArea = geom.area1 * cos(std::abs(geom.beta1rms) * DEG_RAD);

    return geom;
}

ThermoProps filterFileContent_thermo(const std::map<std::string, std::string>& fileContent) {
    try {
        auto key = keyToString(FileKey::FLUID);
        ThermoProps thermo{fileContent.at(key)};

        key = keyToString(FileKey::PRESSURE);
        setPropertyFromFileContent(fileContent, key, thermo.props.P);
        key = keyToString(FileKey::TEMPERATURE);
        setPropertyFromFileContent(fileContent, key, thermo.props.T);

        thermo.set_props("PT", thermo.props.P, thermo.props.T);

        return thermo;
    } catch (std::exception& e) {
        std::cout << "Error setting inlet thermodynamic conditions: " << e.what() << "\n";
        exit(1);
    }
}

double getRMS(const double& val1, const double& val2) { return sqrt((val1 * val1 + val2 * val2) / 2.0); }
