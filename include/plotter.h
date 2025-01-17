#ifndef PLOTTER_H
#define PLOTTER_H

#include <string>

#include "impeller.h"

std::string formatHeaderRow(const std::string& variable);
std::string formatSubHeaderRow(const std::string& variable);
std::string formatVariableRow();
std::string formatTableRow(const std::string& variable, const std::string& units, float value, int precision);
void plotVelocityTriangle(const Impeller& impeller, bool htmlFlag);

#endif  // PLOTTER_H