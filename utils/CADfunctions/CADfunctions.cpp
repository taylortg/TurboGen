#include "CADfunctions.h"

#include <corecrt_math_defines.h>

#include <functional>


double Cfunc(double x) { return sin(x); }

int main() {
    std::function<double(double)> wrapper = Cfunc;

    double L = getArcLength(wrapper, 0.0, M_PI, 100, trapezoidalRule);
    std::cout << "Trapezoidal arc length: " << L << '\n';
    L = getArcLength(wrapper, 0.0, M_PI, 100, simpsonsRule);
    std::cout << "Simpsons arc length: " << L << '\n';
    L = getArcLength(wrapper, 0.0, M_PI, 10, gaussQuadrature);
    std::cout << "Gauss arc length: " << L << '\n';
}