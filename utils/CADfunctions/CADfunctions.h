#ifndef CADFUNCTIONS_H
#define CADFUNCTIONS_H

#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
#include <type_traits>
#include <vector>

double dCdt(std::function<double(double)> Cfunc, double x, double step) {
    double result = (Cfunc(x + step) - Cfunc(x - step)) / (2 * step);
    return result;
}

struct GaussLegendreData {
    std::vector<double> w;  // weights
    std::vector<double> x;  // arguments

    // Constructor to initialize the data for the nth order
    GaussLegendreData(uint8_t n) {
        if (n == 1) {
            w = {2.0};
            x = {0.0};
        } else if (n == 2) {
            w = {1.0, 1.0};
            x = {-0.5773503, 0.5773503};
        } else if (n == 3) {
            w = {0.5555556, 0.8888889, 0.5555556};
            x = {-0.7745967, 0.0, 0.7745967};
        } else if (n == 4) {
            w = {0.347854845, 0.652145155, 0.652145155, 0.347854845};
            x = {-0.861136312, -0.339981044, 0.339981044, 0.861136312};
        } else if (n == 5 || n == 0) {
            w = {0.236926885, 0.478628670, 0.568888889, 0.478628670, 0.236926885};
            x = {-0.906179846, -0.538469310, 0.0, 0.538469310, 0.906179846};
        } else if (n == 10) {
            w = {0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2692667193099963, 0.2190863625159820,
                 0.2190863625159820, 0.1494513491505806, 0.1494513491505806, 0.0666713443086881, 0.0666713443086881};
            x = {-0.1488743389816312, 0.1488743389816312, -0.4333953941292472, 0.4333953941292472,
                 -0.6794095682990244, 0.6794095682990244, -0.8650633666889845, 0.8650633666889845,
                 -0.9739065285171717, 0.9739065285171717};
        } else {
            std::cerr << "Error: " << n << " order guass quadrature intergration not supported\n";
        }
    }
};

// Simpson's Rule Integration
double simpsonsRule(std::function<double(double)> func, double lower, double upper, int n) {
    double h = (upper - lower) / static_cast<double>(n);
    double sum = func(lower) + func(upper);

    for (int i = 1; i < n; i++) {
        double x = lower + i * h;
        if (i % 2 == 0) {
            sum += 4 * func(x);  // Even terms are multiplied by 4
        } else {
            sum += 2 * func(x);  // Odd terms are multiplied by 2
        }
    }

    return sum * h / 3.0;
}

// Trapezoidal Rule Integration
double trapezoidalRule(std::function<double(double)> func, double lower, double upper, int n) {
    double h = (upper - lower) / static_cast<double>(n);
    double sum = func(lower) + func(upper);

    for (int i = 1; i < n; i++) {
        double x = lower + i * h;
        sum += 2 * func(x);
    }

    return sum * h / 2.0;
}

// Gauss Quadrature (with a basic fixed set of points and weights)
double gaussQuadrature(std::function<double(double)> func, double lower, double upper, int n) {
    GaussLegendreData gauss(n);

    double bma = (upper - lower) / 2.0;
    double baa = (upper + lower) / 2.0;
    double sum = 0.0;

    for (int i = 0; i < n; ++i) {
        // Scales the interval in case the interval isn't [-1, 1]
        double xi = bma * gauss.x[i] + baa;
        sum += gauss.w[i] * func(xi);
    }

    return sum * bma;
}

double getArcLength(std::function<double(double)> Cfunc, double lower, double upper, int n,
                    std::function<double(std::function<double(double)>, double, double, int)> integrationMethod) {
    // Define the function to integrate: sqrt(1 + (dC/dx)^2)
    auto arcLengthFunction = [&](double x) {
        double fxPrime = dCdt(Cfunc, x, 1e-3);  // Derivative
        // std::cout << "x: " << x << ", dC/dx: " << fxPrime << "\n";
        return std::sqrt(1 + fxPrime * fxPrime);  // Arc length differential
    };

    return integrationMethod(arcLengthFunction, lower, upper, n);
}

#endif  // CADFUNCTIONS_H