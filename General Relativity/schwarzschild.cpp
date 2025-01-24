#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>

class SchwarzschildMetric {
private:
    const double G = 6.674e-11;
    const double c = 2.998e8;
    double M;
    double rs;

public:
    SchwarzschildMetric(double mass) : M(mass) {
        rs = 2 * G * M / (c * c);
    }

    double g_tt(double r) const {
        return -(1.0 - rs/r);
    }

    double g_rr(double r) const {
        return 1.0 / (1.0 - rs/r);
    }

    double g_theta_theta(double r) const {
        return r * r;
    }

    double g_phi_phi(double r, double theta) const {
        return r * r * sin(theta) * sin(theta);
    }

    bool isInsideHorizon(double r) const {
        return r <= rs;
    }

    double getSchwarzschildRadius() const {
        return rs;
    }
};

int main() {
    double solarMass = 1.989e30;
    SchwarzschildMetric blackHole(solarMass);

    double r = 3 * blackHole.getSchwarzschildRadius();
    double theta = M_PI / 2;

    std::cout << "Metric components at r = " << r << " meters:\n";
    std::cout << "g_tt = " << blackHole.g_tt(r) << "\n";
    std::cout << "g_rr = " << blackHole.g_rr(r) << "\n";
    std::cout << "g_θθ = " << blackHole.g_theta_theta(r) << "\n";
    std::cout << "g_φφ = " << blackHole.g_phi_phi(r, theta) << "\n";

    return 0;
}
