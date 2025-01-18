#include <cmath>
#include <complex>
#include <vector>

class KerrBlackHole {
private:
    double mass;
    double angular_momentum;
    double rs;    // Schwarzschild radius
    double a;     // specific angular momentum

public:
    KerrBlackHole(double M, double J) : mass(M), angular_momentum(J) {
        rs = 2.0 * mass;
        a = angular_momentum / mass;
    }

    struct MetricTensor {
        double g00, g11, g22, g33;
        double g03, g30;
    };

    MetricTensor calculateMetric(double r, double theta) {
        MetricTensor metric;
        double rho2 = r * r + a * a * cos(theta) * cos(theta);
        double delta = r * r - rs * r + a * a;
        double sigma = pow(r * r + a * a, 2) - a * a * delta * sin(theta) * sin(theta);

        metric.g00 = -(1.0 - rs * r / rho2);
        metric.g11 = rho2 / delta;
        metric.g22 = rho2;
        metric.g33 = (r * r + a * a + rs * r * a * a * sin(theta) * sin(theta) / rho2) 
                     * sin(theta) * sin(theta);
        metric.g03 = metric.g30 = -rs * r * a * sin(theta) * sin(theta) / rho2;

        return metric;
    }

    std::pair<double, double> findHorizons() {
        double rplus = mass + sqrt(mass * mass - a * a);
        double rminus = mass - sqrt(mass * mass - a * a);
        return std::make_pair(rplus, rminus);
    }

    bool isInErgoregion(double r, double theta) {
        double rho2 = r * r + a * a * cos(theta) * cos(theta);
        return r < rs && r > 0 && rho2 < rs * r;
    }
};