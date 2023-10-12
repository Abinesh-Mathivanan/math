#include <iostream>
#include <cmath>

const double G = 6.67430e-11;  // Gravitational constant (m^3 kg^-1 s^-2)
const double M = 1.989e30;     // Mass of the central object (e.g., the sun) in kg

// Define a structure to represent a point in 3D space
struct Point {
    double x, y, z;
};

// Define a structure to represent a four-vector (t, x, y, z)
struct FourVector {
    double t, x, y, z;
};

// Function to calculate the Schwarzschild metric
FourVector SchwarzschildMetric(const Point& p) {
    double r = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    double rs = 2.0 * G * M;  // Schwarzschild radius
    double g_tt = 1.0 - rs / r;
    
    FourVector metric;
    metric.t = 1.0 / sqrt(g_tt);
    metric.x = p.x / r;
    metric.y = p.y / r;
    metric.z = p.z / r;
    
    return metric;
}

int main() {
    Point initialPosition;
    double dt;
    int iterations;

    // Get user input for initial position, time step, and number of iterations
    std::cout << "Enter initial x, y, z (in meters): ";
    std::cin >> initialPosition.x >> initialPosition.y >> initialPosition.z;

    std::cout << "Enter time step (in seconds): ";
    std::cin >> dt;

    std::cout << "Enter the number of iterations: ";
    std::cin >> iterations;

    double t = 0.0;  // Initial time

    // Perform numerical integration to follow the geodesic
    for (int i = 0; i < iterations; ++i) {
        // Calculate the Schwarzschild metric at the current position
        FourVector initialFourVector = SchwarzschildMetric(initialPosition);

        // Update the position based on the Schwarzschild metric
        initialPosition.x += initialFourVector.x * dt;
        initialPosition.y += initialFourVector.y * dt;
        initialPosition.z += initialFourVector.z * dt;

        // Update the time
        t += initialFourVector.t * dt;

        // Print the updated position and time
        std::cout << "t: " << t << " x: " << initialPosition.x << " y: " << initialPosition.y << " z: " << initialPosition.z << std::endl;
    }

    return 0;
}