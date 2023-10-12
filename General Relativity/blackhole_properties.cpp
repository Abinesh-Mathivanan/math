#include <iostream>
#include <cmath>

const double h_bar = 1.0545718e-34;  // Reduced Planck constant (J·s)
const double c = 299792458;          // Speed of light (m/s)
const double G = 6.67430e-11;       // Gravitational constant (m^3 kg^-1 s^-2)
const double k = 1.38064852e-23;    // Boltzmann constant (J/K)
const double pi = 3.14159265359;    // π

double calculateHawkingTemperature(double mass) {
    return (h_bar * c * c * c) / (8 * pi * G * mass * k);
}

double calculateMass(double area) {
    double r = sqrt(area / (4 * pi));
    return c * c * r / (2 * G);
}

double calculateEntropy(double area) {
    return area / (4 * h_bar * G);
}

int main() {
    double mass;
    double area;

    std::cout << "Enter the mass of the black hole (in kg): ";
    std::cin >> mass;
    if (mass <= 0) {
        std::cerr << "Mass must be positive." << std::endl;
        return 1;
    }

    std::cout << "Enter the area of the black hole's event horizon (in square meters): ";
    std::cin >> area;
    if (area <= 0) {
        std::cerr << "Area must be positive." << std::endl;
        return 1;
    }

    double hawkingTemperature = calculateHawkingTemperature(mass);
    double entropy = calculateEntropy(area);

    std::cout << "The Hawking temperature of the black hole is: " << hawkingTemperature << " K" << std::endl;
    std::cout << "The mass of the black hole is: " << mass << " kg" << std::endl;
    std::cout << "The entropy of the black hole is: " << entropy << " J/K" << std::endl;

    return 0;
}
