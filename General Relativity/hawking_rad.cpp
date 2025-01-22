#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <random>
#include <vector>

class BlackHole {
private:
    double mass;
    double temperature;
    const double BOLTZMANN = 1.380649e-23;
    const double PLANCK = 6.62607015e-34;
    const double C = 2.998e8;
    const double G = 6.674e-11;

public:
    BlackHole(double initial_mass) : mass(initial_mass) {
        calculateTemperature();
    }

    void calculateTemperature() {
        temperature = (PLANCK * pow(C, 3)) / (8 * M_PI * G * mass * BOLTZMANN);
    }

    double emitRadiation(double time_step) {
        double surface_area = 4 * M_PI * pow(2 * G * mass / pow(C, 2), 2);
        double stefan_boltzmann = 5.67e-8;
        double power = stefan_boltzmann * surface_area * pow(temperature, 4);
        
        double mass_loss = (power * time_step) / pow(C, 2);
        mass -= mass_loss;
        
        calculateTemperature();
        return mass_loss;
    }

    double getMass() const { return mass; }
    double getTemperature() const { return temperature; }
};

int main() {
    BlackHole bh(10);
    
    double time_step = 1e10;
    int steps = 10;

    std::cout << "Simulating Hawking radiation for a " << bh.getMass() 
              << " solar mass black hole\n\n";

    for (int i = 0; i < steps; i++) {
        double mass_loss = bh.emitRadiation(time_step);
        std::cout << "Time step " << i + 1 << ":\n";
        std::cout << "Mass: " << bh.getMass() << " solar masses\n";
        std::cout << "Temperature: " << bh.getTemperature() << " K\n";
        std::cout << "Mass lost: " << mass_loss << " kg\n\n";
    }

    return 0;
}
