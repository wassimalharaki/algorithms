#include <iostream>
#include <functional>
#include <random>
#include <cmath>

// Monte Carlo Approximation Template
template <typename T>
T monteCarloApproximation(int numSamples, const std::function<bool(T, T)>& insideFunction, T rangeMin, T rangeMax) {
    std::random_device rd;  // Random device for seeding
    std::mt19937 gen(rd()); // Mersenne Twister RNG
    std::uniform_real_distribution<T> dist(rangeMin, rangeMax); // Uniform distribution

    int insideCount = 0;

    for (int i = 0; i < numSamples; ++i) {
        T x = dist(gen); // Generate random x
        T y = dist(gen); // Generate random y

        // Check if the point is inside the desired region
        if (insideFunction(x, y)) {
            ++insideCount;
        }
    }

    // Calculate the area of the sampling space
    T area = (rangeMax - rangeMin) * (rangeMax - rangeMin);

    // Approximation: Fraction of points inside * area
    return static_cast<T>(insideCount) / numSamples * area;
}

// Example: Monte Carlo Approximation of Pi
bool isInsideCircle(double x, double y) {
    return x * x + y * y <= 1; // Check if point is inside the unit circle
}

int main() {
    int numSamples = 1000000; // Number of random samples

    // Approximate Pi
    double piApproximation = monteCarloApproximation<double>(
            numSamples,
            isInsideCircle,
            -1.0, // Min range for x and y
            1.0   // Max range for x and y
    );

    std::cout << "Approximated Pi: " << piApproximation << std::endl;
    return 0;
}
