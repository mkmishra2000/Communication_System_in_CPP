#include <iostream>
#include <iterator>
#include <random>
#include <chrono>

int main() {
    // Example data
    std::vector<double> data;

    // Define random generator with Gaussian distribution
    const double mean = 0.0;
    const double stddev = 1;
    
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);

    std::normal_distribution<double> dist(mean, stddev);

    // Add Gaussian noise
    for (int i =0; i<100; i++) {
        data.insert(data.end(),dist(generator));
    }

    // Output the result, for demonstration purposes
    std::copy(begin(data), end(data), std::ostream_iterator<double>(std::cout, "    "));
    std::cout << "\n";

    return 0;
}