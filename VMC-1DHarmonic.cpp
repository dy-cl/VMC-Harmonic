#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

//PDF
double pdf(double r, double a) {
    return std::exp(-2 * a * std::pow(r, 2));
}

//Local Energy
double EL(double a, double r) {
    return -1.0 * a * (2 * a * std::pow(r, 2) - 1) + 0.5 * std::pow(r, 2); //Local energy for 1D harmonic oscillator
}

double thermalise(double a, double h, int Nthermalsteps) {
    std::srand(std::time(0)); //Seed RNG
    double xi = static_cast<double>(std::rand()) / RAND_MAX * 10.0 - 5.0; //Uniform disribution in [-5, 5]
    int accepted_moves = 0;

    for (int i = 0; i < Nthermalsteps; ++i) {
        //Propose step
        double xtrial = xi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; //Add normalised value in [-0.5, 0.5] scaled by h to xi
        //Acceptance criteria
        double r = pdf(xtrial, a) / pdf(xi, a);
        if (r >= 1.0 || static_cast<double>(std::rand()) / RAND_MAX < r) { //Accept move if r >= 1 OR r > Random number
            xi = xtrial;
            accepted_moves += 1;
        }
        //Adaptive step size
        h *= (accepted_moves / static_cast<double>(i + 1) > 0.5) ? 1.005 : 0.995;
    }
    return xi; //Return thermalised position
}

std::vector<double> gen_samples(double a, double xi, double h, int Nsteps) {
    std::srand(std::time(0)); //Seed RNG

    std::vector<double> samples;
    int accepted_moves = 0;

    for (int i = 0; i < Nsteps; ++i) {
        int sign = (std::rand() % 2 == 0) ? 1 : -1;
        //Propose step
        double xtrial = xi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; //Add normalised value in [-0.5, 0.5] scaled by h to xi
        //Acceptance Criteria
        double r = pdf(xtrial, a) / pdf(xi, a);
        if (r >= 1.0 || static_cast<double>(std::rand()) / RAND_MAX < r) { //Accept move if r >= 1 OR r > Random number
            xi = xtrial;
            accepted_moves += 1;
        }
        //Adaptive step size
        h *= (accepted_moves / static_cast<double>(i + 1) > 0.5) ? 1.005 : 0.995;
        samples.push_back(xi);
    }
    double acceptance_rate = static_cast<double>(accepted_moves) / Nsteps;
    return samples;
}

int main() {

    int h = 2; //Step Size
    int Nthermalsteps = 4000; //Thermalisation steps
    int Nsteps = 26000; //VMC steps

    //Generate variational parameter
    double alpha1 = 0.1;
    double alpha2 = 2.0;
    int points = 20;
    double step = 0.1;
    std::vector<double> alpha;
    for (int i = 0; i < points; ++i) {
        double val = alpha1 + i * step;
        alpha.push_back(val);
    }

    int num_walkers = 400; //Number of walkers
    std::vector<std::vector<double>> xi_list(num_walkers);
    std::vector<std::vector<std::vector<double>>> samples(num_walkers);
    std::vector<std::vector<double>> ET(num_walkers, std::vector<double>(points, 0.0));

    for (int w = 0; w < num_walkers; ++w) {
        for (int i = 0; i < points; ++i) {
            xi_list[w].push_back(thermalise(alpha[i], h, Nthermalsteps)); //List of thermalised positions for each alpha
            samples[w].push_back(gen_samples(alpha[i], xi_list[w][i], h, Nsteps)); //List of generated samples for each alpha

            for (int j = 0; j < samples[w][i].size(); ++j) {
                ET[w][i] += EL(alpha[i], samples[w][i][j]); //Sum local energies of each sample for that alpha
            }

            ET[w][i] /= samples[w][i].size(); //Average sum to get local energy
        }
    }

    //Calculate average local energy for all walkers combined
    std::vector<double> average_ET(points, 0.0);
    for (int w = 0; w < num_walkers; ++w) {
        for (int i = 0; i < points; ++i) {
            average_ET[i] += ET[w][i];
        }
    }

    for (int i = 0; i < points; ++i) {
        average_ET[i] /= num_walkers;
    }

    //Output results
    for (int i = 0; i < points; ++i) {
        std::cout << "a: " << alpha[i] << ", Average ET: " << average_ET[i] << std::endl;
    }

    return 0;
}