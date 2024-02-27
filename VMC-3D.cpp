#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

// PDF
double pdf(double r, double a) {
    //return std::exp(-2 * a * std::pow(r, 2)); //PDF for 3D Harmonic oscillator
    return std::exp(-2*a*r);  //PDF for hydrogen atom
}

// Local Energy
double EL(double a, double r) {
    //return 3*a + (0.5 - 2*std::pow(a, 2)) * std::pow(r, 2); //Local energy for 3D harmonic oscillator
    return -(1/r) - (1/2)*a*(a - (2/r)); //Local energy for hydrogen atom
}

std::vector<double> thermalise(double a, double h, int Nthermalsteps) {
    std::srand(std::time(0)); // Seed RNG
    //Initialise positions in 3D
    double xi = static_cast<double>(std::rand()) / RAND_MAX * 10.0 - 5.0; //Uniform distribution in [-5, 5]
    double yi = static_cast<double>(std::rand()) / RAND_MAX * 10.0 - 5.0; 
    double zi = static_cast<double>(std::rand()) / RAND_MAX * 10.0 - 5.0; 

    int accepted_moves = 0;

    for (int i = 0; i < Nthermalsteps; ++i) {
        //Propose step
        double xtrial = xi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; //Add normalised value in [-0.5, 0.5] scaled by h to xi
        double ytrial = yi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; 
        double ztrial = zi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h;

        double ri = sqrt(std::pow(xi, 2) + std::pow(yi, 2) + std::pow(zi, 2));
        double rtrial = sqrt(std::pow(xtrial, 2) + std::pow(ytrial, 2) + std::pow(ztrial, 2));

        //Acceptance criteria
        double p = pdf(rtrial, a) / pdf(ri, a);
        if (p >= 1.0 || static_cast<double>(std::rand()) / RAND_MAX < p) { //Accept move if r >= 1 OR r > Random number
            xi = xtrial;
            yi = ytrial;
            zi = ztrial;
            accepted_moves += 1;
        }

        //Adaptive step size
        h *= (accepted_moves / static_cast<double>(i + 1) > 0.5) ? 1.005 : 0.995;
    }
    return {xi, yi, zi}; //Return thermalised position
}

std::vector<double> gen_samples(double a, double xi, double yi, double zi, double h, int Nsteps) {
    
    std::srand(std::time(0)); //Seed RNG
    std::vector<double> samples;
    int accepted_moves = 0;

    for (int i = 0; i < Nsteps; ++i) {
        //Propose step
        double xtrial = xi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; //Add normalised value in [-0.5, 0.5] scaled by h to xi
        double ytrial = yi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; 
        double ztrial = zi + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h;

        double ri = sqrt(std::pow(xi, 2) + std::pow(yi, 2) + std::pow(zi, 2));
        double rtrial = sqrt(std::pow(xtrial, 2) + std::pow(ytrial, 2) + std::pow(ztrial, 2));

        //Acceptance Criteria
        double p = pdf(rtrial, a) / pdf(ri, a);
        if (p >= 1.0 || static_cast<double>(std::rand()) / RAND_MAX < p) { //Accept move if r >= 1 OR r > Random number
            xi = xtrial;
            yi = ytrial;
            zi = ztrial;
            accepted_moves += 1;
        }
        //Adaptive step size
        h *= (accepted_moves / static_cast<double>(i + 1) > 0.5) ? 1.005 : 0.995;
        samples.push_back(sqrt(std::pow(xi, 2) + std::pow(yi, 2) + std::pow(zi, 2))); //Add r values to vector
    }

    double acceptance_rate = static_cast<double>(accepted_moves) / Nsteps;
    return samples; //Return r values
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

    std::vector<std::vector<double>> starting_pos(num_walkers); //Starting positions in xyz obtained from thermalisation
    std::vector<std::vector<std::vector<double>>> samples(num_walkers); //Positions visited by walkers stored here
    std::vector<std::vector<double>> ET(num_walkers, std::vector<double>(points, 0.0));

    for (int w = 0; w < num_walkers; ++w) {
        for (int i = 0; i < points; ++i) {

            //Thermalize
            std::vector<double> thermalized_pos = thermalise(alpha[i], h, Nthermalsteps);
            starting_pos[w].insert(starting_pos[w].end(), thermalized_pos.begin(), thermalized_pos.end());

            //Generate samples
            std::vector<double> generated_samples = gen_samples(alpha[i], thermalized_pos[0], thermalized_pos[1], thermalized_pos[2], h, Nsteps);
            samples[w].push_back(generated_samples);

            //Calculate local energy
            for (int j = 0; j < generated_samples.size(); ++j) {
                ET[w][i] += EL(alpha[i], generated_samples[j]); //Sum local energies of each sample for that alpha
            }

            ET[w][i] /= generated_samples.size(); //Average sum to get local energy
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