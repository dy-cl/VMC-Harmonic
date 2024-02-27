#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>

//Normalise vector
std::vector<double> normalize(const std::vector<double>& r) {
    /*
    r: 3D vector
    */
    double norm = 0.0;
    //Calculate the sum of squares
    for (int i = 0; i < r.size(); i++) {
        norm += r[i] * r[i];
    }
    //Calculate the square root of the sum of squares
    norm = std::sqrt(norm);
    // Avoid division by zero
    if (norm != 0.0) {
        std::vector<double> result;
        //Normalize each component and push it into the result vector
        for (int i = 0; i < r.size(); i++) {
            result.push_back(r[i] / norm);
        }
        return result;
    } else {
        //Return the original vector if its magnitude is zero
        return r;
    }
}

//Dot product of vectors
double dot_product(const std::vector<double>& r1, const std::vector<double>& r2) {
    /*
    r1 : 3D vector
    r2 : 3D vector
    */
    double product = 0;
    for(int i = 0; i < r1.size(); i++)
        product += r1[i] * r2[i];
    return product;
}

//Difference between vectors
std::vector<double> subtract_vectors(const std::vector<double>& r1, const std::vector<double>& r2) {
    /*
    r1 : 3D vector
    r2 : 3D vector
    */
    std::vector<double> result;
    for (size_t i = 0; i < r1.size(); ++i) {
        result.push_back(r1[i] - r2[i]);
    }

    return result;
}

//Convert 3D cartesian to 1D spherical polar coordinates
double cartesianToPolar(std::vector<double>& r) {
    /*
    r : 3D Vector
    */
    double coord = 0.0;
    for (int i = 0; i < r.size(); i++) {
        coord += r[i] * r[i];
    }
    return std::sqrt(coord);
}

//Distance between vectors
double distance(std::vector<double>& r1, std::vector<double>& r2) {
    /*
    r1 : 3D Vector
    r2 : 3D Vector
    */
    std::vector<double> diff = subtract_vectors(r2, r1);
    return cartesianToPolar(diff);
}

//PDF
double pdf(double r1, double r2, double r12, double a) {
    /*
    a   : variational parameter
    r1  : 1st electron position
    r2  : 2nd electron position
    r12 : distance between electrons
    parameters r1, r2, r12 are spherical polar coordinates
    */
    return std::exp(-4*r1)*std::exp(-4*r2)*std::exp(r12 / (2 * (1 + a*r12)));
}

//Local Energy
double EL(std::vector<double>& r1, std::vector<double>& r2, double r12, double a) {
    /*
    a   : Variational parameter
    r1  : 3D Vector - 1st electron position
    r2  : 3D Vector - 2nd electron position
    r12 : distance between electrons
    parameter r12 is spherical polar coordinates
    */ 
    //Construct parts of local energy expression
    std::vector<double> n_r1 = normalize(r1);
    std::vector<double> n_r2 = normalize(r2);
    std::vector<double> norm_diff = subtract_vectors(n_r1, n_r2);
    double dot_result = dot_product(norm_diff, subtract_vectors(r1, r2));
    double deno1 = r12 * std::pow(1 + a*r12, 2);
    double deno2 = r12 * std::pow(1 + a*r12, 3);
    double deno3 = 4 * std::pow(1 + a*r12, 4);

    //Local energy
    return -4 + (dot_result / deno1) - (1 / deno2) - (1 / deno3) + (1 / r12);
}

//Thermalise
std::pair<std::vector<double>, std::vector<double>> thermalise(double a, double h, int Nthermalsteps) {
    /*
    a : variational parameter
    h : step size
    Nthermalsteps: number of iterations for thermalisation
    */
    std::srand(std::time(0)); //Seed RNG
    //Initialise blank vectors
    std::vector<double> r1(3);
    std::vector<double> r2(3);
    std::vector<double> r1trial(3);
    std::vector<double> r2trial(3);
    int accepted_moves = 0;

    //Initialize positions in 3D 
    for (size_t i = 0; i < 3; ++i) {
        r1[i] = -5.0 + (std::rand() % 1000) / 100.0;  // Random value in [-5, 5]
        r2[i] = -5.0 + (std::rand() % 1000) / 100.0;  // Random value in [-5, 5]
    }
    
    //Propose step
    for (int i = 0; i < Nthermalsteps; ++i) {
        for (int j = 0; j < r1.size(); j++) {
            r1trial[j] = r1[j] + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; //Add normalised value in [-0.5, 0.5] scaled by h
            r2trial[j] = r2[j] + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; 
        }

        double r12 = distance(r1, r2);
        double r12trial = distance(r1trial, r2trial);

        //Acceptance criteria
        double p = pdf(cartesianToPolar(r1trial), cartesianToPolar(r2trial), r12trial, a) / pdf(cartesianToPolar(r1), cartesianToPolar(r2), r12, a);
        if (p >= 1.0 || static_cast<double>(std::rand()) / RAND_MAX < p) { //Accept move if r >= 1 OR r > Random number
                r1 = r1trial;
                r2 = r2trial;
                accepted_moves += 1;
            }
            //Adaptive step size
            h *= (accepted_moves / static_cast<double>(i + 1) > 0.5) ? 1.005 : 0.995;
        
    }

    //Return the pair of vectors
    return std::make_pair(r1, r2);
}

//Generate samples
std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> gen_samples(
    std::vector<double>& r1, std::vector<double>& r2, double a, double h, int Nsteps) {
    /*
    r1      : 3D vector 1 generated via thermalisation process
    r2      : 3D vector 2 generated via thermalisation process
    a       : variational parameter
    h       : step size
    Nsteps  : number of iterations for sample generation
    */

    std::srand(std::time(0)); // Seed RNG
    std::vector<std::vector<double>> positions1; // Initialize vector to store all samples generated for electron 1
    std::vector<std::vector<double>> positions2; // Initialize vector to store all samples generated for electron 2
    std::vector<double> r1trial(3); // Initialize blank vectors
    std::vector<double> r2trial(3);
    int accepted_moves = 0;

    // Propose step
    for (int i = 0; i < Nsteps; ++i) {

        //Store the current positions in the vectors
        positions1.push_back(r1);
        positions2.push_back(r2);

        for (int j = 0; j < r1.size(); j++) {
            r1trial[j] = r1[j] + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h; //Add normalized value in [-0.5, 0.5] scaled by h
            r2trial[j] = r2[j] + (static_cast<double>(std::rand()) / RAND_MAX - 0.5) * h;
        }

        double r12 = distance(r1, r2);
        double r12trial = distance(r1trial, r2trial);

        //Acceptance criteria
        double p = pdf(cartesianToPolar(r1trial), cartesianToPolar(r2trial), r12trial, a) /
                   pdf(cartesianToPolar(r1), cartesianToPolar(r2), r12, a);
        if (p >= 1.0 || static_cast<double>(std::rand()) / RAND_MAX < p) { //Accept move if r >= 1 OR r > Random number
            r1 = r1trial;
            r2 = r2trial;
            accepted_moves += 1;
        }
        //Adaptive step size
        h *= (accepted_moves / static_cast<double>(i + 1) > 0.5) ? 1.005 : 0.995;
    }

    //Return a pair of vectors containing all sampled 3D Vectors for each electron
    return std::make_pair(positions1, positions2);
}


int main() {
  
    int h = 2; //Step Size
    int Nthermalsteps = 4000; //Thermalisation steps
    int Nsteps = 26000; //VMC steps

    //Generate variational parameter alpha
    double alpha1 = 0.05;
    double alpha2 = 0.25;
    int points = 20;
    double step = 0.025;
    std::vector<double> alpha;
    for (int i = 0; i < points; ++i) {
        double val = alpha1 + i * step;
        alpha.push_back(val);
    }

    int num_walkers = 400; //Number of walkers

    //Vector of length num_walkers of 3D vectors to store starting positions for each walker for electron 1
    std::vector<std::vector<double>> starting_pos1(num_walkers);
    //Vector of length num_walkers of 3D vectors to store starting positions for each walker for electron 2
    std::vector<std::vector<double>> starting_pos2(num_walkers);
    //Vector of length num_walkers of vectors length Nsteps of 3D vectors to store each position visited for electron 1
    std::vector<std::vector<std::vector<double>>> samples1(num_walkers); 
    //Vector of length num_walkers of vectors length Nsteps of 3D vectors to store each position visited for electron 2
    std::vector<std::vector<std::vector<double>>> samples2(num_walkers); 

    for (int a = 0; a < alpha.size(); a++) {
        double total_energy = 0.0;

        for (int w = 0; w < num_walkers; w++) {
            //Thermalize
            std::pair<std::vector<double>, std::vector<double>> thermalized_pos = thermalise(alpha[a], h, Nthermalsteps);

            //Separate thermalized positions
            starting_pos1[w] = thermalized_pos.first;
            starting_pos2[w] = thermalized_pos.second;

            //Generate samples
            std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> samples = gen_samples(starting_pos1[w], starting_pos2[w], alpha[a], h, Nsteps);

            //Calculate and sum local energy for each sample
            double walker_energy = 0.0;
            for (int step = 0; step < Nsteps; step++) {
                double r12 = distance(samples.first[step], samples.second[step]);
                double local_energy = EL(samples.first[step], samples.second[step], r12, alpha[a]);

                walker_energy += local_energy;
            }

            //Average over the number of samples
            walker_energy /= Nsteps;

            //Accumulate the total energy for this walker
            total_energy += walker_energy;
        }

        //Average over the number of walkers
        total_energy /= num_walkers;

        //Output the total energy for the current alpha
        std::cout << "Total Energy for Alpha = " << alpha[a] << ": " << total_energy << std::endl;
    }

    return 0;
}