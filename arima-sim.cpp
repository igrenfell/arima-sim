//This program will take the parameters set at the beginning of the main function (around line 230) and output a two column csv, the first being the simulated magnitudues
//and the second being the direction.

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <cmath>




static constexpr double bad_number = -9999.0;
static constexpr double epsilon = 0.000000000000000001;
static constexpr double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679;

bool isFilterOk(double value)
{
    return !(std::abs(value - bad_number) < epsilon);
}



void applyAutoregressiveFilter(std::vector<double>& arimaWindStream, std::vector<double>& autoregressiveParameters)
{
    std::vector<double> out(arimaWindStream.size() + autoregressiveParameters.size());

    int nf = (int)autoregressiveParameters.size();
    int nx = (int)arimaWindStream.size();

    double sum;
    double tmp;

    bool goodStatus = true;
    for (int i = 0; i < nx; i++)
    {
        goodStatus = true;
        sum = arimaWindStream[(unsigned int)i];
        for (int j = 0; j < nf; j++)
        {
            tmp = out[(unsigned int)(nf + i - j - 1)];
            if (isFilterOk(tmp))
            {
                sum += tmp * autoregressiveParameters[(unsigned int)j];
            }
            else
            {
                out[(unsigned int)(nf + i)] = bad_number;
                goodStatus = false;
                break;
            }
        }
        if (goodStatus)
        {
            out[(unsigned int)(nf + i)] = sum;
        }
        else
        {
            continue;
        }
    }

    out.erase(out.begin(), out.begin() + nf); // Delete junk data
    arimaWindStream = out;
}


void applyMovingAverageFilter(std::vector<double>& arimaWindStream, std::vector<double>& movingAverageParameters)
{
    int nshift = 0;
    int nf = (int)movingAverageParameters.size();
    int nx = (int)arimaWindStream.size();

    std::vector<double> out(arimaWindStream.size() + nf);
    double z;
    double tmp;

    bool goodStatus = true;
    for (int i = 0; i < (int)(arimaWindStream.size()); i++)
    {
        goodStatus = true;
        z = 0;
        if (i + nshift - (nf - 1) < 0 || i + nshift >= nx)
        {
            out[(unsigned int)i] = bad_number;
            continue;
        }
        for (int j = std::max(0, nshift + i - nx); j < std::min(nf, i + nshift + 1); j++)
        {
            tmp = arimaWindStream[(unsigned int)(i + nshift - j)];
            if (isFilterOk(tmp))
            {
                z += movingAverageParameters[(unsigned int)j] * tmp;
            }
            else
            {
                out[(unsigned int)i] = bad_number;
                goodStatus = false;
            }
        }
        if (goodStatus)
        {
            out[(unsigned int)i] = z;
        }
        else
        {
            continue;
        }
    }

    out.erase(out.begin(), out.begin() + nf); // Delete junk data
    arimaWindStream = out;
}


std::vector<double> derive_autocorrelation(const std::vector<double>& ar_coeffs, const std::vector<double>& ma_coeffs, size_t k) {
    size_t p = ar_coeffs.size();
    size_t q = ma_coeffs.size();
    size_t max_lag = std::max(p, q + k);

    std::vector<double> acf(max_lag + 1, 0.0);

    // Calculate autocorrelation coefficients using Yule-Walker equations for AR part
    for (size_t l = 1; l <= p; ++l) {
        acf[l] = ar_coeffs[l - 1];
        for (size_t j = 1; j < l; ++j) {
            acf[l] += ar_coeffs[j - 1] * acf[l - j];
        }
    }

    // Calculate autocorrelation coefficients using invertibility conditions for MA part
    for (size_t l = 1; l <= q && l <= k; ++l) {
        for (size_t j = 1; j < l; ++j) {
            acf[l] += ma_coeffs[j - 1] * acf[l - j];
        }
    }

    return acf;
}

double von_mises(double mean, double kappa, std::mt19937& generator) {
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
    std::normal_distribution<double> normal_dist(0.0, 1.0);

    double tau = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    double rho = (tau - std::sqrt(2.0 * tau)) / (2.0 * kappa);
    double r = (1.0 + rho * rho) / (2.0 * rho);

    while (true) {
        double u1 = uniform_dist(generator);
        double z = std::cos(pi * u1);
        double f = (1.0 + r * z) / (r + z);
        double c = kappa * (r - f);

        double u2 = uniform_dist(generator);
        if (u2 < c * (2.0 - c) || u2 <= c * std::exp(1.0 - c)) {
            double u3 = uniform_dist(generator);
            double theta = mean + std::copysign(1.0, u3 - 0.5) * std::acos(f);
            return theta;
        }
    }
}




std::vector<double> ArrayToVector(double* arr, size_t arr_len) {
    return std::vector<double>(arr, arr + arr_len);
}



// Function to simulate ARIMA process
std::vector<double> simulateARIMA(double *phi, double *theta, int num_samples, double sig2, double meanwind) {
    // Initialize variables
    std::vector<double> simulated_series(num_samples, 0.0);
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    // Generate ARIMA process

    // Initialize with random noise
    int k = 30;
    //std::vector<double> acf = derive_autocorrelation(phivec, thetavec, k);
    //std::vector<double> acf = { 1.000000,0.975384,0.941925,0.913408,0.890211,0.872260,0.858886,0.847181,0.834872,0.822290,0.810373,0.797528,0.784385,0.771865,0.759759,0.748905,
    //    0.737786,0.725777,0.713362,0.701186,0.689194,0.676517,0.663835,0.652762,0.642148,0.631309,0.621892,0.612033,0.600785,0.587710,0.575404,0.563168,0.551136,0.541185,0.531647,0.521629 };
    

    std::vector<double> arvec, mavec, acfvec;
    arvec = ArrayToVector(phi, sizeof(*phi));
    mavec = ArrayToVector(theta, sizeof(*theta));
    acfvec = derive_autocorrelation(arvec, mavec, k);
    double cc = 1.0;

    for (int i = 0; i < sizeof(*phi); i++)
        cc -= phi[i] * acfvec[i];
    //cc = sqrt(cc);
    double scfac =  sig2 / cc;
    scfac = sqrt(scfac);

    for (int i = 0; i < num_samples; i++)
    {

        double sample = distribution(generator);
        //sample = sample / sig2;
        sample = sample * scfac;
        simulated_series[i] = sample;
        //file << MagWind << sample << '\n'
       // std::cout << sample << std::endl;
    }

    applyAutoregressiveFilter(simulated_series, arvec);
    //applyMovingAverageFilter(simulated_series, thetavec);


    for (int i = 0; i < simulated_series.size(); i++){
        simulated_series[i] = simulated_series[i]+meanwind;
    }

    return simulated_series;
}

int main() {

    //mean and sd of wind magnitude
    double sig2 = 0.10;
    double meanwind = 5.0;
    // set all ARIMA coefficients
    double phi[] = { 1.17, -0.321, 0.0785, -0.0182, -0.00601, 0.0724 };  // Autoregressive coefficients
    double theta[]{ 0.1 }; // Moving average coefficients
    double vmmean = 1.5*pi; //von mises mean direction in radians counterclockwise from east
    double vmkappa = 20.0; //von mises shape parameter: 1.0 => uniform direction, higher values imply a "pointier" distribution

    // Number of samples to simulate
    int num_samples = 1000;

    // Simulate ARIMA process
    std::vector<double> simulated_series = simulateARIMA(phi, theta, num_samples, sig2, meanwind);
    std::random_device rd;
    std::mt19937 generator(rd());
    
    std::mt19937 gen(rd());
    // Output simulated series
    //for (int i = 0; i < num_samples; ++i) {lear
    //    std::cout << simulated_series[i] << std::endl;
    //}
    std::ofstream outfile("output.txt");

    for (int i = 0; i < num_samples; i++)
    {

        double sample = simulated_series[i];
        double vmsample = von_mises(vmmean, vmkappa, generator);
        double vmdegrees = 90 - (vmsample / pi) * 180;
        if (vmdegrees < 0) {
            vmdegrees += 360.0;
        }        
       // sample = sample * sig2;
        //sample = sample * cc;
        outfile << sample << "," << vmdegrees << std::endl;

        //file << MagWind << sample << '\n';

       // std::cout << sample << std::endl;
    }

    outfile.close();

    return 0;
}



