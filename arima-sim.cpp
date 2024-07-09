#include <iostream>
#include <vector>
#include <random>
#include <fstream>




static constexpr double bad_number = -9999.0;
static constexpr double epsilon = 0.000000000000000001;


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


std::vector<double> ArrayToVector(double* arr, size_t arr_len) {
    return std::vector<double>(arr, arr + arr_len);
}



// Function to simulate ARIMA process
std::vector<double> simulateARIMA(int p, int d, int q, double *phi, double *theta, int num_samples) {
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
    arvec = ArrayToVector(phi, p);
    mavec = ArrayToVector(theta, q);
    acfvec = derive_autocorrelation(arvec, mavec, k);
    double cc = 1.0;

    for (int i = 0; i < p; i++)
        cc -= phi[i] * acfvec[i];
    //cc = sqrt(cc);
    double sig2 = 2;
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


    return simulated_series;
}





int main() {
    // Define ARIMA parameters
    int p = 6; // Autoregressive terms
    int d = 0; // Difference terms
    int q = 0; // Moving average terms

    // ARIMA coefficients
    double phi[] = { 1.17, -0.321, 0.0785, -0.0182, -0.00601, 0.0724 };  // Autoregressive coefficients
    double theta[]{ 0.0 }; // Moving average coefficients

    // Number of samples to simulate
    int num_samples = 1000;

    // Simulate ARIMA process
    std::vector<double> simulated_series = simulateARIMA(p, d, q, phi, theta, num_samples);

    // Output simulated series
    //for (int i = 0; i < num_samples; ++i) {
    //    std::cout << simulated_series[i] << std::endl;
    //}


    std::ofstream outfile("output.txt");


    for (int i = 0; i < num_samples; i++)
    {

        double sample = simulated_series[i];
       // sample = sample * sig2;
        //sample = sample * cc;
        outfile << sample << std::endl;

        //file << MagWind << sample << '\n';

       // std::cout << sample << std::endl;
    }

    outfile.close();


    return 0;
}