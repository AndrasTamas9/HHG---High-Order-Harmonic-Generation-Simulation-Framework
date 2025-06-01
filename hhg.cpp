/**
 * @file hhg.cpp
 * @brief Single-atom high-harmonic generation (HHG) simulation code.
 *
 * Computes the nonlinear dipole response of an atom in a spatially- and temporally-resolved laser field.
 * Includes numerical integration, Fourier transforms, and core physics models.
 */

#define _USE_matH_DEFINES ///< Enable additional math constants

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <complex.h>
#include <omp.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_complex.h>

using namespace std;
using namespace std::complex_literals;

// Physical and numerical constants
const double PI = 3.1415926535;
const int SIZE = 10;       ///< Spatial grid size (radial and longitudinal)
//const int TSIZE = 16384;     ///< Time grid size, initialized later
int TSIZE;

// Main physical variables
double lambda, w0, E0, tau, phi, Ip, eps;

// Derived or secondary variables
double z0, omega;

/**
 * @brief Squares each element of the input vector.
 *
 * @param mat Input vector
 * @return Vector with squared values
 */
vector<double> squareVectorElements(const vector<double>& mat) {
    // Get input vector dimensions
    size_t t_size = mat.size();   // Time dimension
    
    vector<double> squared(t_size, 0.0);
    // Compute the square of each element
    for (size_t t = 0; t < t_size; ++t) {
        squared[t] = mat[t] * mat[t];
    }
    return squared;
}

/**
 * @brief Multiplies each element of the input vector by a scalar.
 *
 * @param mat Input vector
 * @param mul Scalar multiplier
 * @return Vector with scaled values
 */
vector<double> multiplyVectorElements(const vector<double>& mat, double mul) {
    // Get input vector dimensions
    size_t t_size = mat.size();   // Time dimension
    
    vector<double> new_mat(t_size, 0.0);
    // Compute the square of each element
    for (size_t t = 0; t < t_size; ++t) {
        new_mat[t] = mat[t] * mul;
    }
    return new_mat;
}

/**
 * @brief Computes the primitive (indefinite integral) of a function using Simpson's rule.
 *
 * Applies Simpson's rule cumulatively across the time grid to build a vector of primitive values.
 * Also handles the final interval using the trapezoidal rule if the vector size is even.
 *
 * @param f Input function sampled on a uniform time grid
 * @param h Grid spacing (time step)
 * @return Vector of primitive values
 */
vector<double> computePrimitiveSimpson3D(
    const vector<double>& f, double h) { 
    // Get input vector dimensions
    size_t t_size = f.size();   // Time dimension

    // Initialize the output grid with the same dimensions as f
    vector<double> F(t_size, 0.0);

    // Apply Simpson's rule for the time dimension at each spatial point
    for (size_t k = 2; k < t_size; k += 2) {
        F[k] = F[k - 2] + (h / 3.0) * (f[k - 2] + 4.0 * f[k - 1] + f[k]);
    }

    // Handle the last interval if t_size is even
    if (t_size % 2 == 0) {
        F[t_size - 1] = F[t_size - 2] + 0.5 * h * (f[t_size - 2] + f[t_size - 1]);
    }

    // Fill in intermediate values for continuity
    for (size_t k = 1; k < t_size; k += 2) {
        F[k] = F[k - 1]; // Copy the previous cumulative value
    }

    return F;
}

/**
 * @brief Computes the definite integral of a vector using Simpson's rule.
 *
 * Falls back to trapezoidal rule for the last interval if number of subintervals is odd.
 *
 * @param f Function values sampled on a uniform grid
 * @param startIndex Index of integration start
 * @param endIndex Index of integration end
 * @param lower Lower limit of integration (corresponds to startIndex)
 * @param upper Upper limit of integration (corresponds to endIndex)
 * @return Approximated integral value
 */
double simpsonIntegration(const vector<double>& f, int startIndex, int endIndex, double lower, double upper) {
    // Ensure startIndex and endIndex are valid
    if (startIndex < 0 || endIndex >= f.size() || startIndex > endIndex) {
        cout << "Invalid start or end indices." << endl;
        return 0;
    }
    if(startIndex == endIndex) return 0.0;

    // Number of intervals
    int n = endIndex - startIndex;
    double h;
    if(n == 0 || n == 1)
        h = upper - lower;
    else
        h = (upper-lower)/n;
    double integral = 0.0;

    // If n is odd, use Simpson's rule for the first n-1 intervals and trapezoidal rule for the last interval
    if (n % 2 != 0) {
        // Apply Simpson's rule for the first n-1 intervals
        integral += f[startIndex] + f[startIndex + n - 1];
        for (int i = 1; i < n - 1; i++) {
            if (i % 2 == 0) {
                integral += 2.0 * f[startIndex + i]; // Even indices
            } else {
                integral += 4.0 * f[startIndex + i]; // Odd indices
            }
        }
        integral *= (h / 3.0);

        // Apply trapezoidal rule for the last interval
        double lastInterval = (h / 2.0) * (f[startIndex + n - 1] + f[endIndex]);
        integral += lastInterval;
    }
    else {
        // Apply Simpson's rule for all intervals if n is even
        integral += f[startIndex] + f[endIndex];
        for (int i = 1; i < n; i++) {
            if (i % 2 == 0) {
                integral += 2.0 * f[startIndex + i]; // Even indices
            } else {
                integral += 4.0 * f[startIndex + i]; // Odd indices
            }
        }
        integral *= (h / 3.0);
    }

    return integral;
}

/**
 * @brief Computes the spatial phase of the Gaussian laser beam at (r,z).
 *
 * This accounts for the Gouy phase shift and spatial curvature of the beam.
 *
 * @param r Radial coordinate (in a.u.)
 * @param z Longitudinal coordinate (in a.u.)
 * @return Phase value at position (r,z)
 */
double spatialPhase(double r, double z) {
    return atan2(z, z0) - r * r /  (w0 * w0) * z * z0 / (z * z + z0 * z0);
}

/**
 * @brief Computes the spatial amplitude of the Gaussian laser field at (r,z).
 *
 * Includes beam divergence with propagation (Rayleigh length).
 *
 * @param r Radial coordinate (in a.u.)
 * @param z Longitudinal coordinate (in a.u.)
 * @return Field amplitude at position (r,z)
 */
double spatialAmplitude(double r, double z) {
    return E0 / sqrt(1.0 + z * z / (z0 * z0)) * exp(-r * r / (w0 * w0));
}

/**
 * @brief Computes the value of the laser electric field E(r,z,t).
 *
 * Combines spatial amplitude, temporal envelope, and carrier oscillation with phase.
 *
 * @param r Radial coordinate (a.u.)
 * @param z Longitudinal coordinate (a.u.)
 * @param t Time (a.u.)
 * @return Electric field value at (r,z,t)
 */
double laserImpulse(double r, double z, double t) {
    return spatialAmplitude(r, z) * exp(- t * t / (tau * tau)) * cos(-omega * t - spatialPhase(r, z) + phi);
}

/**
 * @brief Writes simulation data to an output file for a fixed (r,z) location.
 *
 * Outputs time-resolved values of a quantity (e.g., E(t), A(t), D(t)).
 *
 * @param mat Vector of time-dependent data
 * @param r Radial position
 * @param z Longitudinal position
 * @param fout Output file stream
 */
void outData(const vector<double>& mat, double r, double z, ofstream &fout) {
    // Get input vector dimensions
    size_t t_size = mat.size();
    fout << r << "\t" << z << "\n";
    for(size_t j = 0; j < t_size; j++)
    {
        fout << mat[j] << "\t";
    }
    fout << "\n";
}

/**
 * @brief Computes the average vector potential over a time interval.
 *
 * This approximates the canonical momentum of the electron.
 *
 * @param A Vector potential A(t)
 * @param startIndex Start index in time
 * @param endIndex End index in time
 * @param lower Start time
 * @param upper End time
 * @return Average vector potential over [lower, upper]
 */
double stationaryMomentum(const vector<double>& A, int startIndex, int endIndex, double lower, double upper) {
    if(startIndex == endIndex)
        return 0.0;
    return 1.0 / (upper - lower) * simpsonIntegration(A, startIndex, endIndex, lower, upper);
}

/**
 * @brief Computes the stationary momentum and action integral for an electron trajectory.
 *
 * These quantities are derived from saddle-point approximation in the Lewenstein model.
 *
 * @param momentum Output: stationary canonical momentum
 * @param action Output: quasiclassical action
 * @param A Vector potential A(t)
 * @param A_sq Squared vector potential A(t)^2
 * @param startIndex Start index in time
 * @param endIndex End index in time
 * @param lower Start time (a.u.)
 * @param upper End time (a.u.)
 */
void stationaryQuantities(double &momentum, double &action, const vector<double>& A,
    const vector<double>& A_sq, int startIndex, int endIndex, double lower, double upper) {
    momentum = stationaryMomentum(A, startIndex, endIndex, lower, upper);
    action = (upper - lower) * (Ip - 0.5 * momentum * momentum) + 0.5 * simpsonIntegration(A_sq, startIndex, endIndex, lower, upper);
}

/**
 * @brief Computes the dipole matrix element for a given electron momentum.
 *
 * This is part of the Lewenstein model and contributes to the nonlinear dipole.
 *
 * @param p Electron momentum
 * @return Value of the dipole matrix element
 */
double dipoleMatrixElement(double p) {
    double c = pow(2.0, 4.75) * pow(Ip, 1.25) / PI;
    double d = p * p + 2.0 * Ip;
    return c * p / (d * d * d);
}

/**
 * @brief Computes the depletion term due to ionization at a given time index.
 *
 * Based on Ammosov–Delone–Krainov (ADK) tunneling ionization rate approximation.
 *
 * @param E Electric field vector
 * @param t Time index
 * @return Depletion contribution at time t
 */
double depletionTerm(const vector<double>& E, int t) {
    return  16.0 * Ip * Ip / (abs(E[t]) / sqrt(2.0 * Ip)) * exp(-4.0 * Ip / (3.0 * abs(E[t]))); // last term in exp should be divided by omega_l
}

/**
 * @brief Computes the nonlinear dipole moment D(t) for an atom using the Lewenstein model.
 *
 * Integrates over past ionization times to calculate the quantum path contributions.
 *
 * @param E Electric field vector E(t)
 * @param A Vector potential A(t)
 * @param A_sq Squared vector potential A(t)^2
 * @param t0 Initial time
 * @param t Evaluation time
 * @param dt Time step
 * @param r Radial coordinate (used for diagnostic output)
 * @param z Longitudinal coordinate (used for diagnostic output)
 * @return Nonlinear dipole moment D(t)
 */
double nonlinearDipoleMoment(const vector<double>& E, const vector<double>& A,
    const vector<double>& A_sq, double t0, double t, double dt, int r, int z) {
    // Get input vector dimensions
    size_t t_size = E.size();   // Time dimension

    vector<double> D(t_size, 0.0);
    vector<double> DepTerm(t_size, 0.0);

    double stmomentum, staction, realcoef;
    double t1, t2;
    double t_curr;
    int index = (t - t0) / dt;

    if (index > t_size)
    {
        cout << "Invalid coodinates.\n";
        return 0;
    }

    for (size_t curr_index = 0; curr_index <= index; ++curr_index) {
        t_curr = t0 + curr_index * dt;

        //calculate the stationary momentum and quasiclassical action
        stationaryQuantities(stmomentum, staction, A, A_sq, curr_index, index, t_curr, t);
        //calculate the real part of complex factors
        complex<double> complexcoef = PI / (eps + 1i * (t - t_curr) / 2.0);
        complex<double> phasefactor = exp(-staction * 1i);
        complex<double> complexnumber = pow(complexcoef, 1.5) * phasefactor * 1i;
        realcoef = real(complexnumber);
        //calculate the dipole matrix elements - product is real!
        double d1 = dipoleMatrixElement(stmomentum - A[index]);
        double d2 = dipoleMatrixElement(stmomentum - A[curr_index]);

        //calculate the dipole integrand for every t_curr from t0 to t
        D[curr_index] = realcoef * d1 * d2 * E[curr_index];
        //calculate the depletion integrand for every t_curr from t0 to t
        DepTerm[curr_index] = depletionTerm(E, curr_index);
    }

    double D_nl = 2.0 * simpsonIntegration(D, 0, index, t0, t);
    double DT = -1.0 * simpsonIntegration(DepTerm, 0, index, t0, t);

    return D_nl * exp(DT);
}

/**
 * @brief Computes the real FFT of the nonlinear dipole moment D(t) using GSL.
 *
 * Converts the result to a magnitude spectrum.
 *
 * @param D_mat Time-domain dipole moment data
 * @param dt Time step size
 * @param spectrum Output vector with frequency spectrum (magnitude only)
 */
void computeFFT_Dmat(const vector<double>& D_mat, double dt, vector<double>& spectrum) {
    size_t N = D_mat.size();
    
    // Copy data into a modifiable array for GSL
    vector<double> data(D_mat); // this will be transformed in-place

    // Create FFT workspace and wavetable
    gsl_fft_real_workspace* work = gsl_fft_real_workspace_alloc(N);
    gsl_fft_real_wavetable* real = gsl_fft_real_wavetable_alloc(N);

    // Perform real FFT
    gsl_fft_real_transform(data.data(), 1, N, real, work);

    // Optional: convert halfcomplex output to complex for magnitude analysis
    vector<gsl_complex> complex_data(N);
    gsl_fft_halfcomplex_unpack(data.data(), reinterpret_cast<double*>(complex_data.data()), 1, N);

    // Prepare output vectors: magnitude spectrum
    spectrum.resize(N / 2);

    for (size_t i = 0; i < N / 2; ++i) {
        double re = GSL_REAL(complex_data[i]);
        double im = GSL_IMAG(complex_data[i]);
        double mag = sqrt(re * re + im * im);
        spectrum[i] = mag;
    }
    // Free resources
    gsl_fft_real_wavetable_free(real);
    gsl_fft_real_workspace_free(work);
}

/**
 * @brief Sets a global stop flag to interrupt long simulations.
 *
 * Meant to be called from external GUI/controller.
 */
extern "C" {
    volatile bool stop_requested = false;
    void requestStop() {
        stop_requested = true;
    }
}


/**
 * @brief Core entry point for HHG simulation from external (e.g. GUI) calls.
 *
 * Accepts physical parameters, initializes grids, computes E(t), A(t), D(t), and spectrum,
 * and saves them to output files.
 *
 * @param wavelength Laser wavelength in nm
 * @param waist Beam waist in microns
 * @param intensity Peak intensity in 10^14 W/cm²
 * @param T Pulse duration in fs
 * @param exponent TSIZE = 2^exponent (temporal resolution)
 */
extern "C" void singleAtomResponse(double wavelength, double waist, double intensity, double T, int exponent)
{
    stop_requested = false;
    
    TSIZE = pow(2, exponent);

    // Conversions
    lambda = wavelength / 0.0529;           // a.u.
    w0 = waist * 1000 / 0.0529;             // a.u.
    E0 = sqrt(intensity / 350.68);          // a.u.
    tau = T * 41.0;                         // a.u.
    z0 = PI * w0 * w0 / lambda;             // a.u.
    omega = 2.0 * PI / lambda * 137.036;    // a.u.
    phi = 0;
    Ip = 0.5;                               // a.u.
    eps = 1E-8;
    double oc = lambda / 137.036;
    
    // Print out computation values
    printf("Wavelength: %f nm = %f a.u.\n", wavelength, lambda);
    printf("Waist: %f mm = %f a.u.\n", waist/1000.0, w0);
    printf("Electric field: %f a.u.\n", E0);
    printf("Tau: %f fs = %f a.u.\n", T, tau);
    printf("Phi: %f rad\n", phi);
    printf("Z0: %f a.u.\n", z0);
    printf("omega: %f a.u.\n", omega);
    printf("oc: %f a.u\n", oc);
    printf("---------------------\n");

    // Define computational grid
    double r0_coord = -16.0;
    double z0_coord = -16.0;
    double t0_coord = -16.0 * oc; // Total time: 32 optical cycles
    double dr = 2.0*fabs(r0_coord)/SIZE;
    double dz = 2.0*fabs(z0_coord)/SIZE;
    double dt = 2.0*fabs(t0_coord)/TSIZE;
    double df = 1.0 / (-2.0 * t0_coord);
    t0_coord += dt/2.0;
    // Print the grid definition
    printf("t=[-%f, %f], dt=%f\n", fabs(t0_coord), fabs(t0_coord), dt);
    printf("r=[-%f, %f], dr=%f\n", fabs(r0_coord), fabs(r0_coord), dr);
    printf("z=[-%f, %f], dz=%f\n", fabs(z0_coord), fabs(z0_coord), dz);
    printf("df = %f\n", df);
    double r_coord, z_coord, t_coord;
    // Output relevant information about the grid to file
    ofstream fout0("results/rzt_data.txt");
    fout0 << r0_coord << "\t" << dr << "\t";
    fout0 << z0_coord << "\t" << dz << "\t";
    fout0 << t0_coord << "\t" << dt << "\t";
    fout0 << SIZE + 1 << "\t" << TSIZE << "\n";
    fout0.close();

    ofstream fout_E("results/E_"+to_string(exponent)+".txt");
    ofstream fout_A("results/A_"+to_string(exponent)+".txt");
    ofstream fout_D("results/D_"+to_string(exponent)+".txt");
    ofstream fout_FFT("results/FFT_"+to_string(exponent)+".txt");

    // Initialize electric field vector
    vector<double> E_mat(TSIZE, 0.0);

    // Initialize the nonlinear dipole moment vector
    vector<double> D_mat(TSIZE, 0.0);
    
    if(stop_requested) cout << "Stop requested\n";

    for(int i = SIZE/2; i <= SIZE; ++i)
    {
        if (stop_requested) return;
        
        r_coord = r0_coord + i * dr;
        for(int j = 0; j <= SIZE; ++j)
        {
            if (stop_requested) return;
            
            z_coord = z0_coord + j * dz;

            // Compute the value of the electric field in each gridpoint
            for(int k = 0; k < TSIZE; ++k)
            {
                if (stop_requested) return;
                
                t_coord = t0_coord + k * dt;
                E_mat[k] = laserImpulse(r_coord, z_coord, t_coord);
            }

            // Output relevant information about the electric field
            outData(E_mat, r_coord, z_coord, fout_E);

            // Calculate the associated vector potential
            vector<double> A = computePrimitiveSimpson3D(E_mat, dt);
            vector<double> A_mat = multiplyVectorElements(A, -1.0);
            // Calculate the squared vector potential
            vector<double> A_mat_squared = squareVectorElements(A_mat);
            // Output relevant information about the vector potential
            outData(A_mat, r_coord, z_coord, fout_A);

            // Calculate the dipole moment along the r=0, z=0 axis - for every t
            #pragma omp parallel for
            for(int k = 0; k < TSIZE; ++k)
            {
                if (stop_requested) continue;
                
                t_coord = t0_coord + k * dt;
                D_mat[k] = nonlinearDipoleMoment(E_mat, A_mat, A_mat_squared, t0_coord, t_coord, dt, r_coord, z_coord);
            }

            // Calculate the fourier transform of the dipole moment
            vector<double> spectrum;
            computeFFT_Dmat(D_mat, dt, spectrum);

            // Output fft of the dipole moment
            fout_FFT << r_coord << "\t" << z_coord << "\n";
            for (size_t i = 0; i < spectrum.size(); ++i) {
                fout_FFT << spectrum[i] << "\t";
            }
            fout_FFT << "\n";

            // Output relevant information about the dipole moment
            outData(D_mat, r_coord, z_coord, fout_D);
        }
    }

    fout_E.close();
    fout_A.close();
    fout_D.close();
}




#ifdef BUILD_EXECUTABLE
/**
 * @brief Standalone CLI mode for testing the HHG simulation.
 *
 * Expects arguments:
 *   wavelength (nm), waist (μm), intensity (10^14 W/cm²), T (fs), exponent (TSIZE=2^n)
 *
 * @param argc Number of arguments
 * @param argv Argument values
 * @return Exit code
 */
int main(int argc, char* argv[])
{   
    if (argc != 6) {
        std::cerr << "Usage: ./hhg <wavelength> <waist> <intensity> <T> <exponent>\n";
        return 1;
    }

    double wavelength = atof(argv[1]);
    double waist = atof(argv[2]);
    double intensity = atof(argv[3]);
    double T = atof(argv[4]);
    int exponent = atoi(argv[5]);

    singleAtomResponse(wavelength, waist, intensity, T, exponent);
    return 0;
}
#endif
