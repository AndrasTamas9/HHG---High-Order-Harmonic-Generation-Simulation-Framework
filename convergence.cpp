/**
 * @file convergence.cpp
 * @brief Performs convergence analysis for HHG dipole outputs (D_n.txt files).
 *
 * This code generates HHG dipole data at various temporal resolutions and compares
 * consecutive resolutions to estimate convergence by computing normalized error ratios.
 */
 
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <regex>

#include <dirent.h>
#include <cstdlib>
#include <cstdio>
#include <sstream>

using namespace std;

/**
 * @brief Flag to allow asynchronous interruption of long computations.
 */
extern "C" {
    volatile bool stop_requested = false;

    /**
     * @brief Signals that a stop has been requested.
     */
    void requestStop() {
        stop_requested = true;
    }
}

/**
 * @brief Calls the external HHG executable to generate D_n.txt files from n = 2 to 8.
 *
 * Each execution uses the provided laser parameters. Stops early if requested.
 *
 * @param wavelength Laser wavelength in nm
 * @param waist Beam waist in μm
 * @param intensity Intensity in 10^14 W/cm²
 * @param impulse Pulse duration in fs
 */
void generate_D_files(double wavelength, double waist, double intensity, double impulse) {
    stop_requested = false;

    for (int n = 2; n <= 10; ++n) {
        if (stop_requested) {
            cout << "Stop requested during generation.\n";
            return;
        }
        
        cout << "Running hhg for n = " << n << "..." << endl;

        // Build the system command to run hhg with laser parameters
        stringstream cmd;
        cmd << "./hhg "   //"./hhg"
            << wavelength << " "
            << waist << " "
            << intensity << " "
            << impulse << " "
            << n;

        int ret = system(cmd.str().c_str());
        if (ret != 0) {
            cerr << "Error running hhg for n = " << n << endl;
            exit(1);
        }
    }
}

/**
 * @brief Linearly interpolates a coarse vector onto a finer grid.
 *
 * Used to rescale coarse dipole vectors to the resolution of the fine reference.
 *
 * @param y_coarse Vector on coarse grid
 * @param n_fine Size of fine grid
 * @return Interpolated vector
 */
vector<double> interpolate(const vector<double>& y_coarse, size_t n_fine) {
    size_t n_coarse = y_coarse.size();
    vector<double> y_interp;
    for (size_t i = 0; i < n_fine; ++i) {
        double t = static_cast<double>(i) / (n_fine - 1);
        double idx = t * (n_coarse - 1);
        size_t i0 = static_cast<size_t>(floor(idx));
        size_t i1 = min(i0 + 1, n_coarse - 1);
        double w = idx - i0;
        double y = (1 - w) * y_coarse[i0] + w * y_coarse[i1];
        y_interp.push_back(y);
    }
    return y_interp;
}

/**
 * @brief Computes the area between two vectors using rectangular integration.
 *
 * @param a First vector (reference)
 * @param b Second vector (interpolated)
 * @return Area between the two vectors
 */
double area_between(const vector<double>& a, const vector<double>& b) {
    size_t N = min(a.size(), b.size());
    double dt = 1.0 / (N - 1);
    double area = 0.0;
    for (size_t i = 0; i < N; ++i) {
        area += fabs(a[i] - b[i]) * dt;
    }
    return area;
}

/**
 * @brief Computes the L1 norm (area under absolute value) of a vector.
 *
 * @param y Input vector
 * @return Area under the curve
 */
double area_under(const vector<double>& y) {
    size_t N = y.size();
    double dt = 1.0 / (N - 1);
    double area = 0.0;
    for (size_t i = 0; i < N; ++i) {
        area += fabs(y[i]) * dt;
    }
    return area;
}


/**
 * @brief Parses a D_n.txt file into a map of (r,z) positions and time-series values.
 *
 * @param filename Path to the file
 * @return Map from (r,z) pairs to time series of D(t)
 */
map<pair<double, double>, vector<double>> read_data(const string& filename) {
    ifstream file(filename);
    map<pair<double, double>, vector<double>> data;
    string line;
    double r = 0.0, z = 0.0;

    while (getline(file, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        vector<double> values;
        double value;
        while (ss >> value) {
            values.push_back(value);
        }

        if (values.size() == 2) {
            r = values[0];
            z = values[1];
        } else if (!values.empty()) {
            data[{r, z}] = values;
        }
    }

    return data;
}

/**
 * @brief Extracts the timestep exponent (n) from a filename of the form "D_n.txt".
 *
 * @param filename File name string
 * @return Integer exponent
 * @throws runtime_error if pattern not matched
 */
int extract_n(const string& filename) {
    regex pattern(R"(D_(\d+)\.txt)");
    smatch match;
    if (regex_search(filename, match, pattern)) {
        return stoi(match[1]);
    }
    throw runtime_error("Could not extract timestep exponent from filename.");
}


/**
 * @brief Comparison function for sorting filenames by extracted timestep exponent.
 */
bool compare_by_n(const string& a, const string& b) {
    return extract_n(a) < extract_n(b);
}

/**
 * @brief Runs convergence analysis across multiple D_n.txt files.
 *
 * Compares each consecutive resolution level, interpolates coarse data,
 * and computes area ratios for all (r,z) points.
 *
 * @param wavelength Laser wavelength in nm
 * @param waist Beam waist in μm
 * @param intensity Peak intensity in 10^14 W/cm²
 * @param impulse Pulse length in fs
 */
void perform_convergence_analysis(double wavelength, double waist, double intensity, double impulse) {
    stop_requested = false;

    generate_D_files(wavelength, waist, intensity, impulse);

    vector<string> files;
    regex file_pattern(R"(D_\d+\.txt)");  // csak a fájlnevet nézzük

    string dir_path = "results";

    DIR* dir;
    struct dirent* ent;
    if ((dir = opendir(dir_path.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            string fname = ent->d_name;
            if (regex_match(fname, file_pattern)) {
                files.push_back(dir_path + "/" + fname);  // elérési út mentés
            }
        }
        closedir(dir);
    } else {
        perror("Could not open directory.");
        return;
    }

    if (files.size() < 2) {
        cerr << "Not enough D_n.txt files to compute convergence.\n";
        return;
    }

    sort(files.begin(), files.end(), compare_by_n);

    ofstream out("results/convergence_ratios.txt");
    out << "# fine_n coarse_n r z area_ratio\n";

    for (size_t i = 1; i < files.size(); ++i) {
    
        if (stop_requested) {
            cout << "Stop requested during generation.\n";
            return;
        }
        
        string fine_file = files[i];
        string coarse_file = files[i - 1];

        int n_fine = extract_n(fine_file);
        int n_coarse = extract_n(coarse_file);

        auto D_fine = read_data(fine_file);
        auto D_coarse = read_data(coarse_file);

        for (const auto& [key, D_fine_vals] : D_fine) {
        
            if (stop_requested) {
                cout << "Stop requested during generation.\n";
                return;
            }
            
            auto it = D_coarse.find(key);
            if (it == D_coarse.end()) continue;

            const auto& D_coarse_vals = it->second;
            auto D_interp = interpolate(D_coarse_vals, D_fine_vals.size());

            double diff_area = area_between(D_fine_vals, D_interp);
            double base_area = area_under(D_fine_vals);

            if (base_area > 1e-20) {
                double ratio = diff_area / base_area;
                out << fixed << setprecision(6)
                    << n_fine << " " << n_coarse << " "
                    << key.first << " " << key.second << " "
                    << ratio << "\n";
            }
        }
    }

    out.close();
    cout << "Wrote convergence data to results/convergence_ratios.txt\n";
}

/**
 * @brief C-compatible wrapper to run convergence analysis.
 *
 * This is designed to be callable from Python or GUI frontends.
 */
extern "C" void runConvergence(double wavelength, double waist, double intensity, double impulse) {
    stop_requested = false;

    perform_convergence_analysis(wavelength, waist, intensity, impulse);
}

#ifdef BUILD_EXECUTABLE
/**
 * @brief CLI entry point for convergence analysis.
 *
 * Usage: ./convergence <wavelength> <waist> <intensity> <impulse>
 */
int main(int argc, char* argv[]) {
    if (argc != 5) {
        cerr << "Usage: " << argv[0] << " wavelength waist intensity impulse\n";
        return 1;
    }

    double wavelength = stod(argv[1]);
    double waist = stod(argv[2]);
    double intensity = stod(argv[3]);
    double impulse = stod(argv[4]);

    perform_convergence_analysis(wavelength, waist, intensity, impulse);
    return 0;
}
#endif


