"""
High-level interface for running HHG simulations.

This module provides:

- run_simulation(): invokes the C++ HHG backend.
- run_convergence_test(): performs convergence diagnostics.
- Automatic loading of shared object libraries (.so) using ctypes.
- Optional plot generation from simulation output data.

Dependencies:

- Shared object files must be built in the parent directory:
  ../libhhg.so and ../libconvergence.so
- Plotting utilities are imported from the `plotting` module.
"""

import os
import ctypes
from .plotting import plot_time_series_from_file, plot_convergence_diagram_from_file

# Load C++ HHG simulation library
_lib_path = os.path.join(os.path.dirname(__file__), "..", "libhhg.so")
hhg_lib = ctypes.CDLL(os.path.abspath(_lib_path))

# Define prototype for singleAtomResponse
hhg_lib.requestStop.argtypes = []
hhg_lib.requestStop.restype = None

hhg_lib.singleAtomResponse.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_int
]
hhg_lib.singleAtomResponse.restype = None


# Load convergence library
_convergence_lib_path = os.path.join(os.path.dirname(__file__), "..", "libconvergence.so")
convergence_lib = ctypes.CDLL(os.path.abspath(_convergence_lib_path))

convergence_lib.runConvergence.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
convergence_lib.runConvergence.restype = None

convergence_lib.requestStop.argtypes = []
convergence_lib.requestStop.restype = None


def run_simulation(wavelength, waist, intensity, pulse_length, exponent, generate_plots=True):
    """
    Run the HHG simulation for a single atom using the C++ backend.

    This function calls the `singleAtomResponse` function from `libhhg.so` with
    specified parameters. Upon success, optional time and spectrum plots are generated.

    Args:
        wavelength (float): Laser wavelength in nm.
        waist (float): Beam waist (FWHM) in μm.
        intensity (float): Peak intensity in 10^14 W/cm².
        pulse_length (float): Pulse duration (FWHM) in fs.
        exponent (int): Time grid resolution exponent (TSIZE = 2^exponent).
        generate_plots (bool): If True, generates D(t) and FFT plots automatically.

    Side Effects:
        - Calls compiled C++ backend for numerical simulation
        - Optionally generates and saves plots in the current working directory
        - Prints diagnostics to console
    """
    
    print("▶ Running HHG simulation with:")
    print(f"  wavelength = {wavelength} nm")
    print(f"  waist      = {waist} µm")
    print(f"  intensity  = {intensity} x10^14 W/cm²")
    print(f"  pulse T    = {pulse_length} fs")
    print(f"  exponent   = {exponent} → TSIZE = 2^{exponent}")

    # Remove old success marker if exists
    try:
        os.remove("output_success.txt")
    except FileNotFoundError:
        pass

    # Run backend simulation
    hhg_lib.singleAtomResponse(wavelength, waist, intensity, pulse_length, exponent)

    print("✔ C++ simulation finished.")
    
    #with open("output_success.txt", "w") as f:
    #    f.write("ok")

    if generate_plots:
        plot_time_series_from_file("D", exponent, ylabel="Dipole moment", title="Dipole moment vs Time")
        plot_time_series_from_file("FFT", exponent, ylabel="Spectral intensity", title="Dipole spectrum")


    print("✔ Simulation complete.")


def run_convergence_test(wavelength, waist, intensity, pulse_length, exponent, generate_plots=True):
    """
    Run convergence analysis across different TSIZE values.

    This function invokes the C++ `runConvergence()` function from `libconvergence.so`
    to compute D(t) at various resolutions. It optionally visualizes results.

    Args:
        wavelength (float): Laser wavelength in nm.
        waist (float): Beam waist (FWHM) in μm.
        intensity (float): Peak intensity in 10^14 W/cm².
        pulse_length (float): Pulse duration in fs.
        exponent (int): Used to locate output files for plotting.
        generate_plots (bool): If True, generates D(t), FFT and convergence ratio plots.

    Side Effects:
        - Executes convergence tests via C++ backend
        - Generates output plots from precomputed data
    """
    
    try:
        os.remove("output_success.txt")
    except FileNotFoundError:
        pass
        
    print("▶ Running HHG convergence test...")
    
    convergence_lib.runConvergence(wavelength, waist, intensity, pulse_length)
    
    print("✔ C++ convergence test finished.")
        
    #with open("output_success.txt", "w") as f:
    #    f.write("ok")
        
    if generate_plots:
        plot_time_series_from_file("D", exponent, ylabel="Dipole moment", title="Dipole moment vs Time")
        plot_time_series_from_file("FFT", exponent, ylabel="Spectral intensity", title="Dipole spectrum")
        
    plot_convergence_diagram_from_file("convergence_ratios", exponent)
