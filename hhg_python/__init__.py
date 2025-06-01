"""
hhg_python

Python bindings and plotting tools for high-harmonic generation (HHG) simulations.

This package provides:
- High-level interface to the compiled C++ HHG backend (`run_simulation`, `run_convergence_test`)
- Plotting utilities for inspecting time-domain and frequency-domain results
- Convergence analysis plotting and automated cleanup of results

Modules:
    - `interface.py`: C++ integration via ctypes
    - `plotting.py`: result visualization and file management

Available functions:
- run_simulation(...)
- run_convergence_test(...)
- plot_time_series_from_file(...)
- plot_convergence_diagram_from_file(...)
- cleanup_result_images()
- cleanup_result_files()
"""

from .plotting import plot_time_series_from_file, cleanup_result_images, cleanup_result_files, plot_convergence_diagram_from_file
from .interface import run_simulation, run_convergence_test

__all__ = ["plot_time_series_from_file", "plot_convergence_diagram_from_file", "cleanup_result_images", "cleanup_result_files", "run_simulation", "run_convergence_test"]

