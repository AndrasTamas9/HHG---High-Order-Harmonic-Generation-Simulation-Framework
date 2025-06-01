"""
Plotting utilities for HHG simulation results.

This module provides:
- Line plotting of time-domain (D) or frequency-domain (FFT) output
- Convergence diagram visualization from tabulated data
- Cleanup of result plots and text files

All data files are expected to be in the `../results/` directory relative to this script.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def plot_time_series_from_file(file_prefix: str, exponent: int, target_coords=(0.0, 0.0), ylabel="Value", title="Time Series"):
    """
    Loads and plots a 1D time or frequency series from a structured text file.

    The file is expected to be located at:
        ../results/<file_prefix>_<exponent>.txt

    Alternating lines are expected:
        - One line with (r, z) coordinates
        - One line with the corresponding time series values

    Args:
        file_prefix (str): File name prefix (e.g., 'D', 'FFT')
        exponent (int): Resolution level (TSIZE = 2^exponent)
        target_coords (tuple): (r, z) coordinates to filter the block
        ylabel (str): Label for the Y-axis
        title (str): Plot title

    Side Effects:
        - Saves a plot image: ../results/<file_prefix>_plot.png
        - Prints messages or warnings if file or data is not found
    """

    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", f"{file_prefix}_{exponent}.txt"))
    
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return

    with open(filename, "r") as f:
        lines = f.readlines()

    series = None

    i = 0
    while i < len(lines) - 1:
        coords = lines[i].strip().split()
        if len(coords) != 2:
            i += 1
            continue

        r, z = float(coords[0]), float(coords[1])
        if r == target_coords[0] and z == target_coords[1]:
            # Found matching coordinates; extract data line
            data_line = lines[i + 1]
            series = [float(val) for val in data_line.strip().split()]
            break
        i += 2

    if series is None:
        print(f"No data found at coordinates {target_coords} in {filename}")
        return

    # Plot
    plt.figure(figsize=(10, 4))
    plt.plot(series)
    plt.xlabel("Time index" if "D_" in file_prefix else "Frequency index")
    plt.ylabel(ylabel)
    plt.title(f"{title} @ r={target_coords[0]}, z={target_coords[1]} (exponent={exponent})")
    plt.tight_layout()
    output_file = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "results", f"{file_prefix.lower()}_plot.png"))
    plt.savefig(output_file)
    plt.close()

    print(f"Plot saved to: {output_file}")
    
    
def plot_convergence_diagram_from_file(file_prefix: str, exponent: int, target_coords=(0.0, 0.0)):
    """
    Reads convergence ratio data and plots area_ratio vs. exponent (fine_n).

    Expects a whitespace-separated file at:
        ../results/<file_prefix>.txt

    Each line should include: fine_n, coarse_n, r, z, area_ratio

    Args:
        file_prefix (str): File name prefix (typically 'convergence_ratios')
        exponent (int): Marks the vertical line at this exponent in the plot
        target_coords (tuple): (r, z) position for filtering convergence data

    Side Effects:
        - Generates convergence_plot.png
        - Displays warning if exponent not present
    """
    
    filename = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", f"{file_prefix}.txt"))

    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return

    df = pd.read_csv(filename, sep=r'\s+', comment='#',
                     names=['fine_n', 'coarse_n', 'r', 'z', 'area_ratio'])

    filtered = df[(df['r'] == target_coords[0]) & (df['z'] == target_coords[1])]
    plot_data = filtered[['fine_n', 'area_ratio']].drop_duplicates().sort_values(by='fine_n')
    
    plt.figure(figsize=(8, 5))
    plt.plot(plot_data['fine_n'], plot_data['area_ratio'], marker='o')

    n0 = exponent
    if n0 in plot_data['fine_n'].values:
        val = plot_data[plot_data['fine_n'] == n0]['area_ratio'].values[0]
        plt.axvline(x=n0, color='red', linestyle='--', label=f'n={n0}')
        plt.text(n0 + 0.1, val, f'{val:.2e}', color='red', verticalalignment='bottom')
    else:
        print(f"n={n0} was not found among the data at (r={target_coords[0]}, z={target_coords[1]})")

    plt.xlabel('n - exponent')
    plt.ylabel(f'area_ratio @ (r={target_coords[0]}, z={target_coords[1]})')
    plt.title('Convergence Test')
    plt.yscale('log')
    plt.grid(True, which='both', linestyle='--')
    plt.legend()
    plt.tight_layout()

    output_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results", "convergence_plot.png"))
    plt.savefig(output_file)
    plt.close()

    print(f"Plot saved to: {output_file}")
    
    
def cleanup_result_images():
    """
    Deletes all PNG plot files from the results/ directory.

    This is typically used to clean up previous simulation output before
    running a new simulation or when closing the GUI.

    Looks for files matching: results/*.png
    """
    
    print("Cleaning up .png images in 'results/' directory...")
    
    results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))
    png_files = glob.glob(os.path.join(results_dir, "*.png"))

    for file in png_files:
        try:
            os.remove(file)
            print(f"Deleted: {file}")
        except Exception as e:
            print(f"Could not delete {file}: {e}")


def cleanup_result_files():
    """
    Deletes all .png and .txt files from the results/ directory.

    This is typically used to clean up previous simulation output before
    running a new simulation or when closing the GUI.

    Looks for files matching: results/*.png and results/*.txt
    """
    
    print("Cleaning up .png and .txt files in 'results/' directory...")
    
    results_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "results"))
    
    # Kiterjesztések listája
    extensions = ["*.png", "*.txt"]
    
    # Fájlok összegyűjtése
    files_to_delete = []
    for ext in extensions:
        files_to_delete.extend(glob.glob(os.path.join(results_dir, ext)))

    # Fájlok törlése
    for file in files_to_delete:
        try:
            os.remove(file)
            print(f"Deleted: {file}")
        except Exception as e:
            print(f"Could not delete {file}: {e}")
