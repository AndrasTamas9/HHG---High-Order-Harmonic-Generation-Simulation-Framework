# HHG---High-Order-Harmonic-Generation-Simulation-Framework

# About the Project

This software simulates high-harmonic generation using numerical methods.
It consists of:
- a C++ core simulation engine,
- a Java-based GUI frontend,
- and a Python integration/management layer.

# Introduction

High-Harmonic Generation (HHG) is a nonlinear optical process that occurs when intense, ultrashort laser pulses interact with atoms or molecules, resulting in the emission of coherent radiation at integer multiples of the driving laser frequency. The simulation of HHG phenomena requires highly accurate numerical models capable of resolving both temporal and spatial variations of the electromagnetic field, as well as the nonlinear atomic response. Due to the inherent complexity of the problem, particularly when modeling at the single-atom level with time-dependent quantum dynamics, a computationally efficient and physically robust simulation framework is required.

The purpose of this project is to develop a modular and user-friendly software platform capable of simulating the nonlinear dipole response of a single atom under a focused laser beam using the Lewenstein model. The system is intended to serve both educational and research purposes, providing intuitive parameter control through a graphical user interface (GUI), while relying on a high-performance C++ backend to carry out the numerically intensive computations. Additionally, the platform should support convergence analysis of the time grid resolution and offer visual feedback in the form of dipole and spectral plots, enabling the user to assess the accuracy and physical validity of the simulation.

# User Documentation

## Installation and Setup

### Repository Download

Download the complete project repository from:

```
https://github.com/AndrasTamas9/HHG---High-Order-Harmonic-Generation-Simulation-Framework
```

This includes all necessary source files, scripts, and documentation.

### Required Dependencies

#### Software requirements:

- Java Runtime Environment (JRE) and Java Development Kit (JDK)
- Python 3
- C++ compiler with OpenMP support (e.g., `g++`)
- GNU Scientific Library (GSL)

#### Python packages:

```bash
pip install jpype1 matplotlib pandas
```

#### System dependencies (Debian/Ubuntu example):

```bash
sudo apt install default-jdk libgsl-dev g++
```

### Compilation Instructions

Run all commands from the **root project directory**, where `hhg_pinGUIn.py` is located.

#### Compile Java GUI:

```bash
javac -d out src/gui/*.java
jar cf gui.jar -C out .
```

#### Compile C++ libraries and executables:

```bash
# Shared libraries used by Python
g++ -fPIC -shared -o libhhg.so hhg.cpp -lgsl -lgslcblas -lm -fopenmp
g++ -fPIC -shared -o libconvergence.so convergence.cpp -lgsl -lgslcblas -lm -fopenmp

# Executables for direct CLI use
g++ -DBUILD_EXECUTABLE -o hhg hhg.cpp -lgsl -lgslcblas -lm -fopenmp
g++ -DBUILD_EXECUTABLE -o convergence convergence.cpp -lgsl -lgslcblas -lm -fopenmp
```

### Running the Application

To launch the GUI:

```bash
python3 hhg_pinGUIn.py
```

---

## Using the Graphical Interface

After launch, a Java Swing-based graphical window appears, allowing users to set physical parameters of the laser pulse:

- Wavelength (nm)
- Beam waist (μm)
- Peak intensity ($10^{14}$ W/cm²)
- Pulse duration (fs)
- Temporal resolution exponent: $TSIZE = 2^n$

Click **"Run"** to start the simulation.

Modes:
- **Standard Mode**: runs a single HHG simulation.
- **Convergence Mode**: runs simulations at increasing resolutions and computes convergence.

---

## Simulation Output

Results are stored in the `results/` folder and include:

- Raw data files: `E(t)`, `A(t)`, `D(t)`, and $|D(ω)|$
- Plots:
  - Dipole moment $D(t)$ in time domain
  - Spectrum $|D(ω)|$ of emitted harmonics
  - Convergence plot (in convergence mode)

Filenames encode the resolution exponent, e.g. `D_13.txt`.

---

## Understanding the Results

- The $D(t)$ plot shows electron dynamics.
- The $|D(ω)|$ plot reveals the harmonic spectrum, plateau, and cutoff.
- The convergence plot shows the area error between resolutions on a log scale.

Important modeling note:

- The underlying model (Lewenstein integral / SFA) is only valid in the **tunneling regime** ($\gamma < 1$).
- Unphysical or invalid results may arise for:
  - Extremely short wavelengths
  - Very low intensities
- Always ensure input parameters lie within the model's domain of validity.

---

## Advanced Features and Termination

- A **Stop** button allows graceful simulation termination.
- Parameter sets can be **saved/loaded**.
- Optional sound feedback is available on simulation start.

The backend handles interrupts safely and prevents data corruption.

---

## Troubleshooting

If simulation output is missing:

- Check Python and Java console output.
- Possible causes:
  - Missing or mislinked shared libraries
  - Invalid input values (e.g., zero intensity)
  - Output files not generated → try cleaning and re-running
