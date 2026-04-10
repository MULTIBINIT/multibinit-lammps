# Multibinit-LAMMPS Tutorial

This tutorial provides a comprehensive guide for installing and using the Multibinit-LAMMPS interface, which enables LAMMPS molecular dynamics simulations using ABINIT's Multibinit effective potential.

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Building ABINIT with Shared Library](#2-building-abinit-with-shared-library)
3. [Building LAMMPS](#3-building-lammps)
4. [Building the Multibinit-LAMMPS Plugin](#4-building-the-multibinit-lammps-plugin)
5. [Basic Usage](#5-basic-usage)
6. [Input File Reference](#6-input-file-reference)
7. [Example: BaHfO3 Simulation](#7-example-bahfo3-simulation)
8. [Troubleshooting](#8-troubleshooting)
9. [Advanced Topics](#9-advanced-topics)

---

## 1. Prerequisites

### 1.1 Required Software

| Software | Version | Purpose |
|----------|---------|---------|
| CMake | ≥ 3.20 | Build system |
| C++ compiler | C++17 compatible | Compile plugin |
| Fortran compiler | gfortran ≥ 9 | Compile ABINIT |
| MPI | OpenMPI or MPICH | Parallel execution |
| BLAS/LAPACK | Any | Linear algebra |
| NetCDF | ≥ 4.0 | I/O library |
| LibXC | ≥ 4.0 | Exchange-correlation |

### 1.2 macOS Installation (Homebrew)

```bash
# Install dependencies
brew install cmake gcc open-mpi netcdf libxc

# Set environment for gfortran (adjust version as needed)
export FC=gfortran-14
export CC=gcc-14
export CXX=g++-14
```

### 1.3 Linux Installation (Ubuntu/Debian)

```bash
# Install dependencies
sudo apt-get install cmake gfortran openmpi-bin libopenmpi-dev \
    libnetcdf-dev libnetcdff-dev libblas-dev liblapack-dev \
    libxc-dev

# Set environment
export FC=gfortran
export CC=gcc
export CXX=g++
```

### 1.4 Directory Structure

The project expects the following directory structure:

```
multibinit_lammps/
├── abinit/           # ABINIT source code (git repository)
├── lammps/           # LAMMPS source code (git repository)
├── mblammps/         # Plugin source code
│   ├── pair_multibinit.h
│   ├── pair_multibinit.cpp
│   ├── multibinitplugin.cpp
│   └── CMakeLists.txt
└── examples/         # Example files
```

---

## 2. Building ABINIT with Shared Library

The Multibinit-LAMMPS plugin requires ABINIT to be built as a shared library with C bindings enabled.

### 2.1 Clone ABINIT (if not already done)

```bash
cd multibinit_lammps
git clone https://github.com/abinit/abinit.git
cd abinit
```

### 2.2 Configure ABINIT

Create a build directory and configure:

```bash
mkdir build_so && cd build_so

# Configure with shared library support
cmake .. \
    -DCMAKE_INSTALL_PREFIX=$PWD/install \
    -DENABLE_SHARED_LIBS=ON \
    -DENABLE_MPI=ON \
    -DENABLE_NETCDF=ON \
    -DENABLE_XC=ON \
    -DCMAKE_BUILD_TYPE=Release
```

**Important CMake options:**

| Option | Value | Description |
|--------|-------|-------------|
| `ENABLE_SHARED_LIBS` | ON | Build shared library (required) |
| `ENABLE_MPI` | ON | Enable MPI support |
| `ENABLE_NETCDF` | ON | Enable NetCDF I/O |
| `ENABLE_XC` | ON | Enable LibXC for exchange-correlation |

### 2.3 Build and Install

```bash
# Build (use multiple cores for faster compilation)
make -j8

# Install to local directory
make install
```

### 2.4 Verify Build

Check that the shared library was created:

```bash
# macOS
ls install/lib/libabinit*.dylib

# Linux
ls install/lib/libabinit*.so
```

Check that the module files exist:

```bash
ls modules/m_mb_cbinding.mod
# Should exist for C bindings
```

---

## 3. Building LAMMPS

### 3.1 Clone LAMMPS (if not already done)

```bash
cd multibinit_lammps
git clone https://github.com/lammps/lammps.git
cd lammps
```

### 3.2 Configure LAMMPS

```bash
mkdir build && cd build

cmake ../cmake \
    -DPKG_PLUGIN=ON \
    -DLAMMPS_EXCEPTIONS=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_BUILD_TYPE=Release
```

**Important CMake options:**

| Option | Value | Description |
|--------|-------|-------------|
| `PKG_PLUGIN` | ON | Enable plugin support (required) |
| `LAMMPS_EXCEPTIONS` | ON | Better error handling |
| `BUILD_SHARED_LIBS` | ON | Build shared library |

### 3.3 Build LAMMPS

```bash
make -j8
```

### 3.4 Verify Build

```bash
./lmp -h
# Should show LAMMPS version and help information
```

---

## 4. Building the Multibinit-LAMMPS Plugin

### 4.1 Configure the Plugin

```bash
cd multibinit_lammps/mblammps
mkdir build && cd build

cmake .. \
    -DABINIT_ROOT=../../abinit/build_so
```

The CMake output should show:

```
-- Found ABINIT includes: .../build_so/modules
-- Found abinit: .../build_so/lib/libabinit.dylib
-- Found netcdf: /path/to/libnetcdf.dylib
-- Found xc: /path/to/libxc.dylib

-- Multibinit LAMMPS Plugin Configuration Summary:
--   LAMMPS source: .../lammps/src
--   ABINIT root: .../abinit/build_so
--   ABINIT includes: .../build_so/modules
--   ABINIT libs: .../libabinit.dylib;...
```

### 4.2 Build the Plugin

```bash
make -j4
```

### 4.3 Verify Build

```bash
ls multibinitplugin.so
# Should exist and be a shared library
```

### 4.4 Set Library Path (if needed)

If you get library loading errors, set the library path:

```bash
# macOS
export DYLD_LIBRARY_PATH=/path/to/abinit/build_so/lib:$DYLD_LIBRARY_PATH

# Linux
export LD_LIBRARY_PATH=/path/to/abinit/build_so/lib:$LD_LIBRARY_PATH
```

---

## 5. Basic Usage

### 5.1 Loading the Plugin

In your LAMMPS input script:

```lammps
# Load the Multibinit plugin
plugin load /path/to/multibinitplugin.so
```

### 5.2 Pair Style Syntax

```lammps
pair_style multibinit nx ny nz keyword value ...

pair_coeff * * system_file [coefficients_file]
```

**Required parameters:**
- `nx ny nz`: Supercell dimensions (number of unit cells in each direction)

**Optional keywords:**

| Keyword | Values | Default | Description |
|---------|--------|---------|-------------|
| `ngqpt` | qx qy qz | 2 2 2 | Q-point grid for dipole-dipole |
| `dipdip` | 0 or 1 | 1 | Enable dipole-dipole interaction |

**File arguments:**
- `system_file`: DDB or XML file containing system definition (required)
- `coefficients_file`: XML file with polynomial coefficients (optional)

### 5.3 Unit System

The plugin uses LAMMPS **metal** units:

| Quantity | Unit |
|----------|------|
| Energy | eV |
| Distance | Angstrom |
| Force | eV/Å |
| Time | ps |
| Pressure | bars |

**Unit Conversions (from atomic units):**

| From | To | Conversion Factor |
|------|-----|--------------------|
| Hartree (Ha) | eV | 27.211386 |
| Bohr | Å | 0.529177 |
| Ha/Bohr | eV/Å | 51.421 |
| Ha/Bohr³ | GPa | 2.9421×10⁴ |

---

## 6. Input File Reference

### 6.1 Complete Input Example

```lammps
# Multibinit-LAMMPS Example
# Perovskite oxide simulation

# ========== System Setup ==========
units           metal
dimension       3
boundary        p p p
atom_style      atomic

# ========== Structure Definition ==========
# Define lattice (must match DDB file)
variable a equal 4.0  # Lattice constant in Angstrom

lattice custom ${a} &
    a1 1.0 0.0 0.0 &
    a2 0.0 1.0 0.0 &
    a3 0.0 0.0 1.0 &
    basis 0.0 0.0 0.0 &    # A-site atom
    basis 0.5 0.5 0.5 &    # B-site atom
    basis 0.5 0.5 0.0 &    # Oxygen
    basis 0.0 0.5 0.5 &    # Oxygen
    basis 0.5 0.0 0.5      # Oxygen

# Create supercell
region box block 0 2 0 2 0 2
create_box 3 box
create_atoms 1 box &
    basis 1 1 &
    basis 2 2 &
    basis 3 3 &
    basis 4 3 &
    basis 5 3

# Set masses
mass 1 137.33  # A-site (e.g., Ba)
mass 2 178.49  # B-site (e.g., Hf)
mass 3 16.00   # Oxygen

# ========== Plugin and Potential ==========
plugin load ../build/multibinitplugin.so

pair_style multibinit 2 2 2 ngqpt 2 2 2 dipdip 1
pair_coeff * * system.DDB coefficients.xml

# ========== Output Settings ==========
thermo 100
thermo_style custom step temp pe ke etotal press vol

# ========== Dynamics ==========
timestep 0.001
fix 1 all nve

run 1000
```

### 6.2 Structure Definition Guidelines

The LAMMPS structure must match what Multibinit expects from the DDB file:

1. **Atom count**: Supercell must have `nx × ny × nz × natom_unit` atoms
2. **Lattice vectors**: Must match the supercell dimensions
3. **Atom positions**: Must be in the same order as Multibinit expects

For a 2×2×2 supercell of a 5-atom perovskite:
- Total atoms = 2 × 2 × 2 × 5 = 40 atoms

### 6.3 Matching DDB File

Check your DDB file for:

```
natom         5           # Atoms per unit cell
acell  0.78411195940000D+01  ...  # Lattice in Bohr
```

Convert `acell` to Angstrom: `7.841 × 0.529177 = 4.149 Å`

---

## 7. Example: BaHfO3 Simulation

### 7.1 Required Files

You need the following files from ABINIT calculations:

1. **DDB file** (`BaHfO3.DDB`): Contains force constants, Born charges, dielectric tensor
2. **XML coefficients file** (`BaHfO3_coeffs.xml`): Contains polynomial coefficients for anharmonic terms

These can be generated using the ABINIT/Multibinit workflow.

### 7.2 Input File

Create `in.bahfo3`:

```lammps
# BaHfO3 molecular dynamics with Multibinit
# 2x2x2 supercell (40 atoms)

units metal
dimension 3
boundary p p p
atom_style atomic

# BaHfO3 perovskite primitive cell
# Lattice: 7.8411 Bohr = 4.1490 Angstrom
# Fractional positions:
#   Ba  (0,0,0)
#   Hf  (0.5,0.5,0.5)
#   O   (0.5,0.5,0), (0,0.5,0.5), (0.5,0,0.5)

variable a equal 4.1490

lattice custom ${a} &
    a1 1.0 0.0 0.0 &
    a2 0.0 1.0 0.0 &
    a3 0.0 0.0 1.0 &
    basis 0.0 0.0 0.0 &
    basis 0.5 0.5 0.5 &
    basis 0.5 0.5 0.0 &
    basis 0.0 0.5 0.5 &
    basis 0.5 0.0 0.5

region box block 0 2 0 2 0 2
create_box 3 box
create_atoms 1 box &
    basis 1 1 &
    basis 2 2 &
    basis 3 3 &
    basis 4 3 &
    basis 5 3

mass 1 137.33  # Ba
mass 2 178.49  # Hf
mass 3 16.00   # O

# Load plugin and define potential
plugin load ../build/multibinitplugin.so

pair_style multibinit 2 2 2 ngqpt 2 2 2 dipdip 1
pair_coeff * * BaHfO3.DDB BaHfO3_coeffs.xml

# Output
thermo 10
thermo_style custom step temp pe ke etotal press vol atoms

# Equilibrate at 300K
timestep 0.001
velocity all create 300.0 12345
fix 1 all nvt temp 300.0 300.0 0.1

run 1000
```

### 7.3 Running the Simulation

```bash
# Single processor
/path/to/lmp -in in.bahfo3

# Parallel (4 MPI processes)
mpirun -np 4 /path/to/lmp -in in.bahfo3
```

---

## 8. Troubleshooting

### 8.1 Plugin Not Found

**Error:**
```
Open of file multibinitplugin.so failed: dlopen(...): symbol not found
```

**Solution:**
The plugin cannot find required libraries. Ensure:

1. ABINIT shared library is in library path:
   ```bash
   # macOS
   export DYLD_LIBRARY_PATH=/path/to/abinit/build_so/lib:$DYLD_LIBRARY_PATH

   # Linux
   export LD_LIBRARY_PATH=/path/to/abinit/build_so/lib:$LD_LIBRARY_PATH
   ```

2. NetCDF and LibXC are installed and accessible:
   ```bash
   # Check if libraries are found
   otool -L multibinitplugin.so  # macOS
   ldd multibinitplugin.so       # Linux
   ```

### 8.2 Atom Count Mismatch

**Error:**
```
Atom count mismatch: LAMMPS has X, Multibinit expects Y
```

**Solution:**
The supercell dimensions don't match. Verify:

1. `pair_style multibinit nx ny nz` matches your supercell
2. Total atoms = nx × ny × nz × natom_unit
3. The DDB file's primitive cell has the correct number of atoms

**Check DDB file:**
```bash
grep "natom" your_file.DDB
# natom should match atoms per unit cell
```

### 8.3 Unrecognized Pair Style

**Error:**
```
ERROR: Unrecognized pair style 'multibinit'
```

**Solution:**
The plugin wasn't loaded successfully. Check:

1. Plugin file exists and is readable
2. All dependencies are available
3. Check for error messages during plugin load:
   ```lammps
   plugin load ../build/multibinitplugin.so
   # Look for "Loading plugin:" message
   ```

### 8.4 High or Unexpected Energy Values

**Symptom:**
Energy is extremely high (e.g., 10⁸ eV) or NaN instead of expected values.

**Possible causes:**
1. Structure mismatch between LAMMPS and DDB
2. Incorrect lattice constant conversion
3. Using wrong DDB/coefficient files
4. Unstable phonon modes in the structure

**Solution:**
1. Verify lattice constant matches DDB file (convert Bohr to Å)
2. Check atom positions are in correct order
3. Ensure using metal units (Angstrom, eV)
4. Try a known stable structure

### 8.5 NaN Energy

**Symptom:**
Energy returns NaN.

**Possible cause:**
The structure has unstable phonon modes (e.g., cubic perovskite BaTiO3 is dynamically unstable).

**Solution:**
Use a stable structure or a different material. Check phonon spectrum with ABINIT first.

---

## 9. Advanced Topics

### 9.1 MPI Parallelization

The plugin supports MPI parallelization through ABINIT's internal MPI:

```bash
# Run with 8 MPI processes
mpirun -np 8 lmp -in input.lammps
```

The plugin automatically uses the LAMMPS MPI communicator.

### 9.2 Generating DDB Files

DDB files are generated from ABINIT DFPT calculations:

1. **Ground state calculation**: Run SCF calculation
2. **DFPT calculation**: Compute phonons, Born charges, dielectric tensor
3. **Merge DDB**: Combine partial DDB files using `mrgddb`

Example ABINIT input for DFPT:

```
# Ground state
tolwfr 1.0d-18
nstep 100
ecut 50.0
...

# DFPT
rfphonon 1
rfdir 1 1 1
nqpt 1
qpt 0.0 0.0 0.0
```

### 9.3 Generating Coefficient Files

Coefficient XML files are generated by Multibinit fitting:

```bash
multibinit < input.files
```

The fitting process determines polynomial coefficients from training data (typically from DFT calculations on displaced structures).

### 9.4 Stress Tensor

The plugin computes the stress tensor (virial) which can be used for:
- NPT simulations (constant pressure)
- Cell optimization
- Elastic constant calculations

```lammps
# Constant pressure simulation
fix 1 all npt temp 300 300 0.1 iso 0 0 1.0
```

### 9.5 Extracting Forces and Stress

```lammps
# Output forces to file
dump 1 all custom 100 dump.forces id type x y z fx fy fz

# Output stress components
variable pxx equal press
variable pyy equal press[2]
variable pzz equal press[3]
fix 2 all print 100 "${pxx} ${pyy} ${pzz}" file stress.dat
```

### 9.6 Different Supercell Sizes

For different supercell sizes, adjust:

```lammps
# For 3x3x3 supercell
pair_style multibinit 3 3 3 ngqpt 2 2 2 dipdip 1

region box block 0 3 0 3 0 3
# Total atoms = 3 × 3 × 3 × 5 = 135 atoms
```

---

## References

- ABINIT Documentation: https://docs.abinit.org/
- Multibinit Documentation: https://docs.abinit.org/topics/multibinit/
- LAMMPS Documentation: https://docs.lammps.org/
- LAMMPS Plugin System: https://docs.lammps.org/Developer_plugins.html
