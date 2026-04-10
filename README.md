# mblammps - Multibinit-LAMMPS Interface

This package provides a LAMMPS pair style for the Multibinit effective potential
from ABINIT.

## Files

- `pair_multibinit.h` - LAMMPS pair style header
- `pair_multibinit.cpp` - LAMMPS pair style implementation
- `CMakeLists.txt` - CMake build configuration
- `Install.sh` - LAMMPS package installation script

## Requirements

- LAMMPS (tested with version >= 2023)
- ABINIT with Multibinit enabled
- Fortran compiler (for ABINIT library)

## Installation

### Method 1: CMake (Recommended)

1. Build and install ABINIT with Multibinit enabled:
   ```bash
   cd abinit
   mkdir build && cd build
   cmake .. -Denable_multibinit=ON -DCMAKE_INSTALL_PREFIX=/path/to/install
   make install
   ```

2. Build LAMMPS with the Multibinit package:
   ```bash
   cd lammps
   mkdir build && cd build
   cmake ../cmake \
     -DPKG_USER-MULTIBINIT=ON \
     -DABINIT_ROOT=/path/to/abinit/install
   make -j
   ```

### Method 2: Traditional Make

1. Copy files to LAMMPS source directory:
   ```bash
   cp pair_multibinit.* /path/to/lammps/src/USER-MULTIBINIT/
   ```

2. Add to LAMMPS Makefile:
   ```makefile
   # In src/MAKE/Makefile.your_machine
   ABINIT_INC = -I/path/to/abinit/include
   ABINIT_PATH = -L/path/to/abinit/lib
   ABINIT_LIB = -labinit -lgfortran
   ```

3. Build LAMMPS:
   ```bash
   cd /path/to/lammps/src
   make yes-user-multibinit
   make your_machine
   ```

## Usage

### LAMMPS Input Script

```lammps
# Use metal units (eV, Angstrom, ps)
units metal

# Define atoms
lattice fcc 4.0
region box block 0 2 0 2 0 2
create_box 1 box
create_atoms 1 box
mass 1 58.93

# Pair style: multibinit ncell nx ny nz [ngqpt qx qy qz] [dipdip flag]
pair_style multibinit 2 2 2 ngqpt 2 2 2 dipdip 1

# Pair coeff: * * sys_file [coeff_file]
# sys_file: DDB or XML file with system definition
# coeff_file: optional XML file with coefficients
pair_coeff * * BaTiO3_DDB BaTiO3_coeff.xml

# Run MD
fix 1 all nve
run 1000
```

### Input Files

The pair style requires the following input files:

1. **System file** (required): DDB or XML file containing:
   - Crystal structure
   - Force constants
   - Effective charges
   - Dielectric tensor

2. **Coefficient file** (optional): XML file with additional coefficients:
   - Anharmonic terms
   - Higher-order force constants

### Parameters

- `ncell nx ny nz`: Supercell dimensions (default: 1 1 1)
- `ngqpt qx qy qz`: q-point grid for dipole-dipole (default: 2 2 2)
- `dipdip flag`: Enable/disable dipole-dipole (0/1, default: 1)

## Unit System

Multibinit uses atomic units internally (Hartree, Bohr). The pair style
converts to LAMMPS metal units (eV, Angstrom):

- Energy: 1 Ha = 27.211 eV
- Length: 1 Bohr = 0.529 Angstrom
- Force: 1 Ha/Bohr = 51.42 eV/Angstrom

## Current Status

**All energy calculations are working!** The plugin successfully computes energies for both harmonic and dipole-dipole calculations.

### Working Features
- Plugin loading and initialization
- Harmonic energy and force calculations (dipdip=0)
- Long-range dipole-dipole via Ewald summation (dipdip=1)
- MPI communicator passing from LAMMPS to Multibinit

### Known Limitations

1. **Atom Ordering (CRITICAL for MD)**: Multibinit expects atoms in a specific order
   that matches its internal supercell generation. LAMMPS `create_atoms` creates atoms
   in a different order, which causes incorrect forces at the ideal positions.
   
   **Symptom**: Large forces (~100-300 eV/Å) at ideal positions instead of ~0.
   
   **Workaround**: Use the `hist2lammps.py` workflow below.

2. **MPI Parallelization**: Requires all atoms gathered to all processors.
   Full domain decomposition support is planned.

3. **Atom Count**: Must match the supercell size specified by `ncell`.

### Getting Correct Atom Ordering for MD

For molecular dynamics simulations, you must use a LAMMPS data file with the
correct atom ordering. The recommended workflow:

1. **Run Multibinit once** to generate a HIST.nc file:
   ```bash
   multibinit < input.files
   ```

2. **Convert to LAMMPS data file**:
   ```bash
   python scripts/hist2lammps.py output_HIST.nc structure.lmp
   ```

3. **Use in LAMMPS**:
   ```lammps
   read_data structure.lmp
   pair_style multibinit 2 2 2 ngqpt 2 2 2 dipdip 1
   pair_coeff * * BaHfO3.DDB
   fix 1 all nve
   run 1000
   ```

**Note**: Energy calculations (`run 0`) work correctly with `create_atoms`, but
MD simulations require the proper atom ordering from the workflow above.



### Quick Test

```bash
cd mblammps/tests
lmp -in in.bho_test
```

Expected output: Energy ~-29241 eV for 2x2x2 BaHfO3 supercell (40 atoms).

### Build and Test from Source

1. **Build ABINIT with shared library**:
   ```bash
   cd abinit/build
   cmake .. -Denable_multibinit=ON -Denable_shared=ON -DCMAKE_INSTALL_PREFIX=$PWD/../install
   make -j4 install
   ```

2. **Build LAMMPS with plugin support**:
   ```bash
   cd lammps/build
   cmake ../cmake -DPKG_PLUGIN=ON -DBUILD_SHARED_LIBS=ON
   make -j4
   ```

3. **Build the Multibinit plugin**:
   ```bash
   cd mblammps/build
   cmake .. -DABINIT_ROOT=$PWD/../../abinit/install
   make -j4
   ```

4. **Run test**:
   ```bash
   cd mblammps/tests
   lmp -in in.bho_test
   ```

## References

- ABINIT: https://www.abinit.org
- Multibinit: https://docs.abinit.org/topics/multibinit/
- LAMMPS: https://www.lammps.org

## License

This package is distributed under the same license as ABINIT (GPL v3+).
