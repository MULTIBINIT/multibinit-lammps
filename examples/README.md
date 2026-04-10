# Multibinit-LAMMPS Examples

This directory contains example input files for using the Multibinit-LAMMPS plugin.

## File Structure

```
examples/
├── README.md              # This file
├── BaHfO3.DDB            # Harmonic force constants (required)
├── BaHfO3.xml            # Anharmonic coefficients (optional)
├── 01_basic_nve/         # Basic NVE example
├── 02_nvt_thermostat/    # NVT with thermostat
├── 03_npt_barostat/      # NPT with pressure control
└── 04_parallel_mpi/      # MPI parallel example
```

## Quick Reference

| Example | Input File | Purpose | Required Files |
|---------|-----------|---------|---------------|
| Basic NVE | `in.bho_test` | Energy conservation test | BaHfO3.DDB |
| NVE with XML | `in.bho_md` | Full potential (harmonic + anharmonic) | BaHfO3.DDB, BaHfO3.xml |
| NVT | `in.nvt_test` | Temperature control at 300K | BaHfO3.DDB, BaHfO3.xml |
| NVE (minimal) | `in.nve_test` | Minimal energy check | BaHfO3.DDB, BaHfO3.xml |
| NPT | `in.npt_harmonic` | Pressure/volume control | BaHfO3.DDB |

## Running Examples

### Basic NVE (Energy Conservation)

```bash
lmp -in in.bho_test
```

**What to check**: Total energy (etotal) should be constant throughout the simulation.

Expected output:
```
Step Temp PotEng KinEng TotEng Press Volume
   0   300 -29240.002 1.512 -29238.49 2826.4 571.37
 100   116 -29239.076 0.587 -29238.49 1092.7 571.37
 200   151 -29239.252 0.763 -29238.49 1421.8 571.37
```

### NVT (Temperature Control)

```bash
lmp -in in.nvt_test
```

**What to check**: Temperature should stabilize around 300K after initial equilibration.

Expected output:
```
Step Temp PotEng KinEng TotEng Press Volume
   0   100 -29240.002 0.504 -29239.50 941.6 571.37
 500   216 -29238.561 1.093 -29237.47 2031.6 571.37
1000   214 -29238.488 1.082 -29237.41 2012.9 571.37
3000   306 -29238.040 1.546 -29236.49 2876.4 571.37
5000   349 -29238.466 1.762 -29236.70 3286.1 571.37
```

Note: Temperature fluctuates around target (300K) - this is normal for a small 40-atom system.

### NPT (Pressure Control)

```bash
lmp -in in.npt_harmonic
```

**What to check**: Pressure fluctuates around target (0), volume adjusts to equilibrium.

Expected output:
```
Step Temp PotEng KinEng TotEng Press Volume
   0   300 -29240.002 1.512 -29238.49 3251.2 571.37
 500   297 -29238.217 1.498 -29236.73 -804.2 571.61
1000   343 -29238.664 1.730 -29236.93 404.4 571.62
1500   279 -29238.634 1.407 -29237.23 453.9 571.39
2000   308 -29238.064 1.553 -29236.51 -504.4 571.50
```

Note: Volume changes slightly as system equilibrates to target pressure.

### MPI Parallel Run

```bash
mpirun -np 4 lmp -in in.bho_xml
```

**What to check**: Same results as serial run, but faster. Load balancing shows in timing breakdown.

## Example Descriptions

### in.bho_test - Basic NVE

Simplest possible test - harmonic potential only, NVE ensemble.
- Verifies energy conservation
- Tests force calculation accuracy
- Minimal dependencies

### in.bho_md - NVE with Anharmonic Terms

Includes XML coefficients for anharmonic contributions.
- Tests anharmonic force evaluation
- Slower than harmonic-only (more terms to compute)
- Good benchmark for accuracy vs speed

### in.bho_xml - MPI-Ready Input

Designed for parallel execution with XML coefficients.
- Includes all potential terms
- Works with `mpirun -np N`
- Demonstrates full capabilities

### in.nvt_test - NVT Thermostat

Canonical ensemble simulation at constant temperature.
- Nosé-Hoover thermostat
- Target: 300K
- Good for equilibrating structures

### in.nve_test - Energy Conservation Check

Microcanonical ensemble - no thermostat or barostat.
- Total energy must be conserved
- Temperature will fluctuate
- Critical validation test

### in.npt_harmonic - NPT Barostat

Isothermal-isobaric ensemble with pressure control.
- Target: 300K, 0 pressure
- Isotropic barostat
- Volume changes to reach equilibrium

## Common Patterns

### Pattern 1: Equilibration followed by Production

```lammps
# 1. Initialize velocities
velocity all create 300.0 12345

# 2. Equilibration (NVT)
fix 1 all nvt temp 300 300 0.1
run 5000

# 3. Switch to NVE for production
unfix 1
fix 2 all nve
run 10000
```

### Pattern 2: Temperature Ramp

```lammps
# Heat from 300K to 1000K
fix 1 all nvt temp 300 1000 0.1
run 10000
```

### Pattern 3: Pressure Scaling

```lammps
# Compress from 0 to 10 kbar
fix 1 all npt temp 300 300 0.1 iso 0 100 10.0
run 5000
```

## Troubleshooting

### "Atom count mismatch"

**Solution**: Ensure your structure matches the DDB supercell dimensions.

### "Pressure explosion in NPT"

**Solution**: Use harmonic potential only (no XML) for NPT, or increase `Pdamp` parameter.

### "Energy drift in NVE"

**Solution**: Check timestep is appropriate (1 fs = 0.001 in metal units). Try smaller timestep.

## Performance Tips

1. **Use MPI**: `mpirun -np 4 lmp -in input` for 4x speedup
2. **Output frequency**: Use `thermo 100` not `thermo 1` to reduce I/O
3. **Timestep**: 1 fs is typical, 2 fs may work for harmonic-only
4. **Thermostat damping**: 0.1 ps (100 fs) is gentle, 0.01 is stiff

## Output Analysis

### Extract Energy
```bash
grep "Loop time" log.lammps
grep -A1 "Step" log.lammps | tail -n +2 > energy.dat
```

### Plot Temperature
```python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('log.lammps', skiprows=2)
step, temp = data[:,0], data[:,1]
plt.plot(step, temp)
plt.xlabel('Step')
plt.ylabel('Temperature (K)')
plt.savefig('temperature.png')
```

## References

- LAMMPS documentation: https://docs.lammps.org/
- ABINIT Multibinit: https://docs.abinit.org/guide/multibinit/
- Examples in this directory mirror the test cases in `../tests/`

---

## Generating DDB Files

DDB files are generated from ABINIT DFPT calculations:

1. Run SCF calculation with ABINIT
2. Run DFPT to compute phonons, Born charges, dielectric tensor
3. Merge partial DDB files using `mrgddb`

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
qpt 0.0  0.0 0.0
```

Then run ABINIT through these files.

## Generating Coefficient Files

Coefficient XML files are generated by Multibinit fitting

```bash
multibinit < input.files
```

The fitting process determines polynomial coefficients through training data (typically from DFT calculations on displaced structures).
