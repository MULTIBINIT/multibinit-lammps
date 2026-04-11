# mblammps - Multibinit-LAMMPS Interface

This package provides a LAMMPS pair style for the Multibinit effective potential
from ABINIT.


## Requirements

- LAMMPS (tested with version >= 2023)
- ABINIT with Multibinit enabled
- Fortran compiler (for ABINIT library)


### Working Features
- Plugin loading and initialization
- Harmonic energy and force calculations (dipdip=0)
- Long-range dipole-dipole via Ewald summation (dipdip=1)
- Anharmonic interaction
- MPI communicator passing from LAMMPS to Multibinit

### Known Limitations


### Getting Correct Atom Ordering for MD

For molecular dynamics simulations, you must use a LAMMPS data file with the
correct atom ordering. The recommended workflow.



## References

- ABINIT: https://www.abinit.org
- Multibinit: https://docs.abinit.org/topics/multibinit/
- LAMMPS: https://www.lammps.org

## License

This package is distributed under the same license as ABINIT (GPL v3+).
