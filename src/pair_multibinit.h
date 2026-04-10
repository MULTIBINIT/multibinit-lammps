/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   Pair style for Multibinit effective potential from ABINIT.
   Supports MPI communicator passing for proper parallel integration.

   References:
   - ABINIT: https://www.abinit.org
   - Multibinit: https://docs.abinit.org/topics/multibinit/
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(multibinit,PairMultibinit);
// clang-format on
#else

#ifndef LMP_PAIR_MULTIBINIT_H
#define LMP_PAIR_MULTIBINIT_H

#include "pair.h"

#include <mpi.h>

namespace LAMMPS_NS {

class PairMultibinit : public Pair {
 public:
  PairMultibinit(class LAMMPS *);
  ~PairMultibinit() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  double init_one(int, int) override;

 protected:
  // Multibinit handle (opaque pointer)
  void *mb_handle;
  bool mb_initialized;

  // Potential file paths
  char *sys_filename;    // DDB or XML system file
  char *coeff_filename;  // Optional coefficients file

  // Supercell dimensions
  int ncell[3];

  // q-point grid for dipole-dipole
  int ngqpt[3];
  int dipdip;

  // MPI communicator (Fortran handle)
  MPI_Fint mpi_comm_f;

  // Unit conversion factors (Multibinit uses atomic units)
  // LAMMPS metal units: energy in eV, distance in Angstrom
  // 1 Ha = 27.211386245988 eV
  // 1 Bohr = 0.529177210903 Angstrom
  static constexpr double HA_TO_EV = 27.211386245988;
  static constexpr double BOHR_TO_ANG = 0.529177210903;
  static constexpr double EV_TO_HA = 1.0 / HA_TO_EV;
  static constexpr double ANG_TO_BOHR = 1.0 / BOHR_TO_ANG;
  // Force: Ha/Bohr to eV/Angstrom
  static constexpr double FORCE_CONV = HA_TO_EV / BOHR_TO_ANG;
  // Stress: Ha/Bohr^3 to eV/Angstrom^3 (then to pressure units)
  static constexpr double STRESS_CONV = HA_TO_EV / (BOHR_TO_ANG * BOHR_TO_ANG * BOHR_TO_ANG);

  // Internal buffers for coordinate conversion
  double *positions_bohr;   // [3*natom] in Bohr
  double *lattice_bohr;     // [9] in Bohr
  double *forces_ha;        // [3*natom] in Ha/Bohr
  double stress_ha[6];      // [6] in Ha/Bohr^3

  int nmax;  // Current allocation size

  // MPI parallelization buffers (for domain decomposition)
  double *all_positions_bohr;  // [3*natoms_total] gathered from all procs
  double *all_forces_ha;       // [3*natoms_total] to scatter back
  int *recv_counts;            // [nprocs] for MPI_Allgatherv
  int *displs;                 // [nprocs] for MPI_Allgatherv
  int nmax_all;                // Allocation size for all-atom arrays

  // Atom mapping for LAMMPS -> Multibinit
  int *lammps_to_mb;             // LAMMPS atom i -> Multibinit atom index
  int *mb_to_lammps;             // Multibinit atom i -> LAMMPS atom index
  double *mb_positions_bohr;    // [3*natom] ideal positions in Bohr
  double *mb_lattice_bohr;     // [9] ideal lattice in Bohr
  bool atom_map_built;            // Flag: mapping has been built

  virtual void allocate();
  void init_multibinit();
  void free_multibinit();
  void create_atoms_from_multibinit();
  void build_atom_map();        // Build the mapping
};

}    // namespace LAMMPS_NS

#endif
#endif
