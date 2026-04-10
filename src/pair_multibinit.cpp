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

#include "pair_multibinit.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "lammps.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cstring>
#include <cmath>
#include <vector>

using namespace LAMMPS_NS;

extern "C" {
void mb_init_potential_simple(const char *sys_filename, int sys_len,
 const char *coeff_filename, int coeff_len,
 const int *ncell, const int *ngqpt, int dipdip,
 void **handle, int *status);

void mb_init_potential_simple_with_comm(const char *sys_filename, int sys_len,
 const char *coeff_filename, int coeff_len,
 const int *ncell, const int *ngqpt, int dipdip,
 int mpi_comm_f, void **handle, int *status);

void mb_evaluate(void *handle, const double *positions, const double *lattice,
                 int natom, double *energy, double *forces, double *stresses,
                 int *status);

void mb_evaluate_with_comm(void *handle, const double *positions, const double *lattice,
                            int natom, int mpi_comm_f, double *energy, double *forces,
                            double *stresses, int *status);

void mb_free_potential(void *handle, int *status);

void mb_get_supercell_structure(void *handle, int *natom, int *species,
                                 double *positions, double *lattice, int *status);
}

/* ---------------------------------------------------------------------- */

PairMultibinit::PairMultibinit(LAMMPS *lmp) :
    Pair(lmp),
    mb_handle(nullptr),
    mb_initialized(false),
    sys_filename(nullptr),
    coeff_filename(nullptr),
    mpi_comm_f(0),
    positions_bohr(nullptr),
    lattice_bohr(nullptr),
    forces_ha(nullptr),
    nmax(0),
    all_positions_bohr(nullptr),
    all_forces_ha(nullptr),
    recv_counts(nullptr),
    displs(nullptr),
    nmax_all(1),
    lammps_to_mb(nullptr),
    mb_to_lammps(nullptr),
    mb_positions_bohr(nullptr),
    mb_lattice_bohr(nullptr),
    atom_map_built(false)
{
  ncell[0] = ncell[1] = ncell[2] = 1;
  ngqpt[0] = ngqpt[1] = ngqpt[2] = 2;
  dipdip = 1;

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  no_virial_fdotr_compute = 1;

  for (int i = 0; i < 6; i++) stress_ha[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

PairMultibinit::~PairMultibinit()
{
  if (copymode) return;

  free_multibinit();

  memory->destroy(positions_bohr);
  memory->destroy(lattice_bohr);
  memory->destroy(forces_ha);
  memory->destroy(setflag);
  memory->destroy(cutsq);
  memory->destroy(all_positions_bohr);
  memory->destroy(all_forces_ha);
  memory->destroy(recv_counts);
  memory->destroy(displs);

  delete[] lammps_to_mb;
  delete[] mb_to_lammps;
  delete[] mb_positions_bohr;
  delete[] mb_lattice_bohr;
  delete[] sys_filename;
  delete[] coeff_filename;
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");

  for (int i = 1; i <= n; i++)
    for (int j = 1; j <= n; j++) setflag[i][j] = 1;
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::settings(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR, "Pair style multibinit requires at least 3 arguments");

  ncell[0] = utils::inumeric(FLERR, arg[0], false, lmp);
  ncell[1] = utils::inumeric(FLERR, arg[1], false, lmp);
  ncell[2] = utils::inumeric(FLERR, arg[2], false, lmp);

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ngqpt") == 0) {
      if (iarg + 3 >= narg)
        error->all(FLERR, "ngqpt requires 3 values");
      ngqpt[0] = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      ngqpt[1] = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      ngqpt[2] = utils::inumeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "dipdip") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "dipdip requires 1 value");
      dipdip = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Unknown pair style setting: {}", arg[iarg]);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::coeff(int narg, char **arg)
{
  if (narg < 3) error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();

  delete[] sys_filename;
  sys_filename = utils::strdup(arg[2]);

  if (narg > 3) {
    delete[] coeff_filename;
    coeff_filename = utils::strdup(arg[3]);
  }

  init_multibinit();
}

/* ---------------------------------------------------------------------- */

double PairMultibinit::init_one(int i, int j)
{
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::init_multibinit()
{
  if (mb_initialized) return;

  int status = 0;
  int sys_len = strlen(sys_filename);
  int coeff_len = coeff_filename ? strlen(coeff_filename) : 0;

  const char *coeff_ptr = coeff_filename ? coeff_filename : "";
  if (!coeff_filename) coeff_len = 0;

  MPI_Comm lammps_comm = lmp->world;
  mpi_comm_f = MPI_Comm_c2f(lammps_comm);

  if (comm->me == 0) {
    utils::logmesg(lmp, "Multibinit: Initializing with LAMMPS MPI communicator\n");
  }

  mb_init_potential_simple_with_comm(sys_filename, sys_len, coeff_ptr, coeff_len,
                                      ncell, ngqpt, dipdip, mpi_comm_f,
                                      &mb_handle, &status);

  if (status != 0 || mb_handle == nullptr) {
    error->all(FLERR, "Failed to initialize Multibinit potential: status {}", status);
  }

  mb_initialized = true;

  int mb_natom = 0;
  mb_get_supercell_structure(mb_handle, &mb_natom, nullptr, nullptr, nullptr, &status);

  if (comm->me == 0) {
    utils::logmesg(lmp, "Multibinit initialized: {} atoms in supercell ({}x{}x{})\n",
                    mb_natom, ncell[0], ncell[1], ncell[2]);
  }
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::free_multibinit()
{
  if (mb_handle) {
    int status = 0;
    mb_free_potential(mb_handle, &status);
    mb_handle = nullptr;
  }
  mb_initialized = false;
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::create_atoms_from_multibinit()
{
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::build_atom_map()
{
  if (atom_map_built) return;

  int me = comm->me;
  int nprocs = comm->nprocs;

  bigint natoms_big = atom->natoms;
  int natoms = static_cast<int>(natoms_big);

  int mb_natom = 0;
  int status = 0;
  mb_get_supercell_structure(mb_handle, &mb_natom, nullptr, nullptr, nullptr, &status);

  if (status != 0) {
    error->all(FLERR, "Failed to get atom count from Multibinit: status {}", status);
  }

  if (mb_natom != natoms) {
    error->all(FLERR, "Atom count mismatch: LAMMPS has {}, Multibinit expects {}", natoms, mb_natom);
  }

  lammps_to_mb = new int[natoms];
  mb_to_lammps = new int[natoms];
  mb_positions_bohr = new double[3 * natoms];
  mb_lattice_bohr = new double[9];

  if (me == 0) {
    int *mb_species = new int[natoms];
    mb_get_supercell_structure(mb_handle, &mb_natom, mb_species, mb_positions_bohr, mb_lattice_bohr, &status);
    if (status != 0) {
      delete[] mb_species;
      error->one(FLERR, "Failed to get structure from Multibinit: status {}", status);
    }
    delete[] mb_species;
  }

  int nlocal = atom->nlocal;
  double **x = atom->x;
  tagint *tag = atom->tag;

  std::vector<double> local_data(4 * nlocal);
  for (int i = 0; i < nlocal; i++) {
    local_data[4*i + 0] = static_cast<double>(tag[i]);
    local_data[4*i + 1] = x[i][0];
    local_data[4*i + 2] = x[i][1];
    local_data[4*i + 3] = x[i][2];
  }

  int *recvcounts = new int[nprocs];
  int *displs_arr = new int[nprocs];
  int mycount = 4 * nlocal;
  MPI_Gather(&mycount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, world);

  int total_count = 0;
  if (me == 0) {
    displs_arr[0] = 0;
    for (int p = 0; p < nprocs; p++) {
      if (p > 0) displs_arr[p] = displs_arr[p-1] + recvcounts[p-1];
      total_count += recvcounts[p];
    }
  }

  std::vector<double> all_data;
  if (me == 0) {
    all_data.resize(total_count);
  }

  MPI_Gatherv(local_data.data(), mycount, MPI_DOUBLE,
              all_data.data(), recvcounts, displs_arr, MPI_DOUBLE,
              0, world);

  delete[] recvcounts;
  delete[] displs_arr;

  if (me == 0) {
    double boxx = mb_lattice_bohr[0] * BOHR_TO_ANG;
    double boxy = mb_lattice_bohr[4] * BOHR_TO_ANG;
    double boxz = mb_lattice_bohr[8] * BOHR_TO_ANG;

    std::vector<bool> mb_used(mb_natom, false);

    int nlocal_total = total_count / 4;
    for (int i = 0; i < nlocal_total; i++) {
      int lammps_tag = static_cast<int>(all_data[4*i + 0]);
      double x_ang = all_data[4*i + 1];
      double y_ang = all_data[4*i + 2];
      double z_ang = all_data[4*i + 3];

      double min_dist_sq = 1e18;
      int best_j = -1;

      for (int j = 0; j < mb_natom; j++) {
        if (mb_used[j]) continue;

        double mb_x_ang = mb_positions_bohr[3 * j] * BOHR_TO_ANG;
        double mb_y_ang = mb_positions_bohr[3 * j + 1] * BOHR_TO_ANG;
        double mb_z_ang = mb_positions_bohr[3 * j + 2] * BOHR_TO_ANG;

        double dx = x_ang - mb_x_ang;
        double dy = y_ang - mb_y_ang;
        double dz = z_ang - mb_z_ang;

        if (boxx > 0) dx -= round(dx / boxx) * boxx;
        if (boxy > 0) dy -= round(dy / boxy) * boxy;
        if (boxz > 0) dz -= round(dz / boxz) * boxz;

        double dist_sq = dx * dx + dy * dy + dz * dz;

        if (dist_sq < min_dist_sq) {
          min_dist_sq = dist_sq;
          best_j = j;
        }
      }

      if (best_j < 0) {
        error->one(FLERR, "Could not find matching Multibinit atom for LAMMPS atom {}", lammps_tag);
      }

      lammps_to_mb[lammps_tag - 1] = best_j;
      mb_to_lammps[best_j] = lammps_tag - 1;
      mb_used[best_j] = true;
    }

    utils::logmesg(lmp, "Atom mapping built: {} atoms matched\n", mb_natom);
  }

  MPI_Bcast(lammps_to_mb, natoms, MPI_INT, 0, world);
  MPI_Bcast(mb_to_lammps, natoms, MPI_INT, 0, world);

  atom_map_built = true;
}

/* ---------------------------------------------------------------------- */

void PairMultibinit::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (!mb_initialized) {
    error->all(FLERR, "Multibinit potential not initialized");
  }

  if (!atom_map_built) {
    build_atom_map();
  }

  int me = comm->me;
  int nlocal = atom->nlocal;
  bigint natoms_big = atom->natoms;
  int natoms = static_cast<int>(natoms_big);
  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;

  double boxx = domain->boxhi[0] - domain->boxlo[0];
  double boxy = domain->boxhi[1] - domain->boxlo[1];
  double boxz = domain->boxhi[2] - domain->boxlo[2];

  if (natoms > nmax_all) {
    memory->destroy(all_positions_bohr);
    memory->destroy(all_forces_ha);
    memory->destroy(lattice_bohr);
    memory->create(all_positions_bohr, 3 * natoms, "pair:all_positions_bohr");
    memory->create(all_forces_ha, 3 * natoms, "pair:all_forces_ha");
    memory->create(lattice_bohr, 9, "pair:lattice_bohr");
    nmax_all = natoms;
  }

  int nprocs = comm->nprocs;
  
  int *recvcounts_pos = new int[nprocs];
  int *displs_pos = new int[nprocs];
  int mycount_pos = 3 * nlocal;
  
  MPI_Allgather(&mycount_pos, 1, MPI_INT, recvcounts_pos, 1, MPI_INT, world);

  displs_pos[0] = 0;
  for (int p = 1; p < nprocs; p++) {
    displs_pos[p] = displs_pos[p-1] + recvcounts_pos[p-1];
  }

  double *gathered_positions = new double[3 * natoms];
  std::vector<double> local_pos(3 * nlocal);
  for (int i = 0; i < nlocal; i++) {
    local_pos[3*i + 0] = x[i][0] * ANG_TO_BOHR;
    local_pos[3*i + 1] = x[i][1] * ANG_TO_BOHR;
    local_pos[3*i + 2] = x[i][2] * ANG_TO_BOHR;
  }
  
  MPI_Allgatherv(local_pos.data(), mycount_pos, MPI_DOUBLE,
                 gathered_positions, recvcounts_pos, displs_pos, MPI_DOUBLE, world);

  int *recvcounts_tag = new int[nprocs];
  int *displs_tag = new int[nprocs];
  MPI_Allgather(&nlocal, 1, MPI_INT, recvcounts_tag, 1, MPI_INT, world);

  displs_tag[0] = 0;
  for (int p = 1; p < nprocs; p++) {
    displs_tag[p] = displs_tag[p-1] + recvcounts_tag[p-1];
  }

  int *gathered_tags = new int[natoms];
  std::vector<int> local_tags(nlocal);
  for (int i = 0; i < nlocal; i++) {
    local_tags[i] = static_cast<int>(tag[i]);
  }
  
  MPI_Allgatherv(local_tags.data(), nlocal, MPI_INT,
                 gathered_tags, recvcounts_tag, displs_tag, MPI_INT, world);

  delete[] recvcounts_pos;
  delete[] displs_pos;
  delete[] recvcounts_tag;
  delete[] displs_tag;

  for (int i = 0; i < natoms; i++) {
    int lammps_tag = gathered_tags[i];
    int mb_idx = lammps_to_mb[lammps_tag - 1];
    all_positions_bohr[3*mb_idx + 0] = gathered_positions[3*i + 0];
    all_positions_bohr[3*mb_idx + 1] = gathered_positions[3*i + 1];
    all_positions_bohr[3*mb_idx + 2] = gathered_positions[3*i + 2];
  }

  delete[] gathered_positions;
  delete[] gathered_tags;

  lattice_bohr[0] = boxx * ANG_TO_BOHR;
  lattice_bohr[1] = 0.0;
  lattice_bohr[2] = 0.0;
  lattice_bohr[3] = 0.0;
  lattice_bohr[4] = boxy * ANG_TO_BOHR;
  lattice_bohr[5] = 0.0;
  lattice_bohr[6] = 0.0;
  lattice_bohr[7] = 0.0;
  lattice_bohr[8] = boxz * ANG_TO_BOHR;

  double energy_ha = 0.0;
  int status = 0;

  mb_evaluate_with_comm(mb_handle, all_positions_bohr, lattice_bohr,
                        natoms, mpi_comm_f, &energy_ha, all_forces_ha,
                        stress_ha, &status);

  if (status != 0) {
    error->all(FLERR, "Multibinit evaluation failed: status {}", status);
  }

  for (int i = 0; i < nlocal; i++) {
    int lammps_tag = static_cast<int>(tag[i]);
    int mb_idx = lammps_to_mb[lammps_tag - 1];

    f[i][0] += all_forces_ha[3*mb_idx + 0] * FORCE_CONV;
    f[i][1] += all_forces_ha[3*mb_idx + 1] * FORCE_CONV;
    f[i][2] += all_forces_ha[3*mb_idx + 2] * FORCE_CONV;
  }

  if (eflag_global && me == 0) {
    eng_vdwl += energy_ha * HA_TO_EV;
  }

  if (vflag_global && me == 0) {
    // Convert stress (Hartree/Bohr^3) to virial (eV)
    // virial = -stress × volume (negative sign for thermodynamic convention)
    double volume = boxx * boxy * boxz;
    for (int i = 0; i < 6; i++) {
      virial[i] -= stress_ha[i] * STRESS_CONV * volume;
    }
  }
}
