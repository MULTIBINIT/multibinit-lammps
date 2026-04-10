/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS Plugin for Multibinit effective potential from ABINIT.

   This file registers the pair_multibinit style as a LAMMPS plugin.
------------------------------------------------------------------------- */

#include "lammpsplugin.h"
#include "version.h"

#include "pair_multibinit.h"

#include <cstring>

using namespace LAMMPS_NS;

static Pair *multibinit_creator(LAMMPS *lmp)
{
  return new PairMultibinit(lmp);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  // Register multibinit pair style
  plugin.version = LAMMPS_VERSION;
  plugin.style = "pair";
  plugin.name = "multibinit";
  plugin.info = "Multibinit effective potential pair style v1.0";
  plugin.author = "ABINIT/LAMMPS integration";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &multibinit_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);
}
