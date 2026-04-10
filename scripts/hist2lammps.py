#!/usr/bin/env python3
"""Convert ABINIT/Multibinit HIST.nc file to LAMMPS data file.

Ensures atom ordering matches what Multibinit expects, which is critical
for the dipole-dipole (Ewald) calculation.

Usage: python hist2lammps.py input_HIST.nc [output.lmp] [timestep]
"""

import sys
import subprocess


def read_hist_nc(filename):
    result = subprocess.run(['ncdump', '-h', filename], capture_output=True, text=True)
    header = result.stdout
    
    dims = {}
    for line in header.split('\n'):
        if '=' in line and ';' in line:
            parts = line.strip().split('=')
            if len(parts) == 2:
                name = parts[0].strip()
                val = parts[1].strip().rstrip(';')
                if val != 'UNLIMITED':
                    try:
                        dims[name] = int(val)
                    except ValueError:
                        pass
    
    vars_to_get = ['xred', 'xcart', 'typat', 'znucl', 'amu', 'rprimd', 'acell']
    data = {}
    
    for var in vars_to_get:
        result = subprocess.run(['ncdump', '-v', var, filename], 
                              capture_output=True, text=True)
        content = result.stdout
        
        in_data = False
        values = []
        for line in content.split('\n'):
            if 'data:' in line:
                in_data = True
                continue
            if in_data:
                if line.strip() == '}':
                    break
                if ' = ' in line and var in line.split(' = ')[0]:
                    line = line.split(' = ', 1)[1]
                nums = line.replace(',', ' ').replace(';', ' ').split()
                for n in nums:
                    try:
                        values.append(float(n))
                    except ValueError:
                        pass
        
        data[var] = values
    
    return dims, data


def write_lammps_data(filename, dims, data, timestep=0):
    natom = dims.get('natom', 40)
    ntypat = dims.get('ntypat', 3)
    
    typat = data.get('typat', [])
    amu = data.get('amu', [])
    
    xred_flat = data.get('xred', [])
    xred = []
    for i in range(natom):
        xred.append([
            xred_flat[timestep * natom * 3 + i * 3],
            xred_flat[timestep * natom * 3 + i * 3 + 1],
            xred_flat[timestep * natom * 3 + i * 3 + 2]
        ])
    
    rprimd_flat = data.get('rprimd', [])
    rprimd = []
    for i in range(3):
        rprimd.append([
            rprimd_flat[timestep * 9 + i * 3],
            rprimd_flat[timestep * 9 + i * 3 + 1],
            rprimd_flat[timestep * 9 + i * 3 + 2]
        ])
    
    bohr_to_ang = 0.529177210903
    
    lx = rprimd[0][0] * bohr_to_ang
    ly = rprimd[1][1] * bohr_to_ang
    lz = rprimd[2][2] * bohr_to_ang
    
    xcart = []
    for i in range(natom):
        x = xred[i][0] * lx
        y = xred[i][1] * ly
        z = xred[i][2] * lz
        xcart.append([x, y, z])
    
    with open(filename, 'w') as f:
        f.write(f"LAMMPS data file from HIST.nc (timestep {timestep})\n\n")
        f.write(f"{natom} atoms\n")
        f.write(f"{ntypat} atom types\n\n")
        f.write(f"0.0 {lx:.10f} xlo xhi\n")
        f.write(f"0.0 {ly:.10f} ylo yhi\n")
        f.write(f"0.0 {lz:.10f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t in range(1, ntypat + 1):
            f.write(f"{t} {amu[t-1]:.10f}\n")
        f.write("\nAtoms\n\n")
        for i in range(natom):
            atom_type = int(typat[i])
            x, y, z = xcart[i]
            f.write(f"{i+1} {atom_type} {x:.10f} {y:.10f} {z:.10f}\n")
    
    print(f"Wrote {filename}")
    print(f"  Atoms: {natom}, Types: {ntypat}")
    print(f"  Box: {lx:.4f} x {ly:.4f} x {lz:.4f} Angstrom")
    print("\nAtom ordering (first 10):")
    for i in range(min(10, natom)):
        atom_type = int(typat[i])
        element = ['Ba', 'Hf', 'O'][atom_type - 1] if atom_type <= 3 else f'T{atom_type}'
        print(f"  {i+1}: {element} at ({xcart[i][0]:.4f}, {xcart[i][1]:.4f}, {xcart[i][2]:.4f})")


def main():
    if len(sys.argv) < 2:
        print("Usage: python hist2lammps.py input_HIST.nc [output.lmp] [timestep]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'data.lmp'
    timestep = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    
    dims, data = read_hist_nc(input_file)
    write_lammps_data(output_file, dims, data, timestep)


if __name__ == '__main__':
    main()
