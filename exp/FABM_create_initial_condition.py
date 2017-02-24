#!/usr/bin/env python

# This script reads a fabm.yaml file, and writes out MOM restart files for all FABM state variables.
# It must be run from the INPUT directory of a MOM setup.
# Requirements:
# - the FABM-Python driver (pyfabm) must have been built and installed.
# - the netCDF Operator toolkit (http://nco.sourceforge.net)

from __future__ import print_function

import argparse
import sys
import os.path
import subprocess

try:
    import pyfabm
except ImportError:
    print('Unable to load pyfabm. See https://github.com/fabm-model/fabm/wiki/python.')
    sys.exit(1)

def run(expdir, rho=1025.):
    model = pyfabm.Model(os.path.join(expdir, 'fabm.yaml'))

    os.chdir(os.path.join(expdir, 'INPUT'))

    print('Writing NCO script file for 3D fields:')
    with open('fabm_init.nco', 'w') as f:
        for variable in model.interior_state_variables:
            # NB the initial value of interior state variable by a reference density 1025 kg m-3,
            # as MOM tracks tracer per seawater mass, rather than tracer per seawater volume as in FABM.
            print('  %s: %s' % (variable.output_name, variable.value))
            f.write('%s=temp*0+%s/%s;%s@long_name="%s";%s@units="%s m3 kg-1";\n' % (variable.output_name, variable.value, rho, variable.output_name, variable.long_path, variable.output_name, variable.units))

    print('Creating 3D restart file...')
    subprocess.check_call(['ncap2', 'ocean_temp_salt.res.nc', 'ocean_fabm.res.nc', '-S', 'fabm_init.nco', '-O', '-v'])
    subprocess.check_call(['ncatted', '-a', 'checksum,,d,,', 'ocean_fabm.res.nc'])

    print('Writing NCO script file for 2D fields:')
    with open('fabm_init_2d.nco', 'w') as f:
        for variable in model.bottom_state_variables:
            print('  %s: %s' % (variable.output_name, variable.value))
            f.write('%s=temp(:,0,:,:)*0+%s;%s@long_name="%s";%s@units="%s";\n' % (variable.output_name, variable.value, variable.output_name, variable.long_path, variable.output_name, variable.units))
        for variable in model.surface_state_variables:
            print('  %s: %s' % (variable.output_name, variable.value))
            f.write('%s=temp(:,0,:,:)*0+%s;%s@long_name="%s";%s@units="%s";\n' % (variable.output_name, variable.value, variable.output_name, variable.long_path, variable.output_name, variable.units))

    print('Creating 2D restart file...')
    subprocess.check_call(['ncap2', 'ocean_temp_salt.res.nc', 'ocean_fabm_2d.res.nc', '-S', 'fabm_init_2d.nco', '-O', '-v'])
    subprocess.check_call(['ncatted', '-a', 'checksum,,d,,', 'ocean_fabm_2d.res.nc'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('expdir', help='path to experiment directory', nargs='?', default='.')
    parser.add_argument('--rho', help='reference density of sea water (kg m-3)', type=float, default=1025.)
    args = parser.parse_args()
    run(args.expdir, args.rho)
