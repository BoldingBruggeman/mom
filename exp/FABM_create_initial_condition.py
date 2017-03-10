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
import netCDF4

try:
    import pyfabm
except ImportError:
    print('Unable to load pyfabm. See https://github.com/fabm-model/fabm/wiki/python.')
    sys.exit(1)

def copy_dimensions(ncin, ncout, dimensions):
    for name in dimensions:
        ncout.createDimension(name, ncin.dimensions[name].size)
        if name in ncin.variables:
            ncvar_in = ncin.variables[name]
            ncvar = ncout.createVariable(name, ncvar_in.dtype, ncvar_in.dimensions)
            ncvar[...] = ncvar_in[...]
            for attname in ncvar_in.ncattrs():
                setattr(ncvar, attname, getattr(ncvar_in, attname))

def run(expdir, rho=1025., custom={}):
    model = pyfabm.Model(os.path.join(expdir, 'fabm.yaml'))

    nc_ts = netCDF4.Dataset(os.path.join(expdir, 'INPUT/ocean_temp_salt.res.nc'), 'r')
    dimensions_3d = nc_ts.variables['temp'].dimensions
    dimensions_2d = dimensions_3d[:-3] + dimensions_3d[-2:]
    dtype = nc_ts.variables['temp'].dtype

    print('Writing 3D restart file: ocean_fabm.res.nc')
    with netCDF4.Dataset(os.path.join(expdir, 'INPUT/ocean_fabm.res.nc'), 'w') as nc3d:
        copy_dimensions(nc_ts, nc3d, dimensions_3d)
        for variable in model.interior_state_variables:
            print('  %s: ' % variable.output_name, end='')
            ncvar = nc3d.createVariable(variable.output_name, dtype, dimensions_3d)
            ncvar.units = '%s m3 kg-1' % variable.units
            ncvar.long_name = variable.long_path
            value = variable.value
            if variable.output_name in custom:
                path, name = custom[variable.output_name]
                with netCDF4.Dataset(path) as nc:
                    value = eval(name, dict([(n, v[...]) for (n, v) in nc.variables.items()]))
                print('[%s read from %s, mean = %s, min = %s, max = %s]' % (name, path, value.mean(), value.min(), value.max()))
            else:
                print(value)

            # Divide the initial value of interior state variable by density (default: 1025 kg m-3),
            # as MOM tracks tracer per seawater mass, rather than tracer per seawater volume as in FABM.
            ncvar[...] = value/rho

    print('Writing 2D restart file: ocean_fabm_2d.res.nc')
    with netCDF4.Dataset(os.path.join(expdir, 'INPUT/ocean_fabm_2d.res.nc'), 'w') as nc2d:
        copy_dimensions(nc_ts, nc2d, dimensions_2d)
        for variable in tuple(model.bottom_state_variables)+tuple(model.surface_state_variables):
            print('  %s: ' % variable.output_name, end='')
            ncvar = nc2d.createVariable(variable.output_name, dtype, dimensions_2d)
            ncvar.units = variable.units
            ncvar.long_name = variable.long_path
            value = variable.value
            if variable.output_name in custom:
                path, name = custom[variable.output_name]
                with netCDF4.Dataset(path) as nc:
                    value = eval(name, dict([(n, v[...]) for (n, v) in nc.variables.items()]))
                print('[%s read from %s, mean = %s, min = %s, max = %s]' % (name, path, value.mean(), value.min(), value.max()))
            else:
                print(value)
            ncvar[...] = value

    print('Writing chlorophyll restart (all 0): ocean_chl.res.nc')
    with netCDF4.Dataset(os.path.join(expdir, 'INPUT/ocean_chl.res.nc'), 'w') as ncchl:
        copy_dimensions(nc_ts, ncchl, dimensions_3d)
        ncvar = ncchl.createVariable('chl', dtype, dimensions_3d)
        ncvar[...] = 0.

    print('Writing irradiance restart (all 0): ocean_irr.res.nc')
    with netCDF4.Dataset(os.path.join(expdir, 'INPUT/ocean_irr.res.nc'), 'w') as ncirr:
        copy_dimensions(nc_ts, ncirr, dimensions_3d)
        ncvar = ncirr.createVariable('irr', dtype, dimensions_3d)
        ncvar[...] = 0.

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('expdir', help='path to experiment directory', nargs='?', default='.')
    parser.add_argument('--rho', help='reference density of sea water (kg m-3)', type=float, default=1025.)
    args = parser.parse_args()
    run(args.expdir, args.rho)
