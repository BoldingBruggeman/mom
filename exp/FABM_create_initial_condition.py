#!/usr/bin/env python

"""Create MOM restart files for FABM state variables.

This script reads ``fabm.yaml`` and writes MOM restart files for all FABM
state variables. It must be pointed to the experiment directory that contains
the ``INPUT`` subdirectory used by MOM.

To use this script, pyfabm (https://fabm.net/python) must be installed.
If you intend to use custom biogeochemical models, pyfabm needs to be built
with those models included. See https://fabm.net/python for details.
"""

import argparse
import sys
import os
from pathlib import Path
from typing import Mapping, Optional, Iterable, Tuple, Union

import netCDF4

try:
    import pyfabm
except ImportError:
    print("Unable to load pyfabm. See https://fabm.net/python.")
    sys.exit(1)


def copy_dimensions(
    ncin: netCDF4.Dataset,
    ncout: netCDF4.Dataset,
    dimensions: Iterable[str],
):
    """Copy selected dimensions and corresponding coordinate variables.

    Args:
        ncin: Source NetCDF dataset.
        ncout: Destination NetCDF dataset.
        dimensions: Names of dimensions to copy from ``ncin`` to ``ncout``.
    """
    for name in dimensions:
        ncout.createDimension(name, ncin.dimensions[name].size)
        if name in ncin.variables:
            ncvar_in = ncin.variables[name]
            ncvar = ncout.createVariable(name, ncvar_in.dtype, ncvar_in.dimensions)
            ncvar[...] = ncvar_in[...]
            for attname in ncvar_in.ncattrs():
                setattr(ncvar, attname, getattr(ncvar_in, attname))


def run(
    expdir: Union[os.PathLike, str],
    rho: float = 1025.0,
    custom: Optional[Mapping[str, Tuple[str, str]]] = {},
) -> None:
    """Generate FABM restart files for a MOM experiment.

    Args:
        expdir: Path to the experiment directory containing ``fabm.yaml`` and
            the ``INPUT`` folder.
        rho: Reference seawater density (kg m-3) used to convert FABM interior
            tracers from per-volume to per-mass units.
        custom: Optional mapping from FABM output variable name to
            ``(path, expression)``. The expression is evaluated against the
            netCDF variables from ``path`` to construct an initial field.
    """
    exp_path = Path(expdir)
    model = pyfabm.Model(os.fspath(exp_path / "fabm.yaml"))

    input_dir = exp_path / "INPUT"
    template_file = input_dir / "ocean_temp_salt.res.nc"
    print(f"Using template file: {template_file}")
    nc_ts = netCDF4.Dataset(template_file)
    dimensions_3d = nc_ts.variables["temp"].dimensions
    dimensions_2d = dimensions_3d[:-3] + dimensions_3d[-2:]
    dtype = nc_ts.variables["temp"].dtype

    state_3d_file = input_dir / "ocean_fabm.res.nc"
    print(f"Writing 3D restart file: {state_3d_file}")
    with netCDF4.Dataset(state_3d_file, "w") as nc3d:
        copy_dimensions(nc_ts, nc3d, dimensions_3d)
        for variable in model.interior_state_variables:
            print(f"  {variable.output_name}: ", end="")
            ncvar = nc3d.createVariable(variable.output_name, dtype, dimensions_3d)
            ncvar.units = f"{variable.units} m3 kg-1"
            ncvar.long_name = variable.long_path
            value = variable.value
            if variable.output_name in custom:
                path, name = custom[variable.output_name]
                with netCDF4.Dataset(path) as nc:
                    value = eval(name, {n: v[...] for (n, v) in nc.variables.items()})
                print(
                    f"[{name} read from {path},"
                    f" mean = {value.mean()}, min = {value.min()}, max = {value.max()}]"
                )
            else:
                print(value)

            # Divide the initial value of interior state variable by density (default: 1025 kg m-3),
            # as MOM tracks tracer per seawater mass, rather than tracer per seawater volume as in FABM.
            ncvar[...] = value / rho

    state_2d_file = input_dir / "ocean_fabm_2d.res.nc"
    print(f"Writing 2D restart file: {state_2d_file}")
    with netCDF4.Dataset(state_2d_file, "w") as nc2d:
        copy_dimensions(nc_ts, nc2d, dimensions_2d)
        for variable in tuple(model.bottom_state_variables) + tuple(
            model.surface_state_variables
        ):
            print(f"  {variable.output_name}: ", end="")
            ncvar = nc2d.createVariable(variable.output_name, dtype, dimensions_2d)
            ncvar.units = variable.units
            ncvar.long_name = variable.long_path
            value = variable.value
            if variable.output_name in custom:
                path, name = custom[variable.output_name]
                with netCDF4.Dataset(path) as nc:
                    value = eval(
                        name, dict([(n, v[...]) for (n, v) in nc.variables.items()])
                    )
                print(
                    f"[{name} read from {path},"
                    f" mean = {value.mean()}, min = {value.min()}, max = {value.max()}]"
                )
            else:
                print(value)
            ncvar[...] = value

    chl_file = input_dir / "ocean_chl.res.nc"
    print(f"Writing chlorophyll restart (all 0): {chl_file}")
    with netCDF4.Dataset(chl_file, "w") as ncchl:
        copy_dimensions(nc_ts, ncchl, dimensions_3d)
        ncvar = ncchl.createVariable("chl", dtype, dimensions_3d)
        ncvar[...] = 0.0

    irr_file = input_dir / "ocean_irr.res.nc"
    print(f"Writing irradiance restart (all 0): {irr_file}")
    with netCDF4.Dataset(irr_file, "w") as ncirr:
        copy_dimensions(nc_ts, ncirr, dimensions_3d)
        ncvar = ncirr.createVariable("irr", dtype, dimensions_3d)
        ncvar[...] = 0.0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument(
        "expdir",
        help="path to experiment directory",
        nargs="?",
        type=Path,
        default=Path("."),
    )
    parser.add_argument(
        "--rho",
        help="reference density of sea water (kg m-3)",
        type=float,
        default=1025.0,
    )
    args = parser.parse_args()
    run(args.expdir, args.rho)
