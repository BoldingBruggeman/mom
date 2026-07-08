#!/usr/bin/env python

"""This script reads a ``fabm.yaml`` file and writes out entries to stdout
that can be pasted into MOM's ``diag_table``.

It must be pointed to the root directory of a MOM setup.

To use this script, pyfabm (https://fabm.net/python) must be installed.
If you intend to use custom biogeochemical models, pyfabm needs to be built
with those models included. See https://fabm.net/python for details.
"""

import sys
import argparse
import os
from pathlib import Path
import logging

try:
    import pyfabm
except ImportError:
    print("Unable to load pyfabm. See https://fabm.net/python.")
    sys.exit(1)

logging.basicConfig(level=logging.INFO)
pyfabm.logger = logging.getLogger("pyfabm")

parser = argparse.ArgumentParser()
parser.add_argument(
    "expdir",
    help="path to experiment directory",
    nargs="?",
    type=Path,
    default=Path("."),
)
parser.add_argument("--file_name", help="name of output file", default="ocean_fabm")
parser.add_argument(
    "--average",
    help="whether to time-average FABM fields",
    action="store_true",
    default=False,
)
parser.add_argument(
    "--precision",
    type=int,
    help="precision of FABM outputs (1: double precision, 2: float, 4: packed 16 bit integers)",
    default=2,
)
parser.add_argument(
    "--show_hidden",
    help="whether to include hidden FABM fields (entries will commented out)",
    action="store_true",
    default=False,
)
args = parser.parse_args()

model = pyfabm.Model(os.fspath(args.expdir / "fabm.yaml"))

print("#")
print("## FABM")
print("# (variables that are not part of FABM's default output are commented out)")
print("#")
print(
    f'"ocean_model","geolat_t","geolat_t","{args.file_name}",'
    f'"all",.false.,"none",{args.precision}'
)
print(
    f'"ocean_model","geolon_t","geolon_t","{args.file_name}",'
    f'"all",.false.,"none",{args.precision}'
)


def writeVariable(variable: pyfabm.Variable, output=None, name=None):
    if output is None:
        output = variable.output
    if name is None:
        name = variable.output_name
    if not (output or args.show_hidden):
        return
    prefix = "" if output else "#"
    print(
        f'{prefix}"ocean_model","{name}","{name}","{args.file_name}","all",{".true." if args.average else ".false."},"none",{args.precision}'
    )


print("#")
print("# State variables")
print("#")
for variable in model.state_variables:
    writeVariable(variable)
print("#")
print("# Diagnostic variables")
print("#")
for variable in model.diagnostic_variables:
    writeVariable(variable)
print("#")
print("# Conserved quantities")
print("#")
for variable in model.conserved_quantities:
    writeVariable(variable, name=variable.name + "_global_int", output=True)
