#!/usr/bin/env python

# This script reads a fabm.yaml file and writes out entries to stdout that can be pasted into MOM's diag_table.
# It must be run from the root directory of a MOM setup.
# Requirements:
# - the FABM-Python driver (pyfabm) must have been built and installed.

from __future__ import print_function

import sys
import os.path
import argparse

try:
    import pyfabm
except ImportError:
    print('Unable to load pyfabm. See https://github.com/fabm-model/code/wiki/python.')
    sys.exit(1)

parser = argparse.ArgumentParser()
parser.add_argument('expdir', help='path to experiment directory', nargs='?', default='.')
parser.add_argument('--file_name', help='name of output file', default='ocean_fabm')
parser.add_argument('--average', help='whether to time-average FABM fields', action='store_true', default=False)
parser.add_argument('--precision', type=int, help='precision of FABM outputs (1: double precision, 2: float, 4: packed 16 bit integers)', default=2)
parser.add_argument('--show_hidden', help='whether to include hidden FABM fields (entries will commented out)', action='store_true', default=False)
args = parser.parse_args()

model = pyfabm.Model(os.path.join(args.expdir, 'fabm.yaml'))

print('#')
print('## FABM')
print('# (variables that are not part of FABM\'s default output are commented out)')
print('#')
print('"ocean_model","geolat_t","geolat_t","%s","all",.false.,"none",%i' % (args.file_name, args.precision))
print('"ocean_model","geolon_t","geolon_t","%s","all",.false.,"none",%i' % (args.file_name, args.precision))

def writeVariable(variable, output=None, name=None):
    if output is None:
        output = variable.output
    if name is None:
        name = variable.output_name
    if not (output or args.show_hidden):
        return
    prefix = '' if output else '#'
    print('%s"ocean_model","%s","%s","%s","all",%s,"none",%i' % (prefix, name, name, args.file_name, '.true.' if args.average else '.false.', args.precision))

print('#')
print('# State variables')
print('#')
for variable in model.state_variables:
    writeVariable(variable)
print('#')
print('# Diagnostic variables')
print('#')
for variable in model.diagnostic_variables:
    writeVariable(variable)
print('#')
print('# Conserved quantities')
print('#')
for variable in model.conserved_quantities:
    writeVariable(variable, name=variable.name+'_global_int', output=True)
