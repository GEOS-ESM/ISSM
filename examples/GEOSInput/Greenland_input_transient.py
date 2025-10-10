import devpath
import os,sys
ISSM_DIR = os.getenv('ISSM_DIR') # for binaries

import numpy as np
from model import *
from loadmodel import loadmodel
from clusters.discover_geos import export_discover
from loadresultsfromdisk import loadresultsfromdisk
from marshall import marshall

md = loadmodel('./Models/Greenland_inversion.nc')
md = loadresultsfromdisk(md,'GreenlandGEOS.outbin')
md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient


# Write the binary input file
# Additional options
md.inversion.iscontrol = 0
md.transient.requested_outputs = ['IceVolume', 'TotalSmb', 'SmbMassBalance']
md.settings.waitonlock = 0
md.private.solution = 'Transient'
md.toolkits = toolkits()
marshall(md) # create .bin file
md.toolkits = toolkits()
export_discover(md,'./Models/Greenland_initialization.nc',delete_rundir=True)

