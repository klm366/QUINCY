#! /usr/bin/env python
# -*- coding:  utf-8 -*-

import sys, os
import numpy as np
from netCDF4 import Dataset


def get_value ( filepath, varname ):
    nc_fid = Dataset ( filepath )
    lat    = nc_fid.variables [ 'lat' ][ : ]
    vardat = nc_fid.variables [ varname ] [ : ]
    nc_fid.close ()
    print vardat  [ :, :, 0 ]
    return vardat [ :, :, 0 ]
    

print "Starting %s" % os.path.realpath(__file__)


#####
##### bench4jsbach library to environment path
#####
try:
    b4j_base = os.environ [ 'BENCH4JSBACH_BASE' ]
except KeyError: 
    print "The environment variable 'BENCH4JSBACH_BASE' is not set."
    sys.exit(1)

bin_dir = os.path.join ( b4j_base, 'bin' )
if os.path.isdir ( bin_dir ):
    sys.path.append ( bin_dir )
else:
    raise Exception ( "Directory %s with source files not found." % bin_dir )
from util import print_item_or_list_rec
from b4jplot import hist_plot_wrapper


#####
##### check input parameter
#####
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.strip ().split ( ':', 1 ) for token in sys.argv [ 1: ])
print kwargs.keys ()
print kwargs
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )


#####
##### argument definition and plot method calls
#####

# varnames = { 'E'   : 'evapotranspiration',
#              'P'   : 'precip_acc',
#              'RD'  : 'disch',
#              'WSC' : 'waterstorage_change' }
varnames = [  'evapotranspiration',
              'precip_acc',
              'disch',
              'waterstorage_change' ]
unit = "mm/year"
can_name = kwargs [ "can_name" ]
canno = kwargs [ "canno" ]
suptitle = "Hydrology mean year sums"
input_can = kwargs [ "H_ysfm" ]
subtitle_can = can_name
print "can_name"
print can_name
print subtitle_can
#subtitle_delta = "Differences of meanyearsums to %s" % kwargs [ "can_name" ]
outpath = kwargs [ "H_ysdfm_plot" ]

filenametemplate = "H_ysdfm_%(refacronym)s" 
reftoken         = "ref"
nchar            = 4
reflist          = [ k for k in kwargs.keys () if k.startswith ( reftoken ) and len ( k ) == nchar ]
input_delta      = []
xlabel_delta     = []
subtitle_delta   = []
for ref in reflist:
    input_delta.append ( kwargs [ filenametemplate % { 'refacronym' : ref }])
    subtitle_delta.append ( kwargs [ ref ])
    xlabel_delta.append ( "$\Delta$ " + ref )

hist_plot_wrapper ( varnames, unit, can_name, canno, suptitle,
                    input_can, subtitle_can,
                    input_delta, subtitle_delta, xlabel_delta,
                    outpath )

#####
##### exiting
#####
print "Exiting %s" % os.path.realpath(__file__)
