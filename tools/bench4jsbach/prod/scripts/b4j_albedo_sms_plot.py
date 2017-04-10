#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
from cdo import Cdo

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
from b4jplot import mapplot, mapplot_seasonal

#####
##### check input parameter
#####

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':', 1 ) for token in sys.argv [ 1: ])
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )

fieldname      = kwargs [ "fieldname" ]

#####
##### plotting
#####

inpath     = kwargs [ 'A_seasmean' ]	
varname    = 'albedo'
outpath    = kwargs [ 'A_seasmean_plot' ]
unit       = '1'
unit_TeX   = unit
title_TeX  = kwargs [ 'plot_title_A_seasmean_plot' ]

#mapplot ( inpath, varname, outpath, unit_TeX, title_TeX, colormap = 'RdBu'  )
#outpath_seas = outpath.replace (  'A_seasmean_plot',  'A_seasmean_plot_test-seasonal' )
mapplot_seasonal ( inpath, varname, outpath, unit_TeX, title_TeX,
                   midpoint = 0.173, colormap = 'RdBu'  )
#####
##### finishing
#####

print "Preprocessing for field %s with '%s' run successfully." % ( fieldname, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
