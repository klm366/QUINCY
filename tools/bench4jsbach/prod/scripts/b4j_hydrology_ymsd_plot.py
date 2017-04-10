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
from b4jplot import mapplot


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

inpath     = kwargs [ 'E_yearsumdiff' ]	
varname    = 'evapotranspiration'
outpath    = kwargs [ 'E_yearsumdiff_plot' ]
unit          = 'mm/year'
unit_TeX   = unit
#title_TeX  = 'Evapotranspiration diff $can - ref$'
title_TeX  = kwargs [ 'plot_title_E_yearsumdiff_plot' ]

mapplot ( inpath, varname, outpath, unit_TeX, title_TeX, midpoint = 0 )

inpath     = kwargs [ 'P_yearsumdiff' ]	
varname    = 'precip_acc'
outpath    = kwargs [ 'P_yearsumdiff_plot' ]
unit       = 'mm/year'
unit_TeX   = unit
#title_TeX  = 'Precipitation diff $can - ref$'
title_TeX  = kwargs [ 'plot_title_P_yearsumdiff_plot' ]

mapplot ( inpath, varname, outpath, unit_TeX, title_TeX, midpoint = 0 )


#####
##### finishing
#####
print "Plotting for field %s with '%s' run successfully." % ( fieldname, __file__ )


print "Exiting %s" % os.path.realpath(__file__)
