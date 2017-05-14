#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import os, sys
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

#####
##### check input parameter
#####

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':' ) for token in sys.argv [ 1: ])
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )

#####
##### constants definition
#####

fieldname      = "hydrology"
unit           = "mm/day"
unit_factor_EP = 86400
sign_adjust_E     = -1
unit_factor_RD = 86400000

#####
##### set unit to mm/day
#####

cdo = Cdo ()

cdo.setname ( "evapotranspiration",
              input = " -setunit,%s -mulc,%d -mulc,%d %s " % (
              unit, unit_factor_EP, sign_adjust_E, kwargs [ 'E_selcode' ]),
              output = kwargs [ 'E_hydrology' ])

cdo.setname ( "precip_acc",
              input = " -setunit,%s -mulc,%d  %s" % ( unit, unit_factor_EP, kwargs [ 'P_selcode' ]),
              output = kwargs [ 'P_hydrology' ])

#####
##### calculate yearly sum from yearly mean from monthly mean, weighted by days per month, days per year
#####

cdo.timmean (  input  = " -muldpy -yearmonmean %s " % ( kwargs [ 'E_hydrology' ]), 
               output = kwargs [ 'E_yearsummean' ])

cdo.timmean (  input  = " -muldpy -yearmonmean %s " % ( kwargs [ 'P_hydrology' ]), 
               output = kwargs [ 'P_yearsummean' ])

#####
#####
##### calculate field mean from yearly sum
#####

cdo.fldmean ( input  = "%s" % ( kwargs [ "E_yearsummean" ]),
              output = kwargs [ "E_ysfm" ])

cdo.fldmean ( input  = "%s" % ( kwargs [ "P_yearsummean" ]),
              output = kwargs [ "P_ysfm" ])

cdo.setname ( 'waterstorage_change',
              input  = " -sub %s %s " %\
                         ( kwargs [ "P_ysfm" ], kwargs [ "E_ysfm" ]),
              output = kwargs [ "WSC_ysfm" ])
cdo.merge ( input  = " %s %s %s " % ( kwargs [ "P_ysfm" ], kwargs [ "E_ysfm" ],
                                      kwargs [ "WSC_ysfm" ]),
            output =  kwargs [ "H_ysfm" ])


#####
##### finishing
#####

print "Preprocessing for field %s with '%s' run successfully." % ( fieldname, __file__ )


print "Exiting %s" % os.path.realpath(__file__)
