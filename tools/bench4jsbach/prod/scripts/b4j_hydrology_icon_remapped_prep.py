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


cdo = Cdo ()

#####
##### calculate zonal mean from yearly sum
#####

cdo.zonmean ( input  = "%s" % ( kwargs [ "E_yearsummean" ]),
              output = kwargs [ "E_ys-zonmean" ])

cdo.zonmean ( input  = "%s" % ( kwargs [ "P_yearsummean" ]),
              output = kwargs [ "P_ys-zonmean" ])

#####
##### finishing
#####

print "Preprocessing for field %s with '%s' run successfully." % ( kwargs [ "fieldname" ], __file__ )


print "Exiting %s" % os.path.realpath(__file__)
