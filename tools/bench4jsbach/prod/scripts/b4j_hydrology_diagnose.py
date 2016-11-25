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
##### constants definition
#####
fieldname = "hydrology"
arglist = []

# arglist   = [ 'startyear',
#               'finalyear',
#               'can.P_hydrology',
#               'P_yseasmean',
#               'P_yearsum_gm' ]
# varname   = "var83"
# varname_2 = "precip_acc"

#####
##### check input parameter
#####
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':' ) for token in sys.argv [ 1: ])

print "subprocess '%s', kwargs = " % __file__
print arglist
print print_item_or_list_rec ( kwargs, 3 )

#assert set ( arglist ).issubset ( set ( kwargs.keys ()))

#####
##### calculation
#####

cdo = Cdo ()

cdo.sub ( input  = "%s %s" % ( kwargs [ "can.E_yearsummean" ], kwargs [ "ref.E_yearsummean" ]),
          output = kwargs ["E_yearsumdiff"  ])

cdo.sub ( input  = "%s %s" % ( kwargs [ "can.P_yearsummean" ], kwargs [ "ref.P_yearsummean" ]),
          output = kwargs ["P_yearsumdiff"  ])



#####
##### finishing
#####
print "Preprocessing for field %s with '%s' run successfully." % ( fieldname, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
