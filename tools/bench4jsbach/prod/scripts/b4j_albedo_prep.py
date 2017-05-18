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
import os

#####
##### check input parameter
#####

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':' ) for token in sys.argv [ 1: ])
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )

#####
##### calculation
#####

fieldname      = "albedo"

cdo = Cdo ()

cdo.setname ( fieldname,
              input=" -setunit,1 -ifthen %s -div %s %s " % (
        kwargs [ 'land_sea_mask' ],
        kwargs [ 'U_selcode' ],
        kwargs [ 'D_selcode' ]),
              output=kwargs [ 'UD_div' ],
              options='-f nc ')

cdo.yseasmean ( input= kwargs [ 'UD_div' ],
                output=kwargs [ 'A_seasmean' ])

#### only references

cdo.yseasmean ( input=" -ifthen %s %s " % (
        kwargs [ 'land_sea_mask' ],
        kwargs [ 'D_selcode' ]),
                output=kwargs [ 'D_seasmean' ])


##### new, only candidates

cdo.fldmean ( input=kwargs [ "A_seasmean" ],
              output=kwargs [ "A_seasmean_gm" ],
              options=' -f %s ' % kwargs [ "suffix" ] )

cdo.timmean ( input=kwargs [ "A_seasmean_gm" ],
              output=kwargs [ "A_seasmean_gtm" ],
              options=' -f %s ' % kwargs [ "suffix" ])

cdo.fldmean ( input=" -abs %s " %  kwargs [ "A_seasmean" ], 
              output=kwargs [ "A_seasmean_gd" ],
              options=' -f %s ' % kwargs [ "suffix" ] )

cdo.timmean ( input=kwargs [ "A_seasmean_gd" ],
              output=kwargs [ "A_seasmean_gtd" ],
              options=' -f %s ' % kwargs [ "suffix" ] )

#####
##### finishing
#####
    
print "\n Preprocessing for field %s with '%s' run successfully." % ( fieldname, __file__ )

print "\nExiting %s" % os.path.realpath(__file__)
