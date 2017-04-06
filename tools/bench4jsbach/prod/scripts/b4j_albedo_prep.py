#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
from cdo import Cdo

print "Starting %s" % os.path.realpath(__file__)


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

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':' ) for token in sys.argv [ 1: ])
print "subprocess '%s', kwargs = " % __file__
if "remap_targets" in kwargs:
    if len ( kwargs [ "remap_targets" ]) > 2:
        kwargs [ "remap_targets" ] = kwargs [ "remap_targets" ][1:-1].split(",")
    else:
        kwargs [ "remap_targets" ] = None

    print kwargs [ "remap_targets" ]
    print "*************************************"
print print_item_or_list_rec ( kwargs, 3 )

### main part

fieldname      = "albedo"
yearliststring = ','.join ( repr ( y ) for y in range ( int ( kwargs [ 'startyear' ]),
                                                        int ( kwargs [ 'finalyear' ]) + 1 ))

cdo = Cdo ()

cdo.setname ( fieldname,
              input="-setunit,1 " +\
              "-setvrange,0.0001,0.99 " +\
              "-div " +\
              ( " -selyear," + yearliststring + \
                " " +  kwargs [ 'U_selcode' ]) +\
              ( " -selyear," + yearliststring + \
                " " +  kwargs [ 'D_selcode' ]),
              output=kwargs [ 'UD_div' ],
              options='-f nc -m 0.00011' )

cdo.yseasmean ( input="-selyear," + yearliststring + \
                " " + kwargs [ 'UD_div' ],
                output=kwargs [ 'A_seasmean' ])

cdo.yseasmean ( input="-selyear," + yearliststring + \
                " " + kwargs  [ 'D_selcode' ],
                output=kwargs [ 'D_seasmean' ])


### closure part, check for output files
	
for token in [ 'UD_div', 'A_seasmean', 'D_seasmean' ]:
    if os.path.isfile ( kwargs [ token ]):
        try:
            filetype = cdo.filedes ( input=kwargs [ token ])
        except:
            filetype = "unsupported filetype"
        else:
            filetype = filetype[1]

        print "\nFile '%s' written at %s" % ( token, kwargs [ token ])
        print "size %s byte, type %s" % ( os.stat ( kwargs [ token ]).st_size, filetype )
    
    else:
        print "\nERROR: File '%s' not written at %s" % ( token, kwargs [ token ])
        sys.exit ( 1 )
    
    print "\n Preprocessing for field %s with '%s' run successfully." % ( fieldname, __file__ )

print "\nExiting %s" % os.path.realpath(__file__)
