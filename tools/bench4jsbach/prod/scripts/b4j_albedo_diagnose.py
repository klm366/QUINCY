#! /usr/bin/env python
# -*- coding:  utf-8 -*-

import sys, os
from cdo import Cdo

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

#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':' ) for token in sys.argv [ 1: ])
#print "subprocess '%s', kwargs = " % __file__
#print print_item_or_list_rec ( kwargs, 3 )


fieldname      = "albedo"
yearliststring = ','.join ( repr ( y ) for y in range ( int ( kwargs [ 'startyear' ]),
                                                        int ( kwargs [ 'finalyear' ]) + 1 ))

cdo = Cdo ()

cdo.setname ( "delta",
              input=" -setunit,'W m-2' -setcode,'9999' -mul -sub " +\
              kwargs [ "can.A_seasmean" ] + " " +\
              kwargs [ "ref.A_seasmean" ] + " " +\
              kwargs [ "ref.D_seasmean" ],
              output=kwargs [ "albedo_wds" ],
              options='-f %s' % kwargs [ "suffix" ])

cdo.fldmean ( input="-setmisstoc,0 -mul " +\
              kwargs [ "albedo_wds" ] + " " +\
              kwargs [ "land_sea_mask" ],
              output=kwargs [ "albedo_wdsgm" ],
              options='-f %s' % kwargs [ "suffix" ] )

cdo.timmean ( input=kwargs [ "albedo_wdsgm" ],
              output=kwargs [ "albedo_wdsgtm" ],
              options='-f %s' % kwargs [ "suffix" ])
        
cdo.fldmean ( input="-abs -setmisstoc,0 -mul " +\
              kwargs [ "albedo_wds" ] + " " +\
              kwargs [ "land_sea_mask" ],
              output=kwargs [ "albedo_wdsgd" ],
              options='-f %s' % kwargs [ "suffix" ] )

cdo.timmean ( input=kwargs [ "albedo_wdsgd" ],
              output=kwargs [ "albedo_wdsgtd" ],
              options='-f %s' % kwargs [ "suffix" ] )
