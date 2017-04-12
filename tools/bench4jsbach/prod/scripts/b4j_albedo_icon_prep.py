#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
from cdo import Cdo
from netCDF4 import Dataset


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
##### calculation
#####

fieldname      = "albedo"

D = "rsds"
D_longname = "Surface Downwelling Solar Radiation"
D_f = kwargs [ 'UD_div' ].replace ( 'UD_div', 'D', 1 )

U = "rsus"
U_f = kwargs [ 'UD_div' ].replace ( 'UD_div', 'U', 1 )
U_longname = "Surface Upwelling Solar Radiation" 

cdo = Cdo ()

cdo.setname ( D, input = " -add %s %s """ % (
    kwargs [ 'Dvis_selcode' ],
    kwargs [ 'Dnir_selcode' ]),
              output = D_f )

nc_fid  = Dataset ( D_f, 'a' )
nc_fid.variables [ D ].setncattr ( 'standard_name', D_longname )
nc_fid.close ()

cdo.setname ( U, input = " -sub %s %s """ % (
    D_f, kwargs [ 'rsns_selcode' ]),
              output = U_f )

nc_fid  = Dataset ( U_f, 'a' )
nc_fid.variables [ U ].setncattr ( 'standard_name', U_longname )
nc_fid.close ()

cdo.setname ( fieldname,
              input=" -setunit,1 " +\
              " -div " + U_f + " " + D_f,
              output=kwargs [ 'UD_div' ],
              options='-f nc ')

nc_fid  = Dataset ( kwargs [ 'UD_div' ], 'a' )
nc_fid.variables [ fieldname ].setncattr ( 'standard_name', fieldname )
nc_fid.close ()

cdo.yseasmean ( input= kwargs [ 'UD_div' ],
                output=kwargs [ 'A_seasmean' ])

cdo.yseasmean ( input= D_f,
                output=kwargs [ 'D_seasmean' ])

#####
##### finishing
#####
    
print "\n Preprocessing for field %s with '%s' run successfully." % ( fieldname, __file__ )

print "\nExiting %s" % os.path.realpath(__file__)
