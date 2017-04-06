#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
from cdo import Cdo
from netCDF4 import Dataset


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
output_items   =  [ 'UD_div', 'A_seasmean', 'D_seasmean' ]

D = "rsds"
D_longname = "Surface Downwelling Solar Radiation"
D_f = kwargs [ 'UD_div' ].replace ( 'UD_div', 'D', 1 )

U = "rsus"
U_f = kwargs [ 'UD_div' ].replace ( 'UD_div', 'U', 1 )
U_longname = "Surface Upwelling Solar Radiation" 

cdo = Cdo ()

cdo.setname ( D, input = " -add %s %s """ % ( kwargs [ 'Dvis_selcode' ],
                                              kwargs [ 'Dnir_selcode' ]),
              output = D_f )

nc_fid  = Dataset ( D_f, 'a' )
nc_fid.variables [ D ].setncattr ( 'standard_name', D_longname )
nc_fid.close ()

cdo.setname ( U, input = " -sub %s %s """ % ( D_f,
                                              kwargs [ 'rsns_selcode' ]),
              output = U_f )

nc_fid  = Dataset ( U_f, 'a' )
nc_fid.variables [ U ].setncattr ( 'standard_name', U_longname )
nc_fid.close ()

cdo.setname ( fieldname,
              input="-setunit,1 " +\
              "-setvrange,0.0001,0.99 " +\
              "-div " + U_f + " " + D_f,
              output=kwargs [ 'UD_div' ],
              options='-f nc -m 0.00011' )

nc_fid  = Dataset ( kwargs [ 'UD_div' ], 'a' )
nc_fid.variables [ fieldname ].setncattr ( 'standard_name', fieldname )
nc_fid.close ()

cdo.yseasmean ( input= kwargs [ 'UD_div' ],
                output=kwargs [ 'A_seasmean' ])

cdo.yseasmean ( input= D_f,
                output=kwargs [ 'D_seasmean' ])

# TODO: THIS PART HAS TO BE MOVED IN THE MAIN TOOL CODE
if "remap_targets" in kwargs and kwargs [ "remap_targets" ]:
    for remap_target in  kwargs [ "remap_targets" ]:
        if remap_target == "MPIESM1":
            gridstring = "t63grid"
            print "WARNING: fixed remapping to gaussian %s." % gridstring
            for output_item in output_items:
                output_item_remap =  output_item + "_remap-" + remap_target
                
                print output_item
                print kwargs [ output_item ]
                print  os.path.isfile ( kwargs [ output_item ])
                print output_item_remap
                print kwargs [ output_item_remap ]
                print "Start remapping os %s to %s. This might take a while ..." % \
                    ( kwargs [ output_item ], gridstring )
                sys.stdout.flush()

                #HACK TO SAVE COMPUTING TIME
                if not os.path.isfile (  kwargs [ output_item_remap ]):
                    cdo.remapycon ( gridstring,
                                    input = " -sellonlatbox,-180,180,-90,90 %s" % kwargs [ output_item ],
                                    output = kwargs [ output_item_remap ])
        else:
            print "ERROR: remapping for %s not implemented yet." % remap_target
            sys.exit ( 1 )


### closure part, check for output files


objective_tokens = output_items
if "remap_targets" in kwargs and kwargs [ "remap_targets" ]:
    for remap_target in  kwargs [ "remap_targets" ]:
        objective_tokens += [ output_item + "_remap-" + remap_target for output_item in output_items ]
	
for token in objective_tokens:
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
