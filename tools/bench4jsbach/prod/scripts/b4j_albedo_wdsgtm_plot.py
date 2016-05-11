#! /usr/bin/env python
# -*- coding:  utf-8 -*-

import sys, os

from netCDF4 import Dataset
from datetime import datetime
from cdo import Cdo

from pycmbs.data import Data
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import matplotlib.ticker as ticker
import numpy as np
from pycmbs.mapping import map_plot

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
from analysis import get_seasonal_period_2d_coords


#print 'Number of arguments:', len(sys.argv), 'arguments.'
#print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':' ) for token in sys.argv [ 1: ])
#print kwargs
#print "subprocess '%s', kwargs = " % __file__
#print print_item_or_list_rec ( kwargs, 3 )


fieldname      = "albedo"
varname        = "delta"
yearliststring = ','.join ( repr ( y ) for y in range ( int ( kwargs [ 'startyear' ]),
                                                        int ( kwargs [ 'finalyear' ]) + 1 ))

seasons =  ['DJF','MAM','JJA','SON']

refstr = 'ref'
refnos = list ( set ([ int ( s.partition ( refstr )[ -1 ]) for s in kwargs.keys() if refstr in s ]))



colors = [ 'r', 'g', 'b', 'c', 'y', 'm' ]
colored_dots = [ 'ro', 'go', 'bo', 'co', 'yo', 'mo' ]
colored_lines = [ 'r-', 'g-', 'b-', 'c-', 'y-', 'm-' ]

hmargin_unitfactor = 1./2.
vmargin_factor = 1./15.


# 4.a.3.b,c)

#globalNtime_mean = self.diagnostic_fields [ A_varname ][ 'globalNtime_mean' ]
#assert globalNtime_mean

fig = plt.figure()

refs = []
xcoords = []
ymin = sys.float_info.max
ymax = - sys.float_info.max

for refno in refnos:

    gm_str = "albedo_wdsgm_ref%d" % refno
    tm_str = "albedo_wdsgtm_ref%d" % refno
    #globalNtime_mean_can [ 'infiles' ].append ( can.B4J_index [ gm_str ])
    nc_fid = Dataset ( kwargs [ gm_str ], 'r' )  
    time = nc_fid.variables [ 'time' ][ : ]
    albedo_gm = nc_fid.variables [ varname ][ : ] 
    nc_fid.close ()

    #globalNtime_mean_can [ 'infiles' ].append ( can.B4J_index [ tm_str ])
    nc_fid = Dataset ( kwargs [ tm_str ], 'r')
    albedo_tm = nc_fid.variables [ varname ][ : ] 
    nc_fid.close ()

    xcoords, ycoords = get_seasonal_period_2d_coords ( time, albedo_gm [ :, 0, 0 ],
                                                       hmargin_unitfactor, vmargin_factor )

    ymin = np.amin ( albedo_gm [ :, 0, 0 ]) if np.amin ( albedo_gm [ :, 0, 0 ]) < ymin  else  ymin 
    ymax = np.amax ( albedo_gm [ :, 0, 0 ]) if np.amax ( albedo_gm [ :, 0, 0 ]) > ymax  else  ymax
    plt.plot ( xcoords, ycoords, colored_lines [ refno ])
    plt.plot ( time, albedo_gm [ :, 0, 0 ],
               colored_dots [ refno ],
               label = r'$gm_{C_{%s}R_{%d}}$'% ( kwargs [ "canno" ], refno ) +\
               r', $R_{%d}$: ' % ( refno ) + kwargs [ "ref%d" % refno ])
    plt.axhline ( albedo_tm [ 0, 0, 0 ],
                  color =  colors [ refno ],
                  label = r'$tm_{C_{%s}R_{%d}}$'% ( kwargs [ "canno" ], refno ))

ymargin = ( ymax - ymin ) * vmargin_factor
plt.axis ([ xcoords [ 0 ], xcoords [ 5 ], ymin - ymargin, ymax + ymargin ])     
plt.title ( r'$gm_{C_{%s}R}$' % kwargs [ "canno" ]  + \
            " and " +\
            r'$tm_{C_{%s}R}$' % kwargs [ "canno" ],
            y = 1.08 )
plt.ylabel ( r'$\frac{W}{m^{2}}$',
             rotation = 0,
             position = ( -0.5, .96 ),
             fontsize = 20 ) 
plt.xlabel ( "season")
plt.xticks ( time, seasons )
plt.grid ( 'on' )
lgd = plt.legend ( numpoints = 1,  
                   loc = 2,                  
                   bbox_to_anchor = ( 1.05, 1 ))
plt.savefig ( kwargs [ 'albedo_wdsgtm_plot' ],
              dpi = 100, 
              bbox_extra_artists = ( lgd, ), 
              bbox_inches = 'tight' )
plt.close ()
