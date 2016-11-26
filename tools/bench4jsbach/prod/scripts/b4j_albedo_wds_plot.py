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

kwargs        =  dict ( token.split ( ':', 1 ) for token in sys.argv [ 1: ])
print "\nsubprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )

fieldname      = kwargs [ "fieldname" ] 


#####
##### constants definition
#####
varname        = "delta"
seasons =  ['DJF','MAM','JJA','SON']


#####
##### plotting
#####
nc_fid = Dataset ( kwargs [ 'albedo_wds' ], 'r')
albedo_delta = nc_fid.variables [ varname ][ : ] 
nc_fid.close ()
vmax = albedo_delta.max ()
vmin = albedo_delta.min ()


# # arguments
nlevs                = 20
tick_interleave_no   = int ( nlevs // 10 )
midval               =  0

# calculate colormap, colorbar, ticks
assert nlevs        >=  5       
unit                 = ( vmax - vmin )  / ( nlevs - 1 )
baselevel_left_val   = midval - unit/2
baselevel_right_val  = midval + unit/2


assert vmin < baselevel_left_val and \
        baselevel_left_val < midval and \
        midval < baselevel_right_val and  \
        baselevel_right_val < vmax
nlevs_left           = int (( baselevel_left_val - vmin ) // unit ) + \
                         ( 1 if int (( baselevel_left_val - vmin ) % unit ) == 0 else 0 )

nlevs_right          = int (( vmax - baselevel_right_val ) // unit ) + \
                         ( 1 if int (( vmax - baselevel_right_val ) % unit ) == 0 else 0 )
nlevs_revised        = nlevs_left + nlevs_right + 1
leftlevel_left_val   = baselevel_left_val -  nlevs_left * unit
rightlevel_right_val = baselevel_right_val + nlevs_right * unit
levels               = np.arange ( leftlevel_left_val, rightlevel_right_val + unit/2, unit )
ticks_left           = np.arange ( midval, leftlevel_left_val, - unit ) [::-1]
ticks_right          = np.delete ( np.arange ( midval, rightlevel_right_val, unit ), 0 )
ticks                = np.concatenate (( ticks_left, ticks_right ), axis = 0 )
ticks_interleave     = np.take ( ticks, np.arange ( nlevs_left % tick_interleave_no,
                                                    len ( ticks ) - 1, tick_interleave_no ))
vals                 = np.interp ( ticks, [ leftlevel_left_val, midval, rightlevel_right_val ],
                                   [ 0, 0.5, 1 ])
colors               = plt.cm.RdBu ( vals )
cmap, norm           = mcolors.from_levels_and_colors ( levels, colors )

cdo = Cdo ()

dates = cdo.showtimestamp ( input =  kwargs [ 'albedo_wds' ])[ 0 ].split ( " ")
dates = [ datetime.strptime ( str ( x ), '%Y-%m-%dT%H:%M:%S') for x in dates if x ]
# dkl: hack, problem with minutes of datetime.datetime in pycmbs data
dates = [ date.replace ( minute=0 ) for date in dates ]

f = plt.figure ( figsize = ( 16, 10 ))
ax = [ f.add_subplot ( 2, 2, i + 1 ) for i in range ( 4 )]
d = [ None ] * 4                
gs1 = gridspec.GridSpec ( 2, 2 )
gs1.update ( wspace = -0.08, hspace = -0.04  )               
for i in range ( 4 ):
    ax [i] = plt.subplot ( gs1 [ i ])
    d [i] = Data ( kwargs [ 'albedo_wds' ],
                   "delta",
                   read       = True,
                   start_time = dates [ i ],
                   label      = seasons [ i ] )    
    x = map_plot ( d[i],
                   ax                  = ax[i],
                   vmin                = leftlevel_left_val,
                   vmax                = rightlevel_right_val,
                   use_basemap         = True,
                   show_colorbar       = False,
                   return_plot_handler = True,
                   cmap_data           = cmap )

cax = f.add_axes ([ 0.17, 0.07, 0.7, 0.03])
cb = f.colorbar ( x.im, cax = cax, orientation = "horizontal",
                  ticks = ticks_interleave,
                  format = '%.2e' )#ticker.FuncFormatter ( '%.2e' ))#'%.3f' )
cb.ax.text ( 0.97, 2.5, r'$\frac{W}{m^{2}}$', size = 32 )
if 'plot_title' in kwargs:
    plt.suptitle ( kwargs [ 'plot_title' ], fontsize = 24 )
plt.savefig ( kwargs [ 'albedo_wds_plot' ])
plt.close ()

# delta_can = delta [ 'cans' ][ '%d' % can.B4J_canno ]
# assert delta_can
# delta_can [ 'infiles' ].append ( can.B4J_index [ delta_str ])
# delta_can [ 'plots' ].append ( can.B4J_index [ plot_delta_str ])


#####
##### finishing
#####
print "Exiting %s" % os.path.realpath(__file__)
