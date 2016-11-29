#! /usr/bin/env python
# -*- coding:  utf-8 -*-
"""The module contains b4j special plotting methods."""


import sys, os

from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs
import numpy as np
    
def mapplot ( inpath, varname, outpath, unit_TeX, title_TeX, midpoint = np.nan ):
    #from datetime import datetime

    print "inpath:", inpath
    print "varname:", varname
    print "outpath:", outpath

    nc_fid  = Dataset ( inpath, 'r')
    lats    = nc_fid.variables [ 'lat' ][ : ]
    lons    = nc_fid.variables [ 'lon' ][ : ]
    time    = nc_fid.variables [ 'time' ][ : ]   
    vardata = nc_fid.variables [ varname ][ 0,:,: ] 
    nc_fid.close ()
    
    vmax    = vardata.max ()
    vmin    = vardata.min ()
    vmean   = vardata.mean () if np.isnan ( midpoint ) else midpoint   
    print "vmin = %0.4g, vmean = %0.4g, vmax = %0.4g" % ( vmin, vmean, vmax )
    
# # arguments
    nlevs                = 50
    tick_interleave_no   = int ( nlevs // 10 )
    midval               =  vmean

# calculate colormap, colorbar, ticks
    assert nlevs        >=  5   
    vrange               = vmax - vmin if vmax > 0 else abs ( vmin - vmax )
    unit                 = vrange  / ( nlevs - 1 )
    baselevel_left_val   = midval - unit/2
    baselevel_right_val  = midval + unit/2   
   # print "vrange = %0.4g, unit = %0.4g, baselevel_left_val = %0.4g, baselevel_right_val = %0.4g" % \
    #    ( vrange, unit, baselevel_left_val, baselevel_right_val )


 #   assert vmin < baselevel_left_val and \
   #     baselevel_left_val < midval and \
    #    midval < baselevel_right_val and  \
     #   baselevel_right_val < vmax
    nlevs_left           = int (( baselevel_left_val - vmin ) // unit ) + 1
    nlevs_right          = int (( vmax - baselevel_right_val ) // unit ) + 1
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
    
 #   print "leftlevel_left_val, midval, rightlevel_right_val ", \
  #      leftlevel_left_val, midval, rightlevel_right_val
    
    colors      = plt.cm.BrBG ( vals )
    cmap, norm  = mcolors.from_levels_and_colors ( levels, colors )
    norm.clip   = True
    norm.vmin   = vmin
    norm.vmax   = vmax
    cmap.set_bad ( color = 'silver', alpha = 1.0 )
    cmap.set_over ( color = 'r', alpha = 1.0 )
    cmap.set_under ( color = 'g', alpha = 1.0 )
    
    f = plt.figure ( figsize = ( 16, 8 ))   
    ax = f.add_subplot ( 111, projection = ccrs.PlateCarree ())
    ax.pcolormesh ( lons, lats, vardata, cmap=cmap, norm=norm,
                    edgecolors='face', transform=ccrs.PlateCarree(), 
                    rasterized=True, zorder=0 )
    ax.coastlines ( linewidth=0.5 )
   
    cax = f.add_axes ([ 0.125, 0.04, 0.775, 0.03 ]) # left, bottom, width, height           
    cb = mpl.colorbar.ColorbarBase ( cax, cmap = cmap, norm = norm,
                                     orientation = 'horizontal', extend = 'neither', 
                                     ticks = ticks_interleave, format = '%.2e' )
    for t in cb.ax.get_xticklabels ():
        t.set_fontsize ( 14 )    
    cb.ax.text ( 0.89, 1.5, "[%s]" % unit_TeX, size = 16 )
    plt.text ( 0.5, 1.04, title_TeX,
               horizontalalignment = 'center',
               fontsize = 20,
               transform = ax.transAxes )
    
    plt.savefig ( outpath ) 
    plt.close ()


def zon_plot ( title, varname, inpathlist, unit, legendlist, outpath ):

    print "start zon_plot"
    print title
    print varname

    colorlist = [ 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w' ]
    print inpathlist
    print type ( inpathlist )
    nrefargs = len ( inpathlist )
    print legendlist
    print outpath
    
    assert len ( legendlist ) == nrefargs and \
        len ( colorlist ) >= nrefargs

    fig = plt.figure ()
    ax = fig.gca()

    for refarg_no in range ( nrefargs ):
        nc_fid = Dataset ( inpathlist [ refarg_no ])  
        lat = nc_fid.variables [ 'lat'][:]
        vardat_ref = nc_fid.variables [ varname ][ : ] 
        nc_fid.close ()
        plt.plot ( lat, vardat_ref [ 0,:, 0 ],
                   colorlist [ refarg_no ],
                   label = legendlist [ refarg_no ])
        plt.axhline ( vardat_ref [ 0,:, 0 ].mean(), color = colorlist [ refarg_no ])

    
    ax.set_xticks ( np.arange ( -90, 90, 15 ))
    ax.yaxis.set_label_coords ( -0.2, 1 ) 
    plt.ylabel ( "[ %s ]" % unit, rotation = 0 )
    ax.xaxis.set_label_coords ( 1.05, - 0.1 ) 
    plt.xlabel ( "latitude [ deg ]")
    plt.title ( title, y = 1.08 )
    plt.grid ( 'on' )
    lgd = plt.legend ( numpoints = 1,  
                       loc = 2,                  
                       bbox_to_anchor = ( 1.05, 1 ))
    plt.savefig ( outpath,
                  dpi = 100, 
                  bbox_extra_artists = ( lgd, ), 
                  bbox_inches = 'tight' ) 
    plt.close ()
    
