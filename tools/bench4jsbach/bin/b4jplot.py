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
import cartopy as cart
import numpy as np
from pandas import DataFrame
import math
from pprint import pprint
from util import print_item_or_list_rec

# from matplotlib import rcParams
# rcParams.update({'figure.autolayout': True})


colorlist = [ 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w' ]
    
darkbrown  = '#c6903f' 
lightbrown = '#e6cd92'
lightgreen = '#ceece8'
darkgreen  = '#6fc1b6'
    
colors_hydrology  = [ darkbrown, lightgreen, lightbrown, darkgreen ]

    
def mapplot ( inpath, varname, outpath, unit_TeX, title_TeX, midpoint = np.nan, vmin = np.nan, vmax = np.nan ):
    #from datetime import datetime

    print "\n\ninpath:", inpath
    print "varname:", varname
    print "outpath:", outpath
    print "CAUTION: tstep_no := 0 hardcoded"

    tstep_no = 0

    nc_fid  = Dataset ( inpath, 'r')

    var       = nc_fid.variables [ varname ]
    units     = var.units
    grid_type = var.grid_type
#    long_name = var.long_name
    time      = nc_fid.variables [ 'time' ][ : ]

    if grid_type == "gaussian":
        lats    = nc_fid.variables [ 'lat' ][ : ]
        lons    = nc_fid.variables [ 'lon' ][ : ] 
        vardata = var [ tstep_no, :, : ]
    elif grid_type == "unstructured":
        
        rad2deg_ = 180 / math.pi
        clon, clat =  nc_fid.variables['clon'][:] * rad2deg_,\
                      nc_fid.variables['clat'][:] * rad2deg_
        clon_vertices, clat_vertices  =  nc_fid.variables['clon_bnds'][:] * rad2deg_,\
                                         nc_fid.variables['clat_bnds'][:] * rad2deg_
        vardata      =  var [ tstep_no, : ] 
    else:
        print "ERROR: grid_type not found."
        sys.exit ( 1 )
        
    nc_fid.close ()


    if grid_type == "unstructured":
        ncells, nv = clon_vertices.shape
        mt_wholeglobe_r = (( np.amax (clon_vertices, axis = 1 ) - np.amin (clon_vertices, axis = 1 )) > 180 ).astype ( int )
        mt_wholeglobe_mask = np.array ([ mt_wholeglobe_r, mt_wholeglobe_r, mt_wholeglobe_r ]).T               
        east = ( clon_vertices < 0 ).astype ( int )
        east_weight = np.sum ( mt_wholeglobe_mask * east, axis=1 )
        east_weight_mask = np.array ([ east_weight,  east_weight,  east_weight ]).T
        shift_eastwest = ( east_weight_mask == 1 ).astype ( int ) * ( east == 1 ).astype ( int ) * 360 + \
                         ( east_weight_mask == 2 ).astype ( int ) * ( east == 0 ).astype ( int ) * ( -360 )
        newxvc = ( clon_vertices + shift_eastwest ).ravel()
        newyvc = clat_vertices.ravel()
        triangles = np.reshape ( np.arange ( ncells * nv ), ( ncells, nv ))

    
    vmax    = vardata.max () if np.isnan ( vmax ) else vmax
    vmin    = vardata.min () if np.isnan ( vmin ) else vmin
    vmean   = vardata.mean () if np.isnan ( midpoint ) else midpoint   
    print "vmin = %0.4g, vmean = %0.4g, vmax = %0.4g" % ( vmin, vmean, vmax )
    
# # arguments
    nlevs                = 20
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

    #print "####### COLORS: VALS, LEVELS: #########"
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
    ax.set_xlim ([ -180, 180 ])
    ax.set_ylim ([ -90, 90 ])
    
    if grid_type == 'gaussian':
        # correct grid with half distance to north and west
        lons_shifted = lons - np.concatenate (( np.zeros ( 1 ), np.abs ( lons[0:-1] - lons [1:] ) / 2 ), axis = 0 )
        lats_shifted = lats + np.concatenate (( np.abs ( lats [0:-1] - lats [1:] ) / 2, np.zeros ( 1 )), axis = 0 )
        #ax.pcolormesh ( lons, lats,
        ax.pcolormesh ( lons_shifted, lats_shifted,
                        vardata, cmap=cmap, norm=norm,
                        edgecolors='face', 
                        transform=ccrs.PlateCarree(), 
                        rasterized=True, zorder=0 )
    elif grid_type == 'unstructured':
        plt.tripcolor( newxvc, newyvc, triangles,
                       facecolors=vardata, cmap=cmap, norm=norm,
                       edgecolor="face",
                       transform=ccrs.PlateCarree ())
        
    ax.coastlines ( linewidth=0.5 )
    ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth=0.5, facecolor="silver" )

    print f
   
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

    # fp = open ( "runtimefigureinfo.txt", 'a' )
    # fp.write ( "\n\n\n########################" )
    # fp.write ( "\ninpath: %s " % inpath )
    # fp.write ( "\noutpath: %s " % outpath )
    # fp.write ( "\nFigure: %s " % type ( f ))
    # fp.write ( "\nsize: %s" % f.get_size_inches ())
    # d = f.__dict__
    # fp.write ("""\n%s""" % print_item_or_list_rec ( d, 3 ))
    # fp.close ()
    
    #f.tight_layout()
    plt.savefig ( outpath, dpi = 300 )
    plt.close ()
  


def zon_plot ( varname, unit, title, inpathlist, legendlist, 
               outpath = None, 
               ax = None,
               art = None,
               mono = False ):  
  
    nrefargs = len ( inpathlist )   
    
    for refarg_no in range ( nrefargs ):
        nc_fid = Dataset ( inpathlist [ refarg_no ])  
        lat = nc_fid.variables [ 'lat'][:]
        vardat_ref = nc_fid.variables [ varname ][ : ] 
        nc_fid.close ()
        ax.plot ( lat, vardat_ref [ 0,:, 0 ],
                   colorlist [ refarg_no ],
                   label = legendlist [ refarg_no ] if mono == False else None )
        ax.axhline ( vardat_ref [ 0,:, 0 ].mean(), color = colorlist [ refarg_no ])   

        
    ax.set_xticks ( np.arange ( -90, 91, 15 ))
    ax.yaxis.set_label_coords ( -.09, .88 ) 
    ax.set_ylabel ( "[ %s ]" % unit,
                    rotation = 90,
                    fontsize = 12 )
    ax.set_xlabel ( "latitude [ deg ]",
                    position = ( 0.9, 0 ),
                    fontsize = 12 )
    ax.set_title ( title, y = 1.04, fontsize = 14 )
    ax.grid ( 'on' )
    lgd = ax.legend ( numpoints = 1,  
                     loc = 2, 
                     bbox_to_anchor = ( 1.05, 1 ) if mono == False else None,
                     bbox_transform=ax.transAxes if mono == False else None )
    if isinstance ( art, list ) and lgd: art.append ( lgd )


def zon_plot_wrapper ( varname, unit, can_name, canno, 
                       input_can, subtitle_can,
                       input_delta, legend_delta, subtitle_delta,
                       outpath ):

    f = plt.figure ( figsize = ( 10, 15 ))
    show_can = True
    supgrid = gridspec.GridSpec( 2 if show_can == True else 1, 1 ) 
    art = []

    ####### first subplot
    subgrid = gridspec.GridSpecFromSubplotSpec ( 1, 1, subplot_spec = supgrid [ 0 ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax )
    
    # first subplot args
    subtitle   = subtitle_can
    inpathlist = [ input_can ]
    legendlist = []

    # first subplot call
    zon_plot ( varname, unit, subtitle, inpathlist, legendlist, 
               ax = ax,
               art = art,
               mono = True )
    
    ####### second subplot
    subgrid = gridspec.GridSpecFromSubplotSpec (1, 1, subplot_spec = supgrid [ 1 ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax ) 

    # first subplot args
    subtitle      = subtitle_delta
    inpathlist = input_delta
    legendlist = legend_delta

    # first subplot call
    zon_plot ( varname, unit, subtitle, inpathlist, legendlist, 
               ax = ax,
               art = art )

    plt.savefig ( outpath,
                  bbox_extra_artists = ( art [ 0 ], ) if art else None,
                  bbox_inches = 'tight' if art else None  )

    plt.close () 


def hist_plot (  colno, colcolors, colnames,
                 alldat, unit, index, title, 
                 outpath = None, 
                 ax = None,
                 art = None,
                 mono = False,
                 solitaire = False,
                 legend = True,
                 ylabel = True ):
    
    df = DataFrame ( alldat, 
                     columns = colnames, 
                     index = index )

    ax_loc = df.plot ( ax = ax, kind = 'bar', color = colcolors, rot = 0,
                       legend = legend )
    ax_loc.set_title ( title, y = 1.02 )
    ax_loc.xaxis.set_label_coords ( 1, -0.1 )
    if ylabel == True:
        ax_loc.yaxis.set_label_coords ( -0.2, .92 ) 
        ax_loc.set_ylabel ( "[ %s ]" % unit, rotation = 90 )
    
    ax_loc.grid ( 'on')
    ax_loc.axhline ( 0, color = 'k' )
    if legend == True:
        lgd = ax_loc.legend ( numpoints = 1,  
                              loc = 2,                  
                              bbox_to_anchor = ( 1.05, 1 ))
        
        if isinstance ( art, list ) and lgd: art.append ( lgd )
        
    if solitaire == True and outpath:
        plt.savefig ( outpath,
                      dpi = 100, 
                      bbox_extra_artists = ( lgd, ), 
                      bbox_inches = 'tight' )
        plt.close ()


def hist_plot_wrapper ( colno, colcolors, colnames,
                        varnames, unit, can_name, canno,
                       input_can, subtitle_can,
                       input_delta, subtitle_delta, xlabel_delta,
                       outpath ):
    print "in hist_plot_wrapper ()"

    show_can = True

    ndeltas = len ( input_delta )

    f        = plt.figure ( figsize = ( 10, 6 ))
    supgrid  = gridspec.GridSpec( 1, ndeltas if show_can == True else 1 ) 
    art      = []

    ####### first subplot
    subgrid = gridspec.GridSpecFromSubplotSpec ( 1, 1, subplot_spec = supgrid [ 0 ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax )

    accum = np.array ([[]])
    nc_fid = Dataset ( input_can )
    vardat = np.concatenate (( nc_fid.variables [ varnames [ 0 ]] [ : ][ 0 ],
                               nc_fid.variables [ varnames [ 1 ]] [ : ][ 0 ],
                               nc_fid.variables [ varnames [ 2 ]] [ : ][ 0 ],
                               nc_fid.variables [ varnames [ 3 ]] [ : ][ 0 ]),
                               axis = 1 )
    nc_fid.close ()
    accum = np.concatenate (( accum, vardat ), axis = 0 ) if accum != [[]] else vardat
    hist_plot (  colno, colcolors, colnames,
                 accum, unit, [ ""], subtitle_can,
                 ax = ax,
                 legend = False )
    

                      

    ####### second subplot
    subgrid = gridspec.GridSpecFromSubplotSpec ( 1, 1, subplot_spec = supgrid [ 1: ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax )
    
 # second subplot args
    nrefs = len ( input_delta )
    accum = np.array ([[]])
    for k in range ( nrefs ):
        nc_fid = Dataset ( input_delta [ k ])
        vardat = np.concatenate (( nc_fid.variables [ varnames [ 0 ]] [ : ][ 0 ],
                                   nc_fid.variables [ varnames [ 1 ]] [ : ][ 0 ],
                                   nc_fid.variables [ varnames [ 2 ]] [ : ][ 0 ],
                                   nc_fid.variables [ varnames [ 3 ]] [ : ][ 0 ]),
                                 axis = 1 )
        nc_fid.close ()
        accum = np.concatenate (( accum, vardat ), axis = 0 ) if accum != [[]] else vardat

    hist_plot (  colno, colcolors, colnames,
                 accum, unit, xlabel_delta, subtitle_delta,
                 ax = ax,
                 art = art,
                 ylabel = not show_can )
    
    plt.savefig ( outpath,
                  bbox_extra_artists = ( art [ 0 ], ) if art else None,
                  bbox_inches = 'tight' if art else None  )

    plt.close ()


def get_globalmean ( filepath, varname ):
    nc_fid = Dataset ( filepath )  
    lat = nc_fid.variables [ 'lat'][:]
    vardat = nc_fid.variables [ varname ][ : ] 
    nc_fid.close ()
    return vardat [ 0,:, 0 ].mean()
