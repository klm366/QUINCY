#! /usr/bin/env python
# -*- coding:  utf-8 -*-
"""The module contains b4j special plotting methods."""


import sys, os

from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import addcyclic, shiftgrid
import cartopy.crs as ccrs
import cartopy as cart
import numpy as np
from pandas import DataFrame
import math
from pprint import pprint

from diagnostics import get_seasonal_period_2d_coords
from util import print_item_or_list_rec


colorlist = [ 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w' ]
    
darkbrown  = '#c6903f' 
lightbrown = '#e6cd92'
lightgreen = '#ceece8'
darkgreen  = '#6fc1b6'
    
colors_hydrology  = [ darkbrown, lightgreen, lightbrown, darkgreen ]

seasonlabel =  ['DJF','MAM','JJA','SON']
nseas       = 4

def triangles_prep ( clon_vertices, clat_vertices, vardata ):
    print "start triangles_prep"
    ncells, nv = clon_vertices.shape

    # shift vertices of cells covering half the globe
    mt_wholeglobe_r = (( np.amax (clon_vertices, axis = 1 ) - np.amin (clon_vertices, axis = 1 )) > 180 ).astype ( int )
    mt_wholeglobe_mask = np.array ([ mt_wholeglobe_r, mt_wholeglobe_r, mt_wholeglobe_r ]).T               
    east = ( clon_vertices < 0 ).astype ( int )
    east_weight = np.sum ( mt_wholeglobe_mask * east, axis=1 )
    east_weight_mask = np.array ([ east_weight,  east_weight,  east_weight ]).T
    shift_eastwest = ( east_weight_mask == 1 ).astype ( int ) * ( east == 1 ).astype ( int ) * 360 + \
                     ( east_weight_mask == 2 ).astype ( int ) * ( east == 0 ).astype ( int ) * ( -360 )
    newxvc = clon_vertices + shift_eastwest

    # copy east/west margin cells
    eastmargincells = ( newxvc > 180  ).astype ( int )
    eastmargincells_rowid, margincells_colid =  eastmargincells.nonzero ()
    eastmargincells_rowid = np.unique ( eastmargincells_rowid )
    eastnnewcells = eastmargincells_rowid.shape [ 0 ]

    westmargincells = ( newxvc < -180 ).astype ( int )
    westmargincells_rowid, margincells_colid =  westmargincells.nonzero ()
    westmargincells_rowid = np.unique ( westmargincells_rowid )
    westnnewcells = westmargincells_rowid.shape [ 0 ]

    nallcells = ncells + eastnnewcells + westnnewcells

    # calculate and concatenate new triangles
    eastnewrows = np.take ( newxvc, eastmargincells_rowid, axis = 0 )
    westnewrows = np.take ( newxvc, westmargincells_rowid, axis = 0 )
    latestxvc = np.concatenate (( newxvc,
                                  eastnewrows - 360,
                                  westnewrows + 360 ), axis = 0 )
    latestyvc = np.concatenate (( clat_vertices, 
                                  np.take ( clat_vertices, eastmargincells_rowid, axis = 0 ), 
                                  np.take ( clat_vertices, westmargincells_rowid, axis = 0)),
                                axis = 0 )
    latestvardata = np.ma.concatenate (( vardata, 
                                         np.take ( vardata, eastmargincells_rowid, axis = 0),
                                         np.take ( vardata, westmargincells_rowid, axis = 0 )), 
                                       axis = 0 )

    newxvc    = latestxvc.ravel ()
    newyvc    = latestyvc.ravel ()
    triangles = np.reshape ( np.arange ( nallcells * nv ), ( nallcells, nv ))

    return newxvc, newyvc, triangles, latestvardata


def colormap_prep ( vmin, vmax, vmean, nlev ):
    
    tick_interleave_no   = int ( nlev // 10 )
    midval               = vmean

# calculate colormap, colorbar, ticks
    assert nlev        >=  5   
    vrange               = vmax - vmin if vmax > 0 else abs ( vmin - vmax )
    unit                 = vrange  / ( nlev - 1 )
    baselevel_left_val   = midval - unit/2
    baselevel_right_val  = midval + unit/2
    
    nlev_left           = int (( baselevel_left_val - vmin ) // unit ) + 1
    nlev_right          = int (( vmax - baselevel_right_val ) // unit ) + 1
    nlev_revised        = nlev_left + nlev_right + 1
    leftlevel_left_val   = baselevel_left_val -  nlev_left * unit
    rightlevel_right_val = baselevel_right_val + nlev_right * unit
    level               = np.arange ( leftlevel_left_val, rightlevel_right_val + unit/2, unit )
    ticks_left           = np.arange ( midval, leftlevel_left_val, - unit ) [::-1]
    ticks_right          = np.delete ( np.arange ( midval, rightlevel_right_val, unit ), 0 )
    ticks                = np.concatenate (( ticks_left, ticks_right ), axis = 0 )
    ticks_interleave     = np.take ( ticks, np.arange ( nlev_left % tick_interleave_no,
                                                    len ( ticks ) - 1, tick_interleave_no ))
    values                 = np.interp ( ticks, [ leftlevel_left_val, midval, rightlevel_right_val ],
                                   [ 0, 0.5, 1 ])

    return values, level, ticks_interleave

    
def mapplot_seasonal ( varname    = None,
                       inpath     = None,
                       outpath    = None,
                       unit       = None,
                       title      = None,
                       midpoint   = None,
                       vmin       = None,
                       vmax       = None,
                       colormap   = "BrBG" ):
    
    midpoint = float ( midpoint ) if midpoint else np.nan
    vmin     = float ( vmin )     if vmin     else np.nan
    vmax     = float ( vmax )     if vmax     else np.nan

    nc_fid  = Dataset ( inpath, 'r')

    var       = nc_fid.variables [ varname ]
    units     = var.units
    grid_type = var.grid_type
    time      = nc_fid.variables [ 'time' ][ : ]

    if grid_type == "gaussian":
        lats          = nc_fid.variables [ 'lat' ][ : ]
        lons_orig     = nc_fid.variables [ 'lon' ][ : ] 
        vardata_orig  = var [ :, :, : ]
        vardata, lons = addcyclic ( vardata_orig, lons_orig )        
        vardata, lons = shiftgrid ( 180., vardata, lons, start = False )
    elif grid_type == "unstructured":
        
        rad2deg_ = 180 / math.pi
        clon, clat =  nc_fid.variables['clon'][:] * rad2deg_,\
                      nc_fid.variables['clat'][:] * rad2deg_
        clon_vertices, clat_vertices  =  nc_fid.variables['clon_bnds'][:] * rad2deg_,\
                                         nc_fid.variables['clat_bnds'][:] * rad2deg_
        vardata      =  var [:, : ]
    else:
        print "ERROR: grid_type not found."
        sys.exit ( 1 )

    nc_fid.close ()

    assert time.shape [ 0 ] == nseas
    
    vmax    = vardata.max () if np.isnan ( vmax ) else vmax
    vmin    = vardata.min () if np.isnan ( vmin ) else vmin
    vmean   = vardata.mean () if np.isnan ( midpoint ) else midpoint
    
    # #print "####### COLORS: VALUES, LEVEL: #########"
    nlev        = 20
    values, level, ticks = colormap_prep ( vmin, vmax, vmean, nlev )
    colors      = plt.cm.get_cmap ( colormap )( values )
    cmap, norm  = mcolors.from_levels_and_colors ( level, colors )
    norm.clip   = True
    norm.vmin   = vmin
    norm.vmax   = vmax
    
    cmap.set_bad   ( color = 'silver', alpha = 1.0 )
    cmap.set_over  ( color = 'r', alpha = 1.0 )
    cmap.set_under ( color = 'g', alpha = 1.0 )
    
    f = plt.figure ( figsize = ( 16, 9 )) 
    ax = [ f.add_subplot ( 2, 2, i + 1, projection = ccrs.PlateCarree () ) \
           for i in range ( nseas )]
    
    for i in range ( 4 ):
        ax [ i ].set_title ( seasonlabel [ i ], fontsize=16 )
        ax [ i ].set_xlim ([ -180, 180 ])
        ax [ i ].set_ylim ([ -90, 90 ])
        ax [ i ].coastlines ( linewidth=0.5 )
        ax [ i ].add_feature ( cart.feature.OCEAN, zorder=100, edgecolor='k',
                               linewidth=0.5, facecolor="silver" )

        if grid_type == 'gaussian':
            # correct grid with half distance to north and west
            lons_shifted = lons - np.concatenate (( np.zeros ( 1 ),
                                                    np.abs ( lons[0:-1] - lons [1:] ) / 2 ),
                                                  axis = 0 )
            lats_shifted = lats + np.concatenate (( np.abs ( lats [0:-1] - lats [1:] ) / 2,
                                                    np.zeros ( 1 )),
                                                  axis = 0 )
            ax [ i ].pcolormesh ( lons_shifted, lats_shifted,
                                  vardata [ i, :, : ], cmap=cmap, norm=norm,
                                  edgecolors='face',  
                                  rasterized=True, zorder=0 )
        elif grid_type == 'unstructured':
            newxvc, newyvc, triangles, latestvardata = triangles_prep (
                clon_vertices, clat_vertices, vardata [ i, : ])
            ax [ i ].tripcolor( newxvc, newyvc, triangles,
                                facecolors=latestvardata, cmap=cmap, norm=norm,
                                edgecolor="face")       
   
    f.subplots_adjust ( wspace = 0.05, hspace = 0.1 )#, bottom= 0.15 )
    cax = f.add_axes ([ 0.125, 0.058, 0.775, 0.03 ]) # left, bottom, width, height
    cb = mpl.colorbar.ColorbarBase ( cax, cmap = cmap, norm = norm,
                                     orientation = 'horizontal', extend = 'neither', 
                                     ticks = ticks, format = '%.2e' )
    for t in cb.ax.get_xticklabels ():
        t.set_fontsize ( 14 )    
    #cb.ax.text ( 0.98, -1.3, "[%s]" % unit, size = 20 )
    cb.ax.text ( 0.89, 1.3, "[%s]" % unit, size = 16 )
    plt.suptitle ( title, fontsize = 24 )
    plt.savefig ( outpath, dpi = 300 )
    plt.close ()


def mapplot ( varname    = None,
              inpath     = None,
              outpath    = None,
              unit       = None,
              title      = None,
              midpoint   = None,
              vmin       = None,
              vmax       = None,
              colormap   = "BrBG",
              tstep_no = 0 ):
    
    midpoint = float ( midpoint ) if midpoint else np.nan
    vmin     = float ( vmin )     if vmin     else np.nan
    vmax     = float ( vmax )     if vmax     else np.nan

    nc_fid  = Dataset ( inpath, 'r')

    var       = nc_fid.variables [ varname ]
    units     = var.units
    grid_type = var.grid_type
    time      = nc_fid.variables [ 'time' ][ : ]

    if grid_type == "gaussian":
        lats          = nc_fid.variables [ 'lat' ][ : ]
        lons_orig     = nc_fid.variables [ 'lon' ][ : ] 
        vardata_orig  = var [ tstep_no, :, : ]
        vardata, lons = addcyclic ( vardata_orig, lons_orig )        
        vardata, lons = shiftgrid ( 180., vardata, lons, start = False )
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
    
    vmax    = vardata.max () if np.isnan ( vmax ) else vmax
    vmin    = vardata.min () if np.isnan ( vmin ) else vmin
    vmean   = vardata.mean () if np.isnan ( midpoint ) else midpoint

    print vmax,vmin,vmean
    
    # #print "####### COLORS: VALUES, LEVEL: #########"
    nlev        = 20
    values, level, ticks = colormap_prep ( vmin, vmax, vmean, nlev )
    colors      = plt.cm.get_cmap ( colormap )( values )
    cmap, norm  = mcolors.from_levels_and_colors ( level, colors )
    norm.clip   = True
    norm.vmin   = vmin
    norm.vmax   = vmax
    
    cmap.set_bad   ( color = 'silver', alpha = 1.0 )
    cmap.set_over  ( color = 'r', alpha = 1.0 )
    cmap.set_under ( color = 'g', alpha = 1.0 )
    
    f = plt.figure ( figsize = ( 16, 8 ))  
    ax = f.add_subplot ( 111, projection = ccrs.PlateCarree ())
    ax.set_xlim ([ -180, 180 ])
    ax.set_ylim ([ -90, 90 ])
    ax.coastlines ( linewidth=0.5 )
    ax.add_feature(cart.feature.OCEAN, zorder=100, edgecolor='k', linewidth=0.5, facecolor="silver" )
    
    if grid_type == 'gaussian':
        # correct grid with half distance to north and west
        lons_shifted = lons - np.concatenate (( np.zeros ( 1 ), np.abs ( lons[0:-1] - lons [1:] ) / 2 ), axis = 0 )
        lats_shifted = lats + np.concatenate (( np.abs ( lats [0:-1] - lats [1:] ) / 2, np.zeros ( 1 )), axis = 0 )
        ax.pcolormesh ( lons_shifted, lats_shifted,
                        vardata, cmap=cmap, norm=norm,
                        edgecolors='face',  
                        rasterized=True, zorder=0 )
    elif grid_type == 'unstructured': 
        newxvc, newyvc, triangles, latestvardata = triangles_prep (
            clon_vertices, clat_vertices, vardata )
        plt.tripcolor( newxvc, newyvc, triangles,
                       facecolors=latestvardata, cmap=cmap, norm=norm,
                       edgecolor="face")       
   
    cax = f.add_axes ([ 0.125, 0.04, 0.775, 0.03 ]) # left, bottom, width, height           
    cb = mpl.colorbar.ColorbarBase ( cax, cmap = cmap, norm = norm,
                                     orientation = 'horizontal', extend = 'neither', 
                                     ticks = ticks, format = '%.2e' )
    for t in cb.ax.get_xticklabels ():
        t.set_fontsize ( 14 )    
    plt.text ( 0.5, 1.04, title,
               horizontalalignment = 'center',
               fontsize = 20,
               transform = ax.transAxes )
    cb.ax.text ( 0.89, 1.5, "[%s]" % unit, size = 16 )
    plt.savefig ( outpath, dpi = 300 )
    plt.close ()



def seas_timemean_plot ( varname       = None,
                         inpathlist_gm = None,
                         inpathlist_tm = None,
                         outpath       = None, 
                         unit          = None,
                         title         = None,
                         legendlist_gm = None,
                         legendlist_tm = None,
                         mono          = False,
                         **kwargs ):
    
    listlen = len ( inpathlist_gm )

    assert listlen > 0

    hmargin_unitfactor = 1./2.
    vmargin_factor = 1./15.
        
    fig     = plt.figure ( figsize = ( 10, 7 ))
    
    ymin    = sys.float_info.max
    ymax    = - sys.float_info.max

    for idx in range ( listlen ):

        nc_fid = Dataset ( inpathlist_gm [ idx ], 'r' )  
        time = nc_fid.variables [ 'time' ][ : ]
        gm = nc_fid.variables [ varname ][ : ] 
        nc_fid.close ()

        nc_fid = Dataset ( inpathlist_tm [ idx ], 'r')
        tm = nc_fid.variables [ varname ][ : ] 
        nc_fid.close ()

        xcoords, ycoords = get_seasonal_period_2d_coords ( time, gm [ :, 0, 0 ],
                                                           hmargin_unitfactor, vmargin_factor )

        ymin = np.amin ( gm [ :, 0, 0 ]) if np.amin ( gm [ :, 0, 0 ]) < ymin  else  ymin 
        ymax = np.amax ( gm [ :, 0, 0 ]) if np.amax ( gm [ :, 0, 0 ]) > ymax  else  ymax
        plt.plot ( xcoords, ycoords, colorlist [ idx ] + '-' )
        plt.plot ( time, gm [ :, 0, 0 ],
                   colorlist [ idx ] + 'o',
                   label = legendlist_gm [ idx ] if mono == False else None )
        plt.axhline ( tm [ 0, 0, 0 ],
                      color =  colorlist [ idx ],
                      label = legendlist_tm [ idx ])
        
    ymargin = ( ymax - ymin ) * vmargin_factor
    plt.axis ([ xcoords [ 0 ], xcoords [ 5 ], ymin - ymargin, ymax + ymargin ])  
    plt.title ( title, y = 1.04, fontsize = 18 )
    plt.ylabel ( "[ %s ]" % unit,
                 rotation = 90,
                 position = ( -0.5, .88 ),
                 fontsize = 14 ) 
    plt.xlabel ( "season")
    plt.xticks ( time, seasonlabel )
    plt.grid ( 'on' )
    lgd = plt.legend ( numpoints = 1,  
                       loc = 2,                  
                       bbox_to_anchor = ( 1.05, 1 ))
    plt.savefig ( outpath,
                  dpi = 100, 
                  bbox_extra_artists = ( lgd, ), 
                  bbox_inches = 'tight' )
    plt.close ()


    
def seas_timemean_subplot ( varname = None,
                            inpathlist_gm = None,
                            inpathlist_tm = None,
                            outpath = None,
                            unit = None,
                            title = None,
                            legendlist_gm = None,
                            legendlist_tm = None,
                            ax = None,
                            art = None,
                            mono = False ):

    time    = [ 1,3,5,7 ]

    listlen = len ( inpathlist_gm )

    hmargin_unitfactor = 1./2.
    vmargin_factor = 1./15.
    
    ymin    = sys.float_info.max
    ymax    = - sys.float_info.max


    for idx in range ( listlen ):

        nc_fid = Dataset ( inpathlist_gm [ idx ], 'r' )  
        gm = nc_fid.variables [ varname ][ : ] 
        nc_fid.close ()

        nc_fid = Dataset ( inpathlist_tm [ idx ], 'r')
        tm = nc_fid.variables [ varname ][ : ] 
        nc_fid.close ()

        xcoords, ycoords = get_seasonal_period_2d_coords ( time, gm [ :, 0, 0 ],
                                                           hmargin_unitfactor, vmargin_factor )
        ymin = np.amin ( gm [ :, 0, 0 ]) if np.amin ( gm [ :, 0, 0 ]) < ymin  else  ymin 
        ymax = np.amax ( gm [ :, 0, 0 ]) if np.amax ( gm [ :, 0, 0 ]) > ymax  else  ymax
        ax.plot ( xcoords , ycoords, colorlist [ idx ] + '-' )
        ax.plot ( time, gm [ :, 0, 0 ],
                  colorlist [ idx ] + 'o',
                  label = legendlist_gm [ idx ] if mono == False else None )
        ax.axhline ( tm [ 0, 0, 0 ],
                     color =  colorlist [ idx ],
                     label = legendlist_tm [ idx ] if mono == False else None )
    
    ymargin = ( ymax - ymin ) * vmargin_factor
    ax.axis ([ xcoords [ 0 ], xcoords [ 5 ], ymin - ymargin, ymax + ymargin ])  
    ax.set_title ( title, y = 1.04, fontsize = 18 )
    ax.set_ylabel ( "[ %s ]" % unit,
                    rotation = 90,
                    position = ( -0.5, .88 ),
                    fontsize = 14 ) 
    ax.set_xlabel ( "season")
    ax.set_xticks ( time )
    ax.set_xticklabels ( seasonlabel )
    ax.grid ( 'on' )
    lgd = ax.legend ( numpoints = 1,  
                     loc = 2, 
                     bbox_to_anchor = ( 1.05, 1 ) if mono == False else None,
                     bbox_transform=ax.transAxes if mono == False else None )
    if isinstance ( art, list ) and lgd: art.append ( lgd )




def seas_timemean_plot_wrapper ( varname = None,
                                 varname_can = None,
                                 inpathlist_gm_can = None,
                                 inpathlist_tm_can = None,
                                 subtitle_can = None,
                                 inpathlist_gm = None,
                                 inpathlist_tm = None,
                                 outpath = None, 
                                 unit = None, 
                                 unit_can = None,
                                 subtitle_delta = None,
                                 legendlist_gm = None,
                                 legendlist_tm = None ):

    f = plt.figure ( figsize = ( 10, 15 ))
    show_can = True
    supgrid  = gridspec.GridSpec( 2 if show_can == True else 1, 1 ) 
    art      = []

    ####### first subplot
    subgrid = gridspec.GridSpecFromSubplotSpec ( 1, 1, subplot_spec = supgrid [ 0 ])
    ax      = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax )

    # first subplot call
    seas_timemean_subplot ( varname       = varname_can,
                            inpathlist_gm = inpathlist_gm_can, 
                            inpathlist_tm = inpathlist_tm_can, 
                            unit          = unit_can, 
                            title         = subtitle_can,
                            ax            = ax,
                            art           = art,
                            mono          = True )
    
    ####### second subplot
    subgrid = gridspec.GridSpecFromSubplotSpec (1, 1, subplot_spec = supgrid [ 1 ])
    ax      = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax ) 

    # second subplot call
    seas_timemean_subplot ( varname       = varname,
                            inpathlist_gm = inpathlist_gm, 
                            inpathlist_tm = inpathlist_tm,
                            unit          = unit,
                            title         = subtitle_delta,
                            legendlist_gm = legendlist_gm,
                            legendlist_tm = legendlist_tm,
                            ax            = ax,
                            art           = art )

    plt.savefig ( outpath,
                  bbox_extra_artists = ( art [ 0 ], ) if art else None,
                  bbox_inches = 'tight' if art else None  )

    plt.close ()



def zon_plot ( varname       = None,
               inpathlist    = None,
               outpath       = None, 
               unit          = None,
               title         = None,
               legendlist    = None,
               ax            = None,
               art           = None,
               mono          = False,
               **kwargs ):

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


def zon_plot_wrapper ( varname_can    = None,
                       inpathlist_can = None,
                       unit_can       = None,
                       subtitle_can   = None,
                       varname        = None,
                       inpathlist     = None,
                       unit           = None, 
                       subtitle       = None,
                       outpath        = None, 
                       legendlist     = None ):
   
    f = plt.figure ( figsize = ( 10, 15 ))
    show_can = True
    supgrid = gridspec.GridSpec( 2 if show_can == True else 1, 1 ) 
    art = []

    ####### first subplot
    subgrid = gridspec.GridSpecFromSubplotSpec ( 1, 1, subplot_spec = supgrid [ 0 ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax )

    zon_plot ( varname       = varname_can,
               inpathlist    = inpathlist_can,
               unit          = unit_can,
               title         = subtitle_can,
               legendlist    = [],
               ax            = ax,
               art           = art,
               mono          = True )
    
    ####### second subplot
    subgrid = gridspec.GridSpecFromSubplotSpec (1, 1, subplot_spec = supgrid [ 1 ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax ) 

    zon_plot ( varname       = varname,
               inpathlist    = inpathlist,
               unit          = unit,
               title         = subtitle,
               legendlist    = legendlist,
               ax            = ax,
               art           = art )

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

        
def hist_plot_wrapper ( nbar          = None,
                        colorlist      = None,
                        legendlist     = None,
                        varnamelist    = None,
                        unit           = None,
                        can_name       = None,
                        canno          = None,
                        input_can      = None,
                        subtitle_can   = None,
                        input_delta    = None,
                        subtitle_delta = None,
                        xticklabellist = None,
                        outpath        = None ):
   
    show_can = True

    ndeltas = len ( input_delta )

    f        = plt.figure ( figsize = ( 10, 6 ))
    supgrid  = gridspec.GridSpec( 1, ndeltas if show_can == True else 1 ) 
    art      = []

    ####### first subplot
    subgrid = gridspec.GridSpecFromSubplotSpec ( 1, 1, subplot_spec = supgrid [ 0 ])
    ax = plt.Subplot ( f, subgrid [ :, : ])
    f.add_subplot ( ax )

    accum = vardat =  np.array ([[]])
    nc_fid = Dataset ( input_can [ 0 ])
    for i in range ( nbar ):
        vardat = np.concatenate (( vardat,
                                   nc_fid.variables [ varnamelist [ i ]] [ : ][ 0 ]),
                                 axis = 1 )
    nc_fid.close ()
    accum = np.concatenate (( accum, vardat ), axis = 0 ) if accum != [[]] else vardat
    
    hist_plot (  nbar, colorlist, legendlist,
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
        vardat = np.array ([[]])
        nc_fid = Dataset ( input_delta [ k ])
        for i in range ( nbar ):
            vardat = np.concatenate (( vardat,
                                       nc_fid.variables [ varnamelist [ i ]] [ : ][ 0 ]),
                                     axis = 1 )
        nc_fid.close ()
        accum = np.concatenate (( accum, vardat ), axis = 0 ) if accum != [[]] else vardat

    hist_plot (  nbar, colorlist, legendlist,
                 accum, unit, xticklabellist, subtitle_delta,
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
