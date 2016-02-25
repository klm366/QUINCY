#! /usr/bin/env python
# -*- coding:  utf-8 -*-
from __future__ import division


import os, sys
import fnmatch
import warnings
from datetime import datetime
from cdo import Cdo
from sets import Set
import copy
from netCDF4 import Dataset

from pycmbs.data import Data
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import numpy as np
from pycmbs.mapping import map_plot

from util import logger as dkl_logger

from report import ( analyses_var_item,
                     analyses_var_test_item,
                     analyses_var_test_can_item,
                     init_analyses_var_test_can_item )


# methods #

def find_files_from_codes ( codes, path, raise_e = False ):
    print codes
    print path
    



def find_files ( pattern, path, single = True, raise_e = False ):
    assert os.path.isdir ( path ) and isinstance ( pattern, str )
    
    paths = []
    for root, dirnames, filenames in os.walk ( path ):
        files   = fnmatch.filter ( filenames, pattern )
        for f in files:
            paths.append ( os.path.join ( root, f ))

    n_files = len ( paths )
    if single:
        if not n_files:
            if raise_e:
                raise Exception ( "no file" )
            else:
                ret = None
        elif n_files > 1:
            if raise_e:
                raise Exception ( "multiple files" )
            else:
                warnings.warn ( "multiple files, take first" )
                ret = paths [ 0 ]
        else:
            ret = paths [ 0 ]
    else:
        if not n_files and raise_e:
            raise Exception ( "no files" )
        ret = paths

    return ret

                    
    # albedo
def get_largest_common_yearlist ( filelist ):
    """Get largest common period for list of datafiles <filelist> given as argument."""
    assert filelist, "No filelist given."
    assert all ( map ( lambda c: os.path.isfile ( c ), filelist )),\
           "Not all files found." 

    cdo = Cdo ()

    yearlists = [ map ( int, ''.join ( cdo.showyear (
        input = f )). encode (). split ())
                  for f in filelist ]
    startyear = max ( l [ 0 ] for l in yearlists )
    endyear = min ( l [ len ( l ) - 1 ] for l in yearlists )
    return range ( startyear, endyear + 1 ) if startyear <= endyear else []


class Analysis ( object ):

    path_template = """%(path)s/%(varname)s_%(classification)s_%(expname)s_%(startdate)s-%(finaldate)s.%(suffix)s"""


    def __init__ ( self, **kwargs ):
        try:
            if self.logger: pass
        except:
            self.logger = dkl_logger


    def  get_seasonal_period_2d_coords ( self, xvals, yvals, xmargin_unitfactor, ymargin_factor ):
        
        xunit         = xvals [ 1 ] - xvals [ 0 ]
        xmargin       = xunit * xmargin_unitfactor
        left_x_coord  = xvals [ 0 ] - xmargin
        right_x_coord = xvals [ 3 ] + xmargin
        
        y_dist        = yvals [ 0 ] - yvals [ 3 ]
        left_y_addon  = y_dist * ( 1 - xmargin_unitfactor ) 
        right_y_addon = y_dist * xmargin_unitfactor
        left_y_coord  = yvals [ 3 ] + left_y_addon
        right_y_coord = yvals [ 3 ] + right_y_addon
        
        x_coords = np.insert ( np.insert ( xvals, 0, left_x_coord ), 5, right_x_coord )
        y_coords = np.insert ( np.insert ( yvals, 0, left_y_coord ), 5, right_y_coord )
        
        return x_coords, y_coords

    def albedo_analyse ( self ):
        """
NAME 
        albedo_analyse

DESCRIPTION
        This method guides the analyses and the benchmarking for the variable albedo. It plots the results and generates texts for 
        the report. 

DETAILS
        for all experiments in benchmark:
    1.  find variable data files:
        - try to find output files for variable albedo                          :  A
        - if not found, calculate albedo from variables rsus and rsds           :  A = U/D
        - if none of the needed output files found                              :  EXIT ON ERROR

    2.  preprocessing, calculate seasonal means for largest common period: 
        - albedo seasonal means                                                 :  A_sm = sm ( A )
        - only for references:  rsds seasonal means                             :  D_sm = sm ( D )

    3.  matrix processing, for all Experiments e and all References r:
        - difference of albedo seasonal means weighted by rsds seasonal means   :  Δ_er  = ( A_sm_e - A_sm_r ) * D_sm_r 
        - global mean                                                           :  gm_er = ∑_wx_x ( Δ_er * land-sea-mask )
        - time mean                                                             :  tm_er = ∑_wt_x ( gm_er )
        - global mean deviation                                                 :  gd_er = ∑_wx_x ( abs ( Δ_er * land-sea-mask )
        - time mean deviation                                                   :  td_er = ∑_wt_x ( gd_er )

    4.  prepare output:
        - store plots for results of calculations in 3.
        - if report is switched on, prepare text blocks
"""

        try:
            from bench import Bench_Object
            if not isinstance ( self, Bench_Object ):
                raise
        except:
            self.logger.warning ( "This method has to be used by a Bench_Object." )

        self.logger.debug ( "Starting." )

        A_varname = "albedo"
        D_varname = "rsds"
        U_varname = "rsus"
        suffix    = "nc"
        
        self.analyses [ A_varname ] = copy.deepcopy ( analyses_var_item )

        # preparation of runtime storage and texts for analyses 
        # 3.a)
        teXt =  """
The difference between candidate and reference experiments albedo seasonal means 
is weighted with the reference  experiments downwelling radiation seasonal means:
\\begin{eqnarray*}
\\Delta_{CR} & = & ( A_{sm\_C} - A_{sm\_R}) * D_{sm\_R}, \\\\
        && A_{sm} = \\texttt{yseasmean}( A ), \\\\
        && D_{sm} = \\texttt{yseasmean} ( D ),
\\end{eqnarray*}
with
\\begin{eqnarray*}
\\Delta & := & \\textnormal{ weighted difference of albedo in } \\frac{W}{m^{2}} \\\\
A & := & \\textnormal{albedo, quotient of upwelling, downwelling radiation in } \% \\\\
D & := & \\textnormal{downwelling solar radiation in } \\frac{W}{m^{2}}
\\end{eqnarray*}

and \\texttt{yseasmean}($X$) defined as
\\begin{eqnarray*}
  y_{SX} & = & n_{S}^{-1}\\sum_{i=0}^{n_{S}} x_{iX}, \\\\
        && t\\mapsto i, t \in S,  \\\\
   && S \in \{\\textnormal{DJF, MAM, JJA, SON} \}.
\\end{eqnarray*}

"""  
        delta = copy.deepcopy ( analyses_var_test_item )
        delta [ 'name' ]   = 'delta'
        delta [ 'teXts' ]. append ( teXt )
        delta [ 'title' ]  = 'Seasonal mean difference'
        delta [ 'matrix' ] = True
        delta [ 'vector' ] = False
        delta [ 'cans' ]   = [ init_analyses_var_test_can_item ( name = can.B4J_name ) \
                               for can in self.experiments if can.B4J_iscan ]
        self.analyses [ A_varname ].append ( delta )

        
        # 3.b,c)
        teXt = """
We calculate the global and time means of the difference of albedo seasonal means 
weighted by downwelling radiation seasonal means from the last section. The global 
means is 
\\begin{eqnarray*}
gm_{CR} & = & \\texttt{fldmean} ( \\Delta_{CR} * \\texttt{land-sea-mask} )
\\end{eqnarray*}
with
\\begin{eqnarray*}
gm & := & \\textnormal{ global mean of weighted difference of albedo  in } \\frac{W}{m^{2}}
\\end{eqnarray*}
and \\texttt{fldmean}($X$) defined as
\\begin{eqnarray*}
  y_{t} & = & \\sum_{i=0}^{n_{X}} w_{Xi} x_{ti} \\\\
\\end{eqnarray*}
The time mean is then calculated with
\\begin{eqnarray*}
tm_{CR} & = & \\texttt{timemean} ( gm_{CR} )
\\end{eqnarray*}
\\begin{eqnarray*}
tm & := & \\textnormal{ time and global mean of weighted difference of seasonal mean of albedo  in } \\frac{W}{m^{2}}
\\end{eqnarray*}
and \\texttt{timemean}($X$) defined as
\\begin{eqnarray*}
  y & = &  n^{-1} \\sum_{i=0}^{n_{T}} x_{i} \\\\
\\end{eqnarray*}

"""
        globalNtime_mean = copy.deepcopy ( analyses_var_test_item )
        globalNtime_mean [ 'name' ]   = 'globalNtime_mean'
        globalNtime_mean [ 'teXts' ]. append ( teXt )
        globalNtime_mean [ 'title' ]  = 'Global and time mean'
        globalNtime_mean [ 'matrix' ] = False
        globalNtime_mean [ 'vector' ] = True
        globalNtime_mean [ 'cans' ]   = [ init_analyses_var_test_can_item ( name = can.B4J_name ) \
                                          for can in self.experiments if can.B4J_iscan ]
        self.analyses [ A_varname ].append ( globalNtime_mean )
        
        
        # 3.d,e)
        teXt = """
We calculate the global and time mean deviation of the difference of albedo seasonal mean 
weighted by downwelling radiation seasonal mean from the last section. The global 
mean deviation is 
\\begin{eqnarray*}
gd_{CR} & = & \\texttt{fldmean} ( \\sqrt{(\\Delta_{CR} * \\texttt{land-sea-mask})^{2}} )
\\end{eqnarray*}
with
\\begin{eqnarray*}
gd & := & \\textnormal{ global mean deviation of weighted difference of albedo  in } \\frac{W}{m^{2}}
\\end{eqnarray*}
and \\texttt{fldmean}($X$) defined as
\\begin{eqnarray*}
  y_{t} & = & \\sum_{i=0}^{n_{X}} w_{Xi} x_{ti} \\\\
\\end{eqnarray*}
The time mean is then calculated with
\\begin{eqnarray*}
td_{CR} & = & \\texttt{timemean} ( gd_{CR} )
\\end{eqnarray*}
\\begin{eqnarray*}
td & := & \\textnormal{time and global mean deviation of weighted difference} \\\\
 & & \\textnormal{seasonal mean of albedo  in } \\frac{W}{m^{2}}
\\end{eqnarray*}
and \\texttt{timemean}($X$) defined as
\\begin{eqnarray*}
  y & = &  n^{-1} \\sum_{i=0}^{n_{T}} x_{i} \\\\
\\end{eqnarray*}

"""
        globalNtime_dev = copy.deepcopy ( analyses_var_test_item )
        globalNtime_dev [ 'name' ]   = 'globalNtime_dev'
        globalNtime_dev [ 'teXts' ]. append ( teXt )
        globalNtime_dev [ 'title' ]  = 'Global and time mean deviation'
        globalNtime_dev [ 'matrix' ] = False
        globalNtime_dev [ 'vector' ] = True
        globalNtime_dev [ 'cans' ]   = [ init_analyses_var_test_can_item ( name = can.B4J_name ) \
                                          for can in self.experiments if can.B4J_iscan ]
        self.analyses [ A_varname ].append ( globalNtime_dev )
        
        for exp in self.experiments:
            exp.B4J_index.clear ()
        in_files = []


        for exp in self.experiments:
            find_files_from_codes ( [21,22], exp.B4J_path )
        

        # 1.  find variable data files
       
        # # 1.a) find variable data files
        # for exp in self.experiments:
        #     pattern = '*%s*%s*.%s' % ( A_varname, exp.B4J_name, suffix )
        #     exp.B4J_index [ "A" ] = A = find_files ( pattern, exp.B4J_path, single = True, raise_e = False )
        #     if A: in_files.append ( A )
        #     if exp.B4J_isref or not A:
        #         pattern = '*%s*%s*.%s' % ( D_varname, exp.B4J_name, suffix )
        #         exp.B4J_index [ "D" ] = D = find_files ( pattern, exp.B4J_path, single = True, raise_e = True )
        #         if D: in_files.append ( D )
        #     if not A:
        #         pattern = '*%s*%s*.%s' % ( U_varname, exp.B4J_name, suffix )
        #         exp.B4J_index [ "U" ]  = U = find_files ( pattern, exp.B4J_path, single = True, raise_e = True )
        #         if U: in_files.append ( U )
        #     self.logger.debug ( exp )


        # yearlist = get_largest_common_yearlist ( in_files )
        # cdo = Cdo ()
        
        # # land sea mask, todo: tests
        # slf_varname = "slf"
        # land_sea_mask = self.data_path + os.sep + \
        #                 slf_varname + "_" + os.path.basename ( self.land_sea_mask )
        # cdo.selname ( slf_varname, input = self.land_sea_mask, output = land_sea_mask )


        # # 1.b) calculate albedo if not found
        # for exp in self.experiments:
        #     if not exp.B4J_index [ "A" ]:
        #         exp.B4J_index [ "A" ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : A_varname,
        #             'classification' : 'divUD',
        #             'expname'        : exp.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.setname ( A_varname,
        #                       input="-setunit,1 " +\
        #                       "-setvrange,0.0001,0.99 " +\
        #                       "-div " +\
        #                       ( " -selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
        #                         " " +  exp.B4J_index [ "U" ] ) +\
        #                       ( " -selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
        #                         " " +  exp.B4J_index [ "D" ] ),
        #                       output=exp.B4J_index [ "A" ],
        #                       options='-f nc -m 0.00011' )
                

        # # 2.  preprocessing, calculate seasonal mean for largest common period
        # for exp in self.experiments:
        #     exp.B4J_index [ "A_sm" ] = self.path_template % { 
        #         'path'           : self.data_path,
        #         'varname'        : A_varname,
        #         'classification' : 'yseasmean',
        #         'expname'        : exp.B4J_name,
        #         'startdate'      : str ( yearlist [ 0 ]),
        #         'finaldate'      : str ( yearlist [ -1 ]),
        #         'suffix'         : suffix }
        #     cdo.yseasmean ( input="-selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
        #                     " " + exp.B4J_index [ "A" ],
        #                     output=exp.B4J_index [ "A_sm" ])
        #     if exp.B4J_isref:
        #         exp.B4J_index [ "D_sm" ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : D_varname,
        #             'classification' : 'yseasmean',
        #             'expname'        : exp.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.yseasmean ( input="-selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
        #                         " " + exp.B4J_index [ "D" ],
        #                         output=exp.B4J_index [ "D_sm" ],
        #                         options='-f %s' % suffix )


        # # 3.  matrix processing, for all Experiments e and all References r
        # # 4.  prepare output

        # seasons =  ['DJF','MAM','JJA','SON']

        # # 4.a) collecting information for common adjustment of colorbar
        # delta_min = sys.float_info.max
        # delta_max = sys.float_info.min
        # can.B4J_index [ "delta_maxs" ] = []
        # can.B4J_index [ "delta_mins" ] = []
    
        
        # for can in self.experiments:
        #     if not can.B4J_iscan: continue
                                   
        #     for ref in self.experiments:
        #         if not ref.B4J_isref:            continue
        #         if ref.B4J_name == can.B4J_name: continue
                
        #         # 3.a) difference of albedo seasonal means weighted by rsds seasonal means
        #         delta_str = "delta-ref%d" % ref.B4J_refno
        #         can.B4J_index [ delta_str ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : A_varname,
        #             'classification' : delta_str,
        #             'expname'        : can.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.mul ( input=" -sub " +\
        #                   can.B4J_index [ "A_sm" ] + " " +\
        #                   ref.B4J_index [ "A_sm" ] + " " +\
        #                   ref.B4J_index [ "D_sm" ],
        #                   output=can.B4J_index [ delta_str ],
        #                   options='-f %s' % suffix )

        #         nc_fid = Dataset ( can.B4J_index [ delta_str ], 'r')
        #         albedo_delta = nc_fid.variables [ A_varname ][ : ] 
        #         nc_fid.close ()
        #         vmax = albedo_delta.max ()
        #         vmin = albedo_delta.min ()
        #         can.B4J_index [ "delta_maxs" ].append ( vmax ) 
        #         can.B4J_index [ "delta_mins" ].append ( vmin )
        #         if vmax > delta_max:
        #             delta_max = vmax
        #         if vmin < delta_min:
        #             delta_min = vmin

                    
        #         # 3.b) global mean
                
        #         gm_str = "gm-ref%d" % ref.B4J_refno
        #         can.B4J_index [ gm_str ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : A_varname,
        #             'classification' : gm_str,
        #             'expname'        : can.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.fldmean ( input=" -mul " +\
        #                       can.B4J_index [ delta_str ] + " " +\
        #                       land_sea_mask,
        #                       output=can.B4J_index [ gm_str ],
        #                       options='-f %s' % suffix )

        #         # 3.c) time mean
        #         tm_str = "tm-ref%d" % ref.B4J_refno
        #         can.B4J_index [ tm_str ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : A_varname,
        #             'classification' : tm_str,
        #             'expname'        : can.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.timmean ( input=can.B4J_index [ gm_str ],
        #                       output=can.B4J_index [ tm_str ],
        #                       options='-f %s' % suffix )

        #         # 3.d) global mean deviation
        #         gd_str = "gd-ref%d" % ref.B4J_refno
        #         can.B4J_index [ gd_str ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : A_varname,
        #             'classification' : gd_str,
        #             'expname'        : can.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.fldmean ( input="-abs -mul " +\
        #                       can.B4J_index [ delta_str ] + " " +\
        #                       land_sea_mask,
        #                       output=can.B4J_index [ gd_str ],
        #                       options='-f %s' % suffix )

        #         # 3.e) time mean deviation
        #         td_str = "td-ref%d" % ref.B4J_refno
        #         can.B4J_index [ td_str ] = self.path_template % { 
        #             'path'           : self.data_path,
        #             'varname'        : A_varname,
        #             'classification' : td_str,
        #             'expname'        : can.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]),
        #             'suffix'         : suffix }
        #         cdo.timmean ( input=can.B4J_index [ gd_str ],
        #                       output=can.B4J_index [ td_str ],
        #                       options='-f %s' % suffix )

                
                
        # # 4.) plotting section
        
        # # 4.a.3.a) colormap adjustment for delta (as serie of mapplots )

        # # arguments
        # nlevs                = 400
        # tick_interleave_no   = int ( nlevs // 16 )
        # midval               =  0
        # vmin                 = delta_min
        # vmax                 = delta_max

        # # calculate colormap, colorbar, ticks
        # assert nlevs        >=  5       
        # unit                 = ( vmax - vmin )  / ( nlevs - 1 )
        # baselevel_left_val   = midval - unit/2
        # baselevel_right_val  = midval + unit/2
        # assert vmin < baselevel_left_val and \
        #         baselevel_left_val < midval and \
        #         midval < baselevel_right_val and  \
        #         baselevel_right_val < vmax
        # nlevs_left           = int (( baselevel_left_val - vmin ) // unit ) + \
        #                          ( 1 if int (( baselevel_left_val - vmin ) % unit ) == 0 else 0 )
        
        # nlevs_right          = int (( vmax - baselevel_right_val ) // unit ) + \
        #                          ( 1 if int (( vmax - baselevel_right_val ) % unit ) == 0 else 0 )
        # nlevs_revised        = nlevs_left + nlevs_right + 1
        # leftlevel_left_val   = baselevel_left_val -  nlevs_left * unit
        # rightlevel_right_val = baselevel_right_val + nlevs_right * unit
        # levels               = np.arange ( leftlevel_left_val, rightlevel_right_val + unit/2, unit )
        # ticks_left           = np.arange ( midval, leftlevel_left_val, - unit ) [::-1]
        # ticks_right          = np.delete ( np.arange ( midval, rightlevel_right_val, unit ), 0 )
        # ticks                = np.concatenate (( ticks_left, ticks_right ), axis = 0 )
        # ticks_interleave     = np.take ( ticks, np.arange ( nlevs_left % tick_interleave_no,
        #                                                     len ( ticks ) - 1, tick_interleave_no ))
        # vals                 = np.interp ( ticks, [ leftlevel_left_val, midval, rightlevel_right_val ],
        #                                    [ 0, 0.5, 1 ])
        # colors               = plt.cm.bwr ( vals )
        # cmap, norm           = mcolors.from_levels_and_colors ( levels, colors )


        # for can in self.experiments:
        #     if not can.B4J_iscan: continue
        #     for ref in self.experiments:
        #         if not ref.B4J_isref: continue
        #         if ref.B4J_name == can.B4J_name: continue
        
                
        #         # 4.a) store plots for results of calculations in 3.

        #         # 4.a.3.a)
                
        #         delta = next (( item for item in self.analyses [ A_varname ] if item [ 'name' ] == "delta" ), None )
        #         assert delta
                
        #         dates = cdo.showtimestamp ( input = can.B4J_index [ delta_str ] )[ 0 ].split ( " ")
        #         dates = [ datetime.strptime ( str ( x ), '%Y-%m-%dT%H:%M:%S') for x in dates if x ]
        #         plot_delta_str = delta_str + "plot"
        #         can.B4J_index [ plot_delta_str ] =  self.path_template % {
        #             'path'           : self.graphic_path,
        #             'varname'        : A_varname,
        #             'classification' : delta_str,
        #             'expname'        : can.B4J_name,
        #             'startdate'      : str ( yearlist [ 0 ]),
        #             'finaldate'      : str ( yearlist [ -1 ]), 
        #             'suffix'         : 'png' }


        #         f = plt.figure ( figsize = ( 16, 10 ))
        #         ax = [ f.add_subplot ( 2, 2, i + 1 ) for i in range ( 4 )]
        #         d = [ None ] * 4                
        #         gs1 = gridspec.GridSpec ( 2, 2 )
        #         gs1.update ( wspace = -0.08, hspace = -0.04  )               
        #         for i in range ( 4 ):
        #             ax [i] = plt.subplot ( gs1 [ i ])
        #             d [i] = Data ( can.B4J_index [ delta_str ],
        #                            'albedo',
        #                            read       = True,
        #                            start_time = dates [ i ],
        #                            label      = seasons [ i ] )    
        #             x = map_plot ( d[i],
        #                            ax                  = ax[i],
        #                            vmin                = leftlevel_left_val,
        #                            vmax                = rightlevel_right_val,
        #                            use_basemap         = True,
        #                            show_colorbar       = False,
        #                            return_plot_handler = True,
        #                            cmap_data           = cmap )
                    
        #         cax = f.add_axes ([ 0.17, 0.07, 0.7, 0.03])
        #         cb = f.colorbar ( x.im, cax=cax, orientation="horizontal", ticks=ticks_interleave, format='%.3f' )
        #         cb.ax.text ( 0.97, 2.5, r'$\frac{W}{m^{2}}$', size = 32 )
        #         plt.suptitle ( r'$\Delta_{C_{%d}R_{%d}}$' % (
        #             can.B4J_canno, ref.B4J_refno ) + \
        #                        r',   $C_{%d}$: %s, $R_{%d}$: %s' % (
        #                            can.B4J_canno, can.B4J_name,
        #                            ref.B4J_refno, ref.B4J_name ),  fontsize = 24 )
        #         plt.savefig ( can.B4J_index [ plot_delta_str ] )
        #         plt.close ()

        #         delta_can = next (( item for item in delta [ 'cans' ] if item [ 'name' ] == can.B4J_name ), None )
        #         assert delta_can
        #         delta_can [ 'infiles' ].append ( can.B4J_index [ delta_str ])
        #         delta_can [ 'plots' ].append ( can.B4J_index [ plot_delta_str ])
                

        # # 4.a.3.b,c,d,e
        
        # colors = [ 'r', 'g', 'b', 'c', 'y', 'm' ]
        # colored_dots = [ 'ro', 'go', 'bo', 'co', 'yo', 'mo' ]
        # colored_lines = [ 'r-', 'g-', 'b-', 'c-', 'y-', 'm-' ]

        # hmargin_unitfactor = 1/2
        # vmargin_factor = 1/15
        
        
        # # 4.a.3.b,c)

        # globalNtime_mean = next (( item for item in self.analyses [ A_varname ] if item [ 'name' ] == "globalNtime_mean" ), None )
        # assert globalNtime_mean
        
        # for can in self.experiments:
        #     if not can.B4J_iscan: continue
            
        #     plot_globalNtime_mean_str =  self.path_template % {
        #         'path'           : self.graphic_path,
        #         'varname'        : A_varname,
        #         'classification' : "gNtmeans",
        #         'expname'        : can.B4J_name,
        #         'startdate'      : str ( yearlist [ 0 ]),
        #         'finaldate'      : str ( yearlist [ -1 ]),  
        #         'suffix'         : 'png' }
        #     globalNtime_mean_can = next (( item for item in globalNtime_mean [ 'cans' ] \
        #                                    if item [ 'name' ] == can.B4J_name ), None )
        #     assert globalNtime_mean_can
        #     globalNtime_mean_can [ 'plots' ].append ( plot_globalNtime_mean_str )
            
        #     fig = plt.figure()
            
        #     refs = []
        #     xcoords = []
        #     ymin = sys.float_info.max
        #     ymax = sys.float_info.min
            
        #     for ref in self.experiments:
        #         if not ref.B4J_isref: continue
        #         if ref.B4J_name == can.B4J_name: continue
                
        #         gm_str = "gm-ref%d" % ref.B4J_refno
        #         tm_str = "tm-ref%d" % ref.B4J_refno
        #         globalNtime_mean_can [ 'infiles' ].append ( can.B4J_index [ gm_str ])
        #         nc_fid = Dataset ( can.B4J_index [ gm_str ])  
        #         time = nc_fid.variables [ 'time' ][ : ]
        #         albedo_gm = nc_fid.variables['albedo'][:] 
        #         nc_fid.close ()

        #         globalNtime_mean_can [ 'infiles' ].append ( can.B4J_index [ tm_str ])
        #         nc_fid = Dataset ( can.B4J_index [ tm_str ], 'r')
        #         albedo_tm = nc_fid.variables [ 'albedo' ][ : ] 
        #         nc_fid.close ()

        #         xcoords, ycoords = self.get_seasonal_period_2d_coords ( time, albedo_gm [ :, 0, 0 ],
        #                                                 hmargin_unitfactor, vmargin_factor )
        #         ymin = np.amin ( albedo_gm [ :, 0, 0 ]) if np.amin ( albedo_gm [ :, 0, 0 ]) < ymin  else  ymin 
        #         ymax = np.amax ( albedo_gm [ :, 0, 0 ]) if np.amax ( albedo_gm [ :, 0, 0 ]) > ymax  else  ymax 
        #         plt.plot ( xcoords, ycoords, colored_lines [ ref.B4J_refno ])
        #         plt.plot ( time, albedo_gm [ :, 0, 0 ],
        #                    colored_dots [ ref.B4J_refno ],
        #                    label = r'$gm_{C_{%s}R_{%s}}$'% ( can.B4J_refno, ref.B4J_refno ) +  r', $R_{%s}$: ' % ref.B4J_refno + ref.B4J_name )
        #         plt.axhline ( albedo_tm [ 0, 0, 0 ],
        #                       color =  colors [ ref.B4J_refno ],
        #                       label = r'$tm_{C_{%s}R_{%s}}$'% ( can.B4J_refno, ref.B4J_refno ))
             

        #     ymargin = ( ymax - ymin ) * vmargin_factor
        #     plt.axis ([ xcoords [ 0 ], xcoords [ 5 ], ymin - ymargin, ymax + ymargin ])     
        #     plt.title ( r'$gm_{C_{%d}R}$' % can.B4J_canno + \
        #                 " and " +\
        #                 r'$tm_{C_{%d}R}$' % can.B4J_canno,
        #                 y = 1.08 )
        #     plt.ylabel ( r'$\frac{W}{m^{2}}$',
        #                  rotation = 0,
        #                  position = ( -0.4, .96 ),
        #                  fontsize = 20 )
        #     plt.xlabel ( "season")
        #     plt.xticks ( time, seasons )
        #     plt.grid ( 'on' )
        #     lgd = plt.legend ( numpoints = 1,  
        #                        loc = 2,                  
        #                        bbox_to_anchor = ( 1.05, 1 ))
        #     plt.savefig ( plot_globalNtime_mean_str, 
        #                   dpi = 100, 
        #                   bbox_extra_artists = ( lgd, ), 
        #                   bbox_inches = 'tight' )
        #     plt.close ()
            
        
        # # 4.a.3.d,e)

        # globalNtime_dev = next (( item for item in self.analyses [ A_varname ] if item [ 'name' ] == "globalNtime_dev" ), None )
        # assert globalNtime_dev
        
        # for can in self.experiments:
        #     if not can.B4J_iscan: continue
            
        #     plot_globalNtime_dev_str =  self.path_template % {
        #         'path'           : self.graphic_path,
        #         'varname'        : A_varname,
        #         'classification' : "gNtdev",
        #         'expname'        : can.B4J_name,
        #         'startdate'      : str ( yearlist [ 0 ]),
        #         'finaldate'      : str ( yearlist [ -1 ]),  
        #         'suffix'         : 'png' }
        #     globalNtime_dev_can = next (( item for item in globalNtime_dev [ 'cans' ] \
        #                                    if item [ 'name' ] == can.B4J_name ), None )
        #     assert globalNtime_dev_can
        #     globalNtime_dev_can [ 'plots' ].append ( plot_globalNtime_dev_str )
            
        #     fig = plt.figure()
            
        #     refs = []
        #     xcoords = []
        #     ymin = sys.float_info.max
        #     ymax = sys.float_info.min
            
        #     for ref in self.experiments:
        #         if not ref.B4J_isref: continue
        #         if ref.B4J_name == can.B4J_name: continue
                
        #         gd_str = "gd-ref%d" % ref.B4J_refno
        #         td_str = "td-ref%d" % ref.B4J_refno
        #         globalNtime_dev_can [ 'infiles' ].append ( can.B4J_index [ gd_str ])
        #         nc_fid = Dataset ( can.B4J_index [ gd_str ])  
        #         time = nc_fid.variables [ 'time' ][ : ]
        #         albedo_gd = nc_fid.variables['albedo'][:] 
        #         nc_fid.close ()

        #         globalNtime_dev_can [ 'infiles' ].append ( can.B4J_index [ td_str ])
        #         nc_fid = Dataset ( can.B4J_index [ td_str ], 'r')
        #         albedo_td = nc_fid.variables [ 'albedo' ][ : ] 
        #         nc_fid.close ()

        #         xcoords, ycoords = self.get_seasonal_period_2d_coords ( time, albedo_gd [ :, 0, 0 ],
        #                                                 hmargin_unitfactor, vmargin_factor )
        #         ymin = np.amin ( albedo_gd [ :, 0, 0 ]) if np.amin ( albedo_gd [ :, 0, 0 ]) < ymin  else  ymin 
        #         ymax = np.amax ( albedo_gd [ :, 0, 0 ]) if np.amax ( albedo_gd [ :, 0, 0 ]) > ymax  else  ymax 
        #         plt.plot ( xcoords, ycoords, colored_lines [ ref.B4J_refno ])
                
        #         plt.plot ( time, albedo_gd [ :, 0, 0 ],
        #                    colored_dots [ ref.B4J_refno ],
        #                    label = r'$gd_{C_{%s}R_{%s}}$'% ( can.B4J_refno, ref.B4J_refno ) +  r', $R_{%s}$: ' % ref.B4J_refno + ref.B4J_name )
        #         plt.axhline ( albedo_td [ 0, 0, 0 ],
        #                       color =  colors [ ref.B4J_refno ],
        #                       label = r'$td_{C_{%s}R_{%s}}$'% ( can.B4J_refno, ref.B4J_refno ))

    
        #     ymargin = ( ymax - ymin ) * vmargin_factor
        #     plt.axis ([ xcoords [ 0 ], xcoords [ 5 ], ymin - ymargin, ymax + ymargin ])               
        #     plt.title ( r'$gd_{C_{%d}R}$' % can.B4J_canno + \
        #                 " and " +\
        #                 r'$td_{C_{%d}R}$' % can.B4J_canno,
        #                 y = 1.08 )
        #     plt.ylabel ( r'$\frac{W}{m^{2}}$',
        #                  rotation = 0,
        #                  position = ( -0.4, .96 ),
        #                  fontsize = 20 )
        #     plt.xlabel ( "season")
        #     plt.xticks ( time, seasons )
        #     plt.grid ( 'on' )
        #     lgd = plt.legend ( numpoints = 1,  
        #                        loc = 2,                  
        #                        bbox_to_anchor=(1.05,1))
        #     plt.savefig ( plot_globalNtime_dev_str, 
        #                   dpi = 100, 
        #                   bbox_extra_artists= ( lgd, ), 
        #                   bbox_inches='tight' )
        #     plt.close ()

        
    
        # # 4.b) if report is switched on, prepare text blocks
            
        self.logger.debug ( "Exiting." )

    
