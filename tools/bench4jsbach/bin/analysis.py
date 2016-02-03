#! /usr/bin/env python
# -*- encoding:  utf-8 -*-


import os
import fnmatch
from cdo import Cdo
from sets import Set

from pycmbs.data import Data
import matplotlib.pyplot as plt
from pycmbs.mapping import map_plot
from util import logger as dkl_logger
import warnings



# methods #
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

    template = """%(path)s/%(varname)s_%(classification)s_%(expname)s_%(startdate)s-%(finaldate)s.%(suffix)s"""


    def __init__ ( self, **kwargs ):
        try:
            if self.logger: pass
        except:
            self.logger = dkl_logger


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
        
        for exp in self.experiments:
            exp.B4J_varpaths.clear ()
        in_files = []
        

        # 1.  find variable data files
       
        # 1.a) find variable data files
        for exp in self.experiments:
            pattern = '*%s*%s*.%s' % ( A_varname, exp.B4J_name, suffix )
            exp.B4J_varpaths [ "A" ] = A = find_files ( pattern, exp.B4J_path, single = True, raise_e = False )
            if A: in_files.append ( A )
            if exp.B4J_isref or not A:
                pattern = '*%s*%s*.%s' % ( D_varname, exp.B4J_name, suffix )
                exp.B4J_varpaths [ "D" ] = D = find_files ( pattern, exp.B4J_path, single = True, raise_e = True )
                if D: in_files.append ( D )
            if not A:
                pattern = '*%s*%s*.%s' % ( U_varname, exp.B4J_name, suffix )
                exp.B4J_varpaths [ "U" ]  = U = find_files ( pattern, exp.B4J_path, single = True, raise_e = True )
                if U: in_files.append ( U )
            self.logger.debug ( exp )


        yearlist = get_largest_common_yearlist ( in_files )
        cdo = Cdo ()


        # 1.b) calculate albedo if not found
        for exp in self.experiments:
            if not exp.B4J_varpaths [ "A" ]:
                exp.B4J_varpaths [ "A" ] = self.template % { 
                    'path'           : self.data_path,
                    'varname'        : A_varname,
                    'classification' : 'divUD',
                    'expname'        : exp.B4J_name,
                    'startdate'      : str ( yearlist [ 0 ]),
                    'finaldate'      : str ( yearlist [ -1 ]),
                    'suffix'         : suffix }
                cdo.setname ( A_varname,
                              input="-setunit,1 " +\
                              "-setvrange,0.0001,0.99 " +\
                              "-div " +\
                              ( " -selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
                                " " +  exp.B4J_varpaths [ "U" ] ) +\
                              ( " -selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
                                " " +  exp.B4J_varpaths [ "D" ] ),
                              output=exp.B4J_varpaths [ "A" ],
                              options='-f nc -m 0.00011' )
                

        # 2.  preprocessing, calculate seasonal means for largest common period
        for exp in self.experiments:
            exp.B4J_varpaths [ "A_sm" ] = self.template % { 
                'path'           : self.data_path,
                'varname'        : A_varname,
                'classification' : 'yseasmean',
                'expname'        : exp.B4J_name,
                'startdate'      : str ( yearlist [ 0 ]),
                'finaldate'      : str ( yearlist [ -1 ]),
                'suffix'         : suffix }
            cdo.yseasmean ( input="-selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
                            " " + exp.B4J_varpaths [ "A" ],
                            output=exp.B4J_varpaths [ "A_sm" ])
            if exp.B4J_isref:
                exp.B4J_varpaths [ "D_sm" ] = self.template % { 
                    'path'           : self.data_path,
                    'varname'        : D_varname,
                    'classification' : 'yseasmean',
                    'expname'        : exp.B4J_name,
                    'startdate'      : str ( yearlist [ 0 ]),
                    'finaldate'      : str ( yearlist [ -1 ]),
                    'suffix'         : suffix }
                cdo.yseasmean ( input="-selyear," + ','.join ( repr ( y ) for y in yearlist ) + \
                                " " + exp.B4J_varpaths [ "D" ],
                                output=exp.B4J_varpaths [ "D_sm" ],
                                options='-f %s' % suffix )


        # 3.  matrix processing, for all Experiments e and all References r
        for can in self.experiments:
            if not can.B4J_iscan: continue
            for ref in self.experiments:
                if not ref.B4J_isref: continue
                delta_str = "delta-ref%d" % ref.B4J_refno
                can.B4J_varpaths [ delta_str ] = self.template % { 
                    'path'           : self.data_path,
                    'varname'        : A_varname,
                    'classification' : delta_str,
                    'expname'        : exp.B4J_name,
                    'startdate'      : str ( yearlist [ 0 ]),
                    'finaldate'      : str ( yearlist [ -1 ]),
                    'suffix'         : suffix }
                cdo.mul ( input=" -sub " +\
                          can.B4J_varpaths [ "A_sm" ] + " " +\
                          ref.B4J_varpaths [ "A_sm" ] + " " +\
                          ref.B4J_varpaths [ "D_sm" ],
                          output=can.B4J_varpaths [ delta_str ],
                          options='-f %s' % suffix )


        
        self.logger.debug ( "Exiting." )

    
