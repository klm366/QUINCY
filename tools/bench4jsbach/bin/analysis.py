#! /usr/bin/env python
# -*- coding:  utf-8 -*-
from __future__ import division

from collections import OrderedDict

import os, sys, shutil
import fnmatch, re
import warnings
from datetime import datetime
import subprocess
import itertools
import configobj
from cdo import Cdo
from sets import Set
import copy
from netCDF4 import Dataset
from util import print_item_or_list_rec

from pycmbs.data import Data
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import matplotlib.ticker as ticker
import numpy as np
from pycmbs.mapping import map_plot

from util import logger as dkl_logger, touch
from experiment import MPIESM1_Experiment, CMIP5_Experiment

from report import teXify

streams = { "echam": [ "echam", "co2" ],
            "jsbach": [ "jsbach", "land", "surf", "veg", "yasso", "accw" ]}


# methods #

def  get_seasonal_period_2d_coords ( xvals, yvals, xmargin_unitfactor, ymargin_factor ):

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



def get_vardes_from_codes ( exp_name, code, path, raise_e = False ):

    ret = None
    for fn in os.listdir ( path + "/work" ):
        stream = re.search ( r'%s_.*_(.*).codes' % exp_name, fn )
        if stream:
            f = open ( path + '/work/' + fn, 'r' )
            for line in f:                
                if line.lstrip().startswith ( "%s " % code ):
                    vardes = { "stream": stream.group ( 1 )}
                    vardes [ "path" ]             = path + '/work/' + fn
                    vardes [ "code" ],       line = line.lstrip().split ( " ", 1 )
                    vardes [ "ntiles" ],     line = line.lstrip().split ( " ", 1 )
                    vardes [ "short_name" ], line = line.lstrip().split ( " ", 1 )
                    vardes [ "attr5" ],      line = line.lstrip().split ( " ", 1 )
                    vardes [ "attr6" ],      line = line.lstrip().split ( " ", 1 )
                    vardes [ "long_name" ]        = line.strip()
                    ret = vardes
                    break
            f.close ()
            if ret: break

    return ret



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
    # TODO: take varname out as it is already in path
    path_template = """%(path)s/%(classification)s_%(expname)s_%(startdate)s-%(finaldate)s.%(suffix)s"""


    def __init__ ( self, **kwargs ):
        super ( Analysis, self ).__init__()
        try:
            if self.logger: pass
        except:
            self.logger = dkl_logger
    

    def diagnostic_generic ( self, fieldname ):
        try:
            from bench import Bench_Object
            if not isinstance ( self, Bench_Object ):
                raise
        except:
            self.logger.warning ( "This method has to be used by a Bench_Object." )

        self.logger.debug ( "Starting." )

        
        fielddict = self.diagnostic_fields [ fieldname ]

        ############ FETCH FILES ###############
        in_files = []
        for exp in self.experiments:
            for fileid, filedict in fielddict [ 'fetch' ].iteritems ():
                exp.B4J_index [ fileid ] = self.path_template % { 
                    'path'           : fielddict [ "data_path" ],
                    'classification' : fileid,
                    'expname'        : exp.B4J_name,
                    'startdate'      : fielddict [ 'startyear' ],
                    'finaldate'      : fielddict [ 'finalyear' ],
                    'suffix'         : fielddict [ 'suffix' ]}
                in_files.append ( exp.fetch_outdata ( filedict [ 'name' ], exp.B4J_index [ fileid ]))
        self.logger.debug ( "Fetched input files for field %s and all experiments." %  fieldname )
        
        ########## CHECK TIME PERIOD, WORKAROUND ON WRONG PLACE ##################
        yearlist = get_largest_common_yearlist ( in_files ) # TODO: put this call before actually fetching the data 
        if int ( yearlist [ 0 ]) > int ( fielddict [ 'startyear' ]) or\
           int ( yearlist [ -1 ]) < int ( fielddict [ 'finalyear' ]):
            raise Exception ( "Error in benchmark configuration, choosen period does not match data files." )
            sys.exit ( 1 )
        startyear =  fielddict [ 'startyear' ]
        finalyear = fielddict [ 'finalyear' ]
        yearlist = range ( int ( startyear ), int ( finalyear ) + 1 )
        
        cdo = Cdo ()
        
        ########## FETCH LAND SEA MASK, TODO: TESTS #######################
        slf_varname = "slf"
        land_sea_mask = fielddict [ "data_path" ] + os.sep + \
                        slf_varname + "_" + os.path.basename ( self.land_sea_mask )
        cdo.selname ( slf_varname, input = self.land_sea_mask, output = land_sea_mask )
        self.logger.debug ( "Fetched land sea mask '%s'." %  land_sea_mask )

        ################ PREPROCESSING ####################tre
        executable = self.get_b4j_path (  fielddict [ 'preprocessing' ][ 'exec' ], "scripts" )
        input_list = fielddict [ 'preprocessing' ][ 'input' ]

        for exp in self.experiments:
            arglist = [ "startyear:" + startyear,  "finalyear:" + finalyear ]
                            
            for key in input_list:
                arglist.append ( key + ':' + exp.B4J_index [ key ])
                            
            for output_item, output_dict in fielddict [ 'preprocessing' ].iteritems ():
                if not isinstance ( output_dict, configobj.Section ): continue
                exp.B4J_index [ output_item ] = self.path_template % { 
                    'path'           : fielddict [ "data_path" ],
                    'classification' : output_item,
                    'expname'        : exp.B4J_name,
                    'startdate'      : fielddict [ 'startyear' ],
                    'finaldate'      : fielddict [ 'finalyear' ],
                    'suffix'         : fielddict [ 'suffix' ]}
                arglist.append ( output_item + ":" + exp.B4J_index [ output_item ])
                                     
            sp = subprocess.Popen( ["python", executable ] +  arglist,
                                   stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
            stdout, stderr = sp.communicate()
            if sp.returncode:
                print "stdout: ", stdout
                print "stderr: ", stderr
                print "returncode: %d" % sp.returncode
                sys.exit ( sp.returncode )
            else:
                print stdout
        self.logger.debug ( "Finished preprocessing with %s for field %s and all experiments successfully." %  \
                            ( executable, fieldname ))
        
        executable = self.get_b4j_path (  fielddict [ 'diagnose' ][ 'exec' ], "scripts" )
        input_list = fielddict [ 'diagnose' ][ 'input' ]

        for can in self.experiments:
            if not can.B4J_iscan: continue
                                   
            for ref in self.experiments:
                if not ref.B4J_isref:            continue
                if ref.B4J_name == can.B4J_name: continue
                
                arglist = [ "startyear:" + startyear,  "finalyear:" + finalyear,
                            "suffix:" +  fielddict [ 'suffix' ],
                            "land_sea_mask:" + land_sea_mask ]
                            
                for key in input_list:
                    arglist.append ( "can.%s" % key + ':' + can.B4J_index [ key ])
                    arglist.append ( "ref.%s" % key + ':' + ref.B4J_index [ key ])
                            
                for output_item, output_dict in fielddict [ 'diagnose' ].iteritems ():
                    if not isinstance ( output_dict, configobj.Section ): continue
                    classification_str = "%s_ref%d" % ( output_item,  ref.B4J_refno )
                    can.B4J_index [ classification_str ] = self.path_template % { 
                        'path'           : fielddict [ "data_path" ],
                        'classification' : classification_str,
                        'expname'        : can.B4J_name,
                        'startdate'      : fielddict [ 'startyear' ],
                        'finaldate'      : fielddict [ 'finalyear' ],
                        'suffix'         : fielddict [ 'suffix' ]}
                    arglist.append ( output_item + ":" + can.B4J_index [ classification_str ])
                                     
                sp = subprocess.Popen( ["python", executable ] +  arglist,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                stdout, stderr = sp.communicate()
                if sp.returncode:
                    print "stdout: ", stdout
                    print "stderr: ", stderr
                    print "returncode: %d" % sp.returncode
                    sys.exit ( sp.returncode )
                else:
                    print stdout
        self.logger.debug ( "Finished diagnostic with %s for field %s and all experiments successfully." %  \
                            ( executable, fieldname ))

        ########### PLOTTING ############################
        
        for p_n, p_dict in fielddict [ 'plotting' ].iteritems ():
            if not isinstance ( p_dict, configobj.Section ): continue
    
            executable = self.get_b4j_path (  p_dict [ 'exec' ], "scripts" )
            input_list = p_dict [ 'input' ]

            if p_dict [ 'nref' ] == "single":
                for can in self.experiments:
                    if not can.B4J_iscan: continue

                    for ref in self.experiments:
                        if not ref.B4J_isref:            continue
                        if ref.B4J_name == can.B4J_name: continue

                        arglist = [ "startyear:" + startyear,  "finalyear:" + finalyear,
                                    "suffix:" +  fielddict [ 'suffix' ],
                                    "land_sea_mask:" + land_sea_mask ]

                        for key in input_list:
                            classification_str =  "%s_ref%d" % ( key,  ref.B4J_refno )
                            arglist.append ( key + ':' + can.B4J_index [ classification_str ])

                        for output_item, output_dict in p_dict.iteritems ():
                            if not isinstance ( output_dict, configobj.Section ): continue
                            classification_str = "%s_ref%d" % ( output_item,  ref.B4J_refno )
                            can.B4J_index [ classification_str ] = self.path_template % { 
                                'path'           : fielddict [ "graphics_path" ],
                                'classification' : classification_str,
                                'expname'        : can.B4J_name,
                                'startdate'      : fielddict [ 'startyear' ],
                                'finaldate'      : fielddict [ 'finalyear' ],
                                'suffix'         : fielddict [ 'plotting' ][ 'suffix' ]}
                            arglist.append ( output_item + ":" + can.B4J_index [ classification_str ])
                            
                        if 'plot_title' in p_dict:
                            try:
                                formats = [ f[0] for f in p_dict [ 'plot_title' ][ 0 ].split ( '%' )[1:]]
                                args = p_dict [ 'plot_title' ][ 1:]
                                if len ( formats ) != len ( args ):
                                    print "ERROR reading plot title"
                                title_args = []
                                for a in args:
                                    exp, mem = a.split ( '_' )
                                    if exp=='can' and mem=='no':
                                        title_args.append ( can.B4J_canno )
                                    elif exp=='can' and mem=='name':
                                        title_args.append ( can.B4J_name )
                                    elif exp=='ref' and mem=='no':
                                        title_args.append ( ref.B4J_canno )
                                    elif exp=='ref' and mem=='name':
                                        title_args.append ( ref.B4J_name )
                                    else:
                                        print "ERROR listing title args"
                                arglist.append ( "plot_title:%s" % p_dict [ 'plot_title' ][ 0 ] % tuple ( title_args ))
                            except:
                                print "ERROR in plotting title."

                        sp = subprocess.Popen( ["python", executable ] +  arglist,
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        stdout, stderr = sp.communicate()
                        print stdout
                        if sp.returncode:
                            print "stdout: ", stdout
                            print "stderr: ", stderr
                            print "returncode: %d" % sp.returncode
                            sys.exit ( sp.returncode )
                        else:
                            print stdout
            elif p_dict [ 'nref' ] == "multiple":
                for can in self.experiments:
                    if not can.B4J_iscan: continue

                    arglist = [ "can_name:" + can.B4J_name, "canno:%d" % can.B4J_canno,
                                "startyear:" + startyear,  "finalyear:" + finalyear ]

                    for ref in self.experiments:
                        if not ref.B4J_isref:            continue
                        if ref.B4J_name == can.B4J_name: continue

                        arglist.append ( "ref%s:" % ( ref.B4J_refno )  + ref.B4J_name  )
                       
                        for key in input_list:
                            classification_str =  "%s_ref%d" % ( key,  ref.B4J_refno )
                            arglist.append ( classification_str + ':' + can.B4J_index [ classification_str ])

                    for output_item, output_dict in p_dict.iteritems ():
                        if not isinstance ( output_dict, configobj.Section ): continue
                        can.B4J_index [ output_item ] = self.path_template % { 
                            'path'           : fielddict [ "graphics_path" ],
                            'classification' : output_item,
                            'expname'        : can.B4J_name,
                            'startdate'      : fielddict [ 'startyear' ],
                            'finaldate'      : fielddict [ 'finalyear' ],
                            'suffix'         : fielddict [ 'plotting' ][ 'suffix' ]}
                        arglist.append ( output_item + ":" + can.B4J_index [ output_item ])

                    sp = subprocess.Popen( ["python", executable ] +  arglist,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    stdout, stderr = sp.communicate()
                    print stdout
                    if sp.returncode:
                        print "stdout: ", stdout
                        print "stderr: ", stderr
                        print "returncode: %d" % sp.returncode
                        sys.exit ( sp.returncode )
                    else:
                        print stdout

            else:
                print ( "ERROR: case '%s' not implemented yet." % p_dict [ 'nref' ])
                sys.exit ( 1 )
                
        self.logger.debug ( "Finished plotting with %s for field %s and all experiments successfully." %  \
                            ( executable, fieldname ))

        self.logger.debug ( "Exiting." )
        return 0
        
