#! /usr/bin/env python
# -*- coding:  utf-8 -*-
"""The module contains diagnostic related classes, functionality and data."""

from __future__ import division

import os, sys
import fnmatch, re
import warnings
from datetime import datetime
import configobj
from sets import Set
from collections import OrderedDict


import numpy as np
from cdo import Cdo

from util import (
    logger as dkl_logger,
    print_item_or_list_rec,
    subprocessing,
    touch )
from experiment import (
    MPIESM1_Experiment,
    CMIP5_Experiment )
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


class Diagnostics ( object ):
    # TODO: take varname out as it is already in path
    path_template = """%(path)s/%(classification)s_%(expname)s_%(startdate)s-%(enddate)s.%(suffix)s"""


    def __init__ ( self, **kwargs ):
        super ( Diagnostics, self ).__init__()
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
        ref_etypes = list ( set ([ e.B4J_etype for e in self.experiments if e.B4J_isref ])) 
        
        cdo = Cdo ()
       
        
        ################ PREPROCESSING ####################

        for exp in self.experiments:
            # TODO: validation of field configuration in preliminary step, EXTEND VALIDATION PROCESS
            assert 'preprocessing' in fielddict and \
                exp.B4J_etype in fielddict [ 'preprocessing' ], \
                "No 'preprocessing'-section defined for %s and %s." % ( fieldname, exp.B4J_etype )
                
            prepdict = fielddict [ 'preprocessing' ][ exp.B4J_etype ]            
            executable = self.get_b4j_path (  prepdict [ 'exec' ], "scripts" )
            input_list = prepdict [ 'input' ]
            
            #arglist = [ "startyear:" + exp.B4J_startyear,  "finalyear:" + exp.B4J_endyear ]
            arglist = [ "fieldname:" + fieldname,
                        "suffix:" +  fielddict [ 'suffix' ],
                        "land_sea_mask:" + exp.land_sea_mask ]
                            
            for key in input_list:
                arglist.append ( key + ':' + exp.B4J_index [ key ])
                
            for output_item in  fielddict [ 'preprocessing' ][ 'output' ]:
                exp.B4J_index [ output_item ] = self.path_template % { 
                    'path'           : fielddict [ "processed_path" ],
                    'classification' : output_item,
                    'expname'        : exp.B4J_name,
                    'startdate'      : exp.B4J_startyear,
                    'enddate'        : exp.B4J_endyear,
                    'suffix'         : fielddict [ 'suffix' ]}
                arglist.append ( output_item + ":" + exp.B4J_index [ output_item ])

            cmd = ["python", executable ] +  arglist
            self.logger.debug ( """Calling childprocess preprocessing: \n%s""" % cmd )
            subprocessing ( cmd, shell_arg = False, stderr_arg = None )

            # Check written files:
            new_indices = [ output_item for output_item in  fielddict [ 'preprocessing' ][ 'output' ]]
            print new_indices
            exp.index_check_output ( new_indices )
        
        self.logger.debug ( "Finished preprocessing with %s for field %s and all experiments successfully." %  \
                            ( executable, fieldname ))

        print "######################## PREPROCESSING FINISHED "

        ########### REMAPPING FOR CANDIDATES #############################

        self.logger.debug ( "Remapping: check necessity and potential processing." )

        for can in self.experiments:
            if not can.B4J_iscan: continue
                
            remap_targets = [ etype for etype in ref_etypes if etype != can.B4J_etype  ]
            
            if remap_targets:
                for output_item in  fielddict [ 'preprocessing' ][ 'output' ]:
                    for remap_target in remap_targets:
                        output_item_remap = output_item + "_remap-" + remap_target
                        can.B4J_index [ output_item_remap ] = self.path_template % { 
                            'path'           : fielddict [ "processed_path" ],
                            'classification' : output_item_remap,
                            'expname'        : can.B4J_name,
                            'startdate'      : can.B4J_startyear,
                            'enddate'        : can.B4J_endyear,
                            'suffix'         : fielddict [ 'suffix' ]}
                        
                        if remap_target == "MPIESM1":
                            gridstring = "t63grid"
                            print "WARNING: fixed remapping to gaussian %s." % gridstring

                            print "Start remapping os %s to %s. This might take a while ..." % \
                                ( can.B4J_index [ output_item ], gridstring )
                            sys.stdout.flush()

                            #HACK TO SAVE COMPUTING TIME
                            if not os.path.isfile (  can.B4J_index [ output_item_remap ]):
                                cdo.remapycon ( gridstring,
                                                input = can.B4J_index [ output_item ],
                                                output = can.B4J_index [ output_item_remap ])
                        else:
                            print "ERROR: remapping for %s not implemented yet." % remap_target
                            sys.exit ( 1 )

                new_indices = [ "%s_remap-%s" % ( output_item, remap_target ) \
                                for output_item in fielddict [ 'preprocessing' ][ 'output' ] \
                                for remap_target in remap_targets ]
                print new_indices
                can.index_check_output ( new_indices )

            
        self.logger.debug ( "Finished remapping." )

        ########### DIAGNOSTIC ############################
        
        self.logger.debug ( "Start diagnostic processing." )
        
        executable = self.get_b4j_path (  fielddict [ 'diagnose' ][ 'exec' ], "scripts" )
        input_list = fielddict [ 'diagnose' ][ 'input' ]

        for can in self.experiments:
            if not can.B4J_iscan: continue
            print
            print print_item_or_list_rec (  can.B4J_index, 3 )
                                   
            for ref in self.experiments:
                if not ref.B4J_isref:            continue
                if ref.B4J_name == can.B4J_name: continue
                print
                print print_item_or_list_rec ( ref.B4J_index, 3 )
             
                arglist = [ "fieldname:" + fieldname,
                            "suffix:" +  fielddict [ 'suffix' ],
                            "land_sea_mask:" + ref.land_sea_mask ]
                            
                for key in input_list:
                    if can.B4J_etype != ref.B4J_etype:
                        gridstring = "t63grid"
                        print "WARNING: fixed remapping to gaussian %s." % gridstring
                        remapped = key + "_remap-" + ref.B4J_etype
                        arglist.append ( "can.%s" % key + ':' + can.B4J_index [ remapped ])
                    else:
                        arglist.append ( "can.%s" % key + ':' + can.B4J_index [ key ])
                    arglist.append ( "ref.%s" % key + ':' + ref.B4J_index [ key ])
                            
                for output_item, output_dict in fielddict [ 'diagnose' ].iteritems ():
                    if not isinstance ( output_dict, configobj.Section ): continue
                    classification_str = "%s_ref%d" % ( output_item,  ref.B4J_refno )
                    can.B4J_index [ classification_str ] = self.path_template % { 
                        'path'           : fielddict [ "processed_path" ],
                        'classification' : classification_str,
                        'expname'        : can.B4J_name,
                        'startdate'      : can.B4J_startyear,
                        'enddate'        : can.B4J_endyear,
                        'suffix'         : fielddict [ 'suffix' ]}
                    arglist.append ( output_item + ":" + can.B4J_index [ classification_str ])

                cmd = ["python", executable ] +  arglist
                subprocessing ( cmd, shell_arg = False, stderr_arg = None )

                new_indices = [ "%s_ref%d" % ( output_item, ref.B4J_refno ) \
                                for output_item, output_dict in fielddict [ 'diagnose' ].iteritems () \
                                if  isinstance ( output_dict, configobj.Section ) \
                                for ref in self.experiments if ref.B4J_isref ]
                print new_indices
                can.index_check_output ( new_indices )
                
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

                        arglist = [ "fieldname:" + fieldname,
                                    "suffix:" +  fielddict [ 'suffix' ],
                                    "land_sea_mask:" + can.land_sea_mask ]
                            
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
                                'startdate'      : can.B4J_startyear,
                                'enddate'        : can.B4J_endyear,
                                'suffix'         : fielddict [ 'plotting' ][ 'suffix' ]}
                        
                            if "plot_title" in output_dict:
                                try:
                                    formats = [ f[0] for f in output_dict [ 'plot_title' ][ 0 ].split ( '%' )[1:]]
                                    args = output_dict [ 'plot_title' ][ 1:]
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
                                            title_args.append ( ref.B4J_refno )
                                        elif exp=='ref' and mem=='name':
                                            title_args.append ( ref.B4J_name )
                                        else:
                                            print "ERROR listing title args"
                                    arglist.append ( "plot_title_%s:%s" % ( output_item,
                                                                            output_dict [ 'plot_title' ][ 0 ] % \
                                                                            tuple ( title_args )))
                                except:
                                    print "ERROR in plotting title."


                            arglist.append ( output_item + ":" + can.B4J_index [ classification_str ])

                        cmd = ["python", executable ] +  arglist
                        subprocessing ( cmd, shell_arg = False, stderr_arg = None )
                        

            elif p_dict [ 'nref' ] == "None":
                for can in self.experiments:
                    if not can.B4J_iscan: continue

                    arglist = [ "fieldname:" + fieldname,
                                "suffix:" +  fielddict [ 'suffix' ],
                                "land_sea_mask:" + can.land_sea_mask,
                                "can_name:" + can.B4J_name,
                                "canno:%d" % can.B4J_canno ]
                    
                    for key in input_list:
                        arglist.append ( key + ':' + can.B4J_index [ key ])
                        
                    for output_item, output_dict in p_dict.iteritems ():
                        if not isinstance ( output_dict, configobj.Section ): continue
                        classification_str = output_item
                        can.B4J_index [ classification_str ] = self.path_template % { 
                            'path'           : fielddict [ "graphics_path" ],
                            'classification' : classification_str,
                            'expname'        : can.B4J_name,
                            'startdate'      : can.B4J_startyear,
                            'enddate'        : can.B4J_endyear,
                            'suffix'         : fielddict [ 'plotting' ][ 'suffix' ]}
                        
                        if "plot_title" in output_dict:                      
                            try:
                                formats = [ f[0] for f in output_dict [ 'plot_title' ][ 0 ].split ( '%' )[1:]]
                                args = output_dict [ 'plot_title' ][ 1:]
                                if len ( formats ) != len ( args ):
                                    print "ERROR reading plot title"
                                title_args = []
                                for a in args:
                                    exp, mem = a.split ( '_' )
                                    if exp=='can' and mem=='no':
                                        title_args.append ( can.B4J_canno )
                                    elif exp=='can' and mem=='name':
                                        title_args.append ( can.B4J_name )
                                    else:
                                        print "ERROR listing title args"

                                arglist.append ( "plot_title_%s:%s" % ( classification_str,
                                                                        output_dict [ 'plot_title' ][ 0 ] % \
                                                                        tuple ( title_args )))
                            except:
                                print "ERROR in plotting title."

                        arglist.append ( output_item + ":" + can.B4J_index [ classification_str ])
                        
                    cmd = ["python", executable ] +  arglist
                    subprocessing ( cmd, shell_arg = False, stderr_arg = None )

                    
            elif p_dict [ 'nref' ] == "multiple":
                for can in self.experiments:
                    if not can.B4J_iscan: continue

                    arglist = [ "can_name:" + can.B4J_name, "canno:%d" % can.B4J_canno ]

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
                            'startdate'      : can.B4J_startyear,
                            'enddate'        : can.B4J_endyear,
                            'suffix'         : fielddict [ 'plotting' ][ 'suffix' ]}
                        
                        print "############## 1 ##############"
                        print output_item
                        print output_dict
                        print "#############################"

                        arglist.append ( output_item + ":" + can.B4J_index [ output_item ])

                    cmd = ["python", executable ] +  arglist
                    subprocessing ( cmd, shell_arg = False, stderr_arg = None )

                    
            elif p_dict [ 'nref' ] == "combined":
                for can in self.experiments:
                    if not can.B4J_iscan: continue

                    arglist = [ "field_name:" + fielddict [ "field_name" ],
                                "can_name:" + can.B4J_name,
                                "canno:%d" % can.B4J_canno ]
                    
                    for classification_str in p_dict [ 'input_can' ]:
                        assert os.path.isfile ( can.B4J_index [ classification_str ])
                        arglist.append ( classification_str + ':' + can.B4J_index [ classification_str ])

                    for ref in self.experiments:
                        if not ref.B4J_isref:            continue
                        if ref.B4J_name == can.B4J_name: continue

                        arglist.append ( "ref%s:" % ( ref.B4J_refno )  + ref.B4J_name  )
                       
                        for key in input_list:
                            classification_str =  "%s_ref%d" % ( key,  ref.B4J_refno )
                            assert os.path.isfile ( can.B4J_index [ classification_str ])
                            arglist.append ( classification_str + ':' + can.B4J_index [ classification_str ])

                    for key in [ 'varname', 'unit', 'varname_can', 'unit_can' ]:
                        if key in p_dict:
                            arglist.append ( '%s:%s' % ( key, p_dict [ key ]))
                        
                    for output_item, output_dict in p_dict.iteritems ():
                        if not isinstance ( output_dict, configobj.Section ): continue
                        classification_str = output_item
                        can.B4J_index [ classification_str ] = self.path_template % { 
                            'path'           : fielddict [ "graphics_path" ],
                            'classification' : classification_str,
                            'expname'        : can.B4J_name,
                            'startdate'      : can.B4J_startyear,
                            'enddate'        : can.B4J_endyear,
                            'suffix'         : fielddict [ 'plotting' ][ 'suffix' ]}
                        
                        if 'plot_title' in output_dict:
                            output_dict.main.interpolation = False
                            arglist.append ( "plot_title_%s:" % ( classification_str ) + \
                                             output_dict [ 'plot_title' ])

                        if 'legend_tmpl' in output_dict:
                            output_dict.main.interpolation = False
                            arglist.append ( "legend_tmpl_%s:" % ( classification_str ) + \
                                             output_dict [ 'legend_tmpl' ])


                        arglist.append ( classification_str + ":" + can.B4J_index [ classification_str ])


                    cmd = ["python", executable ] +  arglist
                    subprocessing ( cmd, shell_arg = False, stderr_arg = None )

            else:
                print ( "ERROR: case '%s' not implemented yet." % p_dict [ 'nref' ])
                sys.exit ( 1 )

            
            # Check written files:
            # new_indices = [ output_item for output_item in  fielddict [ 'preprocessing' ][ 'output' ]]
            # print new_indices
            # exp.index_check_output ( new_indices )
                
        self.logger.debug ( "Finished plotting with %s for field %s and all experiments successfully." %  \
                            ( executable, fieldname ))

        self.logger.debug ( "Exiting." )
        
        return 0
        
    def __repr__ ( self ):
        super ( Diagnostics, self ).__repr__ ()
        rep =  """\n\nSuperclass object {} ()\n\n""".\
              format ( Diagnostics.__name__ )
        
        return rep
