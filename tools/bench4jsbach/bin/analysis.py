#! /usr/bin/env python

import os
import fnmatch
from cdo import Cdo
from sets import Set

from pycmbs.data import Data
import matplotlib.pyplot as plt
from pycmbs.mapping import map_plot


# methods #
def find_files ( pattern, path ):
    assert os.path.isdir ( path ) and isinstance ( pattern, str )
    
    ret = []
    for root, dirnames, filenames in os.walk ( path ):
        files   = fnmatch.filter ( filenames, pattern )
        for f in files:
            ret.append ( os.path.join ( root, f ))

    return ret
                    
    # albedo
def get_period_intersection ( filelist ):
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


def albedo_yseasmean ( rsus_infile, rsds_infile, albedo_outfile, years ):
    assert os.path.isfile ( rsus_infile ), "Input file %s not found." % rsus_infile
    assert os.path.isfile ( rsds_infile ), "Input file %s not found." % rsds_infile

    cdo = Cdo ()
    
    assert Set ( years ).issubset ( Set ( map ( int, ''.join ( cdo.showyear(
        input=rsus_infile )). encode (). split ()))), "Selected years not in " +\
        "range of simulation period."

    cdo.setname ( "albedo",
                  input="-setunit,1 " +\
                  "-setvrange,0.0001,0.99 " +\
                  "-div " +\
                  ( "-yseasmean " + \
                    ( " -selyear," + ','.join ( repr ( y ) for y in years ) + \
                      " " +  rsus_infile )) +\
                  ( " -yseasmean " + \
                    ( " -selyear," + ','.join ( repr ( y ) for y in years ) + \
                      " " +  rsds_infile )),
                  output=albedo_outfile,
                  options='-f nc -m 0.00011' )

    # TO CHECK:
    # cdo setvrange,min,max sets x in (min,max) to missval # ??? the opposite to 
    # what you would expect
    # cdo setrtomiss,min,max sets x not in (min,max) to missval or (???) but also changes
    # values in that range
    
    #  "-setvrange,0.0001,0.99 " +\
    #  "-setrtomiss,0.0001,0.91 " +\
    #  "-setctomiss,1 " +\


class Analysis ( object ):
    def __init__ ( self, **kwargs ):
        pass


    def albedo_analyse ( self ):
        """
        - look for albedo output, if it is not found look for rsus, rsds output and calculate albedo
        - calculate seasonal means for largest common period and store plot with subplots
        - calculate difference to reference experiment and store plot with subplots
        - storage demand over exp loop, temp: subplot, data
        - storage demand for report: paths for plots, report section texts
        - TODO: rename method to 'albedo seasonal means'
        """
        self.logger.debug ( "Starting." )

        from bench import REF_EXP_IDX
        from bench import Bench_Object
        assert isinstance ( self, Bench_Object ), \
               "This method has to be used by a Bench_Object."

        template = """%(path)s/%(varname)s_%(classification)s_%(expname)s_%(startdate)s-%(finaldate)s.%(suffix)s"""

        varname = "albedo"

        path_tokens = [ "rsus","rsds","albedo","DJF","MAM","JJA","SON",
                        "DJF_diff","MAM_diff","JJA_diff","SON_diff" ]
        n_paths = len ( path_tokens )
        
        cdo = Cdo ()

        self.analyses [ varname ] = { "type" : "yseasmean_diff",
                                      "plots" : [],
                                      "diff_plots" : [],
                                      "text" : "" }
        
        n_exp = len ( self.experiments )
        assert n_exp > 1    # TOCHECK: this should be checked while reading bench configfile 
        self.logger.warning ( "Experiments' former filename list will be deleted." )
        
        # TODO: look for albedo file first ...
        # ... no albedo found. Calculate albedo:
        for e in self.experiments:
            e.B4J_buffer = [ None ] * n_paths
            for root, dirnames, filenames in os.walk ( e.B4J_path ):
                for idx in [ path_tokens.index ( "rsus" ),
                             path_tokens.index ( "rsds" )]:
                    pattern = '*%s*%s*.nc' % ( path_tokens [ idx ], e.B4J_name )
                    files   = fnmatch.filter ( filenames, pattern )
                    if len ( files ) != 1:
                        raise Exception ( """Multiple or no files found matching pattern '%s' in '%s':
                                        files = [ %s ]""" %\
                                          ( pattern, e.B4J_path,
                                            ",".join ( "{!s}".format ( f ) for f in files )))
                    e.B4J_buffer [ idx ] =  os.path.join ( root, files [ 0 ])
            self.logger.debug ( e )

        period = get_period_intersection ([ p for e in self.experiments \
                                            for p in e.B4J_buffer if p ])

        seasons      = [ 'DJF', 'MAM', 'JJA', 'SON' ]
        n_seasons    = len ( seasons )
        for e in self.experiments:
            e.B4J_buffer [ path_tokens.index ( varname )] = template % {
                'path'           : self.data_path,
                'varname'        : varname,
                'classification' : 'yseasmean',
                'expname'        : e.B4J_name,
                'startdate'      : str ( period [ 0 ]),
                'finaldate'      : str ( period [ -1 ]),
                'suffix'         : 'nc' }
             
            albedo_yseasmean ( e.B4J_buffer [ path_tokens.index ( "rsus" )],
                               e.B4J_buffer [ path_tokens.index ( "rsds" )],
                               e.B4J_buffer [ path_tokens.index ( "albedo" )],
                               period )


            f            = plt.figure ( figsize = ( 15, 10 ))
            tmp_subplots = [ f.add_subplot ( 2, 2, idx + 1 ) for idx in range ( n_seasons )]
            figurefilepath =  template % { 'path'           : self.graphic_path,
                                           'varname'        : varname,
                                           'classification' : 'yseasmean',
                                           'expname'        : e.B4J_name,
                                           'startdate'      : str ( period [ 0 ]),
                                           'finaldate'      : str ( period [ -1 ]),
                                           'suffix'         : 'png' }
            self.analyses [ varname ][ "plots" ].append ( figurefilepath )
            for season in seasons:
                season_idx = seasons.index ( season )
                buffer_idx = path_tokens.index ( season )
                tmp_data     = [ None ] * n_seasons
                self.logger.debug ( "season_idx=%d, buffer_idx=%d" % (season_idx, buffer_idx ))
                e.B4J_buffer [ buffer_idx ] = template % { 'path'           : self.data_path,
                                                       'varname'        : varname,
                                                       'classification' : 'y%smean' % season,
                                                       'expname'        : e.B4J_name,
                                                       'startdate'      : str ( period [ 0 ]),
                                                       'finaldate'      : str ( period [ -1 ]),
                                                       'suffix'         : 'nc' }
                cdo.selseas ( season,
                              input  = e.B4J_buffer [ path_tokens.index ( "albedo" )],
                              output = e.B4J_buffer [ buffer_idx ])
                tmp_data [ season_idx ] = Data (  e.B4J_buffer [ buffer_idx ],
                                                  'albedo',
                                                  lat_name = 'lat',
                                                  lon_name = 'lon',
                                                  read     = True,
                                                  label    = season )
                
                x = map_plot ( tmp_data [ season_idx ],
                               ax = tmp_subplots [ season_idx ],
                               use_basemap = True )

            plt.show ()
            plt.savefig ( figurefilepath )
            plt.close ()


            ### TODO: external loop
        self.analyses [ varname ][ "diff_plots" ] = [ None ] * n_exp
        for e in self.experiments:
            exp_idx = self.experiments.index ( e )
            if exp_idx == REF_EXP_IDX:
                continue
            
            self.logger.debug ( "Test experiment '%s', make diff plots." % e.B4J_name )
            f            = plt.figure ( figsize = ( 15, 10 ))
            tmp_subplots = [ f.add_subplot ( 2, 2, idx + 1 ) for idx in range ( n_seasons )]
            figurefilepath =  template % { 'path'           : self.graphic_path,
                                           'varname'        : varname,
                                           'classification' : 'yseasmean_diff',
                                           'expname'        : e.B4J_name,
                                           'startdate'      : str ( period [ 0 ]),
                                           'finaldate'      : str ( period [ -1 ]),
                                           'suffix'         : 'png' }
            self.analyses [ varname ][ "diff_plots" ][ exp_idx ] = figurefilepath
            for season in seasons:
                season_idx = seasons.index ( season )
                buffer_in_idx =  path_tokens.index ( season )
                buffer_out_idx = path_tokens.index ( season + "_diff" )
                tmp_data     = [ None ] * n_seasons
                e.B4J_buffer [ buffer_out_idx ] = template % { 'path'           : self.data_path,
                                                           'varname'        : varname,
                                                           'classification' : 'y%smean_diff' % season,
                                                           'expname'        : e.B4J_name,
                                                           'startdate'      : str ( period [ 0 ]),
                                                           'finaldate'      : str ( period [ -1 ]),
                                                           'suffix'         : 'nc' }

                cdo.sub ( input   = e.B4J_buffer [ buffer_in_idx ] + " " +\
                          self.experiments [ REF_EXP_IDX ].B4J_buffer [ buffer_in_idx ],
                          output  = e.B4J_buffer [ buffer_out_idx ],
                          options = '-f nc' )

                self.logger.debug ( "input sub: " + e.B4J_buffer [ season_idx ] + \
                                   " " + self.experiments [ REF_EXP_IDX ].B4J_buffer [ season_idx ] )

                tmp_data [ season_idx ] = Data (  e.B4J_buffer [ buffer_out_idx ],
                                                  'albedo',
                                                  lat_name = 'lat',
                                                  lon_name = 'lon',
                                                  read     = True,
                                                  label    = season + " diff")
                


                x = map_plot ( tmp_data [ season_idx ],
                               ax = tmp_subplots [ season_idx ],
                               use_basemap = True )

            plt.show ()
            plt.savefig ( figurefilepath )
            plt.close ()

        
        self.logger.debug ( "Exiting." )

    
