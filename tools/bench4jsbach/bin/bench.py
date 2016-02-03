#! /usr/bin/env python

import os, subprocess, shutil,getpass,datetime
from experiment import Experiment, CMIP5_Experiment
from util import logger as dkl_logger, setup_logger, print_none, log_message_newline, str2bool
import configobj
from sets import Set
from analysis import Analysis
from report import Report


REF_EXP_IDX       = 0
SWITCHED_ON_IDX   = TYPE_IDX = 0
TEMP_SCALE_IDX    = PATH_IDX = 1
ISREF_IDX         = 2
ISCAN_IDX         = 3
VARNAME_4_CMIP5_METADATA = 'rsds'
bench4jsbach_base =  os.environ [ 'BENCH4JSBACH_BASE' ]
if not bench4jsbach_base: raise Exception ( "Environment variable not set." )
config_template   = os.path.join ( bench4jsbach_base, "configuration", "template.cfg" )
work_base         = os.path.join ( bench4jsbach_base, "work" )
bench_dirnames    = [ "config", "data", "graphic", "report", "logging" ]


class Bench_Object ( Analysis, Report ):
    bench_name        = None
    bench_path        = None
    config_path       = None
    data_path         = None
    graphic_path      = None
    report_path       = None
    logging_path      = None
    logfile_path      = None
    bench_config_path = None
    config            = None
    experiments       = []
    logger            = None
    logger_name       = 'B4J_LOGGER'
    reporting         = True
    
    
    def __init__ ( self, config_arg, **kwargs ):
        dkl_logger.debug ( "Starting." )
        
        if not os.path.exists ( config_arg ):
            self._writeconfig__init__ ( config_arg )
            dkl_logger.debug ( "Exiting." )
            return
        bench_config_path = config_arg
        
        # read config object
        try:
            config = configobj.ConfigObj ( bench_config_path )
        except ( configobj.ConfigObjError, IOError), e:
            print 'Could not read "%s": %s' % ( bench_config_path, e)
        dkl_logger.debug ( "Read configuration file %s." % ( bench_config_path ))
            
        # TODO: validation. First shot:
        try:
            bench_name = config [ "bench_name" ]
        except:
            raise Exception (  "ERROR reading token 'bench_name' from configuration file %s." % \
                               bench_config_path )
        
        # check environment
        bench_config_path = os.path.abspath ( bench_config_path )
        config_path =  os.path.dirname ( bench_config_path ) 
        bench_path = os.path.dirname ( config_path )
        if not os.path.basename ( bench_path ) == bench_name:
            raise Exception ( "Directory structure found differs from assumption." )
        onlydirs = [ d for d in os.listdir ( bench_path ) \
                     if os.path.isdir ( os.path.join ( bench_path, d ))]
        if not Set ( onlydirs ).issuperset ( Set ( bench_dirnames )):
            raise Exception ( "Directory structure found differs from assumption." )

        # all tests succeeded
        self.bench_name        = bench_name
        self.bench_path        = bench_path
        self.bench_config_path = bench_config_path
        self.config            = config
        for n in bench_dirnames:
            path = os.path.join ( self.bench_path, n )
            setattr ( self,  n + "_path",  path )

        self.logfile_path = os.path.join ( self.logging_path,
                                           self.bench_name + '_' +\
                                           str ( datetime.datetime.now ().strftime ( '%Y%m%d%H%M' )) +\
                                           '.out' )

        dkl_logger.debug ( "Logger switched to %s. It will write to file " % self.logger_name + \
                           log_message_newline +\
                           "'%s' of now." % self.logfile_path +\
                           log_message_newline +\
                           "Processing might take a while ..." )

        self.logger = setup_logger ( name = self.logger_name, filename = self.logfile_path )

        
        Analysis.__init__ ( self )
        self.reporting = str2bool ( config [ "reporting" ])
        Report.__init__ ( self ) if self.reporting else None
            
            
        es = self.config [ "experiments" ]
        from analysis import find_files
        n_refs = n_cans = 0
        for en in es:
            name = en
            path = es [ en ][ PATH_IDX ]
            iscan = str2bool ( es [ en ][ ISCAN_IDX ])
            isref = str2bool ( es [ en ][ ISREF_IDX ])
            n_cans += 1 if iscan else 0
            n_refs += 1 if isref else 0
            refno = n_refs - 1 if isref else None
            kwargs = { 'B4J_name'  : name,
                       'B4J_etype' : es [ en ][ TYPE_IDX ],
                       'B4J_path'  : path,
                       'B4J_iscan' : iscan,
                       'B4J_isref' : isref,
                       'B4J_refno' : refno,
                       'B4J_bench' : self }
            if es [ en ][ TYPE_IDX ] == "CMIP5":
                varfile = find_files ( "*%s*%s*%s" % \
                                        ( VARNAME_4_CMIP5_METADATA, name, '.nc' ),
                                        path, raise_e = True )
                kwargs [ 'var_filename' ] = varfile
                e = CMIP5_Experiment ( **kwargs )
            else:
                e = Experiment ( **kwargs )
            self.experiments.append ( e )
        if n_refs < 1 or n_cans < 1: raise Exception ( "No candidate or no reference experiment defined."  )

        self.logger.debug ( "Exiting." )


    def _writeconfig__init__ ( self, config_arg ):
        dkl_logger.debug ( "Starting." )

        self.bench_name = config_arg
        self.bench_path = os.path.join ( work_base, self.bench_name )
        if os.path.exists ( self.bench_path ):
            shutil.rmtree ( self.bench_path )
            #raise Exception ( "Directory %s already exists." % self.bench_path )
        for n in bench_dirnames:
            path = os.path.join ( self.bench_path, n )
            setattr ( self,  n + "_path",  path )
            os.makedirs ( path )

        self.bench_config_path = os.path.join ( self.config_path, self.bench_name + ".cfg" )
        with open ( config_template ) as f:
            template_contents = f.read ()
        config_contents = template_contents % { '<BENCH_NAME>': self.bench_name,
                                                '<AUTHOR_NAME>': getpass.getuser (),
                                                '<TODAY>': str ( datetime.date.today ())}
        with open ( self.bench_config_path, "w" ) as f:
            f.write ( config_contents )

        dkl_logger.debug ( "Exiting." )

        
    def __repr__ ( self ):
        return  """{} ( bench_name = {}, bench_path = {}, config_path = {},
        data_path = {}, graphic_path = {}, report_path = {}, bench_config_path = {}, config = {}, 
        experiments = {}""".format (
            self.__class__.__name__,
            print_none (  self.bench_name        ),
            print_none (  self.bench_path        ),
            print_none (  self.config_path       ),
            print_none (  self.data_path         ),
            print_none (  self.graphic_path      ),
            print_none (  self.report_path       ),
            print_none (  self.bench_config_path ),
            print_none (  self.config            ),
            print_none ( self.experiments )).replace ( "\n", "" )
    

