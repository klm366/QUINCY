#! /usr/bin/env python

import os, sys, subprocess, shutil, getpass, datetime
import copy
import configobj
from validate import Validator
from sets import Set
import subprocess

from wikipage import Bench4JSBach_Wikipage
from experiment import Experiment, CMIP5_Experiment, MPIESM1_Experiment
from util import ( logger as dkl_logger, setup_logger, print_none,
                   log_message_newline, str2bool,print_item_or_list_rec,
                   gethostname )
from diagnostics import Diagnostics, find_files
from report import Report

#from datetime import datetime

LOGSTREAM_2_FILE = False

REF_EXP_IDX              = 0
SWITCHED_ON_IDX          = TYPE_IDX = 0
TEMP_SCALE_IDX           = PATH_IDX = 1
ISREF_IDX                = 2
ISCAN_IDX                = 3
VARNAME_4_CMIP5_METADATA = 'rsds'

bench4jsbach_base          =  os.environ [ 'BENCH4JSBACH_BASE' ]
if not bench4jsbach_base: raise Exception ( "Environment variable not set." )
bench4jsbach_configs       =  os.path.join ( bench4jsbach_base, "configs" )
bench4jsbach_scripts       =  os.path.join ( bench4jsbach_base, "scripts" )
config_template            = os.path.join ( bench4jsbach_configs, "b4j_template.cfg" )
b4j_configspec_fn          = os.path.join ( bench4jsbach_configs, "b4j_config.spec" )
b4j_field_configspec_fn    = os.path.join ( bench4jsbach_configs, "b4j_field_config.spec" )
work_base                  = os.path.join ( bench4jsbach_base, "work" )
bench_dirnames             = [ "configs",
                               "downloads",
                               "graphics",
                               "interpolated_selected",
                               "logging",
                               "processed",
                               "report",
                               "scripts" ]


def get_git_origin ():
    p = subprocess.Popen ([ "git", "remote", "-v" ], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE )
    out, err = p.communicate()

    if not err:
        strl=out.split('\n')
        strm=[ s.split('\t') for s in strl if '\t' in s]
        strdsl=[[l[0]]+l[1].split() for l in strm]
        ret=[l[0]+l[2]+':'+l[1] for l in strdsl]
        ret=','.join(s for s in ret) + ','
    else:
        ret = "error,"

    return ret

    
def get_git_branch ():
    p = subprocess.Popen ([ "git", "branch" ], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE )
    out, err = p.communicate()
    strl=out.split('\n')
    return [ s.strip ( ' *' ) for s in strl if '*' in s ][ 0 ]

    
def get_git_sha ():
    p = subprocess.Popen ([ "git", "rev-parse", "HEAD" ], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE )
    out, err = p.communicate()
    return out.strip() if not err else err.strip()


def get_git_status ():
    p = subprocess.Popen ([ "git", "status", "--short" ], stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE )
    out, err = p.communicate()
    if not err:
        strl=out.strip().split('\n')
        modified_files = [ s.strip ( ' M' ) for s in strl if 'M' in s ]
        ret = "modified" if modified_files else "clean"
    else:
        ret = "error"
    return ret

def get_ip_address ():
    import socket
    return socket.gethostbyname ( socket.gethostname ())


class Bench_Object ( Diagnostics, Report ):
    bench_name         = None
    bench_path         = None
    bench_configs_path = None
    configs_path       = None
    processed_path     = None
    graphics_path      = None
    logging_path       = None
    logfile_path       = None
    report_path        = None
    scripts_path       = None
    download_path      = None
    config             = None
    wikipage           = None
    experiments        = []
    n_cans             = 0
    n_refs             = 0
    logger             = None
    logger_name        = 'B4J_LOGGER'
    reporting          = False
    
    
    def __init__ ( self, config_arg, **kwargs ):
        super ( Bench_Object, self ).__init__ ()
        
        dkl_logger.debug ( "Starting." )
        
        if not os.path.exists ( config_arg ):
            self._writeconfig__init__ ( config_arg )
            dkl_logger.debug ( "Exiting." )
            return
        bench_configs_path = config_arg
 
        # read configspec and config
        if not os.path.isfile ( b4j_configspec_fn ):
            raise Exception ( "\nERROR no file named %s." %  \
                              b4j_configspec_fn )
        try:
            configspec = configobj.ConfigObj ( b4j_configspec_fn,
                                               interpolation=False,
                                               list_values=False,
                                               _inspec=True)
        except ( configobj.ConfigObjError, IOError), e:
            print 'Could not read "%s": %s' % ( b4j_configspec_fn, e )
            
        try:
            config = configobj.ConfigObj ( bench_configs_path,
                                           configspec=configspec )
        except ( configobj.ConfigObjError, IOError), e:
            print '\nCould not read "%s": %s' % ( bench_configs_path, e)
            
        dkl_logger.debug ( "Read configuration file %s." % ( bench_configs_path ))

        val  = Validator()
        test = config.validate ( val )
        if test == True:          
            dkl_logger.debug ( "\nValidation of configuration from %s with specification from %s successfully run." % \
                ( bench_configs_path, b4j_configspec_fn ))
        else:
            raise Exception (  "\nERROR validating configuration %s with specialization %s:\n%s" %\
                               ( bench_configs_path, b4j_configspec_fn,
                                 ', '.join ( "item '%s'" %  k for k,v in test.items () if v == False)))
            sys.exit(1)
            
        bench_name = config [ "bench_name" ]
        
        # check environment
        bench_configs_path = os.path.abspath ( bench_configs_path )
        configs_path =  os.path.dirname ( bench_configs_path ) 
        bench_path = os.path.dirname ( configs_path )
        if not os.path.basename ( bench_path ) == bench_name:
            raise Exception ( "Directory structure found differs from assumption." )
            sys.exit ( 1 )
        onlydirs = [ d for d in os.listdir ( bench_path ) \
                     if os.path.isdir ( os.path.join ( bench_path, d ))]
        if not Set ( onlydirs ).issuperset ( Set ( bench_dirnames )):
            raise Exception ( "Directory structure found differs from assumption." )
            sys.exit ( 1 )

        # all tests succeeded
        self.bench_name        = bench_name
        self.bench_path        = bench_path
        self.bench_configs_path = bench_configs_path
        self.config            = config
        for n in bench_dirnames:
            path = os.path.join ( self.bench_path, n )
            setattr ( self,  n + "_path",  path )

        if LOGSTREAM_2_FILE:
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

        self.reporting     = config [ "reporting" ]
        self.wikipage =  Bench4JSBach_Wikipage ( self.config [ 'redmine_user' ],
                                                 self.config [ 'redmine_key' ])

        es = self.config [ "experiments" ]
        n_refs = n_cans = 0
        for en, attr in es.iteritems ():
            name = en
            name_orig = attr.get ( "name_orig", None )
            etype = attr.get ( "type", "CMIP5" )
            path = attr.get ( "path", None )
            iscan = attr.get ( 'iscan', True )
            isref = attr.get ( 'isref', False )
            startyear = attr.get ( 'startyear', None )
            finalyear = attr.get ( 'finalyear', None )
            n_cans += 1 if iscan else 0
            n_refs += 1 if isref else 0
            canno = n_cans - 1 if iscan else None
            refno = n_refs - 1 if isref else None
            kwargs = { 'B4J_name'         : name,
                       'B4J_name_orig'    : name_orig,
                       'B4J_etype'        : etype,
                       'B4J_path'         : path,
                       'B4J_iscan'        : iscan,
                       'B4J_isref'        : isref,
                       'B4J_startyear'    : startyear,
                       'B4J_finalyear'    : finalyear,
                       'B4J_canno'        : canno,
                       'B4J_refno'        : refno,
                       'B4J_bench'        : self,
                       'redmine_user'     : self.config [ 'redmine_user' ],
                       'redmine_key'      : self.config [ 'redmine_key' ]}
            if etype == "CMIP5":
                varfile = find_files ( "*%s*%s*%s" % \
                                        ( VARNAME_4_CMIP5_METADATA, name, '.nc' ),
                                        path, raise_e = True )
                kwargs [ 'var_filename' ] = varfile
                e = CMIP5_Experiment ( self, **kwargs )
            elif etype == "MPIESM1":
                e = MPIESM1_Experiment ( self, **kwargs )
            else:
                e = Experiment ( self, **kwargs )
            self.experiments.append ( e )
        if n_refs < 1 or n_cans < 1:
            raise Exception ( "No candidate or no reference experiment defined."  )
        self.n_cans = n_cans
        self.n_refs = n_refs

        self.init_diagnostics ()
        
        self.logger.debug ( "Exiting." )

        
    def init_diagnostics ( self ):
        self.logger.debug ( "Starting." )
        
        # read configspec object
        if not os.path.isfile ( b4j_field_configspec_fn ):
            raise Exception ( "\nERROR no file named %s." % \
                              b4j_field_configspec_fn )
            sys.exit ( 1 )
        try:
            field_configspec = configobj.ConfigObj ( b4j_field_configspec_fn,
                                                     interpolation=False,
                                                     list_values=False,
                                                     _inspec=True)
        except ( configobj.ConfigObjError, IOError), e:
            print 'Could not read configspec object from file %s: %s' % \
                ( b4j_field_configspec_fn, e)
            sys.exit ( 1 )

        self.logger.debug (  "Configuration specialization for field configs read." )        

        f_s = self.config [ "fields" ]
        
        for f_n, f_attr in f_s.iteritems ():
            
            if f_attr and isinstance ( f_attr, configobj.Section ):
                field_config_fn = f_attr [ "config_filename" ]
                field_configs_path = self.get_b4j_path ( field_config_fn, "configs" )
                
                try:
                    field_config = configobj.ConfigObj ( field_configs_path,
                                                         configspec=field_configspec )
                except ( configobj.ConfigObjError, IOError), e:
                    print 'Could not read field config object from file %s: %s' % \
                        ( field_configs_path, e)
                    sys.exit ( 1 )

                val  = Validator()
                test = field_config.validate ( val, preserve_errors=True )
                if test == True:
                    self.logger.debug (  """\nValidation of configuration for diagnostics for field '%s' with files 
'%s' and '%s' successfully run.\n""" %  ( f_n, field_configs_path, b4j_field_configspec_fn ))
                else:
                    print test
                    raise Exception (  "\nERROR validating configuration %s with specialization %s.\n%s" %\
                                       ( field_configs_path, b4j_field_configspec_fn,
                                         ', '.join ( "item '%s'" %  k for k,v in test.items () if v == False)))
                    # TODO: use flattenErrors ( test ) to print exact source of error
                    sys.exit ( 1 )

                # validated field configuration now in field_config
                # initializing diagnostic_fields
                
                field_name = field_config [ 'field_name' ]
                if field_name != f_n:
                    raise Exception ( """
\n\nERROR: By convention the value for the attribute 'field_name' in the field 
configuration file '%s' have to accord with the corresponding section name 
 ( here:= '%s' ) in the b4j main configuration file '%s'.\n""" % \
                                      ( field_configs_path, f_n, self.bench_configs_path ))
                    sys.exit ( 1 )

                    
                self.diagnostic_fields [ f_n ] = field_config                
                self.diagnostic_fields [ f_n ][ 'graphics_path' ]  = self.graphics_path + os.sep + f_n
                self.diagnostic_fields [ f_n ][ 'processed_path' ] = self.processed_path + os.sep + f_n
                for d in [ self.diagnostic_fields [ f_n ][ 'graphics_path' ],
                           self.diagnostic_fields [ f_n ][ 'processed_path' ]]:
                    try:
                        os.makedirs ( d )
                    except:
                        pass
                
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
            
        self.bench_configs_path = os.path.join ( self.configs_path, self.bench_name + ".cfg" )
        with open ( config_template ) as f:
            template_contents = f.read ()
        config_contents = template_contents % { '<BENCH_NAME>': self.bench_name,
                                                '<BENCH_CONFIGSPEC_FNAME>': b4j_configspec_fn,
                                                '<AUTHOR_NAME>': getpass.getuser (),
                                                '<GIT_ORIGIN>': '%s' % get_git_origin (),
                                                '<GIT_BRANCH>': get_git_branch (),
                                                '<GIT_SHA>': get_git_sha (),
                                                '<GIT_STATUS>': get_git_status (),
                                                '<HOSTNAME>': gethostname (),
                                                '<IP_ADDRESS>': get_ip_address (),
                                                '<B4J_BASE>': bench4jsbach_base,
                                                '<TODAY>': str ( datetime.date.today ()),
                                                '<TIMESTAMP>': datetime.datetime.utcnow(),
                                                '<FIELDS_CONFIGSPEC_FAME>' : b4j_field_configspec_fn }
        with open ( self.bench_configs_path, "w" ) as f:
            f.write ( config_contents )

        dkl_logger.debug ( "Exiting." )


    def get_b4j_path ( self, filename, b4j_filetype ):
        assert b4j_filetype and b4j_filetype in bench_dirnames

        dirmap = {  "scripts"   : [  self.scripts_path,
                                     bench4jsbach_base + "/prod/scripts" ],
                    "configs"   : [ self.configs_path,
                                    bench4jsbach_base + "/prod/configs" ],
                    "processed" : [ self.processed_path ], 
                    "graphics"  : [ self.graphics_path ],
                    "report"    : [ self.report_path ],
                    "logging"   : [ self.logging_path ]}
        
        dirpath = dirmap [ b4j_filetype ]

        for dp in dirpath:
            filepath = os.path.join ( dp, filename )
            if os.path.isfile ( filepath ):
                break
            elif dirpath.index ( dp ) == len ( dirpath ) -1:
                print ( "ERROR: No file named '%s' found." % filename )
                sys.exit ( 1 )
            else:
                pass

        return filepath
            
        
    def __repr__ ( self ):
        rep =  """\n\n{} ( bench_name = {}, 
\tbench_path        = {},
\tconfigs_path       = {},
\tprocessed_path         = {}, 
\tgraphics_path      = {}, 
\treport_path       = {}, 
\tbench_configs_path = {}, 
\tconfig            = {}, 
\texperiments       = {})\n\n""".format ( self.__class__.__name__,
                print_none ( self.bench_name        ),
                print_none ( self.bench_path        ),
                print_none ( self.configs_path       ),
                print_none ( self.processed_path     ),
                print_none ( self.graphics_path      ),
                print_none ( self.report_path       ),
                print_none ( self.bench_configs_path ),
                print_item_or_list_rec ( self.config, 28 ),
                print_none ( self.experiments ))
        
        rep += Diagnostics.__repr__ ( self )
        rep += Report.__repr__ ( self )
        
        return rep
        
    

