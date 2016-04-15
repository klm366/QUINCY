#! /usr/bin/env python

import os, subprocess, shutil,getpass,datetime
import copy
import configobj
from validate import Validator
from sets import Set


from experiment import Experiment, CMIP5_Experiment, MPIESM1_Experiment
from util import logger as dkl_logger, setup_logger, print_none, log_message_newline, str2bool
from analysis import Analysis
from report import Report
from report import ( analyses_II_analyses,
                     analyses_III_dicts,
                     analyses_IV_cans,
                     analyses_V_benches )


REF_EXP_IDX              = 0
SWITCHED_ON_IDX          = TYPE_IDX = 0
TEMP_SCALE_IDX           = PATH_IDX = 1
ISREF_IDX                = 2
ISCAN_IDX                = 3
VARNAME_4_CMIP5_METADATA = 'rsds'

bench4jsbach_base        =  os.environ [ 'BENCH4JSBACH_BASE' ]
if not bench4jsbach_base: raise Exception ( "Environment variable not set." )
bench4jsbach_configuration =  os.path.join ( bench4jsbach_base, "configuration" )
config_template            = os.path.join ( bench4jsbach_configuration, "template.cfg" )
b4j_configspec_fn          = os.path.join ( bench4jsbach_configuration, "b4j_config.spec" )
b4j_field_configspec_fn    = os.path.join ( bench4jsbach_configuration, "b4j_field_config.spec" )
work_base                  = os.path.join ( bench4jsbach_base, "work" )
bench_dirnames             = [ "config", "data", "graphic", "report", "logging" ]
    


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
    n_cans            = 0
    n_refs            = 0
    logger            = None
    logger_name       = 'B4J_LOGGER'
    reporting         = True
    land_sea_mask     = None
    
    
    def __init__ ( self, config_arg, **kwargs ):
        super ( Bench_Object, self ).__init__ ()
        
        dkl_logger.debug ( "Starting." )
        
        if not os.path.exists ( config_arg ):
            self._writeconfig__init__ ( config_arg )
            dkl_logger.debug ( "Exiting." )
            return
        bench_config_path = config_arg
 
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
            config = configobj.ConfigObj ( bench_config_path,
                                           configspec=configspec )
        except ( configobj.ConfigObjError, IOError), e:
            print '\nCould not read "%s": %s' % ( bench_config_path, e)
            
        dkl_logger.debug ( "Read configuration file %s." % ( bench_config_path ))

        val  = Validator()
        test = config.validate ( val )
        if test == False:
            raise Exception (  "\nERROR validating configuration %s with specialization %s.\n" %\
                               ( bench_config_path, b4j_configspec_fn ))
        else:
            dkl_logger.debug ( "\nValidation of configuration from %s with specification from %s successfully run." % \
                ( bench_config_path, b4j_configspec_fn ))

            
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

        self.land_sea_mask = config [ "land_sea_mask" ]
        if not os.path.isfile ( self.land_sea_mask ):
            raise Exception ( "Error: Path for land sea mask '%s' does not exists." % \
                              self.land_sea_mask )

        es = self.config [ "experiments" ]
        from analysis import find_files
        n_refs = n_cans = 0
        for en, attr in es.iteritems ():
            name = en
            etype = attr.get ( "type", "CMIP5" )
            path = attr [ "path" ]
            iscan = attr.get ( 'iscan', True )
            isref = attr.get ( 'isref', False )
            n_cans += 1 if iscan else 0
            n_refs += 1 if isref else 0
            canno = n_cans - 1 if iscan else None
            refno = n_refs - 1 if isref else None
            kwargs = { 'B4J_name'  : name,
                       'B4J_etype' : etype,
                       'B4J_path'  : path,
                       'B4J_iscan' : iscan,
                       'B4J_isref' : isref,
                       'B4J_canno' : canno,
                       'B4J_refno' : refno,
                       'B4J_bench' : self }
            if etype == "CMIP5":
                varfile = find_files ( "*%s*%s*%s" % \
                                        ( VARNAME_4_CMIP5_METADATA, name, '.nc' ),
                                        path, raise_e = True )
                kwargs [ 'var_filename' ] = varfile
                e = CMIP5_Experiment ( **kwargs )
            elif etype == "MPIESM1":
                e = MPIESM1_Experiment ( **kwargs )
            else:
                e = Experiment ( **kwargs )
            self.experiments.append ( e )
        if n_refs < 1 or n_cans < 1:
            raise Exception ( "No candidate or no reference experiment defined."  )
        self.n_cans = n_cans
        self.n_refs = n_refs

        self.init_analyses ()
        
        self.logger.debug ( "Exiting." )


        
    def init_analyses ( self ):
        self.logger.debug ( "Starting." )
        
        print self.analyses_fields
        
        # read configspec object
        if not os.path.isfile ( b4j_field_configspec_fn ):
            raise Exception ( "\nERROR no file named %s." % \
                              b4j_field_configspec_fn )
        try:
            field_configspec = configobj.ConfigObj ( b4j_field_configspec_fn,
                                                     interpolation=False,
                                                     list_values=False,
                                                     _inspec=True)
        except ( configobj.ConfigObjError, IOError), e:
            print 'Could not read configspec object from file %s: %s' % \
                ( b4j_field_configspec_fn, e)



        print "\nConfiguration specialization for analyses configs read.\n"


        f_s = self.config [ "fields" ]
        for f_n, f_attr in f_s.iteritems ():
            if f_attr and isinstance ( f_attr, configobj.Section ):
                print "in loop for field %s: " % f_n
                field_config_fn = f_attr [ "config_fname" ]
                field_config_path = os.path.join ( self.config_path, field_config_fn )
                if not os.path.isfile ( field_config_path ):
                    raise Exception ( "\nERROR no file named %s." % field_config_path )
                try:
                    field_config = configobj.ConfigObj ( field_config_path,
                                                         configspec=field_configspec )
                except ( configobj.ConfigObjError, IOError), e:
                    print 'Could not read field config object from file %s: %s' % \
                        ( field_config_path, e)
        
                val  = Validator()
                test = field_config.validate ( val )
                if test == False:
                    raise Exception (  "\nERROR validating configuration %s with specialization %s.\n" %\
                                       ( field_config_path, b4j_field_configspec_fn ))
                else:
                    print """\nValidation of configuration for analyses of field '%s' with files 
'%s' and '%s' successfully run.\n""" %  ( f_n, field_config_path, b4j_field_configspec_fn )

                    
                # validated field configuration now in field_config
                # initializing analyses_fields

                field_name = field_config [ 'field_name' ]
                if field_name != f_n:
                    raise Exception ( """
\n\nERROR: By convention the value for the attribute 'field_name' in the field 
configuration file '%s' have to accord with the corresponding section name 
 ( here:= '%s' ) in the b4j main configuration file '%s'.\n""" % \
                                      ( field_config_path, f_n, self.bench_config_path ))

                self.analyses_fields [ f_n ] = copy.deepcopy ( analyses_II_analyses )
                self.analyses_fields [ f_n ][ 'frequency' ] = field_config [ 'frequency' ]
                self.analyses_fields [ f_n ][ 'teXts' ] = field_config [ 'teXts' ]
                a_s = field_config [ 'analyses' ]
                for a_n, a_attr in a_s.iteritems ():
                    if a_attr and isinstance ( a_attr, configobj.Section ):
                        print "in loop for analyse %s: " % a_n
                
                        # frequency, teXts, analyses
                        a_dict = copy.deepcopy ( analyses_III_dicts )
                        a_dict [ 'name' ]     = a_attr [ "name" ]
                        a_dict [ 'title' ]    = a_attr [ "title" ]
                        a_dict [ 'unit' ]     = a_attr [ "unit" ]
                        a_dict [ 'type' ]     = a_attr [ "type" ]
                        a_dict [ 'isbench' ]  = a_attr [ "isbench" ]
                        a_dict [ 'teXts' ]    = a_attr [ "teXts" ]
                        a_dict [ 'cans' ]     = { '%d' % canno : copy.deepcopy ( analyses_IV_cans ) \
                                                 for canno in range ( self.n_cans )}
                        self.analyses_fields [ f_n ][ a_dict [ 'name' ]] = a_dict

                
            
        print self.analyses_fields
        
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
                                                '<BENCH_CONFIGSPEC_FNAME>': b4j_configspec_fn,
                                                '<AUTHOR_NAME>': getpass.getuser (),
                                                '<TODAY>': str ( datetime.date.today ()),
                                                '<FIELDS_CONFIGSPEC_FAME>' : b4j_field_configspec_fn }
        with open ( self.bench_config_path, "w" ) as f:
            f.write ( config_contents )

        dkl_logger.debug ( "Exiting." )

        
    def __repr__ ( self ):
        repr =  """\n\n{} ( bench_name = {}, 
\tbench_path        = {},
\tconfig_path       = {},
\tdata_path         = {}, 
\tgraphic_path      = {}, 
\treport_path       = {}, 
\tbench_config_path = {}, 
\tconfig            = {}, 
\texperiments       = {}\n\n""".format ( self.__class__.__name__,
                print_none ( self.bench_name        ),
                print_none ( self.bench_path        ),
                print_none ( self.config_path       ),
                print_none ( self.data_path         ),
                print_none ( self.graphic_path      ),
                print_none ( self.report_path       ),
                print_none ( self.bench_config_path ),
                print_none ( self.config            ),
                print_none ( self.experiments ))#.replace ( "\n", "" )

        print repr
    
        super ( Bench_Object, self ).__repr__ ()

        return ""
        
    

