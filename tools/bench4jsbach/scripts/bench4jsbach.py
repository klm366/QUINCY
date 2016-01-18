#! /usr/bin/env python
""" BENCH4JSBACH benchmarking script """

import os,sys,shutil

script_dir = os.path.dirname ( os.path.realpath ( __file__ ))
bench4jsbach_base = os.path.dirname ( os.path.realpath ( script_dir ))
bin_dir = os.path.join ( bench4jsbach_base, "bin" )
if os.path.isdir ( bin_dir ):
    sys.path.append ( bin_dir )
else:
    raise Exception ( "Directory %s with source files not found." % bin_dir )
import sys
os.environ [ 'BENCH4JSBACH_BASE' ] = bench4jsbach_base

#from util import logger as dkl_logger
import bench

#dkl_logger.debug ( "hello" )
nargs = len ( sys.argv )
if nargs == 1:
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Usage:" + __doc__
    print "******************************************"
elif not os.path.isfile ( sys.argv [ 1 ]):
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Creating environment and configfile for new bench4jsbach object '%s'." % \
          sys.argv [ 1 ]
    print "******************************************"

    b4j_obj = bench.Bench_Object ( sys.argv [ 1 ])
    b4j_obj.logger.debug ( b4j_obj )
        
    print "\n******************************************"
    print "* BENCH4JSBACH: "
    print "* Created environment \n*\t" +\
          "'%s' \n* and the configuration file \n*\t" % b4j_obj.bench_path +\
          "'%s'.\n* Please edit this file." % b4j_obj.bench_config_path
    print "******************************************"
else:
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Reading configuration file '%s'." %  sys.argv [ 1 ]
    print "******************************************"

    b4j_obj = bench.Bench_Object ( sys.argv [ 1 ])

    # analysis
    for var,var_bench in b4j_obj.config [ 'variables' ].iteritems():

        if not var_bench [ bench.SWITCHED_ON_IDX ]: continue

        for fun_suffix in [ "_analyse", "_section" ]:
            fun_name = var + fun_suffix
            try:
                fun = getattr ( b4j_obj, fun_name )
                b4j_obj.logger.debug ( "Loaded function '%s'." % fun_name )
            except Exception, e:
                raise AttributeError ( "%s. \n Bench object has no attribute '%s'" % 
                                       ( str ( e ), fun_name ))

            ret = fun ()
        
        b4j_obj.logger.debug ( "Finished analysis of '%s' successfully." % var ) if not ret else \
                b4j_obj.logger.debug ( "Analysis of '%s' failed." % var )
    
    b4j_obj.write_report ()
    
    b4j_obj.logger.debug ( "Exiting." )

    print "\n******************************************"
    print "* BENCH4JSBACH: "
    print "* Finished successfully."
    print "******************************************"
