#! /usr/bin/env python
""" BENCH4JSBACH benchmarking script """

import os,sys,shutil
sys.path.append ( "/scratch/mpi/mpiles/m300050/bench4jsbach/bin" )
import sys
from util import logger as dkl_logger
import bench

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
    dkl_logger.debug ( b4j_obj )
        
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

    _old_stdout = sys.stdout
    log_file = open ( b4j_obj.logfile_path, "w")
    sys.stdout = log_file

    # analysis
    for var,var_bench in b4j_obj.config [ 'variables' ].iteritems():

        if not var_bench [ bench.SWITCHED_ON_IDX ]: continue

        for fun_suffix in [ "_analyse", "_section" ]:
            fun_name = var + fun_suffix
            try:
                fun = getattr ( b4j_obj, fun_name )
                dkl_logger.debug ( "Loaded function '%s'." % fun_name )
            except Exception, e:
                raise AttributeError ( "%s. \n Bench object has no attribute '%s'" % 
                                       ( str ( e ), fun_name ))

            ret = fun ()
        
        dkl_logger.debug ( "Finished analysis of '%s' successfully." % var ) if not ret else \
                dkl_logger.debug ( "Analysis of '%s' failed." % var )
    
    b4j_obj.write_report ()
    
    dkl_logger.debug ( "Exiting." )
    sys.stdout = _old_stdout
    log_file.close ()

    print "\n******************************************"
    print "* BENCH4JSBACH: "
    print "* Finished successfully."
    print "******************************************"
