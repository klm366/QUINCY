#! /usr/bin/env python
""" BENCH4JSBACH benchmarking script """

import os,sys,shutil,timeit

### Set environment
script_dir = os.path.dirname ( os.path.realpath ( __file__ ))
bench4jsbach_base = os.path.dirname ( os.path.realpath ( script_dir ))
bin_dir = os.path.join ( bench4jsbach_base, "bin" )
if os.path.isdir ( bin_dir ):
    sys.path.append ( bin_dir )
else:
    raise Exception ( "Directory %s with source files not found." % bin_dir )
os.environ [ 'BENCH4JSBACH_BASE' ] = bench4jsbach_base
import bench


start = timeit.default_timer ()

### Start action:
### Three cases ( say hello / create configuration file and environment / run benchmarking from configuration file ) 
nargs = len ( sys.argv )
if nargs == 1:
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Usage:" + __doc__
    print "******************************************"
elif not os.path.isfile ( sys.argv [ 1 ]):
    # todo: if bench-environment exists exit on error, overwrite only with extra force configuration
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Creating environment and configfile for new bench4jsbach object '%s'." % \
          sys.argv [ 1 ]
    print "******************************************"

    b4j_obj = bench.Bench_Object ( sys.argv [ 1 ])
        
    print "\n******************************************"
    print "* BENCH4JSBACH: "
    print "* Created environment \n*\t" +\
          "'%s' \n* and the configuration file \n*\t" % b4j_obj.bench_path +\
          "'%s'.\n* Please edit this file." % b4j_obj.bench_configs_path
    print "******************************************"
else:
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Reading configuration file '%s'." %  sys.argv [ 1 ]
    print "******************************************"

    b4j_obj = bench.Bench_Object ( sys.argv [ 1 ])
    
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Initialized Bench_Object. Defining Diagnostics."
    print "* Reading configuration file '%s'." %  sys.argv [ 1 ]
    print "******************************************"

    # fetch outdata
    for exp in b4j_obj.experiments:
        exp.fetch_outdata ()

    # diagnosis
    for f_n,f_bench in b4j_obj.config [ 'fields' ].iteritems():
        ret = b4j_obj.diagnostic_generic ( f_n )
        b4j_obj.logger.debug ( "Finished diagnostic for '%s' successfully." % f_n ) if not ret else \
                b4j_obj.logger.debug ( "Diagnostic for '%s' failed." % f_n )

    # reporting
    b4j_obj.write_report () if b4j_obj.reporting else None
    
    b4j_obj.logger.debug ( "Exiting." )

    stop = timeit.default_timer ()

    print "\n******************************************"
    print "* BENCH4JSBACH: "
    print "* Finished successfully. Elapsed time %.4f seconds." % ( stop - start )
    if b4j_obj.reporting: print "* Find the report %s in %s." % ( b4j_obj.bench_name + '.pdf',
                                                                  b4j_obj.report_path )
    print "******************************************"
