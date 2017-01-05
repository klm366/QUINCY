#! /usr/bin/env python
""" BENCH4JSBACH benchmarking run clean script
USAGE:
       >>>./%s <B4JCONFIG>
WITH <B4JCONFIG> being the configuration file of the b4j run in focus."""

import os,sys
from optparse import OptionParser


script_dir = os.path.dirname ( os.path.realpath ( __file__ ))
bench4jsbach_base = os.path.dirname ( os.path.realpath ( script_dir ))
bin_dir = os.path.join ( bench4jsbach_base, "bin" )
if os.path.isdir ( bin_dir ):
    sys.path.append ( bin_dir )
else:
    raise Exception ( "Directory %s with source files not found." % bin_dir )
os.environ [ 'BENCH4JSBACH_BASE' ] = bench4jsbach_base
import bench
from util import print_item_or_list_rec

work_dir = os.path.join ( bench4jsbach_base, "work" )

parser = OptionParser()
(options, args) = parser.parse_args()

nargs = len ( args )
if nargs != 1:
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Usage:" + __doc__
    print "******************************************"
elif not os.path.isfile ( args [ 0 ]):
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* No configuration file found: '%s'." % \
          args [ 0 ]
    print "******************************************"
else:
    print "******************************************"
    print "* BENCH4JSBACH: "
    print "* Found file '%s'." %  args [ 0 ]
    print "******************************************"

    try:
        b4j_obj = bench.Bench_Object ( args [ 0 ])
    except:
        print "******************************************"
        print "* BENCH4JSBACH: "
        print "* File '%s' contains no valid b4j configuration." %  args [ 0 ]
        print "* Exit on error. "
        print "******************************************"
        sys.exit ( 1 )

    dirs2clean = [ b4j_obj.graphics_path,
                   b4j_obj.processed_path,
                   b4j_obj.report_path ]

    count = 0
    for dir in dirs2clean:
        for root, subdir, files in os.walk ( dir ):
            for f in files:
                os.remove ( os.path.join ( root, f ))
                count += 1

    print "\n******************************************"
    print "* BENCH4JSBACH: "
    print "*"
    print "* >> %s %s" % ( __file__, args [ 0 ])
    print "*"
    print "* Finished successfully."
    print "* Deleted %d files out of the directories:" % count
    print print_item_or_list_rec ( dirs2clean, 3 )
    print "******************************************"
