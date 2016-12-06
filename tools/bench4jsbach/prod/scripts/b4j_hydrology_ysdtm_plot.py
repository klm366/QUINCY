#! /usr/bin/env python
# -*- coding:  utf-8 -*-

import sys, os

print "Starting %s" % os.path.realpath(__file__)


#####
##### bench4jsbach library to environment path
#####
try:
    b4j_base = os.environ [ 'BENCH4JSBACH_BASE' ]
except KeyError: 
    print "The environment variable 'BENCH4JSBACH_BASE' is not set."
    sys.exit(1)

bin_dir = os.path.join ( b4j_base, 'bin' )
if os.path.isdir ( bin_dir ):
    sys.path.append ( bin_dir )
else:
    raise Exception ( "Directory %s with source files not found." % bin_dir )
from util import print_item_or_list_rec
from b4jplot import zon_plot


#####
##### check input parameter
#####
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.strip ().split ( ':', 1 ) for token in sys.argv [ 1: ])
print kwargs.keys ()
print kwargs
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )


#####
##### argument definition and plot method calls
#####
varname    = "evapotranspiration"
unit       = "mm/year"
title      = "$ysdzm_{C_{%s}R}$ and $ysdfm_{C_{%s}R}$" % \
             ( kwargs [ "canno" ], kwargs [ "canno" ])

token      = "E_ysd-zonmean_ref"
inpathlist = []
legendlist = []
keylist    =  [ k for k in kwargs.keys () if k.startswith ( token )]
for k in keylist:
    refno = k.split ( token, 1 )[1][0]
    inpathlist.append ( kwargs [ k ])
    legend = "$ysdzm_{C_{%s}R_{%s}}, ysdzm_{C_{%s}R_{%s}}$\n$R_{%s}$: %s" % \
             ( kwargs [ "canno" ], refno, kwargs [ "canno" ], refno,
               refno, kwargs [ "ref%s"  % refno ])
    legendlist.append ( legend )

zon_plot ( title, varname, inpathlist, unit, legendlist, kwargs [ "E_ysd-zonmean_plot" ])


varname    = "precip_acc"
unit       = "mm/day"
title      = "$ysdzm_{C_{%s}R}$ and $ysdfm_{C_{%s}R}$" % \
             ( kwargs [ "canno" ], kwargs [ "canno" ])

token      = "P_ysd-zonmean_ref"
inpathlist = []
legendlist = []
keylist    =  [ k for k in kwargs.keys () if k.startswith ( token )]
for k in keylist:
    refno = k.split ( token, 1 )[1][0]
    inpathlist.append ( kwargs [ k ])
    legend = "$ysdzm_{C_{%s}R_{%s}}, ysdtm_{C_{%s}R_{%s}}$\n$R_{%s}$: %s" % \
             ( kwargs [ "canno" ], refno, kwargs [ "canno" ], refno,
               refno, kwargs [ "ref%s"  % refno ])
    legendlist.append ( legend )

zon_plot ( title, varname, inpathlist, unit, legendlist, kwargs [ "P_ysd-zonmean_plot" ])


#####
##### exiting
#####
print "Exiting %s" % os.path.realpath(__file__)
