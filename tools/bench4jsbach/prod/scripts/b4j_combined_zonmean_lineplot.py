#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
from cdo import Cdo

print "\n\nStarting %s" % os.path.realpath(__file__)

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
from b4jplot import zon_plot_wrapper

#####
##### check input parameter
#####

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':', 1 ) for token in sys.argv [ 1: ])
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )

#####
##### plotting
#####

can_name          = kwargs [ "can_name" ]
canno             = kwargs [ "canno" ]

title_tmpl        = kwargs [ "title_tmpl" ]
legend_tmpl       = kwargs [ "legend_tmpl" ]

field_name        = kwargs [ "field_name" ]

varname           = kwargs [ "varname" ]
unit              = kwargs [ "unit" ]
funstr            = kwargs [ 'funstr' ]
statstr           = kwargs [ 'statstr' ]
statstr_tm        = kwargs.get ( 'statstr_tm', None )

inpathlist_gm_can = [ kwargs [  "input_can" ]]
varname_can       = kwargs [ "varname_can" ]
unit_can          = kwargs [ "unit_can" ]
funstr_can        = kwargs [ 'funstr_can' ]

outpath           = kwargs [ "outpath" ]

subtitle_can      = title_tmpl % { 'funstr'  : funstr_can,
                                   'stat1str' : statstr,
                                   'stat2str' : statstr_tm,
                                   'expabr'   : "C_%s" % ( canno ),
                                   'expl'     : ", $C_{%s}$: %s" % ( canno, can_name )}

subtitle_delta    = title_tmpl % { 'funstr'  : funstr,
                                   'stat1str': statstr,
                                   'stat2str': statstr_tm,
                                   'expabr'  : "C_%sR" % ( canno ),
                                   'expl'    : ", $C_{%s}$: %s" % ( canno, can_name )}

token             = "input_ref"
keylist           = [ k for k in kwargs.keys () if k.startswith ( token )]
inpathlist_gm     = []
legendlist_gm     = []
for k in keylist:
    refno  = k.split ( token, 1 )[1][0]
    inpathlist_gm.append ( kwargs [ k ])
    legend =  legend_tmpl % { 'funstr'  : funstr,
                              'stat1str': statstr,
                              'stat2str': statstr_tm,
                              'expabr'  : "C_%sR_%s" % ( canno, refno ),
                              'expl'    : ",\n $R_{%s}$: %s" % \
                              ( refno,  kwargs [ "ref%s"  % refno ])}
    legendlist_gm.append ( legend )

print inpathlist_gm
print legendlist_gm

zon_plot_wrapper ( varname_can       = varname_can,
                   inpathlist_can    = inpathlist_gm_can, 
                   unit_can          = unit_can,
                   subtitle_can      = subtitle_can,
                   varname           = varname,
                   inpathlist        = inpathlist_gm, 
                   unit              = unit,
                   subtitle          = subtitle_delta,
                   outpath           = outpath, 
                   legendlist        = legendlist_gm )


#####
##### finishing
#####

print "Plotting for field %s with '%s' run successfully." % ( field_name, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
