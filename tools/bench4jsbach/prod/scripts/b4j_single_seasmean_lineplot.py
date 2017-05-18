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
from b4jplot import seas_timemean_plot

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
funstr            = kwargs [ 'funstr' ]
statstr           = kwargs [ 'statstr' ]
statstr_tm        = kwargs [ 'statstr_tm' ]

outpath           = kwargs [ "outpath" ]
field_name        = kwargs [ "field_name" ]
varname           = kwargs [ "varname" ]
unit              = kwargs [ "unit" ]

title             = title_tmpl % { 'funstr'  : funstr,
                                   'stat1str': statstr,
                                   'stat2str': statstr_tm,
                                   'expabr'  : "C_%sR" % ( canno ),
                                   'expl'    : ", $C_{%s}$: %s" % ( canno, can_name )}

token             = "input_ref"
keylist           = [ k for k in kwargs.keys () if k.startswith ( token )]
keylist.sort () if keylist else None
inpathlist_gm     = []
legendlist_gm     = []
for k in keylist:
    refno  = k.split ( token, 1 )[1][0]
    inpathlist_gm.append ( kwargs [ k ])
    legend =  legend_tmpl % { 'funstr'  : funstr,
                              'stat1str': statstr,
                              'expabr'  : "C_%sR_%s" % ( canno, refno ),
                              'expl'    : " $R_{%s}$: %s" % \
                              ( refno,  kwargs [ "ref%s"  % refno ])}
    legendlist_gm.append ( legend )


token             = "input_tm_ref"
keylist           = [ k for k in kwargs.keys () if k.startswith ( token )]
keylist.sort () if keylist else None
inpathlist_tm     = []
legendlist_tm     = []
for k in keylist:
    refno  = k.split ( token, 1 )[1][0]
    inpathlist_tm.append ( kwargs [ k ])
    legend = legend_tmpl % { 'funstr'  : funstr,
                             'stat1str': statstr_tm,
                             'expabr'  : "C_%sR_%s" % ( canno, refno ),
                             'expl'    : "" }
    legendlist_tm.append ( legend )

seas_timemean_plot ( varname           = varname,
                     inpathlist_gm     = inpathlist_gm, 
                     inpathlist_tm     = inpathlist_tm,
                     outpath           = outpath, 
                     unit              = unit,
                     title             = title,
                     legendlist_gm     = legendlist_gm, 
                     legendlist_tm     = legendlist_tm )

#####
##### finishing
#####

print "Plotting for field %s with '%s' run successfully." % ( field_name, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
