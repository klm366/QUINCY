#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
import ast

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
from b4jplot import hist_plot_wrapper, colors_hydrology

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

can_name          = kwargs.get ( "can_name" )
canno             = kwargs.get ( "canno" )

title_tmpl        = kwargs.get ( "title_tmpl" )

field_name        = kwargs.get ( "field_name" )

varnamelist       = ast.literal_eval ( kwargs.get ( "varnamelist" ))
legendlist        = ast.literal_eval ( kwargs.get ( "legendlist" ))
unit              = kwargs.get ( "unit" )
funstr            = kwargs.get ( 'funstr' )
statstr           = kwargs.get ( 'statstr' )

inpathlist_can    = [ kwargs.get (  "input_can" )]
funstr_can        = kwargs.get ( "funstr_can" )

colorlist         = ast.literal_eval ( kwargs.get ( "colorlist" ))

outpath           = kwargs.get ( "outpath" )

nbar              = len ( varnamelist )

subtitle_can      = title_tmpl % { 'funstr'  : funstr_can,
                                   'stat1str' : statstr,
                                   'expabr'   : "C_%s" % ( canno ),
                                   'expl'     : ", field values"}

subtitle_delta    = title_tmpl % { 'funstr'  : funstr,
                                   'stat1str': statstr,
                                   'expabr'  : "C_%sR" % ( canno ),
                                   'expl'    : ", differences" }

token             = "input_ref"
keylist           = [ k for k in kwargs.keys () if k.startswith ( token )]
keylist.sort () if keylist else None
inpathlist        = []
xticklabellist    = []
for k in keylist:
    refno  = k.split ( token, 1 )[1][0]
    inpathlist.append ( kwargs.get ( k ))
    xticklabel = "$R_{%s}$: %s" % ( refno,  kwargs.get ( "ref%s"  % refno ))
    xticklabellist.append ( xticklabel )
    

hist_plot_wrapper ( nbar           = nbar,
                    colorlist      = colorlist,
                    legendlist     = legendlist,
                    varnamelist    = varnamelist,
                    unit           = unit,
                    can_name       = can_name,
                    canno          = canno,
                    input_can      = inpathlist_can,
                    subtitle_can   = subtitle_can,
                    input_delta    = inpathlist,
                    subtitle_delta = subtitle_delta,
                    xticklabellist = xticklabellist,
                    outpath        = outpath )


#####
##### finishing
#####

print "Plotting for field %s with '%s' run successfully." % ( field_name, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
