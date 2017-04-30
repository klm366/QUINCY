#! /usr/bin/env python
# -*- coding:  utf-8 -*-y

import sys, os
from cdo import Cdo

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
from b4jplot import mapplot, mapplot_seasonal

#####
##### check input parameter
#####

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

kwargs        =  dict ( token.split ( ':', 1 ) for token in sys.argv [ 1: ])
print "subprocess '%s', kwargs = " % __file__
print print_item_or_list_rec ( kwargs, 3 )

fieldname      = kwargs [ "fieldname" ]

#####
##### plotting
#####

can_name          = kwargs.get ( "can_name" )
canno             = kwargs.get ( "canno" )

title_tmpl        = kwargs.get ( "title_tmpl" )
funstr            = kwargs.get ( 'funstr', "FUN" )
statstr           = kwargs.get ( 'statstr', 'STAT' )

inpath            = kwargs.get ( 'inpath' )
field_name        = kwargs.get ( "field_name" )
varname           = kwargs.get ( "varname" )
unit              = kwargs.get ( "unit", 'UNKN' )

colormap          = kwargs.get ( 'colormap', 'BrBG' )
midpoint          = kwargs.get ( 'midpoint', None )
vmin              = kwargs.get ( 'vmin', None )
vmax              = kwargs.get ( 'vmax', None )

outpath           = kwargs.get ( "outpath" )

nref              = kwargs.get ( 'nref', None )

if nref == "single":
    refno         = kwargs.get ( 'refno' )
    ref_name      = kwargs.get ( 'ref_name' )
    title         = title_tmpl % { 'funstr': funstr,
                                   'stat1str': statstr,
                                   'expabr': "C_%sR_%s" % ( canno, refno ),
                                   'expl': ", $C_{%s}$: %s, $R_{%s}$: %s" % (
                                       canno, can_name, refno, ref_name )}
else:
    title         = title_tmpl % { 'funstr': funstr,
                                   'stat1str': statstr,
                                   'expabr': "C_%s" % ( canno ),
                                   'expl': ", $C_{%s}$: %s" % ( canno, can_name )}

mapplot_seasonal ( varname    = varname,
                   inpath     = inpath,
                   outpath    = outpath,
                   unit       = unit,
                   title      = title,
                   midpoint   = midpoint,
                   vmin       = vmin,
                   vmax       = vmax,
                   colormap   = colormap )

#####
##### finishing
#####

print "Plotting for field %s with '%s' run successfully." % ( fieldname, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
