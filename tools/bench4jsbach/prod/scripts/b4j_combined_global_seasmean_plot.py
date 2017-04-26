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
from b4jplot import seas_timemean_plot_wrapper

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

field_name  = kwargs [ "field_name" ]
varname     = kwargs [ "varname" ]
varname_can = field_name
unit_can    = kwargs [ "unit_can" ]
unit        = kwargs [ "unit" ]
can_name       = kwargs [ "can_name" ]
canno          = kwargs [ "canno" ]
outpath = kwargs [ "albedo_wdsgtm_plot" ]
inpathlist_gm_can = [ kwargs [  "A_seasmean_gm" ]]
inpathlist_tm_can =  [ kwargs [  "A_seasmean_gtm" ]]
subtitle_can   = kwargs [ "plot_title_albedo_wdsgtm_plot" ] % { 'funstr' : "A",
                                                                'stat1str': "gsm",
                                                                'stat2str': "gtm",
                                                                'expabr'  : "C_%s" % ( canno ),
                                                                'expl'    : ", $C_{%s}$: %s" % ( canno, can_name )}
subtitle_delta = kwargs [ "plot_title_albedo_wdsgtm_plot" ] % { 'funstr' : "w\Delta A",
                                                                'stat1str': "gsm",
                                                                'stat2str': "gtm",
                                                                'expabr'  : "C_%sR" % ( canno ),
                                                                'expl'    : ", $C_{%s}$: %s" % ( canno, can_name )}



                          
token             = "albedo_wdsgm_ref"
keylist           = [ k for k in kwargs.keys () if k.startswith ( token )]
inpathlist_gm     = []
legendlist_gm     = []
for k in keylist:
    refno = k.split ( token, 1 )[1][0]
    inpathlist_gm.append ( kwargs [ k ])
    legend =  kwargs [ "legend_tmpl_albedo_wdsgtm_plot" ] % { 'funstr' : "w \Delta A ",
                                                              'stat1str': "gsm",
                                                              'expabr': "C_%sR_%s" % ( canno, refno ),
                                                              'expl'     : " $R_{%s}$: %s" % \
                                                              ( refno,  kwargs [ "ref%s"  % refno ])}
    legendlist_gm.append ( legend )


token             = "albedo_wdsgtm_ref"
keylist           = [ k for k in kwargs.keys () if k.startswith ( token )]
inpathlist_tm     = []
legendlist_tm     = []
for k in keylist:
    refno = k.split ( token, 1 )[1][0]
    inpathlist_tm.append ( kwargs [ k ])
    legend =  kwargs [ "legend_tmpl_albedo_wdsgtm_plot" ] % { 'funstr'  : "w \Delta A ",
                                                              'stat1str': "gtm",
                                                              'expabr'  : "C_%sR_%s" % ( canno, refno ),
                                                              'expl'    : "" }
    legendlist_tm.append ( legend )

seas_timemean_plot_wrapper ( varname = varname,
                            varname_can = "albedo",
                            inpathlist_gm_can = inpathlist_gm_can, 
                            inpathlist_tm_can = inpathlist_tm_can, 
                            subtitle_can = subtitle_can,
                            inpathlist_gm = inpathlist_gm, 
                            inpathlist_tm = inpathlist_tm,
                            outpath = outpath, 
                            unit = unit,
                            unit_can = unit_can,
                            subtitle_delta = subtitle_delta,
                            legendlist_gm = legendlist_gm, 
                            legendlist_tm = legendlist_tm )


#####
##### finishing
#####

print "Plotting for field %s with '%s' run successfully." % ( field_name, __file__ )

print "Exiting %s" % os.path.realpath(__file__)
