#! /usr/bin/env python
# -*- coding:  utf-8 -*-

import sys, os
import numpy as np
from netCDF4 import Dataset


def get_value ( filepath, varname ):
    nc_fid = Dataset ( filepath )
    lat    = nc_fid.variables [ 'lat' ][ : ]
    vardat = nc_fid.variables [ varname ] [ : ]
    nc_fid.close ()
    print vardat  [ :, :, 0 ]
    return vardat [ :, :, 0 ]
    

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
from b4jplot import hist_plot


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

vartoken = [ 'E', 'P', 'RD' ]
varnames = { 'E': 'evapotranspiration',
             'P': 'precip_acc',
             'RD': 'disch' }
filenametemplate = "%(varacronym)s_ysdfm_%(refacronym)s" 
reftoken = "ref"
nchar = 4
reflist    =  [ k for k in kwargs.keys () if k.startswith ( reftoken ) and len ( k ) == nchar ]
print reflist

# refmap = { k : np.array ([[ get_value ( kwargs [ filenametemplate % { 'varacronym' : 'E',
#                                                                       'refacronym' : k }],
#                                         varnames [ 'E' ])]]) \
#            for k in reflist }
# print refmap
# print print_item_or_list_rec ( refmap, 3 )
# for k,v in refmap:
#     type ( v )
#     print v

acca = np.array ([[]])
for k in reflist:
    accb = np.array ([[]])
    for v in vartoken:
        b = get_value ( kwargs [ filenametemplate % { 'varacronym' : v,
                                                      'refacronym' : k }],
                           varnames [ v ])
        print b
        accb = np.concatenate (( accb, b ), axis = 1 ) 
        print accb
        
    accb = np.concatenate (( accb, np.array ([[ accb [ 0 ][ 1 ] -\
                                                accb [ 0 ][ 0 ] -\
                                                accb [ 0 ][ 2 ]]])), axis = 1 )
    print "accb"
    print accb
    print accb.shape
    print "acca"
    print acca
    print acca.shape
    acca = np.concatenate (( acca, accb ), axis = 0 ) if acca != [[]] else accb
    print acca

hist_plot ( acca, "mm/year", [ "schnick", "schnack", "schuck" ], "hallo" )

#####
##### exiting
#####
print "Exiting %s" % os.path.realpath(__file__)
