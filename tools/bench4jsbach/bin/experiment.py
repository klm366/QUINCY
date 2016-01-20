#! /usr/bin/env python

import os
from netCDF4 import Dataset 

from util import print_none


class Experiment ( object ):
    id    = None
    name  = None
    etype = None
    path  = None
    
    
    def __init__ ( self, **kwargs ):
        self.name     = kwargs.get ( "name" )
        self.etype    = kwargs.get ( "etype", None )
        self.path     = kwargs.get ( "path", None )
        self.buffer   = []
        

    def __repr__ ( self ):
        return "{} ( id = {}, name = {}, etype = {}, path = {}, buffer = {})". \
               format ( self.__class__.__name__,
                        print_none ( self.id ),
                        print_none ( self.name ),
                        print_none ( self.etype ),
                        print_none ( self.path ),
                        print_none ( self.buffer ))


class CMIP5_Experiment ( Experiment ):
    
    def __init__ ( self, **kwargs ):
        var_filename = kwargs.get ( "var_filename", None )
        assert var_filename and os.path.isfile ( var_filename )

        Experiment.__init__ ( self )
        
        nc_fid  = Dataset ( var_filename, 'r' )
        nc_attrs = nc_fid.ncattrs ()
        for nc_attr in nc_attrs:
            self.__dict__[ nc_attr ] = nc_fid.getncattr ( nc_attr )
        nc_fid.close ()


    def __repr__ ( self ):
        return "%s object %s ( %s )" % ( self.__class__.__name__, self.name,
                               ",\n".join ( '{!s} : {!r}'.format ( k, print_none ( v.encode ( 'ascii', 'ignore' )) \
                                                                   if isinstance ( v, unicode ) else \
                                                                   print_none ( v )) \
                                            for (k,v) in self.__dict__.iteritems ()))
                                           
                                   
                                   
if __name__ == "__main__":

    var_filename = "/scratch/mpi/mpiles/m300050/benchmarking/example_data/CMIP5/amip/rsds_Amon_MPI-ESM-LR_amip_r3i1p1_197901-200812.nc"

    assert os.path.isfile ( var_filename )
    
    e = CMIP5_Experiment ( var_filename = var_filename )

    print e.institution

    print e
    
    
