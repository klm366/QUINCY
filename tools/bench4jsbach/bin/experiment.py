#! /usr/bin/env python

import os
from netCDF4 import Dataset 

from util import print_none


class Experiment ( object ):
    B4J_id    = None
    B4J_name  = None
    B4J_etype = None
    B4J_path  = None
    
    
    def __init__ ( self, **kwargs ):
        self.B4J_name     = kwargs.get ( "B4J_name" )
        self.B4J_etype    = kwargs.get ( "B4J_etype", None )
        self.B4J_path     = kwargs.get ( "B4J_path", None )
        self.B4J_buffer   = []
        

    def __repr__ ( self ):
        return "{} ( B4J_id = {}, B4J_name = {}, B4J_etype = {}, B4J_path = {}, B4J_buffer = {})". \
               format ( self.__class__.__B4J_name__,
                        print_none ( self.B4J_id ),
                        print_none ( self.B4J_name ),
                        print_none ( self.B4J_etype ),
                        print_none ( self.B4J_path ),
                        print_none ( self.B4J_buffer ))


class CMIP5_Experiment ( Experiment ):
    
    def __init__ ( self, **kwargs ):
        var_filename = kwargs.pop ( "var_filename", None )
        assert var_filename and os.path.isfile ( var_filename )

        Experiment.__init__ ( self, **kwargs )
        
        nc_fid  = Dataset ( var_filename, 'r' )
        nc_attrs = nc_fid.ncattrs ()
        for nc_attr in nc_attrs:
            self.__dict__[ nc_attr ] = nc_fid.getncattr ( nc_attr )
        nc_fid.close ()


    def __repr__ ( self ):
        return "%s object %s ( %s )" % ( self.__class__.__name__, self.B4J_name,
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
    
    
