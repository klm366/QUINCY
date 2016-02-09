#! /usr/bin/env python

import os
from netCDF4 import Dataset 

from util import print_none


class Experiment ( object ):
    B4J_id       = None
    B4J_name     = None
    B4J_etype    = None
    B4J_path     = None
    B4J_isref    = False
    B4J_iscan    = True
    B4J_varpaths = None
    B4J_canno    = None
    B4J_refno    = None
    
    def __init__ ( self, **kwargs ):
        self.B4J_name     = kwargs.get ( "B4J_name" )
        self.B4J_etype    = kwargs.get ( "B4J_etype", None )
        self.B4J_path     = kwargs.get ( "B4J_path", None )
        self.B4J_isref    = kwargs.get ( "B4J_isref", False )
        self.B4J_iscan    = kwargs.get ( "B4J_iscan", True )
        self.B4J_canno    = kwargs.get ( "B4J_canno", None )
        self.B4J_refno    = kwargs.get ( "B4J_refno", None )
        self.B4J_varpaths = {}
        

    def __repr__ ( self ):
        return "%s object %s ( %s )" % ( self.__class__.__name__, self.B4J_name,
                    ",\n".join ( '{!s} : {!r}'.format ( k, str ( print_none ( v ))) \
                    for (k,v) in sorted ( self.__dict__.iteritems ())))
    

    def teX_repr ( self ):
        from report import lateXify

        def truncate ( s, pl, sl ):
            import re
            assert isinstance ( s, str )
            ret = s
            if len ( s ) > pl + sl:
                ret = re.sub(r'^(.{%d}).*(.{%d})$' % (pl,sl), '\g<1>...\g<2>', s )
            return ret

        def force_breaks ( string, length ):
            assert isinstance ( string, str )
            insertnewline = "\\\\"
            ret = string
            l = len ( string )
            if l > length:
                withnewline = ''
                lastidx = 0
                for idx in range ( length - 1, l, length ):
                    withnewline += string [ lastidx : idx ] + insertnewline
                    lastidx = idx
                ret = '\\makecell[l]{ %s }' % withnewline
                
            return ret
        
        return """
\\begin{tabularx}{\\linewidth}{lL}
\\toprule
\\multicolumn{2}{l}{%s \\textbf{%s object {\\tt %s}}}  \\\\
\\addlinespace
\\cmidrule(r){1-2}
%s
\\\\ \\bottomrule
\\end{tabularx}
""" % ( '$\mathbf{C_{%d}}$:' % self.B4J_canno if self.B4J_iscan else '$\mathbf{R_{%d}}$:' % self.B4J_refno,
        lateXify ( self.__class__.__name__ ),
        lateXify ( self.B4J_name ), 
        """\\\\
""".join ( '{!s} & {!s}'.format ( lateXify ( str ( k )), force_breaks ( lateXify ( str ( v )), 50 )) \
           if not isinstance ( v, dict ) else '&'
           for ( k, v ) in sorted ( self.__dict__.iteritems ())))


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

                                                                                                                
if __name__ == "__main__":

    var_filename = "/scratch/mpi/mpiles/m300050/benchmarking/example_data/CMIP5/amip/rsds_Amon_MPI-ESM-LR_amip_r3i1p1_197901-200812.nc"

    assert os.path.isfile ( var_filename )
    
    e = CMIP5_Experiment ( var_filename = var_filename )

    print e.institution

    print e
    
    
