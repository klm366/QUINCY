#! /usr/bin/env python

import os, sys, subprocess
from netCDF4 import Dataset
import fnmatch, re , shutil
from redmine_read import Wikitable_experiment_description
from util import print_none

ARCHIVE_PREFIX = "/hpss/arch/"



def find_files ( pattern, path, single = True, raise_e = False ):
    assert os.path.isdir ( path ) and isinstance ( pattern, str )
    
    paths = []
    for root, dirnames, filenames in os.walk ( path ):
        files   = fnmatch.filter ( filenames, pattern )
        for f in files:
            paths.append ( os.path.join ( root, f ))
    
    n_files = len ( paths )

    if single:
        if not n_files:
            if raise_e:
                raise Exception ( "no file" )
            else:
                ret = None
        elif n_files > 1:
            if raise_e:
                raise Exception ( "multiple files" )
            else:
                warnings.warn ( "multiple files, take first" )
                ret = paths [ 0 ]
        else:
            ret = paths [ 0 ]
    else:
        if not n_files and raise_e:
            raise Exception ( "no files" )
        ret = paths

    return ret


class Experiment ( object ):
    B4J_parent    = None
    B4J_id        = None
    B4J_name      = None
    B4J_etype     = None
    B4J_startyear = None
    B4J_finalyear = None
    B4J_path      = None
    B4J_isref     = False
    B4J_iscan     = True
    B4J_index     = None
    B4J_canno     = None
    B4J_refno     = None
    
    def __init__ ( self, parent, **kwargs ):
        self.B4J_parent   = parent
        self.B4J_name     = kwargs.get ( "B4J_name" )
        self.B4J_etype    = kwargs.get ( "B4J_etype", None )
        self.B4J_path     = kwargs.get ( "B4J_path", None )
        self.B4J_isref    = kwargs.get ( "B4J_isref", False )
        self.B4J_iscan    = kwargs.get ( "B4J_iscan", True )
        self.B4J_canno    = kwargs.get ( "B4J_canno", None )
        self.B4J_refno    = kwargs.get ( "B4J_refno", None )
        self.B4J_index    = {}

        if not self.B4J_path:
            try:
                redmine_row =  Wikitable_experiment_description (
                    kwargs.get ( 'redmine_url' ),
                    kwargs.get ( 'redmine_user' ),
                    kwargs.get ( 'redmine_key' ),
                    kwargs.get ( 'redmine_project' ),
                    kwargs.get ( 'redmine_wikipage' ),
                    "can" if self.B4J_iscan else "ref",
                    self.B4J_name )
            except:
                print "No experiment description found for %s." % self.B4J_name
                sys.exit ( 1 )
            # print '\nRedmine reference experiment description:  \n%s\n' % \
            #       ''.join('{!s} = {!r},\n'.\
            #               format ( key, val) for ( key, val ) in redmine_row.__dict__.iteritems ())
            self.B4J_path = redmine_row.data_storage
            self.B4J_startyear = redmine_row.simulation_period.split ( '-' )[ 0 ]
            self.B4J_finalyear = redmine_row.simulation_period.split ( '-' )[ 1 ]


    def fetch_outdata_all ( self ):
        from bench import bench4jsbach_scripts

        self.B4J_parent.logger.debug ( "Starting" )

        # if data files are in archive
        if self.B4J_path.startswith ( ARCHIVE_PREFIX ):
            fetchfile_list=[]
            for field_name, field_config in self.B4J_parent.diagnostic_fields.iteritems ():
                print "######### ", field_name
                assert 'fetch' in field_config, "No 'fetch'-section defined for %s." % field_name
                for file_name, file_config in field_config ['fetch' ].iteritems ():
                    print file_name
                    print file_config
                    modeltype = self.__class__.__name__.split ( "_" )[ 0 ]
                    assert modeltype in file_config, "No '%s'-section defined for %s." % file_name
                    print file_config [ modeltype ]
                    if field_config [ "frequence" ] == "monthly":
                        frequence_fix = self.outdata_frequence_fixes [ field_config [ "frequence" ]]
                        finalyearstring = self.outdata_finalyear_tmpl % { 'finalyear': self.B4J_finalyear }
                        print frequence_fix
                        print finalyearstring
                    else:
                        frequence_fix = ""
                        finalyearstring = ""
                    fetchpath = self.B4J_path + os.sep + self.outdata_subpath + os.sep + file_config [ modeltype ][ "submodel" ]

                    fetchfile = self.outdata_fn_tmpl % { "expname"           : self.B4J_name,
                                                         "submodel"          : file_config [ modeltype ][ "submodel" ],
                                                         "stream"            : file_config [ modeltype ][ "stream" ],
                                                         "cond_frequencefix" : frequence_fix,
                                                         "startyear"         : self.B4J_startyear,
                                                         "cond_finalyearstr" : finalyearstring,
                                                         "postfix"           : self.outdata_postfix }
                    fetchfile_list.append (( fetchpath, fetchfile ))
                                                         
            fetchfile_list = list ( set ( fetchfile_list ))
            print fetchfile_list
            print "########### end field fetch file config ############"

            # fetch data files from archive
            workingdir = os.getcwd ()
            os.chdir ( self.B4J_parent.temp_path )
            exe = os.path.join ( bench4jsbach_scripts, "b4j_getdata_pftp" )
            assert os.path.isfile ( exe ) and os.access ( exe, os.X_OK ), "File %s is not executable." % exe
            for fetchfile in fetchfile_list:
                sp = subprocess.Popen( [ exe + " " + fetchfile [ 0 ] + " " + fetchfile [ 1 ]], shell=True,
                                       stdout=subprocess.PIPE,  stderr=subprocess.PIPE )
                stdout, stderr = sp.communicate()
                if sp.returncode:
                    print "stdout: ", stdout
                    print "stderr: ", stderr
                    print "returncode: %d" % sp.returncode
                    sys.exit ( sp.returncode )
                else:
                    print stdout
            os.chdir ( workingdir )
            self.B4J_parent.logger.debug ( "Fetched outdata for %s successfully into %s." % (
                self.B4J_name, self.B4J_parent.temp_path ))

        
        self.B4J_parent.logger.debug ( "Exiting" )


    def fetch_outdata ( self, arg1, arg2 ):
        print "hello from superclass experiments method fetch outdata."

    def __repr__ ( self ):
        return "\n%s object %s ( \n\t%s )" % ( self.__class__.__name__, self.B4J_name,
                    ",\n\t".join ( '{!s}\t: {!r}'.format ( k, str ( print_none ( v ))) \
                    for (k,v) in sorted ( self.__dict__.iteritems ())))

    def teX_short_repr ( self, shortname = 'can'  ):
        from report import teXify
        
        return '%s \\textbf{%s object {\\tt %s}}' % ( teXify ( self.__class__.__name__ ),
                                                      teXify ( self.B4J_name ),
                                                      '($\mathbf{C_{%d}}$)' % self.B4J_canno if shortname == 'can' else \
                                                      '($\mathbf{R_{%d}}$)' % self.B4J_refno )
    

    def teX_repr ( self ):
        from report import teXify

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
\\multicolumn{2}{l}{\\phantomsection\\label{%s}\\hypertarget{%s} %s \\textbf{%s object {\\tt %s}}} \\\\
\\addlinespace
\\cmidrule(r){1-2}
%s
\\\\ \\bottomrule
\\end{tabularx}\-\\\\\\bigskip




""" % ( self.B4J_name,
        self.B4J_name,
        teXify ( self.__class__.__name__ ),
        teXify ( self.B4J_name ),
        '($\mathbf{C_{%d}}$)' % self.B4J_canno if self.B4J_iscan else '($\mathbf{R_{%d}}$)' % self.B4J_refno,
        """\\\\
""".join ( '{!s} & {!s}'.format ( teXify ( str ( k )), force_breaks ( teXify ( str ( v )), 50 )) \
           if not isinstance ( v, dict ) else '&'
           for ( k, v ) in sorted ( self.__dict__.iteritems ())))


class CMIP5_Experiment ( Experiment ):
    suffix = 'nc'
    
    def __init__ ( self,  parent, **kwargs ):
        var_filename = kwargs.pop ( "var_filename", None )
        assert var_filename and os.path.isfile ( var_filename )

        Experiment.__init__ ( self, parent, **kwargs )
        
        nc_fid  = Dataset ( var_filename, 'r' )
        nc_attrs = nc_fid.ncattrs ()
        for nc_attr in nc_attrs:
            self.__dict__[ nc_attr ] = nc_fid.getncattr ( nc_attr )
        nc_fid.close ()

    # fetch_outdata is just a workaround function.
      # TODO: WRITE FUNCTION FETCH_OUTDATA FETCHING VARIABLE FOR CODE AND DATES INDEPENDENT OF B4J DATASTRUCTURE 
    def fetch_outdata ( self, varname, newpath ):
        pattern = '*%s*%s*.%s' % ( varname, self.B4J_name, self.suffix )
        oldpath = find_files ( pattern, self.B4J_path, single = True, raise_e = False )
        if not oldpath:
            raise Exception ( "ERROR: No file matching pattern '%s' found." % pattern )
            sys.exit ( 1 )
        shutil.copyfile ( oldpath, newpath )

        return oldpath



class MPIESM1_Experiment ( Experiment ):
    outdata_fn_tmpl         ="""%(expname)s_%(submodel)s_%(stream)s_%(cond_frequencefix)s%(startyear)s%(cond_finalyearstr)s.%(postfix)s"""
    outdata_finalyear_tmpl  ="""-%(finalyear)s"""
    outdata_frequence_fixes ={ "monthly" : "mm_" }
    outdata_subpath         = None
    outdata_postfix         = "grb"

    def __init__ ( self, parent, **kwargs ):
        print "###################### go on here  ####################"
        Experiment.__init__ ( self, parent, **kwargs )
        #super(Experiment, self).__init__(**kwargs)
        self.outdata_subpath = self.B4J_name + os.sep + "outdata"
        print "MPIESM1_Experiment attributes:"
        print self.__dict__


    def fetch_oudata ( self, varname, newpath ):
        print "fetch_outdata: Starting."
        ret = []
        print "fetch_outdata: Exiting."        

                                                                                                                
if __name__ == "__main__":

    var_filename = "/scratch/mpi/mpiles/m300050/benchmarking/example_data/CMIP5/amip/rsds_Amon_MPI-ESM-LR_amip_r3i1p1_197901-200812.nc"

    assert os.path.isfile ( var_filename )
    
    e = CMIP5_Experiment ( var_filename = var_filename )

    print e.institution

    print e
    
    
