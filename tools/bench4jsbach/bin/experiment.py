#! /usr/bin/env python

import os, sys, subprocess
from netCDF4 import Dataset
import fnmatch, re , shutil
from util import print_none, subprocessing
from cdo import Cdo

ARCHIVE_PREFIX = "/hpss/arch/"
STREAM_SUFFIX  = "_stream"
CODES_SUFFIX   = "_codes"

HOST_AVAIL = { "thunder" : { "login_name" : "thunder4",
                             "names"      : [ "thunder4", "thunder5" ],
                             "cdo_path"   : "/sw/squeeze-x64/cdo-1.7.1-gcc48/bin" },
               "mistral" : { "login_name" : "mistral.dkrz.de",
                             "names"      : [ "mlogin100", "mlogin101", "mlogin102", "mlogin103", "mlogin104", "mlogin105" ],
                             "cdo_path"   : "/sw/rhel6-x64/cdo/cdo-1.7.0-magicsxx-gcc48/bin" },
               "archive" : { "login_name" : None,
                             "names"      : []}}

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
    B4J_obj         = None
    B4J_id          = None
    B4J_name        = None
    B4J_etype       = None
    B4J_atmres      = None
    B4J_oceres      = None
    B4J_startyear   = None
    B4J_finalyear   = None
    B4J_host        = None
    B4J_path        = None
    B4J_landseamask = None
    B4J_isref       = False
    B4J_iscan       = True
    B4J_index       = None
    B4J_canno       = None
    B4J_refno       = None
                    
    def __init__ (   self, b4j_obj, **kwargs ):
        self.B4J_obj   = b4j_obj
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
                exp_metadata = self.B4J_obj.wikipage.fetch_simulation_row (
                    "can" if self.B4J_iscan else "ref", self.B4J_name )
            except:
                print "No experiment description found for %s." % self.B4J_name
                sys.exit ( 1 )


            print exp_metadata
            self.B4J_host      = exp_metadata [ "host" ]
            self.B4J_path      = exp_metadata [ "data_storage" ]
            self.B4J_atmres    = exp_metadata [ "atm_res" ]
            self.B4J_oceres    = exp_metadata [ "oce_res" ]
            self.B4J_startyear = exp_metadata [ "simulation_period" ].split ( '-' )[ 0 ]
            self.B4J_finalyear = exp_metadata [ "simulation_period" ].split ( '-' )[ 1 ]

    def remotehost_processing ( self, cmd ):
        cmd = "ssh %s@%s '%s'" % ( self.B4J_obj.config [ 'author' ],
                                   HOST_AVAIL [ self.B4J_host ][ 'login_name' ],
                                   cmd )
        return subprocessing ( cmd )


    def remotehost_fetching ( self, remote_path, target_path ):
        cmd = "scp %s@%s:%s %s" % ( self.B4J_obj.config [ 'author' ],
                                    HOST_AVAIL [ self.B4J_host ][ 'login_name' ],
                                    remote_path,
                                    target_path )
        return subprocessing ( cmd )

    
    def accumulate_streams ( self, orig_dir, orig_fn_tmpl, startyear, finalyear, selcodes, format, codes_path, target_path, from_remote ):
        self.B4J_obj.logger.debug ( "Starting." )

        sy  = int ( startyear )
        fy  = int ( finalyear )

        if not from_remote:
            cdo = Cdo ()
            try:
                os.remove ( target_path )
            except OSError:
                pass
            cdo.mergetime ( input = [ os.path.join ( orig_dir, orig_fn_tmpl % { 'year' : y }) for y in range ( sy, fy + 1 )],
                            output = target_path )
        else:
            remote_temp_dir  = "B4J_TEMP_DIR"
            remote_temp_fn   = ( orig_fn_tmpl % { "year" : "%s-%s" % ( startyear, finalyear )}).replace ( '.grb', '.nc' )
#            remote_temp_fn = remote_temp_fn.replace ( '.grb', '.nc' )
            remote_temp_path = os.path.join ( remote_temp_dir, remote_temp_fn )

            

            ofnt = orig_fn_tmpl.replace ( '%(year)s', '${y}' )

            remote_cmd = """
            echo 'on server' $(hostname);
            mkdir -p %s;
            for y in `seq %d 1 %d`;do
            ARRAY[$i]=%s/%s;
            echo ${ARRAY[i]}
            i=$i+1;
            done;
            export PATH=%s:$PATH;
            rm -f %s;
            cdo -f %s -t %s select,code=%s ${ARRAY[*]} %s;
            """ % ( remote_temp_dir, sy, fy, orig_dir, ofnt,
            HOST_AVAIL [ self.B4J_host ][ 'cdo_path' ],
            remote_temp_path, format, codes_path, selcodes, remote_temp_path )

            print remote_cmd

            self.remotehost_processing ( remote_cmd )
            self.remotehost_fetching ( remote_temp_path, target_path )
            self.remotehost_processing ( "rm -rf %s" % remote_temp_dir )
        
        self.B4J_obj.logger.debug ( "Exiting." )


    def fetch_outdata_new ( self ):
        from bench import bench4jsbach_scripts, get_hostname
        
        self.B4J_obj.logger.debug ( "Starting" )

        host         = get_hostname ()
        from_archive = self.B4J_host == 'archive'
        from_remote  = host not in HOST_AVAIL [ self.B4J_host ][ 'names' ]
        accumulated  = True if from_archive else False              # ATT: if logic changes, program flow has to be adapted
        modeltype    = self.__class__.__name__.split ( "_" )[ 0 ]

        self.B4J_obj.logger.debug ( """Prepare fetching%s%s output data%s.""" %
                                    (( " remote" if from_remote else "local" ),
                                     ( " acculated" if accumulated else " not accumulated" ),
                                     ( " from %s" % self.B4J_host if from_remote else "" )))
        
        
        if from_remote: fetchfile_list = []
        
        for field_name, field_config in self.B4J_obj.diagnostic_fields.iteritems ():
            assert 'fetch' in field_config, "No 'fetch'-section defined for %s." % field_name
            for file_name, file_config in field_config ['fetch' ].iteritems ():
                assert modeltype in file_config, "No '%s'-section defined for %s." % file_name
                
                print file_name

                if not accumulated:
                    assert "submodel" in file_config [ modeltype ] and\
                           "stream" in file_config [ modeltype ]

                    dest_path_token =  file_name + STREAM_SUFFIX                     
                    fn = self.outdata_fn_tmpl % { "expname"           : self.B4J_name,
                                                  "submodel"          : file_config [ modeltype ][ "submodel" ],
                                                  "stream"            : file_config [ modeltype ][ "stream" ],
                                                  "cond_frequencefix" : self.outdata_frequence_infixes [ 'monthly' ],
                                                  "startyear"         : field_config [ 'startyear' ],
                                                  "cond_finalyearstr" : self.outdata_finalyear_tmpl % { 'finalyear': field_config [ 'finalyear' ]},
                                                  "suffix"            : field_config [ 'suffix' ]}#self.outdata_suffix }
                    self.B4J_index [ dest_path_token ] = os.path.join ( self.B4J_obj.downloads_path, fn )
                    print self.B4J_index [ dest_path_token ]
                    
                    orig_path    = os.path.join ( self.B4J_path, self.outdata_subpath, file_config [ modeltype ][ "submodel" ])
                    orig_fn_tmpl = self.outdata_fn_tmpl % { "expname"           : self.B4J_name,
                                                            "submodel"          : file_config [ modeltype ][ "submodel" ],
                                                            "stream"            : file_config [ modeltype ][ "stream" ],
                                                            "cond_frequencefix" : self.outdata_frequence_infixes [ 'monthly' ],
                                                            "startyear"         : "%(year)s",
                                                            "cond_finalyearstr" : "",
                                                            "suffix"            : self.outdata_suffix }
                    
                    codes_path  = os.path.join ( self.B4J_path, self.B4J_name, self.codesfiles_subpath,
                                                 self.codesfiles_fn_tmpl % { "expname": self.B4J_name,
                                                                             "submodel": file_config [ modeltype ][ "submodel" ],
                                                                             "stream": file_config [ modeltype ][ "stream" ],
                                                                             "suffix": self.codesfiles_suffix })
                    print codes_path
                    self.accumulate_streams ( orig_path,
                                              orig_fn_tmpl,
                                              field_config [ 'startyear' ],
                                              field_config [ 'finalyear' ],
                                              file_config [ modeltype ][ 'code' ],
                                              field_config [ 'suffix' ],
                                              codes_path,
                                              self.B4J_index [ dest_path_token ],
                                              from_remote )
                    
                                    
        self.B4J_obj.logger.debug ( "Exiting" )
        
    

    def fetch_outdata_all ( self ):
        from bench import bench4jsbach_scripts, get_hostname

        self.B4J_obj.logger.debug ( "Starting" )
        
        modeltype    = self.__class__.__name__.split ( "_" )[ 0 ]
        from_archive = self.B4J_path.startswith ( ARCHIVE_PREFIX ) 
        host         = get_hostname ()
        from_remote  = host not in HOST_AVAIL [ self.B4J_host ] or from_archive
        accumulated  = True if from_archive else False

        print "from_remote:  %d" % ( 1 if from_remote else 0 )
        print "accumulated: %d" % ( 1 if accumulated else 0 )
        

        self.B4J_obj.logger.debug ( "Prepare fetching %s data%s." % ( \
            "remote" if from_remote else "local",
            " from archive" if from_archive else "" ))

        if from_remote: fetchfile_list = []
        for field_name, field_config in self.B4J_obj.diagnostic_fields.iteritems ():
            assert 'fetch' in field_config, "No 'fetch'-section defined for %s." % field_name
            for file_name, file_config in field_config ['fetch' ].iteritems ():
                assert modeltype in file_config, "No '%s'-section defined for %s." % file_name
                if field_config [ "frequence" ] == "monthly":
                    frequence_infix = self.outdata_frequence_infixes [ field_config [ "frequence" ]]
                    finalyearstring = self.outdata_finalyear_tmpl % { 'finalyear': self.B4J_finalyear }
                else:
                    frequence_infix = ""
                    finalyearstring = ""
                fetchpath = self.B4J_path + os.sep + self.outdata_subpath + os.sep + file_config [ modeltype ][ "submodel" ]
                print "fetchpath: ", fetchpath

                if accumulated:
                    fetchfile = self.outdata_fn_tmpl % { "expname"           : self.B4J_name,
                                                         "submodel"          : file_config [ modeltype ][ "submodel" ],
                                                         "stream"            : file_config [ modeltype ][ "stream" ],
                                                         "cond_frequencefix" : frequence_infix,
                                                         "startyear"         : self.B4J_startyear,
                                                         "cond_finalyearstr" : finalyearstring,
                                                         "suffix"           : self.outdata_suffix }
                    print "fetchfile: ", fetchfile
                
                if from_archive:
                    self.B4J_index [ file_name + STREAM_SUFFIX ] = os.path.join ( self.B4J_obj.downloads_path, fetchfile )
                    if not os.path.isfile ( self.B4J_index [ file_name + STREAM_SUFFIX ]):
                        fetchfile_list.append (( fetchpath, fetchfile ))
                elif accumulated:
                    self.B4J_index [ file_name + STREAM_SUFFIX ] = os.path.join ( fetchpath, fetchfile )
                else:
                    
                    accumulated_fn = self.outdata_fn_tmpl % { "expname"           : self.B4J_name,
                                                              "submodel"          : file_config [ modeltype ][ "submodel" ],
                                                              "stream"            : file_config [ modeltype ][ "stream" ],
                                                              "cond_frequencefix" : frequence_infix,
                                                              "startyear"         : field_config [ 'startyear' ],
                                                              "cond_finalyearstr" : self.outdata_finalyear_tmpl % { 'finalyear': field_config [ 'finalyear' ]},
                                                              "suffix"            : self.outdata_suffix }
                    self.B4J_index [ file_name + STREAM_SUFFIX ] = os.path.join ( self.B4J_obj.downloads_path, accumulated_fn )
                    origin_template = self.outdata_fn_tmpl % { "expname"           : self.B4J_name,
                                                               "submodel"          : file_config [ modeltype ][ "submodel" ],
                                                               "stream"            : file_config [ modeltype ][ "stream" ],
                                                               "cond_frequencefix" : frequence_infix,
                                                               "startyear"         : "%(year)d",
                                                               "cond_finalyearstr" : "",
                                                               "suffix"            : self.outdata_suffix }
                    self.accumulate_streams ( fetchpath,
                                              origin_template,
                                              field_config [ 'startyear' ],
                                              field_config [ 'finalyear' ],
                                              self.B4J_index [ file_name + STREAM_SUFFIX ] )
                    
                print "self.B4J_index [ file_name + STREAM_SUFFIX ]: ", self.B4J_index [ file_name + STREAM_SUFFIX ]

                fetchcodespath    = self.B4J_path + os.sep + self.subpath 
                fetchcodestarfile = self.codesfiles_tar_tmpl % { "expname" : self.B4J_name }
                fetchcodesfile    = self.codesfiles_fn_tmpl % { "expname"  : self.B4J_name,
                                                        "submodel" : file_config [ modeltype ][ "submodel" ],
                                                        "stream"   : file_config [ modeltype ][ "stream" ],
                                                        "suffix"  : self.codesfiles_suffix }

                print "fetchcodespath: ", fetchcodespath
                print "fetchcodestarfile: ", fetchcodestarfile
                print "fetchcodesfile: ", fetchcodesfile

                if from_archive:
                    self.B4J_index [ file_name + CODES_SUFFIX ] = os.path.join ( self.B4J_obj.downloads_path, fetchcodesfile )
                    if not os.path.isfile ( self.B4J_index [ file_name + CODES_SUFFIX ]) and \
                           not os.path.isfile ( fetchcodestarfile ):
                        fetchfile_list.append (( fetchcodespath, fetchcodestarfile ))
                else:
                    print 
                    # to check: path of codes file in not archived mpiesm1 experiment
                    self.B4J_index [ file_name + CODES_SUFFIX ] = os.path.join (  self.B4J_path, self.codefiles_subpath, fetchcodesfile )
                
                print "self.B4J_index [ file_name + CODES_SUFFIX ]: ", self.B4J_index [ file_name + CODES_SUFFIX ]    


        if from_archive and fetchfile_list:
            fetchfile_list = list ( set ( fetchfile_list ))
            workingdir = os.getcwd ()
            os.chdir ( self.B4J_obj.downloads_path )
            exe = os.path.join ( bench4jsbach_scripts, "b4j_fetch_archived_file.ksh" )
            assert os.path.isfile ( exe ) and os.access ( exe, os.X_OK ), "File %s is not executable." % exe
            self.B4J_obj.logger.debug ( "Fetching outdata and codes files from archive for %s into %s. This might take a while ..." % (
                self.B4J_name, self.B4J_obj.downloads_path ))
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
            self.B4J_obj.logger.debug ( "Fetched outdata from archive for %s successfully into %s." % (
                self.B4J_name, self.B4J_obj.downloads_path ))
          
            # extract codesfiles and delete tarball
            if ( fetchcodespath, fetchcodestarfile ) in fetchfile_list:
                for cmd in [[ "tar -xf %s --wildcards --no-anchored --strip=1 '*.%s'" %
                              ( fetchcodestarfile, self.codesfiles_suffix  )],
                            [ "rm -f %s" % fetchcodestarfile ]]:
                    sp = subprocess.Popen( cmd, shell=True,
                                           stdout=subprocess.PIPE,  stderr=subprocess.PIPE )
                    stdout, stderr = sp.communicate()
                    if sp.returncode:
                        print "stdout: ", stdout
                        print "stderr: ", stderr
                        print "returncode: %d" % sp.returncode
                        sys.exit ( sp.returncode )
                    else:
                        print stdout
                self.B4J_obj.logger.debug ( "Extracted all codesfiles from %s to %s, deleted tarball." %
                                               ( fetchcodestarfile, self.B4J_obj.downloads_path ))
                
            os.chdir ( workingdir )


        # select years and codes

        cdo = Cdo ()

        for fieldname, fielddict in self.B4J_obj.diagnostic_fields.iteritems ():
            assert 'fetch' in fielddict, "No 'fetch'-section defined for %s." % fieldname
            yearlist = range ( int ( fielddict [ 'startyear' ]), int ( fielddict [ 'finalyear' ]) + 1 )
            for fileid, filedict in fielddict ['fetch' ].iteritems ():
                self.B4J_index [ fileid ] = self.B4J_obj.path_template % { 
                    'path'           : fielddict [ "processed_path" ],
                    'classification' : fileid,
                    'expname'        : self.B4J_name,
                    'startdate'      : fielddict [ 'startyear' ],
                    'finaldate'      : fielddict [ 'finalyear' ],
                    'suffix'         : fielddict [ 'suffix' ]}
                assert modeltype in filedict, "No '%s'-section defined for %s." % fileid

                print self.B4J_index [ fileid + CODES_SUFFIX ]
                assert os.path.isfile ( self.B4J_index [ fileid + CODES_SUFFIX ])
                
                cdo.selcode (  filedict [ modeltype ][ 'code' ],
                              input = " -selyear," + ','.join ( repr ( y ) for y in yearlist ) +\
                               " " + self.B4J_index [ fileid + STREAM_SUFFIX ],
                              output = self.B4J_index [ fileid ],
                              options = ' -f %s' % fielddict [ 'suffix' ] +\
                                        ' -t %s' % self.B4J_index [ fileid + CODES_SUFFIX ])

        # fetch land-sea-mask
        lsm_path = self.land_sea_mask_path_tmpl % { "atm_res" : self.B4J_atmres,
                                                    "oce_res" : self.B4J_oceres }
        assert os.path.isfile ( lsm_path ), "No file found for %s." % lsm_path 
        self.land_sea_mask  = self.B4J_obj.path_template % { 
            'path'           : fielddict [ "processed_path" ],
            'classification' : "land_sea_mask",
            'expname'        : self.B4J_name,
            'startdate'      : fielddict [ 'startyear' ],
            'finaldate'      : fielddict [ 'finalyear' ],
            'suffix'         : fielddict [ 'suffix' ]}
        slf_varname = "slf"
        cdo.selname ( slf_varname, input = lsm_path, output = self.land_sea_mask )
        self.B4J_obj.logger.debug ( "Fetched land sea mask '%s'." %  self.land_sea_mask )


        self.B4J_obj.logger.debug ( "Exiting" )


    def fetch_outdata ( self, arg1, arg2 ):
        print "hello from superclass experiments method fetch outdata."

    def __repr__ ( self ):
        return "\n%s object %s ( \n\t%s )" % ( self.__class__.__name__, self.B4J_name,
                    ",\n\t".join ( '{!s}\t: {!r}'.format ( k, str ( print_none ( v ))) \
                    for (k,v) in sorted ( self.__dict__.iteritems ())))

    def teX_short_repr ( self, shortname = 'can'  ):
        from report import teXify
        
        
        return '%s object \\textbf{%s {\\tt %s}}' % ( teXify ( self.__class__.__name__ ),
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
\\multicolumn{2}{l}{\\phantomsection\\label{%s}\\hypertarget{%s} %s object \\textbf{%s {\\tt %s}}} \\\\
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
           if not isinstance ( v, dict ) and not k == "B4J_obj" else '&'
           for ( k, v ) in sorted ( self.__dict__.iteritems ())))


class CMIP5_Experiment ( Experiment ):
    suffix = 'nc'
    
    def __init__ ( self,  b4j_obj, **kwargs ):
        var_filename = kwargs.pop ( "var_filename", None )
        assert var_filename and os.path.isfile ( var_filename )

        Experiment.__init__ ( self, b4j_obj, **kwargs )
        
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


# status: only for archived experiment
class MPIESM1_Experiment ( Experiment ):
    codesfiles_suffix         = "codes"
    codesfiles_subpath        = 'log'
    codesfiles_fn_tmpl        = "%(expname)s_%(submodel)s_%(stream)s.%(suffix)s"
    codesfiles_tar_tmpl       = "%(expname)s_log.tar"
    
    outdata_fn_tmpl           = """%(expname)s_%(submodel)s_%(stream)s_%(cond_frequencefix)s%(startyear)s%(cond_finalyearstr)s.%(suffix)s"""
    outdata_finalyear_tmpl    = """-%(finalyear)s"""
    outdata_frequence_infixes = { "monthly" : "mm_" }
    outdata_subpath           = None
    outdata_suffix            = "grb"

    land_sea_mask_path_tmpl   = "/pool/data/JSBACH/input/r0007/%(atm_res)s/jsbach_%(atm_res)s%(oce_res)s_11tiles_5layers_1850.nc"

    def __init__ ( self, b4j_obj, **kwargs ):
        print "###################### go on here, subclass %s ####################" % self.__class__.__name__
        Experiment.__init__ ( self, b4j_obj, **kwargs )
        #super(Experiment, self).__init__(**kwargs)
        self.subpath           = self.B4J_name
        self.outdata_subpath   = os.path.join ( self.subpath, "outdata" )
        self.codefiles_subpath = os.path.join ( self.subpath, "log" )
        #print "MPIESM1_Experiment attributes:"
        #print self.__dict__


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
    
    
