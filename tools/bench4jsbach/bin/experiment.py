#! /usr/bin/env python
"""Module contains experiment related classes, methods and metadata."""

__author__="Deike Kleberg"

# show all methods of class Experiment
# >>> from experiment import Experiment
# >>> import inspect
# >>> import pprint
# >>> mem = inspect.getmembers ( Experiment, predicate=inspect.ismethod)
# >>> pprint.pprint ([ i for i,j in mem])

            # import pprint
            # pprint.pprint ( dict ( classinstance.__dict__ ))
            # print

# ['__init__',
#  '__repr__',
#  'accumulate_streams',
#  'fetch_archived_streams',
#  'fetch_land_sea_mask',
#  'fetch_outdata',
#  'remotehost_fetching',
#  'remotehost_processing',
#  'teX_repr',
#  'teX_short_repr']


import os, sys, subprocess
from netCDF4 import Dataset
from cdo import Cdo

from util import (
    getdomain,
    gethostname,
    print_none,
    subprocessing,
    stdout2devnull )

#: a test documentation for the archive prefix constant
ARCHIVE_PREFIX = "/hpss/arch/"
"""Property docstring."""

STREAM_SUFFIX  = "_stream"
CODES_SUFFIX   = "_codes"

HOST_AVAIL = { "thunder" : { "login_name" : "thunder4",
                             "names"      : [ "thunder4", "thunder5" ],
                             "domain"     : "zmaw.de",
                             "cdo_path"   : "/sw/squeeze-x64/cdo-1.7.1-gcc48/bin" },
               "mistral" : { "login_name" : "mistral.dkrz.de",
                             "names"      : [ "mlogin100", "mlogin101",
                                              "mlogin102", "mlogin103",
                                              "mlogin104", "mlogin105",
                                              "mistralpp1", "mistralpp2", "mistralpp3", "mistralpp4", "mistralpp5" ],
                             "domain"     : "hpc.dkrz.de",
                             "cdo_path"   : "/sw/rhel6-x64/cdo/cdo-1.7.0-magicsxx-gcc48/bin" },
               "archive" : { "login_name" : None,
                             "names"      : []},
               "scalar"  : { "login_name" : "scalar",
                             "names"      : [ "scalar", "localhost" ]}}

              

def find_files ( pattern, path, single = True, raise_e = False ):
    """Function to find files with a name containing a given infix.
    """
    
    assert os.path.isdir ( path ) and isinstance ( pattern, str )
    
    #: a test documentation for the paths variable
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
    """Base class representing methods and metadata of an experiment.

Specialized classes for model simulations and observational datasets
are derived from it.
    """
    
    B4J_obj         = None
    B4J_id          = None
    B4J_name        = None
    B4J_name_orig   = None
    B4J_etype       = None
    B4J_atmres      = None
    B4J_oceres      = None
    B4J_startyear   = None
    B4J_endyear     = None
    B4J_initialyear = None
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
        """Objects are initialized giving their attributes.

In a bench4jsbach run this is done by the runskript reading a configuration file.      
        """
        self.B4J_obj         = b4j_obj
        self.B4J_name        = kwargs.get ( "B4J_name" )
        self.B4J_name_orig   = kwargs.get ( "B4J_name_orig" )
        self.B4J_etype       = kwargs.get ( "B4J_etype", None )
        #: this is a comment on B4J_path
        self.B4J_path        = kwargs.get ( "B4J_path", None )
        self.B4J_host        = kwargs.get ( "B4J_host", None )
        self.B4J_isref       = kwargs.get ( "B4J_isref", False )
        self.B4J_iscan       = kwargs.get ( "B4J_iscan", True )
        self.B4J_canno       = kwargs.get ( "B4J_canno", None )
        self.B4J_refno       = kwargs.get ( "B4J_refno", None )
        self.B4J_startyear   = kwargs.get ( "B4J_startyear", None )
        self.B4J_endyear     = kwargs.get ( "B4J_endyear", None )
        self.B4J_initialyear = kwargs.get ( "B4J_initialyear", None )
        self.B4J_finalyear   = kwargs.get ( "B4J_finalyear", None )
        self.B4J_index       = {}

        self.B4J_obj.logger.debug ( 'kwargs = { %s }' % ''.join("{!s} = {!r}, ".format(key,val) 
                                                        for (key,val) in kwargs.items()))

        
        if not self.B4J_path:
            exp_metadata = self.B4J_obj.wikipage.fetch_simulation_row_split_period (
                "can" if self.B4J_iscan else "ref", self.B4J_name_orig ) 
            if not exp_metadata:
                raise KeyError ( "No experiment description found for %s." % self.B4J_name_orig )

            self.B4J_obj.logger.debug ( "Read experiment metadata from wikipage :\n %s " % exp_metadata )
            
            self.B4J_host        = exp_metadata [ "host" ]
            self.B4J_path        = exp_metadata [ "data_storage" ]
            self.B4J_initialyear = self.B4J_initialyear or exp_metadata [ "initialyear" ]
            self.B4J_finalyear   = self.B4J_finalyear or exp_metadata [ "finalyear" ]
            self.B4J_atmres      = exp_metadata [ "atm_res" ]
            self.B4J_oceres      = exp_metadata [ "oce_res" ]

        host         = gethostname ()
        self.from_archive = self.B4J_host == 'archive'
        self.from_remote  = host not in HOST_AVAIL [ self.B4J_host ][ 'names' ]
        self.accumulated  = True if self.from_archive else False
        # ATT: if logic changes, program flow has to be adapted


    def check_output ( self, filepath ): 
        cdo = Cdo ()       
        if os.path.isfile ( filepath ):         
            with stdout2devnull ():
                try:
                    filetype = cdo.filedes ( input=filepath )
                except:
                    filetype = "unsupported filetype"
                else:
                    filetype = filetype [ 1 ]
                    
            print "\nFile %s exists:" % ( filepath )
            print "size %s byte, type %s" % ( os.stat ( filepath).st_size, filetype )

        else:
            print "\nERROR: File '%s' does not exist." % ( filepath )
            sys.exit ( 1 )
        self.B4J_obj.logger.debug ( "Exiting." )


    def index_check_output ( self, tokenlist ):
        self.B4J_obj.logger.debug ( "Starting." )
        for token in tokenlist:
            self.check_output ( self.B4J_index [ token ])
        self.B4J_obj.logger.debug ( "Exiting." )
                

    def remotehost_processing ( self, cmd ):
        cmd = "ssh %s@%s '%s'" % ( self.B4J_obj.config [ 'author' ],
                                   HOST_AVAIL [ self.B4J_host ][ 'login_name' ],
                                   cmd )

        self.B4J_obj.logger.debug ( "Exiting" )
        return subprocessing ( cmd )


    def remotehost_fetching ( self, remote_path, target_path ):
        self.B4J_obj.logger.debug ( "Starting." )
        cmd = "scp %s@%s:%s %s" % ( self.B4J_obj.config [ 'author' ],
                                    HOST_AVAIL [ self.B4J_host ][ 'login_name' ],
                                    remote_path,
                                    target_path )
        

        self.B4J_obj.logger.debug ( "Exiting" )
        return subprocessing ( cmd )

    
    def fetch_land_sea_mask ( self, origin_path, target_path, varname ):
        self.B4J_obj.logger.debug ( "Starting." )
        
        target_fn = os.path.basename ( target_path )

        remote_cmd = """
        echo 'on server' $(hostname);
        export PATH=%s:$PATH;
        rm -f %s;
        cdo selname,%s %s %s
        """ % (
            HOST_AVAIL [ self.B4J_host ][ 'cdo_path' ],
            target_fn, varname, origin_path, target_fn )


        self.remotehost_processing ( remote_cmd )
        self.B4J_obj.logger.debug ( "Executed command successfully." )
        self.remotehost_fetching ( target_fn, self.land_sea_mask )
        self.B4J_obj.logger.debug ( "Copied file to local machine successfully." )
        self.remotehost_processing ( "rm -f %s" % target_fn )
        self.B4J_obj.logger.debug ( "Deleted remote file successfully." )

        self.B4J_obj.logger.debug ( "Exiting" )


    def accumulate_streams ( self, orig_dirpath, orig_fn_tmpl, startyear, endyear,
                             selcodes, suffix, codes_path, target_path, from_remote ):
        self.B4J_obj.logger.debug ( "Starting." )
        print "orig_dirpath:", orig_dirpath
        print "orig_fn_tmpl:", orig_fn_tmpl
        print "startyear:", startyear
        print "endyear:", endyear
        print "selcodes:", selcodes
        print "suffix:", suffix
        print "codes_path:", codes_path
        print "target_path:", target_path
        print "from_remote:", from_remote

        sy  = int ( startyear )
        fy  = int ( endyear )
        
        if not os.path.isfile ( target_path ):
            assert 'cdo_path' in HOST_AVAIL [ self.B4J_host ] and \
                'login_name' in  HOST_AVAIL [ self.B4J_host ],\
                "Configuration for host %s incomplete." %  HOST_AVAIL [ self.B4J_host ]
            
            target_fn = os.path.basename ( target_path )
            ofnt     = orig_fn_tmpl.replace ( '%(year)s', '${y}' )

            cmd = """
            echo 'on server' $(hostname);
            for y in `seq %d 1 %d`;do
            ARRAY[$i]=%s/%s;
            echo ${ARRAY[i]}
            i=$i+1;
            done;
            export PATH=%s:$PATH;
            rm -f %s;
            echo %s %s %s %s
            cdo -f %s -t %s select,code=%s ${ARRAY[*]} %s;
            """ % (
                sy, fy, 
                orig_dirpath, ofnt,
                HOST_AVAIL [ self.B4J_host ][ 'cdo_path' ],
                target_fn, 
                suffix, codes_path, selcodes, target_fn,
                suffix, codes_path, selcodes, target_fn )
            
            if from_remote:
                self.B4J_obj.logger.debug ( "Before executing command '%s' on %s." % (
                        cmd, HOST_AVAIL [ self.B4J_host ][ 'login_name' ]))

                print "\n\n clear:\n\n %s \n\n" % cmd

                self.remotehost_processing ( cmd )
                self.B4J_obj.logger.debug ( "Executed command successfully." )
                self.remotehost_fetching ( target_fn, target_path )
                self.B4J_obj.logger.debug ( "Copied file to local machine successfully." )
                self.remotehost_processing ( "rm -f %s" % target_fn )
                self.B4J_obj.logger.debug ( "Deleted remote file successfully." )
            else:
                cmd = "cd %s" % ( self.B4J_obj.interpolated_selected_path ) + cmd

                self.B4J_obj.logger.debug ( "Before executing command '%s' on %s." % (
                        cmd, HOST_AVAIL [ self.B4J_host ][ 'login_name' ]))

                print "\n\n clear:\n\n %s \n\n" % cmd

                sp = subprocess.Popen( cmd, shell=True,
                                       stdout=subprocess.PIPE,  stderr=subprocess.PIPE )
                stdout, stderr = sp.communicate()
                if sp.returncode:
                    print "stdout: ", stdout
                    print "stderr: ", stderr
                    print "returncode: %d" % sp.returncode
                    sys.exit ( sp.returncode )
        else:
            pass
            
        self.B4J_obj.logger.debug ( "Exiting." )


    def fetch_archived_streams ( self, fetchfile_list ): 
        from bench import bench4jsbach_scripts
        
        self.B4J_obj.logger.debug ( "Starting" )

        fetchfile_list = list ( set ( fetchfile_list ))
        workingdir = os.getcwd ()
        os.chdir ( self.B4J_obj.downloads_path )
        exe = os.path.join ( bench4jsbach_scripts, "b4j_fetch_archived_file.ksh" )
        assert os.path.isfile ( exe ) and os.access ( exe, os.X_OK ),\
            "File %s is not executable." % exe
        self.B4J_obj.logger.debug ( """Fetching outdata and codes files from archive 
for %s into %s. This might take a while ...""" % (
            self.B4J_name_orig, self.B4J_obj.downloads_path ))
        
        for origin_dirpath, origin_filename in fetchfile_list:
            sp = subprocess.Popen( [ exe + " " + origin_dirpath + " " + origin_filename ],
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE )
            stdout, stderr = sp.communicate()
            if sp.returncode:
                print "stdout: ", stdout
                print "stderr: ", stderr
                print "returncode: %d" % sp.returncode
                sys.exit ( sp.returncode )
            else:
                print stdout
                
            self.B4J_obj.logger.debug ( "Fetched outdata from archive for %s successfully into %s." % (
                self.B4J_name_orig, self.B4J_obj.downloads_path ))

            # assuming tarball containing codesfiles in work/b4jname/downloads
            codes_tar_fn = self.codesfiles_dict [ 'tar_fn_tmpl' ] % { "expname" : self.B4J_name_orig }
            if origin_filename == codes_tar_fn:
                for cmd in [[ "tar -xf %s --wildcards --no-anchored --strip=1 '*.%s'" %
                              ( codes_tar_fn, self.codesfiles_dict [ "suffix" ])],
                            [ "rm -f %s" % codes_tar_fn ]]:
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
                                            ( codes_tar_fn, os.getcwd ()))

        os.chdir ( workingdir )

        self.B4J_obj.logger.debug ( "Exiting" )
            

    def fetch_outdata ( self ):
        """Fetch outdata.
"""
        from bench import bench4jsbach_scripts
        
        self.B4J_obj.logger.debug ( "Starting" )

        self.B4J_obj.logger.debug ( """Prepare fetching%s%s output data%s.""" %
                                    (( " remote" if self.from_remote else "local" ),
                                     ( " accumulated" if self.accumulated else " not self.accumulated" ),
                                     ( " from %s" % self.B4J_host if self.from_remote else "" )))
  
        if self.from_archive: fetchfile_list = []
        
        for fieldname, fielddict in self.B4J_obj.diagnostic_fields.iteritems ():
            self.B4J_obj.logger.debug ( "Fetch data for field '%s'." % fieldname )

            assert 'fetch' in fielddict and \
                self.B4J_etype in fielddict [ 'fetch' ], \
                "No 'fetch'-section defined for %s and %s." % ( fieldname, self.B4J_etype )
            
            frequence_infix = self.outfiles_dict [ 'frequence_infixes' ][ fielddict [ "frequence" ]]\
                              if fielddict [ "frequence" ] == "monthly" else "" 
            endyearstring = self.outfiles_dict [ 'endyear_tmpl' ] % { 'endyear': self.B4J_endyear }\
                              if self.accumulated else ""
                
            for fileid, filedict in fielddict [ 'fetch' ][ self.B4J_etype ].iteritems ():
                

                assert "submodel" in filedict and \
                       "stream" in filedict,\
                       "No '%s'- or '%s'-section defined for %s." % \
                       ( "submodel", "stream", fileid )
                self.B4J_obj.logger.debug ( "Fetch data for file '%s'." % fileid )

                self.B4J_index [ fileid ] = self.B4J_obj.path_template % { 
                    'path'           : self.B4J_obj.interpolated_selected_path,
                    'classification' : fileid,
                    'expname'        : self.B4J_name,
                    'startdate'      : self.B4J_startyear,
                    'enddate'      : self.B4J_endyear,
                    'suffix'         : fielddict [ 'suffix' ]}

                if os.path.isfile ( self.B4J_index [ fileid ]) and \
                   not self.B4J_obj.fetchremote:
                    continue
                
                orig_dirpath = os.path.join ( self.outdata_dirpath,
                                              filedict [ "submodel" ])

                # fetching outdata, not self.accumulated
                if not self.accumulated:
                    
                    orig_fn_tmpl = self.outfiles_dict [ 'fn_tmpl' ]  % {
                        "expname"             : self.B4J_name_orig,
                        "submodel"            : filedict [ "submodel" ],
                        "stream"              : filedict [ "stream" ],
                        "cond_frequenceinfix" : frequence_infix,
                        "startyear"           : "%(year)s",
                        "cond_endyearstr"   : "",
                        "suffix"              : self.outfiles_dict [ 'suffix' ]}
                    
                    codesfile_path = os.path.join ( self.codesfiles_dirpath,
                                                    self.codesfiles_dict [ "fn_tmpl" ] % {
                                                        "expname" : self.B4J_name_orig,
                                                        "submodel": filedict [ "submodel" ],
                                                        "stream"  : filedict [ "stream" ],
                                                        "suffix"  : self.codesfiles_dict [ 'suffix' ]})
                    
                    self.accumulate_streams ( orig_dirpath,
                                              orig_fn_tmpl,
                                              self.B4J_startyear,
                                              self.B4J_endyear,
                                              filedict [ 'code' ],
                                              fielddict [ 'suffix' ],
                                              codesfile_path,
                                              self.B4J_index [ fileid ],
                                              self.from_remote )
                # finished fetching outdata, not self.accumulated
                else:
                    orig_filename = self.outfiles_dict [ "fn_tmpl" ]  % {
                        "expname"             : self.B4J_name_orig,
                        "submodel"            : filedict [ "submodel" ],
                        "stream"              : filedict [ "stream" ],
                        "cond_frequenceinfix" : frequence_infix,
                        "startyear"           : self.B4J_initialyear,
                        "cond_endyearstr"     : endyearstring,
                        "suffix"              : self.outfiles_dict [ 'suffix' ]}

                    codes_filename = self.codesfiles_dict [  "fn_tmpl" ] % {
                        "expname"  : self.B4J_name_orig,
                        "submodel" : filedict [ "submodel" ],
                        "stream"   : filedict [ "stream" ],
                        "suffix"   : self.codesfiles_dict [ "suffix" ]}

                    if self.from_archive:
                        self.B4J_index [ fileid + STREAM_SUFFIX ] = os.path.join ( self.B4J_obj.downloads_path,
                                                                                   orig_filename )
                        if not os.path.isfile ( self.B4J_index [ fileid + STREAM_SUFFIX ]):
                            fetchfile_list.append (( orig_dirpath, orig_filename ))
                            
                        self.B4J_index [ fileid + CODES_SUFFIX ] = os.path.join ( self.B4J_obj.downloads_path,
                                                                                  codes_filename )
                        if not os.path.isfile ( self.B4J_index [ fileid + CODES_SUFFIX ]):
                            fetchfile_list.append (( self.dirpath,
                                                     self.codesfiles_dict [ 'tar_fn_tmpl' ] % \
                                                     { "expname" : self.B4J_name_orig }))
                    else:
                        self.B4J_index [ fileid + STREAM_SUFFIX ] = os.path.join ( orig_dirpath,
                                                                                   orig_filename )

        self.land_sea_mask  = self.B4J_obj.lsmpath_template % { 
            'path'           : self.B4J_obj.interpolated_selected_path,
            'classification' : "land_sea_mask",
            'expname'        : self.B4J_name,
            'suffix'         : fielddict [ 'suffix' ]}
        
        if os.path.isfile ( self.land_sea_mask ) and \
           not self.B4J_obj.fetchremote:
            return
        
        # TODO : CHECK/IMPLEMENT CASE FOR ACCUMULATED AND NOT FROM_ARCHIVE
        if self.from_archive and fetchfile_list:
            self.fetch_archived_streams ( fetchfile_list )
           
        if self.accumulated:
            cdo = Cdo ()
            
            for fieldname, fielddict in self.B4J_obj.diagnostic_fields.iteritems ():
                assert 'fetch' in fielddict, "No 'fetch'-section defined for %s." % fieldname
                yearlist = range ( int ( self.B4J_startyear ), int ( self.B4J_endyear ) + 1 )
                for fileid, filedict in fielddict ['fetch' ].iteritems ():
                    self.B4J_index [ fileid ] = self.B4J_obj.path_template % { 
                        'path'           : self.B4J_obj.interpolated_selected_path,
                        'classification' : fileid,
                        'expname'        : self.B4J_name,
                        'startdate'      : self.B4J_startyear,
                        'enddate'        : self.B4J_endyear,
                        'suffix'         : fielddict [ 'suffix' ]}
                    assert self.B4J_etype in filedict, "No '%s'-section defined for %s." % fileid
                    assert os.path.isfile ( self.B4J_index [ fileid + CODES_SUFFIX ])
                    try:
                        os.remove ( self.B4J_index [ fileid ])
                    except OSError:
                        pass
                    cdo.selcode (  filedict [ 'code' ],
                                   input = " -selyear," + ','.join ( repr ( y ) for y in yearlist ) +\
                                   " " + self.B4J_index [ fileid + STREAM_SUFFIX ],
                                   output = self.B4J_index [ fileid ],
                                   options = ' -f %s' % fielddict [ 'suffix' ] +\
                                   ' -t %s' % self.B4J_index [ fileid + CODES_SUFFIX ])

        
        # fetching land-sea-mask
        lsm_path = self.pool_dict [ "lsm_path_tmpl" ] % { "atm_res" : self.B4J_atmres,
                                                          "oce_res" : self.B4J_oceres }
        slf_varname = "slf"
        self.B4J_obj.logger.debug ( "Fetching land_sea_mask '%s' from %s." % ( slf_varname, lsm_path ))
        b4j_domain = getdomain ()

        if b4j_domain in self.pool_dict [ 'domains' ]:
            assert os.path.isfile ( lsm_path ), "No file found for %s." % lsm_path
            cdo = Cdo ()
            try:
                os.remove ( self.land_sea_mask )
            except OSError:
                pass
            cdo.selname ( slf_varname, input = lsm_path, output = self.land_sea_mask )
        else:
            self.fetch_land_sea_mask ( lsm_path, self.land_sea_mask, slf_varname )

        self.B4J_obj.logger.debug ( "Fetched land sea mask '%s'." %  self.land_sea_mask )
        # finished fetching land-sea-mask
                                    
        self.B4J_obj.logger.debug ( "Exiting" )
        

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


class ICON_Experiment ( Experiment ):
    """Derived class specialized on simulations done for Reiners atm_amip ICON.
"""
    gridtype = "unstructured"
    suffix = 'nc'
    outfiles_dict = { "suffix"            : "nc",
                      "fn_tmpl"           : "%(exp_name)s%(submodel)s_%(year)s%(mon_int)s01T000000Z.nc" }
    # lsm_dict = { "fn" : "r2b4_slf_bin.nc",
    #              "var": "notsea",
    #              "bin": True }
    lsm_dict = { "fn"  : "bc_land_frac.nc",
                 "varn": "notsea" }
    
    def __init__ ( self,  b4j_obj, **kwargs ):

        Experiment.__init__ ( self, b4j_obj, **kwargs )
        self.B4J_etype = "ICON"
        
        
    def fetch_outdata ( self ):
        """Fetch outdata.
"""
        from bench import bench4jsbach_scripts
        
        self.B4J_obj.logger.debug ( "Starting" )
        
        for fieldname, fielddict in self.B4J_obj.diagnostic_fields.iteritems ():
            self.B4J_obj.logger.debug ( "Fetch data for field '%s'." % fieldname )

            assert 'fetch' in fielddict and \
                self.B4J_etype in fielddict [ 'fetch' ], \
                "No 'fetch'-section defined for %s and %s." % ( fieldname, self.B4J_etype )

            cdo = Cdo ()
                
            for fileid, filedict in fielddict [ 'fetch' ][ self.B4J_etype ].iteritems ():
                self.B4J_index [ fileid ] = self.B4J_obj.path_template % { 
                    'path'           : self.B4J_obj.interpolated_selected_path,
                    'classification' : fileid,
                    'expname'        : self.B4J_name,
                    'startdate'      : self.B4J_startyear,
                    'enddate'        : self.B4J_endyear,
                    'suffix'         : fielddict [ 'suffix' ]}

                #if not self.B4J_obj.fetchremote: continue
                if os.path.isfile ( self.B4J_index [ fileid ]): continue


                orig_dirpath =  self.B4J_path
                fn_tmpl = "%(exp_name)s%(submodel)s_%(year)s%(mon_int)s01T000000Z.nc"
                orig_fn_tmpl = self.outfiles_dict [ "fn_tmpl" ] % { "exp_name" : self.B4J_name_orig,
                                                                    "submodel" : "_%s" % filedict [ 'submodel' ],
                                                                    "year"     : "${x}",
                                                                    "mon_int"  : "${y}" }
                suffix = 'nc'
                
                target_path = self.B4J_index [ fileid ]
                target_fn = os.path.basename ( target_path )

                cmd = """
                echo 'on server' $(hostname);
                echo $(pwd -P);
                logstring=$(date +%%Y%%m%%d%%H%%M)
                tempfn=temp_${logstring};
                echo ${tempfn};
                for x in `seq %s 1 %s`;do
                for y in {01..12};do
                ARRAY[$i]=%s/%s
                echo ${ARRAY[i]};
                i=$i+1;
                done;
                done;
                export PATH=%s:$PATH;
                echo ${PATH};
                echo %s %s %s;
                rm -f %s;
                cdo -P 8 -f %s select,name=%s ${ARRAY[*]} ${tempfn};
                cdo monmean ${tempfn} %s;
                rm -r ${tempfn}
                """ % ( self.B4J_startyear, self.B4J_endyear,
                        self.B4J_path, orig_fn_tmpl,
                        HOST_AVAIL [ self.B4J_host ][ 'cdo_path' ],
                        suffix, filedict [ 'name' ],target_fn,
                        target_fn,
                        suffix, filedict [ 'name' ], 
                        target_fn )
            
                if self.from_remote:
                    self.B4J_obj.logger.debug ( "Before executing command '%s' on %s." % (
                            cmd, HOST_AVAIL [ self.B4J_host ][ 'login_name' ]))

                    print "\n\n clear:\n\n %s \n\n" % cmd
                
                    self.remotehost_processing ( cmd )
                    self.B4J_obj.logger.debug ( "Executed command successfully." )
                    self.remotehost_fetching ( target_fn, target_path )
                    self.B4J_obj.logger.debug ( "Copied file to local machine successfully." )
                    self.remotehost_processing ( "rm -f %s" % target_fn )
                    self.B4J_obj.logger.debug ( "Deleted remote file successfully." )
                else:
                    cmd = "cd %s" % ( self.B4J_obj.interpolated_selected_path ) + cmd

                    self.B4J_obj.logger.debug ( "Before executing command '%s' on %s." % (
                            cmd, HOST_AVAIL [ self.B4J_host ][ 'login_name' ]))

                    print "\n\n clear:\n\n %s \n\n" % cmd

                    sp = subprocess.Popen( cmd, shell=True,
                                           stdout=subprocess.PIPE,  stderr=subprocess.PIPE )
                    stdout, stderr = sp.communicate()
                    if sp.returncode:
                        print "stdout: ", stdout
                        print "stderr: ", stderr
                        print "returncode: %d" % sp.returncode
                        sys.exit ( sp.returncode )

        
        self.land_sea_mask  = self.B4J_obj.lsmpath_template % { 
            'path'           : self.B4J_obj.interpolated_selected_path,
            'classification' : "land_sea_mask",
            'expname'        : self.B4J_name_orig,
            'suffix'         : "nc" }

        if not os.path.isfile ( self.land_sea_mask ):
            # fetching land-sea-mask
            lsm_path = """%(path)s/%(fn)s""" % {
                'path' : self.B4J_path,
                'fn': self.lsm_dict [ 'fn' ]}

            self.B4J_obj.logger.debug ( "Fetching land_sea_mask '%s' from %s." % ( self.lsm_dict [ "varn" ], 
                                                                                   lsm_path ))
            if not self.from_remote:
                assert os.path.isfile ( lsm_path ), "No file found for %s." % lsm_path
                cdo = Cdo ()
                cdo.gtc ( 0,
                          input = " -select,name=%s %s " % ( self.lsm_dict [ "varn" ],
                                                             lsm_path ),
                          output = self.land_sea_mask )
            else:
                print "Fetching land sea mask: remote case not implemented jet. EXITING."

                
                cmd = """
                echo 'on server' $(hostname);
                echo $(pwd -P);
                export PATH=%s:$PATH;
                echo ${PATH};
                echo %s %s %s %s;
                rm -f %s;
                cdo  -f %s -gtc,0 -select,name=%s %s %s;
                """ % ( HOST_AVAIL [ self.B4J_host ][ 'cdo_path' ],
                        fielddict [ 'suffix' ], self.lsm_dict [ "varn" ],
                        lsm_path, os.path.basename ( self.land_sea_mask ),
                        os.path.basename ( self.land_sea_mask ),
                        fielddict [ 'suffix' ], self.lsm_dict [ "varn" ],
                        lsm_path, os.path.basename ( self.land_sea_mask ))
                
                self.B4J_obj.logger.debug ( "Before executing command '%s' on %s." % (
                    cmd, HOST_AVAIL [ self.B4J_host ][ 'login_name' ]))

                print "\n\n clear:\n\n %s \n\n" % cmd

                self.remotehost_processing ( cmd )
                self.B4J_obj.logger.debug ( "Executed command successfully." )
                self.remotehost_fetching (  os.path.basename ( self.land_sea_mask ),
                                            self.land_sea_mask )
                self.B4J_obj.logger.debug ( "Copied file to local machine successfully." )
                self.remotehost_processing ( "rm -f %s" % os.path.basename ( self.land_sea_mask ))
                self.B4J_obj.logger.debug ( "Deleted remote file successfully." )

            self.check_output ( self.land_sea_mask )
            self.B4J_obj.logger.debug ( "Fetched land sea mask '%s'." %  self.land_sea_mask )
            # finished fetching land-sea-mask
                
        self.B4J_obj.logger.debug ( "Exiting" )


class CMIP5_Experiment ( Experiment ):
    """Derived class specialized on simulations done for CMIP5.
"""
    gridtype = "gaussian"
    suffix = 'nc'
    
    def __init__ ( self,  b4j_obj, **kwargs ):
        var_filename = kwargs.pop ( "var_filename", None )
        assert var_filename and os.path.isfile ( var_filename )

        Experiment.__init__ ( self, b4j_obj, **kwargs )
        self.B4J_etype = "CMIP5"
        
        nc_fid  = Dataset ( var_filename, 'r' )
        nc_attrs = nc_fid.ncattrs ()
        for nc_attr in nc_attrs:
            self.__dict__[ nc_attr ] = nc_fid.getncattr ( nc_attr )
        nc_fid.close ()


# status: only for archived experiment
class MPIESM1_Experiment ( Experiment ):
    """Derived class specialized on simulations with MPIESM1.
    """
    gridtype = "gaussian"
    
    codesfiles_dict = { "suffix"        : "codes",
                        "dirname"       : "log",
                        "fn_tmpl"       : "%(expname)s_%(submodel)s_%(stream)s.%(suffix)s",
                        "is_tared"      : False,
                        "tar_fn_tmpl"   : "%(expname)s_log.tar" }

    outfiles_dict = { "suffix"            : "grb",
                      "dirname"           : "outdata",
                      "fn_tmpl"           : """%(expname)s_%(submodel)s_%(stream)s_%(cond_frequenceinfix)s%(startyear)s%(cond_endyearstr)s.%(suffix)s""",
                      "frequence_infixes" : { "monthly" : "mm_" },
                      "endyear_tmpl"    : "-%(endyear)s" }
                      
    pool_dict = { "domains"       : [ "zmaw.de", "hpc.dkrz.de" ],
                  "lsm_path_tmpl" : "/pool/data/JSBACH/input/r0007/%(atm_res)s/jsbach_%(atm_res)s%(oce_res)s_11tiles_5layers_1850.nc" }

    def __init__ ( self, b4j_obj, **kwargs ):
        Experiment.__init__ ( self, b4j_obj, **kwargs )
        self.dirpath            = os.path.join ( self.B4J_path, self.B4J_name_orig )
        self.outdata_dirpath    = os.path.join ( self.dirpath, self.outfiles_dict [ "dirname" ])
        self.codesfiles_dirpath = os.path.join ( self.dirpath, self.codesfiles_dict [ 'dirname' ])
        self.B4J_etype          = "MPIESM1"
        self.B4J_atmres         = kwargs.get ( "B4J_atmres" )
        self.B4J_oceres         = kwargs.get ( "B4J_oceres" )

        self.codesfiles_dict [ "is_tared" ] = self.from_archive
        
        #print "MPIESM1_Experiment attributes:"
        #print self.__dict__       

                                                                                                                
if __name__ == "__main__":

    var_filename = "/scratch/mpi/mpiles/m300050/benchmarking/example_data/CMIP5/amip/rsds_Amon_MPI-ESM-LR_amip_r3i1p1_197901-200812.nc"

    assert os.path.isfile ( var_filename )
    
    e = CMIP5_Experiment ( var_filename = var_filename )

    print e.institution

    print e
    
    
