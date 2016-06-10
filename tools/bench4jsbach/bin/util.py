#! /usr/bin/env python

import os,sys,logging,textwrap
import subprocess


default_logger_name = 'DKL_LOGGER'
log_message_newline = "\n  #\t"

class WrappedFixedIndentingLog ( logging.Formatter ):
    def __init__ ( self, fmt = None, datefmt = None, width = 120,
                   indent = 2 ):
        super ( WrappedFixedIndentingLog, self ).\
              __init__ ( fmt = fmt, datefmt = datefmt )
        self.wrapper = textwrap.TextWrapper ( width = width,
                                              subsequent_indent = ' '*indent + '# ',
                                              replace_whitespace = False )

    def format ( self, record ):
        return self.wrapper.fill ( super ( WrappedFixedIndentingLog,
                                           self ).format ( record))


def setup_logger ( *args, **kwargs ):

    name     = kwargs.get ( 'name', default_logger_name )
    filename = kwargs.get ( 'filename', None )

    logger   = logging.getLogger ( name )
    logger.propagate = False
    
    if not len ( logger.handlers ):
        if filename:
            handler = logging.FileHandler ( filename )
        else:
            handler = logging.StreamHandler ( sys.stdout )
            
        handler.setFormatter ( WrappedFixedIndentingLog (
            """\n### %(name)s - %(levelname)s -  %(asctime)-8s - MODULE %(module)s - FUNCTION %(funcName)s - LINE %(lineno)d :
            
            %(message)s\n###""",
            "%Y/%m/%d - %H:%M:%S",
            indent=2 )) 

        logger.addHandler ( handler )
        logger.setLevel ( logging.DEBUG )
    return logger


logger = setup_logger ()


def str2bool ( s ) :
  return s.lower() in ( "true", "t", "1", "yes", "y" )

def print_none ( item ):
    return item if item is not None else "None"

def print_item_or_list_rec ( something, indent ):
    import collections

    nextindent = indent + 3
    
    baseind = " " *   indent
    nextind = " " *   nextindent
    joinind = " " * ( nextindent + 20 )
    sepind  = " " * ( nextindent + 20 )
    if not something:
        return "None"
    elif isinstance ( something, basestring ):
        return something
    elif not isinstance ( something, collections.Iterable ):
        return something
    elif isinstance ( something, list ):
        return "[\n%s]" % ( ",\n".join ( '{!s}{!s}'.format ( nextind, print_item_or_list_rec ( i, nextindent + 3 ))\
                                       for i in something ))
    elif isinstance ( something, dict ):
        return "{\n%s}" % ( ",\n".join ( '{!s}{!s}\t: {!s}'.format ( nextind, k, print_item_or_list_rec ( v, nextindent + 3 )) \
                                       for k,v in something.iteritems ()))
    else:
        print "It's somthing different."
        


def pprint_class_member_methods ( modulename, classname ):
    """python -c "import sys;sys.path.append ( '/scratch/mpi/mpiles/m300050/bench4jsbach/bin' );import util; util.pprint_class_member_methods ( 'bench', 'Bench_Object' )"
"""
    import inspect
    import pprint
    try:
        module = __import__ ( modulename )
    except Exception,e:
        raise ImportError ( "%s. \nNo module named '%s'." % ( str ( e ), modulename ))
    try:
        class_ = getattr ( module, classname )
    except Exception,e:
        raise AttributeError ( "%s. \nNo class named '%s' in '%s'." % ( str ( e ),
                                                                        classname,
                                                                        modulename ))
    mem = inspect.getmembers ( class_, predicate=inspect.ismethod)
    pprint.pprint ([ i for i,j in mem])


def touch ( filename ):
    if os.path.exists ( filename ):
        os.utime ( filename, None )
    else:
        open ( filename, 'a' ).close ()



def subprocessing ( cmd ):
    sp = subprocess.Popen ( cmd, shell=True,
                            stdout=subprocess.PIPE,  stderr=subprocess.PIPE )
    stdout, stderr = sp.communicate()
    if sp.returncode:
        print "stdout: ", stdout
        print "stderr: ", stderr
        print "returncode: %d" % sp.returncode
        sys.exit ( sp.returncode )
    else:
        print stdout


def getdomain ():
    import socket
    
    name, hlist, alist = socket.gethostbyaddr (
        socket.gethostbyname ( socket.gethostname ())) 
    dot = '.'
    domainname = None

    for h in hlist + [ name ]:
        if h.count ( dot ) >= 2:
            d = h.split ( dot, 1 )[ 1 ]
            if not domainname:
                domainname = d
            else:
                if not d == domainname:
                    raise Exception ( "Got different domainnames for host '%s'." % name )

    return domainname



def gethostname ():
    import socket
    name, hlist, alist = socket.gethostbyaddr ( 
        socket.gethostbyname ( socket.gethostname ())) 
    return name.split ( '.' )[ 0 ]


if __name__ == '__main__':
    message = "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
    s = "0123456789"*100
    
    logger.debug ( "Starting" )
    logger.debug ( message )
    logger.error ( "This is a test." )
    logger.info  ( "This is a test." )
    logger.debug ( "Exiting" )



