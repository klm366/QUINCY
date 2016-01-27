#! /usr/bin/env python

import sys,logging,textwrap

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


if __name__ == '__main__':
    message = "Lorem ipsum dolor sit amet, consectetur adipisicing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum."
    s = "0123456789"*100
    
    logger.debug ( "Starting" )
    logger.debug ( message )
    logger.error ( "This is a test." )
    logger.info  ( "This is a test." )
    logger.debug ( "Exiting" )



