#! /usr/bin/env python

from util import log_message_newline,logger as dkl_logger, print_none

class Experiment ( object ):
    id    = None
    name  = None
    etype = None
    path  = None
    
    def __init__ ( self, **kwargs ):
        dkl_logger.debug ( "Starting" )
        self.name     = kwargs.get ( "name" )
        self.etype    = kwargs.get ( "etype", None )
        self.path     = kwargs.get ( "path", None )
        self.buffer   = []
        dkl_logger.debug ( "Exiting" )
        

    def __repr__ ( self ):
        return "{} ( id = {}, name = {}, etype = {}, path = {}, buffer = {})". \
               format ( self.__class__.__name__,
                        print_none ( self.id ),
                        print_none ( self.name ),
                        print_none ( self.etype ),
                        print_none ( self.path ),
                        print_none ( self.buffer ))


    def __str__ ( self ):
        return self.__repr__ ()
