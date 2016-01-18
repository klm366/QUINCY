#! /usr/bin/env python

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


    def __str__ ( self ):
        return self.__repr__ ()
