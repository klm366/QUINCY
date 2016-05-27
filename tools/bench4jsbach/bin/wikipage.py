#! /usr/bin/env python

from redmine import Redmine
import re, os, sys

roles                         = [ "ref", "can" ]
ref_table_starttoken          = '\n<notextile></notextile>\r\n\r\n'
ref_table_endtoken            = '\r\n\r\n<notextile></notextile>\r\n'
can_table_starttoken          = '\n<notextile></notextile><notextile></notextile>\r\n\r\n'
can_table_endtoken            = '\r\n\r\n<notextile></notextile><notextile></notextile>'
b4j_table_starttoken          = '\n<notextile></notextile><notextile></notextile><notextile></notextile>\r\n\r\n'
b4j_table_endtoken            = '\r\n\r\n<notextile></notextile><notextile></notextile><notextile></notextile>'

column_header_simulation_name = 'simulation_name'
column_header_data_storage    = 'data_storage'


class Bench4JSBach_Wikipage ( object ):
    codezmaw_url     = "https://code.zmaw.de"
    project_name     = "jsbach-bench"
    wikipage_name    = "Bench4jsbach_simulations"
    redmine          = None
    redmine_project  = None
    redmine_wikipage = None

    def __init__ ( self, redmine_user, redmine_key ):
        try:
            redmine = Redmine ( self.codezmaw_url, username = redmine_user,
                                key = redmine_key )
        except:
            print "ERROR: Could not connect to redmine %s for user %s with key %s." %\
                  ( self.codezmaw_url, redmine_user, redmine_key )
            sys.exit ( 1 )
        self.redmine              = redmine
        self.redmine_project      = redmine.project.get ( self.project_name )
        self.redmine_wikipage     = redmine.wiki_page.get ( self.wikipage_name,
                                                    project_id = self.redmine_project.id )


    def fetch_simulation_row ( self, experiment_role, experiment_key ):
        assert experiment_role in roles, "No role for experiment given."
        table_starttoken = ref_table_starttoken if experiment_role == "ref" else can_table_starttoken
        table_endtoken = ref_table_endtoken if experiment_role == "ref" else can_table_endtoken
        
        tmp = re.search ( r'.*%s(.*)%s.*' % ( table_starttoken, table_endtoken ),
                          self.redmine_wikipage.text, re.DOTALL )
        try:
            table_raw = tmp.group(1).strip()
        except:
            raise Exception (  'ERROR reading experiment table from wiki page.' )
    
        rows = [ v.strip ().strip ( '|' ) for v in table_raw.split ( '\n' )]
        
        try:
            headings =  [ v.split ( '.' ).pop ().strip ().encode ().replace ( ' ', '_' )
                          for v in rows [ 0 ].split ( '|' )]
        except:
            raise Exception ( 'ERROR reading headings from experiment table.' )
        assert column_header_simulation_name in headings and\
               column_header_data_storage in headings ,\
               "Table headers are not as exspected."

        dicts = []   
        if len ( rows ) > 1:
            for r in rows[1:]:
                values = [ v.strip ().encode ()  for v in r.split ( '|' )]
                expdict = dict ( zip ( headings, values ))
                dicts.append ( expdict )
        else:
            print ( 'No experiments found in experiment table.' )

        try:
            return  next (( d for d in dicts if d [ column_header_simulation_name ] ==
                            experiment_key ), {} )
        except:
            raise Exception ( 'ERROR searching table columns with key "%s".' %
                              column_header_simulation_name )


    def write_b4j_row ( self ):
        
        tmp = re.search ( r'.*%s(.*)%s.*' % ( b4j_table_starttoken, b4j_table_endtoken ),
                          self.redmine_wikipage.text, re.DOTALL )
        try:
            table_raw = tmp.group(1).strip()
        except:
            raise Exception (  'ERROR reading experiment table from wiki page.' )

        print table_raw
    
        rows = [ v.strip ().strip ( '|' ) for v in table_raw.split ( '\n' )]
        
        try:
            headings =  [ v.split ( '.' ).pop ().strip ().encode ().replace ( ' ', '_' )
                          for v in rows [ 0 ].split ( '|' )]
        except:
            raise Exception ( 'ERROR reading headings from experiment table.' )

        print headings


    def upload_file ( self, filepath ):

        assert filepath and os.path.isfile ( filepath )

        print 
        print self.redmine.__dict__
        print
        print self.redmine_wikipage.__dict__
        print
        print self.redmine_wikipage.manager.__dict__
        print
        
        #token = self.redmine_wikipage.upload ( filepath )
        #self.redmine.wiki_page.update ("Bench4jsbach_simulations", project_id=79, uploads=[{'path': filepath}])
        b = self.redmine_wikipage.uploads = [{'path': filepath}]
        a = self.redmine_wikipage.comments = "Some interesting comment."
        print a
        print b
        #self.redmine_wikipage.save ()

    def test ( self ):
        print self.redmine.project.all ().__dict__
    


if __name__ == '__main__':
    b4j_wikipage = Bench4JSBach_Wikipage ( 'deike.kleberg', '676af888885116cd94e9c36cd0a703826aa706e8' )
    print b4j_wikipage.fetch_simulation_row ( "ref", "vga0188" )
    print b4j_wikipage.fetch_simulation_row ( "can", "vga0186" )
    b4j_wikipage.write_b4j_row ()
    b4j_wikipage.upload_file ( "/scratch/mpi/mpiles/m300050/git/jsbach/tools/bench4jsbach/bin/testfile.txt" )
    b4j_wikipage.test ()
