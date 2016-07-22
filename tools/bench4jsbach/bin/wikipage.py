#! /usr/bin/env python

from redmine import Redmine
import re, os, sys

TABLE_SEP                     = "|"
TABLE_NCOL                    = 10
TABLE_SHIFT                   = 1
TABLE_MARKER                  = "_{background:lightgrey}."
TABLE_ERROR_MESSAGE           = "Wikipage table for experiments is not readable."

roles                         = [ "ref", "can" ]
ref_table_starttoken          = '\r\n<notextile></notextile>\r\n\r\n'
ref_table_endtoken            = '\r\n\r\n<notextile></notextile>\r\n'
can_table_starttoken          = '\r\n<notextile></notextile><notextile></notextile>\r\n\r\n'
can_table_endtoken            = '\r\n\r\n<notextile></notextile><notextile></notextile>\r\n'
b4j_table_starttoken          = '\r\n<notextile></notextile><notextile></notextile><notextile></notextile>\r\n\r\n'
b4j_table_endtoken            = '\r\n\r\n<notextile></notextile><notextile></notextile><notextile></notextile>\r\n'

column_header_simulation_name   = 'simulation_name'
column_header_data_storage      = 'data_storage'
column_header_simulation_period = 'simulation_period'
detail_initialyear              = 'initialyear'
detail_finalyear                = 'finalyear'


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

    def fetch_simulation_row_split_period ( self, experiment_role, experiment_key  ):
        def is_number ( s ):
            try:
                int ( s )
                return True
            except ValueError:
                return False
            
        tablerow = self.fetch_simulation_row ( experiment_role, experiment_key )
        if tablerow:
            simulation_period = tablerow.get ( 'simulation_period' )
            if simulation_period:
                years = filter ( None, simulation_period.split ( '-' ))
                assert len ( years ) == 2 and \
                       all ( is_number ( s ) for s in years ) and \
                       years [ 0 ] < years [ 1 ]
            tablerow [ detail_initialyear ] = years [ 0 ]
            tablerow [ detail_finalyear ]   = years [ 1 ]
            tablerow.pop ( column_header_simulation_period, None )
            
        return tablerow            


    def fetch_simulation_row ( self, experiment_role, experiment_key ):
        assert experiment_role in roles, "No role for experiment given."
        table_starttoken = ref_table_starttoken if experiment_role == "ref" else can_table_starttoken
        table_endtoken = ref_table_endtoken if experiment_role == "ref" else can_table_endtoken
        
        sre_match_obj = re.search ( r'.*%s(.*)%s.*' % ( table_starttoken, table_endtoken ),
                          self.redmine_wikipage.text, re.DOTALL )       
        try:
            textblock = sre_match_obj.group(1).strip()
        except:
            raise Exception (  'ERROR reading experiment table from wiki page.' )
        
        cells = textblock.split ( TABLE_SEP )
        ncells = len ( cells )
        # SHIFT are induced due to the seperators in front of the rows in the tables on the wikipage.
        rawrowlist = [ cells [ i + TABLE_SHIFT : i + TABLE_NCOL + TABLE_SHIFT ]  \
                       for i in range ( 0, ncells, TABLE_NCOL + TABLE_SHIFT ) \
                       if cells [ i + 1 : i + TABLE_NCOL + TABLE_SHIFT ]]

        # CONFIDENCE TEST BLOCK
        assert all ( map ( lambda c: len ( c ) == TABLE_NCOL, rawrowlist )), TABLE_ERROR_MESSAGE
        assert all ( map ( lambda c: c.startswith ( TABLE_MARKER ), rawrowlist [ 0 ])), TABLE_ERROR_MESSAGE
        for row in rawrowlist [ 1: ]:
            assert all ( map ( lambda c: not c.startswith ( TABLE_MARKER ), row )), TABLE_ERROR_MESSAGE

        headings =   [ v.split ( '.' ).pop ().strip ().encode ().replace ( ' ', '_' )
                          for v in rawrowlist [ 0 ]]
        values = [[ v.strip ().encode ().replace ( "\r\n", "," ) for v in row ] for row in rawrowlist [ 1: ]]
        table = dict (( row [ 0 ], dict ( zip ( headings, row ))) for row in values ) 

        return table [ experiment_key ] if experiment_key in table else None


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
    print b4j_wikipage.fetch_simulation_row ( "can", "slo0022thuO" )
    # b4j_wikipage.write_b4j_row ()
    # b4j_wikipage.upload_file ( "/scratch/mpi/mpiles/m300050/git/jsbach/tools/bench4jsbach/bin/testfile.txt" )
    # b4j_wikipage.test ()
