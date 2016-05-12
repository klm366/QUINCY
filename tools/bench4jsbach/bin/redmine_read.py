#! /usr/bin/env python

from redmine import Redmine
import re, os

redmine_url                   = 'https://code.zmaw.de'
redmine_user                  = 'deike.kleberg'
redmine_key                   = '676af888885116cd94e9c36cd0a703826aa706e8'
redmine_project               = 'jsbach-bench'
redmine_wikipage              = 'Bench4jsbach_simulations'

roles                         = [ "ref", "can" ]
ref_table_starttoken          = '\n<notextile></notextile>\r\n\r\n'
ref_table_endtoken            = '\r\n\r\n<notextile></notextile>\r\n'
can_table_starttoken          = '\n<notextile></notextile><notextile></notextile>\r\n\r\n'
can_table_endtoken            = '\r\n\r\n<notextile></notextile><notextile></notextile>'

column_header_simulation_name = 'simulation_name'
column_header_data_storage    = 'data_storage'

this_key                      = 'vga0187'
this_role                     = 'can'


class Wikitable_experiment_description ( object ):

    def __init__ ( self, redmine_url, redmine_user, redmine_key,
                   redmine_project, redmine_wikipage,
                   experiment_role, experiment_key ):
        self.__dict__.update ( self.redmine_wiki_table_read ( redmine_url,
                                                              redmine_user,
                                                              redmine_key,
                                                              redmine_project,
                                                              redmine_wikipage,
                                                              experiment_role,
                                                              experiment_key ))


    def redmine_wiki_table_read ( self, redmine_url, redmine_user,
                                  redmine_key, redmine_project, redmine_wikipage,
                                  experiment_role, experiment_key ):

        assert experiment_role in roles, "No role for experiment given."
        table_starttoken = ref_table_starttoken if experiment_role == "ref" else can_table_starttoken
        table_endtoken = ref_table_endtoken if experiment_role == "ref" else can_table_endtoken

        redmine   = Redmine ( redmine_url, username = redmine_user,
                              key = redmine_key )
        project   = redmine.project.get ( redmine_project )
        wikipage = redmine.wiki_page.get ( redmine_wikipage,
                                            project_id = project.id )
        
        tmp = re.search ( r'.*%s(.*)%s.*' % ( table_starttoken, table_endtoken ),
                          wikipage.text, re.DOTALL )
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

    
    def check_for_local_data ( self ):
        try:
            path = getattr ( self, column_header_data_storage )
        except:
            raise Exception ( 'ERROR trying to get object attribute "%s".' %
                              column_header_data_storage )

        return os.path.isdir ( path )



if __name__ == '__main__':

    name = this_key
    role = this_role
    
    red = Wikitable_experiment_description ( redmine_url, redmine_user, redmine_key,
                   redmine_project, redmine_wikipage, role, name )
    assert red.__dict__, 'No experiment "%s" found.' % name

    print '\nRedmine reference experiment description:  \n%s\n' % \
          ''.join('{!s} = {!r},\n'.\
                  format ( key, val) for ( key, val ) in red.__dict__.iteritems ())
    print red.data_storage
    print 'Data directory of', this_key, 'exists: ',  red.check_for_local_data ()
    
