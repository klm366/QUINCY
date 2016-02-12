#! /usr/bin/env python

import os, subprocess, shutil, re
import json
import copy

from util import log_message_newline

doc_template = """
\\documentclass[]{article}
\\usepackage[english]{babel}
\\usepackage{graphicx}
\\usepackage{color}

\\usepackage{booktabs,tabularx,ragged2e,makecell}
\\newcolumntype{Y}{>{\\raggedright\\arraybackslash}X}
\\newcolumntype{L}{>{\\arraybackslash}X}
\
\\usepackage{titling}
\\usepackage{lipsum}
\\renewcommand{\\familydefault}{\\sfdefault}
\\newcommand{\graphicpath}[1]{%(graphicpath)s/#1}

\\setlength{\parindent}{0pt}


\\begin{document}

\\title{%(title)s}
\\author{%(author)s, %(institute)s, %(department)s}
\maketitle
\\thispagestyle{empty}

\section{Benchmark description}
%(abstract)s

%(body)s

\section{Summary}
%(epilog)s \lipsum[1]

\end{document}
"""
section_template = """\section{%(section)s}
%(body)s"""

subsection_template = """\subsection{%(subsection)s}"""
subsubsection_template = """\subsubsection{%(subsubsection)s}"""

report_dirname = "report"
source_dirname = "sources"


analyses_var_item = []

analyses_var_test_item = { 'name'     : None,
                           'infiles'  : [],
                           'outfiles' : [],
                           'plots'    : [],
                           'teXts'    : [],
                           'matrix'   : False,
                           'vector'   : False,
                           'cans'     : [],
                           'title'    : None }

analyses_var_test_can_item = { 'name'     : None,
                               'infiles'  : [],
                               'outfiles' : [],
                               'plots'    : []}

def init_analyses_var_test_can_item ( name = None ):
    ret =  copy.deepcopy ( analyses_var_test_can_item )
    ret [ 'name' ] = name 
    return ret
                   


def includegraphics_template ( path, width = '1', hspace = '0' ):
    
    ret = "\\includegraphics[width=%s\\textwidth]{%s}" % ( width, path )
    if hspace:
        ret = """

\\hspace*{%scm}%s

""" % ( hspace, ret )
    else:
        ret = """ 

\\begin{centering}
%s
\\end{centering}

    """ % ret

    return ret
    

def teXify ( string ):
    assert isinstance ( string, str )
    
    ret = string
    ret = string.replace ( '_', '\_' )
    ret = re.sub ( "\"(.*?)\"",r"``\1''",ret )
    ret = re.sub ( "'(.*?)'",r"`\1'",ret )

    return ret
    


class Report ( object ):
    body = ""
    graphic = None
    analyses = {}

    def __init__ ( self ):
        self.doc_template = doc_template

    def write_description ( self ):

        analysis_table = [ ]

        description = """
In this benchmark report, the candidate experiments
\\begin{itemize}
%s
\\end{itemize}
are compared with the reference experiments
\\begin{itemize}
%s
\\end{itemize} 
\\ \\\\ \\\\

The variables tested and the analyses and plots made are listed in the 
table below.
\\ \\\\

\\begin{tabularx}{\\linewidth}{lLL}
\\toprule
\\textbf{variable} & \\textbf{analyis method} & \\textbf{shortname} \\\\
\\addlinespace
\\cmidrule(r){1-3}
%s
\\\\ \\bottomrule
\\end{tabularx}

""" %(
    ''.join ( '\\item {!s} ($C_{{{!s}}}$)'.format ( teXify ( can.B4J_name ),
                                                    str ( can.B4J_canno )) \
              if can.B4J_iscan else '' \
                for can in self.experiments ),
    ''.join ( '\\item {!s} ($R_{{{!s}}}$)'.format ( teXify ( ref.B4J_name ),
                                                   str ( ref.B4J_refno )) \
              if ref.B4J_isref else '' \
                for ref in self.experiments ),
    '\\\\\n'.join ( '%s & %s & %s' %\
                    (( teXify ( varname ),
                      teXify ( str ( analysesitemdict [ 'title' ])),
                      teXify ( str ( analysesitemdict [ 'name' ]))) \
                    if self.analyses [ varname ].index ( analysesitemdict ) == 0 \
                    else ( '',
                           teXify ( str ( analysesitemdict [ 'title' ])),
                           teXify ( str ( analysesitemdict [ 'name' ])))) \
                    for varname in self.analyses \
                    for analysesitemdict in self.analyses [ varname ] ))

        description += subsection_template % { 'subsection': 'Candidate experiments to benchmark' }
    
        for can in self.experiments:
            if not can.B4J_iscan: continue
            description +=  can.teX_repr ()

        description += subsection_template % { 'subsection': 'Reference experiments' }
        for ref in self.experiments:
            if not ref.B4J_isref: continue
            if ref.B4J_name == can.B4J_name: continue
            description +=  ref.teX_repr ()

        return description 
    

    def albedo_section ( self ):
        self.logger.debug ( "Starting." )
        varname = "albedo"
        albedo_body = "" 
        self.logger.debug ( "self.analyses buffer: " + log_message_newline +\
                            json.dumps ( self.analyses ))

        print self.analyses

        for analyse in self.analyses [ varname ]:
            
            albedo_body += subsection_template % { 'subsection': teXify ( analyse [ 'title' ])}
            albedo_body +=  analyse [ 'teXts' ][ 0 ]
            if analyse [ 'matrix' ]:
                for can in self.experiments:
                    if not can.B4J_iscan: continue
                    can_name = teXify ( can.B4J_name )
                    title = "Candidate experiment {\\tt %s} ($\mathbf{C_{%d}}$)" %\
                            ( can_name, can.B4J_canno ) 
                    albedo_body += subsubsection_template % { 'subsubsection': title }
                    analyse_can = next (( item for item in analyse [ 'cans' ] \
                                          if item [ 'name' ] == can.B4J_name ), None )
                    if len ( analyse_can [ 'plots' ]):
                        no = 0
                        for ref in self.experiments:
                            if not ref.B4J_isref: continue
                            if ref.B4J_name == can.B4J_name: continue
                            plot = analyse_can [ 'plots' ][ no ]
                            albedo_body += includegraphics_template (  plot,
                                                                       width      = '1.5',
                                                                       hspace     = '-3.3' )
                            no += 1
                            
            elif analyse [ 'vector' ]:
                for can in self.experiments:
                    if not can.B4J_iscan: continue
                    can_name = teXify ( can.B4J_name )
                    title = "Candidate experiment {\\tt %s} ($\mathbf{C_{%d}}$)" %\
                            ( can_name, can.B4J_canno ) 
                    albedo_body += subsubsection_template % { 'subsubsection': title }
                    analyse_can = next (( item for item in analyse [ 'cans' ] \
                                          if item [ 'name' ] == can.B4J_name ), None )
                    if len ( analyse_can [ 'plots' ]):
                        plot = analyse_can [ 'plots' ][ 0 ]
                        albedo_body += includegraphics_template ( plot,
                                                                  width    = '1',
                                                                  hspace   = '0' )


        section_intro = """The the output data of the analyses can be found in \\\\ 
\\bigskip

  \\-\ \\hspace{.5cm} \\texttt{%s},\\\\
\\medskip
the plots\\\\
\\medskip
   \\-\ \\hspace{.5cm}  \\texttt{%s}.
\\medskip""" % ( teXify ( self.data_path ), teXify ( self.graphic_path ))

        section = section_template % { 'section' : "%s: analyses and benchmarks for variable" % \
                                       varname.title (),
                                       'body' : section_intro + albedo_body }
        self.body += section
    

        
    def write_report ( self ):
        self.logger.debug ( "Starting." )
        source_path  = os.path.join ( self.report_path, "sources" )
        tmp_path     = os.path.join ( self.report_path, "tmp" )
        if not os.path.exists ( source_path ):
            os.makedirs ( source_path )
        if not os.path.exists ( tmp_path ):
            os.makedirs ( tmp_path )
        os.chdir ( source_path )

        title_second = "%s" % ''.join ( '\\\\{!s} ({!s})'. format ( teXify ( can.B4J_name ),  teXify ( can.B4J_etype )) \
                if can.B4J_iscan else '' \
                for can in self.experiments )

        page = doc_template % { 'title'      : """
        {\\normalsize \\textbf{Bench4JSBach}} {\\Large Report for Experiments} 
        {\\LARGE %s}
        """ % title_second, 
                                'author'     : "%s" % self.config [ 'author' ],
                                'institute'  : "Max-Planck-Institute for Meteorology",
                                'department' : "Land Department",
                                'graphicpath': self.graphic_path,
                                'abstract'   : self.write_description (),
                                'body'       : self.body,
                                'epilog'     : """This chapter provides a brief overview of the analyses and 
                                benchmarks being done."""}

        self.logger.debug ( "Latex source code: " + log_message_newline +\
                            page )

        print page 
        filename = self.bench_name + ".tex"
        with open ( filename, "w" ) as text_file:
            text_file.write ( page )
        args = ['pdflatex',
                '-output-directory',
                '%s'%tmp_path,
                '%s'%filename]
        
        self.logger.debug (  "subprocess pdflatex, args =" + log_message_newline +\
                             "[ %s ]" %','.join ("{!s}".format ( l ) for l in args ))
        

        try:
            subprocess.check_call ( args )
        except:
            raise Exception ( "Subprocess call 'pdflatex' failed" )

        pdf_filename = self.bench_name + '.pdf'
        try:
            shutil.move ( os.path.join ( tmp_path, pdf_filename ),
                          os.path.join ( self.report_path, pdf_filename ))
        except:
            raise Exception ( "Moving of report-pdf failed." )

        self.logger.debug ( "Exiting." )


if __name__ == '__main__':

   
    s = "asdksdf_sdjfsd_3 'fsdfsf' sdf_"

    print s

    s2 = teXify ( s )

    print s2
