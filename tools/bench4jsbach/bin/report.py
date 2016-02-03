#! /usr/bin/env python

import os, subprocess, shutil, re
import json

from util import log_message_newline

doc_template = """
\\documentclass[]{article}
\\usepackage[english]{babel}
\\usepackage{graphicx}
\\usepackage{color}

\\usepackage{booktabs,tabularx,ragged2e}
\\newcolumntype{Y}{>{\\raggedright\\arraybackslash}X}
\\newcolumntype{L}{>{\\arraybackslash}X}
\
\\usepackage{titling}
\\usepackage{lipsum}
\\renewcommand{\\familydefault}{\\sfdefault}
\\newcommand{\graphicpath}[1]{%(graphicpath)s/#1}

\\begin{document}

\\title{%(title)s}
\\author{%(author)s, %(institute)s, %(department)s}
\maketitle
\\thispagestyle{empty}

\section{Benchmark description}
%(abstract)s

%(body)s

\section{Outlook}
%(epilog)s \lipsum[1]

\end{document}
"""
section_template = """\section{%(section)s}
%(body)s"""

subsection_template = """\subsection{%(subsection)s}"""

includegraphics_template = """
\\begin{center}
\includegraphics[width=\\textwidth]{%(graphicpath)s}
\\end{center}
"""

report_dirname = "report"
source_dirname = "sources"


def lateXify ( string ):
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

        description = ""

        description += subsection_template % { 'subsection': 'Experiments to benchmark' }
        
        from bench import REF_EXP_IDX
        for e in self.experiments:
            if self.experiments.index ( e ) == REF_EXP_IDX:
                continue
            description += lateXify ( e.teX_repr ())

        description += subsection_template % { 'subsection': 'Reference experiments' }
        description += lateXify ( self.experiments [ REF_EXP_IDX ] .teX_repr ())

        print description

        return description 
    

    def albedo_section ( self ):
        self.logger.debug ( "Starting." )
        varname = "albedo"
        albedo_body = "" 
        self.logger.debug ( "self.analyses buffer: " + log_message_newline +\
                            json.dumps ( self.analyses ))

        # from bench import REF_EXP_IDX
        
        # for e in self.experiments:
        #     e_idx = self.experiments.index ( e )
        #     ref = True if e_idx == REF_EXP_IDX else False
        #     name =  e.B4J_name.replace ( '_','\_' )
        #     title = "Reference experiment {\\tt %s}" % name if ref else \
        #             "Experiment {\\tt %s}" % name
        #     albedo_body += subsection_template % { 'subsection': title }
        #     plot =  self.analyses [ varname ][ 'plots' ][ e_idx ]
        #     if plot:
        #         albedo_body += includegraphics_template % { 'graphicpath': plot }
        #     diff_plot =  self.analyses [ varname ][ 'diff_plots' ][ e_idx ]
        #     if diff_plot:
        #         albedo_body += """
        #     %s difference plot to reference experiment {\\tt %s} :\n
        #     %s """ % ( varname.title (), e.B4J_name.replace ( '_','\_' ),
        #     includegraphics_template %  { 'graphicpath': diff_plot })
        #     #print "e_idx", e_idx,"albedo_body", albedo_body

        # section = section_template % { 'section' : varname.title (),
        #                                'body' : albedo_body }
        # #print "section", section
        # self.body += section
    

        
    def write_report ( self ):
        self.logger.debug ( "Starting." )
        source_path  = os.path.join ( self.report_path, "sources" )
        tmp_path     = os.path.join ( self.report_path, "tmp" )
        if not os.path.exists ( source_path ):
            os.makedirs ( source_path )
        if not os.path.exists ( tmp_path ):
            os.makedirs ( tmp_path )
        os.chdir ( source_path )

        page = doc_template % { 'title'      : 'Bench4JSBACH Report', 
                                'author'     : "%s" % self.config [ 'author' ],
                                'institute'  : "Max-Planck-Institute for Meteorology",
                                'department' : "Land Department",
                                'graphicpath': self.graphic_path,
                                'abstract'   : self.write_description (),
                                'body'       : self.body,
                                'epilog'     : "This is an epilog."}
        self.logger.debug ( "Latex source code: " + log_message_newline +\
                            page )
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

    s2 = lateXify ( s )

    print s2
