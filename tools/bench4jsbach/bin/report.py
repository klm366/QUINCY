#! /usr/bin/env python

import os,subprocess
from util import logger as dkl_logger, log_message_newline
import shutil

template = """\documentclass[]{article}
\\usepackage[english]{babel}
\\usepackage{graphicx}
\\usepackage{color}
\\usepackage{titling}
\\usepackage{lipsum}
\\renewcommand{\\familydefault}{\\sfdefault}
\\newcommand{\graphicpath}[1]{%(graphicpath)s/#1}
\\begin{document}
\\title{%(title)s}
\\author{%(author)s, %(institute)s, %(department)s}
\maketitle
\\thispagestyle{empty}

\\begin{abstract}
%(abstract)s \lipsum[1]
\end{abstract}

%(body)s

\section{outlook}
%(epilog)s \lipsum[1]
\end{document}
"""
section_template = """\section{%(section)s}
%(body)s"""

report_dirname = "report"
source_dirname = "sources"

class Report ( object ):
    body = ""
    graphic = None
    analyses = {}

    def __init__ ( self ):
        self.template = template

    def albedo_section ( self ):
        dkl_logger.debug ( "Starting." )
        varname = "albedo"
        print self.analyses

        section = section_template % { 'section' : varname,
                                       'body'    : """
This text is about analyzing '%s'.
%s
""" % ( varname,
        ' \n\n '.join ( """
        \\begin{{center}}
        \includegraphics[width=\\textwidth]{{{!s}}}
        \\end{{center}}""".format ( p ) for p in \
                    self.analyses [ varname ][ 'plots' ]))}

        dkl_logger.debug ( section )
        self.body += section
        dkl_logger.debug ( "Exiting." )
        
    def write_report ( self ):
        dkl_logger.debug ( "Starting." )

        source_path  = os.path.join ( self.report_path, "sources" )
        tmp_path     = os.path.join ( self.report_path, "tmp" )
        if not os.path.exists ( source_path ):
            os.makedirs ( source_path )
        if not os.path.exists ( tmp_path ):
            os.makedirs ( tmp_path )
        os.chdir ( source_path )

        page = template % { 'title'      : 'Bench4JSBACH Report', 
                            'author'     : "%s" % self.config [ 'author' ],
                            'institute'  : "Max-Planck-Institute for Meteorology",
                            'department' : "Land Department",
                            'graphicpath': "/scratch/mpi/mpiles/m300050/bench4jsbach/work/yab/graphic",
                            'abstract'   : "This is an abstract.",
                            'body'       : self.body,
                            'epilog'     : "This is an epilog."}
        print page
        filename = self.bench_name + ".tex"
        with open ( filename, "w" ) as text_file:
            text_file.write ( page )
        args = ['pdflatex',
                '-output-directory',
                '%s'%tmp_path,
                '%s'%filename]
        
        dkl_logger.debug (  "subprocess pdflatex, args =" + log_message_newline +\
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

        dkl_logger.debug ( "Exiting." )
