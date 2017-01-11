#! /usr/bin/env python
"""The module contains report related classes, functionality and data."""

import sys, os, subprocess, shutil, re
from collections import OrderedDict
import configobj
from netCDF4 import Dataset

from util import (
    get_filesystem_tree_dirtree,
    get_filesystem_tree_tikz,
    log_message_newline,
    print_none,
    print_item_or_list_rec,
    subprocessing )


doc_template = """
\\documentclass[]{article}
\\usepackage[english]{babel}
\\usepackage{graphicx}
\\usepackage{color}

\\usepackage{booktabs,tabularx,ragged2e,makecell, rotating}
\\newcolumntype{Y}{>{\\raggedright\\arraybackslash}X}
\\newcolumntype{L}{>{\\arraybackslash}X}
\
\\usepackage{titling}
\\usepackage{lipsum}

\\usepackage{lmodern}

\\renewcommand{\\familydefault}{\\sfdefault}
\\newcommand{\graphicpath}[1]{%(graphicpath)s/#1}
\\newcommand{\\ttbf}[1]{\\textnormal{\\ttfamily \\bfseries #1}}

\\setlength{\parindent}{0pt}

\\usepackage{geometry}
 \\geometry{
 a4paper,
 total={140mm,247mm},
 left=30mm,
 top=40mm,
 bottom=30mm
 }

%% dina4 : 21 cm x 29,7 cm 

%(use_filesystem_display)s

\\usepackage{hyperref}
\\hypersetup{colorlinks=true, linkcolor=blue, urlcolor=blue}

\\usepackage[printwatermark]{xwatermark}
\\newsavebox\\mybox
\\savebox\\mybox{\\tikz[color=gray,font=\\sffamily,opacity=0.3]\\node{\\textbf{TESTRUN}};}
\\newwatermark*[
  allpages,
  angle=55,
  scale=16,
  xpos=-65,
  ypos=40
]{\\usebox\\mybox}


\\usepackage{xspace}
\\def\\b4j{Bench4JSBACH\\xspace}

\\usepackage[toc,section=section,style=super]{glossaries}
\\makeglossaries

\\newglossaryentry{yearsummean}{name={\\ttbf{yearsummean}}, description={\\texttt{yearsummean}($V_{TX}$) is defined as
\\begin{eqnarray*}
  y_{X} & = & n_{y}^{-1}\\sum_{i=0}^{n_{d}-1} x_{iX},\\\\
&& n_{d} \\textnormal{ number of days, }\\\\
&& n_{y} \\textnormal{ number of years in period (consisting of whole years)}
\\end{eqnarray*}}}

\\newglossaryentry{ysm}{name={\\ttbf{ysm}}, description={\\emph{See} 
\\gls{yearsummean}}}

\\newglossaryentry{zonmean}{name={\\ttbf{zonmean}}, 
description={\\texttt{zonmean}($V_{X}$) is defined as
\\begin{eqnarray*}
  y_{L} & = & n_{L}^{-1}\sum_{i=0}^{n_{L}-1} x_{iL},\\\\
&& n_{L} \\textnormal{ number of points on latitude } L
\\end{eqnarray*}}}

\\newglossaryentry{zm}{name={\\ttbf{zm}}, description={\\emph{See} 
\\gls{zonmean}}}

\\newglossaryentry{fieldmean}{name={\\ttbf{fieldmean}}, description={\\texttt{fieldmean}($V_{TX}$) is defined as
\\begin{eqnarray*}
  y_{tX} & = & n_{X}^{-1}\sum_{i=0}^{n_{X}-1} x_{ti},\\\\
&& n_{X} \\textnormal{ number of gridpoints of field } X
\\end{eqnarray*}}}

\\newglossaryentry{fm}{name={\\ttbf{fm}}, description={\\emph{See} 
\\gls{fieldmean}}}


\\begin{document}

\\title{%(title)s}
\\author{%(author)s, %(institute)s, %(department)s}

\\maketitle

\\tableofcontents

\\thispagestyle{empty}

\\section{Benchmark description}
%(abstract)s

%(body)s

\\section{Results and Summary}
%(epilog)s

\\addcontentsline{toc}{section}{Glossary}
\\printglossaries

\\end{document}
"""
section_template       = """\\section{%(section)s}"""
subsection_template    = """\\subsection{%(subsection)s}"""
subsubsection_template = """\\subsubsection{%(subsubsection)s}"""

use_forest_style = """
\usepackage{tikz}
\usepackage{forest}
\usetikzlibrary{arrows.meta}
\\forestset{
dir tree/.style={
for tree={
parent anchor=south west,
child anchor=west,
anchor=mid west,
inner ysep=-3.5pt,
grow'=0,
align=left,
edge path={
\\noexpand\\path [draw, \\forestoption{edge}] (!u.parent anchor) ++(1em,0) |- (.child anchor)\\forestoption{edge label};
},
font=\\sffamily,
if n children=0{}{
delay={
prepend={[,phantom, calign with current]}
}
},
fit=band,
before computing xy={
l=2em
}
},
}
}
"""

use_dirtree = "\usepackage{dirtree}\n"

report_dirname = "report"
source_dirname = "sources"


def includegraphics_template ( path, width = None, hspace = None, scale = None ):
    
    if width:
        param =  "width=%s\\textwidth" % width
    elif scale:
        param = "scale=%s" % scale
    else:
        param = ""
        
    ret = "\\includegraphics[%s]{%s}" % ( param, path )
    
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
    
    ret = re.sub ( "\"(.*?)\"",r"``\1''", ret )
    ret = re.sub ( "'(.*?)'",r"`\1'", ret )
    ret = ret.replace ( '_', '\_' )
    ret = ret.replace ( '#', '\#' )

    return ret


def lib_avail ( libname ):
    cmd_tmpl = "kpsewhich %(lib)s"
    ret = subprocessing ( cmd_tmpl % { "lib" : libname }, exit_on_error = False )
    if not ret:
        print "laTeX library %s found." % \
            ( libname )
    else:
        print "laTeX library %s not found. Returncode: %d" % \
            ( libname, ret )
        
    return True if ret == 0 else False
    


class Report ( object ):
    diagnostic_fields = OrderedDict ()
    TeXlib_forest  = False
    TeXlib_dirtree = False

    def __init__ ( self ):
        super ( Report, self ).__init__ ()
        self.doc_template = doc_template
        self.TeXlib_forest = lib_avail ( "tikz.sty" ) and lib_avail ( "forest.sty" )
        self.TeXlib_dirtree = lib_avail ( "dirtree.sty" )

    def write_description ( self ):
        self.logger.debug ( "Starting." )

        def writeRepositoryNLocalisation_old ():

# tree ~/git/jsbach/tools/bench4jsbach/ > workingdir_bench4jsbach.txt
# paps workingdir_bench4jsbach.txt | ps2pdf - workingdir_bench4jsbach.pdf


            redmine_url = "https://code.zmaw.de/projects/jsbach/repository/jsbach/show/branches/jsbach-dkl/tools/bench4jsbach"
            
            localtree_rootpath = self.graphics_path + os.sep + "workingdir_tree"
            localtree_path_txt = localtree_rootpath + ".txt"
            localtree_path = localtree_rootpath + ".pdf"
            cmd_localtree_1 = [ "tree -d %s > %s" % ( self.bench_path, localtree_path_txt )]
            cmd_localtree_2 = [ "paps %s | ps2pdf - %s" % ( localtree_path_txt, localtree_path )]
            
            retval = os.getcwd ()
            os.chdir ( os.environ [ 'BENCH4JSBACH_BASE' ] )
            check = subprocess.check_output ( cmd_localtree_1, shell=True )
            check = subprocess.check_output ( cmd_localtree_2, shell=True )
            os.chdir ( retval )

            ret = """\\-\\bigskip
\-\\bigskip


This report was generated by \\textbf{Bench4JSBach} from \\url{%s} with commit {\\tt %s}, see also on the ZMAW (Zentrum f{\\"u}r Marine und Atmosph{\\"a}rische Wissenschaften) redmine web site \\url{%s}. After runtime the working directory of this diagnostic run was located as follows:\\\\\\bigskip

  %% crop/trim a pdf when including it using the trim=left botm right top
\\hspace*{0cm}\\includegraphics[clip, trim=1.cm 22.cm 1.cm 1.cm, width=\\textwidth]{%s}

""" % (
    self.config [ "git_repository" ][ 0 ].split ( ':' )[ 1 ] ,
    self.config [ "git_sha" ],
    redmine_url,
    localtree_path )


            self.logger.debug ( "Exiting." )
            return ret

        def writeRepositoryNLocalisation ():
            self.logger.debug ( "Starting." )
            
            redmine_url = "https://code.zmaw.de/projects/jsbach/repository/jsbach/show/branches/jsbach-dkl/tools/bench4jsbach"


            

            ret = """
\\-\\bigskip
\\-\\bigskip

This report was generated by \\textbf{Bench4JSBach} from \\url{%s} with commit {\\tt %s}, see also on the ZMAW (Zentrum f{\\"u}r Marine und Atmosph{\\"a}rische Wissenschaften) redmine web site \\url{%s}. After runtime the working directory of this diagnostic run was located as follows:\\\\\\bigskip
%s
\-\\\\\\bigskip
\-\\\\\\bigskip
""" % (
    self.config [ "git_repository" ][ 0 ].split ( ':' )[ 1 ] ,
    self.config [ "git_sha" ],
    redmine_url,
    teXify ( get_filesystem_tree_tikz (  self.bench_path ) if self.TeXlib_forest else \
             get_filesystem_tree_dirtree ( self.bench_path ) if self.TeXlib_dirtree else "" )
)

            self.logger.debug ( "Exiting." )
            return ret
        

        def writeDiagnosticsTabularx ():
            ret = ""
            for field in self.diagnostic_fields.keys ():
                fieldrow = 0
                ret += "\\midrule\n"
                for bench, _dict in { k : v for k, v in self.diagnostic_fields [ field ][ 'report' ].iteritems () \
                                      if isinstance ( v, configobj.Section )}.iteritems () :
                    ret += "&" if fieldrow > 0 else "%s & " % field 
                    ret += ' %s & %s & %s & %s \\\\' % \
                           ( teXify ( _dict [ 'title' ]),
                             teXify ( _dict [ 'name' ]),
                             teXify ( _dict [ 'type' ]),
                              "p.\\pageref{%s}\n" % _dict [ 'name' ])

                    fieldrow += 1
            
            return ret

        description = """
In this benchmark report, the candidate experiments
\\begin{itemize}
%s
\\end{itemize}
""" % ''.join ( '\\item {!s} ($C_{{{!s}}}$), selected period {!s}-{!s}'.\
                format ( teXify ( can.B4J_name ),
                         str ( can.B4J_canno ),
                         can.B4J_startyear,
                         can.B4J_endyear ) \
                if can.B4J_iscan else '' for can in self.experiments ) +\
"""
are compared with the reference experiments
\\begin{itemize}
%s
\\end{itemize} 
\\ \\\\ \\\\
""" % ''.join ( '\\item {!s} ($R_{{{!s}}}$), selected period {!s}-{!s}'.\
                format ( teXify ( ref.B4J_name ),
                         str ( ref.B4J_refno ),
                         ref.B4J_startyear,
                         ref.B4J_endyear ) \
              if ref.B4J_isref else '' for ref in self.experiments ) +\
"""
The fields tested and the diagnostics and plots made are listed in the 
table below. You may find detailed information on each method following the appropriate link.
\\ \\\\
\\begin{center}
\\begin{tabular}{rrrrr}
\\toprule
\\textbf{field} & \\textbf{diagnostic} & \\textbf{name} & \\textbf{type} & \\textbf{page} \\\\
\\addlinespace
%s
\\bottomrule
\\end{tabular}
\\end{center}

""" %  writeDiagnosticsTabularx ()

        description += \
        writeRepositoryNLocalisation () +\
            subsection_template % { 'subsection': 'Candidate experiments to benchmark' }
    
        for can in self.experiments:
            if not can.B4J_iscan: continue
            description +=  can.teX_repr ()

        description += subsection_template % { 'subsection': 'Reference experiments' }
        for ref in self.experiments:
            if not ref.B4J_isref: continue
            if ref.B4J_iscan:
                description +=  """
%s %s
\\\\\\medskip
""" % ( ref.teX_short_repr ( shortname = 'ref' ),
         ' see p.\\pageref{%s} ({\\tt $C_{%d}$})' % (  ref.B4J_name, ref.B4J_canno ))
            else: 
                description +=  ref.teX_repr ()
                
        return description


    def report_diagnostic ( self ):
        self.logger.debug ( "Starting." )

        dss = ""

        for f_n, f_dict in self.diagnostic_fields.iteritems ():
            fdr =  f_dict [ 'report' ]
            dss += section_template % { 'section' : f_n.capitalize () + \
                                        ( ": " + fdr [ 'title' ] if 'title' in fdr else "" )}
            dss += fdr [ 'abstract' ] if 'abstract' in fdr else ""
            
            sss = [ k for k,v in fdr.iteritems () if isinstance ( v, configobj.Section )]
            for subsection in sss:
                ss_dict = fdr [ subsection ]
                dss += "" if not "name" in ss_dict else \
                       "\\label{%s}\\hypertarget{%s}{}" % ( ss_dict [ 'name' ], ss_dict [ 'name' ] )
                dss += subsection_template % { 'subsection' : ss_dict [ 'title' ].capitalize ()}
                dss += teXify ( ss_dict [ 'abstract' ])
                dss += """\\phantom \\\\\\medskip
"""
                dss += ss_dict [ 'technical_teXed' ]
                
                for can in self.experiments:
                    if not can.B4J_iscan: continue
                    can_name = teXify ( can.B4J_name )
                    title = "Candidate experiment {\\tt %s} ($\mathbf{C_{%d}}$)" %\
                            ( can_name, can.B4J_canno ) 
                    dss += subsubsection_template % { 'subsubsection': title }

                    if ss_dict [ 'ismatrix' ] == True:
                        
                        if "canplot" in ss_dict:
                            plot = can.B4J_index [ ss_dict [ "canplot" ]]
                            dss += includegraphics_template (  plot,
                                                               width      = '1.5',
                                                               hspace     = '-3.3' )
                        
                        for ref in self.experiments:
                            if not ref.B4J_isref:            continue
                            if ref.B4J_name == can.B4J_name: continue
                            plot = can.B4J_index [ '%s_plot_ref%d' % ( subsection, ref.B4J_refno )]
                            dss += includegraphics_template (  plot,
                                                               width      = '1.5',
                                                               hspace     = '-3.3' )
                    elif ss_dict [ 'isvector' ] == True:
                        
                        if "canplot" in ss_dict:
                            plot = can.B4J_index [ ss_dict [ "canplot" ]]
                            dss += "\\raggedright{%s}" % \
                                   includegraphics_template (  plot,
                                                               scale = .5,
                                                               hspace = '-.5' )
                            plot = can.B4J_index [ '%s_plot' %  subsection ]
                            dss += "\\raggedright{%s}" % \
                                   includegraphics_template (  plot,
                                                               scale = .5,
                                                               hspace = '-.5' )
                        else:
                            plot = can.B4J_index [ '%s_plot' %  subsection ]
                            dss += includegraphics_template ( plot,
                                                              scale = .56,
                                                              hspace = '-1' )              
                
        self.logger.debug ( "Exiting." )
        return dss
    

    def write_summary ( self ):
        self.logger.debug ( "Starting." )
        
        def writeBenchmarkTabularX ():
            n_col_static = 4
            n_col = n_col_static + self.n_cans
            ret = ""

            tablerow = 0
            for field in self.diagnostic_fields.keys ():
                fieldrow = 0
                for bench, _dict in { k : v for k, v in self.diagnostic_fields [ field ][ 'report' ].iteritems () \
                                      if isinstance ( v, configobj.Section ) and \
                                         'isbench' in v and v [ 'isbench' ] == True and \
                                         'benchvar' in v }.iteritems () :
                    benchrow = 0
                    for ref in self.experiments:
                        if not ref.B4J_isref: continue
                        
                        ret += "\\midrule" if tablerow > 0 and \
                               fieldrow == 0 else ""
                        ret += "\\cmidrule(r){2-%d}" % n_col if fieldrow > 0 and \
                               benchrow == 0 else ""
                        ret += ( "&" if fieldrow > 0 else " %s & " % field )
                        ret += ( "&&" if benchrow > 0 else "\\hyperlink{%s}{%s} & %s &" % (
                            _dict [ 'name' ], teXify ( _dict [ 'name' ]), _dict [ 'benchunit' ]))
                        ret += ( '\\hyperlink{%s}{$\mathbf{R_{%d}}$}' % ( ref.B4J_name, ref.B4J_refno ))
                        for can in self.experiments:
                            if not can.B4J_iscan: continue
                            
                            if ref.B4J_name == can.B4J_name:
                                ret += "& -"
                            else:                            
                                classification_str = "%s_ref%d" % ( bench,  ref.B4J_refno )
                                path = can.B4J_index [ classification_str ]
                                if not os.path.exists ( path ):
                                    print "File not found."
                                    sys.exit ( 1 )
                                nc_fid   = Dataset (  can.B4J_index [ classification_str ], 'r')
                                benchval = nc_fid.variables [ _dict [ 'benchvar' ]][ : ][ 0 ][ 0 ][ 0 ]
                                unit     = nc_fid.variables [ _dict [ 'benchvar' ]].units
                                nc_fid.close ()
                                ret += ( '& %.4e' % benchval )
                        ret += "\\\\\n"
                        benchrow += 1
                        fieldrow += 1
                        tablerow += 1
                        
            return ret
      
        self.logger.debug ( "Exiting." )
        return """This chapter provides an overview of the benchmarks being done. For evaluation of the values follow the link of the diagnostics' qualification string. The associated chapter should supply evaluation criteria. \\\\\\medskip
        
        \\begin{tabularx}{\\linewidth}{rrrr%s}
        """ % ( self.n_cans * "r" ) +\
        """
        \\toprule
        \\textbf{field} & \\textbf{analyis} & \\textbf{unit} & \\textbf{reference} &
        \\multicolumn{%d}{p{%.1fcm}}{\\RaggedLeft \\textbf{candidate}} \\\\
        """ % ( self.n_cans,  self.n_cans * 1.5 ) +\
        """
        \\addlinespace
        & & & & %s \\\\
        \\addlinespace
        \\midrule
        """ % ( ' & '.join ( ' \\textbf{%s}' % teXify ( can.B4J_name ) if len ( can.B4J_name ) <= 10 else \
        '\\hyperlink{%s}{$\mathbf{C_{%d}}$}' % (can.B4J_name, can.B4J_canno ) for can in self.experiments if can.B4J_iscan )) +\
        """
        %s
        \\bottomrule
        \\end{tabularx}
        
        """ % writeBenchmarkTabularX ()
    
            
    def write_report ( self ):
        self.logger.debug ( "Starting." )
        source_path  = os.path.join ( self.report_path, "sources" )
        tmp_path     = os.path.join ( self.report_path, "temp" )
        if not os.path.exists ( source_path ):
            os.makedirs ( source_path )
        if not os.path.exists ( tmp_path ):
            os.makedirs ( tmp_path )
        os.chdir ( source_path )

        title_second = "%s" % ''.join ( '\\\\{!s} ({!s})'. format ( teXify ( can.B4J_name ),  teXify ( can.B4J_etype )) \
                if can.B4J_iscan else '' \
                for can in self.experiments )

        use_filesystem_display = use_forest_style if self.TeXlib_forest else \
                                 use_dirtree if self.TeXlib_dirtree else ""

        page = doc_template % {
            'use_filesystem_display' : use_filesystem_display,
            'title'      : """{\\Large Benchmark Report for Experiments} 
{\\LARGE %s}\\\\
{\\normalsize \\textbf{Bench4JSBach}}""" % title_second, 
            'author'     : "%s" % self.config [ 'author' ],
            'institute'  : "%s" % self.config [ 'institute' ],
            'department' : "%s" % self.config [ 'department' ],
            'graphicpath': self.graphics_path,
            'abstract'   : self.write_description (),
            'body'       : self.report_diagnostic (),
            'epilog'     : self.write_summary ()}

        filename = self.bench_name + ".tex"
        with open ( filename, "w" ) as text_file:
            text_file.write ( page )
        # args = ['pdflatex',
        #         '-output-directory',
        #         '%s'%tmp_path,
        #         '%s'%filename]
        args = ['pdflatex',
                '%s'%filename]
        
        self.logger.debug (  "subprocess pdflatex, args =" + log_message_newline +\
                             "[ %s ]" %','.join ("{!s}".format ( l ) for l in args ))
        

        for i in range ( 2 ):
            try:
                subprocess.check_call ( args )
            except:
                raise Exception ( "Subprocess call 'pdflatex' failed" )

           # makeindex -l -s xtimeperiod.ist -o xtimeperiod.gls xtimeperiod.glo
        args_makeindex = [ 'makeindex',
                           '-l',
                           '-s',
                           '%s.ist' % self.bench_name,
                           '-o',
                           '%s.gls' % self.bench_name,
                           '%s.glo' % self.bench_name ]

        print os.getcwd()
        print args_makeindex

        try:
            subprocess.check_call ( args_makeindex )
        except:
            raise Exception ( "Subprocess call 'makeindex' failed" )

        try:
            subprocess.check_call ( args )
        except:
            raise Exception ( "Subprocess call 'pdflatex after' failed" )

        pdf_filename = self.bench_name + '.pdf'
        try:
            shutil.move ( pdf_filename,
                          os.path.join ( self.report_path, pdf_filename ))
        except:
            raise Exception ( "Moving of report-pdf failed." )

        self.logger.debug ( "Exiting." )
        
    def __repr__ ( self ):
        super ( Report, self ).__repr__ ()
        rep =  """\n\nSuperclass object {} ()\n\n""".\
              format ( Report.__name__ )
                
        return rep

if __name__ == '__main__':

   
    s = "asdksdf_sdjfsd_3 'fsdfsf' sdf_"

    print s

    s2 = teXify ( s )

    print s2
