#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  dsl4jsb.py
#
#  Fortran preprocessor for JSBACH DSL language
#
#  Call "dsl4jsb.py -h" for usage.
#  To test, call "dsl4jsb.py" without options in the source directory (containing "dsl4jsb.py" and "test.f90")
#    
#  Nov 10, 2016: Reiner Schnur, MPIM
#  

import sys, getopt, re, os, fnmatch

# Search pattern to find lines with dsl4jsb macros and extract macro name
#re_function = re.compile(r'.*(dsl4jsb_(.*?)( |\())+', re.IGNORECASE)
re_function = re.compile(r'dsl4jsb_(.*?)( |\()', re.IGNORECASE)

# Table of macros
#   Each entry is of the form
#       macro name: [re.compile(match pattern, re.I), replacement pattern]
#   (macro name without the "dsl4jsb_" and in lower case)
#
macros = {
    'use_processes' : [re.compile(
        r'dsl4jsb_Use_processes (.*)',                                   re.IGNORECASE),
        r'USE mo_jsb_process_factory, ONLY: Get_process_memory, Get_process_config, \1'],
    'use_memory' : [re.compile(
        r'dsl4jsb_Use_memory\((.*)\)',                                   re.IGNORECASE),
        r'USE mo_\1_memory_class, ONLY: t_\1_memory'],
    'def_memory' : [re.compile(
        r'dsl4jsb_Def_memory\((.*)\)',                                   re.IGNORECASE),
        r'TYPE(t_\1_memory), POINTER :: \1__mem'],
    'get_memory' : [re.compile(
        r'dsl4jsb_Get_memory\((.*)\)',                                   re.IGNORECASE),
        r'CALL Get_process_memory(tile%mem, \1, \1__mem)'],
    'use_config' : [re.compile(
        r'dsl4jsb_Use_config\((.*)\)',                                   re.IGNORECASE),
        r'USE mo_\1_config_class, ONLY: t_\1_config'],
    'def_config' : [re.compile(
        r'dsl4jsb_Def_config\((.*)\)',                                   re.IGNORECASE),
        r'TYPE(t_\1_config), POINTER :: \1__conf'],
    'get_config' : [re.compile(
        r'dsl4jsb_Get_config\((.*)\)',                                   re.IGNORECASE),
        r'CALL Get_process_config(model%processes(\1)%p%config, \1__conf)'],
    'config' : [re.compile(
        r'dsl4jsb_Config\((.*?)\)',                                      re.IGNORECASE),
        r'\1__conf'],
    'real2d_onchunk' : [re.compile(
        r'dsl4jsb_Real2D_onChunk',                                       re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:)'],
    'real3d_onchunk' : [re.compile(
        r'dsl4jsb_Real3D_onChunk',                                       re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:,:)'],
    'get_var2d_onchunk' : [re.compile(
        r'dsl4jsb_Get_var2D_onChunk\((.*?), *(.*?), *(.*?)\)',           re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(ics:ice,iblk)'],
    'get_var3d_onchunk' : [re.compile(
        r'dsl4jsb_Get_var3D_onChunk\((.*?), *(.*?), *(.*?)\)',           re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(ics:ice,:,iblk)'],
    'var2d_onchunk' : [re.compile(
        r'dsl4jsb_var2D_onChunk\((.*?), *(.*?)\)',                       re.IGNORECASE),
        r'\1__mem%\2%ptr(ics:ice,iblk)'],
    'var3d_onchunk' : [re.compile(
        r'dsl4jsb_var3D_onChunk\((.*?), *(.*?)\)',                       re.IGNORECASE),
        r'\1__mem%\2%ptr(ics:ice,:,iblk)'],
    'real2d_ondomain' : [re.compile(
        r'dsl4jsb_Real2D_onDomain',                                      re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:,:)'],
    'real3d_ondomain' : [re.compile(
        r'dsl4jsb_Real3D_onDomain',                                      re.IGNORECASE),
        r'REAL(wp), POINTER, DIMENSION(:,:,:)'],
    'get_var2d_ondomain' : [re.compile(
        r'dsl4jsb_Get_var2D_onDomain\((.*?), *(.*?), *(.*?)\)',          re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(:,:)'],
    'get_var3d_ondomain' : [re.compile(
        r'dsl4jsb_Get_var3D_onDomain\((.*?), *(.*?), *(.*?)\)',          re.IGNORECASE),
        r'\2 => \1__mem%\2%ptr(:,:,:)'],
    'var2d_ondomain' : [re.compile(
        r'dsl4jsb_var2D_onDomain\((.*?), *(.*?)\)',                      re.IGNORECASE),
        r'\1__mem%\2%ptr(:,:)'],
    'var3d_ondomain' : [re.compile(
        r'dsl4jsb_var3D_onDomain\((.*?), *(.*?)\)',                      re.IGNORECASE),
        r'\1__mem%\2%ptr(:,:,:)'],
    'var2d_ptr' : [re.compile(
        r'dsl4jsb_var2D_ptr\((.*?), *(.*?)\)',                           re.IGNORECASE),
        r'\1__mem%\2%ptr'],
    'var3d_ptr' : [re.compile(
        r'dsl4jsb_var3D_ptr\((.*?), *(.*?)\)',                           re.IGNORECASE),
        r'\1__mem%\2%ptr'],
    }

help_text = """
    dsl4jsb.py [-h] [OPTION...]

        -h                usage
        -v, --verbose     print list of conversions
        -i, --ifile=FILE  input file (test.f90 if missing and input directory (-d) not given)
        -o, --ofile=FILE  output file (stdout if missing)
        -d, --idir=DIR    directory with source files
        -t, --odir=DIR    target directory for output files (must exist)
        -k, --keep-dirs   keep directory tree structure of files in input directory (if -d and -t given)
        -p, --pp-add=STR  add string STR before .suffix in processed files (pp = pre-processed)
    """

def get_input_files(root):
    for dir, subdirs, files in os.walk(root):
        for file in fnmatch.filter(files, '*.f90'):
            if dir == '.':
                yield(file)
            else:
                yield os.path.join(dir.replace('./',''), file)

def process_file(inputfile, outputfile, verbose=False, cwd=None):
    with open(inputfile, 'r') as infile:
        with (open(outputfile, 'w') if outputfile != '' else sys.stdout) as outfile:
            outfile.write("!\n! This file was automatically generated from:\n")
            outfile.write("!         " + os.path.realpath(os.path.abspath(inputfile)) + "\n!\n")
            if verbose:
                if cwd:
                    print '   Pre-processing: ' + os.path.relpath(os.path.realpath(os.path.abspath(inputfile)), cwd) + ' => ' + os.path.relpath(outputfile, cwd)
                else:
                    print '   Pre-processing: ' + os.path.realpath(inputfile) + ' => ' + outputfile
            for line in infile:
                if outputfile == '':
                    outfile.write('! ' + line)
                matches = re.finditer(re_function, line)
                if matches:
                    for match in matches:         # match.groupt(1) contains the macro name without "dsl4jsb_"
                        try:
                            key = match.group(1).lower()
                            pat, replace = macros[key]
                            line = pat.sub(replace, line)
                        except:
                            print 'ERROR in', inputfile, ': dsl4jsb macro not found:', match.group(1).lower()
                            sys.exit()
                outfile.write(line)

def main(argv):
    inputfile = ''
    outputfile = ''
    inputdir = ''
    outputdir = ''
    verbose = False
    keep = False
    pp_add = None
    try:
        opts, args = getopt.getopt(argv,"hvki:o:d:t:p:",["ifile=","ofile=","idir=","odir=","pp-add="])
    except getopt.GetoptError:
        print help_text
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print help_text
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-k", "--keep-dirs"):
            keep = True
        elif opt in ("-p", "--pp-add"):
            pp_add = arg
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg
        elif opt in ("-d", "--idir"):
            inputdir = os.path.abspath(arg)     # Make inputdir absolute
        elif opt in ("-t", "--odir"):
            outputdir = os.path.abspath(arg)    # Make outputdir absolute
            if not os.path.exists(outputdir):
                print 'Output directory ' + outputdir + ' does not exist.'
                sys.exit()
    if inputdir != '' and outputdir == '':
        print 'Input directory specified, but no output directory'
        sys.exit()
    if inputfile != '' and inputdir != '':
        print 'Specify either input file or input directory but not both!'
        sys.exit()
    if inputdir != '' and outputdir != '':
        inputfile = ''
        outputfile = ''
    if inputfile != '' and outputfile == '' and outputdir != '':
        if pp_add:
            root,ext = os.path.splitext(os.path.basename(inputfile))
            outputfile = os.path.join(outputdir, root + pp_add + ext)
        else:
            outputfile = os.path.join(outputdir, os.path.basename(inputfile))
    if inputfile == '' and inputdir == '':
        inputfile = 'test.f90'

    cwd = os.getcwd()    # Save current directory
    if inputdir != '':
        os.chdir(inputdir)   # Change to input directory in order to get relative directory names
        for ifile in get_input_files('.'):     # gets input files including (relative to inputdir) directory name
            if keep:
                if pp_add:
                    root,ext = os.path.splitext(ifile)
                    ofile = os.path.join(outputdir, root + pp_add + ext)
                else:
                    ofile = os.path.join(outputdir, ifile)
                dir = os.path.dirname(ofile)
                if not os.path.exists(dir):
                    os.makedirs(dir)
            else:
                if pp_add:
                    root,ext = os.path.splitext(os.path.basename(ifile))
                    ofile = os.path.join(outputdir, root + pp_add + ext)
                else:
                    ofile = os.path.join(outputdir, os.path.basename(ifile))
            process_file(os.path.join(inputdir, ifile), ofile, verbose, cwd)
        os.chdir(cwd)
    else:
        process_file(inputfile, outputfile, verbose, cwd)

    return 0

if __name__ == '__main__':
     sys.exit(main(sys.argv[1:]))
