#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  dsl4jsb.py
#
#  Copyright 2016 Reiner Schnur <reiner@reinerpc-vb>
#

import sys, getopt, re

# Search pattern to find lines with dsl4jsb macros and extract macro name
re_function = re.compile(r'.*dsl4jsb_(.*?)( |\()', re.IGNORECASE)

# Table of macros
#   Each entry is of the form
#       macro name: [re.compile(match pattern, re.I), replacement pattern]
#   (macro name without the "dsl4jsb_")
#
macros = {
    'use_processes' :
        [re.compile(r'dsl4jsb_Use_processes (.*)', re.IGNORECASE),
         r'USE mo_jsb_process_factory, ONLY: Get_processes_memory, Get_process_config, \1'],
    'use_memory' :
        [re.compile(r'dsl4jsb_Use_memory\((.*)\)', re.IGNORECASE),
         r'USE mo_\1_memory_class, ONLY: t_\1_memory'],
    'def_memory' :
        [re.compile(r'dsl4jsb_Def_memory\((.*)\)', re.IGNORECASE),
         r'TYPE(t_\1_memory) :: \1__mem'],
    'get_memory' :
        [re.compile(r'dsl4jsb_Get_memory\((.*)\)', re.IGNORECASE),
         r'CALL Get_process_memory(tile%mem, \1, \1__mem)'],
    'use_config' :
        [re.compile(r'dsl4jsb_Use_config\((.*)\)', re.IGNORECASE),
         r'USE mo_\1_config_class, ONLY: t_\1_config'],
    'def_config' :
        [re.compile(r'dsl4jsb_Def_config\((.*)\)', re.IGNORECASE),
         r'TYPE(t_\1_config) :: \1__conf'],
    'get_config' :
        [re.compile(r'dsl4jsb_Get_config\((.*)\)', re.IGNORECASE),
         r'CALL Get_process_config(model%processes(\1)%p%config, \1__conf)'],
    'config' :
        [re.compile(r'dsl4jsb_Config\((.*?)\)', re.IGNORECASE),
         r'\1__conf'],
    'real2d_onchunk' :
        [re.compile(r'dsl4jsb_Real2D_onChunk', re.IGNORECASE),
         r'REAL(wp), POINTER, DIMENSION(:)'],
    'real3d_onchunk' :
        [re.compile(r'dsl4jsb_Real3D_onChunk', re.IGNORECASE),
         r'REAL(wp), POINTER, DIMENSION(:,:)'],
    'get_var2d_onchunk' :
        [re.compile(r'dsl4jsb_Get_var2D_onChunk\((.*?), *(.*?), *(.*?)\)', re.IGNORECASE),
         r'\2 => \1__mem%\2%ptr(ics:ice,iblk); CALL \1__mem%\2%Set_inout("\3")'],
    'get_var3d_onchunk' :
        [re.compile(r'dsl4jsb_Get_var3D_onChunk\((.*?), *(.*?), *(.*?)\)', re.IGNORECASE),
         r'\2 => \1__mem%\2%ptr(ics:ice,:,iblk); CALL \1__mem%\2%Set_inout("\3")'],
    'var2d_onchunk' :
        [re.compile(r'dsl4jsb_var2D_onChunk\((.*?), *(.*?)\)', re.IGNORECASE),
         r'\1__mem%\2%ptr(ics:ice,iblk)'],
    'var3d_onchunk' :
        [re.compile(r'dsl4jsb_var3D_onChunk\((.*?), *(.*?)\)', re.IGNORECASE),
         r'\1__mem%\2%ptr(ics:ice,:,iblk)'],
    'real2d_ondomain' :
        [re.compile(r'dsl4jsb_Real2D_onDomain', re.IGNORECASE),
         r'REAL(wp), POINTER, DIMENSION(:,:)'],
    'real3d_ondomain' :
        [re.compile(r'dsl4jsb_Real3D_onDomain', re.IGNORECASE),
         r'REAL(wp), POINTER, DIMENSION(:,:,:)'],
    'get_var2d_ondomain' :
        [re.compile(r'dsl4jsb_Get_var2D_onDomain\((.*?), *(.*?), *(.*?)\)', re.IGNORECASE),
         r'\2 => \1__mem%\2%ptr(:,:); CALL \1__mem%\2%Set_inout("\3")'],
    'get_var3d_ondomain' :
        [re.compile(r'dsl4jsb_Get_var3D_onDomain\((.*?), *(.*?), *(.*?)\)', re.IGNORECASE),
         r'\2 => \1__mem%\2%ptr(:,:,:); CALL \1__mem%\2%Set_inout("\3")'],
    'var2d_ondomain' :
        [re.compile(r'dsl4jsb_var2D_onDomain\((.*?), *(.*?)\)', re.IGNORECASE),
         r'\1__mem%\2%ptr(:,:)'],
    'var3d_ondomain' :
        [re.compile(r'dsl4jsb_var3D_onDomain\((.*?), *(.*?)\)', re.IGNORECASE),
         r'\1__mem%\2%ptr(:,:,:)'],
    'var2d_ptr' :
        [re.compile(r'dsl4jsb_var2D_ptr\((.*?), *(.*?)\)', re.IGNORECASE),
         r'\1__mem%\2%ptr'],
    'var3d_ptr' :
        [re.compile(r'dsl4jsb_var3D_ptr\((.*?), *(.*?)\)', re.IGNORECASE),
         r'\1__mem%\2%ptr'],
    }

def main(argv):
    inputfile = 'test.f90'
    outputfile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
        print 'dsl4py.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'dsl4jsb.py -i <inputfile> -o <outputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    with open(inputfile, 'r') as infile, (open(outputfile, 'w') if outputfile != '' else sys.stdout) as outfile:
        outfile.write("!\n! This file was automatically generated from:\n")
        if outputfile != '':
            outfile.write("!    dsl4jsb.py -i " + inputfile + " -o " + outputfile + "\n!\n")
        else:
            outfile.write("!    dsl4jsb.py -i " + inputfile + "\n!\n")
        for line in infile:
            match = re.match(re_function, line)  # match.groupt(1) contains the macro name without "dsl4jsb_"
            if match:
                try:
                    pat, replace = macros[match.group(1).lower()]
                    line = pat.sub(replace, line)
                except:
                    print 'ERROR in', inputfile, ': dsl4jsb macro not found:', match.group(1)
                    sys.exit()
            outfile.write(line)

    return 0

if __name__ == '__main__':
     sys.exit(main(sys.argv[1:]))
