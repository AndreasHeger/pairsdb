####
####
##
## Project PairsDB
##
## Copyright (C) 2002 Andreas Heger All rights reserved
##
## Author: Andreas Heger <heger@ebi.ac.uk>
##
## $Id: Nrdb.py,v 1.1.1.1 2002/07/02 10:46:57 heger Exp $
##
##
####
####

USAGE="""python nrdb2rsdb.py [OPTIONS] nrdb.fasta

create an rsdb table from an nrdb fasta file using
cd-hit.

CMD_NRDB90=$(DIR_PROGS)cd-hit -M 2000 -n 5 -c 0.90 -d 40
CMD_NRDB70=$(DIR_PROGS)cd-hit -M 2000 -n 5 -c 0.70 -d 40
"""

# cd-hit commandline options
LEVELS= {
    90 : "-n 5 -c 0.90 -d 40",
    70 : "-n 5 -c 0.70 -d 40" }

import sys, re, string, optparse, time, os, tempfile, shutil, subprocess

import Experiment

if __name__ == '__main__':

    parser = optparse.OptionParser( version = "%prog version: $Id$", usage=USAGE )

    parser.add_option( "-t", "--tempdir", dest="tempdir", type="string" ,
                       help="temporary directory use. Default is system default [/tmp].")

    parser.add_option( "--filename-rsdb-table", dest="output_filename_table", type="string" ,
                       help="OUTPUT filename with rsdb sequences in table format.")

    parser.add_option( "--filename-rsdb-fasta", dest="output_filename_fasta", type="string" ,
                       help="OUTPUT filename with rsdb sequences in fasta format.")

    parser.add_option( "--filename-pairsdb-table", dest="output_pairsdb_table", type="string" ,
                       help="OUTPUT filename with pairsdb entries - groupies of each cluster.")

    parser.add_option( "--level", dest="level", type="choice",
                       choices=('90', '70'),
                       help="level of target rsdb.")

    parser.add_option( "--memory", dest="memory", type="int",
                       help="memory reserved for cd-hit use." )

    parser.set_defaults( executable="cd-hit",
                         memory = 2000,
                         level = '90',
                         output_filename_fasta = "rsdb.fasta",
                         output_filename_pairsdb = "pairsdb.table",
                         output_filename_table = "rsdb.table",                         
                         tempdir = None)

    (options, args) = Experiment.Start( parser, add_pipe_options = True )

    options.level = int(options.level)

    if len(args) != 1:
        raise "please supply style fasta file."

    if options.tempdir:
        tempdir = tempfile.mkdtemp( dir = options.tempdir )
    else:
        tempdir = tempfile.mkdtemp()        

    
    if options.loglevel >= 2:
        options.stdlog.write("# working diretory: %s\n" % tempdir )

    if options.level not in LEVELS:
        raise "unsupported level %i: implemented are: %s" % (options.level, "".join( map(str, LEVELS.keys())))

    input_filename = os.path.abspath( args[0] )
    output_filename = "result"
    
    statement = "%s -M %i -i %s -o %s %s" % (options.executable,
                                             options.memory,
                                             input_filename,
                                             output_filename,
                                             LEVELS[options.level])
    
    if options.loglevel >= 1:
        options.stdlog.write("# executing statement: %s\n" % statement )
        options.stdlog.flush()

    time_t0 = time.time()

    s = subprocess.Popen( statement,
                          shell = True,
                          stdout = subprocess.PIPE,
                          stderr = subprocess.PIPE,
                          cwd = tempdir,
                          close_fds = True)                              

    (out, err) = s.communicate()

    if s.returncode != 0:
        raise "Execution error\n%s\n%s\nTemporary directory in %s" % (err, out, tempdir)
    
    if options.loglevel >= 1:
        options.stdlog.write("# finished execution in %i seconds.\n" % (time.time() - time_t0 ))
        options.stdlog.write("# cd-hit output:\n%s\n" % out )
        options.stdlog.flush()
    
    infile = open( tempdir + "/" + output_filename + ".clstr" )

    outfile_pairsdb = open( options.output_filename_pairsdb, "w" )
    outfile_table = open( options.output_filename_table, "w" )
    
    rx = re.compile( ">(\d+)" )

    def __print( rep_nid, mem_nids ):
        if rep_nid:
            outfile_table.write( rep_nid + "\n" )                     
            for mem_nid in mem_nids:
                if mem_nid != rep_nid:
                    outfile_pairsdb.write( rep_nid + "\t" + mem_nid + "\n")               

    rep_nid, mem_nids = 0, []
    
    for line in infile:
        if line[0] == ">":
            __print(rep_nid, mem_nids )
            rep_nid, mem_nids = 0, []
            continue

        nid = rx.search( line ).groups()[0]
        
        if line[-2] == "*":
            rep_nid = nid
        else:
            mem_nids.append(nid)

    __print(rep_nid, mem_nids )            

    infile.close()
    outfile_pairsdb.close()
    outfile_table.close()
    
    os.system( "mv %s %s" % (tempdir + "/" + output_filename,
                             options.output_filename_fasta ) )
        
    shutil.rmtree( tempdir )
    
    Experiment.Stop()
