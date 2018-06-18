#!/usr/bin/python
import sys
import os

sys.path.append( "." )

import argparse
from src import align
from src import score
from src import merge

EPILOG = "Commands:\n\
  align      Align the FASTQ sequencing files into bait and prey sequences to generate interaction matrix\n\
  score      Usging two interaction matries to generate an interaction score matrix\n\
  merge      Merge several interaction score matries and apply quartile correction to generate a final average interaction score matrix\n\n\
Run 'recYnH.py COMMAND --help' for more information on a command."

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='recYnH program',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog = EPILOG )
    
    #parser.add_argument('cmd', metavar='COMMAND', help="set COMMAND ('align'|'merge')", default="" )
    subparsers = parser.add_subparsers(dest='cmd', metavar="COMMAND", help='sub-command help')
    align_parser = subparsers.add_parser('align', help='a help for align')
    align_parser.add_argument('-p', '--program', default='Y2H', help="set the experiments type ('Y2H'|'Y3H') (default 'Y2H')" ) # Y2H or Y3H
    align_parser.add_argument('-r', '--relaxed', action='store_true', help="set relaxed alignment mode" ) # Y2H or Y3H
    align_parser.add_argument('-i1', '--fasta1', required=True, help="set the sequence of baits and preys; if i2 is set then it is for baits" ) # Y2H or Y3H
    align_parser.add_argument('-i2', '--fasta2', required=False, help="set the sequence of baits (OPTIONAL)" ) # Y2H or Y3H
    align_parser.add_argument('-l1', '--lastnt1', default=150, type=int, required=False, help="set the sequence of length of baits (default 150nt)" ) # Y2H or Y3H
    align_parser.add_argument('-l2', '--lastnt2', default=150, type=int, required=False, help="set the sequence of length of preys (default 150nt)" ) # Y2H or Y3H
    align_parser.add_argument('-f1', '--fastq1', required=True, help="set the FASTQ file (read 1 = baits)" ) # Y2H or Y3H
    align_parser.add_argument('-f2', '--fastq2', required=True, help="set the FASTQ file (read 2 = preys)" ) # Y2H or Y3H
    align_parser.add_argument('-o', '--output', required=False, help="set the output directory path (default = same folder as FASTQ file 1)" ) # Y2H or Y3H
    align_parser.add_argument('-n', '--name', default='recYnH.raw', required=False, help="set the output filename (default 'recYnH.raw')" ) # Y2H or Y3H
    
    score_parser = subparsers.add_parser('score', help='a help for score') #, epilog = "Run 'recYnH.py score --help' for more information on a command.")
    score_parser.add_argument('-p', '--program', default='Y2H', help="set the experiments type ('Y2H'|'Y3H') (default 'Y2H')" ) # Y2H or Y3H
    score_parser.add_argument('-m1', '--matrix1', required=True, help="set the interaction matrix of non-selection condition" ) # Y2H or Y3H
    score_parser.add_argument('-m2', '--matrix2', required=True, help="set the interaction matrix of selection condition" ) # Y2H or Y3H
    score_parser.add_argument('-o', '--output', required=False, help="set the output folder name (default = same folder as interaction matrix file)" ) # Y2H or Y3H
    score_parser.add_argument('-n', '--name', default='recYnH.IS', required=False, help="set the output filename (default 'recYnH.IS')" ) # Y2H or Y3H   
    
    merge_parser = subparsers.add_parser('merge', help='a help for merge') #, epilog = "Run 'recYnH.py merge --help' for more information on a command.")
    merge_parser.add_argument('-i', '--input', required = True, nargs = '*', help="list of interaction score matrices" ) # Y2H or Y3H
    merge_parser.add_argument('-q', '--quartile', default=0.75, required=False, help="set the percentile auto-activation signal correction value (default 0.75)" ) # Y2H or Y3H
    merge_parser.add_argument('-o', '--output', required=False, help="set the output folder name (default = same folder as interaction matrix file 1)" ) # Y2H or Y3H
    merge_parser.add_argument('-n', '--name', default='recYnH.avgIS', required=False, help="set the output filename (default 'recYnH.avgIS')" ) # Y2H or Y3H   
    args = parser.parse_args()

    if args.cmd not in [ "align", "score", "merge" ]:
        exit(0 )

    print "[ Starting rec-YnH Pipeline ]"

    if ( args.cmd == "align" ):
        align.run( args )
    elif ( args.cmd == "score" ):
        score.run( args )
    elif ( args.cmd == "merge" ):
        merge.run( args )

    print "[ Finishing rec-YnH Pipeline ]"