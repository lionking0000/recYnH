#!/usr/bin/python
import sys
import argparse
import align
import merge

EPILOG = "Commands:\n\
  align      Align the FASTQ sequencing files into bait and prey sequences to generate interaction matrix\n\
  merge      Merge two interaction matries to generate an interaction score matrix\n\n\
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
    align_parser.add_argument('-f1', '--fastq1', required=True, help="set the FASTQ file (read 1 = baits)" ) # Y2H or Y3H
    align_parser.add_argument('-f2', '--fastq2', required=True, help="set the FASTQ file (read 2 = preys)" ) # Y2H or Y3H
    align_parser.add_argument('-o', '--output', required=False, help="set the output directory path (default = same folder as FASTQ file 1)" ) # Y2H or Y3H
    align_parser.add_argument('-n', '--name', default='recYnH', required=False, help="set the output filename (default 'recYnH')" ) # Y2H or Y3H
    
    merge_parser = subparsers.add_parser('merge', help='a help for merge') #, epilog = "Run 'recYnH.py merge --help' for more information on a command.")
    merge_parser.add_argument('-p', '--program', default='Y2H', help="set the experiments type ('Y2H'|'Y3H') (default 'Y2H')" ) # Y2H or Y3H
    merge_parser.add_argument('-m1', '--matrix1', required=True, help="set the interaction matrix of non-selection condition" ) # Y2H or Y3H
    merge_parser.add_argument('-m2', '--matrix2', required=True, help="set the interaction matrix of selection condition" ) # Y2H or Y3H
    merge_parser.add_argument('-o', '--output', required=False, help="set the output folder name (default = same folder as interaction matrix file 1)" ) # Y2H or Y3H
    merge_parser.add_argument('-n', '--name', default='recYnH', required=False, help="set the output filename (default 'recYnH')" ) # Y2H or Y3H   
    
    args = parser.parse_args()

    print args.cmd

    if ( args.cmd == "align" ):
        align.run( args )
    elif ( args.cmd == "merge" ):
        merge.run( args )
    else:
        assert( False )  

    '''
    if "align" not in sys.argv and "merge" not in sys.argv:
        args = parser.parse_args()
    
    if ( "align" in sys.argv ): # or args.cmd == "align"):
        align_parser = argparse.ArgumentParser(prog='recYnH.py align', description='A recYnH interaction matrix generator')
        
        align_parser.add_argument('-p', '--program', default='Y2H', help="set the experiments type ('Y2H'|'Y3H') (default 'Y2H')" ) # Y2H or Y3H
        align_parser.add_argument('-i1', '--fasta1', required=True, help="set the sequence of baits and preys; if i2 is set then it is for baits" ) # Y2H or Y3H
        align_parser.add_argument('-i2', '--fasta2', required=False, help="set the sequence of baits (OPTIONAL)" ) # Y2H or Y3H
        align_parser.add_argument('-f1', '--fastq1', required=True, help="set the FASTQ file (read 1 = baits)" ) # Y2H or Y3H
        align_parser.add_argument('-f2', '--fastq2', required=True, help="set the FASTQ file (read 2 = preys)" ) # Y2H or Y3H
        align_parser.add_argument('-o', '--output', required=False, help="set the output folder name (default 'recYnH_output')" ) # Y2H or Y3H
        
        align_args = align_parser.parse_args()

    if ( "merge" in sys.argv ): # or args.cmd == "merge"):
        merge_parser = argparse.ArgumentParser(prog='recYnH.py merge', description='A recYnH interaction score generator')
        
        merge_parser.add_argument('-p', '--program', default='Y2H', help="set the experiments type ('Y2H'|'Y3H') (default 'Y2H')" ) # Y2H or Y3H
        merge_parser.add_argument('-m1', '--matrix1', required=True, help="set the interaction matrix of selection condition" ) # Y2H or Y3H
        merge_parser.add_argument('-m2', '--matrix2', required=False, help="set the interaction matrix of non-selection condition (OPTIONAL)" ) # Y2H or Y3H
        merge_parser.add_argument('-o', '--output', required=False, help="set the output folder name (default 'recYnH_output')" ) # Y2H or Y3H
        
        merge_args = merge_parser.parse_args()
        print merge_args.program
        print merge_args.matrix1
        print merge_args.matrix2
        print merge_args.output
    '''
        

    #parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                help='an integer for the accumulator')
    #parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                const=sum, default=max,
    #                help='sum the integers (default: find the max)')
    #
    #print(args.accumulate(args.integers))

    #print args.program


    
'''
A self-sufficient runtime for containers

Options:
      --config string      Location of client config files (default "/Users/jyang/.docker")
  -D, --debug              Enable debug mode
  -H, --host list          Daemon socket(s) to connect to
  -l, --log-level string   Set the logging level ("debug"|"info"|"warn"|"error"|"fatal") (default "info")
      --tls                Use TLS; implied by --tlsverify
      --tlscacert string   Trust certs signed only by this CA (default "/Users/jyang/.docker/ca.pem")
      --tlscert string     Path to TLS certificate file (default "/Users/jyang/.docker/cert.pem")
      --tlskey string      Path to TLS key file (default "/Users/jyang/.docker/key.pem")
      --tlsverify          Use TLS and verify the remote
  -v, --version            Print version information and quit

Management Commands:
  checkpoint  Manage checkpoints
  config      Manage Docker configs
  container   Manage containers
  image       Manage images
  network     Manage networks
  node        Manage Swarm nodes
  plugin      Manage plugins
  secret      Manage Docker secrets
  service     Manage services
  swarm       Manage Swarm
  system      Manage Docker
  trust       Manage trust on Docker images
  volume      Manage volumes

Commands:
  attach      Attach local standard input, output, and error streams to a running container
  build       Build an image from a Dockerfile
  commit      Create a new image from a container's changes
  cp          Copy files/folders between a container and the local filesystem

'''