#!/usr/bin/python
import os
import threading

VERBOSE = True # False

NUM_THREADS = 1 #

def run_cmd( cmd ):
    if VERBOSE: print cmd
    os.system( cmd )

def align_subprocess( original_fasta, fasta, fastq  ):
    # convert fastq to fasta ( Temporarily now using fastq file generated in Friedrich folder; since it is the same fastq. But We need to change it to Blastn for general cases )
    #cmd = "python main.py fastq_to_fasta %s > %s" % ( fastq, fasta )
    cmd = "main.py fastq_to_fasta %s > %s" % ( fastq, fasta )
    run_cmd( cmd )
    
    # blastn-short search                                                                         (20) 
    cmd = "blastn -num_threads %d -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 5 -evalue 1e-8 > %s.blastn" % ( NUM_THREADS, original_fasta, fasta, fasta )
    run_cmd( cmd )
   

def run( args ):
    #print args.program
    #print args.fasta1
    #print args.fasta2
    #print args.fastq1
    #print args.fastq2
    #print args.output   

    if args.fasta2 == None:
        args.fasta2 = args.fasta1

    [ dirname1, fq1 ] = os.path.split( args.fastq1 )
    [ dirname2, fq2 ] = os.path.split( args.fastq2 )

    if args.output == None:
        args.output = dirname1
    
    if fq1[-3:] == ".gz":
        fq1 = fq1[:-3]
    if fq2[-3:] == ".gz":
        fq2 = fq2[:-3]
    
    # make output folder
    if os.path.exists( args.output ) == False:
        os.makedirs( args.output )
    
    temp_dir = os.path.join( args.output, "tmp" )
    if os.path.exists( temp_dir ) == False:
        os.makedirs( temp_dir )
    
    fq1 = os.path.join( temp_dir, fq1 )
    fq2 = os.path.join( temp_dir, fq2 )
    fa1 = fq1 + ".fa"
    fa2 = fq2 + ".fa"

    # remove 5' end adaptor sequences in both fastq file; remove reads if any of read does not have an adaptor sequence ( Temporarily now using fastq files in Friedrich folder which are already trimmed )
    # -m 15 : Discard trimmed reads that are shorter than LENGTH.
    # --discard-untrimmed : Discard reads that do not contain the adapter.
    #
    #	bait      --CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT
    #	prey      GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT
    #	             ***  ** **** ****************************************
    #
    #cmd = "cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/$1/$6/$2.fastq -p output/$1/$6/$3.fastq -m 15 --discard-untrimmed ./fastq/$2.fastq.gz ./fastq/$3.fastq.gz"
    cmd = "cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o %s -p %s -m 15 --discard-untrimmed %s %s" % ( fq1, fq2, args.fastq1, args.fastq2 )
    run_cmd( cmd )
    
    '''
    # convert fastq to fasta ( Temporarily now using fastq file generated in Friedrich folder; since it is the same fastq. But We need to change it to Blastn for general cases )
    cmd = "python main.py fastq_to_fasta %s > %s" % ( fq1, fa1 )
    run_cmd( cmd )
    
    cmd = "python main.py fastq_to_fasta %s > %s" % ( fq2, fa2 )
    run_cmd( cmd )
    
    # blastn-short search
    cmd = "blastn -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > %s.blastn" % ( args.fasta1, fa1, fa1 )
    run_cmd( cmd )

    cmd = "blastn -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > %s.blastn" % ( args.fasta2, fa2, fa2 )
    run_cmd( cmd )
    '''

    # multi-threading
    #th1 = thread.start_new_thread( align_subprocess, ( args.fasta1, fa1, fq1 ) )
    th1 = threading.Thread(target=align_subprocess, args = ( args.fasta1, fa1, fq1 ) )
    th1.start()

    #th2 = thread.start_new_thread( align_subprocess, ( args.fasta2, fa2, fq2 ) )
    th2 = threading.Thread(target=align_subprocess, args = ( args.fasta2, fa2, fq2 ) )
    th2.start()

    th1.join()
    th2.join()

    # python parse blastn output and make ppi map
    # currently stringent case, no restriction for position cases
    # maybe ignoring orientation could be added in the future
    if args.relaxed == True:
        # no restriction for aligned position
        cmd = "main.py BLASTN_RELAXED %s %s.blastn %s.blastn > %s/%s.ppi.txt" % ( args.fasta1, fa1, fa2, args.output, args.name )
        run_cmd( cmd )  
    else:
        # very stringent case
        cmd = "main.py BLASTN_NEW %s %s.blastn %s.blastn > %s/%s.ppi.txt" % ( args.fasta1, fa1, fa2, args.output, args.name )
        run_cmd( cmd )