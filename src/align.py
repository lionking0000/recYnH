#!/usr/bin/python
import os
import threading
import fasta

VERBOSE = False # True

NUM_THREADS = 1 #
BLASTN_E_VALUE = "1e-8"

def run_cmd( cmd ):
    if VERBOSE: print cmd
    os.system( cmd )

def FastqToFasta( fastq, fasta ):
    '''
    @M03766:53:000000000-B63MG:1:1101:12798:1849 1:N:0:4
    AACGTAAAATGATATAAATATCAATATATTAAATTATATTTTGCATAAAAAACAGTCTACATAATACTGTAAATCACAACATATCCTGTCACT
    +
    AFGFHHHFHHFHH5FDGBGEGGBD55B4FG444BFG444FFF434B43433///?43BBD4B4?F44FGH44B443B33?/B1B?22?B1?12
    '''
    f = open( fastq )
    fout = open( fasta, "w" )
    for line in f.xreadlines():
        print >> fout, ">%s" % line[:-1]
        print >> fout, f.next()[:-1]
        f.next()
        f.next()
    fout.close()
    f.close()

def GenerateLastNts( fasta_file, length = 150 ):
    output_file = fasta_file
    output_file += ".-%d" % length
    fa = fasta.read_fasta( fasta_file )

    lastXnt_dic = {}

    fout = open( output_file, "w" )
    for id in fa:
        lastXnt = fa[id][-length:].upper()
        print >> fout, ">%s\n%s" % ( id, lastXnt ) #, lastXnt in lastXnt_dic
        lastXnt_dic[ lastXnt ] = id
    fout.close()

    return output_file

def CheckBlastnDB( fasta_file, lastnt_length ):
    print "[ Checking BlastN Database ]", fasta_file, lastnt_length
    bChecked = True
    if os.path.exists( fasta_file + ".-%d.nhr" % lastnt_length ) == False:
        bChecked = False
    if os.path.exists( fasta_file + ".-%d.nin" % lastnt_length ) == False:
        bChecked = False
    if os.path.exists( fasta_file + ".-%d.nsq" % lastnt_length ) == False:
        bChecked = False

    return bChecked

def MakeBlastnDB( fasta_file, lastnt_length ):
    print "[ Making BlastN Database ]", fasta_file, lastnt_length
    command = "makeblastdb -in %s.-%d -dbtype nucl" % ( fasta_file, lastnt_length )
    run_cmd( command )


def BLASTN_NEW( fasta_file, lastnt_length, filepath1, filepath2, output_file ):
    #print ( fasta_file, filepath1, filepath2, output_file )
    '''
        parse blastn output and make ppi
    '''

    ppi_cnt_dic = {}
    read1_cnt_dic = {}
    read2_cnt_dic = {}
    read1_dic = {}
    read2_dic = {}
    '''
    ['M03766:33:000000000-AT3T3:1:1101:21081:6509', '113', 'NMI', '880', '36', '54S41M', 'IGF2', '452', '0', 'ATTTTGATCATATGACTGCTCTGTTTCATTTTTTTCAATAAACCCTTTACAATTAAGTGTTCTCTAGGTCAACCTCACATAGCATACTTTGAAGA', 'HHFFHHHHFDHHHGHHHHHHHHEHHHHHGGHFHGBHGHHGHHEG4GHHHHHHHHHHHHHFFFG3GEBGBFHHHHGHHHHHGHHFHFHGHHGHHHH', 'AS:i:82', 'XN:i:0', 'XM:i:0', 'XO:i:0', 'XG:i:0', 'NM:i:0', 'MD:Z:41', 'YS:i:174', 'YT:Z:DP']
    ['M03766:33:000000000-AT3T3:1:1101:21081:6509', '177', 'IGF2', '452', '36', '5S87M', 'NMI', '880', '0', 'TCTCTAGGCCAAACGTCACCGTCCCCTGATTGCTCTACCCACCCAAGACCCCGCCCACGGGGGCGCCCCCCCAGAGATGGCCAGCAATCGGA', '/BBB/BBBFFFEFFFEEFAFB?FFFFBFFFFFFFEB;@-DFFFFFFD@FFFFEFFFFAFFFFDAFGCGGHGGHHHHHHHFFHHHGFEGFHHH', 'AS:i:174', 'XN:i:0', 'XM:i:0', 'XO:i:0', 'XG:i:0', 'NM:i:0', 'MD:Z:87', 'YS:i:82', 'YT:Z:DP']^C
    '''
    #if len( sys.argv ) < 2:
    #    print "python SAM.py ../data/roth2016_control_set_plus_control.fa output/2016-12-22_MiSeq/Friedrich/17543_S1.sam"
    #    sys.exit(0)
    total_cnt = 0
    #fa = fasta.read_fasta_file( sys.argv[2] )
    #filepath1 = sys.argv[3]
    #filepath2 = sys.argv[4]
    fa = fasta.read_fasta_file( "%s.-%d" % ( fasta_file, lastnt_length ) )
    #filepath1 = sys.argv[3]
    #filepath2 = sys.argv[4]


    PREV_QNAME = ""
    f = open( filepath1 )
    read_cnt = 0
    for line in f.xreadlines():
        #if read_cnt % 10000 == 0: print read_cnt
        read_cnt += 1
        ## READ 1
        # @M03766:53:000000000-B63MG:1:1101:13982:1738  cask_p142       98.969  97      1       0       1       97      99      3       3.06e-50        184
        [ QNAME, TARGET, PERCENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE ] = line[:-1].split("\t")
        if QNAME == PREV_QNAME: continue
        if int( SEND ) > int ( SSTART ): continue # DEFAULT CONDITION. ==> 1968532 for 2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt
        if int( QSTART ) > 10 or int( SSTART ) < 90: continue ## NEW CONDITION. 2017-10-13 => 1878423
        #if int( QSTART ) >= 5 or int( SSTART ) <= 95: continue ## NEW3 CONDITION. 2017-10-21
        #if int( QSTART ) + int( SSTART ) < 90: continue # NEW4
        read1_dic[ QNAME ] = TARGET
        read1_cnt_dic[ TARGET ] = read1_cnt_dic.get( TARGET, 0 ) + 1
        PREV_QNAME = QNAME
    f.close()

    PREV_QNAME = ""
    read_cnt = 0
    f = open( filepath2 )
    for line in f.xreadlines():
        #if read_cnt % 10000 == 0: print read_cnt
        read_cnt += 1
        ## READ 2
        # @M03766:53:000000000-B63MG:1:1101:13982:1738  cask_p142       98.969  97      1       0       1       97      99      3       3.06e-50        184
        [ QNAME, TARGET, PERCENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE ] = line[:-1].split("\t")
        if QNAME == PREV_QNAME: continue
        if int( SEND ) > int ( SSTART ): continue # DEFAULT CONDITION. ==> 1968532 for 2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt
        if int( QSTART ) > 10 or int( SSTART ) < 90: continue ## NEW CONDITION. 2017-10-13
        #if int( QSTART ) >= 5 or int( SSTART ) <= 95: continue ## NEW3 CONDITION. 2017-10-21
        #if int( QSTART ) + int( SSTART ) < 90: continue # NEW4
        read2_dic[ QNAME ] = TARGET
        read2_cnt_dic[ TARGET ] = read2_cnt_dic.get( TARGET, 0 ) + 1
        PREV_QNAME = QNAME
    f.close()

    for QNAME in read1_dic:
        TARGET2 = read2_dic.get( QNAME, "" )
        if TARGET2 == "": continue
        TARGET1 = read1_dic[ QNAME ]
        ppi_cnt_dic[ ( TARGET1, TARGET2 ) ] = ppi_cnt_dic.get( ( TARGET1, TARGET2 ), 0 ) + 1
        total_cnt += 1
    f.close()

    id_list = fa.keys()
    id_list.sort()

    fout = open( output_file, "w" )
    print >> fout, "# This file is generated by BLASTN_NEW"
    print >> fout, "DB(Read 1) \ AD(Read 2)\t"+"\t".join( id_list )
    for id1 in id_list:
        output = id1
        for id2 in id_list:
            cnt = ppi_cnt_dic.get( (id1, id2 ), 0 )
            output += "\t%d" % cnt
        print >> fout, output
    fout.close()
    # print total_cnt # only for debug

    # output/$1/Blastn/$2.blastn
    fout1 = open( "%s.cnt.txt" % filepath1, "w" )
    fout2 = open( "%s.cnt.txt" % filepath2, "w" )
    for id in id_list:
        print >> fout1, "%s\t%d" % ( id, read1_cnt_dic.get( id, 0 ) )
        print >> fout2, "%s\t%d" % ( id, read2_cnt_dic.get( id, 0 ) )
    fout1.close()
    fout2.close()

def align_subprocess( original_fasta, lastnt_length, fasta_file, fastq_file  ):
    # convert fastq to fasta ( Temporarily now using fastq file generated in Friedrich folder; since it is the same fastq. But We need to change it to Blastn for general cases )
    #cmd = "python main.py fastq_to_fasta %s > %s" % ( fastq, fasta )
    #cmd = "main.py fastq_to_fasta %s > %s" % ( fastq, fasta )
    #run_cmd( cmd )
    FastqToFasta( fastq_file, fasta_file )
    
    # blastn-short search                                                                         (20) 
    cmd = "blastn -num_threads %d -db %s.-%d  -query %s -task blastn-short -outfmt 6 -max_target_seqs 5 -evalue %s > %s.blastn" % ( NUM_THREADS, original_fasta, lastnt_length, fasta_file, BLASTN_E_VALUE, fasta_file )
    run_cmd( cmd )
   


def run( args ):
    if args.fasta2 == None:
        args.fasta2 = args.fasta1

    try:
        args.lastnt1 = int( args.lastnt1 )
        args.lastnt2 = int( args.lastnt2 )
    except:
        exit(0)

    # Blastn check
    if CheckBlastnDB( args.fasta1, args.lastnt1 ) == False:
        GenerateLastNts( args.fasta1, args.lastnt1 )
        MakeBlastnDB( args.fasta1, args.lastnt1 )
    if CheckBlastnDB( args.fasta2, args.lastnt2 ) == False:
        GenerateLastNts( args.fasta2, args.lastnt2 )
        MakeBlastnDB( args.fasta2, args.lastnt2 )

    [ dirname1, fq1 ] = os.path.split( args.fastq1 )
    [ dirname2, fq2 ] = os.path.split( args.fastq2 )

    if args.output == None:
        args.output = dirname1

    #if os.path.exists( args.output ) == None:
    #    sys.exit(0)
    #    #args.output = dirname1
    
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
    print "[ Removing adaptor sequences ]", fq1, fq2 
    cmd = "cutadapt --quiet -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o %s -p %s -m 15 --discard-untrimmed %s %s" % ( fq1, fq2, args.fastq1, args.fastq2 )
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
    print "[ Mapping sequences into reference baits and preys sequences ]", fq1, fq2 
    th1 = threading.Thread(target=align_subprocess, args = ( args.fasta1, args.lastnt1, fa1, fq1 ) )
    th1.start()

    th2 = threading.Thread(target=align_subprocess, args = ( args.fasta2, args.lastnt2, fa2, fq2 ) )
    th2.start()

    th1.join()
    th2.join()

    # python parse blastn output and make ppi map
    # currently stringent case, no restriction for position cases
    # maybe ignoring orientation could be added in the future
    print "[ Making a pair matrix ]", os.path.join( args.output, args.name ) 
    if args.relaxed == True:
        # no restriction for aligned position
        cmd = "main.py BLASTN_RELAXED %s.-%d %s.blastn %s.blastn > %s/%s" % ( args.fasta1, args.lastnt1, fa1, fa2, args.output, args.name )
        run_cmd( cmd )  
    else:
        # very stringent case
        #cmd = "main.py BLASTN_NEW %s %s.blastn %s.blastn > %s/%s" % ( args.fasta1, fa1, fa2, args.output, args.name )
        #run_cmd( cmd )
        BLASTN_NEW( args.fasta1, args.lastnt1, "%s.blastn" % fa1, "%s.blastn" % fa2, "%s/%s" % ( args.output, args.name ) )
