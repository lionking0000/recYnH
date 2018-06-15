#!/opt/conda/bin/python2.7
'''
 *
 *  recY2H CommandLine v1.0
 *
 *  Created by Jae-Seong Yang on 01/01/18.
 *  Copyright 2018 CRG. All rights reserved.
 *
 *   
'''

import fasta
#import sys
#import time
#from jslib import object
#from jslib import util
#from exp import Experiment
import command_center
#from jslib import fileutil

#======================================================================================
# This program calculate DAMRatio from the ELMSeq data
# If it is necessary, user should change corresponding file paths or prefix sequence
#======================================================================================

# User needs to setup output file folder path
OUTPUT_FOLDER = "./output"

# User needs to setup corresponding fastq file path
# Example sequencing data for transcription study (promoter)
# For EXPERIMENT_TYPE = "promoter"  
NO_SELECT_FASTQ_PATH = "./data/dam_screen_lib_2_12466_CGTGAT.fastq.gz"
SELECT_FASTQ_PATH = "./data/dam_screen_lib_2_12466_ACATCG.fastq.gz"

# User needs to specify which experiment it is among followings: ( promoter / utr_with_strong_promoter / utr_with_weak_promoter )
EXPERIMENT_TYPE = "promoter" # "promoter", "utr_with_strong_promoter", "utr_with_weak_promoter" 
READ_COUNT_CUT_OFF = 100
ALPHA = 1.0 # Used when calculating CPM

# If it is necessary, user needs to change prefix sequence filter
PROMOTER_RANDOM_SEQ_LEN = 38
PROMOTER_SEQ_PATTERN = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)GACCGGAACTTCTATGATCGAGATCGAGATCGAGATCGCGGCCGCAAC\w*"
USE_PROMOTER_TATTAT_FILTER = True  # True or False

UTR_RANDOM_SEQ_LEN = 25
UTR_STRONG_PROMOTER_SEQ_PATTERN = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)AGTTTATATTATAACACTTTAACCTATGGC\w+" 
UTR_WEAK_PROMOTER_SEQ_PATTERN = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)TGCAATTATTCTAACAAACCCCAAACTTATTTCAA\w+"

WRITE_ALL_DAMRATIOS = False  # True or False

#======================================================================================

import os
#import sh
import subprocess
import time
import re
#import psutil


def read_summary_file( input_filepath="./output/exact_match/summary/Translation_Weak_Promoter_llmp200.result.txt" ):
    '''
    # uncut_all_cnt	uncut_match_cnt	DpnI_all_cnt	DpnI_match_cnt	MboI_all_cnt	MboI_match_cnt
    # 86070864	14926650	95232908	16485908	81139744	13929837
    # DAMRatio	motif	uncut_cnt	DpnI_cnt	MboI_cnt	uncut_CPM	DpnI_CPM	MboI_CPM
    5.562432e+02	AGTAGATTTTTCGATTACTCAATTTTATAATCATTTAA	51	0	469	3.483702	0.060658	33.740524
    4.852334e+02	TCATGCAACAGATAAAGGCTGTTGTTATAATTAACAAT	42	0	409	2.880754	0.060658	29.433223
    4.497285e+02	GATGAAAGAAAAATTGGTTTTATAATATAATATTTTTA	42	0	379	2.880754	0.060658	27.279573
    '''
    data = []
    f = open( input_filepath )
    title = f.next()
    total_cnt = f.next()
    [ uncut_all_cnt, uncut_match_cnt, DpnI_all_cnt, DpnI_match_cnt, MboI_all_cnt, MboI_match_cnt ] = total_cnt[2:].split()

    cnt_infos = [ int(uncut_all_cnt), int(uncut_match_cnt), int(DpnI_all_cnt), int(DpnI_match_cnt), int(MboI_all_cnt), int(MboI_match_cnt) ]

    for line in f.xreadlines():
        if line[0] == "#": continue
        [ DAMRatio, motif, uncut_cnt, DpnI_cnt, MboI_cnt, uncut_CPM, DpnI_CPM, MboI_CPM ] = line[:-1].split( "\t" )
        if "N" in motif: continue
        if motif[25:25+6] != "TATAAT": continue
        data.append( [ float(DAMRatio), motif, int(uncut_cnt), int(DpnI_cnt), int(MboI_cnt), float(uncut_CPM), float(DpnI_CPM), float(MboI_CPM) ] )
    f.close()

    print len( data )
    return title, cnt_infos, data

'''
def MonitorMemory():
    process = psutil.Process(os.getpid())
    print(process.memory_info().rss/1024.0/1024.0)
'''

def getCommandOutput(command, exedir = None):
    if exedir != None:
        os.chdir( exedir )
    task = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    return task

def ReadSequenceFile( filepath ):
    if filepath[-3:] == ".gz":
        cmd = "zcat %s" % ( filepath )
    else:
        cmd = "cat %s" % ( filepath )
    task = getCommandOutput( cmd )
    return task


def revcomp(dna, reverse):
    bases = 'ATGCNTACGN'
    complement_dict = {bases[i]:bases[i+5] for i in range(5)}
    if reverse:
        dna = reversed(dna)
    result = [complement_dict[base] for base in dna]
    return ''.join(result)

def sort_by_cnt( aDic, _reverse ):
    output = []
    for key in aDic:
        output.append( [ aDic[ key ], key ] )
    output.sort( reverse=_reverse )
    return output


def match_cnt( ref, query ):
    match = 0
    for i in range(len(ref)):
        if ref[i] == query[i]: match += 1
    return match

def motif_cnt( f, pattern = "\w+TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)AGTTTATATTATAACACTTTAACCTATGGC\w+" ):
    motif_cnt_dic = {}
    all_cnt = 0
    match_cnt = 0

    for seq in f.xreadlines():
        all_cnt += 1
        m = re.match( pattern, seq )
        if m == None: continue
        match_cnt += 1
        motif_rc = m.group( "motif" )
        motif = revcomp( motif_rc, True )
        motif_cnt_dic[ motif ] = motif_cnt_dic.get( motif, 0 ) + 1

    return motif_cnt_dic, match_cnt, all_cnt


def ReadSequenceFiles( output_folder, uncut_fastq_path, dpnI_fastq_path, mboI_fastq_path, pattern = "TGCCCACTTCAAAAAAGCGCGATTTTTCTTCAT(?P<motif>\w+)GACCGGAACTTCTATGATCGAGATCGAGATCGAGATCGCGGCCGCAAC\w*" ):
    uncut_motif_cnt_dic = {}
    DpnI_motif_cnt_dic = {}
    MboI_motif_cnt_dic = {}

    # 1. Read No enzyme treated file
    print "[ Read No Enzyme Treated Sequence File ]", uncut_fastq_path
    fastq = ReadSequenceFile( uncut_fastq_path )
    [ uncut_motif_cnt_dic, uncut_match_cnt, uncut_all_cnt ] = motif_cnt( fastq.stdout, pattern )
    fastq.stdout.close()
    fastq.kill()

    # 2. Read DpnI treated file
    print "[ Read DpnI Treated Sequence File ]", dpnI_fastq_path
    fastq = ReadSequenceFile( dpnI_fastq_path )
    [ DpnI_motif_cnt_dic, DpnI_match_cnt, DpnI_all_cnt ] = motif_cnt( fastq.stdout, pattern )
    fastq.stdout.close()
    fastq.kill()

    # 3. Read MboI treated file
    print "[ Read MboI Treated Sequence File ]", mboI_fastq_path
    fastq = ReadSequenceFile( mboI_fastq_path )
    [ MboI_motif_cnt_dic, MboI_match_cnt, MboI_all_cnt ] = motif_cnt( fastq.stdout, pattern )
    fastq.stdout.close()
    fastq.kill()

    return uncut_motif_cnt_dic, uncut_match_cnt, uncut_all_cnt, DpnI_motif_cnt_dic, DpnI_match_cnt, DpnI_all_cnt, MboI_motif_cnt_dic, MboI_match_cnt, MboI_all_cnt


def CalculateDAMRatio( experiment_type, output_folder, uncut_fastq_path, dpnI_fastq_path, mboI_fastq_path ):
    # User needs to specify which experiment it is among followings: ( promoter / utr_with_strong_promoter / utr_with_weak_promoter )
    print "[ Starting ELMSeq Analysis ]", "It might take a long time (more than hours) depending on sequencing size"

    assert( experiment_type in [ "promoter", "utr_with_strong_promoter", "utr_with_weak_promoter" ] )
    if experiment_type == "promoter": 
        pattern = PROMOTER_SEQ_PATTERN
    if experiment_type == "utr_with_strong_promoter": 
        pattern = UTR_STRONG_PROMOTER_SEQ_PATTERN
    if experiment_type == "utr_with_weak_promoter": 
        pattern = UTR_WEAK_PROMOTER_SEQ_PATTERN

    [ uncut_motif_cnt_dic, uncut_match_cnt, uncut_all_cnt, DpnI_motif_cnt_dic, DpnI_match_cnt, DpnI_all_cnt, MboI_motif_cnt_dic, MboI_match_cnt, MboI_all_cnt ] = ReadSequenceFiles( output_folder, uncut_fastq_path, dpnI_fastq_path, mboI_fastq_path, pattern )
    _motifs = set( uncut_motif_cnt_dic.keys() )
    _motifs.update( DpnI_motif_cnt_dic.keys() )
    _motifs.update( MboI_motif_cnt_dic.keys() )

    print "[ Calculating DAMRatios ]"

    ################################################
    # Random sequence length and TATAAT motif filter
    ################################################
    motifs = set()
    if experiment_type == "promoter": 
        for motif in _motifs:
            if len( motif ) != PROMOTER_RANDOM_SEQ_LEN: continue
            if USE_PROMOTER_TATTAT_FILTER and motif[25:25+6] != "TATAAT": continue
            motifs.add( motif )
    else:
        for motif in _motifs:
            if len( motif ) != UTR_RANDOM_SEQ_LEN: continue
            motifs.add( motif )

    ################################################
    # Calculate CPM and DAMRatio
    ################################################
    data_all = []
    for motif in motifs:
        uncut_cnt = uncut_motif_cnt_dic.get( motif, 0 )
        DpnI_cnt = DpnI_motif_cnt_dic.get( motif, 0 )
        MboI_cnt = MboI_motif_cnt_dic.get( motif, 0 )
        CPM_uncut = ( float( uncut_cnt ) + ALPHA ) / uncut_match_cnt * (10**6)
        CPM_DpnI = ( float( DpnI_cnt ) + ALPHA) / DpnI_match_cnt * (10**6)
        CPM_MboI = ( float( MboI_cnt ) + ALPHA ) / MboI_match_cnt * (10**6)
        data_all.append( [ CPM_MboI / CPM_DpnI, motif, uncut_cnt, DpnI_cnt, MboI_cnt, CPM_uncut, CPM_DpnI, CPM_MboI ] )
    data_all.sort( reverse = True )

    ################################################
    # Write output files 
    ################################################
    output_all_filepath = "%s/%s.all.txt" % ( output_folder, experiment_type )
    output_filepath = "%s/%s.txt" % ( output_folder, experiment_type )
    fout = open( output_filepath, "w" )
    if WRITE_ALL_DAMRATIOS: fout_all = open( output_all_filepath, "w" )
    print >> fout, "# uncut_all_cnt\tuncut_match_cnt\tDpnI_all_cnt\tDpnI_match_cnt\tMboI_all_cnt\tMboI_match_cnt"
    print >> fout, "# %d\t%d\t%d\t%d\t%d\t%d" % ( uncut_all_cnt, uncut_match_cnt, DpnI_all_cnt, DpnI_match_cnt, MboI_all_cnt, MboI_match_cnt )
    if WRITE_ALL_DAMRATIOS: print >> fout_all, "# uncut_all_cnt\tuncut_match_cnt\tDpnI_all_cnt\tDpnI_match_cnt\tMboI_all_cnt\tMboI_match_cnt"
    if WRITE_ALL_DAMRATIOS: print >> fout_all, "# %d\t%d\t%d\t%d\t%d\t%d" % ( uncut_all_cnt, uncut_match_cnt, DpnI_all_cnt, DpnI_match_cnt, MboI_all_cnt, MboI_match_cnt )
    print >> fout, "# DAMRatio\tmotif\tuncut_cnt\tDpnI_cnt\tMboI_cnt\tuncut_CPM\tDpnI_CPM\tMboI_CPM"
    if WRITE_ALL_DAMRATIOS: print >> fout_all, "# DAMRatio\tmotif\tuncut_cnt\tDpnI_cnt\tMboI_cnt\tuncut_CPM\tDpnI_CPM\tMboI_CPM"
    for DAMRatio, motif, uncut_cnt, DpnI_cnt, MboI_cnt, CPM_uncut, CPM_DpnI, CPM_MboI in data_all:
        if WRITE_ALL_DAMRATIOS: print >> fout_all, "\t".join( [ "%e"%DAMRatio, motif, "%d"%uncut_cnt, "%d"%DpnI_cnt, "%d"%MboI_cnt, "%f"%CPM_uncut, "%f"%CPM_DpnI, "%f"%CPM_MboI ] )
        if uncut_cnt >= READ_COUNT_CUT_OFF or DpnI_cnt >= READ_COUNT_CUT_OFF or MboI_cnt >= READ_COUNT_CUT_OFF:
            print >> fout, "\t".join( [ "%e"%DAMRatio, motif, "%d"%uncut_cnt, "%d"%DpnI_cnt, "%d"%MboI_cnt, "%f"%CPM_uncut, "%f"%CPM_DpnI, "%f"%CPM_MboI ] )
    fout.close()
    if WRITE_ALL_DAMRATIOS: fout_all.close()

    print "[ Finishing ELMSeq Analysis ]"

def GenerateLastNts( blast_db_path, filepath, length = 30 ):
    fa = fasta.read_fasta( filepath )

    lastXnt_dic = {}

    fout = open( blast_db_path, "w" )
    for id in fa:
        lastXnt = fa[id][-length:].upper()
        print >> fout, ">%s\n%s" % ( id, lastXnt ) #, lastXnt in lastXnt_dic
        lastXnt_dic[ lastXnt ] = id
    fout.close()

    return lastXnt_dic

def ReadLastNts( blast_db_path, length = 30 ):
    fa = fasta.read_fasta( blast_db_path )

    lastXnt_dic = {}

    for id in fa:
        lastXnt = fa[id][-length:].upper()
        lastXnt_dic[ lastXnt ] = id

    return lastXnt_dic

def FastqToFasta( fastq, fasta ):
    '''
    @M03766:53:000000000-B63MG:1:1101:12798:1849 1:N:0:4
    AACGTAAAATGATATAAATATCAATATATTAAATTATATTTTGCATAAAAAACAGTCTACATAATACTGTAAATCACAACATATCCTGTCACT
    +
    AFGFHHHFHHFHH5FDGBGEGGBD55B4FG444BFG444FFF434B43433///?43BBD4B4?F44FGH44B443B33?/B1B?22?B1?12
    '''
    fout = open( fasta, "w" )
    f = open( fastq )
    for line in f.xreadlines():
        print >> fout, ">%s" % line[:-1]
        print >> fout, f.next()[:-1]
        f.next()
        f.next()
    f.close()
    fout.close()

def BLASTN( input_fasta, noselect_blastn_read1_name, noselect_blastn_read2_name, noselect_ppi_file ):  
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
    fa = fasta.read_fasta( input_fasta )
    filepath1 = noselect_blastn_read1_name
    filepath2 = noselect_blastn_read2_name

    PREV_QNAME = ""
    f = open( filepath1 )
    read_cnt = 0
    for line in f.xreadlines():
        #if read_cnt % 10000 == 0: print read_cnt
        read_cnt += 1
        ## READ 1
        # @M03766:53:000000000-B63MG:1:1101:13982:1738	cask_p142	98.969	97	1	0	1	97	99	3	3.06e-50	184
        [ QNAME, TARGET, PERCENT, LENGTH, MISMATCH, GAPOPEN, QSTART, QEND, SSTART, SEND, EVALUE, BITSCORE ] = line[:-1].split("\t")
        if QNAME == PREV_QNAME: continue
        if int( SEND ) > int ( SSTART ): continue # DEFAULT CONDITION. ==> 1968532 for 2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt
        #if int( QSTART ) > 10 or int( SSTART ) < 90: continue ## NEW CONDITION. 2017-10-13 => 1878423
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
        #if int( QSTART ) > 10 or int( SSTART ) < 90: continue ## NEW CONDITION. 2017-10-13
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

    fout = open( noselect_ppi_file, "w" )
    print >> fout, "# This file is generated by recY2H"
    print >> fout, "DB(Read 1) \ AD(Read 2)\t"+"\t".join( id_list )
    for id1 in id_list:
        output = id1
        for id2 in id_list:
            cnt = ppi_cnt_dic.get( (id1, id2 ), 0 )
            output += "\t%d" % cnt
        print >> fout, output
    fout.close()

    # print total_cnt # only for debug
    #noselect_ppi_file

    # output/$1/Blastn/$2.blastn
    fout1 = open( "%s.cnt.txt" % filepath1, "w" )
    fout2 = open( "%s.cnt.txt" % filepath2, "w" )
    for id in id_list:
        print >> fout1, "%s\t%d" % ( id, read1_cnt_dic.get( id, 0 ) )
        print >> fout2, "%s\t%d" % ( id, read2_cnt_dic.get( id, 0 ) )
    fout1.close()
    fout2.close()


def CheckBlastDB( blast_db_path ):
    if os.path.exists( "%s.nhr" % blast_db_path ) == False: return False
    if os.path.exists( "%s.nin" % blast_db_path ) == False: return False
    if os.path.exists( "%s.nsq" % blast_db_path ) == False: return False
    return True
    
# CalculateRecY2H( EXPERIMENT_TYPE, OUTPUT_FOLDER, NOSELCT_FASTQ_PATH, SELECT_FASTQ_PATH )
def CalculateRecY2H( argv ):
    # experiment_type, output_folder, noselect_fastq_path, select_fastq_path, input_fasta
    input_fasta = "./data/roth2016_control_set_plus_control.fa"
    nt_len = 150
    blast_db_path = "./data/roth2016_control_set_plus_control.-150.fa"
    output_folder = "./output/RothExp1"
    noselect_fastq_path_read1 = "../../2017-03-03_MiSeq/S1_R1.fastq.gz"
    noselect_fastq_path_read2 = "../../2017-03-03_MiSeq/S1_R2.fastq.gz"
    select_fastq_path_read1 = "../../2017-03-03_MiSeq/S2_R1.fastq.gz"
    select_fastq_path_read2 = "../../2017-03-03_MiSeq/S2_R2.fastq.gz"
    noselect_ppi_file = os.path.join( output_folder, "nonselect.ppi.txt" )
    select_ppi_file = os.path.join( output_folder, "select.ppi.txt" )
    
    
    # 1. Make Blast DB
    # 1.1 Generate last 150nt db 
    # python main.py generate_last_nts ../data/A463-MGj69.RBP-MAP.fa 150 > ../data/A463-MGj69.RBP-MAP.-150.fa 
    if os.path.exists( blast_db_path ) == False:
        lastXnt_dic = GenerateLastNts( blast_db_path, input_fasta, nt_len )
    else:
        lastXnt_dic = ReadLastNts( blast_db_path, nt_len )


    # 1.2 Blast DB
    # makeblastdb -in ../data/A463-MGj69.RBP-MAP.-150.fa -dbtype nucl
    if CheckBlastDB( blast_db_path ) == False:
        command = "makeblastdb -in %s -dbtype nucl" % blast_db_path
        os.system( command )
    
    # 2. Run Short Blast
    # 2.1 Make output folder
    if os.path.exists( output_folder ) == False:
        os.makedirs( output_folder )

    # 2.2 Remove Adaptor
    noselect_fastq_read1_name = os.path.join( output_folder,"1.fastq" ) #os.path.split( noselect_fastq_path_read1 )[-1]
    noselect_fastq_read2_name = os.path.join( output_folder,"2.fastq" ) #os.path.split( noselect_fastq_path_read2 )[-1]
    select_fastq_read1_name = os.path.join( output_folder,"3.fastq" ) #os.path.split( select_fastq_path_read1 )[-1]
    select_fastq_read2_name = os.path.join( output_folder,"4.fastq" ) #os.path.split( select_fastq_path_read2 )[-1]
    
    Adaptor1 = "CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT"
    Adaptor2 = "GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT"

    if os.path.exists( noselect_fastq_read1_name ) == False or os.path.exists( noselect_fastq_read2_name ) == False:
        print "[ Trimming Adapter sequences for non-selection - it will take some time. ]"
        command = "cutadapt -g %s -G %s -o %s -p %s -m 15 --discard-untrimmed %s %s" % ( Adaptor1, Adaptor2, noselect_fastq_read1_name, noselect_fastq_read2_name, noselect_fastq_path_read1, noselect_fastq_path_read2 )

    if os.path.exists( select_fastq_read1_name ) == False or os.path.exists( select_fastq_read2_name ) == False:
        print "[ Trimming Adapter sequences for selection - it will take some time. ]"
        command = "cutadapt -g %s -G %s -o %s -p %s -m 15 --discard-untrimmed %s %s" % ( Adaptor1, Adaptor2, select_fastq_read1_name, select_fastq_read2_name, select_fastq_path_read1, select_fastq_path_read2 )
        print command
        os.system( command )

    # 2.3 Convert Fastq file into Fasta file
    #python main.py fastq_to_fasta output/2017-07-03_MiSeq/Friedrich/60_Q_R1.fastq > output/2017-07-03_MiSeq/Blastn/60_Q_R1.fa
    #python main.py fastq_to_fasta output/2017-07-03_MiSeq/Friedrich/60_Q_R2.fastq > output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa
    
    noselect_fasta_read1_name = os.path.join( output_folder,"1.fa" ) #os.path.split( noselect_fastq_path_read1 )[-1]
    noselect_fasta_read2_name = os.path.join( output_folder,"2.fa" ) #os.path.split( noselect_fastq_path_read2 )[-1]
    select_fasta_read1_name = os.path.join( output_folder,"3.fa" ) #os.path.split( select_fastq_path_read1 )[-1]
    select_fasta_read2_name = os.path.join( output_folder,"4.fa" ) #os.path.split( select_fastq_path_read2 )[-1]
    if os.path.exists( noselect_fastq_read1_name ) == False:
        print "[ Convert Fastq file into Fasta file ]"
        FastqToFasta( noselect_fastq_read1_name, noselect_fasta_read1_name )
    if os.path.exists( noselect_fastq_read2_name ) == False:
        print "[ Convert Fastq file into Fasta file ]"
        FastqToFasta( noselect_fastq_read2_name, noselect_fasta_read2_name )
    if os.path.exists( select_fastq_read1_name ) == False:
        print "[ Convert Fastq file into Fasta file ]"
        FastqToFasta( select_fastq_read1_name, select_fasta_read1_name )
    if os.path.exists( select_fastq_read2_name ) == False:
        print "[ Convert Fastq file into Fasta file ]"
        FastqToFasta( select_fastq_read2_name, select_fasta_read2_name )

    # 2.4 Run Short Blast
    noselect_blastn_read1_name = os.path.join( output_folder,"1.blastn" ) #os.path.split( noselect_fastq_path_read1 )[-1]
    noselect_blastn_read2_name = os.path.join( output_folder,"2.blastn" ) #os.path.split( noselect_fastq_path_read2 )[-1]
    select_blastn_read1_name = os.path.join( output_folder,"3.blastn" ) #os.path.split( select_fastq_path_read1 )[-1]
    select_blastn_read2_name = os.path.join( output_folder,"4.blastnfa" ) #os.path.split( select_fastq_path_read2 )[-1]
    if os.path.exists( noselect_blastn_read1_name ) == False:
        print "[ Blastn into target sequence with Fasta file 1 ]"
        command = "blastn -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > %s" % ( blast_db_path, noselect_fasta_read1_name, noselect_blastn_read1_name )
        os.system( command ) 
    if os.path.exists( noselect_blastn_read2_name ) == False:
        print "[ Blastn into target sequence with Fasta file 2 ]"
        command = "blastn -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > %s" % ( blast_db_path, noselect_fasta_read2_name, noselect_blastn_read2_name )
        os.system( command ) 
    if os.path.exists( select_blastn_read1_name ) == False:
        print "[ Blastn into target sequence with Fasta file 3 ]"
        command = "blastn -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > %s" % ( blast_db_path, select_fasta_read1_name, select_blastn_read1_name )
        os.system( command ) 
    if os.path.exists( select_blastn_read2_name ) == False:
        print "[ Blastn into target sequence with Fasta file 4 ]"
        command = "blastn -db %s  -query %s -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > %s" % ( blast_db_path, select_fasta_read2_name, select_blastn_read2_name )
        os.system( command ) 

    # 3. PPI mapping
    # python parse blastn output and make ppi map
    #noselect_ppi_file = os.path.join( output_folder, "nonselect.ppi.txt" )
    #select_ppi_file = os.path.join( output_folder, "select.ppi.txt" )

    #echo "python main.py BLASTN $4.fa output/$1/$6/$2.blastn output/$1/$6/$3.blastn > output/$1/$6/$5.ppi.txt"
    if os.path.exists( noselect_ppi_file ) == False:
        BLASTN( input_fasta, noselect_blastn_read1_name, noselect_blastn_read2_name, noselect_ppi_file )

    if os.path.exists( select_ppi_file ) == False:
        BLASTN( input_fasta, select_blastn_read1_name, select_blastn_read2_name, select_ppi_file )
    # cat output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn | awk '{print $1}' | uniq | wc    ## ==> 3209686
    # cat output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn | awk '{print $1}' | uniq | wc    ## ==> 3029067
    # python main.py BLASTN ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn    ##  ==> 2507307



'''
#===============================================================================
# Y2H_Blastn.sh
# created 2017.07.03
# written by Jae-Seong Yang
#
# Blastn short-read based ppi map
#===============================================================================

# READ: https://www.ncbi.nlm.nih.gov/books/NBK279668/

# $1 : directory
# $2 : File1 (without .fastq.gz)
# $3 : File2 (without .fastq.gz)
# $4 : Fasta (without .fa)
# $5 : Outname (prefix only)
# $6 : Output Folder

#===============================================================================
# Need to make blastn database (only -100 bp containing) as following instructions before running this script
#
# python main.py generate_last_nts ../data/P170_4_library_MGj5615-Jun-2017122014.fa 100 > ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa
# makeblastdb -in P170_4_library_MGj5615-Jun-2017122014.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/E375pLT-20170706.fa 100 > ../data/E375pLT-20170706.-100.fa 
# makeblastdb -in ../data/E375pLT-20170706.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/roth2016_control_set_plus_control.fa 100 > ../data/roth2016_control_set_plus_control.-100.fa 
# makeblastdb -in ../data/roth2016_control_set_plus_control.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/EGFR_entry_clones.fa  100 > ../data/EGFR_entry_clones.-100.fa
# makeblastdb -in EGFR_entry_clones.-100.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/roth2016_control_set_plus_control.fa 150 > ../data/roth2016_control_set_plus_control.-150.fa 
# makeblastdb -in ../data/roth2016_control_set_plus_control.-150.fa -dbtype nucl
#
# python main.py generate_last_nts ../data/A463-MGj69.RBP-MAP.fa 150 > ../data/A463-MGj69.RBP-MAP.-150.fa 
# makeblastdb -in ../data/A463-MGj69.RBP-MAP.-150.fa -dbtype nucl
#===============================================================================


#===============================================================================
# JOB_LIST (what have been done)
#
### output/2016-12-22_MiSeq/Friedrich/
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17543_S1_R1 17543_S1_R2 ../data/roth2016_control_set_plus_control.-100 17543_S1 > qjobs/qjob_2016-12-22_MiSeq_S1.sh # ==> 1968532 (Roth -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17544_S2_R1 17544_S2_R2 ../data/roth2016_control_set_plus_control.-100 17544_S2 > qjobs/qjob_2016-12-22_MiSeq_S2.sh # ==> 1541502 (Roth -A) from 2418986 (trimed seq) from xxx
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17545_S3_R1 17545_S3_R2 ../data/roth2016_control_set_plus_control.-100 17545_S3 > qjobs/qjob_2016-12-22_MiSeq_S3.sh # ==> 2190445 (Roth Seaprep -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17546_S4_R1 17546_S4_R2 ../data/roth2016_control_set_plus_control.-100 17546_S4 > qjobs/qjob_2016-12-22_MiSeq_S4.sh # ==> 1877953  (Roth Seaprep -A)
#
### output/2017-02-04_MiSeq/Friedrich/ ## (not use because of contamination)
#
#
### output/2017-02-22_MiSeq/Friedrich/
#
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control.-100 S1 > qjobs/qjob_2017-02-22_MiSeq_S1.sh # (Roth -W) # ==> 960245
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control.-100 S2 > qjobs/qjob_2017-02-22_MiSeq_S2.sh # (Roth -A) # ==> 2024754
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control.-100 S3 > qjobs/qjob_2017-02-22_MiSeq_S3.sh # (Roth -Q) # ==> 813912
#
### output/2017-03-03_MiSeq/Friedrich (Roth)
#
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control.-100 S1 > qjobs/qjob_2017-03-03_MiSeq_S1.sh # (Roth -W) # ==> 1414435
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control.-100 S2 > qjobs/qjob_2017-03-03_MiSeq_S2.sh # (Roth -A) # ==> 969159
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control.-100 S3 > qjobs/qjob_2017-03-03_MiSeq_S3.sh # (Roth -Q) # ==> 1255002
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S4_R1 S4_R2 ../data/roth2016_control_set_plus_control.-100 S4 > qjobs/qjob_2017-03-03_MiSeq_S4.sh # (Roth Aonly -W) to check toxicity # ==> 23
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S5_R1 S5_R2 ../data/roth2016_control_set_plus_control.-100 S5 > qjobs/qjob_2017-03-03_MiSeq_S5.sh # (Roth Bonly -W) to check toxicity # ==> 169
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S6_R1 S6_R2 ../data/roth2016_control_set_plus_control.-100 S6 > qjobs/qjob_2017-03-03_MiSeq_S6.sh # (Roth Bonly -A) to check auto-activation # ==> 1296
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S7_R1 S7_R2 ../data/roth2016_control_set_plus_control.-100 S7 > qjobs/qjob_2017-03-03_MiSeq_S7.sh # (Roth Seaprep -Q) # ==> 734766
# 
# 2017-06-08_MiSeq (P170 toxicity test and Roth Seaprep)
#  
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S49_R1 S49_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S49 > qjobs/qjob_2017-06-08_MiSeq_S49.sh # (P170 Bonly -W) to check toxicity           # ==> 9299 
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S50_R1 S50_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S50 > qjobs/qjob_2017-06-08_MiSeq_S50.sh # (P170 Bonly -W-A) to check auto-activation  # ==> 7991
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S51_R1 S51_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S51 > qjobs/qjob_2017-06-08_MiSeq_S51.sh # (P170 Aonly -W) to check toxicity           # ==> 3965
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S52_R1 S52_R2 ../data/roth2016_control_set_plus_control.-100 S52 > qjobs/qjob_2017-06-08_MiSeq_S52.sh # (Roth Seaprep -W) to check complexcity          # ==> 2801344
#
#
# 2017-06-12_MiSeq (P170)
# sh Y2H_Blastn.sh 2017-06-12_MiSeq S53_R1 S53_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S53 > qjobs/qjob_2017-06-08_MiSeq_S53.sh # (P170 -W) no selection   # ==> 3250858
# sh Y2H_Blastn.sh 2017-06-12_MiSeq S54_R1 S54_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S54 > qjobs/qjob_2017-06-08_MiSeq_S54.sh # (P170 -W-A) selection    # ==> 3553164
#
### 2017-07-03_MiSeq : 60
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 60_W_R1 60_W_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 60_W # ==> 2629509
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 60_Q_R1 60_Q_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 60_Q # ==> 2507307
#
### 2017-07-03_MiSeq : 58
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 58_AW_R1 58_AW_R2 ../data/E375pLT-20170706.-100 58_AW # ==> 219
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 58_BW_R1 58_BW_R2 ../data/E375pLT-20170706.-100 58_BW # ==> 219
# sh Y2H_Blastn.sh 2017-07-03_MiSeq 58_BQ_R1 58_BQ_R2 ../data/E375pLT-20170706.-100 58_BQ" # ==> 306
#
# 2017-08-15_MiSeq (Test for Seaprep again)
# sh Y2H_Blastn.sh 2017-08-15_MiSeq S61_SWD_R1 S61_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S61_SWD > qjobs/qjob_2017-08-15_MiSeq_S61_SWD.sh # SWD (P170 Seaprep -W)       # ==> 2909222 out of 6012622 (48.4%)
# sh Y2H_Blastn.sh 2017-08-15_MiSeq S61_PWD_R1 S61_PWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S61_PWD > qjobs/qjob_2017-08-15_MiSeq_S61_PWD.sh # SWD (P170 Plate -W)         # ==> 1956131 out of 4473498 (43.7%)
# sh Y2H_Blastn.sh 2017-08-15_MiSeq S61_PQD_R1 S61_PQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S61_PQD > qjobs/qjob_2017-08-15_MiSeq_S61_PQD.sh # SWD (P170 Plate -Q)         # ==> 1684777 out of 3553019 (47.4%)
#
# 2017-08-22_MiSeq (Roth - Seaprep with selection)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SWD_S1_R1 S64-SWD_S1_R2 ../data/roth2016_control_set_plus_control.-100 S64_SWD > qjobs/qjob_2017-08-22_MiSeq_S64_SWD.sh # SWD (Roth Seaprep -W)  # ==> 1491017 out of 2562522 (58.2%)           ==> 1450453 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA4D_S2_R1 S64-SA4D_S2_R2 ../data/roth2016_control_set_plus_control.-100 S64_SA4D > qjobs/qjob_2017-08-22_MiSeq_S64_SA4D.sh # SWD (Roth Seaprep +A1/4) # ==> 1177404 out of 2321645 (50.7%)     ==> 1145754 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA8D_S3_R1 S64-SA8D_S3_R2 ../data/roth2016_control_set_plus_control.-100 S64_SA8D > qjobs/qjob_2017-08-22_MiSeq_S64_SA8D.sh # SWD (Roth Seaprep +A1/8) # ==> 1125083 out of 2076630 (54.2%)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SQD_S4_R1 S64-SQD_S4_R2 ../data/roth2016_control_set_plus_control.-100 S64_SQD > qjobs/qjob_2017-08-22_MiSeq_S64_SQD.sh # SWD (Roth Seaprep -Q)    # ==> 618177 out of 1437540 (43.0%)
#
# 2017-08-28_MiSeq (P170 - Seaprep)
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SWD_R1 S68_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SWD > qjobs/qjob_2017-08-28_MiSeq_S68_SWD.sh # SWD (P170 Seaprep -W)            # 1999636
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SA4D_R1 S68_SA4D_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SA4D > qjobs/qjob_2017-08-28_MiSeq_S68_SA4D.sh # SA4D (P170 Seaprep +A1/4)    # 1943064
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SA8D_R1 S68_SA8D_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SA8D > qjobs/qjob_2017-08-28_MiSeq_S68_SA8D.sh # SA8D (P170 Seaprep +A1/8)    # 1372704
# sh Y2H_Blastn.sh 2017-08-28_MiSeq S68_SQD_R1 S68_SQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S68_SQD > qjobs/qjob_2017-08-28_MiSeq_S68_SQD.sh # SQD (P170 Seaprep -Q)            # 1247668
#
# 2017-08-30_MiSeq (P170 - Seaprep)
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SWD_R1 S65_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SWD > qjobs/qjob_2017-08-30_MiSeq_S65_SWD.sh # SWD (P170 Seaprep -W)            # 1836282
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SAD_R1 S65_SAD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SAD > qjobs/qjob_2017-08-30_MiSeq_S65_SAD.sh # SAD (P170 Seaprep +A?)           # 2108754
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SQD_R1 S65_SQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SQD > qjobs/qjob_2017-08-30_MiSeq_S65_SQD.sh # SQD (P170 Seaprep -Q)            # 1870105
#
# 2017-10-09_MiSeq (EGFR - 1)
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S1_WD_R1 S1_WD_R2 ../data/EGFR_entry_clones.-100 S1_WD > qjobs/qjob_2017-10-09_MiSeq_S1_WD.sh # WD
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S2_WD_R1 S2_WD_R2 ../data/EGFR_entry_clones.-100 S2_WD > qjobs/qjob_2017-10-09_MiSeq_S2_WD.sh
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S3_A8D_R1 S3_A8D_R2 ../data/EGFR_entry_clones.-100 S3_A8D > qjobs/qjob_2017-10-09_MiSeq_S3_A8D.sh
# sh Y2H_Blastn.sh 2017-10-09_MiSeq S4_QD_R1 S4_QD_R2 ../data/EGFR_entry_clones.-100 S4_QD > qjobs/qjob_2017-10-09_MiSeq_S4_QD.sh
#
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S1_WD_R1 S1_WD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S1_WD > qjobs/qjob_2017-10-16_MiSeq_S1_WD.sh # WD (P170 Seaprep -W) 
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S2_WD_R1 S2_WD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S2_WD > qjobs/qjob_2017-10-16_MiSeq_S2_WD.sh # WD (P170 Seaprep -W) 
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S3_AD_R1 S3_AD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S3_AD > qjobs/qjob_2017-10-16_MiSeq_S3_AD.sh # WD (P170 Seaprep -W) 
# sh Y2H_Blastn.sh 2017-10-16_MiSeq S4_AD_R1 S4_AD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S4_AD > qjobs/qjob_2017-10-16_MiSeq_S4_AD.sh # WD (P170 Seaprep -W) 
#
# 2017-10-30_MiSeq (R75 auto-activator)
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S1_BWD_R1 S1_BWD_R2 ../data/roth2016_control_set_plus_control.-100 S1_BWD > qjobs/qjob_2017-10-30_MiSeq_S1_BWD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S2_BA2D_R1 S2_BA2D_R2 ../data/roth2016_control_set_plus_control.-100 S2_BA2D > qjobs/qjob_2017-10-30_MiSeq_S2_BA2D.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S3_BQD_R1 S3_BQD_R2 ../data/roth2016_control_set_plus_control.-100 S3_BQD > qjobs/qjob_2017-10-30_MiSeq_S3_BQD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S4_AWD_R1 S4_AWD_R2 ../data/roth2016_control_set_plus_control.-100 S4_AWD > qjobs/qjob_2017-10-30_MiSeq_S4_AWD.sh
#
# 2017-11-03_MiSeq (R75 technical repeat)
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S1_W_R1 S1_W_R2 ../data/roth2016_control_set_plus_control.-100 S1_W > qjobs/qjob_2017-11-03_MiSeq_S1_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S2_Q_R1 S2_Q_R2 ../data/roth2016_control_set_plus_control.-100 S2_Q > qjobs/qjob_2017-11-03_MiSeq_S2_Q.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S3_W_R1 S3_W_R2 ../data/roth2016_control_set_plus_control.-100 S3_W > qjobs/qjob_2017-11-03_MiSeq_S3_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S4_Q_R1 S4_Q_R2 ../data/roth2016_control_set_plus_control.-100 S4_Q > qjobs/qjob_2017-11-03_MiSeq_S4_Q.sh
#
#===============================================================================


#===============================================================================
# Use all reference sequences to map
### 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 ) 
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17543_S1_R1 17543_S1_R2 ../data/roth2016_control_set_plus_control 17543_S1 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S1.sh # ==> 1968532 (Roth -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17544_S2_R1 17544_S2_R2 ../data/roth2016_control_set_plus_control 17544_S2 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S2.sh # ==> 1541502 (Roth -A) from 2418986 (trimed seq) from xxx
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17545_S3_R1 17545_S3_R2 ../data/roth2016_control_set_plus_control 17545_S3 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S3.sh # ==> 2190445 (Roth Seaprep -W)
# sh Y2H_Blastn.sh 2016-12-22_MiSeq 17546_S4_R1 17546_S4_R2 ../data/roth2016_control_set_plus_control 17546_S4 Blastn_All_Ref > qjobs/all.qjob_2016-12-22_MiSeq_S4.sh # ==> 1877953  (Roth Seaprep -A)

### 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA ) 
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control S1 Blastn_All_Ref > qjobs/all.qjob_2017-02-22_MiSeq_S1.sh # (Roth -W) # ==> 960245
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control S2 Blastn_All_Ref > qjobs/all.qjob_2017-02-22_MiSeq_S2.sh # (Roth -A) # ==> 2024754
# sh Y2H_Blastn.sh 2017-02-22_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control S3 Blastn_All_Ref > qjobs/all.qjob_2017-02-22_MiSeq_S3.sh # (Roth -Q) # ==> 813912

### 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S1_R1 S1_R2 ../data/roth2016_control_set_plus_control S1 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S1.sh # (Roth -W) # ==> 1414435
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S2_R1 S2_R2 ../data/roth2016_control_set_plus_control S2 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S2.sh # (Roth -A) # ==> 969159
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S3_R1 S3_R2 ../data/roth2016_control_set_plus_control S3 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S3.sh # (Roth -Q) # ==> 1255002
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S4_R1 S4_R2 ../data/roth2016_control_set_plus_control S4 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S4.sh # (Roth Aonly -W) to check toxicity # ==> 23
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S5_R1 S5_R2 ../data/roth2016_control_set_plus_control S5 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S5.sh # (Roth Bonly -W) to check toxicity # ==> 169
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S6_R1 S6_R2 ../data/roth2016_control_set_plus_control S6 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S6.sh # (Roth Bonly -A) to check auto-activation # ==> 1296
# sh Y2H_Blastn.sh 2017-03-03_MiSeq S7_R1 S7_R2 ../data/roth2016_control_set_plus_control S7 Blastn_All_Ref > qjobs/all.qjob_2017-03-03_MiSeq_S7.sh # (Roth Seaprep -Q) # ==> 734766

### 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SWD_S1_R1 S64-SWD_S1_R2 ../data/roth2016_control_set_plus_control S64_SWD Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SWD.sh # SWD (Roth Seaprep -W)  # ==> 1491017 out of 2562522 (58.2%)           ==> 1450453 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA4D_S2_R1 S64-SA4D_S2_R2 ../data/roth2016_control_set_plus_control S64_SA4D Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SA4D.sh # SWD (Roth Seaprep +A1/4) # ==> 1177404 out of 2321645 (50.7%)     ==> 1145754 (For New)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SA8D_S3_R1 S64-SA8D_S3_R2 ../data/roth2016_control_set_plus_control S64_SA8D Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SA8D.sh # SWD (Roth Seaprep +A1/8) # ==> 1125083 out of 2076630 (54.2%)
# sh Y2H_Blastn.sh 2017-08-22_MiSeq S64-SQD_S4_R1 S64-SQD_S4_R2 ../data/roth2016_control_set_plus_control S64_SQD Blastn_All_Ref > qjobs/all.qjob_2017-08-22_MiSeq_S64_SQD.sh # SWD (Roth Seaprep -Q)    # ==> 618177 out of 1437540 (43.0%)

### 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W; 
# sh Y2H_Blastn.sh 2017-06-08_MiSeq S52_R1 S52_R2 ../data/roth2016_control_set_plus_control S52 Blastn_All_Ref > qjobs/all.qjob_2017-06-08_MiSeq_S52.sh # (Roth Seaprep -W) to check complexcity   

### 2017-10-30_MiSeq (R75 auto-activator)
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S1_BWD_R1 S1_BWD_R2 ../data/roth2016_control_set_plus_control S1_BWD Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S1_BWD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S2_BA2D_R1 S2_BA2D_R2 ../data/roth2016_control_set_plus_control S2_BA2D Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S2_BA2D.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S3_BQD_R1 S3_BQD_R2 ../data/roth2016_control_set_plus_control S3_BQD Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S3_BQD.sh
# sh Y2H_Blastn.sh 2017-10-30_MiSeq S4_AWD_R1 S4_AWD_R2 ../data/roth2016_control_set_plus_control S4_AWD Blastn_All_Ref > qjobs/all.qjob_2017-10-30_MiSeq_S4_AWD.sh

# 2017-11-03_MiSeq; Roth75 - technical repeat  
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S1_W_R1 S1_W_R2 ../data/roth2016_control_set_plus_control S1_W Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S1_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S2_Q_R1 S2_Q_R2 ../data/roth2016_control_set_plus_control S2_Q Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S2_Q.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S3_W_R1 S3_W_R2 ../data/roth2016_control_set_plus_control S3_W Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S3_W.sh
# sh Y2H_Blastn.sh 2017-11-03_MiSeq S4_Q_R1 S4_Q_R2 ../data/roth2016_control_set_plus_control S4_Q Blastn_All_Ref > qjobs/all.qjob_2017-11-03_MiSeq_S4_Q.sh   


#===============================================================================


#===============================================================================
# A463 library (all)
#sh Y2H_Blastn.sh 2017-11-17_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2017-11-17_MiSeq_S1_W2.sh
#sh Y2H_Blastn.sh 2017-11-17_MiSeq S2_QD_R1 S2_QD_R2 ../data/A463-MGj69.RBP-MAP.-150 S2_Q > qjobs/qjob_Sebastian_2017-11-17_MiSeq_S2_Q2.sh


# A463 library (w/o auto activator)
#sh Y2H_Blastn.sh 2017-12-07_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2017-12-07_MiSeq_S1_W.sh                
#sh Y2H_Blastn.sh 2017-12-07_MiSeq S2_QD_R1 S2_QD_R2 ../data/A463-MGj69.RBP-MAP.-150 S2_Q > qjobs/qjob_Sebastian_2017-12-07_MiSeq_S2_Q.sh                

# A463 library (w/o auto activator)
#sh Y2H_Blastn.sh 2018-02-21_MiSeq S1_WD_R1 S1_WD_R2 ../data/A463-MGj69.RBP-MAP.-150 S1_W > qjobs/qjob_Sebastian_2018-02-21_MiSeq_S1_W.sh                
#sh Y2H_Blastn.sh 2018-02-21_MiSeq S2_QD_R1 S2_QD_R2 ../data/A463-MGj69.RBP-MAP.-150 S2_Q > qjobs/qjob_Sebastian_2018-02-21_MiSeq_S2_Q.sh                
#===============================================================================


# make output folder
echo "cd /users/lserrano/jyang/work/Mireia/src"
#echo "mkdir output/$1/Blastn"
echo "mkdir output/$1/$6"

# remove 5' end adaptor sequences in both fastq file; remove reads if any of read does not have an adaptor sequence ( Temporarily now using fastq files in Friedrich folder which are already trimmed )
# -m 15 : Discard trimmed reads that are shorter than LENGTH.
# --discard-untrimmed : Discard reads that do not contain the adapter.
#
#	bait      --CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT
#	prey      GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT
#	             ***  ** **** ****************************************
#
echo "cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/$1/$6/$2.fastq -p output/$1/$6/$3.fastq -m 15 --discard-untrimmed ../$1/$2.fastq.gz ../$1/$3.fastq.gz"
#cutadapt -g CGCTGCAGGTCGACGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -G GCAGCTCGAGCTCGATGGATCTTAGTTACTTACCACTTTGTACAAGAAAGCTGGGT -o output/$1/Blastn/$2.fastq -p output/$1/Blastn/$3.fastq -m 15 --discard-untrimmed ../$1/$2.fastq.gz ../$1/$3.fastq.gz

# convert fastq to fasta ( Temporarily now using fastq file generated in Friedrich folder; since it is the same fastq. But We need to change it to Blastn for general cases )
#echo "python main.py fastq_to_fasta output/$1/Friedrich/$2.fastq > output/$1/Blastn/$2.fa"
#echo "python main.py fastq_to_fasta output/$1/Friedrich/$3.fastq > output/$1/Blastn/$3.fa"
echo "python main.py fastq_to_fasta output/$1/$6/$2.fastq > output/$1/$6/$2.fa"
echo "python main.py fastq_to_fasta output/$1/$6/$3.fastq > output/$1/$6/$3.fa"
#python main.py fastq_to_fasta output/2017-07-03_MiSeq/Friedrich/60_Q_R1.fastq > output/2017-07-03_MiSeq/Blastn/60_Q_R1.fa
#python main.py fastq_to_fasta output/2017-07-03_MiSeq/Friedrich/60_Q_R2.fastq > output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa


# blastn-short search
echo "blastn -db $4.fa  -query output/$1/$6/$2.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/$1/$6/$2.blastn"
#blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/2017-07-03_MiSeq/Blastn/60_Q_R1.fa -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn
#echo "blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn"
echo "blastn -db $4.fa  -query output/$1/$6/$3.fa -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8 > output/$1/$6/$3.blastn"
#blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/2017-07-03_MiSeq/Blastn/60_Q_R2.fa -task blastn-short -outfmt 6 -max_target_seqs 10 -evalue 1e-8 > output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn
#echo "blastn -db ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa  -query output/$1/Friedrich/$2.fastq -task blastn-short -outfmt 6 -max_target_seqs 20 -evalue 1e-8"


# python parse blastn output and make ppi map
echo "python main.py BLASTN $4.fa output/$1/$6/$2.blastn output/$1/$6/$3.blastn > output/$1/$6/$5.ppi.txt"
# cat output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn | awk '{print $1}' | uniq | wc    ## ==> 3209686
# cat output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn | awk '{print $1}' | uniq | wc    ## ==> 3029067
# python main.py BLASTN ../data/P170_4_library_MGj5615-Jun-2017122014.-100.fa output/2017-07-03_MiSeq/Blastn/60_Q_R1.blastn output/2017-07-03_MiSeq/Blastn/60_Q_R2.blastn    ##  ==> 2507307


#if __name__ == "__main__":
#    #CalculateDAMRatio( EXPERIMENT_TYPE, OUTPUT_FOLDER, UNCUT_FASTQ_PATH, DPNI_FASTQ_PATH, MBOI_FASTQ_PATH )
#    CalculateRecY2H( EXPERIMENT_TYPE, OUTPUT_FOLDER, NOSELCT_FASTQ_PATH, SELECT_FASTQ_PATH )
'''

if __name__ == "__main__":
    aCommand = command_center.Command()
    aCommand.AddCommand( "recY2H", CalculateRecY2H )
    
    aCommand.Run()

