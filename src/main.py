#!/usr/bin/python
import fasta
import sys
import time
import object
import command_center
import library_checker
#from exp import Experiment

'''
def ReadConfig( filepath ):
    aExp = Experiment()

    f = open( filepath )
    for line in f.xreadlines():        
        line = line.strip()
        if len(line) == 0: continue
        if "#######" in line: break
        if line[0] == "#": continue
        
        line = line.split( "#" )[0]
        [key,value] = line.split("=")
        key = key.strip().upper()
        value = value.strip()
        if key == "MOTIF_LEN": aExp.conditions[ "MOTIF_LEN" ] = int( value )
        elif key == "COUNT_CUT": aExp.conditions[ "COUNT_CUT" ] = int( value )
        elif key == "IS_REVERSE_COMPLEMENT": aExp.conditions[ "IS_REVERSE_COMPLEMENT" ] = eval( value ) 
        elif key == "SOURCE_FASTA": aExp.conditions[ "SOURCE_FASTA" ] = value # "Kd_Pair.fasta" # "40_40.fasta"
        elif key == "MAX_3PRIME_LENGTH": aExp.conditions[ "MAX_3PRIME_LENGTH" ] = int( value ) # MAX_3PRIME_LENGTH
        elif key == "READ1": aExp.conditions[ "READ1" ] = value # sys.argv[2] + "1.fastq.gz" #1A_12153_GCTCCA_read1.fastq.gz" 
        elif key == "READ2": aExp.conditions[ "READ2" ] = value # sys.argv[2] + "2.fastq.gz" #1A_12153_GCTCCA_read2.fastq.gz"  # 40%
        elif key == "READ1_FILTER": aExp.conditions[ "READ1_FILTER" ] = value #READ1_FILTER
        elif key == "READ2_FILTER": aExp.conditions[ "READ2_FILTER" ] = value #READ2_FILTER 
        elif key == "USE_ONLY_UNIQUE_SEQ": aExp.conditions[ "USE_ONLY_UNIQUE_SEQ" ] = eval( value ) # USE_ONLY_UNIQUE_SEQ
        elif key == "IGNORE_ORIENTATION": aExp.conditions[ "IGNORE_ORIENTATION" ] = eval( value ) #IGNORE_ORIENTATION
        elif key == "INFO_OUTPUT_PATH" : aExp.conditions[ "INFO_OUTPUT_PATH" ] = value 
        elif key == "PPI_OUTPUT_PATH" : aExp.conditions[ "PPI_OUTPUT_PATH" ] = value  
        elif key == "READ_LEN": aExp.conditions[ "READ_LEN" ] = int( value ) # 151
        elif key == "MOTIF_START": aExp.conditions[ "MOTIF_START" ] = int( value ) # 151
        elif key == "SAVE_INFO": aExp.conditions[ "SAVE_INFO" ] = eval( value ) # True / False # default = False
        elif key == "STOP_CNT": aExp.conditions[ "STOP_CNT" ] = int( value ) # 0
        else: aExp.conditions[ key ] = value

    f.close()

    return aExp
'''
def Run( argv ):
    aExp = ReadConfig( sys.argv[2] )
    # test
    #COUNT_CUT = 1
    #MOTIF_LEN = 10
    #aExp.conditions[ "COUNT_CUT" ] = COUNT_CUT
    #aExp.conditions[ "MOTIF_LEN" ] = MOTIF_LEN
    #aExp.conditions[ "PPI_OUTPUT_PATH" ] = aExp.conditions[ "PPI_OUTPUT_PATH" ] + ".%d.%d.test" % ( COUNT_CUT, MOTIF_LEN )
    aExp.run()

def ParameterScan( argv ):
    global pos_read1_dic
    global pos_read2_dic
    for COUNT_CUT in range(1,11):
        #for MOTIF_LEN in range(10,32,2):    
        for MOTIF_LEN in range(11,33,2):    
            pos_read1_dic = {}           
            pos_read2_dic = {}
            aExp = ReadConfig( sys.argv[2] )
            aExp.conditions[ "COUNT_CUT" ] = COUNT_CUT
            aExp.conditions[ "MOTIF_LEN" ] = MOTIF_LEN
            aExp.conditions[ "PPI_OUTPUT_PATH" ] = aExp.conditions[ "PPI_OUTPUT_PATH" ] + ".%d.%d" % ( COUNT_CUT, MOTIF_LEN )
            aExp.conditions[ "DEBUG" ] = False
            aExp.conditions[ "SAVE_INFO" ] = False
            aExp.run()

# input: ppi output
# output: real ppi matrix


def ExtractMatrix( postfix, filepath, title ):
    f = open( filepath )
    fout = open( filepath + postfix, "w" )
    bOK = False
    for line in f.xreadlines():
        if title in line:
            bOK = True
            continue
        if len( line.strip() ) == 0:
            bOK = False
            continue
        if bOK == True:
            print >> fout, line[:-1]
    f.close()
    fout.close()

def ExtractRealPPI( argv ):
    filepath = argv[2]
    ExtractMatrix( ".real.txt", filepath, "###################### print PRC_ vs PRC_ PPI matrix (Real PPI) #############################" )
    ExtractMatrix( ".prc_p.txt", filepath, "###################### print PRC vs P_ PPI matrix (Self PPI) #############################" )
    ExtractMatrix( ".p_prc.txt", filepath, "###################### print P_ vs PRC PPI matrix (Self PPI) #############################" )  
    ExtractMatrix( ".p_p.txt", filepath, "###################### print P_ vs P_ PPI matrix (Weired case) #############################" ) 
    ExtractMatrix( ".prc_v.txt", filepath, "###################### print PRC_ vs V_ PPI matrix #############################" ) 
    ExtractMatrix( ".prc_vrc.txt", filepath, "###################### print PRC_ vs VRC_ PPI matrix #############################" )
    ExtractMatrix( ".p_v.txt", filepath, "###################### print P_ vs V_ PPI matrix #############################" )       
    ExtractMatrix( ".p_vrc.txt", filepath, "###################### print P_ vs VRC_ PPI matrix #############################" )  
    ExtractMatrix( ".v_prc.txt", filepath, "###################### print V_ vs PRC_ PPI matrix #############################" )
    ExtractMatrix( ".vrc_prc.txt", filepath, "###################### print VRC_ vs PRC_ PPI matrix #############################" ) 
    ExtractMatrix( ".v_p.txt", filepath, "###################### print V_ vs P_ PPI matrix #############################" )
    ExtractMatrix( ".vrc_p.txt", filepath, "###################### print VRC_ vs P_ PPI matrix #############################" )
    ExtractMatrix( ".vrc_vrc.txt", filepath, "###################### print VRC_ vs VRC_ PPI matrix (Real PPI) #############################" )
    ExtractMatrix( ".vrc_v.txt", filepath, "###################### print VRC_ vs V_ PPI matrix (Real PPI) #############################" )
    ExtractMatrix( ".v_vrc.txt", filepath, "###################### print V_ vs VRC_ PPI matrix (Real PPI) #############################" )
    ExtractMatrix( ".v_v.txt", filepath, "###################### print V_ vs V_ PPI matrix (Real PPI) #############################" ) 

def PKAMotifAnalysis( argv ):
    from jslib import utr_designer
    
    fa = fasta.read_fasta( "../data/P170_4_library_MGj5615-Jun-2017122014.fa" )
    codon_dic = utr_designer.GetCodonDic()
    codon_dic[ "TAA" ] = "*"
    codon_dic[ "TGA" ] = "*"
    codon_dic[ "TAG" ] = "*"

    for id in fa:
        seq = fa[id].upper()
        aas = ""
        for i in range( 0, len(seq), 3 ):
            codon = seq[i:i+3]
            if len( codon ) < 3: break
            aas += codon_dic[ codon ]
        print "%s\t%s" % ( id, aas )
        

def RepSeq( profile, alphabet ):
    seq = ""
    i = 0
    for p in profile.array:
        print i, "\t", p, "\t", alphabet[ p.argmax() ]
        i += 1
        freq = float( p.max() ) / sum( p )
        if freq >= 0.8:
            seq += alphabet[ p.argmax() ]
        elif freq >= 0.4:
            seq += alphabet[ p.argmax() ].lower()
        else:
            seq += "n"
    return seq 



def ShowRefSeq( argv ):
    '''
        python main.py show_ref_seq ../2017-07-03_MiSeq/60_W_R1.fastq.gz
        python main.py show_ref_seq ../2017-07-03_MiSeq/60_W_R2.fastq.gz
    '''
    from weblogolib import corebio
    from fusion_ppi import ReadSequenceFile

    filepath1 = argv[2]
    read1 = ReadSequenceFile( filepath1 )

    seqs = []  

    index = 0 
    for line in read1.stdout.xreadlines():
        #print line # id
        seq = read1.stdout.next().strip()
        xx = read1.stdout.next()
        qual = read1.stdout.next()
        index += 1
        if index < 10000: continue

        seqs.append( seq )
        if index > 200000: break
    read1.stdout.close()
    read1.kill()  

    aSeqList = corebio.seq.SeqList(seqs, corebio.seq.unambiguous_dna_alphabet)   
    profile = aSeqList.profile()
    alphabet = profile.alphabet

    print RepSeq( profile, alphabet )

def ReadPPIMap( filepath, only_name = True ):
    output = {}
    f = open( filepath )
    titles = f.next()[:-1].split("\t")
    for line in f.xreadlines():
        fields = line[:-1].split("\t")
        bait = fields[0]
        if only_name == True: bait = bait.split("_")[0]
        for i in range( 1, len(fields) ):
            ppi_read = int( fields[i] )
            prey = titles[i]
            if only_name == True: prey = prey.split("_")[0]
            output[ (bait,prey) ] = ppi_read
    f.close()
    return output

def PerformanceCheck( argv ):
    if len( argv ) < 3:
        ppi_file = "./output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt"
    else:
        ppi_file = argv[2] # output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt
    print ppi_file
    ppi_map = ReadPPIMap( ppi_file )
    neg_map = {}
    pos_map = {} 
    '''
    Positive and Negative set from Sebastian Mauer's suggestions
    '''
    negC=[('dmrtb1','reep2'),('topbp1','dapk1'),('lam','larget'),('gaa','ccl14')]
    posC=[('dmrtb1','ctbp1'),('dmrtb1','rbfox2'),('cep55','cep55'),('p53','larget')]
    print negC
    print "bait \t prey \t ppi_read"
    for bait, prey in negC:
        print bait, "\t", prey, "\t", ppi_map[ (bait,prey) ]
        neg_map[(bait,prey)] = ppi_map[ (bait,prey) ]
        if bait != prey:
            print prey, "\t", bait, "\t", ppi_map[ (prey,bait) ]
            neg_map[(prey,bait)] = ppi_map[ (prey,bait) ]
    print posC
    print "bait \t prey \t ppi_read"
    for bait, prey in posC:
        print bait, "\t", prey, "\t", ppi_map[ (bait,prey) ]
        pos_map[(bait,prey)] = ppi_map[ (bait,prey) ]
        if bait != prey:
            print prey, "\t", bait, "\t", ppi_map[ (prey,bait) ]
            pos_map[(prey,bait)] = ppi_map[ (prey,bait) ]
    return ppi_map, neg_map, pos_map

def PerformanceP170Check( argv ):
    import numpy
    from scipy import stats
    '''
        python main.py performance_P170
    '''
    ## check spiked in positive/negative set
    ## For P170 x P170 experiments
    #python main.py performance ./output/2017-06-12_MiSeq/Blastn/S53.ppi.txt # BD P170-4 x AD P170-4; -W;
    #python main.py performance ./output/2017-06-12_MiSeq/Blastn/S54.ppi.txt # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2;
    #python main.py performance ./output/2017-07-03_MiSeq/Blastn/60_W.ppi.txt # BD P170-4 x AD P170-4; -W; 
    #python main.py performance ./output/2017-07-03_MiSeq/Blastn/60_Q.ppi.txt # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; 
    #python main.py performance ./output/2017-08-15_MiSeq/Blastn/S61_SWD.ppi.txt # BD P170-4 x AD P170-4; -W; (Seaprep)
    #python main.py performance ./output/2017-08-15_MiSeq/Blastn/S61_PWD.ppi.txt # BD P170-4 x AD P170-4; -W; 
    #python main.py performance ./output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; 

    no_selection_filelist = [ "./output/2017-06-12_MiSeq/Blastn/S53.ppi.txt", "./output/2017-07-03_MiSeq/Blastn/60_W.ppi.txt", "./output/2017-08-15_MiSeq/Blastn/S61_SWD.ppi.txt", "./output/2017-08-15_MiSeq/Blastn/S61_PWD.ppi.txt" ]
    selection_filelist = [ "./output/2017-06-12_MiSeq/Blastn/S54.ppi.txt", "./output/2017-07-03_MiSeq/Blastn/60_Q.ppi.txt", "./output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt" ]
  
    neg_all = []
    pos_all = []
    for filepath in no_selection_filelist:
        ppi_map, neg_map, pos_map = PerformanceCheck( [None,None,filepath] )
        total_reads = float( sum( ppi_map.values() ) )
        for bait, prey in neg_map:
            print "NEG:", bait, prey, neg_map[(bait,prey)]/total_reads*(len(ppi_map))
            if neg_map[(bait,prey)] != 0: neg_all.append( neg_map[(bait,prey)]/total_reads*(len(ppi_map)) )
        for bait, prey in pos_map:
            print "POS:", bait, prey, pos_map[(bait,prey)]/total_reads*(len(ppi_map))
            if pos_map[(bait,prey)] != 0: pos_all.append( pos_map[(bait,prey)]/total_reads*(len(ppi_map)) )
    print "No Selection:", numpy.mean( neg_all ), len(neg_all), numpy.mean( pos_all ), len( pos_all ), stats.mannwhitneyu( neg_all, pos_all )
       
    for filepath in selection_filelist:
        ppi_map, neg_map, pos_map = PerformanceCheck( [None,None,filepath] )
        total_reads = float( sum( ppi_map.values() ) )
        for bait, prey in neg_map:
            print "NEG:", bait, prey, neg_map[(bait,prey)]/total_reads*(len(ppi_map))
            if neg_map[(bait,prey)] != 0: neg_all.append( neg_map[(bait,prey)]/total_reads*(len(ppi_map)) )
        for bait, prey in pos_map:
            print "POS:", bait, prey, pos_map[(bait,prey)]/total_reads*(len(ppi_map))
            if pos_map[(bait,prey)] != 0: pos_all.append( pos_map[(bait,prey)]/total_reads*(len(ppi_map)) )
 
    print "Selection:", numpy.mean( neg_all ), len(neg_all), numpy.mean( pos_all ), len( pos_all ), stats.mannwhitneyu( neg_all, pos_all )

    '''
    Output:
        No Selection: 0.320275342749 18 0.782570789647 13 (85.5, 0.10718855347264772)
        Selection: 0.316370718779 26 22.4757085384 28 (176.5, 0.00060264500742833926)
    '''

    
if __name__ == "__main__":
    aCommand = command_center.Command()
    aCommand.AddCommand( "generate_last_nts", library_checker.GenerateLastNts )
    aCommand.AddCommand( "compare_last_nts", library_checker.CompareLastNts )
    aCommand.AddCommand( "check_last_nts", library_checker.CheckLastNts )
    aCommand.AddCommand( "fastq_to_fasta", library_checker.FastqToFasta )
    aCommand.AddCommand( "trim_fastq_to_fasta", library_checker.TrimFastqToFasta )
    aCommand.AddCommand( "trim_fastq_to_fastq", library_checker.TrimFastqToFastq )
    aCommand.AddCommand( "sam_to_matrix", library_checker.SamToMatrix )
    aCommand.AddCommand( "performance", PerformanceCheck )
    aCommand.AddCommand( "performance_P170", PerformanceP170Check )
    
    aCommand.AddCommand( "run", Run )

    aCommand.AddCommand( "parameter_scan", ParameterScan )
    aCommand.AddCommand( "extract_real_ppi", ExtractRealPPI )

    aCommand.AddCommand( "PKA_motif_analysis", PKAMotifAnalysis )
    
    aCommand.AddCommand( "show_ref_seq", ShowRefSeq )
    aCommand.AddCommand( "BLASTN", library_checker.BLASTN )
    aCommand.AddCommand( "BLASTN_NEW", library_checker.BLASTN_NEW )
    aCommand.AddCommand( "BLASTN_BARCODE", library_checker.BLASTN_BARCODE )
    aCommand.AddCommand( "BLASTN_RELAXED", library_checker.BLASTN_RELAXED )
    aCommand.Run()


    
    
    
    
    
    
    



