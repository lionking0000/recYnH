######################################################################################################################################  
#
# Version 1.0
# visualization.R
#
#
######################################################################################################################################


source( "/usr/local/bin/init.R" )

#######################
# Read Argument
#######################
print("[recYnH] Visualization : Read Argument")
argv=commandArgs(T)
print( argv[1] )	# args.program # the experiments type ('Y2H'|'Y3H')
print( argv[2] )	# args.matrix1 # the interaction matrix of selection condition
print( argv[3] )	# args.matrix2 # the interaction matrix of non-selection condition
print( argv[4] )	# args.output # the output folder name
print( argv[5] )	# args.name # output name

IS_output_path <- sprintf( "%s/%s.is.txt", argv[4], argv[5] )
NIS_output_path <- sprintf( "%s/%s.nis.txt", argv[4], argv[5] )

PNG_is_output_path <- sprintf( "%s/%s.is.png", argv[4], argv[5] )
PNG_nis_output_path <- sprintf( "%s/%s.nis.png", argv[4], argv[5] )

#m1 = drawMatrix("../share/output/test/recYnH.ppi.txt",0.5) 
m1 = drawMatrix( argv[2], 0.5 ) 
#m2 = drawMatrix("../share/output/test/recYnH.ppi.txt",0.5) 
m2 = drawMatrix( argv[3], 0.5) 
is2 = InteractionScores(m1,m2,1.0)  
nis2 = NewInteractionScores(m1,m2,1.0)  
CorMatrix(is2,nis2)
#SaveMatrixPlot(is2, 1.0, "../share/output/test/recYnH.ppi.is.png", 4 ) 
SaveMatrixPlot(is2, 1.0, PNG_is_output_path, 4 ) 
#SaveMatrixPlot(nis2, 1.0, "../share/output/test/recYnH.ppi.nis.png", 4 ) 
SaveMatrixPlot(nis2, 1.0, PNG_nis_output_path, 4 ) 

write.table(is2,file=IS_output_path,sep="\t")
write.table(nis2,file=NIS_output_path,sep="\t")

quit()
###################################################################################################################################### 
## Initilize functions and load data
source( "/Volumes/users/lserrano/jyang/work/Mireia/src/init.R" )


######################################################################################################################################


       

## NOT USED ##   
# 2016-05-16_MiSeq; 33x54 RNA binding proteins
m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-05-16_MiSeq/14075.ppi.txt.refine.txt",0.5) 
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-05-16_MiSeq/14076.ppi.txt.refine.txt",0.5) 
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-05-16_MiSeq/14077.ppi.txt.refine.txt",0.5) 
  
# 2016-06-27_MiSeq; 33x54 RNA binding proteins   
m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-06-27_MiSeq/14477.ppi.txt.refine.txt",0.5) 
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-06-27_MiSeq/14478.ppi.txt.refine.txt",0.5) 
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-06-27_MiSeq/14479.ppi.txt.refine.txt",0.5)
                                                 
# 2016-07-06_MiSeq; 33x54 RNA binding proteins                                               
m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-06_MiSeq/14813.ppi.txt.refine.txt",0.5)
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-06_MiSeq/14814.ppi.txt.refine.txt",0.5)  
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-06_MiSeq/14815.ppi.txt.refine.txt",0.5)
m4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-06_MiSeq/14816.ppi.txt.refine.txt",0.5)
m5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-06_MiSeq/14817.ppi.txt.refine.txt",0.5)
m6 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-06_MiSeq/14818.ppi.txt.refine.txt",0.5)
                                                 
# 2016-07-08_MiSeq; 33x54 RNA binding proteins   
m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-08_MiSeq/14813.ppi.txt.refine.txt",0.5)
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-08_MiSeq/14814.ppi.txt.refine.txt",0.5)  
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-08_MiSeq/14815.ppi.txt.refine.txt",0.5)
m4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-08_MiSeq/14816.ppi.txt.refine.txt",0.5)
m5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-08_MiSeq/14818.ppi.txt.refine.txt",0.5)
m6 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2016-07-08_MiSeq/14819.ppi.txt.refine.txt",0.5)
is2 = InteractionScores(m1,m2,1.0)  
is3 = InteractionScores(m1,m3,1.0)  
CorMatrix(is2,is3)
nis2 = NewInteractionScores(m1,m2,1.0)  
nis3 = NewInteractionScores(m1,m3,1.0)  
CorMatrix(nis2,nis3)
  
#======================================================================================================================================================================================
#                  
# 2016-12-22_MiSeq; Roth75; MGj46  
#
# These are never cloned into pENTR223: 
# R75_41,WDR7,WD repeat domain 7
# R75_72,PLEKHG7,pleckstrin homology domain containing, family G (with RhoGef domain) member 7
#
m1_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt.3.15.refine.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1947500"
m2_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt.3.15.refine.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1533027" 
m3_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17545.ppi.txt.3.15.refine.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2121357" 
m4_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17546.ppi.txt.3.15.refine.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1819108"    

# from blastn-short
#bm1_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
#bm2_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
#bm3_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
#bm4_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"    
                                               
PercentageIncreaseByBlast = ( c( sum(bm1_1)/sum(m1_1), sum(bm2_1)/sum(m2_1), sum(bm3_1)/sum(m3_1), sum(bm4_1)/sum(m4_1) ) - 1 ) * 100  # 2.03% increase

is2_1 = InteractionScores(m1_1,m2_1,1.0)  
is4_1 = InteractionScores(m3_1,m4_1,1.0) 
                                             
nis2_1 = NewInteractionScores(m1_1,m2_1,1.0)  
nis4_1 = NewInteractionScores(m3_1,m4_1,1.0) 


bis2_1 = InteractionScores(bm1_W,bm1_A,1.0,7)  
bis4_1 = InteractionScores(bm1_RW,bm1_RA,1.0,7) 

bnis2_1 = NewInteractionScores(bm1_W,bm1_A,1.0,7)  
bnis4_1 = NewInteractionScores(bm1_RW,bm1_RA,1.0,7) 

# comparison with BioGrid
x_BioGrid = MaskMatrix(bis2_1,BioGrid)
simple_plotMatrix(x_BioGrid>3,1,7)
y_BioGrid = MaskMatrix(bis2_1,BioGrid, 1)
simple_plotMatrix(y_BioGrid>3,1,7) 

x_HIPPIE = MaskMatrix(bis2_1,HIPPIE)
simple_plotMatrix(x_HIPPIE>3,1,7)
y_HIPPIE = MaskMatrix(bis2_1,HIPPIE, 1)
simple_plotMatrix(y_HIPPIE>3,1,7) 




CorMatrix(is2_1,is4_1,method="pearson")    # 0.9120527 --> 0.9227895 (proper one?)  
CorMatrix(nis2_1,nis4_1)    # 0.9747  --> 0.9628 (proper one)    
CorMatrix(bis2_1,bis4_1)  # 0.9217131     
CorMatrix(bnis2_1,bnis4_1)  # 0.9610908      

     
#sort(rowSums(m1_1))    
# PLEKHG7     WDR7    UBE2I    AP2B1      BIK    PSMC1   DMRTB1    RING1   TSG101 
#       0        1     2859     3253     3338     4653     5782     6531     7730
#sort(colSums(m1_1))
#    WDR7  PLEKHG7      p53    CTBP1   BCL2L2    VAMP3    IKZF5    VAMP2     PCNA 
#       0        1       53      147      156      157      263      267      373
rm1_1 = RefineMatrix(m1_1) # "% detected cnt = 86.4" : 580400  ?? --> 1876156         
BasicStats(rm1_1)
rm2_1 = RefineMatrix(m2_1) 
rm3_1 = RefineMatrix(m3_1) # "% detected cnt = 85.0" : 2057781         
BasicStats(rm3_1)
rm4_1 = RefineMatrix(m4_1)
SaveMatrixPlot(rm1_1, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2016-12-22_MiSeq_R75_MGj46/R71_76.S1.png" )
SaveMatrixPlot(rm2_1, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2016-12-22_MiSeq_R75_MGj46/R71_76.S2.png" )
SaveMatrixPlot(rm3_1, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2016-12-22_MiSeq_R75_MGj46/R71_76.S3.png" )
SaveMatrixPlot(rm4_1, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2016-12-22_MiSeq_R75_MGj46/R71_76.S4.png" )
ris2_1 = InteractionScores(rm1_1,rm2_1,1.0)    
SaveMatrixPlot(2^ris2_1, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2016-12-22_MiSeq_R75_MGj46/R71_76.IS2.png" )                         
ris4_1 = InteractionScores(rm3_1,rm4_1,1.0)    
SaveMatrixPlot(2^ris4_1, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2016-12-22_MiSeq_R75_MGj46/R71_76.IS4.png" )


# 2017-02-04_MiSeq; Roth75; MGj48   ==> seems not work
m1_x = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2017-02-04_MiSeq/17861.ppi.txt.refine.txt",0.5) # no selection  --> It seems there is selection ; WCC R10D10 
m2_x = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2017-02-04_MiSeq/17862.ppi.txt.refine.txt",0.5) # selection                                     ; ACC (1/2 Aba) R10D10
m3_x = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2017-02-04_MiSeq/17863.ppi.txt.refine.txt",0.5) # selection                                     ; A1CC (1 Aba) R10D10
m4_x = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2017-02-04_MiSeq/17864.ppi.txt.refine.txt",0.5) # no selection                                  ; WCC R6D6
# no DNA: m5_x = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2017-02-04_MiSeq/17865.ppi.txt.refine.txt",0.5) # selection                           ; ACC (1/2 Aba) R6D6
# no DNA: m6_x = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src_backup_20170622/output/2017-02-04_MiSeq/17866.ppi.txt.refine.txt",0.5) # selection                           ; A1CC (1 Aba) R6D6
is2_x = InteractionScores(m1_x,m2_x,1.0)  
is3_x = InteractionScores(m1_x,m3_x,1.0)  
CorMatrix(is2_x,is3_x)     # 0.4334 
nis2_x = NewInteractionScores(m1_x,m2_x,1.0)  
nis3_x = NewInteractionScores(m1_x,m3_x,1.0)  
CorMatrix(nis2_x,nis3_x)   # 0.7182
rm1_x = RefineMatrix(m1_x) 
rm2_x = RefineMatrix(m2_x) 
rm3_x = RefineMatrix(m3_x)
rm4_x = RefineMatrix(m4_x) 
SaveMatrixPlot(rm1_x, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-04_MiSeq_R75_MGj48/R71_76.S1.png" )
SaveMatrixPlot(rm2_x, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-04_MiSeq_R75_MGj48/R71_76.S2.png" )
SaveMatrixPlot(rm3_x, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-04_MiSeq_R75_MGj48/R71_76.S3.png" )
SaveMatrixPlot(rm4_x, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-04_MiSeq_R75_MGj48/R71_76.S4.png" )
ris2_x = InteractionScores(rm1_x,rm2_x,1.0)    
SaveMatrixPlot(2^ris2_x, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-04_MiSeq_R75_MGj48/R71_76.IS2.png" ) 
ris3_x = InteractionScores(rm1_x,rm3_x,1.0)    
SaveMatrixPlot(2^ris3_x, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-04_MiSeq_R75_MGj48/R71_76.IS3.png" ) 


#
# 2017-02-22_MiSeq; Roth75; R75_MGj51
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
# R75_11,AP2B1,adaptor-related protein complex 2, beta 1 subunit
# R75_30,CTBP1,C-terminal binding protein 1
# R75_52,RBBP8,retinoblastoma binding protein 8
# R75_58,FHL3,four and a half LIM domains 3 
# R75_70,TEX11,testis expressed 11
#
# These are never cloned into pENTR223: 
# R75_41,WDR7,WD repeat domain 7
# R75_72,PLEKHG7,pleckstrin homology domain containing, family G (with RhoGef domain) member 7
#
m1_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/S1.ppi.txt.3.15.refine.txt",0.5) # no selection; R75_MGj51  "sum = 919439"
m2_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/S2.ppi.txt.3.15.refine.txt",0.5) # selection; ACCRD         "sum = 1956718"
m3_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/S3.ppi.txt.3.15.refine.txt",0.5) # selection; QCCRD         "sum = 791387"   

bm1_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; R75_MGj51  "sum = 960245"
bm2_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; ACCRD         "sum = 2024754"
bm3_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK


is2_2 = InteractionScores(m1_2,m2_2,1.0)  
is3_2 = InteractionScores(m1_2,m3_2,1.0)   
CorMatrix(is2_2,is3_2)		# 0.6129

nis2_2 = NewInteractionScores(m1_2,m2_2,1.0)  
nis3_2 = NewInteractionScores(m1_2,m3_2,1.0)   
CorMatrix(nis2_2,nis3_2)	# 0.8309  -->  0.8425  

#nis2_2 = NewInteractionScores(m1_2,m2_2,1.0,9,2)  
#nis3_2 = NewInteractionScores(m1_2,m3_2,1.0,9,2)   
#CorMatrix(nis2_2,nis3_2)	# 0.8474
 
                                                     
#nis2_2 = NewInteractionScores(m1_2,m2_2,1.0,9,2)  
#nis3_2 = NewInteractionScores(m1_2,m3_2,1.0,9,2)


bis2_A = InteractionScores(bm2_W,bm2_A,1.0)  
bis2_Q = InteractionScores(bm2_W,bm2_Q,1.0)   
CorMatrix(bis2_A,bis2_Q)		# 0.6835282

nbis2_A = NewInteractionScores(bm2_W,bm2_A,1.0)  
nbis2_Q = NewInteractionScores(bm2_W,bm2_Q,1.0)   
CorMatrix(nbis2_A,nbis2_Q)	# 0.8438049

# comparison with BioGrid
x_KnownPPI = MaskMatrix(bis2_A,KnownPPI,FALSE)
simple_plotMatrix(x_KnownPPI,1,7) 
simple_plotMatrix(x_KnownPPI>3,1,7)
y_KnownPPI = MaskMatrix(bis2_A,KnownPPI,TRUE)   
simple_plotMatrix(y_KnownPPI,1,7) 
simple_plotMatrix(y_KnownPPI>3,1,7) 

x_KnownPPI = MaskMatrix(nbis2_A,KnownPPI,FALSE)
simple_plotMatrix(x_KnownPPI,1,7) 
simple_plotMatrix(x_KnownPPI>3,1,7)
y_KnownPPI = MaskMatrix(nbis2_A,KnownPPI,TRUE)   
simple_plotMatrix(y_KnownPPI,1,7) 
simple_plotMatrix(y_KnownPPI>3,1,7)



#sort(rowSums(m1)) # BD-vector
#   CTBP1  PLEKHG7    AP2B1     WDR7    RBBP8    TEX11     FHL3    UBE2I   TOPBP1 
#       0        0        1        1        9       62       71     1543     1738
#sort(colSums(m1)) # AD-vector
# PLEKHG7     WDR7   BCL2L2    VAMP2    CTBP1    IKZF5      MME      p53    VAMP3 
#       0        0       51       53       57       76      123      161      202
rm1_2 = RefineMatrix(m1_2) # "% detected cnt = 87.2"  
BasicStats(rm1_2)
rm2_2 = RefineMatrix(m2_2) 
rm3_2 = RefineMatrix(m3_2) 
SaveMatrixPlot(rm1_2, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-22_MiSeq_R75_MGj51/R71_76.S1.png" )
SaveMatrixPlot(rm2_2, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-22_MiSeq_R75_MGj51/R71_76.S2.png" )  
SaveMatrixPlot(rm3_2, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-22_MiSeq_R75_MGj51/R71_76.S3.png" )
ris2_2 = InteractionScores(rm1_2,rm2_2,1.0)    
SaveMatrixPlot(2^ris2_2, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-22_MiSeq_R75_MGj51/R71_76.IS2.png" )  
ris3_2 = InteractionScores(rm1_2,rm3_2,1.0)    
SaveMatrixPlot(2^ris3_2, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-02-22_MiSeq_R75_MGj51/R71_76.IS3.png" )

                                                                                           
#
# 2017-03-03_MiSeq; Roth75
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
# R75_11,AP2B1,adaptor-related protein complex 2, beta 1 subunit
# R75_30,CTBP1,C-terminal binding protein 1
# R75_52,RBBP8,retinoblastoma binding protein 8
# R75_58,FHL3,four and a half LIM domains 3
# R75_70,TEX11,testis expressed 11
#     
# 3.15                                                                                                                             
m1_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S1.ppi.txt.3.15.refine.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = 1369865"
m2_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S2.ppi.txt.3.15.refine.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = 947672"
m3_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S3.ppi.txt.3.15.refine.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1218835"
m4_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S4.ppi.txt.3.15.refine.txt",0.5) # no selection; 53_A73e_WRD	18712	0,30%	"sum = 72" # For toxicity
m5_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S5.ppi.txt.3.15.refine.txt",0.5) # no selection; 53_eB68_WRD	18713	0,30%	"sum = 619"    # For auto-activation
m6_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S6.ppi.txt.3.15.refine.txt",0.5) # selection; 53_eB68_ARD	18714	0,30%	"sum = 1192"   # For auto-activation
m7_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S7.ppi.txt.3.15.refine.txt",0.5) # selection; 53_QL4RD	18715	24,7%	"sum = 715800"    

bm1_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = 1414435"
bm2_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = 969159"
bm3_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"
bm4_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S4.ppi.txt",0.5) # no selection; 53_A73e_WRD	18712	0,30%	"sum = 23" # For toxicity     
bm5_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S5.ppi.txt",0.5) # no selection; 53_eB68_WRD	18713	0,30%	"sum = 169"    # For auto-activation 
bm6_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S6.ppi.txt",0.5) # selection; 53_eB68_ARD	18714	0,30%	"sum = 1296"   # For auto-activation  
bm7_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt",0.5) # selection; 53_QL4RD	18715	24,7%	"sum = 734766"    

# more output/2017-03-03_MiSeq/Blastn/S7_R1.blastn | grep BCL2 | sort -k4 | awk '{print $4}' | uniq -c 

is2_3 = InteractionScores(m1_3,m2_3,1.0)
is3_3 = InteractionScores(m1_3,m3_3,1.0)
is7_3 = InteractionScores(m1_3,m7_3,1.0)  
CorMatrix(is2_3,is3_3)	# 0.7122323
CorMatrix(is2_3,is7_3)	# 0.5847621

nis2_3 = NewInteractionScores(m1_3,m2_3,1.0)
nis3_3 = NewInteractionScores(m1_3,m3_3,1.0)
nis7_3 = NewInteractionScores(m1_3,m7_3,1.0) 
CorMatrix(nis2_3,nis3_3)	# 0.7752 	--> 0.782
CorMatrix(nis2_3,nis7_3)	# 0.5841 	--> 0.5893
simple_plotLog2Matrix(nis2_3*nis3_3*nis7_3,1)

bis2_3 = InteractionScores(bm1_3,bm2_3,1.0)
bis3_3 = InteractionScores(bm1_3,bm3_3,1.0)
bis7_3 = InteractionScores(bm1_3,bm7_3,1.0) 
CorMatrix(bis2_3,bis3_3)	# 0.722062
CorMatrix(bis2_3,bis7_3)	# 0.5871871

nbis2_3 = NewInteractionScores(bm1_3,bm2_3,1.0)
nbis3_3 = NewInteractionScores(bm1_3,bm3_3,1.0)
nbis7_3 = NewInteractionScores(bm1_3,bm7_3,1.0) 
CorMatrix(nbis2_3,nbis3_3)	# 0.7849249
CorMatrix(nbis2_3,nbis7_3)	# 0.5897572


rm1_3 = RefineMatrix(m1_3) # "% detected cnt = 89.5"      
BasicStats(rm1_3)
rm2_3 = RefineMatrix(m2_3) 
rm3_3 = RefineMatrix(m3_3)
rm4_3 = RefineMatrix(m4_3)
rm5_3 = RefineMatrix(m5_3)
rm6_3 = RefineMatrix(m6_3)
rm7_3 = RefineMatrix(m7_3)
SaveMatrixPlot(rm1_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S1.png" )
SaveMatrixPlot(rm2_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S2.png" )  
SaveMatrixPlot(rm3_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S3.png" )  
SaveMatrixPlot(rm4_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S4.png" )  
SaveMatrixPlot(rm5_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S5.png" )  
SaveMatrixPlot(rm6_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S6.png" )  
SaveMatrixPlot(rm7_3, 1, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S7.png" )  
ris2_3 = InteractionScores(rm1_3,rm2_3,1.0)    
SaveMatrixPlot(2^ris2_3, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.IS2.png" )  
ris3_3 = InteractionScores(rm1_3,rm3_3,1.0)    
SaveMatrixPlot(2^ris3_3, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.IS3.png" )
ris7_3 = InteractionScores(rm1_3,rm7_3,1.0)    
SaveMatrixPlot(2^ris7_3, 0.0, "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.IS7.png" )
melted_rm1_3 = melt(rm1_3)
write.csv(melted_rm1_3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_R71_76.S1.txt")
melted_ris2_3 = melt(ris2_3)
write.csv(melted_ris2_3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_R71_76.IS2.txt")
melted_ris3_3 = melt(ris3_3)
write.csv(melted_ris3_3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_R71_76.IS3.txt")                
melted_ris7_3 = melt(ris7_3)
write.csv(melted_ris7_3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_R71_76.IS7.txt")
     
# Check Bias
barplot( sort(rowSums(rm1_3)) ) # with DNA binding domain
barplot( sort(colSums(rm1_3)) ) # with Activation domain              
                        
col_sum_1_1 = colSums(rm1_1)
col_sum_1_3 = colSums(rm1_3)
cor(R75_gene_length,col_sum_1_1[match( R75_gene_name, names(col_sum_1_1))],use = "na")
cor(R75_gene_length,col_sum_1_3[match( R75_gene_name, names(col_sum_1_3))],use = "na") 
plot(R75_gene_length,col_sum_1_1[match( R75_gene_name, names(col_sum_1_1))])
plot(R75_gene_length,col_sum_1_3[match( R75_gene_name, names(col_sum_1_3))])

row_sum_1_1 = rowSums(rm1_1)
row_sum_1_3 = rowSums(rm1_3)
cor(R75_gene_length,row_sum_1_1[match( R75_gene_name, names(row_sum_1_1))],use = "na")  
cor(R75_gene_length,row_sum_1_3[match( R75_gene_name, names(row_sum_1_3))],use = "na")
plot(R75_gene_length,row_sum_1_1[match( R75_gene_name, names(row_sum_1_1))])  
plot(R75_gene_length,row_sum_1_3[match( R75_gene_name, names(row_sum_1_3))])
                                                                                                                                                             











#================================================================================================================================================================
# Read BFG_Y2H (Only R75 Set we used 71 x 76 )
ccc1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_1_+His.tsv.refined.txt",0.5)
ccc1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_1_-His.tsv.refined.txt",0.5)
ccc1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_1_+3AT.tsv.refined.txt",0.5)            
ccc1_Row_Index = c( 3,4,6,8,9,11,13,14,15,19,24,29,30,31,32,35,36,45,48,49,50,51,56,57,59,60,61,64,65,66,67,70,71,72,74 );
ccc1_Col_Index = c( 1,2,5,7,12,13,14,16,17,18,25,26,29,33,35,37,38,39,41,42,46,47,49,51,53,55,57,58,62,65,67,68,71,73 );
simple_plotLog2Matrix(ccc1_m1[ccc1_Row_Index,ccc1_Col_Index], 1) # read sum = 7042
simple_plotLog2Matrix(ccc1_m2[ccc1_Row_Index,ccc1_Col_Index], 1) # read sum = 4574
simple_plotLog2Matrix(ccc1_m3[ccc1_Row_Index,ccc1_Col_Index], 1) # read sum = 4433                             
CCC_1_is2 = InteractionScores(ccc1_m1,ccc1_m2,1.0,7)  
CCC_1_is3 = InteractionScores(ccc1_m1,ccc1_m3,1.0,7)  
CorMatrix(CCC_1_is2,CCC_1_is3) 		# 0.717    
nCCC_1_is2 = NewInteractionScores(ccc1_m1,ccc1_m2,1.0)  
nCCC_1_is3 = NewInteractionScores(ccc1_m1,ccc1_m3,1.0)  
CorMatrix(nCCC_1_is2,nCCC_1_is3)	# 0.7138


ccc2_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_2_+His.tsv.refined.txt",0.5)
ccc2_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_2_-His.tsv.refined.txt",0.5)
ccc2_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_2_+3AT.tsv.refined.txt",0.5)                                                                                                                             
ccc2_Row_Index = c( 3,4,6,8,9,11,13,14,15,19,24,29,30,31,32,35,36,45,48,49,50,51,56,57,59,60,61,64,65,66,67,70,71,72,74 );
ccc2_Col_Index = c( 1,2,5,7,12,13,14,16,17,18,25,26,29,33,35,37,38,39,41,42,46,47,49,51,53,55,57,58,62,65,67,68,71,73 );
simple_plotLog2Matrix(ccc2_m1[ccc2_Row_Index,ccc2_Col_Index], 1) # read sum = 10105 
simple_plotLog2Matrix(ccc2_m2[ccc2_Row_Index,ccc2_Col_Index], 1) # read sum = 3611 
simple_plotLog2Matrix(ccc2_m3[ccc2_Row_Index,ccc2_Col_Index], 1) # read sum = 5088                               
CCC_2_is2 = InteractionScores(ccc2_m1,ccc2_m2,1.0)  
CCC_2_is3 = InteractionScores(ccc2_m1,ccc2_m3,1.0)  
CorMatrix(CCC_2_is2,CCC_2_is3)		# 0.774 
nCCC_2_is2 = NewInteractionScores(ccc2_m1,ccc2_m2,1.0)  
nCCC_2_is3 = NewInteractionScores(ccc2_m1,ccc2_m3,1.0)  
CorMatrix(nCCC_2_is2,nCCC_2_is3)	# 0.7699
                                                                   
                                                               
cent_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_-AA_+His.tsv.refined.txt",0.5)  # read sum = 515092
cent_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_-AA_-His.tsv.refined.txt",0.5)  # read sum = 3743381
cent_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_-AA_+3AT.tsv.refined.txt",0.5)  # read sum = 11983468     
cent_Row_Index = c( 3,4,6,8,9,11,13,14,15,19,24,30,31,32,35,36,41,45,48,49,50,56,57,59,60,61,65,70,74 );
cent_Col_Index = c( 2,5,7,12,13,16,17,18,25,26,29,33,35,37,38,39,41,46,47,51,53,55,57,58,62,65,71,73 );
simple_plotLog2Matrix(cent_m1[cent_Row_Index,cent_Col_Index], 1) # read sum = 515092
simple_plotLog2Matrix(cent_m2[cent_Row_Index,cent_Col_Index], 1) # read sum = 3743381
simple_plotLog2Matrix(cent_m3[cent_Row_Index,cent_Col_Index], 1) # read sum = 11983468
CENT_AA_is2 = InteractionScores(cent_m1,cent_m2,1.0)  
CENT_AA_is3 = InteractionScores(cent_m1,cent_m3,1.0)  
CorMatrix(CENT_AA_is2,CENT_AA_is3)
                                                       

centall_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_ALL_+His.tsv.refined.txt",0.5)  # read sum = 531690  
centall_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_ALL_-His.tsv.refined.txt",0.5)  # read sum = 566750  
centall_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_ALL_+3AT.tsv.refined.txt",0.5)  # read sum = 5736220         
centall_Row_Index = c( 3,4,6,8,9,11,13,14,15,19,24,30,31,32,35,36,41,45,48,49,50,56,57,59,60,61,65,70,74 );
centall_Col_Index = c( 2,5,7,12,13,16,17,18,25,26,29,33,35,37,38,39,41,46,47,51,53,55,57,58,62,65,71,73 );
Null_centall = centall_m1 * 0 
centall_Row_Not_Index = (1:78)[-centall_Row_Index]
centall_Col_Not_Index = (1:78)[-centall_Col_Index]   
Null_centall[centall_Row_Not_Index,] = 1
Null_centall[,centall_Col_Not_Index] = 1  
simple_plotMatrix(Null_centall, 1)     
BasicStats(Null_centall)

simple_plotLog2Matrix(centall_m1[centall_Row_Index,centall_Col_Index], 1)
simple_plotLog2Matrix(centall_m2[centall_Row_Index,centall_Col_Index], 1)
simple_plotLog2Matrix(centall_m3[centall_Row_Index,centall_Col_Index], 1)
CENT_ALL_is2 = InteractionScores(centall_m1,centall_m2,1.0)  
CENT_ALL_is3 = InteractionScores(centall_m1,centall_m3,1.0)  
CorMatrix(CENT_ALL_is2,CENT_ALL_is3)
nCENT_ALL_is2 = NewInteractionScores(centall_m1,centall_m2,1.0)  
nCENT_ALL_is3 = NewInteractionScores(centall_m1,centall_m3,1.0)  
CorMatrix(nCENT_ALL_is2,nCENT_ALL_is3)
                                             

cva1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+His.tsv.refined.txt",0.5,7)   	 # read sum = 8597 
cva1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_-His.tsv.refined.txt",0.5,7)     # read sum = 5588 
cva1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+3AT.tsv.refined.txt",0.5,7)     # read sum = 6009 
cva1_Row_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,26,28,29,30,31,32,33,35,36,37,39,41,42,45,46,47,48,50,51,52,53,54,56,57,58,59,60,61,62,64,66,67,70,71,72,73,74 );
cva1_Col_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,24,26,28,29,31,32,33,35,36,37,38,41,42,45,46,47,48,50,51,52,53,56,57,58,60,61,62,63,64,66,67,68,70,71,72,74 );    
Null_cva1 = cva1_m1 * 0 #matrix(0,nrow=78,ncol=78)
cva1_Row_Not_Index = (1:78)[-cva1_Row_Index]
cva1_Col_Not_Index = (1:78)[-cva1_Col_Index]   
Null_cva1[cva1_Row_Not_Index,] = 1     
Null_cva1[,cva1_Col_Not_Index] = 1  
simple_plotMatrix(Null_cva1, 1, 7)   
BasicStats(Null_cva1)      
simple_plotLog2Matrix(cva1_m1[cva1_Row_Index,cva1_Col_Index], 1)
simple_plotLog2Matrix(cva1_m2[cva1_Row_Index,cva1_Col_Index], 1)
simple_plotLog2Matrix(cva1_m3[cva1_Row_Index,cva1_Col_Index], 1)
CVA_1_is2 = InteractionScores(cva1_m1,cva1_m2,1.0,filename="cva_1_is2.pdf",cellwidth = 10,cellheight = 10)     
CVA_1_is3 = InteractionScores(cva1_m1,cva1_m3,1.0,filename="cva_1_is3.pdf",cellwidth = 10,cellheight = 10)     
CorMatrix(CVA_1_is2,CVA_1_is3)		# 0.6172 
nCVA_1_is2 = NewInteractionScores(m1,m2,1.0)
nCVA_1_is3 = NewInteractionScores(m1,m3,1.0)
CorMatrix(nCVA_1_is2,nCVA_1_is3)    # 0.7195
ocva1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+His.tsv.original.txt",0.5)
ocva1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_-His.tsv.original.txt",0.5)
ocva1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+3AT.tsv.original.txt",0.5)
OCVA_1_is2 = InteractionScores(ocva1_m1,ocva1_m2,1.0)  
OCVA_1_is3 = InteractionScores(ocva1_m1,ocva1_m3,1.0)  
CorMatrix(OCVA_1_is2,OCVA_1_is3)		#                          
nOCVA_1_is2 = NewInteractionScores(ocva1_m1,ocva1_m2,1.0)  
nOCVA_1_is3 = NewInteractionScores(ocva1_m1,ocva1_m3,1.0)  
CorMatrix(nOCVA_1_is2,nOCVA_1_is3)	#



cva2_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_+His.tsv.refined.txt",0.5) # read sum = 10703
cva2_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_-His.tsv.refined.txt",0.5) # read sum = 5858
cva2_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_+3AT.tsv.refined.txt",0.5) # read sum = 5769
cva2_Row_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,26,28,29,30,31,32,33,35,36,37,39,41,42,45,46,47,48,50,51,52,53,54,56,57,58,59,60,61,62,64,66,67,70,71,72,73,74 );
cva2_Col_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,24,26,28,29,31,32,33,35,36,37,38,41,42,45,46,47,48,50,51,52,53,56,57,58,60,61,62,63,64,66,67,68,70,71,72,74 );
Null_cva2 = cva2_m1 * 0 #matrix(0,nrow=78,ncol=78)
cva2_Row_Not_Index = (1:78)[-cva2_Row_Index]
cva2_Col_Not_Index = (1:78)[-cva2_Col_Index]   
Null_cva2[cva2_Row_Not_Index,] = 1     
Null_cva2[,cva2_Col_Not_Index] = 1  
simple_plotMatrix(Null_cva2, 1)   
BasicStats(Null_cva2)
simple_plotLog2Matrix(cva2_m1[cva2_Row_Index,cva2_Col_Index], 1)
simple_plotLog2Matrix(cva2_m2[cva2_Row_Index,cva2_Col_Index], 1)
simple_plotLog2Matrix(cva2_m3[cva2_Row_Index,cva2_Col_Index], 1)
CVA_2_is2 = InteractionScores(cva2_m1,cva2_m2,1.0,filename="cva_2_is2.pdf",cellwidth = 10,cellheight = 10) 
CVA_2_is3 = InteractionScores(cva2_m1,cva2_m3,1.0,filename="cva_2_is3.pdf",cellwidth = 10,cellheight = 10) 
CorMatrix(CVA_2_is2,CVA_2_is3)		# 0.6565  
nCVA_2_is2 = NewInteractionScores(m1,m2,1.0)
nCVA_2_is3 = NewInteractionScores(m1,m3,1.0)
CorMatrix(nCVA_2_is2,nCVA_2_is3)	# 0.7527  
ocva2_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_+His.tsv.original.txt",0.5)
ocva2_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_-His.tsv.original.txt",0.5)
ocva2_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_+3AT.tsv.original.txt",0.5)
OCVA_2_is2 = InteractionScores(ocva2_m1,ocva2_m2,1.0)  
OCVA_2_is3 = InteractionScores(ocva2_m1,ocva2_m3,1.0)  
CorMatrix(OCVA_2_is2,OCVA_2_is3)		#                          
nOCVA_2_is2 = NewInteractionScores(ocva2_m1,ocva2_m2,1.0)  
nOCVA_2_is3 = NewInteractionScores(ocva2_m1,ocva2_m3,1.0)  
CorMatrix(nOCVA_2_is2,nOCVA_2_is3)	#
                                                                    


#===========================
# CV_1
# 1145 X 1345
#
cv1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_+His.tsv.refined.txt",0.5) # read sum = 43064  
cv1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_-His.tsv.refined.txt",0.5) # read sum = 34318  
cv1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_+3AT.tsv.refined.txt",0.5) # read sum = 22932    
simple_plotLog2Matrix(cv1_m1,1,7,filename="cv1_m1.pdf");

cv1_Row_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,26,28,29,30,31,32,33,35,36,37,39,41,42,45,46,47,48,50,51,52,53,54,56,57,58,59,60,61,62,64,66,67,70,71,72,73,74 );
cv1_Col_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,24,26,28,29,31,32,33,35,36,37,38,41,42,45,46,47,48,50,51,52,53,56,57,58,60,61,62,63,64,66,67,68,70,71,72,74 );
Null_cv1 = cv1_m1 * 0 #matrix(0,nrow=78,ncol=78)
cv1_Row_Not_Index = (1:78)[-cv1_Row_Index]
cv1_Col_Not_Index = (1:78)[-cv1_Col_Index]   
Null_cv1[cv1_Row_Not_Index,] = 1     
Null_cv1[,cv1_Col_Not_Index] = 1  
simple_plotMatrix(Null_cv1, 1,7,filename="cv1_m1_null.pdf");         
BasicStats(Null_cv1)
simple_plotLog2Matrix(cv1_m1[cv1_Row_Index,cv1_Col_Index], 1)
simple_plotLog2Matrix(cv1_m2[cv1_Row_Index,cv1_Col_Index], 1)
simple_plotLog2Matrix(cv1_m3[cv1_Row_Index,cv1_Col_Index], 1)
CV_1_is2 = InteractionScores(cv1_m1,cv1_m2,1.0)  
CV_1_is3 = InteractionScores(cv1_m1,cv1_m3,1.0)  
CV_1_is2 = InteractionScores(cv1_m1,cv1_m2,1.0,7,filename="cv_1_is2.pdf",cellwidth = 10,cellheight = 10)  
CV_1_is3 = InteractionScores(cv1_m1,cv1_m3,1.0,7,filename="cv_1_is3.pdf",cellwidth = 10,cellheight = 10)
CorMatrix(CV_1_is2,CV_1_is3)		# 0.7204                         
nCV_1_is2 = NewInteractionScores(cv1_m1,cv1_m2,1.0)  
nCV_1_is3 = NewInteractionScores(cv1_m1,cv1_m3,1.0)  
CorMatrix(nCV_1_is2,nCV_1_is3)		# 0.8083                              
ocv1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_+His.tsv.original.txt",0.5)
ocv1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_-His.tsv.original.txt",0.5)
ocv1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_+3AT.tsv.original.txt",0.5)
OCV_1_is2 = InteractionScores(ocv1_m1,ocv1_m2,1.0)  
OCV_1_is3 = InteractionScores(ocv1_m1,ocv1_m3,1.0)  
CorMatrix(OCV_1_is2,OCV_1_is3)		#                          
nOCV_1_is2 = NewInteractionScores(ocv1_m1,ocv1_m2,1.0)  
nOCV_1_is3 = NewInteractionScores(ocv1_m1,ocv1_m3,1.0)  
CorMatrix(nOCV_1_is2,nOCV_1_is3)	#                              


#########################
# Read BFG_Y2H (All)
#########################

# 325x316 53.1% detection for non-selection
ccc1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_1_+His.tsv.original.txt",0.5)
ccc1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_1_-His.tsv.original.txt",0.5)
ccc1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_1_+3AT.tsv.original.txt",0.5)
ccc1_is2 = InteractionScores(ccc1_m1,ccc1_m2,1.0)  
ccc1_is3 = InteractionScores(ccc1_m1,ccc1_m3,1.0)  
CorMatrix(is2,is3)		# 0.6613                    
ccc1_nis2 = NewInteractionScores(ccc1_m1,ccc1_m2,1.0)  
ccc1_nis3 = NewInteractionScores(ccc1_m1,ccc1_m3,1.0)  
CorMatrix(ccc1_nis2,ccc1_nis3)		# 0.7027                             

                                                                          
# 325x316
ccc2_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_2_+His.tsv.original.txt",0.5)
ccc2_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_2_-His.tsv.original.txt",0.5)
ccc2_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CCC_2_+3AT.tsv.original.txt",0.5)                                                                                                                             
ccc2_is2 = InteractionScores(ccc2_m1,ccc2_m2,1.0)  
ccc2_is3 = InteractionScores(ccc2_m1,ccc2_m3,1.0)  
CorMatrix(ccc2_is2,ccc2_is3)		# 0.673                   
ccc2_nis2 = NewInteractionScores(ccc2_m1,ccc2_m2,1.0)  
ccc2_nis3 = NewInteractionScores(ccc2_m1,ccc2_m3,1.0)  
CorMatrix(ccc2_nis2,ccc2_nis3)		# 0.7188
         
# comparison between replicate experiments (CCC_1 vs CCC_2)
CorMatrix(ccc1_is2, ccc2_is2)	# 0.677
CorMatrix(ccc1_is3, ccc2_is3)	# 0.7115
CorMatrix(ccc1_nis2, ccc2_nis2)	# 0.7303
CorMatrix(ccc1_nis3, ccc2_nis3)	# 0.7594

                                        
# 129x139                                                                                          
cent_aa_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_-AA_+His.tsv.original.txt",0.5)
cent_aa_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_-AA_-His.tsv.original.txt",0.5)
cent_aa_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_-AA_+3AT.tsv.original.txt",0.5)
cent_aa_is2 = InteractionScores(cent_aa_m1,cent_aa_m2,1.0)  
cent_aa_is3 = InteractionScores(cent_aa_m1,cent_aa_m3,1.0)  
CorMatrix(cent_aa_is2,cent_aa_is3)		# 0.5301
cent_aa_nis2 = NewInteractionScores(cent_aa_m1,cent_aa_m2,1.0)  
cent_aa_nis3 = NewInteractionScores(cent_aa_m1,cent_aa_m3,1.0)  
CorMatrix(cent_aa_nis2,cent_aa_nis3)		# 0.4507

# 129x139
cent_all_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_ALL_+His.tsv.original.txt",0.5)
cent_all_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_ALL_-His.tsv.original.txt",0.5)
cent_all_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CENT_ALL_+3AT.tsv.original.txt",0.5)
cent_all_is2 = InteractionScores(cent_all_m1,cent_all_m2,1.0)  
cent_all_is3 = InteractionScores(cent_all_m1,cent_all_m3,1.0)  
CorMatrix(cent_all_is2,cent_all_is3)		# 0.7456
cent_all_nis2 = NewInteractionScores(cent_all_m1,cent_all_m2,1.0)  
cent_all_nis3 = NewInteractionScores(cent_all_m1,cent_all_m3,1.0)  
CorMatrix(cent_all_nis2,cent_all_nis3)		# 0.6797             
                                                                                 
# 1600x1632 32.2% detection for non-selection 
cva1_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+His.tsv.original.txt",0.5)
cva1_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_-His.tsv.original.txt",0.5)
cva1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+3AT.tsv.original.txt",0.5)

cva1_is2 = InteractionScores(cva1_m1,cva1_m2,1.0)  
cva1_is3 = InteractionScores(cva1_m1,cva1_m3,1.0)  
CorMatrix(cva1_is2,cva1_is3)		# 0.673                   
cva1_nis2 = NewInteractionScores(cva1_m1,cva1_m2,1.0)  
cva1_nis3 = NewInteractionScores(cva1_m1,cva1_m3,1.0)  
CorMatrix(cva1_nis2,cva1_nis3)		# 0.7188

cva2_m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_+His.tsv.original.txt",0.5)
cva2_m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_-His.tsv.original.txt",0.5)
cva2_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_2_+3AT.tsv.original.txt",0.5)

cva2_is2 = InteractionScores(cva2_m1,cva2_m2,1.0)  
cva2_is3 = InteractionScores(cva2_m1,cva2_m3,1.0)  
CorMatrix(cva2_is2,cva2_is3)		# 0.673                   
cva2_nis2 = NewInteractionScores(cva2_m1,cva2_m2,1.0)  
cva2_nis3 = NewInteractionScores(cva2_m1,cva2_m3,1.0)  
CorMatrix(cva2_nis2,cva2_nis3)		# 0.7188

m1=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_+His.tsv.original.txt",0.5)
m2=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_-His.tsv.original.txt",0.5)
m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CV_1_+3AT.tsv.original.txt",0.5)


../BFG_Y2H/Yachie_Petsalaki_Data_S3/CCC_1_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CCC_1_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CCC_1_-His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CCC_2_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CCC_2_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CCC_2_-His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CENT_-AA_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CENT_-AA_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CENT_-AA_-His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CENT_ALL_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CENT_ALL_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CENT_ALL_-His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CVA_1_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CVA_1_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CVA_1_-His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CVA_2_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CVA_2_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CVA_2_-His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CV_1_+3AT.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CV_1_+His.tsv
../BFG_Y2H/Yachie_Petsalaki_Data_S3/CV_1_-His.tsv
          

# Performance Test 
boxplot( is7_3[ which( BioGrid == 1 ) ], is7_3[ which( BioGrid == 0 ) ] ) 
t.test( is7_3[ which( BioGrid == 1 ) ], is7_3[ which( BioGrid == 0 ) ] )

boxplot( is2_3[ which( BioGrid == 1 ) ], is2_3[ which( BioGrid == 0 ) ] ) 
t.test( is2_3[ which( BioGrid == 1 ) ], is2_3[ which( BioGrid == 0 ) ] )

boxplot( is3_3[ which( BioGrid == 1 ) ], is3_3[ which( BioGrid == 0 ) ] ) 
t.test( is3_3[ which( BioGrid == 1 ) ], is3_3[ which( BioGrid == 0 ) ] )


boxplot( is7_3[ which( HIPPIE == 1 ) ], is7_3[ which( HIPPIE == 0 ) ] ) 
t.test( is7_3[ which( HIPPIE == 1 ) ], is7_3[ which( HIPPIE == 0 ) ] )

boxplot( is2_3[ which( HIPPIE == 1 ) ], is2_3[ which( HIPPIE == 0 ) ] ) 
t.test( is2_3[ which( HIPPIE == 1 ) ], is2_3[ which( HIPPIE == 0 ) ] )    

boxplot( is3_3[ which( HIPPIE == 1 ) ], is3_3[ which( HIPPIE == 0 ) ] ) 
t.test( is3_3[ which( HIPPIE == 1 ) ], is3_3[ which( HIPPIE == 0 ) ] )
   

is2_sum = is2_1 + is2_2 + is2_3
boxplot( is2_sum[ which( BioGrid == 1 ) ], is2_sum[ which( BioGrid == 0 ) ] ) 
t.test( is2_sum[ which( BioGrid == 1 ) ], is2_sum[ which( BioGrid == 0 ) ] )

boxplot( is2_sum[ which( HIPPIE == 1 ) ], is2_sum[ which( HIPPIE == 0 ) ] ) 
t.test( is2_sum[ which( HIPPIE == 1 ) ], is2_sum[ which( HIPPIE == 0 ) ] )



boxplot( is2[ which( BioGrid == 1 ) ], is2[ which( BioGrid == 0 ) ] ) 
t.test( is2[ which( BioGrid == 1 ) ], is2[ which( BioGrid == 0 ) ] )

boxplot( is2[ which( HIPPIE == 1 ) ], is2[ which( HIPPIE == 0 ) ] ) 
t.test( is2[ which( HIPPIE == 1 ) ], is2[ which( HIPPIE == 0 ) ] )   

boxplot( is3[ which( BioGrid == 1 ) ], is3[ which( BioGrid == 0 ) ] ) 
t.test( is3[ which( BioGrid == 1 ) ], is3[ which( BioGrid == 0 ) ] )

boxplot( is3[ which( HIPPIE == 1 ) ], is3[ which( HIPPIE == 0 ) ] ) 
t.test( is3[ which( HIPPIE == 1 ) ], is3[ which( HIPPIE == 0 ) ] )
                                                                       

## Comparison between data
data1 = data.frame(Exp1=as.vector(log2(rm1_1+1)),Exp2=as.vector(log2(rm1_2+1)),Exp3=as.vector(log2(rm1_3+1)))
png( "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S1.comparison.png" )   
#pairs(data1)                                                                                                    
pairs(data1,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (no selection media)")  
dev.off()
cor(data1)
#           first    second     third
#first  1.0000000 0.8333677 0.8255604
#second 0.8333677 1.0000000 0.8872547
#third  0.8255604 0.8872547 1.0000000
data2 = data.frame(first=as.vector(log2(rm2_1+1)),second=as.vector(log2(rm2_2+1)),third=as.vector(log2(rm2_3+1)))       
png( "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S2.comparison.png" ) 
pairs(data2)  
dev.off()   
cor(data2)
#           first    second     third
#first  1.0000000 0.5758140 0.6732435
#second 0.5758140 1.0000000 0.7681709
#third  0.6732435 0.7681709 1.0000000
data3 = data.frame(Exp1A=as.vector(log2(rm2_1+1)),Exp2A=as.vector(log2(rm2_2+1)),Exp3A=as.vector(log2(rm2_3+1)),Exp2Q=as.vector(log2(rm3_2+1)),Exp3Q=as.vector(log2(rm3_3+1)))    
png( "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S2_3.comparison.png" )     
#pairs(data3)
pairs(data3,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (selection media)")   
cor(data3) 
dev.off()  



roth=read.csv("~/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/src/roth2016_ref.txt",sep="\t",header=TRUE, comment.char = "#")
View(roth)

roth_pos = roth[which(roth$Interaction.re.annotation == "Interaction"),] 
roth_neg = roth[which(roth$Interaction.re.annotation == "No Interaction"),]              

dim(roth_pos)
dim(roth_neg)

rm1_1[roth_pos$rowindex,roth_pos$colindex]

for (i in 1:35){                              
	#print( roth_pos[i,1] )
	#print( roth_pos[i,2] )
	print( rm1_1[roth_pos[i,6],roth_pos[i,7]])
}

for (i in 1:23){                              
	#print( roth_pos[i,1] )
	#print( roth_pos[i,2] )
	print( rm1_1[roth_neg[i,6],roth_neg[i,7]])
}                         

rm1_1[ cbind( roth_pos[,6],roth_pos[,7] ) ]
rm1_1[ cbind( roth_neg[,6],roth_neg[,7] ) ]

par(mfrow=c(2,4)) 
# Roth2016_ref_read_counts                                                                                                                  
# Exp1
boxplot( log2(rm1_1[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm1_1[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp1 no selection" ) # no selection
boxplot( log2(rm2_1[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm2_1[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp1 selection 1" ) # Aba
View( cbind( roth_neg, rm2_1[ cbind( roth_neg[,6],roth_neg[,7] ) ]) )    
      
# Exp2
boxplot( log2(rm1_2[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm1_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp2 no selection" ) # no selection    
boxplot( log2(rm2_2[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm2_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp2 selection 1" ) # Aba 
boxplot( log2(rm3_2[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm3_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp2 selection 2" ) # Q 
View( cbind( roth_neg, rm3_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) )    

# Exp3
boxplot( log2(rm1_3[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm1_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp3 no selection" ) # no selection          
boxplot( log2(rm2_3[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm2_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp3 selection 1" ) # Aba 
boxplot( log2(rm3_3[ cbind( roth_pos[,6],roth_pos[,7] ) ]+1), log2(rm3_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]+1), main="Exp3 selection 2" ) # Q 
View( cbind( roth_neg, rm3_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]) )   


par(mfrow=c(2,4)) 
# Roth2016_ref_IS                                                                                                                  
# Exp1
boxplot( ris2_1[ cbind( roth_pos[,6],roth_pos[,7] ) ], ris2_1[ cbind( roth_neg[,6],roth_neg[,7] ) ], main="Exp1 selection 1" ) # Aba
cut2_1 = mean(ris2_1[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 5.0            
cut2_1 = max(ris2_1[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 2.0   
par(mfrow=c(1,1))
visweb(ris2_1> cut2_1, type="none" )                                    

# Exp2
boxplot( ris2_2[ cbind( roth_pos[,6],roth_pos[,7] ) ], ris2_2[ cbind( roth_neg[,6],roth_neg[,7] ) ], main="Exp2 selection 1" ) # Aba 
boxplot( ris3_2[ cbind( roth_pos[,6],roth_pos[,7] ) ], ris3_2[ cbind( roth_neg[,6],roth_neg[,7] ) ], main="Exp2 selection 2" ) # Q  
cut2_2 = mean(ris2_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 5.0         
cut3_2 = mean(ris3_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 5.0 
cut2_2 = max(ris2_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 2.0         
cut3_2 = max(ris3_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 2.0 
par(mfrow=c(1,2))
visweb(ris2_2> cut2_2, type="none" )  
visweb(ris3_2> cut3_2, type="none" )
View( cbind( roth_neg, ris2_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) )       
View( cbind( roth_neg, ris3_2[ cbind( roth_neg[,6],roth_neg[,7] ) ]) )       


# Exp3
boxplot( ris2_3[ cbind( roth_pos[,6],roth_pos[,7] ) ], ris2_3[ cbind( roth_neg[,6],roth_neg[,7] ) ], main="Exp3 selection 1" ) # Aba 
boxplot( ris3_3[ cbind( roth_pos[,6],roth_pos[,7] ) ], ris3_3[ cbind( roth_neg[,6],roth_neg[,7] ) ], main="Exp3 selection 2" ) # Q
cut2_3 = mean(ris2_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 5.0         
cut3_3 = mean(ris3_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 5.0   
cut2_3 = max(ris2_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 2.0         
cut3_3 = max(ris3_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]) * 2.0
par(mfrow=c(1,2))
visweb(ris2_3> cut2_3, type="none" )  
visweb(ris3_3> cut3_3, type="none" )  
View( cbind( roth_neg, ris3_3[ cbind( roth_neg[,6],roth_neg[,7] ) ]) )


# Count times detected
ris_call_sum = (ris2_1>cut2_1) + (ris2_2>cut2_2) + (ris3_2>cut3_2) + (ris2_3>cut2_3) + (ris3_3>cut3_3)

pheatmap(ris_call_sum, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))    

bsum = baseline(sum)                                                                      
pheatmap(bsum, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))    
             
melted_ris_call_sum = melt(ris_call_sum) 
write.csv(melted_ris_call_sum,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/ris_call_sum.txt")        

                       
rHIPPIE=RefineMatrix(HIPPIE)
rBioGrid=RefineMatrix(BioGrid)
melted_HIPPIE = melt(rHIPPIE)
melted_BioGrid = melt(rBioGrid)    

write.csv(melted_HIPPIE,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_HIPPIE.txt")  
write.csv(melted_BioGrid,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_BioGrid.txt")     

melted_CCC_1_is2 = melt(CCC_1_is2)
write.csv(melted_CCC_1_is2,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CCC_1_is2.txt")   

melted_CCC_1_is3 = melt(CCC_1_is3)
write.csv(melted_CCC_1_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CCC_1_is3.txt")
                                  
melted_CCC_2_is2 = melt(CCC_2_is2)
write.csv(melted_CCC_2_is2,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CCC_2_is2.txt")   

melted_CCC_2_is3 = melt(CCC_2_is3)
write.csv(melted_CCC_2_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CCC_2_is3.txt")

melted_CENT_AA_is2 = melt(CENT_AA_is2)
write.csv(melted_CENT_AA_is2,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CENT_AA_is2.txt")

melted_CENT_AA_is3 = melt(CENT_AA_is3)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CENT_AA_is3.txt")

melted_CENT_ALL_is2 = melt(CENT_ALL_is2)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CENT_ALL_is2.txt")

melted_CENT_ALL_is3 = melt(CENT_ALL_is3)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CENT_ALL_is3.txt")

melted_CVA_1_is2 = melt(CVA_1_is2)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CVA_1_is2.txt")

melted_CVA_1_is3 = melt(CVA_1_is3)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CVA_1_is3.txt")

melted_CVA_2_is2 = melt(CVA_2_is2)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CVA_2_is2.txt")

melted_CVA_2_is3 = melt(CVA_2_is3)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CVA_2_is3.txt")

melted_CV_1_is2 = melt(CV_1_is2)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CV_1_is2.txt")

melted_CV_1_is3 = melt(CV_1_is3)
write.csv(melted_CENT_AA_is3,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/melted_CV_1_is3.txt")






############# codon usage ####################         
# actually human genes we are using are sub-optimal for yeast. For human average codon usage they are 0.774 and for yeast 0.603 (relative maximum 1.0) and its significant difference (P = 7.53E-32).
# However, there were no difference between protein pairs found in Y2H and didn't find in Y2H

library(readxl)                    
SpotTestMGj54_JS_R75 <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "R75_Final_List") 
R75_Human_CAI = SpotTestMGj54_JS_R75[1:76,16]
R75_Yeast_CAI = SpotTestMGj54_JS_R75[1:76,17]
R75_Ecoli_CAI = SpotTestMGj54_JS_R75[1:76,18]
boxplot( R75_Human_CAI, R75_Yeast_CAI, R75_Ecoli_CAI, names = c("CAI (Human)","CAI (Yeast)", "CAI (E.coli)"), ylab="Codon adaptation index (CAI)", main="CAI for tested proteins" )

R75_gene_length = SpotTestMGj54_JS_R75[1:76,12]   
R75_gene_name = SpotTestMGj54_JS_R75[1:76,14]

> SpotTestMGj54_JS <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "CAI_yeast_summary")

par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS[1:65,5], SpotTestMGj54_JS[66:185,5], names = c("Found in Y2H","No found in Y2H"), ylab="CAI (yeast)", main="CAI for protein 1" )
boxplot( SpotTestMGj54_JS[1:65,6], SpotTestMGj54_JS[66:185,6], names = c("Found in Y2H","No found in Y2H"), ylab="CAI (yeast)", main="CAI for protein 2" )
boxplot( SpotTestMGj54_JS[1:65,7], SpotTestMGj54_JS[66:185,7], names = c("Found in Y2H","No found in Y2H"), ylab="CAI (yeast)", main="CAI p1 x CAI p2" )
boxplot( SpotTestMGj54_JS[1:65,8], SpotTestMGj54_JS[66:185,8], names = c("Found in Y2H","No found in Y2H"), ylab="CAI (yeast)", main="Average CAI for p1 and p2" )
boxplot( SpotTestMGj54_JS[1:65,9], SpotTestMGj54_JS[66:185,9], names = c("Found in Y2H","No found in Y2H"), ylab="CAI (yeast)", main="Min( CAI for p1 and p2)" )
boxplot( SpotTestMGj54_JS[1:65,10], SpotTestMGj54_JS[66:185,10], names = c("Found in Y2H","No found in Y2H"), ylab="CAI (yeast)", main="Max( CAI for p1 and p2)" )         


SpotTestMGj54_JS_K_diff <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "AA_K_differences")
par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_K_diff[1:65,5], SpotTestMGj54_JS_K_diff[66:185,5], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of K", main="K (Freq) for protein 1" )
boxplot( SpotTestMGj54_JS_K_diff[1:65,6], SpotTestMGj54_JS_K_diff[66:185,6], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of K", main="K (Freq) for protein 2" )
boxplot( SpotTestMGj54_JS_K_diff[1:65,7], SpotTestMGj54_JS_K_diff[66:185,7], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of K", main="K (Freq) p1 x K (Freq) p2" )
boxplot( SpotTestMGj54_JS_K_diff[1:65,8], SpotTestMGj54_JS_K_diff[66:185,8], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of K", main="Average K (Freq) for p1 and p2" )
boxplot( SpotTestMGj54_JS_K_diff[1:65,9], SpotTestMGj54_JS_K_diff[66:185,9], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of K", main="Min( K (Freq) for p1 and p2)" )
boxplot( SpotTestMGj54_JS_K_diff[1:65,10], SpotTestMGj54_JS_K_diff[66:185,10], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of K", main="Max( K (Freq) for p1 and p2)" )
                                                                                      

SpotTestMGj54_JS_Y_diff <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "AA_Y_differences")
par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_Y_diff[1:65,5], SpotTestMGj54_JS_Y_diff[66:185,5], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Y (Freq) for protein 1" )
boxplot( SpotTestMGj54_JS_Y_diff[1:65,6], SpotTestMGj54_JS_Y_diff[66:185,6], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Y (Freq) for protein 2" )
boxplot( SpotTestMGj54_JS_Y_diff[1:65,7], SpotTestMGj54_JS_Y_diff[66:185,7], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Y (Freq) p1 x Y (Freq) p2" )
boxplot( SpotTestMGj54_JS_Y_diff[1:65,8], SpotTestMGj54_JS_Y_diff[66:185,8], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Average Y (Freq) for p1 and p2" )
boxplot( SpotTestMGj54_JS_Y_diff[1:65,9], SpotTestMGj54_JS_Y_diff[66:185,9], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Min( Y (Freq) for p1 and p2)" )
boxplot( SpotTestMGj54_JS_Y_diff[1:65,10], SpotTestMGj54_JS_Y_diff[66:185,10], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Max( Y (Freq) for p1 and p2)" )

par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_Y_diff[1:25,5], SpotTestMGj54_JS_Y_diff[26:50,5], SpotTestMGj54_JS_Y_diff[51:65,5], SpotTestMGj54_JS_Y_diff[66:79,5], SpotTestMGj54_JS_Y_diff[80:185,5], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Y (Freq) for protein 1" ) 
boxplot( SpotTestMGj54_JS_Y_diff[1:25,6], SpotTestMGj54_JS_Y_diff[26:50,6], SpotTestMGj54_JS_Y_diff[51:65,6], SpotTestMGj54_JS_Y_diff[66:79,6], SpotTestMGj54_JS_Y_diff[80:185,6], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Y (Freq) for protein 2" )               
boxplot( SpotTestMGj54_JS_Y_diff[1:25,7], SpotTestMGj54_JS_Y_diff[26:50,7], SpotTestMGj54_JS_Y_diff[51:65,7], SpotTestMGj54_JS_Y_diff[66:79,7], SpotTestMGj54_JS_Y_diff[80:185,7], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Y (Freq) p1 x Y (Freq) p2" )   
boxplot( SpotTestMGj54_JS_Y_diff[1:25,8], SpotTestMGj54_JS_Y_diff[26:50,8], SpotTestMGj54_JS_Y_diff[51:65,8], SpotTestMGj54_JS_Y_diff[66:79,8], SpotTestMGj54_JS_Y_diff[80:185,8], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Average Y (Freq) for p1 and p2" )    
boxplot( SpotTestMGj54_JS_Y_diff[1:25,9], SpotTestMGj54_JS_Y_diff[26:50,9], SpotTestMGj54_JS_Y_diff[51:65,9], SpotTestMGj54_JS_Y_diff[66:79,9], SpotTestMGj54_JS_Y_diff[80:185,9], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Min( Y (Freq) for p1 and p2)" )    
boxplot( SpotTestMGj54_JS_Y_diff[1:25,10], SpotTestMGj54_JS_Y_diff[26:50,10], SpotTestMGj54_JS_Y_diff[51:65,10], SpotTestMGj54_JS_Y_diff[66:79,10], SpotTestMGj54_JS_Y_diff[80:185,10], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Max( Y (Freq) for p1 and p2)" )           


SpotTestMGj54_JS_P_diff <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "AA_P_differences")
par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_P_diff[1:65,5], SpotTestMGj54_JS_P_diff[66:185,5], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="P (Freq) for protein 1" )
boxplot( SpotTestMGj54_JS_P_diff[1:65,6], SpotTestMGj54_JS_P_diff[66:185,6], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="P (Freq) for protein 2" )
boxplot( SpotTestMGj54_JS_P_diff[1:65,7], SpotTestMGj54_JS_P_diff[66:185,7], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="P (Freq) p1 x P (Freq) p2" )
boxplot( SpotTestMGj54_JS_P_diff[1:65,8], SpotTestMGj54_JS_P_diff[66:185,8], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Average P (Freq) for p1 and p2" )
boxplot( SpotTestMGj54_JS_P_diff[1:65,9], SpotTestMGj54_JS_P_diff[66:185,9], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Min( P (Freq) for p1 and p2)" )
boxplot( SpotTestMGj54_JS_P_diff[1:65,10], SpotTestMGj54_JS_P_diff[66:185,10], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Max( P (Freq) for p1 and p2)" )

par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_P_diff[1:25,5], SpotTestMGj54_JS_P_diff[26:50,5], SpotTestMGj54_JS_P_diff[51:65,5], SpotTestMGj54_JS_P_diff[66:79,5], SpotTestMGj54_JS_P_diff[80:185,5], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="P (Freq) for protein 1" ) 
boxplot( SpotTestMGj54_JS_P_diff[1:25,6], SpotTestMGj54_JS_P_diff[26:50,6], SpotTestMGj54_JS_P_diff[51:65,6], SpotTestMGj54_JS_P_diff[66:79,6], SpotTestMGj54_JS_P_diff[80:185,6], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="P (Freq) for protein 2" )               
boxplot( SpotTestMGj54_JS_P_diff[1:25,7], SpotTestMGj54_JS_P_diff[26:50,7], SpotTestMGj54_JS_P_diff[51:65,7], SpotTestMGj54_JS_P_diff[66:79,7], SpotTestMGj54_JS_P_diff[80:185,7], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="P (Freq) p1 x P (Freq) p2" )   
boxplot( SpotTestMGj54_JS_P_diff[1:25,8], SpotTestMGj54_JS_P_diff[26:50,8], SpotTestMGj54_JS_P_diff[51:65,8], SpotTestMGj54_JS_P_diff[66:79,8], SpotTestMGj54_JS_P_diff[80:185,8], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Average P (Freq) for p1 and p2" )    
boxplot( SpotTestMGj54_JS_P_diff[1:25,9], SpotTestMGj54_JS_P_diff[26:50,9], SpotTestMGj54_JS_P_diff[51:65,9], SpotTestMGj54_JS_P_diff[66:79,9], SpotTestMGj54_JS_P_diff[80:185,9], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Min( P (Freq) for p1 and p2)" )    
boxplot( SpotTestMGj54_JS_P_diff[1:25,10], SpotTestMGj54_JS_P_diff[26:50,10], SpotTestMGj54_JS_P_diff[51:65,10], SpotTestMGj54_JS_P_diff[66:79,10], SpotTestMGj54_JS_P_diff[80:185,10], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Max( P (Freq) for p1 and p2)" )           


SpotTestMGj54_JS_D_diff <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "AA_D_differences")
par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_D_diff[1:65,5], SpotTestMGj54_JS_D_diff[66:185,5], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="D (Freq) for protein 1" )
boxplot( SpotTestMGj54_JS_D_diff[1:65,6], SpotTestMGj54_JS_D_diff[66:185,6], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="D (Freq) for protein 2" )
boxplot( SpotTestMGj54_JS_D_diff[1:65,7], SpotTestMGj54_JS_D_diff[66:185,7], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="D (Freq) p1 x D (Freq) p2" )
boxplot( SpotTestMGj54_JS_D_diff[1:65,8], SpotTestMGj54_JS_D_diff[66:185,8], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Average D (Freq) for p1 and p2" )
boxplot( SpotTestMGj54_JS_D_diff[1:65,9], SpotTestMGj54_JS_D_diff[66:185,9], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Min( D (Freq) for p1 and p2)" )
boxplot( SpotTestMGj54_JS_D_diff[1:65,10], SpotTestMGj54_JS_D_diff[66:185,10], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Max( D (Freq) for p1 and p2)" )

par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_D_diff[1:25,5], SpotTestMGj54_JS_D_diff[26:50,5], SpotTestMGj54_JS_D_diff[51:65,5], SpotTestMGj54_JS_D_diff[66:79,5], SpotTestMGj54_JS_D_diff[80:185,5], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="D (Freq) for protein 1" ) 
boxplot( SpotTestMGj54_JS_D_diff[1:25,6], SpotTestMGj54_JS_D_diff[26:50,6], SpotTestMGj54_JS_D_diff[51:65,6], SpotTestMGj54_JS_D_diff[66:79,6], SpotTestMGj54_JS_D_diff[80:185,6], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="D (Freq) for protein 2" )               
boxplot( SpotTestMGj54_JS_D_diff[1:25,7], SpotTestMGj54_JS_D_diff[26:50,7], SpotTestMGj54_JS_D_diff[51:65,7], SpotTestMGj54_JS_D_diff[66:79,7], SpotTestMGj54_JS_D_diff[80:185,7], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="D (Freq) p1 x D (Freq) p2" )   
boxplot( SpotTestMGj54_JS_D_diff[1:25,8], SpotTestMGj54_JS_D_diff[26:50,8], SpotTestMGj54_JS_D_diff[51:65,8], SpotTestMGj54_JS_D_diff[66:79,8], SpotTestMGj54_JS_D_diff[80:185,8], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Average D (Freq) for p1 and p2" )    
boxplot( SpotTestMGj54_JS_D_diff[1:25,9], SpotTestMGj54_JS_D_diff[26:50,9], SpotTestMGj54_JS_D_diff[51:65,9], SpotTestMGj54_JS_D_diff[66:79,9], SpotTestMGj54_JS_D_diff[80:185,9], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Min( D (Freq) for p1 and p2)" )    
boxplot( SpotTestMGj54_JS_D_diff[1:25,10], SpotTestMGj54_JS_D_diff[26:50,10], SpotTestMGj54_JS_D_diff[51:65,10], SpotTestMGj54_JS_D_diff[66:79,10], SpotTestMGj54_JS_D_diff[80:185,10], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Max( D (Freq) for p1 and p2)" )           


SpotTestMGj54_JS_KRDE_diff <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "AA_KRDE_differences")
par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_KRDE_diff[1:65,5], SpotTestMGj54_JS_KRDE_diff[66:185,5], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Charged AA (Freq) for protein 1" )
boxplot( SpotTestMGj54_JS_KRDE_diff[1:65,6], SpotTestMGj54_JS_KRDE_diff[66:185,6], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Charged AA (Freq) for protein 2" )
boxplot( SpotTestMGj54_JS_KRDE_diff[1:65,7], SpotTestMGj54_JS_KRDE_diff[66:185,7], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Charged AA (Freq) p1 x Charged AA (Freq) p2" )
boxplot( SpotTestMGj54_JS_KRDE_diff[1:65,8], SpotTestMGj54_JS_KRDE_diff[66:185,8], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Average Charged AA (Freq) for p1 and p2" )
boxplot( SpotTestMGj54_JS_KRDE_diff[1:65,9], SpotTestMGj54_JS_KRDE_diff[66:185,9], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Min( Charged AA (Freq) for p1 and p2)" )
boxplot( SpotTestMGj54_JS_KRDE_diff[1:65,10], SpotTestMGj54_JS_KRDE_diff[66:185,10], names = c("Found in Y2H","No found in Y2H"), ylab="Frequency of Y", main="Max( Charged AA (Freq) for p1 and p2)" )

par(mfrow=c(2,3))
boxplot( SpotTestMGj54_JS_KRDE_diff[1:25,5], SpotTestMGj54_JS_KRDE_diff[26:50,5], SpotTestMGj54_JS_KRDE_diff[51:65,5], SpotTestMGj54_JS_KRDE_diff[66:79,5], SpotTestMGj54_JS_KRDE_diff[80:185,5], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Charged AA (Freq) for protein 1" ) 
boxplot( SpotTestMGj54_JS_KRDE_diff[1:25,6], SpotTestMGj54_JS_KRDE_diff[26:50,6], SpotTestMGj54_JS_KRDE_diff[51:65,6], SpotTestMGj54_JS_KRDE_diff[66:79,6], SpotTestMGj54_JS_KRDE_diff[80:185,6], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Charged AA (Freq) for protein 2" )               
boxplot( SpotTestMGj54_JS_KRDE_diff[1:25,7], SpotTestMGj54_JS_KRDE_diff[26:50,7], SpotTestMGj54_JS_KRDE_diff[51:65,7], SpotTestMGj54_JS_KRDE_diff[66:79,7], SpotTestMGj54_JS_KRDE_diff[80:185,7], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Charged AA (Freq) p1 x Charged AA (Freq) p2" )   
boxplot( SpotTestMGj54_JS_KRDE_diff[1:25,8], SpotTestMGj54_JS_KRDE_diff[26:50,8], SpotTestMGj54_JS_KRDE_diff[51:65,8], SpotTestMGj54_JS_KRDE_diff[66:79,8], SpotTestMGj54_JS_KRDE_diff[80:185,8], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Average Charged AA (Freq) for p1 and p2" )    
boxplot( SpotTestMGj54_JS_KRDE_diff[1:25,9], SpotTestMGj54_JS_KRDE_diff[26:50,9], SpotTestMGj54_JS_KRDE_diff[51:65,9], SpotTestMGj54_JS_KRDE_diff[66:79,9], SpotTestMGj54_JS_KRDE_diff[80:185,9], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Min( Charged AA (Freq) for p1 and p2)" )    
boxplot( SpotTestMGj54_JS_KRDE_diff[1:25,10], SpotTestMGj54_JS_KRDE_diff[26:50,10], SpotTestMGj54_JS_KRDE_diff[51:65,10], SpotTestMGj54_JS_KRDE_diff[66:79,10], SpotTestMGj54_JS_KRDE_diff[80:185,10], names = c("5 times","4 times","3 times", "2 ~ 1 times", "not detected"), ylab="Frequency of Y", main="Max( Charged AA (Freq) for p1 and p2)" )


#############################################################################
#
# 2017-06-08

m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.refine.txt",0.5) # BD P170-4 x AD empty; -W; R75_01 ~ R75_07 spike in                      "sum = 6489"
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.refine.txt",0.5) # BD P170-4 x AD empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in        "sum = 7337"
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.refine.txt",0.5) # BD Empty x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "sum = 2192"      

m4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S52.ppi.txt.refine.txt",0.5) # BD R68 x AD R73; -W; Seaprep gel; R75_01 ~ R75_07 spike in              "sum = 2656366"

m5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                   "sum = 3118284"
m6 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in     "sum = 3426672"
is6 = InteractionScores(m5,m6,1.0,4) 
nis6 = NewInteractionScores(m5,m6,1.0,4) 
                                                                 
                                                                                                                   
## DATA generated by blastn-short
bm1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S49.ppi.txt",0.5)  	# ==> 9299 
bm2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S50.ppi.txt",0.5)    # ==> 7991
bm3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S51.ppi.txt",0.5)    # ==> 3965
bm4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt",0.5)    # ==> 2801344

bm5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S53.ppi.txt",0.5)    # ==> 3250858
bm6 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S54.ppi.txt",0.5)    # ==> 3553164
bis6 = InteractionScores(bm5,bm6,1.0,4) 
bnis6 = NewInteractionScores(bm5,bm6,1.0,4)                                 

# # BD R68 x AD R73; -W; Seaprep gel; R75_01 ~ R75_07 spike in  
> BasicStats(m4)
[1] "BD cnt = 78"
[1] "AD cnt = 78"
[1] "ALL cnt = 6084"
[1] "ZERO cnt = 1665"
[1] "% detected cnt = 72.6"
[1] "sum = 2656366"
[1] "max read = 7472"
[1] "min read = 0"
[1] "max AD = 169266"
[1] "min AD = 0"
[1] "max BD = 94331"
[1] "max BD = 0"   

rm4 = RefineMatrix(m4) # "% detected cnt = 81.7"
BasicStats(rm4)


> BasicStats(m5)
[1] "BD cnt = 173"
[1] "AD cnt = 173"
[1] "ALL cnt = 29929"
[1] "ZERO cnt = 10270"
[1] "% detected cnt = 65.7"
[1] "sum = 3118284"
[1] "max read = 2601"
[1] "min read = 0"
[1] "max AD = 88144"
[1] "min AD = 0"
[1] "max BD = 59913"
[1] "max BD = 0"         

## Comparison between data2
data1 = data.frame(Exp1=as.vector(log2(rm1_1+1)),Exp2=as.vector(log2(rm1_2+1)),Exp3=as.vector(log2(rm1_3+1)),Exp4=as.vector(log2(rm4+1)))
#png( "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S1.comparison.png" )   
#pairs(data1)                                                                                                    
pairs(data1,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (no selection media)")  
#dev.off()
cor(data1)
#           first    second     third
#first  1.0000000 0.8333677 0.8255604
#second 0.8333677 1.0000000 0.8872547
#third  0.8255604 0.8872547 1.0000000
data2 = data.frame(first=as.vector(log2(rm2_1+1)),second=as.vector(log2(rm2_2+1)),third=as.vector(log2(rm2_3+1)))       
png( "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S2.comparison.png" ) 
pairs(data2)  
dev.off()   
cor(data2)
#           first    second     third
#first  1.0000000 0.5758140 0.6732435
#second 0.5758140 1.0000000 0.7681709
#third  0.6732435 0.7681709 1.0000000
data3 = data.frame(Exp1A=as.vector(log2(rm2_1+1)),Exp2A=as.vector(log2(rm2_2+1)),Exp3A=as.vector(log2(rm2_3+1)),Exp2Q=as.vector(log2(rm3_2+1)),Exp3Q=as.vector(log2(rm3_3+1)))    
png( "/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-03-03_MiSeq_R73x68_MGj53/R71_76.S2_3.comparison.png" )     
#pairs(data3)
pairs(data3,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (selection media)")   
cor(data3) 
dev.off()



# Save file
EB_1 = InteractionScores(m5,m6,1.0,4)     
raw_is_EB = RawInteractionScores(m5,m6,1.0) 
raw_is_EB[which(raw_is_EB=="NaN")] = 0
raw_is_EB[which(raw_is_EB==Inf)] = 0
pheatmap(log2(raw_is_EB+1.0), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"),fontsize=4)                      
melted_raw_is_EB = melt(raw_is_EB)
write.table(melted_raw_is_EB,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-06-08_MiSeq_EB_Par_MGj56/melted_EB.IS.txt",sep="\t")
melted_S53 = melt(m5)
write.table(melted_S53,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-06-08_MiSeq_EB_Par_MGj56/melted_S53.count.txt",sep="\t")
melted_S54 = melt(m6)
write.table(melted_S54,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-06-08_MiSeq_EB_Par_MGj56/melted_S54.count.txt",sep="\t")
                                         
NEB_1 = NewInteractionScores(m5,m6,1.0,4)     


rm4 = RefineMatrix(m4)
simple_plotLog2Matrix(rm4,0.5)
BasicStats(rm4) # --> S52_BasicStats.pdf 81.7% detected  
[1] "BD cnt = 71"
[1] "AD cnt = 76"
[1] "ALL cnt = 5396"
[1] "ZERO cnt = 985"
[1] "% detected cnt = 81.7"
[1] "sum = 2655845"
[1] "max read = 7472"
[1] "min read = 0"
[1] "max AD = 169266"
[1] "min AD = 76"
[1] "max BD = 94331"
[1] "max BD = 3558"

rrm4 = RefineMatrix2(m4)
simple_plotLog2Matrix(rrm4,0.5)
BasicStats(rrm4) # --> S52_BasicStats.pdf "% detected cnt = 94.2"  
                                       
rrm1_1 = RefineMatrix2(m1_1) 
simple_plotLog2Matrix(rrm1_1,0.5)
BasicStats(rrm1_1) # --> S52_BasicStats.pdf "% detected cnt = 96.9"

rrm1_2 = RefineMatrix2(m1_2) 
simple_plotLog2Matrix(rrm1_2,0.5)
BasicStats(rrm1_2) # --> S52_BasicStats.pdf "% detected cnt = 98.0"

rrm1_3 = RefineMatrix2(m1_3) 
simple_plotLog2Matrix(rrm1_3,0.5)
BasicStats(rrm1_3) # --> S52_BasicStats.pdf "% detected cnt = 98.5"
                          

############## Elementwise multi
EB_MUL = EB_1 * EB_2
simple_plotLog2Matrix(EB_MUL,1.0,4)

NEB_MUL = NEB_1 * NEB_2
simple_plotLog2Matrix(NEB_MUL,1.0,4)  

#
# BD P170-4 x AD empty; -W; R75_01 ~ R75_07 spike in    
# Without selection media to define auto-activator  
#
S49.prc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.real.txt",0.5)
S49.prc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.prc_p.txt",0.5)
S49.p_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.p_prc.txt",0.5)
S49.p_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.p_p.txt",0.5)
S49.prc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.prc_v.txt",0.5)
S49.prc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.prc_vrc.txt",0.5)
S49.p_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.p_v.txt",0.5)
S49.p_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.p_vrc.txt",0.5)    
S49.v_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.v_prc.txt",0.5)  
S49.vrc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.vrc_prc.txt",0.5)  
S49.v_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.v_p.txt",0.5)       
S49.vrc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.vrc_p.txt",0.5)   
S49.vrc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.vrc_vrc.txt",0.5)
S49.vrc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.vrc_v.txt",0.5)
S49.v_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.v_vrc.txt",0.5)
S49.v_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S49.ppi.txt.v_v.txt",0.5)

row_sum.S49.prc_prc = rowSums(S49.prc_prc) 
row_sum.S49.prc_p = rowSums(S49.prc_p)
row_sum.S49.prc_v = rowSums(S49.prc_v)   
row_sum.S49.prc_vrc = rowSums(S49.prc_vrc)   

plot(row_sum.S49.prc_prc) # 7337                                            
plot(row_sum.S49.prc_p)   # 257276
plot(row_sum.S49.prc_v)   # 103287
plot(row_sum.S49.prc_vrc) # 945066     
                                    
plot(row_sum.S49.prc_prc/sum(row_sum.S50.prc_prc)*100) # 7337                                            
plot(row_sum.S49.prc_p)   # 257276
plot(row_sum.S49.prc_v)   # 103287
plot(row_sum.S49.prc_vrc) # 945066        

length(row_sum.S49.prc_prc) # 7337                                            
length(row_sum.S49.prc_p)   # 257276
length(row_sum.S49.prc_v)   # 103287
length(row_sum.S49.prc_vrc) # 945066    

PRC	PRC 	6489
PRC	P 	390660
P	PRC 	240
P	P 	83
PRC	V 	130961
PRC	VRC 	1190839
P	V 	1160
P	VRC 	1479
V	PRC 	2124
VRC	PRC 	385
V	P 	22737
VRC	P 	128
VRC	VRC 	35693
VRC	V 	98028
V	VRC 	527442
V	V 	165691
S49.stat = read.clipboard(header = FALSE)

#
# BD P170-4 x AD empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in     
# Selection media to define auto-activator  
#
S50.prc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.real.txt",0.5)
S50.prc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.prc_p.txt",0.5)
S50.p_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.p_prc.txt",0.5)
S50.p_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.p_p.txt",0.5)
S50.prc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.prc_v.txt",0.5)
S50.prc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.prc_vrc.txt",0.5)
S50.p_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.p_v.txt",0.5)
S50.p_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.p_vrc.txt",0.5)    
S50.v_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.v_prc.txt",0.5)  
S50.vrc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.vrc_prc.txt",0.5)  
S50.v_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.v_p.txt",0.5)       
S50.vrc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.vrc_p.txt",0.5)   
S50.vrc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.vrc_vrc.txt",0.5)
S50.vrc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.vrc_v.txt",0.5)
S50.v_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.v_vrc.txt",0.5)
S50.v_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S50.ppi.txt.v_v.txt",0.5)

# S50 final stat
PRC	PRC 	7337
PRC	P 	257276
P	PRC 	143
P	P 	19
PRC	V 	103287
PRC	VRC 	945066
P	V 	335
P	VRC 	296
V	PRC 	1432
VRC	PRC 	175
V	P 	9517
VRC	P 	51
VRC	VRC 	21565
VRC	V 	84193
V	VRC 	381116
V	V 	104481
S50.stat = read.clipboard(header = FALSE)

row_sum.S50.prc_prc = rowSums(S50.prc_prc) 
row_sum.S50.prc_p = rowSums(S50.prc_p)
row_sum.S50.prc_v = rowSums(S50.prc_v)   
row_sum.S50.prc_vrc = rowSums(S50.prc_vrc)   

plot(row_sum.S50.prc_prc) # 7337                                            
plot(row_sum.S50.prc_p)   # 257276
plot(row_sum.S50.prc_v)   # 103287
plot(row_sum.S50.prc_vrc) # 945066     
                                    
plot(row_sum.S50.prc_prc/sum(row_sum.S50.prc_prc)*100) # 7337                                            
plot(row_sum.S50.prc_p)   # 257276
plot(row_sum.S50.prc_v)   # 103287
plot(row_sum.S50.prc_vrc) # 945066        

length(row_sum.S50.prc_prc) # 7337                                            
length(row_sum.S50.prc_p)   # 257276
length(row_sum.S50.prc_v)   # 103287
length(row_sum.S50.prc_vrc) # 945066

plot(row_sum.S50.prc_p,row_sum.S50.prc_vrc)
> cor(row_sum.S50.prc_p,row_sum.S50.prc_vrc)
[1] 0.9999333
> cor(log(row_sum.S50.prc_p,2),log(row_sum.S50.prc_vrc,2))
[1] NaN
> cor(log(row_sum.S50.prc_p+1,2),log(row_sum.S50.prc_vrc+1,2))
[1] 0.9426546
> plot(log(row_sum.S50.prc_p+1,2),log(row_sum.S50.prc_vrc+1,2))

 
row_sum.S50.prc_prc[as.vector(row_sum.S50.prc_prc/sum(row_sum.S50.prc_prc)*100)>0.6097561]
row_sum.S50.prc_p[as.vector(row_sum.S50.prc_p/sum(row_sum.S50.prc_p)*100)>0.6097561]
row_sum.S50.prc_v[as.vector(row_sum.S50.prc_v/sum(row_sum.S50.prc_v)*100)>0.6097561]
row_sum.S50.prc_vrc[as.vector(row_sum.S50.prc_vrc/sum(row_sum.S50.prc_vrc)*100)>0.6097561]
                                                                  

pdf("AA_from_BD_x_Empty.@1.pdf")                             
qqnorm(row_sum.S50.prc_vrc, main = "Normal QQ Plot - x")        
qqline(row_sum.S50.prc_vrc, col = "red")
dev.off()
          
x= ((rowSums(bP170_2_1))/sum(bP170_2_1))/((rowSums(bP170_1_1))/sum(bP170_1_1))       
x[is.na(x)==TRUE] = 0
cor(x[P166],row_sum.S50.prc_vrc[P166])  # [1] 0.9320228
plot(x,row_sum.S50.prc_vrc)  
               
pdf("P170.AA_from_BD_x_Empty.@2.pdf")
plot(row_sum.S49.prc_vrc[P166],row_sum.S50.prc_vrc[P166], main = "", xlab = "Number of Reads in No-Selection Media", ylab = "Number of Reads in Selection Media" ) 
dev.off() 
png("P170.AA_from_BD_x_Empty.@2.png")
plot(row_sum.S49.prc_vrc[P166],row_sum.S50.prc_vrc[P166], main = "", xlab = "Number of Reads in No-Selection Media", ylab = "Number of Reads in Selection Media" ) 
dev.off()


# S51 empty x AD
S51.prc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.real.txt",0.5)
S51.prc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.prc_p.txt",0.5)
S51.p_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.p_prc.txt",0.5)
S51.p_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.p_p.txt",0.5)
S51.prc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.prc_v.txt",0.5)
S51.prc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.prc_vrc.txt",0.5)
S51.p_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.p_v.txt",0.5)
S51.p_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.p_vrc.txt",0.5)    
S51.v_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.v_prc.txt",0.5)  
S51.vrc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.vrc_prc.txt",0.5)  
S51.v_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.v_p.txt",0.5)       
S51.vrc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.vrc_p.txt",0.5)   
S51.vrc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.vrc_vrc.txt",0.5)
S51.vrc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.vrc_v.txt",0.5)
S51.v_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.v_vrc.txt",0.5)
S51.v_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/S51.ppi.txt.v_v.txt",0.5)

col_sum.S51.prc_prc = colSums(S51.prc_prc) 
col_sum.S51.prc_p = colSums(S51.prc_p)
col_sum.S51.prc_v = colSums(S51.prc_v)   
col_sum.S51.prc_vrc = colSums(S51.prc_vrc)   

plot(col_sum.S51.prc_prc) # 7337                                            
plot(col_sum.S51.prc_p)   # 257276
plot(col_sum.S51.prc_v)   # 103287
plot(col_sum.S51.prc_vrc) # 945066     
                                    
plot(col_sum.S51.prc_prc/sum(row_sum.S50.prc_prc)*100) # 7337                                            
plot(col_sum.S51.prc_p)   # 257276
plot(col_sum.S51.prc_v)   # 103287
plot(col_sum.S51.prc_vrc) # 945066        

length(col_sum.S51.prc_prc) # 7337                                            
length(col_sum.S51.prc_p)   # 257276
length(col_sum.S51.prc_v)   # 103287
length(col_sum.S51.prc_vrc) # 945066




## without selection


S53.prc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.real.txt",0.5)
S53.prc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.prc_p.txt",0.5)
S53.p_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.p_prc.txt",0.5)
S53.p_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.p_p.txt",0.5)
S53.prc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.prc_v.txt",0.5)
S53.prc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.prc_vrc.txt",0.5)
S53.p_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.p_v.txt",0.5)
S53.p_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.p_vrc.txt",0.5)    
S53.v_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.v_prc.txt",0.5)  
S53.vrc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.vrc_prc.txt",0.5)  
S53.v_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.v_p.txt",0.5)       
S53.vrc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.vrc_p.txt",0.5)   
S53.vrc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.vrc_vrc.txt",0.5)
S53.vrc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.vrc_v.txt",0.5)
S53.v_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.v_vrc.txt",0.5)
S53.v_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.v_v.txt",0.5)


row_sum.S53.prc_prc = rowSums(S53.prc_prc) 
row_sum.S53.prc_p = rowSums(S53.prc_p)
row_sum.S53.prc_v = rowSums(S53.prc_v)   
row_sum.S53.prc_vrc = rowSums(S53.prc_vrc)

plot(row_sum.S53.prc_prc)  
plot(row_sum.S53.prc_p)
plot(row_sum.S53.prc_v)
plot(row_sum.S53.prc_vrc)
                         
PRC	PRC 	3118284
PRC	P 	1126015
P	PRC 	960386
P	P 	980
PRC	V 	652716
PRC	VRC 	120506
P	V 	102059
P	VRC 	155
V	PRC 	413661
VRC	PRC 	85210
V	P 	80244
VRC	P 	166
VRC	VRC 	7659
VRC	V 	251380
V	VRC 	200929
V	V 	991349



## with selection
                                                                                                                    
S54.prc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.real.txt",0.5)
S54.prc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.prc_p.txt",0.5)
S54.p_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.p_prc.txt",0.5)
S54.p_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.p_p.txt",0.5)
S54.prc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.prc_v.txt",0.5)
S54.prc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.prc_vrc.txt",0.5)
S54.p_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.p_v.txt",0.5)
S54.p_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.p_vrc.txt",0.5)    
S54.v_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.v_prc.txt",0.5)  
S54.vrc_prc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.vrc_prc.txt",0.5)  
S54.v_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.v_p.txt",0.5)       
S54.vrc_p = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.vrc_p.txt",0.5)   
S54.vrc_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.vrc_vrc.txt",0.5)
S54.vrc_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.vrc_v.txt",0.5)
S54.v_vrc = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.v_vrc.txt",0.5)
S54.v_v = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.v_v.txt",0.5)


row_sum.S54.prc_prc = rowSums(S54.prc_prc) 
row_sum.S54.prc_p = rowSums(S54.prc_p)
row_sum.S54.prc_v = rowSums(S54.prc_v)   
row_sum.S54.prc_vrc = rowSums(S54.prc_vrc)

plot(row_sum.S54.prc_prc)  
plot(row_sum.S54.prc_p)
plot(row_sum.S54.prc_v)
plot(row_sum.S54.prc_vrc)                   


sum(row_sum.S54.prc_prc)  # 3426672
sum(row_sum.S54.prc_p)    # 1064236
sum(row_sum.S54.prc_v)    # 796500
sum(row_sum.S54.prc_vrc)  # 70666


plot(row_sum.S50.prc_vrc, row_sum.S54.prc_vrc) 

row_sum.S54.prc_prc[as.vector(row_sum.S54.prc_prc/sum(row_sum.S54.prc_prc)*100)>0.6097561]
row_sum.S54.prc_p[as.vector(row_sum.S54.prc_p/sum(row_sum.S54.prc_p)*100)>0.6097561]
row_sum.S54.prc_v[as.vector(row_sum.S54.prc_v/sum(row_sum.S54.prc_v)*100)>0.6097561]
row_sum.S54.prc_vrc[as.vector(row_sum.S54.prc_vrc/sum(row_sum.S54.prc_vrc)*100)>0.6097561]      


length(row_sum.S54.prc_prc) # 7337                                            
length(row_sum.S54.prc_p)   # 257276
length(row_sum.S54.prc_v)   # 103287
length(row_sum.S54.prc_vrc) # 945066
                 


plot(-log(sort(1/1:164),2),-log(sort( row_sum.S50.prc_vrc/sum(row_sum.S50.prc_vrc) ),2) )

                                                                                 
# outlier detection
http://eurekastatistics.com/using-the-median-absolute-deviation-to-find-outliers/

require(BHH2)
#x <- c(1, 2, 3, 3, 4, 4, 4, 5, 5.5, 6, 6, 6.5, 7, 7, 7.5, 8, 9, 12, 52, 90) 
x=log2(row_sum.S54.prc_prc+1)
x=log2(row_sum.S50.prc_vrc+1) 
x=(row_sum.S50.prc_vrc)    
dotPlot(x)
mean.x <- mean(x)
sd.x   <- sd(x)
lines(rep(mean.x         , 2), c(0.2, 0.25)) # Vertical line at mean
lines(rep(mean.x + 2*sd.x, 2), c(0.2, 0.25)) # Vertical line at mean + 2 SD
text(mean.x         , 0.3, expression(bar(x)))
text(mean.x + 2*sd.x, 0.3, expression(paste(bar(x), " + 2s")))

                                                                     
PRC_nde1_p44	0
PRC_kif1b_p120	0
PRC_pard3b_p159	0
PRC_pard3b_p160	0
PRC_pard6g_p163	0
PRC_cdc42_p168	0
PRC_ywhae_p170	0


PRC_cdc42_p168	0	0
PRC_nde1_p44	0	0
PRC_pard3b_p159	0	0
  





#
# 2017-07-03_MiSeq; P170
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
#                                                                                                                                   
m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/60_W.ppi.txt.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                  "sum = 2538094"
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/60_Q.ppi.txt.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in    "sum = 2431908"
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/58_BW.ppi.txt.refine.txt",0.5) # BD E375pLT x AD Empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in     "sum = 144"
m4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/58_BQ.ppi.txt.refine.txt",0.5) # BD E375pLT x AD Empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in     "sum = 244"
m5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/58_AW.ppi.txt.refine.txt",0.5) # BD Empty x AD E375pLT; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in     "sum = 126"
EB_2 = InteractionScores(m1,m2,1.0,4)  
NEB_2 = NewInteractionScores(m1,m2,1.0,4)  

# DATA from blastn-short
bm1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_W.ppi.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "sum = 2629509"
bm2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_Q.ppi.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in        "sum = 2507307"
bm3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/58_BW.ppi.txt",0.5) # BD E375pLT x AD Empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in         "sum = 219"
bm4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/58_BQ.ppi.txt",0.5) # BD E375pLT x AD Empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in         "sum = 306"
bm5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/58_AW.ppi.txt",0.5) # BD Empty x AD E375pLT; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in         "sum = 219"


# Save file
melted_raw_is_EB2 = melt(EB_2)
write.table(melted_raw_is_EB2,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-07-03_MiSeq/melted_EB.IS.txt",sep="\t")
melted_60_W = melt(m1)
write.table(melted_60_W,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-07-03_MiSeq/melted_60_W.count.txt",sep="\t")
melted_60_Q = melt(m2)
write.table(melted_60_Q,file="/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/2017-07-03_MiSeq/melted_60_Q.count.txt",sep="\t")

       



BW_58 = detailed_ppi( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/58_BW.ppi.txt" )    
AW_58 = detailed_ppi( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/58_BQ.ppi.txt" )     ## File is changed with AW
BQ_58 = detailed_ppi( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/58_AW.ppi.txt" )     ## File is changed with BQ

# FIGURE - E375 set to check auto-activator          
plot( BW_58$prc_vrc[,7], BQ_58$prc_vrc[,7], xlab="Read counts from no selection", ylab="Read counts from selection", main="E375 set with Binding domain"  )

x = cbind( BW_58$prc_vrc[,7], BQ_58$prc_vrc[,7] )
write_clip(x)
write_clip(row.names(x))


# BW_58
[1] "prc_prc	144"
[1] "prc_p	18174"
[1] "p_prc	14"
[1] "p_p	6"
[1] "prc_v	48"
[1] "prc_vrc	83069"
[1] "p_v	0"
[1] "p_vrc	31"
[1] "v_prc	3"
[1] "vrc_prc	4"
[1] "v_p	5"
[1] "vrc_p	1"
[1] "vrc_vrc	1132"
[1] "vrc_v	4265"
[1] "v_vrc	18790"
[1] "v_v	533"

# BQ_58  ==> seems like AW_58
[1] "prc_prc	244"
[1] "prc_p	270"
[1] "p_prc	15517"
[1] "p_p	5"
[1] "prc_v	20"
[1] "prc_vrc	54"
[1] "p_v	4"
[1] "p_vrc	3"
[1] "v_prc	121"
[1] "vrc_prc	73411"
[1] "v_p	1"
[1] "vrc_p	85"
[1] "vrc_vrc	1108"
[1] "vrc_v	21441"
[1] "v_vrc	5555"
[1] "v_v	1532"
        
# AW_58  ==> seems like BQ_58
[1] "prc_prc	126"
[1] "prc_p	24969"
[1] "p_prc	28"
[1] "p_p	19"
[1] "prc_v	125"
[1] "prc_vrc	91306"
[1] "p_v	0"
[1] "p_vrc	101"
[1] "v_prc	13"
[1] "vrc_prc	5"
[1] "v_p	13"
[1] "vrc_p	5"
[1] "vrc_vrc	1233"
[1] "vrc_v	5253"
[1] "v_vrc	22641"
[1] "v_v	1549"
  
##
## Plot for auto-activator  (FIGURE!!)
##     
plot(BW_58$prc_vrc[,7], BQ_58$prc_vrc[,7],xlab="- Selection",ylab="+ Selection") 
#names(BW_58$prc_vrc[,7])
index_list = which(BQ_58$prc_vrc[,7]>500)        
#BQ_58$prc_vrc[index_list,7]    
#BW_58$prc_vrc[index_list,7]       
#names(BW_58$prc_vrc[index_list,7]) 
text(BW_58$prc_vrc[index_list,7]+50, BQ_58$prc_vrc[index_list,7]+500, names(BW_58$prc_vrc[index_list,7]),offset=10000, cex=0.5 )
length(index_list)
                  
##
## Plot for toxicity     (FIGURE!!)
##     
#plot(BW_58$prc_vrc[,7], t(BQ_58$vrc_prc)[,7],xlab="Binding domain",ylab="Activation domain")   
plot(BW_58$prc_vrc[,7], BQ_58$vrc_prc[7,],xlab="Binding domain",ylab="Activation domain")       
# plot( rowSums( BQ_58$p_prc ) )
# with AD
index_list = which( BQ_58$vrc_prc[7,] == 0 )
length( index_list )
BQ_58$vrc_prc[7,index_list]                      
# with BD   
index_list = which( BW_58$prc_vrc[7,] == 0 )
length( index_list )
BW_58$vrc_prc[7,index_list]
 


#
# 2017-08-15_MiSeq; P170
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
#                                                                                                                                   
m1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/S61_SWD.ppi.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                 # "sum = 2538094"
m2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/S61_PWD.ppi.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in   # "sum = 2431908"
m3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/S61_PQD.ppi.refine.txt",0.5) # BD E375pLT x AD Empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in     #"sum = 144"
EB_3 = InteractionScores(m2,m3,1.0,4)  
NEB_3 = NewInteractionScores(m2,m3,1.0,4)   

# DATA from blastn-short
bm1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_SWD.ppi.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "# ==> 2909222 out of 6012622 (48.4%)
bm2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PWD.ppi.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in        "# ==> 1956131 out of 4473498 (43.7%)
bm3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt",0.5) # BD E375pLT x AD Empty; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in          "# ==> 1684777 out of 3553019 (47.4%)

EB_3_1 = InteractionScores(bm1,bm3,1.0,4)
EB_3_2 = InteractionScores(bm2,bm3,1.0,4)     
CorMatrix(EB_3_1,EB_3_2)	# 0.936

NEB_3_1 = NewInteractionScores(bm1,bm3,1.0,4)     
NEB_3_2 = NewInteractionScores(bm2,bm3,1.0,4)     
CorMatrix(NEB_3_1,NEB_3_2)	#0.9816  --> 0.975



#
# 2017-08-22_MiSeq; Roth with Seaprep
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
m1_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt",0.5) # Roth; -W; R75_01 ~ R75_07 spike in                  "sum = 2538094"
m2_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt",0.5) # Roth; -W/-H/-A/Aba 1/4; R75_01 ~ R75_07 spike in    "sum = 2431908"
m3_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt",0.5) # Roth; -W/-H/-A/Aba 1/8; R75_01 ~ R75_07 spike in     "sum = 144"
m4_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt",0.5) # Roth; -W/-H/-A/Aba ?; R75_01 ~ R75_07 spike in     "sum = 144"
is2_4 = InteractionScores(m1_4,m2_4,1.0)
is3_4 = InteractionScores(m1_4,m3_4,1.0)
is4_4 = InteractionScores(m1_4,m4_4,1.0)
CorMatrix(is2_4, is3_4)		#[1] 0.6794
CorMatrix(is2_4, is4_4) 	#[1] 0.6523
CorMatrix(is3_4, is4_4)     #[1] 0.5395
                            
nis2_4 = NewInteractionScores(m1_4,m2_4,1.0)
nis3_4 = NewInteractionScores(m1_4,m3_4,1.0)
nis4_4 = NewInteractionScores(m1_4,m4_4,1.0)
CorMatrix(nis2_4, nis3_4)	#[1] 0.9227 	--> 0.9286
CorMatrix(nis2_4, nis4_4) 	#[1] 0.7908     --> 0.7935
CorMatrix(nis3_4, nis4_4)   #[1] 0.7402     --> 0.7412



#############################################################################################################################
# Parameter Scan Result - 2017-07-04

file1 = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt"
file2 = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt"

ParameterScanLogP( file1, file2, PPI_Data = HIPPIE, len_start=10, len_end=32, len_step=1, count_start=1, count_end=11, count_step=1 )
#ParameterScanLogP <-function( file1 = "", file2 = "", len_start=10, len_end=30, len_step=2, count_start=3, count_end=10, count_step=1, PPI_Data = HIPPIE)     

#
# ParameterScanLogP results for HIPPIE 
# file1 = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt"
# file2 = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt"
# PPI_Data = HIPPIE 
# [wrong]
#	10	12	14	16	18	20	22	24	26	28
#3	12.119535	12.61457	12.995637	13.568711	13.555425	13.173335	13.172121	13.068685	13.057941	13.160244
#4	11.962086	12.539481	12.915856	13.438762	13.466109	12.97254	13.099795	12.921178	12.930902	13.094324
#5	11.816555	12.567495	12.746132	13.315637	12.9815	12.898217	12.924346	12.843103	12.823358	12.533881
#6	11.867978	12.480877	12.583509	13.228927	12.767213	12.836097	12.779727	12.721758	12.741005	12.424757
#7	11.773232	12.255916	12.158692	12.734643	12.699663	12.646758	12.695861	12.610918	12.178405	12.345744
#8	11.697342	12.072075	12.050032	12.595044	12.63795	12.537751	12.574609	12.512774	12.079857	12.268284
#9	11.506608	11.946044	11.834187	12.545103	12.449787	12.430694	12.484846	11.96775	12.023242	12.127461  
#
# [new one]
#	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31
#1	10.802748	12.119876	13.059394	12.516753	12.679873	13.174251	13.358734	13.388644	13.57586	13.611758	13.556464	13.553254	13.24803	13.225925	13.176139	13.17876	13.064775	13.053777	13.050122	13.142196	13.168991	13.242222
#2	12.028187	12.903751	12.701286	12.645538	12.838988	13.371585	13.386635	13.57257	13.610713	13.555993	13.564166	13.16492	13.228682	13.174743	13.177468	13.064745	13.053391	13.048987	13.141561	13.165478	13.239575	12.671434
#3	12.119535	12.663326	12.61457	12.812952	12.995637	13.389265	13.568711	13.601736	13.555425	13.563776	13.173335	13.146376	13.172121	13.175241	13.068685	13.053229	13.057941	13.135784	13.160244	13.226631	12.67565	12.677124
#4	12.072279	12.469607	12.739373	12.965911	13.042091	13.556264	13.592613	13.558071	13.564456	13.168683	13.155822	13.087448	13.171926	13.066824	13.051383	13.061278	13.13338	13.155135	13.216382	12.679304	12.673967	12.645071
#5	11.970603	12.728077	12.893139	13.01701	13.204044	13.576312	13.549096	13.555541	13.173143	13.151897	13.095193	13.089115	13.061447	13.056215	13.056906	13.139007	13.148219	13.205375	12.665492	12.679329	12.644034	12.592768
#6	12.062684	12.770879	12.929889	13.173612	13.209984	13.536477	13.55951	13.165151	13.153011	13.092724	13.098619	12.984299	13.05286	13.054162	13.140107	13.142862	13.201517	12.654203	12.681026	12.640963	12.590432	12.487664
#7	11.989712	12.795174	13.098632	13.164469	12.900135	13.547202	13.161114	13.150605	13.088555	13.099987	12.993479	12.967755	13.056882	13.133865	13.148103	13.195974	12.641698	12.665015	12.642248	12.588336	12.506169	12.656675
#8	11.917396	12.952679	13.086192	12.850607	12.918858	13.157995	13.147218	13.090899	13.097322	12.992114	12.975215	12.96538	13.129625	13.153036	13.187509	12.637184	12.655507	12.643698	12.587812	12.501444	12.654235	12.648473
#9	11.885813	12.983968	13.050536	12.857467	12.78083	13.143426	13.098946	13.088311	12.993779	12.972183	12.966887	13.038405	13.145921	13.192334	12.631345	12.643534	12.627494	12.586145	12.49991	12.675291	12.644911	12.594152
#10	11.801935	12.949871	12.758453	12.729774	12.783486	13.095665	13.098971	12.99578	12.964226	12.968021	13.040553	13.060533	13.198094	12.622922	12.638572	12.617795	12.586961	12.496354	12.670401	12.642901	12.599481	12.751867


ParameterScanLogP( file1, file2, PPI_Data = BioGrid )

#
# ParameterScanLogP results for BioGrid 
# file1 = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt"
# file2 = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt"
# PPI_Data = BioGrid
#
#	10	12	14	16	18	20	22	24	26	28
#3	14.495799	15.029193	15.347186	15.948191	15.956298	15.564222	15.537226	15.408823	15.439875	15.391427
#4	14.420761	14.923497	15.285795	15.88855	15.880093	15.336893	15.483688	15.334507	15.315759	15.391104
#5	14.238185	14.972133	15.123764	15.741253	15.392489	15.300403	15.298795	15.285252	15.220273	14.800905
#6	14.267287	14.90115	15.020496	15.675457	15.149908	15.258833	15.225915	15.167845	15.212719	14.723557
#7	14.191203	14.673798	14.541356	15.173834	15.115809	15.045039	15.170441	15.074197	14.61398	14.680292
#8	14.11844	14.532967	14.447463	15.011128	15.067272	15.021188	15.052324	15.053545	14.55032	14.644636
#9	13.922697	14.380653	14.229211	14.99126	14.856728	14.942121	14.983127	14.462983	14.524834	14.474399
#





#############################################################################################################################
# Comparison with Bowtie mapping method (Y2H_Friedrich.sh) - 2017.07.11
#
# Friedrich method was done by Y2H_Friedrich.sh file with manual modification
# After running script, run SAM.py file to get ppi matrix as following
#                                                        
# python SAM.py ../data/roth2016_control_set_plus_control.fa output/2016-12-22_MiSeq/Friedrich/17543_S1.sam >  output/2016-12-22_MiSeq/Friedrich/17543_S1.ppi.txt
# python SAM.py ../data/roth2016_control_set_plus_control.fa output/2016-12-22_MiSeq/Friedrich/17544_S2.sam >  output/2016-12-22_MiSeq/Friedrich/17544_S2.ppi.txt                     
#
sm1 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Friedrich/17543_S1.ppi.txt", 0.5 )      # ==> it gives diagonal signal (un-reliable) and less number of PPI map
js1 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt.3.10.refine.txt", 0.5 )

sm2 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Friedrich/17544_S2.ppi.txt", 0.5 )       #==> it gives diagonal signal (un-reliable) and less number of PPI map
js2 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt.3.10.refine.txt", 0.5 )

sm_i = InteractionScores(sm1,sm2,1.0)
js_i = InteractionScores(js1,js2,1.0)

PPI_PerformanceCheck( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Friedrich/17543_S1.ppi.txt", "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Friedrich/17544_S2.ppi.txt", HIPPIE )
# "8.064136-1345436-1141479"
PPI_PerformanceCheck( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt.3.10.refine.txt", "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt.3.10.refine.txt", HIPPIE )
# "12.119535-1666813-1472273"
PPI_PerformanceCheck( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt.3.16.refine.txt", "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt.3.16.refine.txt", HIPPIE )
# "13.568711-1911898-1505797"


sm3 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Friedrich/S53.ppi.txt", 0.5, 4 )       ==> it gives diagonal signal (un-reliable) and less number of PPI map
js1 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt.3.10.refine.txt", 0.5 )

sm4 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Friedrich/S54.ppi.txt", 0.5, 4 )       ==> it gives diagonal signal (un-reliable) and less number of PPI map
js2 = drawMatrix( "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt.3.10.refine.txt", 0.5 )

sm3_i = InteractionScores(sm3,sm4,1.0)
js_i = InteractionScores(js1,js2,1.0)





                            
                           
########
# check uniform dist in column
   
x = ToxicityTest(m1_3,1,78)

y=uniform.test( hist(m4[,]) )
x1 = c()
x2 = c()
for (i in 1:78)
{
	print(i)
	y=uniform.test( hist(m1_3[,i]) )
	x1[i] = y$statistic
	x2[i] = sum(m1_3[,i])
}
plot(x1,x2,xlab="Uniform-test statistics",ylab="PPI read count")


	y=uniform.test( hist(m4[,]) )

x1 = c()
x2 = c()
for (i in 1:173)
{        
	zeros = 1:173 * 0
	print(i)             
	#if (all.equal(m1[,i]) && m1[1,i]==0){  
	if ( sum( as.vector( m1[,i] != zeros ) ) == 0 ){  
		x1[i] = 600.0
	}else{
		y=uniform.test( hist(m1[,i]) )
		x1[i] = y$statistic
	}
	x2[i] = sum(m1[,i])
}
plot(x1,x2,xlab="Uniform-test statistics",ylab="PPI read count")
                                                                           


















































#===========================================================================            
# For Figures
#===========================================================================            

### Load Data
# from blastn-short                                                          

#setwd("/Users/jyang/Dropbox_CRG/Collaboration/mireia/Paper/v1/Figures")         
setwd("/Users/jyang/Dropbox_CRG/matrixppi/figures")   

#                  
# 2016-12-22_MiSeq; Roth75-exp1; MGj46  
#
# These are never cloned into pENTR223: 
# R75_41,WDR7,WD repeat domain 7
# R75_72,PLEKHG7,pleckstrin homology domain containing, family G (with RhoGef domain) member 7 

m1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17543.ppi.txt.3.15.refine.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1947500"
m1_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17544.ppi.txt.3.15.refine.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1533027" 
m1_RW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17545.ppi.txt.3.15.refine.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2121357" 
m1_RA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/17546.ppi.txt.3.15.refine.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1819108"

bm1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
bm1_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
bm1_RW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
bm1_RA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"       

jb1_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2016-12-22.S1.ppi.txt", 0.5)										# no selection; sum = 1968532  == bm1_1
jb1_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2016-12-22.S2.ppi.txt", 0.5)										# selection; sum = 1541502     == bm2_1
se1_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_17543_S1_sorted_flag113tidyformat.xlsx.txt",0.5)	# no selection; sum = 1167422
se1_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_17544_S2_sorted_flag113tidyformat.xlsx.txt",0.5)	# selection; sum = 899571
sl1_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_17543_S1_sorted_flag113tidyformat.xlsx.txt",0.5)	# no selection; sum = 1487814
sl1_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_17544_S2_sorted_flag113tidyformat.xlsx.txt",0.5)	# selection; sum = 1183181

is_j_1 = InteractionScores(m1_1,m2_1,1.0,4) 
is_jb_1 = InteractionScores(jb1_1,jb2_1,1.0,4) 
is_se_1 = InteractionScores(se1_1,se2_1,1.0,4) 
is_sl_1 = InteractionScores(sl1_1,sl2_1,1.0,4) 

nis_j_1 = NewInteractionScores(m1_1,m2_1,1.0,4)  
nis_jb_1 = NewInteractionScores(jb1_1,jb2_1,1.0,4) 
nis_se_1 = NewInteractionScores(se1_1,se2_1,1.0,4) 
nis_sl_1 = NewInteractionScores(sl1_1,sl2_1,1.0,4) 

bisR75_1 = InteractionScores(bm1_1,bm2_1,1.0,4) 
bnisR75_1 = NewInteractionScores(bm1_1,bm2_1,1.0,4) 
bisR75_1b = InteractionScores(bm3_1,bm4_1,1.0,4) 
bnisR75_1b = NewInteractionScores(bm3_1,bm4_1,1.0,4)
                                                       
simple_plotLog2Matrix( jb1_W, 0.5, 7, filename="./jb1_W.pdf" )
simple_plotLog2Matrix( se1_W, 0.5, 7, filename="./se1_W.pdf" )
simple_plotLog2Matrix( sl1_W, 0.5, 7, filename="./sl1_W.pdf" )

simple_plotLog2Matrix( jb1_A, 0.5, 7, filename="./jb1_A.pdf" )
simple_plotLog2Matrix( se1_A, 0.5, 7, filename="./se1_A.pdf" )
simple_plotLog2Matrix( sl1_A, 0.5, 7, filename="./sl1_A.pdf" )




#
# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
# R75_11,AP2B1,adaptor-related protein complex 2, beta 1 subunit
# R75_30,CTBP1,C-terminal binding protein 1
# R75_52,RBBP8,retinoblastoma binding protein 8
# R75_58,FHL3,four and a half LIM domains 3
# R75_70,TEX11,testis expressed 11
#
# These are never cloned into pENTR223: 
# R75_41,WDR7,WD repeat domain 7
# R75_72,PLEKHG7,pleckstrin homology domain containing, family G (with RhoGef domain) member 7       

m2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/S1.ppi.txt.3.15.refine.txt",0.5) # no selection; R75_MGj51  "sum = 919439"
m2_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/S2.ppi.txt.3.15.refine.txt",0.5) # selection; ACCRD         "sum = 1956718"
m2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/S3.ppi.txt.3.15.refine.txt",0.5) # selection; QCCRD         "sum = 791387"

bm2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; R75_MGj51  "sum = 960245"
bm2_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; ACCRD         "sum = 2024754"
bm2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK

bisR75_2 = InteractionScores(bm1_2,bm2_2,1.0,4) 
bnisR75_2 = NewInteractionScores(bm1_2,bm2_2,1.0,4)
bisR75_2b = InteractionScores(bm1_2,bm3_2,1.0,4) 
bnisR75_2b = NewInteractionScores(bm1_2,bm3_2,1.0,4)

jb2_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//jae-seong/blast/X75/2017-02-22.S1.ppi.txt", 0.5, 7)										# no selection; sum = 1968532
jb2_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//jae-seong/blast/X75/2017-02-22.S2.ppi.txt", 0.5, 7)										# selection; sum = 1541502
jb2_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//jae-seong/blast/X75/2017-02-22.S3.ppi.txt", 0.5, 7)										# selection; sum = 1541502
se2_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//sebastian_old/endon/X75/Ref_roth75_seq_51-WCCRD_S1_sorted_flag113tidyformat.xlsx.txt",0.5,7)
se2_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//sebastian_old/endon/X75/Ref_roth75_seq_51-A1CCRD_S2_sorted_flag113tidyformat.xlsx.txt",0.5,7)
se2_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//sebastian_old/endon/X75/Ref_roth75_seq_51-QCCRD_S3_sorted_flag113tidyformat.xlsx.txt",0.5,7)
sl2_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//sebastian_old/local/X75/Ref_roth75_seq_51-WCCRD_S1_sorted_flag113tidyformat.xlsx.txt",0.5,7)
sl2_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//sebastian_old/local/X75/Ref_roth75_seq_51-A1CCRD_S2_sorted_flag113tidyformat.xlsx.txt",0.5,7)
sl2_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data//sebastian_old/local/X75/Ref_roth75_seq_51-QCCRD_S3_sorted_flag113tidyformat.xlsx.txt",0.5,7)

simple_plotLog2Matrix( jb2_W, 0.5, 7, filename="./jb2_W.pdf" )
simple_plotLog2Matrix( se2_W, 0.5, 7, filename="./se2_W.pdf" )
simple_plotLog2Matrix( sl2_W, 0.5, 7, filename="./sl2_W.pdf" )        

simple_plotLog2Matrix( jb2_A, 0.5, 7, filename="./jb2_A.pdf" )
simple_plotLog2Matrix( se2_A, 0.5, 7, filename="./se2_A.pdf" )
simple_plotLog2Matrix( sl2_A, 0.5, 7, filename="./sl2_A.pdf" )
                                                                  
simple_plotLog2Matrix( jb2_Q, 0.5, 7, filename="./jb2_Q.pdf" )
simple_plotLog2Matrix( se2_Q, 0.5, 7, filename="./se2_Q.pdf" )
simple_plotLog2Matrix( sl2_Q, 0.5, 7, filename="./sl2_Q.pdf" )

is_j_2 = InteractionScores(m1_2,m2_2,1.0,4) 
is_jb_2 = InteractionScores(jb1_2,jb2_2,1.0,4) 
is_se_2 = InteractionScores(se1_2,se2_2,1.0,4) 
is_sl_2 = InteractionScores(sl1_2,sl2_2,1.0,4) 
nis_j_2 = NewInteractionScores(m1_2,m2_2,1.0,4) 
nis_jb_2 = NewInteractionScores(jb1_2,jb2_2,1.0,4) 
nis_se_2 = NewInteractionScores(se1_2,se2_2,1.0,4) 
nis_sl_2 = NewInteractionScores(sl1_2,sl2_2,1.0,4)
        
#
# 2017-03-03_MiSeq; Roth75-exp3; [53]
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
# R75_11,AP2B1,adaptor-related protein complex 2, beta 1 subunit
# R75_30,CTBP1,C-terminal binding protein 1
# R75_52,RBBP8,retinoblastoma binding protein 8
# R75_58,FHL3,four and a half LIM domains 3
# R75_70,TEX11,testis expressed 11
                                                             
m3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S1.ppi.txt.3.15.refine.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = 1369865"
m3_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S2.ppi.txt.3.15.refine.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = 947672"
m3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S3.ppi.txt.3.15.refine.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1218835"
m3_QL = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/S7.ppi.txt.3.15.refine.txt",0.5) # selection; 53_QL4RD	18715	 "sum = 715800"     ## NOT USE (LOW QUALITY) 

bm3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
bm3_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
bm3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
bm3_QL = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      

jb3_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-03-03.S1.ppi.txt", 0.5)										# no selection; sum = 1968532
jb3_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-03-03.S2.ppi.txt", 0.5)										# selection; sum = 1541502
jb3_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-03-03.S3.ppi.txt", 0.5)										# selection; sum = 1541502
se3_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_18689_S1_sorted_flag113tidyformat.xlsx.txt",0.5)
se3_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_18690_S2_sorted_flag113tidyformat.xlsx.txt",0.5)
se3_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_18691_S3_sorted_flag113tidyformat.xlsx.txt",0.5)
sl3_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_18689_S1_sorted_flag113tidyformat.xlsx.txt",0.5)
sl3_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_18690_S2_sorted_flag113tidyformat.xlsx.txt",0.5)
sl3_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_18691_S3_sorted_flag113tidyformat.xlsx.txt",0.5)     

sa3_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18689_S1_sorted_flag113tidyformat.xlsx.txt",0.5) 
sa3_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18690_S2_sorted_flag113tidyformat.xlsx.txt",0.5) 
sa3_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18691_S3_sorted_flag113tidyformat.xlsx.txt",0.5)  
is_sa3_A = InteractionScores(sa3_W,sa3_A,1.0,4) 
is_sa3_Q = InteractionScores(sa3_W,sa3_Q,1.0,4)          

 


is_jb_3 = InteractionScores(jb3_W,jb3_Q,1.0,4) 
is_se_3 = InteractionScores(se3_W,se3_Q,1.0,4) 
is_sl_3 = InteractionScores(sl3_W,sl3_Q,1.0,4) 
nis_jb_3 = NewInteractionScores(jb3_W,jb3_Q,1.0,4) 
nis_se_3 = NewInteractionScores(se3_W,se3_Q,1.0,4) 
nis_sl_3 = NewInteractionScores(sl3_W,sl3_Q,1.0,4)


bisR75_3 = InteractionScores(bm3_W,bm3_A,1.0,4) 
bnisR75_3 = NewInteractionScores(bm3_W,bm3_A,1.0,4)
bisR75_3b = InteractionScores(bm3_W,m3_Q,1.0,4) 
bnisR75_3b = NewInteractionScores(bm3_W,m3_Q,1.0,4)  
bisR75_3c = InteractionScores(bm3_W,m3_QL,1.0,4)      
bnisR75_3c = NewInteractionScores(bm3_W,m3_QL,1.0,4)   

#
# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61]
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator       
#
# 250ml seaprep 0.5% 
m1_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/S64_SWD.ppi.txt.3.15.refine.txt",0.5) # Roth; -W;                "sum = "
m2_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/S64_SA4D.ppi.txt.3.15.refine.txt",0.5) # Roth; -W/Aba 1/4;    "sum = "
m3_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/S64_SA8D.ppi.txt.3.15.refine.txt",0.5) # Roth; -W/Aba 1/8;    "sum = "
m4_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/S64_SQD.ppi.txt.3.15.refine.txt",0.5) # Roth; -W/-H/-A;      "sum = "

bm4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt",0.5) # Roth; -W;                 "sum = 1491017"
bm4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
bm4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
bm4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt",0.5) # Roth; -W/-H/-A;      "sum = 618177"

bm4_SW_new = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt.new",0.5) # Roth; -W;                 "sum = 1450453"
bm4_S4A_new = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt.new",0.5) # Roth; -W/Aba 1/4;   "sum = 1145754"
bisR75_4 = InteractionScores(bm4_SW,bm4_S4A,1.0,4) 
bnisR75_4 = NewInteractionScores(bm4_SW,bm4_S4A,1.0,4)
bisR75_4b = InteractionScores(bm4_SW,bm4_S8A,1.0,4) 
bnisR75_4b = NewInteractionScores(bm4_SW,bm4_S8A,1.0,4)
bisR75_4c = InteractionScores(bm4_SW,bm4_SQ,1.0,4) 
bnisR75_4c = NewInteractionScores(bm4_SW,bm4_SQ,1.0,4)          

jb4_SW = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-08-22.S64_SWD.ppi.txt", 0.5, 7)										# no selection; sum = 1968532
jb4_S8A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-08-22.S64_SA8D.ppi.txt", 0.5, 7)										# selection; sum = 1541502
jb4_S4A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-08-22.S64_SA4D.ppi.txt", 0.5, 7)										# selection; sum = 1541502
jb4_SQ = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/jae-seong/blast/X75/2017-08-22.S64_SQD.ppi.txt", 0.5, 7)										# selection; sum = 1541502
jb1_4 = jb4_SW; jb2_4 = jb4_S8A; jb3_4 = jb4_S4A; jb4_4 = jb4_SQ;                          

se4_SW = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_21387-64-SWD_S1_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
se4_S8A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_21389-64-SA8D_S3_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
se4_S4A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_21388-64-SA4D_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
se4_SQ = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/endon/X75/Ref_roth75_seq_21390-64-SQD_S4_sorted_flag113tidyformat.xlsx.txt",0.5, 7)

sl4_SW = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_21387-64-SWD_S1_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
sl4_S8A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_21389-64-SA8D_S3_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
sl4_S4A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_21388-64-SA4D_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
sl4_SQ = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X75/Ref_roth75_seq_21390-64-SQD_S4_sorted_flag113tidyformat.xlsx.txt",0.5, 7)
                                                 
sa4_SW = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21387-64-SWD_S1_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
sa4_S8A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21389-64-SA8D_S3_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
sa4_S4A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21388-64-SA4D_S2_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
sa4_SQ = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21390-64-SQD_S4_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
is_sa4_S8A = InteractionScores(sa4_SW,sa4_S8A,1.0,4)    
is_sa4_S4A = InteractionScores(sa4_SW,sa4_S4A,1.0,4)    
is_sa4_SQ = InteractionScores(sa4_SW,sa4_SQ,1.0,4)    

simple_plotLog2Matrix( jb4_SW, 0.5, 7, filename="./figures/jb4_SW.pdf" )
simple_plotLog2Matrix( se4_SW, 0.5, 7, filename="./figures/se4_SW.pdf" )
simple_plotLog2Matrix( sl4_SW, 0.5, 7, filename="./figures/sl4_SW.pdf" )
                       
is_b_4 = InteractionScores(m1_4,m2_4,1.0,4) 
is_jb_4 = InteractionScores(jb4_SW,jb4_S4A,1.0,4) 
is_se_4 = InteractionScores(se4_SW, se4_S4A,1.0,4) 
is_sl_4 = InteractionScores(sl4_SW,sl4_S4A,1.0,4) 

nis_b_4 = NewInteractionScores(m1_4,jb2_4,1.0,4) 
nis_jb_4 = NewInteractionScores(jb4_SW,jb4_S4A,1.0,4) 
nis_se_4 = NewInteractionScores(se4_SW,se4_S4A,1.0,4) 
nis_sl_4 = NewInteractionScores(sl4_SW,sl4_S4A,1.0,4)


                        
jb1_W = BasicStats( jb1_1 )
jb2_W = BasicStats( jb1_2 )
jb3_W = BasicStats( jb1_3 )
jb4_SW = BasicStats( jb4_SW )
se1_W = BasicStats( se1_1 )
se2_W = BasicStats( se1_2 )
se3_W = BasicStats( se1_3 )
se4_SW = BasicStats( se4_SW )
sl1_W = BasicStats( sl1_1 )
sl2_W = BasicStats( sl1_2 )
sl3_W = BasicStats( sl1_3 )
sl4_SW = BasicStats( sl4_SW )

        
#
# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;  []
#
# Removed these genes in the BD-vector
# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
                                                     
bm5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
jb5_W = bm5_SW           



#================================================================================================================
# P170  - EB-Par interacting proteins
#================================================================================================================           

# 4 positive controls from R75 set 
# 47, 26 'dmrtb1','ctbp1' 
# 47, 142 'dmrtb1','rbfox2' 
# 20, 20 'cep55','cep55' 
# 120, 93 'p53','larget' 
# 4 negative controls from R75 set 
# 47, 149 'dmrtb1','reep2' 
# 163, 28 'topbp1','dapk1' 
# 92, 93 'lam','larget' 
# 61, 16 'gaa','ccl14'          

# posC={'dmrtb1','ctbp1';'dmrtb1','rbfox2';'cep55','cep55';'p53','larget'}
# negC={'dmrtb1','reep2';'topbp1','dapk1';'lam','larget';'gaa','ccl14'}

NULL_MATRIX = matrix(0, nrow=nrow(bnisP170_2_1),ncol=ncol(bnisP170_2_1) ) 
rownames(NULL_MATRIX) = rownames(bnisP170_2_1)                       
colnames(NULL_MATRIX) = colnames(bnisP170_2_1)                       
                  

CALIBRATION_POS_PAIRS = matrix(1:16, nrow=8,ncol=2 )
CALIBRATION_POS_PAIRS[1,] = c(47,26)
CALIBRATION_POS_PAIRS[2,] = c(47,142)     
CALIBRATION_POS_PAIRS[3,] = c(20,20)    
CALIBRATION_POS_PAIRS[4,] = c(120,93)
CALIBRATION_POS_PAIRS[5,] = c(26,47)
CALIBRATION_POS_PAIRS[6,] = c(142,47)     
CALIBRATION_POS_PAIRS[7,] = c(20,20)    
CALIBRATION_POS_PAIRS[8,] = c(93,120)

CALIBRATION_NEG_PAIRS = matrix(1:16, nrow=8,ncol=2 )
CALIBRATION_NEG_PAIRS[1,] = c(47,149)
CALIBRATION_NEG_PAIRS[2,] = c(163,28)     
CALIBRATION_NEG_PAIRS[3,] = c(92,93)    
CALIBRATION_NEG_PAIRS[4,] = c(61,16)
CALIBRATION_NEG_PAIRS[5,] = c(149,47)
CALIBRATION_NEG_PAIRS[6,] = c(28,163)     
CALIBRATION_NEG_PAIRS[7,] = c(93,92)    
CALIBRATION_NEG_PAIRS[8,] = c(16,61)
                                                                                 
NULL_MATRIX[CALIBRATION_POS_PAIRS] = 100
NULL_MATRIX[CALIBRATION_NEG_PAIRS] = 50
simple_plotMatrix( NULL_MATRIX[CALIBRATION,CALIBRATION],0,filename= "P170_Spikded_set.pdf");

#simple_plotMatrix( cbind( bnisP170_2_1[CALIBRATION_POS_PAIRS], bnisP170_2_1[CALIBRATION_NEG_PAIRS] ),0 )  
#simple_plotMatrix( cbind( bnisP170_2_2[CALIBRATION_POS_PAIRS], bnisP170_2_2[CALIBRATION_NEG_PAIRS] ),0 )  
#simple_plotMatrix( cbind( bnisP170_2_3[CALIBRATION_POS_PAIRS], bnisP170_2_3[CALIBRATION_NEG_PAIRS] ),0 )  
#simple_plotMatrix( cbind( bnisP170_3_3[CALIBRATION_POS_PAIRS], bnisP170_3_3[CALIBRATION_NEG_PAIRS] ),0 )  
t.test( bnisP170_2_1[CALIBRATION_POS_PAIRS], bnisP170_2_1[CALIBRATION_NEG_PAIRS]  )  
t.test( bnisP170_2_2[CALIBRATION_POS_PAIRS], bnisP170_2_2[CALIBRATION_NEG_PAIRS]  )  
t.test( bnisP170_3_3[CALIBRATION_POS_PAIRS], bnisP170_3_3[CALIBRATION_NEG_PAIRS]  )  


P166 = 8:173    
CALIBRATION = c(16,20,26,28,47,61,92,93,120,142,149,163)        

# 16 "ccl14_p85"  
# 20 "cep55_p86"  
# 26 "ctbp1_p87"  
# 28 "dapk1_p88"  
# 47 "dmrtb1_p89" 
# 61 "gaa_p90"    
# 92 "lam_p91"    
# 93 "larget_p92"
# 120 "p53_p96"    
# 142 "rbfox2_p97" 
# 149 "reep2_p98"  
# 163 "topbp1_p99"
      

# [P170 - Exp1]                      
# ALL PROTEIN INCLUDED 
# R1 - R7 spike-in
# 2017-06-12_MiSeq; P170-exp1    ## 3250858 + 3553164 = 6804022       
#
P170_1_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S53.ppi.txt.3.15.refine.txt",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 
P170_2_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/S54.ppi.txt.3.15.refine.txt",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 

bP170_1_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S53.ppi.txt",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3250858
bP170_2_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S54.ppi.txt",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3553164

slP170_1_1 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X163/Ref_X163nostp_seq_20053-56-4W_S1_sorted_flag113tidyformat.xlsx.txt",0.5, 4)
slP170_2_1  = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian_old/local/X163/Ref_X163nostp_seq_20054-56-4Q_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 4)

bisP170_2_1 = InteractionScores(bP170_1_1,bP170_2_1,1.0,4) 
bnisP170_2_1 = NewInteractionScores(bP170_1_1,bP170_2_1,1.0,4)                              

sl_isP170_2_1 = InteractionScores(slP170_1_1,slP170_2_1,1.0,4) 

BasicStats(P170_1_1[P166,P166])
BasicStats(P170_2_1[P166,P166])
BasicStats(bP170_1_1[P166,P166])
BasicStats(bP170_2_1[P166,P166])
                               

# [P170 - Exp2]                      
# REMOVED AUTO-ACTIVATOR: TACC3_p108, FEZ1_p143, FEZ2_p144
# 2017-07-03_MiSeq; P170-exp2     ## 2629509 + 2507307 = 5136816   
# B: R1,2,3,5 spike-in
# A: R1,2,3,4 spike-in
# 60_W = 8x square plate  -W (No selection)
# 60_Q = 8x square plate  -W/-H/-A/Aba 1/2; (Strong selection)          
# 
P170_1_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/60_W.ppi.txt.3.15.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "sum = 2629509"
P170_2_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/60_Q.ppi.txt.3.15.refine.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2 ? ; R75_01 ~ R75_07 spike in        "sum = 2507307"

bP170_1_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_W.ppi.txt",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "sum = 2629509"
bP170_2_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_Q.ppi.txt",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in        "sum = 2507307"
bisP170_2_2 = InteractionScores(bP170_1_2,bP170_2_2,1.0,4) 
bnisP170_2_2 = NewInteractionScores(bP170_1_2,bP170_2_2,1.0,4)

slP170_1_2 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_20959-60-W_S1_sorted_flag113tidyformat.xlsx.txt",0.5, 4)
slP170_2_2  = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_20960-60-Q_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 4)
sl_isP170_2_2 = InteractionScores(slP170_1_2,slP170_2_2,1.0,4)   
simple_plotMatrix(sl_isP170_2_2[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_2_2[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_2_2[60:120,60:120],0,7)

   
# [P170 - Exp3]      
# 2017-08-15_MiSeq; P170-exp3; [61]     ## 2909222 + 1956131 + 1684777 = 6550130   
#
# 61_SWD = 2x250 ml seaprep 0.5%                    -W (No selection)
# 61_PWD = 250 ml seaprep 0.5% --> 8x square plate  -W (No selection)
# 61_QWD = 250 ml seaprep 0.5% --> 8x square plate  -W/-H/-A/Aba 1/2; (Strong selection)
#
P170_1_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/S61_SWD.ppi.txt.3.15.refine.txt",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
P170_2_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/S61_PWD.ppi.txt.3.15.refine.txt",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
P170_3_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/S61_PQD.ppi.txt.3.15.refine.txt",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)

bP170_1_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_SWD.ppi.txt",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
bP170_2_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PWD.ppi.txt",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
bP170_3_3 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)

bisP170_3_3 = InteractionScores(bP170_1_3,bP170_3_3,1.0,4) 
bnisP170_3_3 = NewInteractionScores(bP170_1_3,bP170_3_3,1.0,4)        
bisP170_3_3_ = InteractionScores(bP170_2_3,bP170_3_3,1.0,4) 
bnisP170_3_3_ = NewInteractionScores(bP170_2_3,bP170_3_3,1.0,4)                              

# correlation between liquid culture and solid culture (after liquid culture)
CorMatrix(P170_1_3[P166,P166], P170_2_3[P166,P166],LogScale = TRUE)  #  0.8412698        
                

slP170_2_3 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21098-61-PWD_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
slP170_3_3 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21099-61-PQD_S3_sorted_flag113tidyformat.xlsx.txt",0.5, 4)# BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
sl_isP170_2_3 = InteractionScores(slP170_2_3,slP170_3_3,1.0,4)    
simple_plotMatrix(sl_isP170_2_3[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_2_3[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_2_3[120:161,120:161],0,7)

# [P170 - Exp4]      
# 2017-08-28_MiSeq; [68]; P170-exp4 (J68 samples = submitted on 18/08/2017 = P168xP168 library)   ## 1999636 + 1943064 + 1372704 + 1247668 = 6563072         
#
# S68_SWD  = 250 ml seaprep 0.5%      -W 
# S68_SA4D = 250 ml seaprep 0.5%      -W/ Aba 1/4
# S68_SA8D = 250 ml seaprep 0.5%      -W/ Aba 1/8
# S68_SQD  = 250 ml seaprep 0.5%      -W/-H/-A
#
bP170_1_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
bP170_2_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SA4D.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/4; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
bP170_3_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SA8D.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/8; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
bP170_4_4 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/-H/-A; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)

bisP170_2_4 = InteractionScores(bP170_1_4,bP170_2_4,1.0,4) 
bnisP170_2_4 = NewInteractionScores(bP170_1_4,bP170_2_4,1.0,4)        

bisP170_3_4 = InteractionScores(bP170_1_4,bP170_3_4,1.0,4) 
bnisP170_3_4 = NewInteractionScores(bP170_1_4,bP170_3_4,1.0,4)

bisP170_4_4 = InteractionScores(bP170_1_4,bP170_4_4,1.0,4) 
bnisP170_4_4 = NewInteractionScores(bP170_1_4,bP170_4_4,1.0,4)   
                                                                   
slP170_1_4 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21762-68-WD_S1_sorted_flag113tidyformat.xlsx.txt",0.5, 4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
slP170_2_4 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21763-68-A4D_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 4)# BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
slP170_3_4 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21764-68-A8D_S3_sorted_flag113tidyformat.xlsx.txt",0.5, 4)# BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
slP170_4_4 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21765-68-QD_S4_sorted_flag113tidyformat.xlsx.txt",0.5, 4)# BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
sl_isP170_2_4 = InteractionScores(slP170_1_4,slP170_2_4,1.0,4) 
simple_plotMatrix(sl_isP170_2_4[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_2_4[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_2_4[120:161,120:161],0,7)   
sl_isP170_3_4 = InteractionScores(slP170_1_4,slP170_3_4,1.0,4)
simple_plotMatrix(sl_isP170_3_4[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_3_4[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_3_4[120:161,120:161],0,7)    
sl_isP170_4_4 = InteractionScores(slP170_1_4,slP170_4_4,1.0,4)
simple_plotMatrix(sl_isP170_4_4[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_4_4[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_4_4[120:161,120:161],0,7)    

                                                             



# [P170 - Exp5]      
# 2017-08-30_MiSeq; [65] P170-exp5 (J65 samples = submitted on 18/08/2017 = P168xP168 library)   ## 1836282 + 2108754 + 1870105 = 5815141   
#
# S65_SWD  = 250 ml seaprep 0.5%      -W 
# S65_SAD = 250 ml seaprep 0.5%      -W/ Aba 1/8
# S65_SQD  = 250 ml seaprep 0.5%      -W/-H/-A
#              
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SWD_R1 S65_SWD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SWD > qjobs/qjob_2017-08-30_MiSeq_S65_SWD.sh # SWD (P170 Seaprep -W)            # ?
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SAD_R1 S65_SAD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SAD > qjobs/qjob_2017-08-30_MiSeq_S65_SAD.sh # SAD (P170 Seaprep +A?)           # ?
# sh Y2H_Blastn.sh 2017-08-30_MiSeq S65_SQD_R1 S65_SQD_R2 ../data/P170_4_library_MGj5615-Jun-2017122014.-100 S65_SQD > qjobs/qjob_2017-08-30_MiSeq_S65_SQD.sh # SQD (P170 Seaprep -Q)            # ?
# qsub -cwd -q long-sl65 -l virtual_free=16G,h_rt=72:00:00 -M jae-seong.yang@crg.es -m abe -v CUSTOM_EMAIL=yes qjobs/qjob_2017-08-30_MiSeq_S65_SWD.sh 
# qsub -cwd -q long-sl65 -l virtual_free=16G,h_rt=72:00:00 -M jae-seong.yang@crg.es -m abe -v CUSTOM_EMAIL=yes qjobs/qjob_2017-08-30_MiSeq_S65_SAD.sh 
# qsub -cwd -q long-sl65 -l virtual_free=16G,h_rt=72:00:00 -M jae-seong.yang@crg.es -m abe -v CUSTOM_EMAIL=yes qjobs/qjob_2017-08-30_MiSeq_S65_SQD.sh
bP170_1_5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W;             # ==> 1836282 out of ? (48.4%)
bP170_2_5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SAD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/8;  # ==>  out of  (43.7%)
bP170_3_5 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/-H/-A; R75_01 # ==>  out of  (47.4%)

bisP170_2_5 = InteractionScores(bP170_1_5,bP170_2_5,1.0,4)
bnisP170_2_5 = NewInteractionScores(bP170_1_5,bP170_2_5,1.0,4)

bisP170_3_5 = InteractionScores(bP170_1_5,bP170_3_5,1.0,4) 
bnisP170_3_5 = NewInteractionScores(bP170_1_5,bP170_3_5,1.0,4)
                                                                

slP170_1_5 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21766-65-WD_S1_sorted_flag113tidyformat.xlsx.txt",0.5, 4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
slP170_2_5 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21767-65-AD_S2_sorted_flag113tidyformat.xlsx.txt",0.5, 4)# BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
slP170_3_5 = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X163/Ref_X163nostp_seq_21768-65-QD_S3_sorted_flag113tidyformat.xlsx.txt",0.5, 4)# BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
sl_isP170_2_5 = InteractionScores(slP170_1_5,slP170_2_5,1.0,4) 
simple_plotMatrix(sl_isP170_2_5[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_2_5[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_2_5[120:161,120:161],0,7)   
sl_isP170_3_5 = InteractionScores(slP170_1_5,slP170_3_5,1.0,4)
simple_plotMatrix(sl_isP170_3_5[1:60,1:60],0,7)
simple_plotMatrix(sl_isP170_3_5[60:120,60:120],0,7)
simple_plotMatrix(sl_isP170_3_5[120:161,120:161],0,7)    



                                               

1/4A
CorMatrix(bisP170_2_5, bisP170_2_4)     # 0.4260798
CorMatrix(bisP170_2_5, bisP170_3_4)   # 0.4127579

Qudraple cases:
CorMatrix(bnisP170_3_5, bnisP170_4_4)  # 0.8386545
CorMatrix(bnisP170_3_5, bnisP170_3_3)  # 0.6671311  
CorMatrix(bnisP170_3_5, bnisP170_3_3_)  # 0.6683618 
CorMatrix(bnisP170_4_4, bnisP170_3_3_)  # 0.6661674
CorMatrix(bnisP170_3_5, bnisP170_2_2)  # 0.7096627

CorMatrix(bisP170_3_5, bisP170_4_4)  # 0.7132838
CorMatrix(bisP170_3_5, bisP170_3_3)  # 0.5947822
CorMatrix(bisP170_3_5, bisP170_3_3_)  # 0.607885
CorMatrix(bisP170_4_4, bisP170_3_3_) # 0.6019529           
CorMatrix(bisP170_3_5, bisP170_2_2)  # 0.5174159 


### Calibration set                                                               
simple_plotMatrix( bnisP170_2_1[CALIBRATION,CALIBRATION],0 )   # -W/-H/-A/Aba 1/2
simple_plotMatrix( bnisP170_2_2[CALIBRATION,CALIBRATION],0 )   # -W/-H/-A/Aba 1/2 
simple_plotMatrix( bnisP170_3_3[CALIBRATION,CALIBRATION],0 )   # -W/-H/-A/Aba 1/2
simple_plotMatrix( bnisP170_3_3_[CALIBRATION,CALIBRATION],0 )  # -W/-H/-A/Aba 1/2

simple_plotMatrix( bnisP170_2_4[CALIBRATION,CALIBRATION],0 )  # -W/-H/-A/Aba 1/2
simple_plotMatrix( bnisP170_3_4[CALIBRATION,CALIBRATION],0 )  # -W/-H/-A/Aba 1/2
simple_plotMatrix( bnisP170_4_4[CALIBRATION,CALIBRATION],0 )  # -W/-H/-A/Aba 1/2
simple_plotMatrix( bnisP170_2_5[CALIBRATION,CALIBRATION],0 )  # -W/-H/-A/Aba 1/2
simple_plotMatrix( bnisP170_3_5[CALIBRATION,CALIBRATION],0 )  # -W/-H/-A/Aba 1/2



simple_plotMatrix( bnisP170_2_5[CALIBRATION,CALIBRATION],0 ) 

simple_plotMatrix( bisP170_2_1[CALIBRATION,CALIBRATION],0 )

                                    

### For Figure ###                   

setwd("/Users/jyang/Dropbox_CRG/matrixppi/figures")                             

bait21 = c(45,24,13,11,36,56,15,35,48,6,61,59,4,32,23,63,22,27,40,28,69)
prey21 = c(47,25,16,13,35,55,41,37,46,7,62,57,5,33,34,44,54,10,43,75,20)
bait19 = c(45,24,13,11,36,56,15,35,48,6,61,59,4,32,23,63,22,27,69) # without PLEKHG7, WDR7
prey19 = c(47,25,16,13,35,55,41,37,46,7,62,57,5,33,34,44,54,10,20) # without PLEKHG7, WDR7

# #75 = R75_41,WDR7,WD repeat domain 7
# #43 = R75_72,PLEKHG7,pleckstrin homology domain containing, family G (with RhoGef domain) member 7
bait76 = setdiff( 1:78, c(43,75) ) 
prey76 = setdiff( 1:78, c(43,75) )

# #2 = R75_11,AP2B1,adaptor-related protein complex 2, beta 1 subunit
# #19 = R75_30,CTBP1,C-terminal binding protein 1
# #51 = R75_52,RBBP8,retinoblastoma binding protein 8
# #26 = R75_58,FHL3,four and a half LIM domains 3
# #68 = R75_70,TEX11,testis expressed 11 
bait71 = setdiff( 1:78, c(2,19,26,43,51,68,75) ) 

simple_plotLog2Matrix(bm1_1[bait76,prey76],1) 
simple_plotLog2Matrix(bm1_2[bait76,prey76],1) 
simple_plotLog2Matrix(bm1_3[bait76,prey76],1)   

simple_plotLog2Matrix(bm1_1[bait71,prey76],1) 
simple_plotLog2Matrix(bm1_2[bait71,prey76],1) 
simple_plotLog2Matrix(bm1_3[bait71,prey76],1)
simple_plotLog2Matrix(bm1_4[bait71,prey76],1)
simple_plotLog2Matrix((bm1_1+bm1_2+bm1_3+bm1_4)[bait71,prey76],1)   

pdf("46_WC.pdf"); simple_plotLog2Matrix(bm1_1,1.0,7); dev.off();
pdf("46_AC.pdf"); simple_plotLog2Matrix(bm2_1,1.0,7); dev.off();
pdf("46_WRC.pdf"); simple_plotLog2Matrix(bm3_1,1.0,7); dev.off();
pdf("46_WAC.pdf"); simple_plotLog2Matrix(bm4_1,1.0,7); dev.off();
pdf("46.IS1.pdf"); bis2_1 = InteractionScores(bm1_1,bm2_1,1.0,7); dev.off();
pdf("46.IS2.pdf"); bis4_1 = InteractionScores(bm3_1,bm4_1,1.0,7); dev.off();   
pdf("46.NIS1.pdf"); bnis2_1 = NewInteractionScores(bm1_1,bm2_1,1.0,7); dev.off();
pdf("46.NIS2.pdf"); bnis4_1 = NewInteractionScores(bm3_1,bm4_1,1.0,7); dev.off();   
CorMatrix(bis2_1,bis4_1,"Interaction Score 1","Interaction Score 2")                   # 0.9217131
CorMatrix(bnis2_1,bnis4_1,"Interaction Score 1","Interaction Score 2")                 # 0.9610908
simple_plotMatrix(bnis2_1,1,4,TRUE,TRUE)   
simple_plotMatrix(bis2_1[bait21,prey21],1,9)   
simple_plotMatrix(bnis2_1[bait21,prey21],1,9)    
                                               


png("51_WC.png"); simple_plotLog2Matrix(bm1_2,1.0,7); dev.off();
png("51_AC.png"); simple_plotLog2Matrix(bm2_2,1.0,7); dev.off();
png("51_QC.png"); simple_plotLog2Matrix(bm3_2,1.0,7); dev.off();
png("51.IS1.png"); bis2_2 = InteractionScores(bm1_2,bm2_2,1.0,7); dev.off();
png("51.IS2.png"); bis3_2 = InteractionScores(bm1_2,bm3_2,1.0,7); dev.off();   
png("51.NIS1.png"); bnis2_2 = NewInteractionScores(bm1_2,bm2_2,1.0,7); dev.off();
png("51.NIS2.png"); bnis3_2 = NewInteractionScores(bm1_2,bm3_2,1.0,7); dev.off();
CorMatrix(bis2_2,bis3_2,"Interaction Score 1","Interaction Score 2")				# 0.6835282
CorMatrix(bnis2_2,bnis3_2,"Interaction Score 1","Interaction Score 2")              # 0.8438049
simple_plotMatrix(bnis2_2,1,4,TRUE,TRUE)   
simple_plotMatrix(bis2_2[bait21,prey21],1,9)   
simple_plotMatrix(bnis2_2[bait21,prey21],1,9)  
                     
par(mfrow=c(1,3))
simple_plotMatrix(bm1_2[bait21,prey21],1,9)
simple_plotMatrix(bm2_2[bait21,prey21],1,9)
simple_plotMatrix(bm3_2[bait21,prey21],1,9)
par(mfrow=c(1,1))    

png("53_WC.png"); simple_plotLog2Matrix(bm1_3,1.0,7); dev.off();
png("53_AC.png"); simple_plotLog2Matrix(bm2_3,1.0,7); dev.off();
png("53_QC.png"); simple_plotLog2Matrix(bm3_3,1.0,7); dev.off();
png("53.IS1.png"); bis2_3 = InteractionScores(bm1_3,bm2_3,1.0,7); dev.off();
png("53.IS2.png"); bis3_3 = InteractionScores(bm1_3,bm3_3,1.0,7); dev.off();   
png("53.NIS1.png"); bnis2_3 = NewInteractionScores(bm1_3,bm2_3,1.0,7); dev.off();
png("53.NIS2.png"); bnis3_3 = NewInteractionScores(bm1_3,bm3_3,1.0,7); dev.off();
CorMatrix(bis2_3,bis3_3,"Interaction Score 1","Interaction Score 2")     			# 0.722062
CorMatrix(bnis2_3,bnis3_3,"Interaction Score 1","Interaction Score 2")                          
simple_plotMatrix(bis2_3,1,4,TRUE,TRUE)  
simple_plotMatrix(bnis2_3,1,4,TRUE,TRUE)   

par(mfrow=c(1,3))
simple_plotMatrix(log2(bm1_3[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm2_3[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm3_3[bait21,prey21]+1),1,9)
par(mfrow=c(1,1))


png("64_SWC.png"); simple_plotLog2Matrix(bm1_4,1.0,7); dev.off();
png("64_SA4C.png"); simple_plotLog2Matrix(bm2_4,1.0,7); dev.off();
png("64_SA8C.png"); simple_plotLog2Matrix(bm3_4,1.0,7); dev.off();
png("64_SQC.png"); simple_plotLog2Matrix(bm4_4,1.0,7); dev.off();
png("64.IS1.png"); bis2_4 = InteractionScores(bm1_4,bm2_4,1.0,7); dev.off();
png("64.IS2.png"); bis3_4 = InteractionScores(bm1_4,bm3_4,1.0,7); dev.off();   
png("64.IS3.png"); bis4_4 = InteractionScores(bm1_4,bm4_4,1.0,7); dev.off();   
png("64.NIS1.png"); bnis2_4 = NewInteractionScores(bm1_4,bm2_4,1.0,7); dev.off();
png("64.NIS2.png"); bnis3_4 = NewInteractionScores(bm1_4,bm3_4,1.0,7); dev.off();
png("64.NIS3.png"); bnis4_4 = NewInteractionScores(bm1_4,bm4_4,1.0,7); dev.off();
CorMatrix(bis2_4,bis3_4,"Interaction Score 1","Interaction Score 2")     			# 0.7223925
CorMatrix(bis2_4,bis4_4,"Interaction Score 1","Interaction Score 2")     			# 0.689696
CorMatrix(bnis2_4,bnis3_4,"Interaction Score 1","Interaction Score 2")				# 0.9287979                          
CorMatrix(bnis2_4,bnis4_4,"Interaction Score 1","Interaction Score 2") 				# 0.793468                         

par(mfrow=c(1,3))
simple_plotMatrix(log2(bm1_4[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm2_4[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm3_4[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm4_4[bait21,prey21]+1),1,9)  

simple_plotMatrix(log2(bis2_4[bait21,prey21]+1),1,9) 
simple_plotMatrix(bis4_4[bait21,prey21],1,9) 
simple_plotMatrix(bnis4_4[bait21,prey21],1,9) 
par(mfrow=c(1,1))

                                                                                               
par(mfrow=c(3,2))
CorMatrix(bis2_1,bis4_1,"Interaction Score 1","Interaction Score 2")                # 0.9217131
CorMatrix(bnis2_1,bnis4_1,"Interaction Score 1","Interaction Score 2")              # 0.9610908
CorMatrix(bis2_2,bis3_2,"Interaction Score 1","Interaction Score 2")				# 0.6835282
CorMatrix(bnis2_2,bnis3_2,"Interaction Score 1","Interaction Score 2")              # 0.8438049
CorMatrix(bis2_3,bis3_3,"Interaction Score 1","Interaction Score 2")     			# 0.722062
CorMatrix(bnis2_3,bnis3_3,"Interaction Score 1","Interaction Score 2")
             

                            
## Plot Raw PPI-Pair reads for 21 Reference set - No selection
# NOTE; PLEKHG7, WDR7 were not cloned..
pdf("21pair.bm1_1.pdf"); simple_plotMatrix(log2(bm1_1[bait21,prey21]+1),1,9); dev.off()
pdf("21pair.bm1_2.pdf"); simple_plotMatrix(log2(bm1_2[bait21,prey21]+1),1,9); dev.off() 
pdf("21pair.bm1_3.pdf"); simple_plotMatrix(log2(bm1_3[bait21,prey21]+1),1,9); dev.off()  
pdf("21pair.bm1_4.pdf"); simple_plotMatrix(log2(bm1_4[bait21,prey21]+1),1,9); dev.off()  

## Plot Raw PPI-Pair reads for 21 Reference set - Selection   
# NOTE; PLEKHG7, WDR7 were not cloned.. 
simple_plotMatrix(log2(bm2_1[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm2_2[bait21,prey21]+1),1,9)
simple_plotMatrix(log2(bm2_3[bait21,prey21]+1),1,9) 
simple_plotMatrix(log2(bm2_4[bait21,prey21]+1),1,9) 

## Plot Raw PPI-Pair reads for 21 Reference set - Interaction   
# NOTE; PLEKHG7, WDR7 were not cloned.. 
InteractionScores(bm1_1[bait21,prey21], bm2_1[bait21,prey21], 1 , 9 )
InteractionScores(bm1_2[bait21,prey21], bm2_2[bait21,prey21], 1 , 9 )
InteractionScores(bm1_3[bait21,prey21], bm2_3[bait21,prey21], 1 , 9 )
InteractionScores(bm1_4[bait21,prey21], bm2_4[bait21,prey21], 1 , 9 )
     
## Plot Raw PPI-Pair reads for 21 Reference set - NewInteraction   
# NOTE; PLEKHG7, WDR7 were not cloned.. 
NewInteractionScores(bm1_1[bait21,prey21], bm2_1[bait21,prey21], 1 , 9 )
NewInteractionScores(bm1_2[bait21,prey21], bm2_2[bait21,prey21], 1 , 9 )
NewInteractionScores(bm1_3[bait21,prey21], bm2_3[bait21,prey21], 1 , 9 )
NewInteractionScores(bm1_4[bait21,prey21], bm2_4[bait21,prey21], 1 , 9 )
           
NewInteractionScores(m1_1[bait21,prey21], m2_1[bait21,prey21], 1 , 9 )
NewInteractionScores(m1_2[bait21,prey21], m2_2[bait21,prey21], 1 , 9 )
NewInteractionScores(m1_3[bait21,prey21], m2_3[bait21,prey21], 1 , 9 )
NewInteractionScores(m1_4[bait21,prey21], m2_4[bait21,prey21], 1 , 9 )


## Plot Deviation from expected matrix PPI-Pair reads for 21 Reference set -    
# NOTE; PLEKHG7, WDR7 were not cloned.. 
ep1_1 = Epistasis(bm1_1[bait21,prey21],0.0) 
ep1_2 = Epistasis(bm1_2[bait21,prey21],0.0) 
ep1_3 = Epistasis(bm1_3[bait21,prey21],0.0) 
ep1_4 = Epistasis(bm1_4[bait21,prey21],0.0) 

ep2_1 = Epistasis(bm2_1[bait21,prey21],0.1) 
ep2_2 = Epistasis(bm2_2[bait21,prey21],0.1) 
ep2_3 = Epistasis(bm2_3[bait21,prey21],0.1) 
ep2_4 = Epistasis(bm2_4[bait21,prey21],0.1)
 
ep1_1 = Epistasis(bm1_1,0.0)   
                             
CorMatrix( ep1_1, ep1_2 )
CorMatrix( ep1_1, ep1_3 )
CorMatrix( ep1_1, ep1_4 )
CorMatrix( ep1_2, ep1_3 )

CorMatrix( ep2_1, ep2_2 )
CorMatrix( ep2_1, ep2_3 )
CorMatrix( ep2_1, ep2_4 )
CorMatrix( ep2_2, ep2_3 )


                           

## Plot Deviation from epistatic interaction PPI-Pair reads for 21 Reference set -    
# NOTE; PLEKHG7, WDR7 were not cloned.. 
ep1_1 = Epistasis(bm2_1[bait21,prey21],0.0) > 2
simple_plotMatrix(ep1_1,1)
ep1_1 = Epistasis(bm2_2[bait21,prey21],0.0) > 2
simple_plotMatrix(ep1_1,1)
ep1_1 = Epistasis(bm2_3[bait21,prey21],0.0) > 2
simple_plotMatrix(ep1_1,1)
ep1_1 = Epistasis(bm2_4[bait21,prey21],0.0) > 2
simple_plotMatrix(ep1_1,1)


## Plot Deviation from expected matrix PPI-Pair reads for 21 Reference set -    
# NOTE; PLEKHG7, WDR7 were not cloned.. 
ep1_1 = Epistasis(bm1_1[bait19,prey19],0.0) 
ep1_2 = Epistasis(bm1_2[bait19,prey19],0.0) 
ep1_3 = Epistasis(bm1_3[bait19,prey19],0.0) 
ep1_4 = Epistasis(bm1_4[bait19,prey19],0.0) 

ep2_1 = Epistasis(bm2_1[bait19,prey19],0.1) 
ep2_2 = Epistasis(bm2_2[bait19,prey19],0.1) 
ep2_3 = Epistasis(bm2_3[bait19,prey19],0.1) 
ep2_4 = Epistasis(bm2_4[bait19,prey19],0.1)

ep1_1 = Epistasis(bm2_1[bait19,prey19],0.0) > 1.645  
simple_plotMatrix(ep1_1,1)
ep1_1 = Epistasis(bm2_2[bait19,prey19],0.0) > 1.645  
simple_plotMatrix(ep1_1,1)
ep1_1 = Epistasis(bm2_3[bait19,prey19],0.0) > 1.645  
simple_plotMatrix(ep1_1,1)
ep1_1 = Epistasis(bm2_4[bait19,prey19],0.0) > 1.645  
simple_plotMatrix(ep1_1,1)

# For all proteins
ep2_1 = Epistasis(bm2_1,0.0) >  1.645     
simple_plotMatrix(ep2_1,1)            
ep4_1 = Epistasis(bm4_1,0.0) >  1.645 
simple_plotMatrix(ep4_1,1)            
ep2_2 = Epistasis(bm2_2,0.0) >  1.645     
simple_plotMatrix(ep2_2,1)  
ep2_3 = Epistasis(bm2_3,0.0) >  1.645     
simple_plotMatrix(ep2_3,1)  
ep2_4 = Epistasis(bm2_4,0.0) >  1.645     
simple_plotMatrix(ep2_4,1)  
simple_plotMatrix(ep2_1+ep4_1+ep2_2+ep2_3+ep2_4,1)

ep1_1 = Epistasis(bm2_1,0.0) >  2     
simple_plotMatrix(ep1_1,1)  
ep1_2 = Epistasis(bm2_2,0.0) >  2     
simple_plotMatrix(ep1_2,1)  
ep1_3 = Epistasis(bm2_3,0.0) >  2     
simple_plotMatrix(ep1_3,1)  
ep1_4 = Epistasis(bm2_4,0.0) >  2     
simple_plotMatrix(ep1_4,1)  
simple_plotMatrix(ep1_1+ep1_2+ep1_3+ep1_4,1)

epP170_2_1 = Epistasis(bP170_2_1,0.0) >  1.645     
simple_plotMatrix(epP170_2_1,1)  
epP170_2_2 = Epistasis(bP170_2_2,0.0) >  1.645     
simple_plotMatrix(epP170_2_2,1)  
epP170_3_3 = Epistasis(bP170_3_3,0.0) >  1.645     
simple_plotMatrix(epP170_3_3,1)  

simple_plotMatrix(epP170_2_1+epP170_2_2+epP170_3_3,1)   
simple_plotMatrix((epP170_2_1+epP170_2_2+epP170_3_3)>=2,1,4)  



#===================================================================

plot(log2(bm1_1[bait21,prey21]+1),log2(bm1_2[bait21,prey21]+1))
plot(log2(bm1_3[bait21,prey21]+1),log2(bm1_4[bait21,prey21]+1))
plot(log2(bm1_1[bait21,prey21]+1),log2(bm1_2[bait21,prey21]+1))  





par(mfrow=c(1,1))



par(mfrow=c(1,1))
   
# one tail Z-score: Z > 1.645 for alpha = 0.05
                                              
                                                   
hist(bnis2_4,breaks =100)
hist(bnis2_4,breaks =100,xlim=c(0,7),ylim=c(0,100))  

sum(bnis2_4>1.645) # 183 PPI pair
dim(bnis2_4) # 78 * 78 = 6084
183/6084.0   # 0.0300789
                                    
par(mfrow=c(1,2))
bnis_overlap = (bnis2_1>1.645) + (bnis2_2>1.645) + (bnis2_3>1.645) + (bnis2_4>1.645) + (bnis3_4>1.645)  + (bnis4_4>1.645)          
hist(bnis_overlap,ylim=c(0,300))

bis_overlap = (bis2_1>1.645) + (bis2_2>1.645) + (bis2_3>1.645) + (bis2_4>1.645) + (bis3_4>1.645)  + (bis4_4>1.645)          
hist(bis_overlap,ylim=c(0,300))     
                              
simple_plotMatrix(bnis_overlap,1)
simple_plotMatrix(bis_overlap,1)

                                  



## P170

bnisP170_overlap = (bnisP170_1>1.645) + (bnisP170_2>1.645) + (bnisP170_3>1.645) + (bnisP170_4>1.645)   
bnisP170_overlap = bnisP170_1 * bnisP170_2 * bnisP170_3 * bnisP170_4        
simple_plotMatrix(bnisP170_overlap,1)                           
simple_plotMatrix(bnisP170_overlap,1,4,TRUE,TRUE)      
hist(bnisP170_overlap,xlim=c(0.5,3), ylim=c(0,1000))  
simple_plotLog2Matrix(bnisP170_overlap,1,4,TRUE,TRUE)
#===========================================================================


                                                                                                 
#===========================================================================
### RNA - PROTEIN INTERACTION                                               
#===========================================================================

rpi1_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-08-24_MiSeq/Blastn/S1_W.rpi.txt",0.5,4) # no selection ; 		    "sum = 3199976"     
rpi1_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-08-24_MiSeq/Blastn/S2_WH.rpi.txt",0.5,4) #  selection ; 		    "sum = 3988896"     
         
rpi1_1m = rpi1_1[1:350,351:451]
rpi1_2m = rpi1_2[1:350,351:451]
     
View( rpi1_2[1:350,351:451] )   

ep_rpi1_1 = Epistasis(rpi1_1,0.0) >  1.645     
simple_plotMatrix(ep_rpi1_1,1,1)
                                    
ep_rpi1_2 = Epistasis(rpi1_2,0.0) >  1.645     
simple_plotMatrix(ep_rpi1_2,1,1)
                                   

> simple_plotMatrix(ep_rpi1_2[1:350,351:451],1,1,FALSE,FALSE)
> simple_plotMatrix(ep_rpi1_1[1:350,351:451],1,1,FALSE,FALSE)
                
simple_plotLog2Matrix(ep_rpi1_1[1:350,351:451],1,1,FALSE,FALSE) 
simple_plotLog2Matrix(ep_rpi1_2[1:350,351:451],1,1,FALSE,FALSE)
                                                             

bnis_rpi1_2 = NewInteractionScores(rpi1_1,rpi1_2,1.0,2)                

bnis_rpi1_2 = NewInteractionScores(rpi1_1[1:350,351:451],rpi1_2[1:350,351:451],1.0,2)      


simple_plotMatrix(bnis_rpi1_2[1:350,351:451],1,1,FALSE,FALSE)


       
RNA
SECIS RNA		[2PJP:B]  	controls
LSL RNA		[2L3C:B]  	controls
1MMS:C		[1MMS:C]  	controls     #
25 nt RNA hairpin		[2ANN:B]  	controls
pre rRNA 		[1RKJ:B]  	controls
pMS22KH empty			controls

Protein
Z97	Adarb1	25367	Rat	2L3C:A 
Z98	SelB 	948103	E.coli	2PJP:A
Z99	rplK	897484	Thermotoga maritima	1MMS:A   #
Z100	NCL	4691	Mesocricetus auratus	1RKJ:A
Z101	Nova-1	4857	H.sapiens	2ANN:A




# 2017-09-07_MiSeq; (For Nele)
N6_W_R1.fastq.gz
N6_W_R2.fastq.gz
N6_WH_R1.fastq.gz
N6_WH_R2.fastq.gz                 

# sh Y3H_Blastn.sh 2017-09-07_MiSeq N6_W_R1 N6_W_R2 ../data/1st_set.-100 N6_W > qjobs/N6_W_2017-09-07_MiSeq.sh             # rpi reads = ?
# sh Y3H_Blastn.sh 2017-09-07_MiSeq N6_WH_R1 N6_WH_R2 ../data/1st_set.-100 N6_WH > qjobs/N6_WH_2017-09-07_MiSeq.sh          # rpi reads = ? , when allow both direction = ?
# qsub -cwd -q long-sl65 -l virtual_free=16G,h_rt=72:00:00 -M jae-seong.yang@crg.es -m abe -v CUSTOM_EMAIL=yes qjobs/N6_W_2017-09-07_MiSeq.sh 
# qsub -cwd -q long-sl65 -l virtual_free=16G,h_rt=72:00:00 -M jae-seong.yang@crg.es -m abe -v CUSTOM_EMAIL=yes qjobs/N6_WH_2017-09-07_MiSeq.sh

rpi2_1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-09-07_MiSeq/Blastn/N6_W.rpi.txt",0.5,4) # no selection ; 		    "sum = 3413429"     
rpi2_2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-09-07_MiSeq/Blastn/N6_WH.rpi.txt",0.5,4) #  selection ; 		    "sum = 4052587 ? 4044055"     
m1 = rpi2_1[1:350,351:451]
m2 = rpi2_2[1:350,351:451]
bis_rpi2_2 = InteractionScores(rpi2_1[1:350,351:451],rpi2_2[1:350,351:451],1.0,2) 
bis_rpi2_2_reorder = InteractionScores(rpi2_1[1:350,351:451][RNA_INDEX,RBP_INDEX],rpi2_2[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2)   
bnis_rpi2_2 = NewInteractionScores(rpi2_1[1:350,351:451],rpi2_2[1:350,351:451],1.0,2)       
bnis_rpi2_2_reorder = NewInteractionScores(rpi2_1[1:350,351:451][RNA_INDEX,RBP_INDEX],rpi2_2[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2)       

cis_rpi2_2 = ColumnWiseInteractionScores(rpi2_1[1:350,351:451],rpi2_2[1:350,351:451],1.0,2)   
cis_rpi2_2_reorder = ColumnWiseInteractionScores(rpi2_1[1:350,351:451][RNA_INDEX,RBP_INDEX] ,rpi2_2[1:350,351:451][RNA_INDEX,RBP_INDEX] ,1,2)  
    
simple_plotMatrix(rpi2_2[1:350,351:451][RNA_INDEX,RBP_INDEX][,1:19],0,4)

                
> CorMatrix(bnis_rpi1_2, bnis_rpi2_2)
[1] 0.3813725        

SRSF5 ==> R127_28s_rn28s1_681_720

rplK (50S ribosomal protein L11) ==> R139_28s_rn28s1_921_960    


RPL22
LARP7
                                                         
                                                                                                                         
### New experiment

bS1_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S1_W.rpi.txt",0.5,4)
bS1_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S2_WH.rpi.txt",0.5,4)   
bS1_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4)
bS1_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4)
is_S1_RPI = InteractionScores(bS1_W_rpi[1:350,351:451],bS1_WH_rpi[1:350,351:451],1,4)             
is_S1_RPI_reorder = InteractionScores(bS1_W_rpi[1:350,351:451][RNA_INDEX,RBP_INDEX],bS1_WH_rpi[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2) 
               
simple_plotMatrix(bS1_WH_rpi[1:350,351:451][RNA_INDEX,RBP_INDEX][,1:19],0,4)

m1 = bS1_W_rpi[1:350,351:451]
m2 = bS1_WH_rpi[1:350,351:451]
cis_S1_RPI = ColumnWiseInteractionScores(bS1_W_rpi[1:350,351:451],bS1_WH_rpi[1:350,351:451],1,2)                       
cis_S1_RPI = ColumnWiseInteractionScores(bS1_W_rpi[1:350,351:451][RNA_INDEX,RBP_INDEX] ,bS1_WH_rpi[1:350,351:451][RNA_INDEX,RBP_INDEX] ,1,2) 
simple_plotMatrix(cis_S1_RPI>50,0,2)   
                                              
# over CIS=50
simple_plotMatrix(cis_S1_RPI[1:60,1:50]>50,0,2)   
                                                
RBP_INDEX = c(1,59,85,19,27,47,68,15,34,9,18,37,46,48,54,57,39,38,87,3,4,5,6,7,8,10,11,12,13,14,16,17,20,21,22,23,24,25,26,28,29,30,31,32,33,35,36,40,41,42,43,44,45,49,50,51,52,53,55,56,58,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,88,89,90,91,92,93,94,95,96,97,101,2,98,99,100)
RNA_INDEX = c(111,222,284,295,306,317,328,339,350,11,22,33,44,55,66,77,88,99,110,122,133,144,155,166,177,188,199,210,221,233,244,255,266,277,279,280,281,282,283,285,286,287,288,289,290,291,292,293,294,296,297,298,299,300,301,302,303,304,305,307,308,309,310,311,312,313,314,315,316,318,319,320,321,322,323,324,325,326,327,329,330,331,332,333,334,335,336,337,338,340,341,342,343,344,345,346,347,348,349,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,45,46,47,48,49,50,51,52,53,54,56,57,58,59,60,61,62,63,64,65,67,68,69,70,71,72,73,74,75,76,78,79,80,81,82,83,84,85,86,87,89,90,91,92,93,94,95,96,97,98,100,101,102,103,104,105,106,107,108,109,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,134,135,136,137,138,139,140,141,142,143,145,146,147,148,149,150,151,152,153,154,156,157,158,159,160,161,162,163,164,165,167,168,169,170,171,172,173,174,175,176,178,179,180,181,182,183,184,185,186,187,189,190,191,192,193,194,195,196,197,198,200,201,202,203,204,205,206,207,208,209,211,212,213,214,215,216,217,218,219,220,223,224,225,226,227,228,229,230,231,232,234,235,236,237,238,239,240,241,242,243,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,267,268,269,270,271,273,278,276,272,274,275)                                       

cis_S1_RPI_order = cis_S1_RPI[RNA_INDEX,RBP_INDEX]                               

simple_plotMatrix(cis_S1_RPI_order[1:60,1:50]>50,0,7,filename="cis_S1_RPI.over50.1-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[61:120,1:50]>50,0,7,filename="cis_S1_RPI.over50.2-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[121:180,1:50]>50,0,7,filename="cis_S1_RPI.over50.3-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[181:240,1:50]>50,0,7,filename="cis_S1_RPI.over50.4-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[241:300,1:50]>50,0,7,filename="cis_S1_RPI.over50.5-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[301:350,1:50]>50,0,7,filename="cis_S1_RPI.over50.6-1.pdf")

simple_plotMatrix(cis_S1_RPI_order[1:60,51:101]>50,0,7,filename="cis_S1_RPI.over50.1-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[61:120,51:101]>50,0,7,filename="cis_S1_RPI.over50.2-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[121:180,51:101]>50,0,7,filename="cis_S1_RPI.over50.3-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[181:240,51:101]>50,0,7,filename="cis_S1_RPI.over50.4-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[241:300,51:101]>50,0,7,filename="cis_S1_RPI.over50.5-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[301:350,51:101]>50,0,7,filename="cis_S1_RPI.over50.6-2.pdf")


simple_plotMatrix(cis_S1_RPI_order[1:60,1:50],0,7,filename="cis_S1_RPI.1-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[61:120,1:50],0,7,filename="cis_S1_RPI.2-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[121:180,1:50],0,7,filename="cis_S1_RPI.3-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[181:240,1:50],0,7,filename="cis_S1_RPI.4-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[241:300,1:50],0,7,filename="cis_S1_RPI.5-1.pdf")
simple_plotMatrix(cis_S1_RPI_order[301:350,1:50],0,7,filename="cis_S1_RPI.6-1.pdf")

simple_plotMatrix(cis_S1_RPI_order[1:60,51:101],0,7,filename="cis_S1_RPI.1-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[61:120,51:101],0,7,filename="cis_S1_RPI.2-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[121:180,51:101],0,7,filename="cis_S1_RPI.3-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[181:240,51:101],0,7,filename="cis_S1_RPI.4-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[241:300,51:101],0,7,filename="cis_S1_RPI.5-2.pdf")
simple_plotMatrix(cis_S1_RPI_order[301:350,51:101],0,7,filename="cis_S1_RPI.6-2.pdf")


                                                                                                                                                                     
KPR_1_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-19_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4)  # I think the samples were switched.
KPR_1_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-19_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4)  # I think the samples were switched.
is_KPR_WH_rpi = InteractionScores(log(KPR_1_W_rpi[c(1:16,18),c(17,19:36)]),log(KPR_1_WH_rpi[c(1:16,18),c(17,19:36)]),1,4)   

is_KPR_WH_rpi = InteractionScores(KPR_1_W_rpi[c(1:16,18),c(17,19:36)],KPR_1_WH_rpi[c(1:16,18),c(17,19:36)],1,4)
m1 = KPR_1_W_rpi[c(1:16,18),c(17,19:36)]
m2 = KPR_1_WH_rpi[c(1:16,18),c(17,19:36)]
nis_KPR_WH_rpi = NewInteractionScores(KPR_1_W_rpi[c(1:16,18),c(17,19:36)],KPR_1_WH_rpi[c(1:16,18),c(17,19:36)],1,4)        
pis_KPR_WH_rpi = PairInteractionScores(KPR_1_W_rpi[c(1:16,18),c(17,19:36)],KPR_1_WH_rpi[c(1:16,18),c(17,19:36)],1,4)        
cis_KPR_WH_rpi = ColumnWiseInteractionScores(KPR_1_W_rpi[c(1:16,18),c(17,19:36)],KPR_1_WH_rpi[c(1:16,18),c(17,19:36)],1,4)  

ColumnWiseInteractionScores <- function(m1,m2,alpha, fontsize=9){      
	rcnt = dim(m1)[1]
	colSumsMat1 = ( (1:rcnt * 0 + 1) %*% t(colSums(m1)) )
	colNorm_m1 = m1 /  colSumsMat1                          
	colNorm_m1[ which(colNorm_m1=="NaN")] = 0                    
	rowFreq = rowSums(colNorm_m1) / length(colnames(colNorm_m1))
	
	colSumsMat2 = ( (1:rcnt * 0 + 1) %*% t(colSums(m2)) )
	colNorm_m2 = m2 /  colSumsMat2	 
	simple_plotMatrix(colNorm_m2,0,fontsize)  
	
	refined = colNorm_m2/ rowFreq  
	refined[ which(refined == "NaN") ] = 0 
	refined[ is.na(refined) ] = 0 
	refined[ which(refined == Inf) ] = max( refined[ which(refined != Inf) ] ) + 1.0
	  
	simple_plotMatrix(refined,0,fontsize)               
	return (refined)
}

RawInteractionScores <- function(m1,m2, alpha,show=TRUE, filename=NA, cellwidth = NA, cellheight = NA){
	# m1 : no selection matrix
	# m2 : selection matrix  
	                 
	NullMatrix_matrix1 <- data.matrix(m1) + alpha  
	total_reads1 = sum(NullMatrix_matrix1)  
	freq_matrix1 = NullMatrix_matrix1 / total_reads1   
	rows1 = rowSums(freq_matrix1)
	cols1 = colSums(freq_matrix1)
	reconst_matrix1 = rows1%*%t(cols1)  # null
	
	total_reads2 = sum(m2)
	freq_matrix2 = ( m2 + alpha ) / total_reads2   
	                  
	is_matrix = freq_matrix2 / reconst_matrix1
	log2_is_matrix = log2( is_matrix )
	log2_is_matrix[ which(log2_is_matrix==-Inf)] = NA
	log2_is_matrix[ which(log2_is_matrix==Inf)] = NA   
	log2_is_matrix[ which(log2_is_matrix=="NaN")] = NA
	                                                    
	if (show == TRUE){
		heatmap.2(log2_is_matrix, Rowv = F, Colv = F, trace = "none", na.color = "Black",col = colorRampPalette(c("grey","white","blue"))(21), filename=filename, cellwidth = cellwidth, cellheight = cellwidth) 
	}
	return( is_matrix )
}
raw_is = RawInteractionScores(m1,m2,alpha,FALSE)
raw_is[which(raw_is=="NaN")] = 0
raw_is[which(raw_is==Inf)] = 0  
raw_is[which(raw_is==-Inf)] = 0 
pheatmap(log2(raw_is+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth)  
return(log2(raw_is+alpha))
 


#=======================================================










## Supplement Figures   



#=======================================================
# [@S1] Transformation efficiency - Transformation_efficiency.@S1.png
# Data from the table in crg intranet
# For -W condition           
# J46 (small scale): 43200 = 432000 / 20 * 2
# J48 (small scale): 39107.14 = 365000 / 28 * 3
# J51 (small scale): 51525 = 549600/32*3
# J53 (small scale): 104500 = 1254000/36*3
# J56 (small scale): 37872.91667 = 908950/48*2
# J60 (large scale): 808000 = 4848000/12*2


small_scale = c(43200,	39107.14, 51525,	104500,	37872.91667) 
large_scale = c(808000)
png("Transformation_efficiency.@S1.png"); boxplot(small_scale, large_scale, names=c("Small scale","Large scale"), xlab="Method", ylab="# of colonies"); dev.off();  
mean(small_scale) #  [1] 55241.01
large_scale/mean(small_scale) #   [1] 14.62681



# [@S2] PPI - Read comparison - No_selection.@S2.png
# 
# For no selection condition 
# https://gist.githubusercontent.com/stephenturner/3492773/raw/18741cd67e4ffdbf1eb0cc5c032ba661bcd0a045/explore-correlations.r          
no_selection_1 = bm1_1[bait71,prey76] # 2016-12-22_MiSeq; no selection
no_selection_2 = bm3_1[bait71,prey76] # 2016-12-22_MiSeq; no selection; RCA
no_selection_3 = bm1_2[bait71,prey76] # 2017-02-22_MiSeq; no selection 
no_selection_4 = bm1_3[bait71,prey76] # 2017-03-03_MiSeq; no selection 
no_selection_5 = bm1_4[bait71,prey76] # 2017-08-22_MiSeq; no selection; Seaprep  

no_selection_data = data.frame( Exp1N=as.vector(log2(no_selection_1+1)),Exp2N=as.vector(log2(no_selection_2+1)),Exp3N=as.vector(log2(no_selection_3+1)),Exp4N=as.vector(log2(no_selection_4+1)),Exp5N=as.vector(log2(no_selection_5+1)))        
#pairs(no_selection_data,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (No selection media)")   
png("No_selection.@S2.png",width = 600, height = 600); chart.Correlation(no_selection_data); dev.off()
                                                          
                         
# For no selection condition     
# Exp1 = 2017_02_22 -W condition
# Exp2 = 2017_03_03 -W condition
# Exp3 = 2017_08_22 -W condition
# Exp4 = 2017_06_08 -W / Seaprep condition
no_selection_data = data.frame( Exp1=as.vector(log2(bm1_2+1)),Exp2=as.vector(log2(bm1_3+1)),Exp3=as.vector(log2(bm1_4+1)),Exp4=as.vector(log2(bm1_5+1)))        
pdf("No_selection.@S2.V2.pdf",width = 600, height = 600); chart.Correlation(no_selection_data,method="spearman"); dev.off()
png("No_selection.@S2.V2.png",width = 600, height = 600); chart.Correlation(no_selection_data,method="spearman"); dev.off()



# [@S3] PPI - Read comparison - Selection.@S3.png
# For selection condition 
# https://gist.githubusercontent.com/stephenturner/3492773/raw/18741cd67e4ffdbf1eb0cc5c032ba661bcd0a045/explore-correlations.r          
selection_1 = bm2_1[bait71,prey76] # 2016-12-22_MiSeq; selection                    # Exp1A
selection_2 = bm4_1[bait71,prey76] # 2016-12-22_MiSeq; selection; RCA               # Exp2A
selection_3 = bm2_2[bait71,prey76] # 2017-02-22_MiSeq; selection; A                 # Exp3A
selection_4 = bm3_2[bait71,prey76] # 2017-02-22_MiSeq; selection; Q                 # Exp3Q
selection_5 = bm2_3[bait71,prey76] # 2017-03-03_MiSeq; selection; A                 # Exp4A
selection_6 = bm3_3[bait71,prey76] # 2017-03-03_MiSeq; selection; Q                 # Exp4Q
selection_7 = bm2_4[bait71,prey76] # 2017-08-22_MiSeq; selection; Seaprep A4        # Exp5A4
selection_8 = bm3_4[bait71,prey76] # 2017-08-22_MiSeq; selection; Seaprep A8        # Exp5A8
selection_9 = bm4_4[bait71,prey76] # 2017-08-22_MiSeq; selection; Seaprep Q         # Exp5Q

selection_data = data.frame( Exp1A=as.vector(log2(selection_1+1)),Exp2A=as.vector(log2(selection_2+1)),Exp3A=as.vector(log2(selection_3+1)),Exp3Q=as.vector(log2(selection_4+1)),Exp4A=as.vector(log2(selection_5+1)),Exp4Q=as.vector(log2(selection_6+1)),Exp5A4=as.vector(log2(selection_7+1)),Exp5A8=as.vector(log2(selection_8+1)),Exp5Q=as.vector(log2(selection_9+1)))        
pairs(selection_data,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (Selection media)")   
png("Selection.@S3.png",width = 800, height = 800); chart.Correlation(selection_data); dev.off()



# [@S4] PPI - Read comparison - All.@S4.png
# For selection condition 
# https://gist.githubusercontent.com/stephenturner/3492773/raw/18741cd67e4ffdbf1eb0cc5c032ba661bcd0a045/explore-correlations.r

all_data = data.frame( Exp1N=as.vector(log2(no_selection_1+1)),Exp2N=as.vector(log2(no_selection_2+1)),Exp3N=as.vector(log2(no_selection_3+1)),Exp4N=as.vector(log2(no_selection_4+1)),Exp5N=as.vector(log2(no_selection_5+1)),Exp1A=as.vector(log2(selection_1+1)),Exp2A=as.vector(log2(selection_2+1)),Exp3A=as.vector(log2(selection_3+1)),Exp3Q=as.vector(log2(selection_4+1)),Exp4A=as.vector(log2(selection_5+1)),Exp4Q=as.vector(log2(selection_6+1)),Exp5A4=as.vector(log2(selection_7+1)),Exp5A8=as.vector(log2(selection_8+1)),Exp5Q=as.vector(log2(selection_9+1)) ) 
pairs(all_data,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of pair read counts (All)")   
png("All.@S4.png",width = 800, height = 800); chart.Correlation(all_data); dev.off()




#=======================================================


## Figures
#=======================================================
# [@1] Toxicity and Auto-Activation -                    
drawVector <- function( ppi_output_path, alpha,fontsize=9, breaks = 30){
	#ppi_output_path = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S6_R2.blastn.cnt.txt"
	NullMatrix <- read.delim( ppi_output_path, header=FALSE )   
	#hist(log10(NullMatrix[,2]+1),breaks = breaks)
	hist(NullMatrix[,2],breaks = breaks)
	return (NullMatrix)
}
                                          
par(mfrow=c(1,3))   
# Toxicity
#A73e_WRD_R1 = drawVector("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S4_R1.blastn.cnt.txt",0.5, breaks= 10)
A73e_WRD_R2 = drawVector("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S4_R2.blastn.cnt.txt",0.5, breaks= 10)

# Auto-Activation (No Selection)
eB68_WRD_R1 = drawVector("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S6_R1.blastn.cnt.txt",0.5, breaks= 10)
#eB68_WRD_R2 = drawVector("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S6_R2.blastn.cnt.txt",0.5, breaks= 10)

# Auto-Activation (Selection)
eB68_ARD_R1 = drawVector("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S5_R1.blastn.cnt.txt",0.5, breaks= 10)     # eB68_ARD <-> eB68_WRD  # samples were switched
#eB68_ARD_R2 = drawVector("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S5_R2.blastn.cnt.txt",0.5, breaks= 10)     

plot( eB68_WRD_R1$V2, eB68_ARD_R1$V2 ) 
ix = sort( eB68_ARD_R1$V2, index.return = TRUE )$ix
plot( eB68_ARD_R1$V2[ix] / eB68_WRD_R1$V2[ix] )
eB68_ARD_R1$V1[ix]      
                                          
# Auto-activator ranking:
PIAS1 (transcritipn co-regulation; 3842), CNOT7 (transcription complex subunit 7; 1861), CDKN2D (Cyclin-dependent kinase 4 inhibitor D; 335), CCND3 (G1/S-specific cyclin-D3; 41), CDKN1A (Cyclin-dependent kinase inhibitor 1; 33)


# [@2] InteractionScore Correlation -   
# bisR75_1, bisR75_1b, bisR75_2, bisR75_2b, bisR75_3, bisR75_3b, bisR75_4, bisR75_4b, bisR75_4c                 
is_data = data.frame( Exp1A=as.vector(bisR75_1[bait71,prey76]), Exp2A=as.vector(bisR75_1b[bait71,prey76]), Exp3A=as.vector(bisR75_2[bait71,prey76]), Exp3Q=as.vector(bisR75_2b[bait71,prey76]), Exp4A=as.vector(bisR75_3[bait71,prey76]), Exp4Q=as.vector(bisR75_3b[bait71,prey76]), Exp5A4 =as.vector(bisR75_4[bait71,prey76]), Exp5A8=as.vector(bisR75_4b[bait71,prey76]), Exp5Q=as.vector(bisR75_4c[bait71,prey76])  )
pairs(is_data,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of Interaction Scores")   
png("IS.@2.png",width = 800, height = 800); chart.Correlation(is_data); dev.off()                                 

nis_data = data.frame( Exp1A=as.vector(bnisR75_1[bait71,prey76]), Exp2A=as.vector(bnisR75_1b[bait71,prey76]), Exp3A=as.vector(bnisR75_2[bait71,prey76]), Exp3Q=as.vector(bnisR75_2b[bait71,prey76]), Exp4A=as.vector(bnisR75_3[bait71,prey76]), Exp4Q=as.vector(bnisR75_3b[bait71,prey76]), Exp5A4 =as.vector(bnisR75_4[bait71,prey76]), Exp5A8=as.vector(bnisR75_4b[bait71,prey76]), Exp5Q=as.vector(bnisR75_4c[bait71,prey76])  )
pairs(nis_data,lower.panel=panel.smooth, upper.panel=panel.cor,pch=20, main="Correlation of New Interaction Scores")   
png("NIS.@2.png",width = 800, height = 800); chart.Correlation(nis_data); dev.off()                                 

all_nis = bnisR75_1 + bnisR75_1b + bnisR75_2 + bnisR75_2b + bnisR75_3 + bnisR75_3b + bnisR75_4 + bnisR75_4b + bnisR75_4c
simple_plotMatrix( all_nis, 1 )

threshold = 2.0
all_nis2 = (bnisR75_1>threshold) + (bnisR75_1b>threshold) + (bnisR75_2>threshold) + (bnisR75_2b>threshold) + (bnisR75_3>threshold) + (bnisR75_3b>threshold) + (bnisR75_4>threshold) + (bnisR75_4b>threshold) + (bnisR75_4c>threshold)     
simple_plotMatrix( all_nis2[bait71,prey76], 1, 7, filename="IS.@2.pdf" )

                
chart.Correlation( data.frame( Exp2=as.vector(is_jb_2), Exp3=as.vector(is_jb_3), Exp4=as.vector(is_jb_4) ), method="spearman" )    # 0.60, 0.58, 0.66 (spearman)    # 0.73, 0.70, 0.77
chart.Correlation( data.frame( Exp2=as.vector(is_se_2), Exp3=as.vector(is_se_3), Exp4=as.vector(is_se_4) ), method="spearman" )    # 0.53, 0.51, 0.64 (spearman)    # 0.69, 0.65, 0.75
chart.Correlation( data.frame( Exp2=as.vector(is_sl_2), Exp3=as.vector(is_sl_3), Exp4=as.vector(is_sl_4) ), method="spearman" )    # 0.59, 0.57, 0.66 (spearman)    # 0.77, 0.74, 0.81

chart.Correlation( data.frame( Exp2=as.vector(nis_jb_2), Exp3=as.vector(nis_jb_3), Exp4=as.vector(nis_jb_4) ), method="spearman" )    # 0.50, 0.45, 0.58 (spearman)    # 0.86, 0.83, 0.87
chart.Correlation( data.frame( Exp2=as.vector(nis_se_2), Exp3=as.vector(nis_se_3), Exp4=as.vector(nis_se_4) ), method="spearman" )    # 0.46, 0.41, 0.55 (spearman)    # 0.85, 0.82, 0.87
chart.Correlation( data.frame( Exp2=as.vector(nis_sl_2), Exp3=as.vector(nis_sl_3), Exp4=as.vector(nis_sl_4) ), method="spearman" )    # 0.52, 0.47, 0.60 (spearman)    # 0.86, 0.84, 0.88


simple_plotMatrix( bnisR75_4b[bait71,prey76]>threshold, 1, 7 )            
                 

                        
# [@3] Performance   
diag(HIPPIE) = 2
t_j_2 = t.test( nis_j_2[ which( HIPPIE == 1 ) ], nis_j_2[ which( HIPPIE == 0 ) ] )  
t_jb_2 = t.test( nis_jb_2[ which( HIPPIE == 1 ) ], nis_jb_2[ which( HIPPIE == 0 ) ] )  
t_se_2 = t.test( nis_se_2[ which( HIPPIE == 1 ) ], nis_se_2[ which( HIPPIE == 0 ) ] )  
t_sl_2 = t.test( nis_sl_2[ which( HIPPIE == 1 ) ], nis_sl_2[ which( HIPPIE == 0 ) ] )   

t_j_3 = t.test( nis_j_3[ which( HIPPIE == 1 ) ], nis_j_3[ which( HIPPIE == 0 ) ] )  
t_jb_3 = t.test( nis_jb_3[ which( HIPPIE == 1 ) ], nis_jb_3[ which( HIPPIE == 0 ) ] )  
t_se_3 = t.test( nis_se_3[ which( HIPPIE == 1 ) ], nis_se_3[ which( HIPPIE == 0 ) ] )  
t_sl_3 = t.test( nis_sl_3[ which( HIPPIE == 1 ) ], nis_sl_3[ which( HIPPIE == 0 ) ] )   

t_j_4 = t.test( nis_j_4[ which( HIPPIE == 1 ) ], nis_j_4[ which( HIPPIE == 0 ) ] )  
t_jb_4 = t.test( nis_jb_4[ which( HIPPIE == 1 ) ], nis_jb_4[ which( HIPPIE == 0 ) ] )  
t_se_4 = t.test( nis_se_4[ which( HIPPIE == 1 ) ], nis_se_4[ which( HIPPIE == 0 ) ] )  
t_sl_4 = t.test( nis_sl_4[ which( HIPPIE == 1 ) ], nis_sl_4[ which( HIPPIE == 0 ) ] )   

boxplot( c(t_jb_2$statistic,t_jb_3$statistic,t_jb_4$statistic), c(t_se_2$statistic,t_se_3$statistic,t_se_4$statistic), c(t_sl_2$statistic,t_sl_3$statistic,t_sl_4$statistic) )  

diag(HIPPIE) = 2
t_j_2 = t.test( is_j_2[ which( HIPPIE == 1 ) ], is_j_2[ which( HIPPIE == 0 ) ] )  
t_jb_2 = t.test( is_jb_2[ which( HIPPIE == 1 ) ], is_jb_2[ which( HIPPIE == 0 ) ] )  
t_se_2 = t.test( is_se_2[ which( HIPPIE == 1 ) ], is_se_2[ which( HIPPIE == 0 ) ] )  
t_sl_2 = t.test( is_sl_2[ which( HIPPIE == 1 ) ], is_sl_2[ which( HIPPIE == 0 ) ] )   

t_j_3 = t.test( is_j_3[ which( HIPPIE == 1 ) ], is_j_3[ which( HIPPIE == 0 ) ] )  
t_jb_3 = t.test( is_jb_3[ which( HIPPIE == 1 ) ], is_jb_3[ which( HIPPIE == 0 ) ] )  
t_se_3 = t.test( is_se_3[ which( HIPPIE == 1 ) ], is_se_3[ which( HIPPIE == 0 ) ] )  
t_sl_3 = t.test( is_sl_3[ which( HIPPIE == 1 ) ], is_sl_3[ which( HIPPIE == 0 ) ] )   

t_j_4 = t.test( is_j_4[ which( HIPPIE == 1 ) ], is_j_4[ which( HIPPIE == 0 ) ] )  
t_jb_4 = t.test( is_jb_4[ which( HIPPIE == 1 ) ], is_jb_4[ which( HIPPIE == 0 ) ] )  
t_se_4 = t.test( is_se_4[ which( HIPPIE == 1 ) ], is_se_4[ which( HIPPIE == 0 ) ] )  
t_sl_4 = t.test( is_sl_4[ which( HIPPIE == 1 ) ], is_sl_4[ which( HIPPIE == 0 ) ] )   

HIPPIE = BioGrid

boxplot( c(t_jb_2$statistic,t_jb_3$statistic,t_jb_4$statistic), c(t_se_2$statistic,t_se_3$statistic,t_se_4$statistic), c(t_sl_2$statistic,t_sl_3$statistic,t_sl_4$statistic) )  


          




Ones <- function( len ){
	output = 1:len / 1:len
	return( output )
}                   

Zeros <- function( len ){
	output = 1:len * 0
	return (output )
}      

ComparisonRothRef <-function( IS_matrix ){     
	# library(pROC)
	roth=read.csv("~/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/src/roth2016_ref.txt",sep="\t",header=TRUE, comment.char = "#")
	#View(roth)
	roth_pos = roth[which(roth$Interaction.re.annotation == "Interaction"),] 
   	roth_neg = roth[which(roth$Interaction.re.annotation == "No Interaction"),]       
    labels = c( Ones(dim(roth_pos)[1]), Zeros(dim(roth_neg)[1]) )
    
	a = RefineMatrix(IS_matrix)    
	rnames = rownames( a )  
	cnames = colnames( a )     
	print( cbind( rnames[ roth_pos[,6] ], cnames[ roth_pos[,7] ], a[cbind( roth_pos[,6],roth_pos[,7] )] ) )
	
	print( cbind( rnames[ roth_neg[,6] ], cnames[ roth_neg[,7] ], a[cbind( roth_neg[,6],roth_neg[,7] )] ) )
	
	boxplot( a[cbind( roth_pos[,6],roth_pos[,7] )], a[cbind( roth_neg[,6],roth_neg[,7] )] )
	print( t.test( a[cbind( roth_pos[,6],roth_pos[,7] )], a[cbind( roth_neg[,6],roth_neg[,7] )] ) )
	predictions = c( a[cbind( roth_pos[,6],roth_pos[,7] )], a[cbind( roth_neg[,6],roth_neg[,7] )] )     
	
	roc_obj = roc( labels, predictions  )  
	plot(roc_obj)
	print ( auc(roc_obj) )
	wilcox.test( predictions, labels )     
}                          

ComparisonRothRef( bisR75_1 ) 		# t = 5.3204, df = 37.461, p-value = 5.052e-06
ComparisonRothRef( bisR75_1b )
ComparisonRothRef( bisR75_2 )
ComparisonRothRef( bisR75_2b )
ComparisonRothRef( bisR75_3 )
ComparisonRothRef( bisR75_3b )
ComparisonRothRef( bisR75_4 )
ComparisonRothRef( bisR75_4b )
ComparisonRothRef( bisR75_4c )    

bis_all = bisR75_1 + bisR75_1b + bisR75_2 + bisR75_2b + bisR75_3 + bisR75_3b + bisR75_4 + bisR75_4b + bisR75_4c       
simple_plotMatrix( bis_all[bait71,prey76], 1 )      
                              
ComparisonRothRef( bnisR75_1 ) 		# t = 4.4644, df = 34.213, p-value = 8.311e-05
ComparisonRothRef( bnisR75_1b )
ComparisonRothRef( bnisR75_2 )
ComparisonRothRef( bnisR75_2b )
ComparisonRothRef( bnisR75_3 )
ComparisonRothRef( bnisR75_3b )
ComparisonRothRef( bnisR75_4 )
ComparisonRothRef( bnisR75_4b )
ComparisonRothRef( bnisR75_4c )

bnis_all = bnisR75_1 + bnisR75_1b + bnisR75_2 + bnisR75_2b + bnisR75_3 + bnisR75_3b + bnisR75_4 + bnisR75_4b + bnisR75_4c       
simple_plotMatrix( bnis_all[bait71,prey76], 1 )




