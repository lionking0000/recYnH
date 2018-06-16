######################################################################################################################################  
#
# Version 1.0
# figure_EGFR.R
# 
# This for analysis of EGFR - protein-protein interaction set
#
######################################################################################################################################

###################################################################################################################################### 
## Initilize functions and load data
source( "/Volumes/users/lserrano/jyang/work/Mireia/src/init.R" )    

######################################################################################################################################
    
#	output/2017-10-09_MiSeq/Blastn/S1_WD.ppi.txt
#	output/2017-10-09_MiSeq/Blastn/S2_WD.ppi.txt
#	output/2017-10-09_MiSeq/Blastn/S3_A8D.ppi.txt
#	output/2017-10-09_MiSeq/Blastn/S4_QD.ppi.txt
                                                                                   
#	../2017-10-09_MiSeq/S1_WD_R1.fastq.gz	4147358 reads
#	../2017-10-09_MiSeq/S2_WD_R1.fastq.gz	2810205 reads 
#	../2017-10-09_MiSeq/S3_A8D_R1.fastq.gz	4077660 reads
#	../2017-10-09_MiSeq/S4_QD_R1.fastq.gz	3857744 reads


E1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S1_WD.ppi.txt",0.5, 4)   				# 1175807 (28.3%) 
E1_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S2_WD.ppi.txt",0.5, 4)   				# 987387  (37.3%)
E1_A8D = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S3_A8D.ppi.txt",0.5, 4)   			# 1522180 (38.7%) 
E1_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S4_QD.ppi.txt",0.5, 4)   				# 1496051 (35.1%)

is2_E1_WH <<- InteractionScores(E1_W,NormalMixture2( E1_WH ),1,4)  
is2_E1_A8D <<- InteractionScores(E1_W,NormalMixture2( E1_A8D ),1,4)
is2_E1_Q <<- InteractionScores(E1_W,NormalMixture2( E1_Q ),1,4)  
                 
rownames(is2_E1_WH)     
KRAS_INDEX = c(95,141)   
KRAS_PARTNER_INDEX = which( colMaxs(is2_E1_WH[KRAS_INDEX,]) > 2.0 )
pheatmap(is2_E2B_WHA1[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE)     

KRAS_INDEX = c(95,141)   
KRAS_PARTNER_INDEX = which( colMaxs(is2_E1_A8D[KRAS_INDEX,]) > 2.0 )
pheatmap(is2_E1_A8D[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE)     

KRAS_INDEX = c(95,141)   
KRAS_PARTNER_INDEX = which( colMaxs(is2_E1_A8D[KRAS_INDEX,]) > 2.0 )
pheatmap(is2_E1_A8D[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE)

rownames(is2_E2_WHA1)     
KRAS_INDEX = c(95,141)   
KRAS_PARTNER_INDEX = which( colMaxs(is2_E1_Q[KRAS_INDEX,]) > 2.0 )
pheatmap(is2_E1_Q[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE, cluster_cols = FALSE)        

pheatmap(is2_E2B_WHA1[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE, cluster_cols = FALSE)  

is2_E1_WH[95,]
                       
is2_E1_ABI1_Bait = rbind(is2_E1_A8D[95,],is2_E1_WH[95,],is2_E1_Q[95,])      
rownames( is2_E1_ABI1_Bait ) = c("ABI1-A8D","ABI1-WH","ABI1-Q")
ABI_PARTNER_INDEX = which( colMaxs(is2_E1_ABI1_Bait) > 2.0 )   
pheatmap(is2_E1_ABI1_Bait[,ABI_PARTNER_INDEX], cluster_rows = FALSE, cluster_cols = FALSE)     
pheatmap(is2_E1_ABI1_Bait[,ABI_PARTNER_INDEX], cluster_rows = FALSE)  
pheatmap(is2_E1_ABI1_Bait[,ABI_PARTNER_INDEX])            

is2_E1_ABI1_Prey = cbind(is2_E1_A8D[,95],is2_E1_WH[,95],is2_E1_Q[,95])     
colnames( is2_E1_ABI1_Prey ) = c("ABI1-A8D","ABI1-WH","ABI1-Q")    
ABI_PARTNER_INDEX2 = which( rowMaxs(is2_E1_ABI1_Prey) > 2.0 )
pheatmap(is2_E1_ABI1_Prey[ABI_PARTNER_INDEX2,], cluster_cols = FALSE)     
pheatmap(is2_E1_ABI1_Prey[ABI_PARTNER_INDEX2,])   

# KRAS as Prey                                                                                  
sort(is2_E1_WH[,141])
#  DMRTB1         GAA      RBFOX2       REEP2      TOPBP1         lam      largeT         p53    178_DOK1      42_FOS     9_BCAR1  109_PLSCR1     124_PXN 
#0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    2.851541    5.483724    5.504017    6.055119    6.524412

# KRAS as Bait
sort(is2_E1_WH[141,])
#    REEP2       TOPBP1          p53          lam       largeT   X41_ERRFI1    X180_RIN2   X159_SPRY2   X96_PIK3CB    X175_WNK1  X201_RASSF5   X133_REPS1 
#0.0000000    0.0000000    0.0000000    0.2205286    1.0315000    1.2421353    1.4727650    1.7880413    2.0764461    2.4176489    2.6939521    2.6941753 
#X187_BRAF X198_RAPGEF4 
#2.9783049    3.1957314
 
# KRAS as Prey                                                                                  
sort(is2_E1_Q[,141])
  DMRTB1         GAA      RBFOX2       REEP2      TOPBP1         lam      largeT         p53    178_DOK1      42_FOS  109_PLSCR1     9_BCAR1     124_PXN 
0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    3.699496    3.936039    5.391909    6.450903    6.830877

# KRAS as Bait
sort(is2_E1_Q[141,])
#      GAA       RBFOX2        REEP2       TOPBP1          lam       largeT          p53      X22_CRK   X41_ERRFI1 X143_SH3BGRL   X99_PIK3R1     X13_CAV2 
# 0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     0.000000     1.046404     1.153740     1.202685     1.416681     1.603941 
#X179_RIN1  X201_RASSF5 
# 1.973320     5.685164
 


is2_E1_ABI1_Bait = rbind(is2_E1_A8D[141,],is2_E1_WH[141,],is2_E1_Q[141,])      
rownames( is2_E1_ABI1_Bait ) = c("KRAS-A8D","KRAS-WH","KRAS-Q")
KRAS_PARTNER_INDEX = which( colMaxs(is2_E1_ABI1_Bait) > 2.0 )   
pheatmap(is2_E1_ABI1_Bait[,KRAS_PARTNER_INDEX], cluster_rows = FALSE, cluster_cols = FALSE)     
pheatmap(is2_E1_ABI1_Bait[,KRAS_PARTNER_INDEX], cluster_rows = FALSE)  
pheatmap(is2_E1_ABI1_Bait[,KRAS_PARTNER_INDEX])            

is2_E1_ABI1_Prey = cbind(is2_E1_A8D[,141],is2_E1_WH[,141],is2_E1_Q[,141])     
colnames( is2_E1_ABI1_Prey ) = c("KRAS-A8D","KRAS-WH","KRAS-Q")    
KRAS_PARTNER_INDEX2 = which( rowMaxs(is2_E1_ABI1_Prey) > 2.0 )
pheatmap(is2_E1_ABI1_Prey[KRAS_PARTNER_INDEX2,], cluster_cols = FALSE)     
pheatmap(is2_E1_ABI1_Prey[KRAS_PARTNER_INDEX2,])




###################################################################################################################################### 
# P25674
# S1: 7557993 total reads (7.5M) = 54.4%
# S2: 6218655 total reads (6.2M) = 44.7%
# S3: 125640 total reads (0.1M) = 0.9% (60.15595 fold less than S1)==> 39450 PPI (31.4%)
# Total: 13902288 reads (13.9M)
#                                 
# .final files only contains bait (protein) X prey (protein)
#  
# E2_W did 12PCR + 12PCR and E2_W2 did 9PCR + 9PCR, so E2_W should have 2^6=64 times more amplification  
# ==> 10PCR + 10PCR could be fine..
# 
#
# Check briefly (with EGFR.fa) ==> need to update fasta file
E2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-03-26_MiSeq/Blastn_All_Ref/S1_W12D.ppi.txt",0.5, 4)   		# 2350895  
E2_WHA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-03-26_MiSeq/Blastn_All_Ref/S2_WHA12D.ppi.txt",0.5, 4)   	# 2050405  
E2_W2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-03-26_MiSeq/Blastn_All_Ref/S3_W9D.ppi.txt",0.5, 4)   		# 39450  (59.59176 fold less)  
                                                                      
 
E2B_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-03-26_MiSeq/Blastn_All_Ref/S1_W12D.barcode.ppi.txt",0.5, 4)   		# 2351542   (31.1% usable reads)
E2B_WHA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-03-26_MiSeq/Blastn_All_Ref/S2_WHA12D.barcode.ppi.txt",0.5, 4)   	# 2050919   (32.9% usable reads)
E2B_W2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-03-26_MiSeq/Blastn_All_Ref/S3_W9D.barcode.ppi.txt",0.5, 4)   		# 39455  (59.59176 fold less)  (31.4% usable reads)


is2_E2_WHA1 <<- InteractionScores(E2_W,NormalMixture2( E2_WHA ),1,4)  
is2_E2_WHA2 <<- InteractionScores(E2_W2,NormalMixture2( E2_WHA ),1,4)  
                                                                        
is2_E2B_WHA1 <<- InteractionScores(E2B_W,NormalMixture2( E2B_WHA ),1,4)  
is2_E2B_WHA2 <<- InteractionScores(E2B_W2,NormalMixture2( E2B_WHA ),1,4)
                                  

#===============================================================================================
# Micro-exon altered bait
#===============================================================================================
                       
sort(is2_E2B_WHA1[95,])		# 1_ABI1 
         p53      X90_NDUFA13        X118_PTK6      X201_RASSF5       X98_PIK3CG        X67_KRT18         ABI1_SK1        X83_MAPK7      X109_PLSCR1 
   0.0000000        0.5932329        1.5164355        1.6200937        1.6621957        2.5055204        2.6355827        2.6483578        3.3600963 
X202_APBB1IP        X174_WASL        X181_RIN3        X66_KRT17         X88_NCK1      X145_SH3GL3 
   3.5428836        4.4779186        4.5249310        4.6414867        5.0828617        5.4015564


sort(is2_E2B_WHA1[171,])	# ABI1_SK1  
#X153_SNRPD2     X139_RPS6KA2         X68_KRT7     X202_APBB1IP      X109_PLSCR1         X88_NCK1 
#   2.103799         2.445119         3.122090         4.349392         5.156770         6.221699

sort(is2_E2B_WHA1[172,])	# ABI1_SK2    
#X18_CEACAM1		X174_WASL      X109_PLSCR1      X145_SH3GL3     X202_APBB1IP        X181_RIN3         X88_NCK1 
#2.575729	    2.765225         3.713297         4.029224         4.604442         4.913382         6.183000

sort(is2_E2B_WHA1[173,])	# ABI1_WT     
#X10_CAMK2A		X92_PAK1		X83_MAPK7        X187_BRAF      X109_PLSCR1        X174_WASL        X181_RIN3      X145_SH3GL3 
#1.470405		1.880597        2.059028         2.431756         3.480401         3.819830         5.088359         5.405397
                                        

ABI1_INDEX = c(95,173,171,172)   
ABI1_PARTNER_INDEX = which( colMaxs(is2_E2B_WHA1[ABI1_INDEX,]) > 2.0 )
pheatmap(is2_E2B_WHA1[ABI1_INDEX,ABI1_PARTNER_INDEX], cluster_rows = FALSE)  
## Save "is2_E2B_WHA1.ABI1_Bait.pdf"

E2B_WHA[ABI1_INDEX,]  
pheatmap(E2B_WHA[ABI1_INDEX,])
pheatmap(E2B_WHA[,ABI1_INDEX])     
pheatmap(E2B_WHA[ABI1_INDEX,])
pheatmap(E2B_WHA[,ABI1_INDEX])     

pheatmap(is2_E2B_WHA1[ABI1_INDEX,ABI1_PARTNER_INDEX])
pheatmap(is2_E2B_WHA1[,ABI1_INDEX])                     

                     

pheatmap(is2_E2B_WHA1[ABI1_INDEX,])
pheatmap(is2_E2B_WHA1[,ABI1_INDEX]) 

plot( E2B_WHA[95,],E2B_WHA[171,] )                                           


ABI1_INDEX = c(95,173,171,172)   
ABI1_PARTNER_INDEX2 = which( rowMaxs(is2_E2B_WHA1[,ABI1_INDEX]) > 2.0 )
pheatmap(is2_E2B_WHA1[ABI1_PARTNER_INDEX2,ABI1_INDEX], cluster_cols = FALSE)       
View( E2B_WHA[ABI1_PARTNER_INDEX2,ABI1_INDEX] )  
## Save "is2_E2B_WHA1.ABI1_Prey.pdf"       

sort(is2_E2B_WHA1[141,])	# 65_KRAS                                       
#   largeT        X179_RIN1        X183_RGL2        X187_BRAF        X182_RGL1      X201_RASSF5 
#0.5160358        1.5843251        1.6891226        2.3156491        2.5513723        4.9445770


sort(is2_E2B_WHA1[180,])	# KRAS_D38A_1  
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0

sort(is2_E2B_WHA1[181,])	# KRAS_D38A_2       
#   REEP2           TOPBP1              lam           largeT              p53        X29_DUSP1 
#0.000000         0.000000         0.000000         0.000000         0.000000         2.230141


sort(is2_E2B_WHA1[182,])	# KRAS_G12D-I21G_1
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0
            
sort(is2_E2B_WHA1[183,])	# KRAS_G12D-I21G_2
#  TOPBP1              lam           largeT              p53     X138_RPS6KA1      X201_RASSF5 
#0.000000         0.000000         0.000000         0.000000         3.935559         5.605779
 
                                   
sort(is2_E2B_WHA1[184,])	# KRAS_G12D-I36F_1
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0

sort(is2_E2B_WHA1[185,])	# KRAS_G12D-I36F_2
RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
     0                0                0                0                0                0


sort(is2_E2B_WHA1[186,])	# KRAS_G12D_1          
#     lam           largeT              p53        X183_RGL2       X132_RBBP7      X201_RASSF5 
#0.000000         0.000000         0.000000         3.119532         3.141956         5.105464

sort(is2_E2B_WHA1[187,])	# KRAS_G12D_2    
#   REEP2           TOPBP1              lam           largeT              p53      X201_RASSF5 
#0.000000         0.000000         0.000000         0.000000         0.000000         4.897635

               

sort(is2_E2B_WHA1[188,])	# KRAS_WT_1
#   REEP2           TOPBP1              lam           largeT              p53      X201_RASSF5 
#0.000000         0.000000         0.000000         0.000000         0.000000         7.231998

sort(is2_E2B_WHA1[189,])	# KRAS_WT_2    
#   REEP2           TOPBP1              lam           largeT              p53        X183_RGL2 
#0.000000         0.000000         0.000000         0.000000         0.000000         3.878032

      
KRAS_INDEX = c(141,188,189,180,181,182,183,184,185,186,187)   
KRAS_PARTNER_INDEX = which( colMaxs(is2_E2B_WHA1[KRAS_INDEX,]) > 2.0 )
pheatmap(is2_E2B_WHA1[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE)    
View(E2B_WHA[KRAS_INDEX,KRAS_PARTNER_INDEX])         
#pheatmap(E2B_WHA[KRAS_INDEX,KRAS_PARTNER_INDEX], cluster_rows = FALSE)         
## Save "is2_E2B_WHA1.KRAS_Bait.pdf"
       
KRAS_INDEX = c(141,188,189,180,181,182,183,184,185,186,187)   
KRAS_PARTNER_INDEX2 = which( rowMaxs(is2_E2B_WHA1[,KRAS_INDEX]) > 2.0 )
pheatmap(is2_E2B_WHA1[KRAS_PARTNER_INDEX2,KRAS_INDEX], cluster_cols = FALSE)  
View(E2B_WHA[KRAS_PARTNER_INDEX2,KRAS_INDEX])       
#pheatmap(E2B_WHA[KRAS_PARTNER_INDEX2,KRAS_INDEX], cluster_cols = FALSE)       
## Save "is2_E2B_WHA1.KRAS_Prey.pdf"   

pheatmap(is2_E2B_WHA1[KRAS_PARTNER_INDEX2,KRAS_PARTNER_INDEX])    




#===============================================================================================
# Micro-exon altered prey
#===============================================================================================
# Appears to function in the signal transduction from Ras activation to actin cytoskeletal remodeling.                                  

sort(is2_E2B_WHA1[,95])		# 1_ABI1     
#       lam           largeT              p53         178_DOK1      202_APBB1IP        117_PTK2B         66_KRT17           42_FOS          50_GRB2 
#  0.000000         0.000000         0.000000         1.102549         3.454692         3.525416         4.596553         4.644379         5.189339 
#109_PLSCR1          124_PXN         174_WASL          88_NCK1           22_CRK          9_BCAR1 
#  5.630853         5.636423         5.837341         6.127684         6.444814         6.635316
                                  



sort(is2_E2B_WHA1[,171])	# ABI1_SK1  
#    RBFOX2           1_ABI1        117_PTK2B         66_KRT17      202_APBB1IP         174_WASL           42_FOS          50_GRB2          9_BCAR1 
# 0.7388753        2.6355827        3.2090263        3.7308942        3.8448107        4.3902802        4.6291291        5.1708952        5.2696092 
#109_PLSCR1          124_PXN           22_CRK          15_CBLB          88_NCK1         38_EPS15 
# 5.3929841        5.6162516        5.7519069        5.8587988        7.2709317        8.0076963

sort(is2_E2B_WHA1[,172])	# ABI1_SK2    
#       lam           largeT              p53            CCL14        117_PTK2B      202_APBB1IP         66_KRT17           42_FOS         174_WASL 
# 0.0000000        0.0000000        0.0000000        0.9065791        3.1118870        3.4719272        3.8802583        4.5490172        5.1000207 
#109_PLSCR1          50_GRB2          124_PXN          9_BCAR1           22_CRK          88_NCK1 
# 5.5432346        5.6571045        6.2777241        6.4426012        6.8960603        7.6531540

sort(is2_E2B_WHA1[,173])	# ABI1_WT     
#       lam           largeT              p53            CEP55           DMRTB1        99_PIK3R1           42_FOS         174_WASL         66_KRT17 
#  0.000000         0.000000         0.000000         2.158618         2.578776         2.888718         4.382267         4.909099         5.358077 
#109_PLSCR1          50_GRB2          124_PXN          9_BCAR1           22_CRK          88_NCK1 
#  5.726413         6.028973         6.187529         6.862809         7.049371         7.502358

                                          



sort(is2_E2B_WHA1[,141])	# 65_KRAS   
#  TOPBP1              lam           largeT              p53         178_DOK1           42_FOS 
#0.000000         0.000000         0.000000         0.000000         3.457911         4.474492


sort(is2_E2B_WHA1[,180])	# KRAS_D38A_1  
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0

sort(is2_E2B_WHA1[,181])	# KRAS_D38A_2       
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0


sort(is2_E2B_WHA1[,182])	# KRAS_G12D-I21G_1
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0
            
sort(is2_E2B_WHA1[,183])	# KRAS_G12D-I21G_2
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0
 
                                   
sort(is2_E2B_WHA1[,184])	# KRAS_G12D-I36F_1
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0

sort(is2_E2B_WHA1[,185])	# KRAS_G12D-I36F_2
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0


sort(is2_E2B_WHA1[,186])	# KRAS_G12D_1          
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0

sort(is2_E2B_WHA1[,187])	# KRAS_G12D_2    
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0
               

sort(is2_E2B_WHA1[,188])	# KRAS_WT_1
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0

sort(is2_E2B_WHA1[,189])	# KRAS_WT_2    
#RBFOX2            REEP2           TOPBP1              lam           largeT              p53 
#     0                0                0                0                0                0





E2_WN = convertNullMatrix(E2_W,1)
E2_W2N = convertNullMatrix(E2_W2,1)
par(mfrow=c(1,1))
CorMatrix(E2_WN, E2_W2N)	# Null Matrix Correlation 0.9488831

E2_WN = convertNullMatrix(E2_W,0.1)
E2_W2N = convertNullMatrix(E2_W2,0.1)
par(mfrow=c(1,1))
CorMatrix(E2_WN, E2_W2N)	# Null Matrix Correlation 0.9898084

E2_WN = convertNullMatrix(E2_W,0.01)
E2_W2N = convertNullMatrix(E2_W2,0.01)
par(mfrow=c(1,1))
CorMatrix(E2_WN, E2_W2N)	# Null Matrix Correlation 0.9912469

E2_WN = convertNullMatrix(E2_W,0.001)
E2_W2N = convertNullMatrix(E2_W2,0.001)
par(mfrow=c(1,1))
CorMatrix(E2_WN, E2_W2N)	# Null Matrix Correlation 0.9913261

is2_E2_WHA1 <<- InteractionScores(E2_W,NormalMixture2( E2_WHA ),0.001,4)  
is2_E2_WHA2 <<- InteractionScores(E2_W2,NormalMixture2( E2_WHA ),0.001,4)   

CorMatrix(is2_E2_WHA1, is2_E2_WHA2)
# [1] 0.9996754



dim(E2_W2)
# [1] 182 182
182*182
# [1] 33124
39450/182/182
# [1] 1.190979  # 1.2x depth

# E2_W detect 19866 pairs 
# sum(E2_W>=1)

# E2_W2 detect 10649 pairs 
# sum(E2_W>=1)

# E2_W only detect 9901 pairs
# sum(E2_W2==0 & E2_W>=1)

# E2_W2 only detect 684 pairs
# sum(E2_W2>=1 & E2_W==0)
                               
# Commonly 9965 pairs detected
# sum(E2_W2>=1 & E2_W>=1)    

# BasicStats(E2_W)  
# [1] "============== m1 =============="
# [1] "BD cnt = 182, 92.9 % detected (169)"
# [1] "AD cnt = 182, 92.9 % detected (169)"
# [1] "ALL cnt = 33124"
# [1] "ZERO cnt = 13258"
# [1] "Detected cnt = 19866"
# [1] "% detected cnt = 60.0"
# [1] "sum = 2350895"
# [1] "max read = 6480"
# [1] "min read = 0"
# [1] "max AD = 240147"
# [1] "min AD = 0"
# [1] "max BD = 127095"
# [1] "max BD = 0"

# BasicStats(E2_W2)
# [1] "BD cnt = 182, 92.3 % detected (168)"
# [1] "AD cnt = 182, 92.3 % detected (168)"
# [1] "ALL cnt = 33124"
# [1] "ZERO cnt = 22475"
# [1] "Detected cnt = 10649"
# [1] "% detected cnt = 32.1"
# [1] "sum = 39450"
# [1] "max read = 94"
# [1] "min read = 0"
# [1] "max AD = 3945" (60.87376 fold less)
# [1] "min AD = 0"
# [1] "max BD = 2075" (61.2506 fold less)
# [1] "max BD = 0"



CorMatrix(E2_W,E2_W2)   
# [1] 0.830946

CorMatrix(is2_E2_WHA1,is2_E2_WHA2)
# [1] 0.9966708

# KRAS Index
which( colnames( E2_W2 ) == "X65_KRAS" )

# KRAS as Prey                                                                                  
sort(is2_E2_WHA1[,141])
#CEP55       CTBP1       DAPK1      DMRTB1         GAA      RBFOX2       REEP2      TOPBP1         lam      largeT         p53    178_DOK1      42_FOS 
#0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    0.000000    3.768731    4.656115
                       
# KRAS as Bait
sort(is2_E2_WHA1[141,])
#       lam       largeT     X92_PAK1   X10_CAMK2A    X11_CASP9    X29_DUSP1    X179_RIN1 X138_RPS6KA1   X113_PRKCG    X187_BRAF    X182_RGL1   X132_RBBP7 
#0.08188503   0.22587940   0.39857214   0.47895685   0.51482800   0.53609339   0.75420426   0.97029509   1.09369096   1.21283621   1.23908619   1.25711201 
# X183_RGL2  X201_RASSF5 
#2.08304924   4.33211332


simple_plotMatrix( (is2_E1_Q>2) + (is2_E2_WHA1>2), 0, 4 )
simple_plotMatrix( (is2_E1_Q>1) + (is2_E2_WHA1>1), 0, 4 )          
                                                                     

###################################################################################################################################### 
# Save Overlap
simple_plotMatrix( (is2_E2_WHA1>2) + (is2_E1_Q>2),0,4,  )
simple_plotMatrix( (is2_E2_WHA1>2) + (is2_E1_Q>2),0,3, filename = "is2_EGFR.Over2Q.pdf" )    
is2_EGFR.Over2Q =  (is2_E2_WHA1>2) + (is2_E1_Q>2)           
#pheatmap( is2_EGFR.Over2Q[1:60,1:60], border_color = "black", fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.1-1.pdf" )   

pheatmap( is2_EGFR.Over2Q[1:60,1:60], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.1-1.pdf" )    
pheatmap( is2_EGFR.Over2Q[1:60,61:120], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.1-2.pdf" )   
pheatmap( is2_EGFR.Over2Q[1:60,121:182], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.1-3.pdf" )    
pheatmap( is2_EGFR.Over2Q[61:120,1:60], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.2-1.pdf" )    
pheatmap( is2_EGFR.Over2Q[61:120,61:120], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.2-2.pdf" )    
pheatmap( is2_EGFR.Over2Q[61:120,121:182], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.2-3.pdf" )    
pheatmap( is2_EGFR.Over2Q[121:182,1:60], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.3-1.pdf" )    
pheatmap( is2_EGFR.Over2Q[121:182,61:120], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.3-2.pdf" )    
pheatmap( is2_EGFR.Over2Q[121:182,121:182], fontsize=6, cluster_rows=FALSE, cluster_cols=FALSE, filename = "is2_EGFR.Over2Q.3-3.pdf" )    


###################################################################################################################################### 
# P25674
# S1: 7557993 total reads (7.5M)
# S2: 6218655 total reads (6.2M)
# S3: 125640 total reads (0.1M) ==> 39450 PPI (31.4%)
# Total: 13902288 reads (13.9M)
#                                 
# .final files only contains bait (protein) X prey (protein)
#  
# Check briefly (with EGFR.fa) ==> need to update fasta file
                                                                           



                          
###################################################################################################################################### 
# KRAS Entry barcode mapping
            
      BARCODE
      ######
GGATCCnnnnnnTTACATAATTACACACTTTGTCTTT.... 
GGATCCatcgat                    
GGATCCGATAGCTTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCACATTGATGGTTCCAAACATCTTCCAT
Query  13   TTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCA  71
            |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
Sbjct  567  TTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCA  509           



>@M03766:102:000000000-BL377:1:1109:12272:27786 1:N:0:1    
GGATCCnnnnnnTTACATAATTACACACTTTGTCTTT.... 
GGATCCGATAGCTTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCACATTGATGGTTCCAAACATCTTCCAT
            TTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCA

>@M03766:102:000000000-BL377:1:1109:10250:28054 1:N:0:1          
GGATCCnnnnnnTTACATAATTACACACTTTGTCTTT.... 
CGTACCAGGAGCTTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCATCTTTTCTTTATGTTTTCGAATTTCT 
            TTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCA


>@M03766:102:000000000-BL377:1:1109:20030:27741 1:N:0:1    
GGATCCnnnnnnTTACATAATTACACACTTTGTCTTT.... 
GGATCCTCTTCCTTACATAATTACACATTATGTCTTTGACTTCTTTTTCTTTTTTTCTCCGTCCTTGCCCATCTTTTCTTTATGTTTTCGAATTTCT
            TTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCGTCTTTGCTCA
            CTACATAATTACACACTTTGTCTTTGACTTCTTTTTCTTCTTTTTACCATCTTTGCTCATCTTTTCTTTATGTTTT

                                                                                                              


# Checking Barcode KRAS count.     
[Read 1]
BARCODE      	NAME      	E2B_W	E2B_WHA	E2B_W2
GGATCCATCGAT 	KRAS_D38A_1 	4182	458	91
GGATCCAGGGTT 	KRAS_D38A_2 	4874	499	66
GGATCCGATAGC 	KRAS_G12D_1 	7479	2534	92
GGATCCGGAGGA 	KRAS_G12D_2 	4246	1499	62
GGATCCGGGGGT 	KRAS_G12D-I21G_1 	556	162	13
GGATCCTCTTCC 	KRAS_G12D-I21G_2 	3159	1423	59
GGATCCCGTTAG 	KRAS_G12D-I36F_1 	4571	60	64
GGATCCTGAAAC 	KRAS_G12D-I36F_2 	3328	183	50
GGATCCGAGTAA 	KRAS_WT_1 	372	315	11
GGATCCAGGGGC 	KRAS_WT_2 	2469	291	52
                                          
[Read 1]  
CCTACATAATTA 	CCTACA 	TAATTA 	12765             
CCTACATAATTA 	CCTACA 	TAATTA 	4572     
CCTACATAATTA 	CCTACA 	TAATTA 	164


[Read 2]                               
CCTACATAATTA 	CCTACA 	TAATTA 	525      
CCTACATAATTA 	CCTACA 	TAATTA 	945 
CCTACATAATTA 	CCTACA 	TAATTA 	5

###################################################################################################################################### 
# ABI1 Entry barcode mapping     - see CheckABI1BarcodedSamples() in Barcode.py   
GGATCC nnnnnn ......
                ATTAATCAGTATAGTGCATGATTGATTCAACATAGTTCCCAGGGAACAGACCAGTC     

              TCATTAATCAGTATAGTGCATGATTGATTCAACATAGTTCCCAGGGAACAGACCAGTCGGGTGTCATGAAGGCATTGTATTTGGG                  # (ABI1 from Violeta's ORF)        
                            HPTGLFPGNYVESIMHYTD-   
                CCTAATCAGTATAGTGCATGATTGATTCAACATAGTTCCCAGGGAACAGACCAGTCACTCGATTGCAGACTCCTTCATACCAGCCATCATCATTCTT    # (ABI1 from Christina's ORF)                  
                KNDDGWYEGVCNRVTGLFPGNYVESIMHYTD-

[Read 1]                                   
BARCODE      	NAME      	E2B_W	E2B_WHA	E2B_W2
GGATCCCTCTCT 	ABI1_WT 	10698 	4858 	160
GGATCCAGGATG 	ABI1_SK1 	7296 	2608 	87
GGATCCTCTATT 	ABI1_SK2 	10060 	4284 	156


[Read 2]
BARCODE      	NAME      	E2B_W	E2B_WHA	E2B_W2
GGATCCCTCTCT 	ABI1_WT 	5567 	14887 	90      
GGATCCAGGATG 	ABI1_SK1 	6128 	10955 	87
GGATCCTCTATT 	ABI1_SK2 	5672 	13343 	104

         
[Read 1]                       
CCTAATCAGTAT 	CCTAAT 	CAGTAT 	15077     # (ABI1 from Christina's ORF) - E2B_W     
CCTAATCAGTAT 	CCTAAT 	CAGTAT 	11334     # (ABI1 from Christina's ORF) - E2B_WHA    
CCTAATCAGTAT 	CCTAAT 	CAGTAT 	265       # (ABI1 from Christina's ORF) - E2B_W2 (56.9 fold less)  

[Read 2]                                
CCTAATCAGTAT 	CCTAAT 	CAGTAT 	3461      # (ABI1 from Christina's ORF) - E2B_W         
CCTAATCAGTAT 	CCTAAT 	CAGTAT 	10529     # (ABI1 from Christina's ORF) - E2B_WHA        
CCTAATCAGTAT 	CCTAAT 	CAGTAT 	69        # (ABI1 from Christina's ORF) - E2B_W2 (50.2 fold less)        
                                           

###################################################################################################################################### 
# coefficient of variation, C.V.
sd(E2B_W) / mean(E2B_W)
# [1] 3.053972

sd(E2B_W2) / mean(E2B_W2)
# [1] 3.110317

sd(E2B_WHA) / mean(E2B_WHA)
# [1] 8.047536


CoV(E1_W)
# [1] 2.257432
                           
CoV(E1_A8D)
# [1] 4.605787

CoV(E1_WH)
# [1] 6.692903

CoV(E1_Q)
# [1] 7.036585


###################################################################################################################################### 
# Name - Index
.
> rownames(is2_E2B_WHA1)
  [1] "101_PIK3R3"       "102_PITPNA"       "105_PLCG2"        "106_PLD1"         "107_PLD2"         "109_PLSCR1"       "10_CAMK2A"       
  [8] "110_PRKAR1A"      "111_PRKCA"        "112_PRKCB1"       "113_PRKCG"        "114_PRKCI"        "115_PRKCZ"        "116_PRKD1"       
 [15] "117_PTK2B"        "118_PTK6"         "119_PTPN11"       "11_CASP9"         "120_PTPN12"       "121_PTPN5"        "122_PTPN6"       
 [22] "123_PTPRR"        "124_PXN"          "125_RAB5A"        "126_RAC1"         "127_RAF1"         "128_RALB"         "129_RALBP1"      
 [29] "12_CAV1"          "130_RALGDS"       "131_RASA1"        "132_RBBP7"        "133_REPS1"        "135_RFXANK"       "136_RGS16"       
 [36] "137_RIPK1"        "138_RPS6KA1"      "139_RPS6KA2"      "13_CAV2"          "140_RPS6KA3"      "142_SH2D3C"       "143_SH3BGRL"     
 [43] "145_SH3GL3"       "146_SH3KBP1"      "147_SHC1"         "148_SHOC2"        "149_SIN3A"        "150_SMAD2"        "151_SMAD3"       
 [50] "152_SNCA"         "153_SNRPD2"       "154_SOCS1"        "155_SOCS3"        "156_SOS1"         "158_SP1"          "159_SPRY2"       
 [57] "15_CBLB"          "160_SRC"          "161_STAT1"        "163_STAT3"        "164_STAT5A"       "165_STAT5B"       "167_TGIF1"       
 [64] "168_TNIP1"        "169_TNK2"         "16_CBLC"          "170_USP6NL"       "171_VAV1"         "174_WASL"         "175_WNK1"        
 [71] "176_YWHAB"        "177_ZNF259"       "178_DOK1"         "179_RIN1"         "17_CDC42"         "180_RIN2"         "181_RIN3"        
 [78] "182_RGL1"         "183_RGL2"         "184_RALA"         "185_DUSP4"        "186_DUSP6"        "187_BRAF"         "189_KSR2"        
 [85] "18_CEACAM1"       "190_PDPK1"        "191_PRKCE"        "192_UBE2L3"       "193_BAD"          "194_LIMK2"        "195_RAPGEF3"     
 [92] "196_RAPGEF1"      "198_RAPGEF4"      "199_RAPGEF5"      "1_ABI1"           "201_RASSF5"       "202_APBB1IP"      "205_SNX1"        
 [99] "207_RGS12"        "208_RGS14"        "209_MYLK"         "211_RASGRP3"      "216_EXOC8"        "217_EXOC2"        "218_RAP2A"       
[106] "21_CREB1"         "221_ARFIP2"       "222_RHOA"         "223_ZFYVE20"      "226_RASSF1"       "22_CRK"           "24_CSK"          
[113] "25_CTNND1"        "27_DNM1"          "28_DOK2"          "29_DUSP1"         "2_AKT1"           "30_EEF1A1"        "32_EGFR"         
[120] "33_ELF3"          "34_ELK1"          "35_ELK4"          "36_EPN1"          "38_EPS15"         "3_AP2A1"          "41_ERRFI1"       
[127] "42_FOS"           "44_GAB1"          "45_GAB2"          "47_GJA1"          "48_GRB10"         "4_APPL1"          "50_GRB2"         
[134] "51_GRB7"          "54_HDAC1"         "55_HIP1"          "56_HIST3H3"       "5_APPL2"          "62_JUN"           "64_KLF11"        
[141] "65_KRAS"          "66_KRT17"         "67_KRT18"         "68_KRT7"          "69_KRT8"          "6_ARAF"           "71_MAP2K2"       
[148] "72_MAP2K3"        "73_MAP2K5"        "76_MAP3K14"       "7_ARF4"           "80_MAPK1"         "81_MAPK14"        "82_MAPK3"        
[155] "83_MAPK7"         "84_MAPK8"         "85_MCF2"          "86_MTA2"          "87_MYC"           "88_NCK1"          "8_ATF1"          
[162] "90_NDUFA13"       "91_NRAS"          "92_PAK1"          "93_PEBP1"         "95_PIK3CA"        "96_PIK3CB"        "98_PIK3CG"       
[169] "99_PIK3R1"        "9_BCAR1"          "ABI1_SK1"         "ABI1_SK2"         "ABI1_WT"          "CCL14"            "CEP55"           
[176] "CTBP1"            "DAPK1"            "DMRTB1"           "GAA"              "KRAS_D38A_1"      "KRAS_D38A_2"      "KRAS_G12D-I21G_1"
[183] "KRAS_G12D-I21G_2" "KRAS_G12D-I36F_1" "KRAS_G12D-I36F_2" "KRAS_G12D_1"      "KRAS_G12D_2"      "KRAS_WT_1"        "KRAS_WT_2"       
[190] "RBFOX2"           "REEP2"            "TOPBP1"           "lam"              "largeT"           "p53"

                                                                               




.
GGATCCGAGTAA 	315
GGATCCATCGAT 	458
GGATCCAGGGTT 	499
GGTTCCGATAGC 	17
GGATCCTGAAAC 	183
CGATCCGATAGC 	17
CCGATAGCTTAC 	20
GGTTCCGGAGGA 	19
GGATCCAGGGGC 	291
GGTTCCTCTTCC 	12
CGATCCTCTTCC 	11
GGATCCTCTTCC 	
GGATCCGATAGC 	
GGATCCGGAGGA 	
GGATCCCGTTAG 	60
GGATCCGCTAGC 	18
GGATCCGGGGGT 	162


GGATCCGAGTAA 	11
GGATCCAGGGGC 	52
GGATCCGATAGC 	92
GGATCCATCGAT 	91
GGATCCGGAGGA 	62
GGATCCAGGGTT 	66
GGATCCTCTTCC 	59
GGATCCCGTTAG 	64
GGATCCTGAAAC 	50
GGATCCGGGGGT 	13

KRAS WT clone2	GCCCCT	AGGGGC	2469	291
KRAS WT clone1	TTACTC	GAGTAA	372	315
KRAS G12D-I36F clone2	GTTTCA	TGAAAC	3328	183
KRAS G12D-I36F clone1	CTAACG	CGTTAG	4571	60
KRAS G12D-I21G clone2	GGAAGA	TCTTCC	3159	1423
KRAS G12D-I21G clone1	ACCCCC	GGGGGT	556
KRAS G12D clone2	TCCTCC	GGAGGA	4246	1499
KRAS G12D clone1	GCTATC	GATAGC	7479	2534
KRAS D38A clone2	AACCCT	AGGGTT	4874
KRAS D38A clone1	ATCGAT	ATCGAT	4182




GGATCC GAGTAA




###################################################################################################################################### 
# STRING - ABI1 - APBB1IP

https://string-db.org/cgi/network.pl?taskId=QUPWs5soAU8f

#node1	node2	node1_string_internal_id	node2_string_internal_id	node1_external_id	node2_external_id	neighborhood_on_chromosome	gene_fusion	phylogenetic_cooccurrence	homology	coexpression	experimentally_determined_interaction	database_annotated	automated_textmining	combined_score
CYFIP2	NCKAP1L	1850773	1847410	9606.ENSP00000325817	9606.ENSP00000293373	0	0	0	0	0.571	0.967	0.900	0.642	0.999
NCKAP1	CYFIP2	1854096	1850773	9606.ENSP00000354251	9606.ENSP00000325817	0	0	0	0	0.549	0.967	0.900	0.883	0.999
CYFIP1	NCKAP1L	1850651	1847410	9606.ENSP00000324549	9606.ENSP00000293373	0	0	0	0	0.549	0.967	0.900	0.678	0.999
NCKAP1	CYFIP1	1854096	1850651	9606.ENSP00000354251	9606.ENSP00000324549	0	0	0	0	0.549	0.967	0.900	0.811	0.999
NCKAP1	WASF1	1854096	1853845	9606.ENSP00000354251	9606.ENSP00000352425	0	0	0	0	0.100	0.921	0.900	0.763	0.998
WASF2	NCKAP1	1859921	1854096	9606.ENSP00000396211	9606.ENSP00000354251	0	0	0	0	0.072	0.912	0.900	0.802	0.998
WASF1	CYFIP2	1853845	1850773	9606.ENSP00000352425	9606.ENSP00000325817	0	0	0	0	0.154	0.921	0.900	0.780	0.998
SOS1	GRB2	1858968	1852164	9606.ENSP00000384675	9606.ENSP00000339007	0	0	0	0	0.053	0.825	0.900	0.931	0.998
WASF2	CYFIP2	1859921	1850773	9606.ENSP00000396211	9606.ENSP00000325817	0	0	0	0	0.102	0.912	0.900	0.800	0.998
WASF2	CYFIP1	1859921	1850651	9606.ENSP00000396211	9606.ENSP00000324549	0	0	0	0	0.102	0.912	0.900	0.730	0.997
WASF3	CYFIP2	1851801	1850773	9606.ENSP00000335055	9606.ENSP00000325817	0	0	0	0	0.154	0.913	0.900	0.660	0.997
WASF1	CYFIP1	1853845	1850651	9606.ENSP00000352425	9606.ENSP00000324549	0	0	0	0	0.154	0.921	0.900	0.526	0.996
BRK1	CYFIP1	1861910	1850651	9606.ENSP00000432472	9606.ENSP00000324549	0	0	0	0	0	0.905	0.900	0.693	0.996
WIPF1	WASL	1853898	1842952	9606.ENSP00000352802	9606.ENSP00000223023	0	0	0	0	0.074	0.685	0.900	0.880	0.996
BRK1	NCKAP1	1861910	1854096	9606.ENSP00000432472	9606.ENSP00000354251	0	0	0	0	0	0.825	0.900	0.836	0.996
WASF2	BAIAP2	1859921	1849864	9606.ENSP00000396211	9606.ENSP00000316338	0	0	0	0	0.042	0.774	0.900	0.798	0.995
WASF3	CYFIP1	1851801	1850651	9606.ENSP00000335055	9606.ENSP00000324549	0	0	0	0	0.154	0.858	0.900	0.641	0.995
NCKAP1	WASF3	1854096	1851801	9606.ENSP00000354251	9606.ENSP00000335055	0	0	0	0	0.127	0.861	0.900	0.651	0.995
WASF3	NCKAP1L	1851801	1847410	9606.ENSP00000335055	9606.ENSP00000293373	0	0	0	0	0.100	0.861	0.900	0.594	0.994
SOS1	HRAS	1858968	1849235	9606.ENSP00000384675	9606.ENSP00000309845	0	0	0	0	0.053	0.884	0.900	0.583	0.994
CDC42	WASL	1849707	1842952	9606.ENSP00000314458	9606.ENSP00000223023	0	0	0	0	0	0.375	0.900	0.915	0.994
WASF1	WASL	1853845	1842952	9606.ENSP00000352425	9606.ENSP00000223023	0	0	0	0	0	0.774	0.900	0.768	0.994
WASF2	ABI1	1859921	1856444	9606.ENSP00000396211	9606.ENSP00000365312	0	0	0	0	0	0.774	0.900	0.757	0.994
WASF2	WASL	1859921	1842952	9606.ENSP00000396211	9606.ENSP00000223023	0	0	0	0	0	0.774	0.900	0.751	0.993
ACTR2	WASL	1856793	1842952	9606.ENSP00000367220	9606.ENSP00000223023	0	0	0	0	0	0.769	0.900	0.744	0.993
BRK1	CYFIP2	1861910	1850773	9606.ENSP00000432472	9606.ENSP00000325817	0	0	0	0	0.049	0.608	0.900	0.843	0.993
SOS1	EPS8	1858968	1846713	9606.ENSP00000384675	9606.ENSP00000281172	0	0	0	0	0	0.403	0.900	0.883	0.992
ABI1	NCKAP1	1856444	1854096	9606.ENSP00000365312	9606.ENSP00000354251	0	0	0	0	0	0.774	0.900	0.713	0.992
BRK1	WASF1	1861910	1853845	9606.ENSP00000432472	9606.ENSP00000352425	0	0	0	0	0	0.779	0.900	0.675	0.992
ABI1	EPS8	1856444	1846713	9606.ENSP00000365312	9606.ENSP00000281172	0	0	0	0	0	0.360	0.900	0.884	0.991
TRIP10	WASL	1850251	1842952	9606.ENSP00000320493	9606.ENSP00000223023	0	0	0	0	0.053	0.774	0.900	0.658	0.991
SOS1	KRAS	1858968	1844404	9606.ENSP00000384675	9606.ENSP00000256078	0	0	0	0	0.104	0.448	0.900	0.839	0.991
SOS1	ABI1	1858968	1856444	9606.ENSP00000384675	9606.ENSP00000365312	0	0	0	0	0.072	0.360	0.900	0.861	0.990
CTTN	WASL	1856519	1842952	9606.ENSP00000365745	9606.ENSP00000223023	0	0	0	0	0.104	0.399	0.900	0.850	0.990
WASF1	NCKAP1L	1853845	1847410	9606.ENSP00000352425	9606.ENSP00000293373	0	0	0	0	0.100	0.667	0.900	0.701	0.989
BRK1	NCKAP1L	1861910	1847410	9606.ENSP00000432472	9606.ENSP00000293373	0	0	0	0	0	0.744	0.900	0.629	0.989
ABI1	CYFIP1	1856444	1850651	9606.ENSP00000365312	9606.ENSP00000324549	0	0	0	0	0	0.774	0.900	0.522	0.988
ABI1	BAIAP2	1856444	1849864	9606.ENSP00000365312	9606.ENSP00000316338	0	0	0	0	0	0.774	0.900	0.480	0.987
BAIAP2	EPS8	1849864	1846713	9606.ENSP00000316338	9606.ENSP00000281172	0	0	0	0	0	0.622	0.900	0.695	0.987
GRB2	HRAS	1852164	1849235	9606.ENSP00000339007	9606.ENSP00000309845	0	0	0	0	0.100	0.379	0.900	0.792	0.986
BRK1	WASF2	1861910	1859921	9606.ENSP00000432472	9606.ENSP00000396211	0	0	0	0	0	0.605	0.900	0.677	0.986
ABI1	WASL	1856444	1842952	9606.ENSP00000365312	9606.ENSP00000223023	0	0	0	0	0	0.774	0.900	0.441	0.986
ACTR2	CTTN	1856793	1856519	9606.ENSP00000367220	9606.ENSP00000365745	0	0	0	0	0.072	0.740	0.900	0.472	0.985
ABI1	CYFIP2	1856444	1850773	9606.ENSP00000365312	9606.ENSP00000325817	0	0	0	0	0	0.774	0.900	0.417	0.985
ABI1	WASF1	1856444	1853845	9606.ENSP00000365312	9606.ENSP00000352425	0	0	0	0	0	0.774	0.900	0.408	0.985
BRK1	WASF3	1861910	1851801	9606.ENSP00000432472	9606.ENSP00000335055	0	0	0	0	0	0.678	0.900	0.555	0.984
CYFIP2	BAIAP2	1850773	1849864	9606.ENSP00000325817	9606.ENSP00000316338	0	0	0	0	0.053	0.774	0.900	0.342	0.984
NCKAP1	BAIAP2	1854096	1849864	9606.ENSP00000354251	9606.ENSP00000316338	0	0	0	0	0.049	0.774	0.900	0.350	0.984
CYFIP1	BAIAP2	1850651	1849864	9606.ENSP00000324549	9606.ENSP00000316338	0	0	0	0	0.046	0.774	0.900	0.312	0.983
ABI1	NCKAP1L	1856444	1847410	9606.ENSP00000365312	9606.ENSP00000293373	0	0	0	0	0	0.774	0.900	0.286	0.982
TRIP10	CDC42	1850251	1849707	9606.ENSP00000320493	9606.ENSP00000314458	0	0	0	0	0.053	0.377	0.900	0.741	0.982
WASF2	NCKAP1L	1859921	1847410	9606.ENSP00000396211	9606.ENSP00000293373	0	0	0	0	0.094	0.629	0.900	0.515	0.981
BAIAP2	CDC42	1849864	1849707	9606.ENSP00000316338	9606.ENSP00000314458	0	0	0	0	0	0.360	0.900	0.732	0.981
GRB2	KRAS	1852164	1844404	9606.ENSP00000339007	9606.ENSP00000256078	0	0	0	0	0.100	0.379	0.900	0.678	0.979
WASF3	WASL	1851801	1842952	9606.ENSP00000335055	9606.ENSP00000223023	0	0	0	0	0.067	0.417	0.900	0.654	0.978
ACTR2	WIPF1	1856793	1853898	9606.ENSP00000367220	9606.ENSP00000352802	0	0	0	0	0.106	0.162	0.900	0.729	0.977
ACTR2	CDC42	1856793	1849707	9606.ENSP00000367220	9606.ENSP00000314458	0	0	0	0	0.088	0.100	0.900	0.678	0.970
CTTN	CDC42	1856519	1849707	9606.ENSP00000365745	9606.ENSP00000314458	0	0	0	0	0.067	0	0.900	0.686	0.968
WASF1	CDC42	1853845	1849707	9606.ENSP00000352425	9606.ENSP00000314458	0	0	0	0	0	0.380	0.900	0.537	0.968
BAIAP2	WASL	1849864	1842952	9606.ENSP00000316338	9606.ENSP00000223023	0	0	0	0	0.064	0.153	0.900	0.630	0.966
WASF2	CDC42	1859921	1849707	9606.ENSP00000396211	9606.ENSP00000314458	0	0	0	0	0	0.360	0.900	0.510	0.965
WASF1	BAIAP2	1853845	1849864	9606.ENSP00000352425	9606.ENSP00000316338	0	0	0	0	0.047	0.329	0.900	0.491	0.963
SOS1	CDC42	1858968	1849707	9606.ENSP00000384675	9606.ENSP00000314458	0	0	0	0	0.054	0.174	0.900	0.564	0.961
WASF2	ACTR2	1859921	1856793	9606.ENSP00000396211	9606.ENSP00000367220	0	0	0	0	0	0.285	0.900	0.497	0.960
CDC42	EPS8	1849707	1846713	9606.ENSP00000314458	9606.ENSP00000281172	0	0	0	0	0	0	0.900	0.618	0.960
ACTR2	WASF1	1856793	1853845	9606.ENSP00000367220	9606.ENSP00000352425	0	0	0	0	0.047	0.360	0.900	0.409	0.959
WIPF1	CDC42	1853898	1849707	9606.ENSP00000352802	9606.ENSP00000314458	0	0	0	0	0.049	0.305	0.900	0.449	0.958
BRK1	ACTR2	1861910	1856793	9606.ENSP00000432472	9606.ENSP00000367220	0	0	0	0	0.049	0.097	0.900	0.565	0.957
ACTR2	CYFIP2	1856793	1850773	9606.ENSP00000367220	9606.ENSP00000325817	0	0	0	0	0.094	0	0.900	0.544	0.955
BRK1	ABI1	1861910	1856444	9606.ENSP00000432472	9606.ENSP00000365312	0	0	0	0	0	0.146	0.900	0.493	0.952
ACTR2	CYFIP1	1856793	1850651	9606.ENSP00000367220	9606.ENSP00000324549	0	0	0	0	0.064	0	0.900	0.526	0.951
CYFIP2	WASL	1850773	1842952	9606.ENSP00000325817	9606.ENSP00000223023	0	0	0	0	0	0	0.900	0.526	0.950
CYFIP2	CDC42	1850773	1849707	9606.ENSP00000325817	9606.ENSP00000314458	0	0	0	0	0.088	0.095	0.900	0.461	0.949
ACTR2	WASF3	1856793	1851801	9606.ENSP00000367220	9606.ENSP00000335055	0	0	0	0	0.047	0.093	0.900	0.482	0.949
NCKAP1	CDC42	1854096	1849707	9606.ENSP00000354251	9606.ENSP00000314458	0	0	0	0	0.088	0.277	0.900	0.316	0.948
NCKAP1	WASL	1854096	1842952	9606.ENSP00000354251	9606.ENSP00000223023	0	0	0	0	0.181	0	0.900	0.409	0.947
BRK1	WASL	1861910	1842952	9606.ENSP00000432472	9606.ENSP00000223023	0	0	0	0	0	0	0.900	0.496	0.947
WASF3	CDC42	1851801	1849707	9606.ENSP00000335055	9606.ENSP00000314458	0	0	0	0	0.052	0.072	0.900	0.465	0.946
ACTR2	NCKAP1	1856793	1854096	9606.ENSP00000367220	9606.ENSP00000354251	0	0	0	0	0.052	0.095	0.900	0.434	0.944
ACTR2	NCKAP1L	1856793	1847410	9606.ENSP00000367220	9606.ENSP00000293373	0	0	0	0	0.140	0.095	0.900	0.370	0.944
CYFIP1	CDC42	1850651	1849707	9606.ENSP00000324549	9606.ENSP00000314458	0	0	0	0	0.088	0.095	0.900	0.396	0.943
CYFIP1	WASL	1850651	1842952	9606.ENSP00000324549	9606.ENSP00000223023	0	0	0	0	0	0	0.900	0.431	0.940
BRK1	BAIAP2	1861910	1849864	9606.ENSP00000432472	9606.ENSP00000316338	0	0	0	0	0	0	0.900	0.408	0.938
CDC42	NCKAP1L	1849707	1847410	9606.ENSP00000314458	9606.ENSP00000293373	0	0	0	0	0.088	0.277	0.900	0.166	0.937
BRK1	CDC42	1861910	1849707	9606.ENSP00000432472	9606.ENSP00000314458	0	0	0	0	0.054	0	0.900	0.368	0.935
ABI1	CDC42	1856444	1849707	9606.ENSP00000365312	9606.ENSP00000314458	0	0	0	0	0.054	0	0.900	0.364	0.934
WASF3	BAIAP2	1851801	1849864	9606.ENSP00000335055	9606.ENSP00000316338	0	0	0	0	0.042	0	0.900	0.365	0.933
WIPF1	WASF1	1853898	1853845	9606.ENSP00000352802	9606.ENSP00000352425	0	0	0	0	0	0	0.900	0.355	0.932
CYFIP2	CYFIP1	1850773	1850651	9606.ENSP00000325817	9606.ENSP00000324549	0	0	0	0.988	0	0.336	0.900	0.596	0.931
ABI1	WASF3	1856444	1851801	9606.ENSP00000365312	9606.ENSP00000335055	0	0	0	0	0	0	0.900	0.304	0.927
ACTR2	ABI1	1856793	1856444	9606.ENSP00000367220	9606.ENSP00000365312	0	0	0	0	0.108	0	0.900	0.241	0.926
WASF2	WIPF1	1859921	1853898	9606.ENSP00000396211	9606.ENSP00000352802	0	0	0	0	0.049	0	0.900	0.277	0.925
ACTR2	BAIAP2	1856793	1849864	9606.ENSP00000367220	9606.ENSP00000316338	0	0	0	0	0	0	0.900	0.275	0.924
CTTN	TRIP10	1856519	1850251	9606.ENSP00000365745	9606.ENSP00000320493	0	0	0	0	0.052	0	0.900	0.258	0.923
WIPF1	WASF3	1853898	1851801	9606.ENSP00000352802	9606.ENSP00000335055	0	0	0	0	0	0	0.900	0.213	0.918
BRK1	WIPF1	1861910	1853898	9606.ENSP00000432472	9606.ENSP00000352802	0	0	0	0	0	0	0.900	0.193	0.915
NCKAP1L	WASL	1847410	1842952	9606.ENSP00000293373	9606.ENSP00000223023	0	0	0	0	0.056	0	0.900	0.177	0.915
WASF1	WASF3	1853845	1851801	9606.ENSP00000352425	9606.ENSP00000335055	0	0	0	0.853	0.063	0	0.900	0.930	0.915
WASF2	WASF3	1859921	1851801	9606.ENSP00000396211	9606.ENSP00000335055	0	0	0	0.847	0	0	0.900	0.930	0.914
WIPF1	BAIAP2	1853898	1849864	9606.ENSP00000352802	9606.ENSP00000316338	0	0	0	0	0	0	0.900	0.170	0.913
WIPF1	CYFIP2	1853898	1850773	9606.ENSP00000352802	9606.ENSP00000325817	0	0	0	0	0.110	0	0.900	0.094	0.912
WASF2	WASF1	1859921	1853845	9606.ENSP00000396211	9606.ENSP00000352425	0	0	0	0.878	0	0	0.900	0.955	0.911
SOS1	APBB1IP	1858968	1856462	9606.ENSP00000384675	9606.ENSP00000365411	0	0	0	0	0.091	0	0.900	0.088	0.910
WIPF1	NCKAP1L	1853898	1847410	9606.ENSP00000352802	9606.ENSP00000293373	0	0	0	0	0.114	0	0.900	0.046	0.908
APBB1IP	HRAS	1856462	1849235	9606.ENSP00000365411	9606.ENSP00000309845	0	0	0	0	0	0.114	0.900	0.050	0.908
APBB1IP	KRAS	1856462	1844404	9606.ENSP00000365411	9606.ENSP00000256078	0	0	0	0	0	0.114	0.900	0	0.907
ABI1	WIPF1	1856444	1853898	9606.ENSP00000365312	9606.ENSP00000352802	0	0	0	0	0	0	0.900	0.090	0.905
ACTR2	TRIP10	1856793	1850251	9606.ENSP00000367220	9606.ENSP00000320493	0	0	0	0	0.054	0	0.900	0.079	0.905
BAIAP2	NCKAP1L	1849864	1847410	9606.ENSP00000316338	9606.ENSP00000293373	0	0	0	0	0	0	0.900	0.082	0.904
HRAS	KRAS	1849235	1844404	9606.ENSP00000309845	9606.ENSP00000256078	0	0	0.527	0.983	0	0	0.900	0.736	0.902
NCKAP1	NCKAP1L	1854096	1847410	9606.ENSP00000354251	9606.ENSP00000293373	0	0	0	0.971	0	0	0.900	0.684	0.901
NCKAP1	WIPF1	1854096	1853898	9606.ENSP00000354251	9606.ENSP00000352802	0	0	0	0	0	0	0.900	0.058	0.901
WIPF1	CYFIP1	1853898	1850651	9606.ENSP00000352802	9606.ENSP00000324549	0	0	0	0	0	0	0.900	0.041	0.900
WASF1	TRIP10	1853845	1850251	9606.ENSP00000352425	9606.ENSP00000320493	0	0	0	0	0.053	0.780	0	0.245	0.829
GRB2	WASL	1852164	1842952	9606.ENSP00000339007	9606.ENSP00000223023	0	0	0	0	0	0.360	0	0.683	0.789
CTTN	WIPF1	1856519	1853898	9606.ENSP00000365745	9606.ENSP00000352802	0	0	0	0	0	0.407	0	0.651	0.785
GRB2	CDC42	1852164	1849707	9606.ENSP00000339007	9606.ENSP00000314458	0	0	0	0	0.055	0.140	0	0.662	0.702
CTTN	GRB2	1856519	1852164	9606.ENSP00000365745	9606.ENSP00000339007	0	0	0	0	0.057	0.360	0	0.466	0.650
CTTN	ABI1	1856519	1856444	9606.ENSP00000365745	9606.ENSP00000365312	0	0	0	0	0	0.503	0	0.290	0.632
WIPF1	GRB2	1853898	1852164	9606.ENSP00000352802	9606.ENSP00000339007	0	0	0	0	0.068	0.238	0	0.489	0.606
WASF2	CTTN	1859921	1856519	9606.ENSP00000396211	9606.ENSP00000365745	0	0	0	0	0	0	0	0.470	0.470
CYFIP1	TRIP10	1850651	1850251	9606.ENSP00000324549	9606.ENSP00000320493	0	0	0	0	0	0.417	0	0.086	0.444
EPS8	WASL	1846713	1842952	9606.ENSP00000281172	9606.ENSP00000223023	0	0	0	0	0.054	0	0	0.432	0.439
CTTN	BAIAP2	1856519	1849864	9606.ENSP00000365745	9606.ENSP00000316338	0	0	0	0	0	0	0	0.437	0.437
ACTR2	KRAS	1856793	1844404	9606.ENSP00000367220	9606.ENSP00000256078	0	0	0	0	0.220	0.191	0	0.177	0.435
WASF2	GRB2	1859921	1852164	9606.ENSP00000396211	9606.ENSP00000339007	0	0	0	0	0	0.285	0	0.234	0.428
ACTR2	GRB2	1856793	1852164	9606.ENSP00000367220	9606.ENSP00000339007	0	0	0	0	0.090	0.094	0	0.354	0.422
CTTN	WASF1	1856519	1853845	9606.ENSP00000365745	9606.ENSP00000352425	0	0	0	0	0	0	0	0.413	0.413
WASF1	GRB2	1853845	1852164	9606.ENSP00000352425	9606.ENSP00000339007	0	0	0	0	0	0.285	0	0.204	0.406























