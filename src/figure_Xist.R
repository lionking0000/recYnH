######################################################################################################################################  
#
# Version 1.0
# figure_Xist.R
# 
# This for analysis of Xist RNA - protein interaction set
#
######################################################################################################################################

###################################################################################################################################### 
## Initilize functions and load data
source( "/Volumes/users/lserrano/jyang/work/Mireia/src/init.R" )    

######################################################################################################################################
 
#                                 
# .final files only contains bait (RNA) X prey (protein)
#
#X1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/S3_W.rpi.txt",0.5)   # 1341307   
#X1_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/S4_WH.rpi.txt",0.5)  # 1840167
X1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/S3_W.rpi.txt.final",0.5)   # 1341307   
X1_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/S4_WH.rpi.txt.final",0.5)  # 1840167  

eX1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/S3_W.exact.rpi.txt",0.5)   # 1341307     
eX1_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/S4_WH.exact.rpi.txt",0.5)  # 1840167

#X2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S1_W.rpi.txt",0.5)  # 1319368
#X2_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S3_WH.rpi.txt",0.5) 
X2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S1_W.rpi.txt.final",0.5)  # 1319368
X2_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S3_WH.rpi.txt.final",0.5)
eX2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5)  # 1319368
eX2_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S3_WH.exact.rpi.txt",0.5)
 
#X3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S2_W.rpi.txt",0.5)  # 1319860      
#X3_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S4_WH.rpi.txt",0.5) # 2737776      
X3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S2_W.rpi.txt.final",0.5)  # 1319860      
X3_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-14_MiSeq/Blastn/S4_WH.rpi.txt.final",0.5) # 2737776      
          
CorMatrix(X1_W,X2_W)  # 0.6761619
CorMatrix(X1_W,X3_W)  # 0.6336603
CorMatrix(X2_W,X3_W)  # 0.6588208

BasicStats(X1_W)    # 53.1%
BasicStats(X2_W)    # 22.1%
BasicStats(X3_W)    # 19.1%

CorMatrix(X1_WH,X2_WH)  # 0.7291064
CorMatrix(X1_WH,X3_WH)  # 0.4994048
CorMatrix(X2_WH,X3_WH)  # 0.5347829

BasicStats(X1_WH)     # 13.5
BasicStats(X2_WH)     # 12.2
BasicStats(X3_WH)     # 11.6
                                                                                  

#=================================================================================

is2_X1 <<- InteractionScores(X1_W,NormalMixture2( X1_WH ),1,4)      
is2_X2 <<- InteractionScores(X2_W,NormalMixture2( X2_WH ),1,4)      
is2_X3 <<- InteractionScores(X3_W,NormalMixture2( X3_WH ),1,4)      

CorMatrix(is2_X1,is2_X2)  # 0.127369
CorMatrix(is2_X1,is2_X3)  # 0.1143734
CorMatrix(is2_X2,is2_X3)  # 0.09590551


CorMatrixWithoutZeros(is2_X1,is2_X2)  # 0.7331705
CorMatrixWithoutZeros(is2_X1,is2_X3)  # 0.6923892
CorMatrixWithoutZeros(is2_X2,is2_X3)  # 0.6017814

View( (is2_X1>2)+(is2_X2>2)+(is2_X3>2) )

simple_plotMatrix( (is2_X1>2)+(is2_X2>2)+(is2_X3>2),0,4 )  

Overlap = ( ((is2_X1>4)+(is2_X2>4)+(is2_X3>4)) >=2 )   
simple_plotMatrix( Overlap[1:111,], 0, 4 )
simple_plotMatrix( Overlap[112:222,], 0, 4 )
simple_plotMatrix( Overlap[223:335,], 0, 4 )

              ctctttagaata
AGAAGAGTCTCTGGCTCTTTAGAATACTGATCCCATTGAAGATACCACGCTGCAT
              
NR_001564.2      AAGAGTCTCTGGCTCTTTAGAATACTGATCCCAT-------------------------T
Mouse_XIST       --cagactctggctgtttagactacaggatgaatttggagtctgttttgtgctcctgcct
                    ** ******** ****** *** *     **                         *
HG19
#=================================================================================
                              
is3_X1 <<- InteractionScores(X1_W+5,NormalMixture2( X1_WH ),1,4)    
is3_X2 <<- InteractionScores(X2_W+5,NormalMixture2( X2_WH ),1,4)    
is3_X3 <<- InteractionScores(X3_W+5,NormalMixture2( X3_WH ),1,4)             

CorMatrix(is3_X1,is3_X2)  # 0.1656337
CorMatrix(is3_X1,is3_X3)  # 0.1495879
CorMatrix(is3_X2,is3_X3)  # 0.1559151

CorMatrixWithoutZeros(is3_X1,is3_X2)  # 0.6659246
CorMatrixWithoutZeros(is3_X1,is3_X3)  # 0.6260491
CorMatrixWithoutZeros(is3_X2,is3_X3)  # 0.4993193





View( X1_W[,100:120] )
View( X1_WH[,100:120] )
View( X1_WH[,100:120]/(X1_W[,100:120]+1) )
simple_plotMatrix( X1_WH[,100:120]/(X1_W[,100:120]+1),0.5 )
simple_plotMatrix( X2_WH[,100:120]/(X2_W[,100:120]+1),0.5 )
simple_plotMatrix( X3_WH[,100:120]/(X3_W[,100:120]+1),0.5 )     

XIS1 = X1_WH[,100:120]/(X1_W[,100:120]+1)
XIS2 = X2_WH[,100:120]/(X2_W[,100:120]+1)
XIS3 = X3_WH[,100:120]/(X3_W[,100:120]+1)

XIS1 = X1_WH[,100:120]
XIS2 = X2_WH[,100:120]
XIS3 = X3_WH[,100:120]    

CorMatrix(XIS1,XIS2)
CorMatrix(XIS1,XIS3)
CorMatrix(XIS2,XIS3)      
                        
simple_plotMatrix( (XIS1>1000) + (XIS2>1000) + (XIS3>1000), 0, 1 )

#===========================================================                                                         
index = sort( (melt(eX1_WH))[,3], index.return=TRUE, decreasing = TRUE )$ix        

X1_WH_index = sort( (melt(X1_WH))[,3], index.return=TRUE, decreasing = TRUE )$ix        
X1_W_index = sort( (melt(X1_W))[,3], index.return=TRUE, decreasing = TRUE )$ix   
     
GetRowColName(X1_W, X1_W_index[1:40])   
GetRowColName(X1_WH, X1_WH_index[1:40])   

        
#  890  5330  8439 26313 32412
GetRowColName(XIS1, 890)
# [1] "1HQ1:A"
# [1] "X1HQ1.B"

GetRowColName(XIS1, 5330)  
# [1] "1S03:H"
# [1] "X1S03_H.A"

GetRowColName(XIS1, 8439)    
# [1] "2L2K:B"
# [1] "X2L3C.B"

GetRowColName(XIS1, 26313)  
# [1] "HNRNPA1"
# [1] "Xist_12370_12483"

GetRowColName(XIS1, 32412)    
# [1] "IRP"
# [1] "IRE"






CheckRPIControls( X1_W, X1_WH )   
#	[1] 376249
#	[1] 476
#	[1] 790.4391
CheckRPIControls( X2_W, X2_WH )  
#	[1] 31786
#	[1] 90
#	[1] 353.1778
CheckRPIControls( X3_W, X3_WH )    
#	[1] 45749
#	[1] 19
#	[1] 2407.842

XIS1 = X1_WH / ( X1_W + 1 )
XIS1[ire,irp]
XIS2 = X2_WH / ( X2_W + 1 )
XIS2[ire,irp]
XIS3 = X3_WH / ( X3_W + 1 )
XIS3[ire,irp]

simple_plotMatrix(XIS1>300,0)
simple_plotMatrix(XIS2>300,0)
simple_plotMatrix(XIS3>300,0)

simple_plotMatrix( ((XIS1>300) + (XIS2>300) + (XIS3>300))>=3, 0 )      



1EC6:A	2L2K:B	2LBS:B	3NMR:A


1S03:H  
2XNR:A 
1HQ1:A 
1MMS:A  
2L2K:B
2LBS:B
4ED5:A  
3R9W:A 
1EC6:A
2PJP:A  
2ANN:A 
2L3C:A 
3NMR:A
1NYB:A  
1RLG:A 
1RKJ:A


PrintTopPairs <-function(m, range=1:10){
	index = sort( (melt(m))[,3], index.return=TRUE, decreasing = TRUE )$ix   
	return( GetRowColName(m, index[range]) )  
}

f_X1 <<- InteractionScores(X1_W+5,FilterMatrix( X1_WH, 1000 ),1,4) 
f_X2 <<- InteractionScores(X2_W+5,FilterMatrix( X2_WH, 1000 ),1,4) 
f_X3 <<- InteractionScores(X3_W+5,FilterMatrix( X3_WH, 5000 ),1,4) 

f_X1_index = sort( (melt(f_X1))[,3], index.return=TRUE, decreasing = TRUE )$ix        
f_X2_index = sort( (melt(f_X2))[,3], index.return=TRUE, decreasing = TRUE )$ix  
f_X3_index = sort( (melt(f_X3))[,3], index.return=TRUE, decreasing = TRUE )$ix  
                      
CorMatrix(f_X1,f_X2)
CorMatrix(f_X1,f_X3)     
CorMatrixWithoutZeros(f_X1,f_X2)
CorMatrixWithoutZeros(f_X1,f_X3)
                                                 



###################################################################################################################################### 
## EmulsionPCRCheck For Nele

EmulsionPCRCheck <- read.delim("/Volumes/users/lserrano/jyang/work/Nele/src/output/2018-03-12_MiSeq/Blastn/Emulsion_PCR.check.txt", header=FALSE, comment.char="#")
> pcr = EmulsionPCRCheck$V2 / (1189564/4.0) * 10^6
> epcr = EmulsionPCRCheck$V4 / (1116916/4.0) * 10^6
> epcr_check = cbind( pcr, epcr)
> View( epcr_check)
> wilcox.test( pcr, epcr)

wilcox.test( pcr, epcr, paired=TRUE)

plot( density(pcr) )
plot( density(epcr) )





