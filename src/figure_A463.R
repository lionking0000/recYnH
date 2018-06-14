######################################################################################################################################  
#
# Version 1.0
# figure_A463.R
# 
# This for analysis of A463 set
#
######################################################################################################################################

RemoveAA_A463 <- function(m1){          
	# Need to check again!
	
	#protein        reads
	#'cask'        205930
	#'fam103a1' 53325
	#'grip1' 185311
	#'kif24' 369656
	#'krt31' 363843
	#'map4' 111397
	#'mkrn1' 40985
    
	#Is krr1 in AD very important to remove in AD?
	
	RowIndex = 1:dim(m1w)[1]   
	ColIndex = 1:dim(m1w)[2]   
	
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="CASK" ) )
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="FAM103A1" ) )
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="GRIP1" ) )
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="KIF24" ) )
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="KRT31" ) )
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="MAP4" ) )
	RowIndex = setdiff( RowIndex, which( rownames(m1q)=="MKRN1" ) )
    
	ColIndex = setdiff( ColIndex, which( colnames(m1q)=="KRR1" ) )                   
	
	return( m1[RowIndex,ColIndex] )      # 454 x 460  ==> 456 x 462 (expected)
}

###################################################################################################################################### 
## Initilize functions and load data
source( "/Volumes/users/lserrano/jyang/work/Mireia/src/init.R" )

######################################################################################################################################

## A463 ##  
# 2017-11-17_MiSeq; 461 x 461 RNA binding proteins  
m1w = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-17_MiSeq/S1_W.ppi.txt",0.5,2) 
m1q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-17_MiSeq/S2_Q.ppi.txt",0.5,2)    
m1q_nm2 = NormalMixture2(m1q)
is1 = InteractionScores(m1w,m1q,1.0,2) 
is1_nm2 = InteractionScores(m1w,m1q_nm2,1.0,2)   

# 2017-12-07_MiSeq; 461 x 461 RNA binding proteins (After auto-activator removal)     
m2w = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-12-07_MiSeq/S1_W.ppi.txt",0.5,2) 
m2q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-12-07_MiSeq/S2_Q.ppi.txt",0.5,2)    
m2q_nm2 = NormalMixture2(m2q)
is2 = InteractionScores(m2w,m2q,1.0,2) 
is2_nm2 = InteractionScores(m2w,m2q_nm2,1.0,2) 

# 2018-02-21_MiSeq; 461 x 461 RNA binding proteins (After auto-activator removal)        
m3w = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-02-21_MiSeq/S1_W.ppi.txt",0.5,2) 
m3q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2018-02-21_MiSeq/S2_Q.ppi.txt",0.5,2)     
m3q_nm2 = NormalMixture2(m3q)    
is3 = InteractionScores(m3w,m3q,1.0,2)                            
is3_nm2 = InteractionScores(m3w,m3q_nm2,1.0,2)           


# Remove Auto-activator   
aa_m1w = RemoveAA_A463(m1w)          
aa_m2w = RemoveAA_A463(m2w)             
aa_m3w = RemoveAA_A463(m3w)
aa_m1q = RemoveAA_A463(m1q)          
aa_m2q = RemoveAA_A463(m2q)             
aa_m3q = RemoveAA_A463(m3q)
             
aa_is1 = RemoveAA_A463(is1)
aa_is2 = RemoveAA_A463(is2)
aa_is3 = RemoveAA_A463(is3)     

CorMatrixWithoutZeros(aa_is1, aa_is2)	# 0.854107  
CorMatrixWithoutZeros(aa_is1, aa_is3)	# 0.8452934     
CorMatrixWithoutZeros(aa_is2, aa_is3)	# 0.8730909      

aa_is1_nm2 = RemoveAA_A463(is1_nm2)
aa_is2_nm2 = RemoveAA_A463(is2_nm2)
aa_is3_nm2 = RemoveAA_A463(is3_nm2)
                 
CorMatrixWithoutZeros(aa_is1_nm2, aa_is2_nm2)	# 0.8583654  
CorMatrixWithoutZeros(aa_is1_nm2, aa_is3_nm2)	# 0.8537419    
CorMatrixWithoutZeros(aa_is2_nm2, aa_is3_nm2)	# 0.8830641       




#=======================================================================
# Exp-1 vs Exp-2
#=======================================================================

# Auto-activator
plot( (rowSums(m1q)/sum(m1q))/(rowSums(m1w)/sum(m1w)) )
plot( (colSums(m1q)/sum(m1q))/(colSums(m1w)/sum(m1w)) )
           
CorMatrix( m1w, m2w )    		# 0.5475902
CorMatrix( m1q, m2q )    		# 0.2961716
CorMatrix( m1q_nm2, m2q_nm2 )   # 0.2961637
CorMatrix( is1, is2 )    		# 0.4280532
CorMatrix( is1_nm2, is2_nm2 )   # 0.421677
                                                    
CorMatrixWithoutZeros( m1w, m2w )   		# 0.5076568
CorMatrixWithoutZeros( m1q, m2q )   		# 0.5116938
CorMatrixWithoutZeros( m1q_nm2, m2q_nm2 )   # 0.5081825
CorMatrixWithoutZeros( is1, is2 )   		# 0.8520726
CorMatrixWithoutZeros( is1_nm2, is2_nm2 )   # 0.8289255

CorMatrixWithoutZeros( is1, is2, cutoff=1.0 )   		# 0.7208608
CorMatrixWithoutZeros( is1_nm2, is2_nm2, cutoff=1.0 )   # 0.7534212

CorMatrixWithoutZeros( is1, is2, cutoff=5.0 )   		# 0.5432292
CorMatrixWithoutZeros( is1_nm2, is2_nm2, cutoff=5.0 )   # 0.5449886
                                                                       
CorMatrixWithoutZeros( is1, is2, cutoff=10.0 )   		# 0.36464
CorMatrixWithoutZeros( is1_nm2, is2_nm2, cutoff=10.0 )   # 0.4717686


## Change INT_CUTOFF 2 ~ 7               
INT_CUTOFF = 3.0
iCount1 = (is1_nm2>=INT_CUTOFF)
iCount2 = (is2_nm2>=INT_CUTOFF)
iCount = iCount1 + iCount2

sum(iCount1)	#3:2754
sum(iCount2)    #3:2563
x = hist(iCount)
x 
pheatmap(iCount,fontsize = 1)  # save as iCount_cutoff_3.pdf etc...
                                 
# 208091  3543     887

#=======================================================================
# Exp-2 vs Exp-3
#=======================================================================
      
# Auto-activator
plot( (rowSums(m2q)/sum(m2q))/(rowSums(m2w)/sum(m2w)) )
plot( (colSums(m2q)/sum(m2q))/(colSums(m2w)/sum(m2w)) )

CorMatrix( m2w, m3w )    		# 0.7592146
CorMatrix( m2q, m3q )    		# 0.8321926
CorMatrix( m2q_nm2, m3q_nm2 )   # 0.8321494
CorMatrix( is2, is3 )    		# 0.5569398
CorMatrix( is2_nm2, is3_nm2 )   # 0.5514135
                                                    
CorMatrixWithoutZeros( m2w, m3w )   		# 0.7161943
CorMatrixWithoutZeros( m2q, m3q )   		# 0.8249402
CorMatrixWithoutZeros( m2q_nm2, m3q_nm2 )   # 0.8012406
CorMatrixWithoutZeros( is2, is3 )   		# 0.8778643
CorMatrixWithoutZeros( is2_nm2, is3_nm2 )   # 0.8820979   

CorMatrixWithoutZeros( is2, is3, cutoff=1.0 )   		# 0.8230502
CorMatrixWithoutZeros( is2_nm2, is3_nm2, cutoff=1.0 )   # 0.8732925

CorMatrixWithoutZeros( is2, is3, cutoff=5.0 )   		# 0.7747705
CorMatrixWithoutZeros( is2_nm2, is3_nm2, cutoff=5.0 )   # 0.8386607
                                                                       
CorMatrixWithoutZeros( is2, is3, cutoff=10.0 )   		# 0.489933
CorMatrixWithoutZeros( is2_nm2, is3_nm2, cutoff=10.0 )   # 0.7970224

## Change INT_CUTOFF 2 ~ 7               
INT_CUTOFF = 7.0

iCount2 = (is2_nm2>=INT_CUTOFF)
iCount3 = (is3_nm2>=INT_CUTOFF)
iCount = iCount2 + iCount3

sum(iCount2)
sum(iCount3)
x = hist(iCount)
x 
#simple_plotMatrix(iCount,0.0,2)   
pheatmap(iCount,fontsize = 1)  # save as iCount_cutoff_3.pdf etc...

# cutoff = 2, 3425, 1800, 208439, 2939, 1143 
# cutoff = 3, 2563, 1380, 209535, 2029, 957 
# cutoff = 4, 1856, 1036, 210358, 1434, 729   
# cutoff = 5, 1321, 696, 211029, 967, 525   
# cutoff = 6, 842, 444, 211537, 682, 302
# cutoff = 7, 529, 230, 211939, 405, 177 

overlap_proportion = c(1143/(2939+1143), 957/(2029+957), 729/(1436+729), 525/(967+525), 302/(682+302), 177/(405+177))
plot( overlap_proportion )      # ==> Detected_Proportion_with_Cutoff.pdf

min_proportion = c(1143/(1800), 957/(1380), 729/(1036), 525/(696), 302/(444), 177/(230))
plot( min_proportion )      # ==> Detected_Proportion_with_Cutoff.pdf                



#=======================================================================
# Exp-1 vs Exp-3
#=======================================================================
      
# Auto-activator
plot( (rowSums(m3q)/sum(m3q))/(rowSums(m3w)/sum(m3w)) )
plot( (colSums(m3q)/sum(m3q))/(colSums(m3w)/sum(m3w)) )
                                                                     


CorMatrix( m1w, m3w )    		# 0.6015055
CorMatrix( m1q, m3q )    		# 0.256507
CorMatrix( m1q_nm2, m3q_nm2 )   # 0.2564931
CorMatrix( is1, is3 )    		# 0.416474
CorMatrix( is1_nm2, is3_nm2 )   # 0.3914876
                                                    
CorMatrixWithoutZeros( m1w, m3w )   		# 0.5631287
CorMatrixWithoutZeros( m1q, m3q )   		# 0.7187181
CorMatrixWithoutZeros( m1q_nm2, m3q_nm2 )   # 0.7299034
CorMatrixWithoutZeros( is1, is3 )   		# 0.8498881
CorMatrixWithoutZeros( is1_nm2, is3_nm2 )   # 0.8519958

CorMatrixWithoutZeros( is1, is3, cutoff=1.0 )   		# 0.7637179
CorMatrixWithoutZeros( is1_nm2, is3_nm2, cutoff=1.0 )   # 0.8108409

CorMatrixWithoutZeros( is1, is3, cutoff=5.0 )   		# 0.601232
CorMatrixWithoutZeros( is1_nm2, is3_nm2, cutoff=5.0 )   # 0.6489626
                                                                       
CorMatrixWithoutZeros( is1, is3, cutoff=10.0 )   		# 0.285764
CorMatrixWithoutZeros( is1_nm2, is3_nm2, cutoff=10.0 )   # 0.2397119      


## Change INT_CUTOFF 2 ~ 7               
INT_CUTOFF = 3.0
iCount1 = (is1_nm2>=INT_CUTOFF)
iCount3 = (is3_nm2>=INT_CUTOFF)
iCount = iCount1 + iCount3

sum(iCount1)	#3:2754
sum(iCount3)    #3:1380
x = hist(iCount)
x 
pheatmap(iCount,fontsize = 1)  # save as iCount_cutoff_3.pdf etc...
                                 
# 209010  2888     623

#=======================================================================
# Exp-1 & Exp-2 & Exp-3
#=======================================================================

# union analysis
## Change INT_CUTOFF 2 ~ 7               
INT_CUTOFF = 3.0

iCount1 = (is1_nm2>=INT_CUTOFF)
iCount2 = (is2_nm2>=INT_CUTOFF)
iCount3 = (is3_nm2>=INT_CUTOFF)
iCount = iCount1 + iCount2 + iCount3

sum(iCount1)	#3:2754
sum(iCount2)    #3:2563
sum(iCount3)    #3:1380
x = hist(iCount)
x 
pheatmap(iCount,fontsize = 1)  # save as iCount_cutoff_3.pdf etc...
# 207703	# 3527	# 703	# 588
# 3527+703+588 ==> 4818
