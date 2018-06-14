#================================================================================================================================================  
#
#	Functions for Fig3
#            
#================================================================================================================================================  

Fig3D_ACC2 <- function(average, KnownPPI, UseSymmetryScore=FALSE, filename = FALSE){
    acc2_average = AutoActivationCorrection2( average )	
    if (UseSymmetryScore){
        acc2_average = SymmetricSum(acc2_average,0)
    }   
    melted_acc2_average = melt(acc2_average)
    if (filename != FALSE){
        write.table(melted_acc2_average,file=filename,sep="\t")
    }

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    for( i in 0:100 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrix(KnownPPI),RefineMatrix(acc2_average), i/10) ) }
    f1_list[is.na(f1_list)] = 0.0

    plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - acc2_average"  )
    plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2, main="f1 - acc2_average"  )
    plot( 0:100 /10, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - acc2_average" )
    plot( 0:100 /10, precision_list, col = "green", pch = 15, cex = 2, main="precision - acc2_average"  )
    
    plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2, main="acc - acc2_average"  )
    plot( 0:100 /10, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - acc2_average"  )

    print( roc(KnownPPI, acc2_average) )
    print( max(f1_list) )
    is_cutoff = (1:100/10)[order(f1_list, decreasing = TRUE)][1]
    print ( is_cutoff )
    print ( sum(RefineMatrix(acc2_average)>is_cutoff) )
    View( cbind(0:100 /10, f1_list) )

    return( acc2_average )
}


Fig3D_ACC3 <- function(average, KnownPPI, UseSymmetryScore=FALSE, filename = FALSE){
    acc3_average = AutoActivationCorrection3( average )	
    if (UseSymmetryScore){
        acc3_average = SymmetricSum(acc3_average,0)
    }   
    melted_acc3_average = melt(acc3_average)
    if (filename != FALSE){
        write.table(melted_acc3_average,file=filename,sep="\t")
    }

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    for( i in 0:100 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrix(KnownPPI),RefineMatrix(acc3_average), i/10) ) }
    f1_list[is.na(f1_list)] = 0.0

    plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - acc3_average"  )
    plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2, main="f1 - acc3_average"  )
    plot( 0:100 /10, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - acc3_average" )
    plot( 0:100 /10, precision_list, col = "green", pch = 15, cex = 2, main="precision - acc3_average"  )
    
    plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2, main="acc - acc3_average"  )
    plot( 0:100 /10, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - acc3_average"  )

    print( roc(KnownPPI, acc3_average) )
    print( max(f1_list) )
    is_cutoff = (1:100/10)[order(f1_list,decreasing = TRUE)][1]
    print ( is_cutoff )
    print ( sum(RefineMatrix(acc3_average)>is_cutoff) )
    View( cbind(0:100 /10, f1_list) )

    return( acc3_average )
}


GetACC <- function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    negatives = labels[which(predictions<cutoff)]
    TP = sum(positives)
    FN = sum(negatives)
    FP = length(positives) - TP
    TN = length(negatives) - FN
    return ( ( TP + TN ) / ( TP + TN + FP + FN ) )
}


GetF1Score <- function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    negatives = labels[which(predictions<cutoff)]
    TP = sum(positives)
    FN = sum(negatives)
    FP = length(positives) - TP
    TN = length(negatives) - FN
    return ( ( 2 * TP ) / ( 2* TP + FP + FN ) )
}

GetMCC <- function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    negatives = labels[which(predictions<cutoff)]
    TP = sum(positives)
    FN = sum(negatives)
    FP = length(positives) - TP
    TN = length(negatives) - FN
    return ( ( TP * TN - FP * FN ) / ( sqrt(TP+FP) * sqrt(TP+FN) * sqrt(TN+FP) * sqrt(TN+FN) ) )
}

GetPrecision <-function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    return ( sum(positives)/length(positives) ) # TP / (TP+FP)
}	

#GetRecall <-function( labels, predictions, cutoff ){
#    positives = labels[which(predictions>=cutoff)]
#    return ( sum(positives)/length(positives) ) # TP / (TP+FP)
#}	

GetSensitivity <-function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    tot_positives = sum(labels)
    return ( sum(positives)/tot_positives ) # TP / P = TP / (TP+FN)
}	

GetSpecificity <-function( labels, predictions, cutoff ){
    negatives = labels[which(predictions<cutoff)]
	FN = sum(negatives)
	TN = length(negatives) - FN
	tot_positives = sum(labels)
    tot_negatives = length(labels) - tot_positives
    return ( TN /tot_negatives ) # TN / N = TN / (TN+FP)
}	

GetFDR <- function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    TP = sum(positives)
    FP = length(positives) - TP
	return( FP / ( FP + TP ) ) 
}

GetInformedness <- function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    negatives = labels[which(predictions<cutoff)]
    TP = sum(positives)
    FN = sum(negatives)
    FP = length(positives) - TP
    TN = length(negatives) - FN
	TPR = TP / ( TP + FN )
	TNR = TN / ( TN + FP )
    return ( TPR + TNR - 1 )
}

GetInformation <- function( labels, predictions, cutoff ){
    positives = labels[which(predictions>=cutoff)]
    negatives = labels[which(predictions<cutoff)]
    TP = sum(positives)
    FN = sum(negatives)
    FP = length(positives) - TP
    TN = length(negatives) - FN
    return (TP, TN, FP, FN )
}

            
#================================================================================================================================================    







# 192 Spot test
if (TRUE){  
	#SpotTestMGj54_JS_R75 <- read_excel("~/Dropbox (CRG)/Collaboration/Mireia/Spot_test/MGj54/SpotTestMGj54_JS.xlsx", sheet = "R75_Final_List") 
	SpotTest_192 <- read_excel("~/Dropbox (CRG)/matrixppi/figures/spot_test_temp.xlsx", sheet = "192_spot_test")   
	
	pos_weak_index = SpotTest_192$`SD/-W/AbA1` == 1
	pos_strong_index = SpotTest_192$`SD/-W/AbA1/-H/-A` == 1
	
	pos_weak_TRUE = SpotTest_192[which(pos_weak_index==TRUE),5:6]    
	pos_weak_FALSE = SpotTest_192[which(pos_weak_index==FALSE),5:6]    
	pos_weak_ALL = SpotTest_192[which(is.na(pos_weak_index)==FALSE),5:6]  
	                                                     
	bm1_A[ ((pos_weak_TRUE[,2]-1)*78+pos_weak_TRUE[,1])[,1] ]    
	> bm1_A[ ((pos_weak_TRUE[,2]-1)*78+pos_weak_TRUE[,1])[,1] ]         
	 [1]    22   341   135 22535  1357   213   822  3230   294  1045   175 11762  2349  2304     0  2671   476   918  1159  4553  1570   144
	[23]   333  2211   731   970  5814  8678   233 30446    71  2756  2483   194   182    19  2075   390   620   264  2239  1554  3806    59
	[45]   129     0  6642  1631   952    28   108   241  1512  2904  1387   348  5571  1438  3395     0  4662   388    81   744   309   809
	[67]     9     5  4209    32    30   276   351    18   440    27  3560  1490  2135  1237   452  2925   161  4708    54    84   516   670
	[89]     5  2144   158    16   257     6
	is2_bm2_A[ ((pos_weak_TRUE[,2]-1)*78+pos_weak_TRUE[,1])[,1] ]    
	  
	> bm1_A[ ((pos_weak_FALSE[,2]-1)*78+pos_weak_FALSE[,1])[,1] ]  
	 [1]    0    0   90    0   24   55    3    0    0    4    2  170   51    4   29  616   41    0    1 1443    0    0 1136  895    8  814
	[27]  120  770   26   68   73   20    8    0    7    2   69   34   10   80   52   30  317    3   13   22  146   10   78   14   24    2
	[53]   26    0   20 2433   16    0    2 1670    0    9    0   20    0   68   20    9   44    1    1    0   20   47    0    0  240    0
	[79]   43    0   33    0    0    0   57    7   36    0 2453    0   45    0   16    4    0   16  256   47
	is2_bm2_A[ ((pos_weak_FALSE[,2]-1)*78+pos_weak_FALSE[,1])[,1] ]           
	
	
	HIPPIE[ ((pos_weak_FALSE[,2]-1)*78+pos_weak_FALSE[,1])[,1] ] 
	
	BioGrid[ ((pos_weak_FALSE[,2]-1)*78+pos_weak_FALSE[,1])[,1] ] 
	
	is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ]           
	is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ]   
	order( is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ] )




### 192 Spot test
#========================================================================             
# Firgure3A - Final
#
# is2_order ( cut = 3)
#
#========================================================================
    ## Fig 3a
pdf("192spot.is2_overlap2.pdf")
par(mfrow=c(4,1)) 
is2_order = order( is2_overlap2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ], decreasing = TRUE )
barplot( is2_overlap2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
barplot( pos_weak_index[is2_order] )        
barplot( HIPPIE[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
barplot( BioGrid[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] )
dev.off()


### 192 Spot test
#========================================================================             
# Firgure3C - Final
#
# AUC graphs
#
#========================================================================

	
	par(mfrow=c(4,1))
	is2_order = order( is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ], decreasing = TRUE )
	barplot( is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
	barplot( pos_weak_index[is2_order] )        
	barplot( HIPPIE[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
	barplot( BioGrid[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] )
	
	boxplot( bm1_A[ ((pos_weak_TRUE[,2]-1)*78+pos_weak_TRUE[,1])[,1] ], bm1_A[ ((pos_weak_FALSE[,2]-1)*78+pos_weak_FALSE[,1])[,1] ])


	par(mfrow=c(4,1)) 
	is2_order = order( is2_bm7_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ], decreasing = TRUE )
	barplot( is2_bm7_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
	barplot( pos_weak_index[is2_order] )        
	barplot( HIPPIE[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
	barplot( BioGrid[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] )

    pos_weak_ALL_info = cbind( pos_weak_ALL, pos_weak_index[1:192] )
    View( pos_weak_ALL_info[ is2_order, ] )
 
	is2_bm1_A = InteractionScores(bm1_W,NormalMixture2(bm1_A),1,9 )	        # auc = 0.8077
	is2_bm1_RA = InteractionScores(bm1_W,NormalMixture2(bm1_RA),1,9 )	    # auc = 0.8077
	is2_bm2_A = InteractionScores(bm2_W,NormalMixture2(bm2_A),1,9 )			# auc = 0.8077	  
	is2_bm2_Q = InteractionScores(bm2_W,NormalMixture2(bm2_Q),1,9 )	        # auc = 0.7959
	is2_bm3_A = InteractionScores(bm3_W,NormalMixture2(bm3_A),1,9 )	        # auc = 0.8077
	is2_bm3_Q = InteractionScores(bm3_W,NormalMixture2(bm3_Q),1,9 )	        # auc = 0.8077
	is2_bm3_QL = InteractionScores(bm3_W,NormalMixture2(bm3_QL),1,9 )	    # auc = 0.8077
	is2_bm4_S4A = InteractionScores(bm4_SW,NormalMixture2(bm4_S4A),1,9 )	# auc = 0.8077  
	is2_bm4_S8A = InteractionScores(bm4_SW,NormalMixture2(bm4_S8A),1,9 )	# auc = 0.8077  
	is2_bm4_SQ = InteractionScores(bm4_SW,NormalMixture2(bm4_SQ),1,9 )	    # auc = 0.8077
	is2_bm6_SQ = InteractionScores(bm6_SW,NormalMixture2(bm6_SQ),1,9 )	    # auc = 0.8077
	is2_bm7_SQ = InteractionScores(bm7_SW,NormalMixture2(bm7_SQ),1,9 )      # auc = 0.8077



	pdf("192spot.is2_overlap2.auc.pdf")     

	v = 3
	is2_overlap2 = (is2_bm1_A>v)+(is2_bm1_RA>v)+(is2_bm2_A>v)+(is2_bm2_Q>v)+(is2_bm3_A>v)+(is2_bm3_Q>v)+(is2_bm3_QL>v)+(is2_bm4_S4A>v)+(is2_bm4_S8A>v)+(is2_bm4_SQ>v)+(is2_bm6_SQ>v)+(is2_bm7_SQ>v)   
    
	par(mfrow=c(1,1)) 
	labels = pos_weak_index[is2_order]
	predictions = is2_overlap2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] 
	roc_obj = roc( labels, predictions  )
	if( TRUE ){plot(roc_obj) }   


    predictions = is2_bm1_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_1 = roc( labels, predictions  )  # 0.782
	predictions = is2_bm1_RA[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_2 = roc( labels, predictions  )   # 0.8239
	predictions = is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_3 = roc( labels, predictions  )    # 0.8077
	predictions = is2_bm2_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_4 =  roc( labels, predictions  )    # 0.7959
	predictions = is2_bm3_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_5 = roc( labels, predictions  )    # 0.8571
	predictions = is2_bm3_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_6 =  roc( labels, predictions  )    # 0.8035
	predictions = is2_bm3_QL[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_7 =  roc( labels, predictions  )    # 0.7787
	predictions = is2_bm4_S4A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_8 = roc( labels, predictions  )   # 0.8472
	predictions = is2_bm4_S8A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_9 = roc( labels, predictions  )   # 0.7682
	predictions = is2_bm4_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_10 = roc( labels, predictions  )   # 0.7765
	predictions = is2_bm6_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_11 = roc( labels, predictions  )  # 0.7993
	predictions = is2_bm7_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_12 = roc( labels, predictions  ) # 0.8061

    lines(roc_1$specificities,roc_1$sensitivities)
    lines(roc_2$specificities,roc_2$sensitivities)
    lines(roc_3$specificities,roc_3$sensitivities)
    lines(roc_4$specificities,roc_4$sensitivities)
    lines(roc_5$specificities,roc_5$sensitivities)
    lines(roc_6$specificities,roc_6$sensitivities)
    lines(roc_7$specificities,roc_7$sensitivities)
    lines(roc_8$specificities,roc_8$sensitivities)
    lines(roc_9$specificities,roc_9$sensitivities)
    lines(roc_10$specificities,roc_10$sensitivities)
    lines(roc_11$specificities,roc_11$sensitivities)
    lines(roc_12$specificities,roc_12$sensitivities)

	dev.off()

	     
	#pdf("192spot.is2_bm2_A.auc.pdf")          
	#par(mfrow=c(1,1)) 
	#labels = pos_weak_index[is2_order]
	#predictions = is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] 
	#roc_obj = roc( labels, predictions  )
	#if( TRUE ){plot(roc_obj) }   
	#dev.off()   


	#boxplot( c(0.782,0.8239,0.8077,0.7959,0.8571,0.8035,0.7787,0.8472,0.7682,0.7765,0.7993,0.8061) )  
	pdf("192spot.is2.auc.box.pdf")      
	boxplot( c(0.782,0.8239,0.8077,0.7959,0.8571,0.8035,0.7787,0.8472,0.7682,0.7765,0.7993,0.8061),ylim=c(0.7,0.9) )
	dev.off()








### 192 Spot test
#========================================================================             
# Firgure3D - Final
#
# AUC graphs
#
#========================================================================


predictions_1 = is2_bm1_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_1 = roc( labels, predictions_1  )  # 0.782
predictions_2 = is2_bm1_RA[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_2 = roc( labels, predictions_2  )   # 0.8239
predictions_3 = is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_3 = roc( labels, predictions_3  )    # 0.8077
predictions_4 = is2_bm2_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_4 =  roc( labels, predictions_4  )    # 0.7959
predictions_5 = is2_bm3_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_5 = roc( labels, predictions_5  )    # 0.8571
predictions_6 = is2_bm3_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_6 =  roc( labels, predictions_6  )    # 0.8035
predictions_7 = is2_bm3_QL[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_7 =  roc( labels, predictions_7  )    # 0.7787
predictions_8 = is2_bm4_S4A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_8 = roc( labels, predictions_8  )   # 0.8472
predictions_9 = is2_bm4_S8A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_9 = roc( labels, predictions_9  )   # 0.7682
predictions_10 = is2_bm4_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_10 = roc( labels, predictions_10  )   # 0.7765
predictions_11 = is2_bm6_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_11 = roc( labels, predictions_11  )  # 0.7993
predictions_12 = is2_bm7_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_12 = roc( labels, predictions_12  ) # 0.8061

pdf("figure3D.pdf")
# single exp
precision_list = c()
sensitivity_list = c()
specificity_list = c()
for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(labels,predictions_2, i/2) ) }
for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,predictions_2, i/2) ) }
for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(labels,predictions_2, i/2) ) }
#plot( 0:20 /2, precision_list )
plot( 0:20 /2, sensitivity_list, col = "red", pch = 15, cex = 4  )
lines( 0:20 /2, sensitivity_list, col = "red" )
#points( 0:20 /2, sensitivity_list, col = "red" )
points( 0:20 /2, specificity_list, col = "green", pch = 15, cex = 4 )
lines( 0:20 /2, specificity_list, col = "green" )

# double exp
SUM_MASK = 1:length(predictions_1) * 0 
SUM_MASK[which(predictions_1>0.0)] = 1
SUM_MASK[which(predictions_2>0.0)] = SUM_MASK[which(predictions_2>0.0)] + 1
average = predictions_1 + predictions_2
average = average / SUM_MASK
average[is.na(average)] = 0.0

precision_list = c()
sensitivity_list = c()
specificity_list = c()
for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(labels,average, i/2) ) }
for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,average, i/2) ) }
for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(labels,average, i/2) ) }
#plot( 0:20 /2, precision_list )
points( 0:20 /2, sensitivity_list, col = "red", pch = 16, cex = 4  )
lines( 0:20 /2, sensitivity_list, col = "red"  )
#points( 0:20 /2, sensitivity_list, col = "red" )
points( 0:20 /2, specificity_list, col = "green", pch = 16, cex = 4 )
lines( 0:20 /2, specificity_list, col = "green"  )


# 4 exp
SUM_MASK = 1:length(predictions_1) * 0 
SUM_MASK[which(predictions_1>0.0)] = 1
SUM_MASK[which(predictions_2>0.0)] = SUM_MASK[which(predictions_2>0.0)] + 1
SUM_MASK[which(predictions_3>0.0)] = SUM_MASK[which(predictions_3>0.0)] + 1
SUM_MASK[which(predictions_4>0.0)] = SUM_MASK[which(predictions_4>0.0)] + 1
average = predictions_1 + predictions_2 + predictions_3 + predictions_4
average = average / SUM_MASK
average[is.na(average)] = 0.0

precision_list = c()
sensitivity_list = c()
specificity_list = c()
for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(labels,average, i/2) ) }
for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,average, i/2) ) }
for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(labels,average, i/2) ) }
#plot( 0:20 /2, precision_list )
points( 0:20 /2, sensitivity_list, col = "red", pch = 17, cex = 4  )
lines( 0:20 /2, sensitivity_list, col = "red"  )

#points( 0:20 /2, sensitivity_list, col = "red" )
points( 0:20 /2, specificity_list, col = "green", pch = 17, cex = 4 )
lines( 0:20 /2, specificity_list, col = "green"  )


# 8 exp
SUM_MASK = 1:length(predictions_1) * 0 
SUM_MASK[which(predictions_1>0.0)] = 1
SUM_MASK[which(predictions_2>0.0)] = SUM_MASK[which(predictions_2>0.0)] + 1
SUM_MASK[which(predictions_3>0.0)] = SUM_MASK[which(predictions_3>0.0)] + 1
SUM_MASK[which(predictions_4>0.0)] = SUM_MASK[which(predictions_4>0.0)] + 1
SUM_MASK[which(predictions_5>0.0)] = SUM_MASK[which(predictions_5>0.0)] + 1
SUM_MASK[which(predictions_6>0.0)] = SUM_MASK[which(predictions_6>0.0)] + 1
SUM_MASK[which(predictions_7>0.0)] = SUM_MASK[which(predictions_7>0.0)] + 1
SUM_MASK[which(predictions_8>0.0)] = SUM_MASK[which(predictions_8>0.0)] + 1
average = predictions_1 + predictions_2 + predictions_3 + predictions_4 + predictions_5 + predictions_6 + predictions_7 + predictions_8
average = average / SUM_MASK
average[is.na(average)] = 0.0

precision_list = c()
sensitivity_list = c()
specificity_list = c()
for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(labels,average, i/2) ) }
for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,average, i/2) ) }
for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(labels,average, i/2) ) }
#plot( 0:20 /2, precision_list )
points( 0:20 /2, sensitivity_list, col = "red", pch = 18, cex = 4  )
lines( 0:20 /2, sensitivity_list, col = "red"  )

#points( 0:20 /2, sensitivity_list, col = "red" )
points( 0:20 /2, specificity_list, col = "green", pch = 18, cex = 4  )
lines( 0:20 /2, specificity_list, col = "green"  )


# 12 exp [Weighted average]
SUM_MASK = 1:length(predictions_1) * 0 
SUM_MASK[which(predictions_1>0.0)] = 1
SUM_MASK[which(predictions_2>0.0)] = SUM_MASK[which(predictions_2>0.0)] + 1
SUM_MASK[which(predictions_3>0.0)] = SUM_MASK[which(predictions_3>0.0)] + 1
SUM_MASK[which(predictions_4>0.0)] = SUM_MASK[which(predictions_4>0.0)] + 1
SUM_MASK[which(predictions_5>0.0)] = SUM_MASK[which(predictions_5>0.0)] + 1
SUM_MASK[which(predictions_6>0.0)] = SUM_MASK[which(predictions_6>0.0)] + 1
SUM_MASK[which(predictions_7>0.0)] = SUM_MASK[which(predictions_7>0.0)] + 1
SUM_MASK[which(predictions_8>0.0)] = SUM_MASK[which(predictions_8>0.0)] + 1
SUM_MASK[which(predictions_9>0.0)] = SUM_MASK[which(predictions_9>0.0)] + 1
SUM_MASK[which(predictions_10>0.0)] = SUM_MASK[which(predictions_10>0.0)] + 1
SUM_MASK[which(predictions_11>0.0)] = SUM_MASK[which(predictions_11>0.0)] + 1
SUM_MASK[which(predictions_12>0.0)] = SUM_MASK[which(predictions_12>0.0)] + 1
average = predictions_1 + predictions_2 + predictions_3 + predictions_4 + predictions_5 + predictions_6 + predictions_7 + predictions_8 + predictions_9 + predictions_10 + predictions_11 + predictions_12
average = average / SUM_MASK
average[is.na(average)] = 0.0

precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(labels,average, i/2) ) }
for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,average, i/2) ) }
for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(labels,average, i/2) ) }
for( i in 0:20 ){ mcc_list = c( mcc_list, GetMCC(labels,average, i/2) ) }
for( i in 0:20 ){ f1_list = c( f1_list, GetF1Score(labels,average, i/2) ) }
for( i in 0:20 ){ acc_list = c( acc_list, GetACC(labels,average, i/2) ) }
#plot( 0:20 /2, precision_list )
points( 0:20 /2, sensitivity_list, col = "red", pch = 8, cex = 4  )
lines( 0:20 /2, sensitivity_list, col = "red"  )

#points( 0:20 /2, sensitivity_list, col = "red" )
points( 0:20 /2, specificity_list, col = "green", pch = 8, cex = 4 )
lines( 0:20 /2, specificity_list, col = "green"  )

dev.off()


plot( 0:20 /2, mcc_list, col = "red", pch = 15, cex = 4  )
plot( 0:20 /2, f1_list, col = "red", pch = 15, cex = 4  )
plot( 0:20 /2, acc_list, col = "red", pch = 15, cex = 4  )



# 12 exp [Normal average]
average2 = predictions_1 + predictions_2 + predictions_3 + predictions_4 + predictions_5 + predictions_6 + predictions_7 + predictions_8 + predictions_9 + predictions_10 + predictions_11 + predictions_12
average2 = average2 / 12

precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(labels,average2, i/2) ) }
for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,average2, i/2) ) }
for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(labels,average2, i/2) ) }
for( i in 0:20 ){ mcc_list = c( mcc_list, GetMCC(labels,average2, i/2) ) }
for( i in 0:20 ){ f1_list = c( f1_list, GetF1Score(labels,average2, i/2) ) }
for( i in 0:20 ){ acc_list = c( acc_list, GetACC(labels,average, i/2) ) }

plot( 0:20 /2, mcc_list, col = "red", pch = 15, cex = 2  )
plot( 0:20 /2, f1_list, col = "red", pch = 15, cex = 2  )
plot( 0:20 /2, acc_list, col = "red", pch = 15, cex = 2  )
 



#========================================================================             
# Firgure3E - no Final
#
# Top 400 rank list (with Average IS of 12 experiments)
#
#========================================================================

is2_Roth_Average2_old =  is2_bm1_A + is2_bm1_RA + is2_bm2_A + is2_bm2_Q + is2_bm3_A + is2_bm3_Q + is2_bm3_QL + is2_bm4_S4A + is2_bm4_S8A + is2_bm4_SQ + is2_bm6_SQ + is2_bm7_SQ 
is2_Roth_Average2_old = is2_Roth_Average2_old / 12
is2acc_Roth_Average2_old <<- AutoActivationCorrection2( is2_Roth_Average2_old )	
is2acc3_Roth_Average2_old <<- AutoActivationCorrection3( is2_Roth_Average2_old )	

Fig3D_ACC2( is2_Roth_Average2_old, KnownPPI )
Fig3D_ACC3( is2_Roth_Average2_old, KnownPPI )	# FINAL?

simple_plotMatrix(is2acc_Roth_Average2,0,4)
simple_plotMatrix(is2acc3_Roth_Average2,0,4)

par(mfrow=c(4,1))	
is2acc_order = order( is2acc_Roth_Average2, decreasing = TRUE )
barplot(is2acc_Roth_Average2[is2acc_order][1:400])	
barplot(HIPPIE[is2acc_order][1:400])	
barplot(BioGrid[is2acc_order][1:400])	
barplot(KnownPPI[is2acc_order][1:400])	

pdf( "Fig3E.melted_is2acc_Roth_Average2.pdf" )
par(mfrow=c(4,1))
melted_is2acc_Roth_Average2 = melt(is2acc_Roth_Average2)
is2acc_order2 = order( melted_is2acc_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc_Roth_Average2$value[is2acc_order2][1:400] ) 
barplot(HIPPIE[is2acc_order2][1:400])	
barplot(BioGrid[is2acc_order2][1:400])	
barplot(KnownPPI[is2acc_order2][1:400])	
write.table(melted_is2acc_Roth_Average2,file="/Users/jyang/melted_is2acc_Roth_Average2.txt",sep="\t")
dev.off()

pdf( "Fig3E.melted_is2acc_Roth_Average2.detail.pdf" )
par(mfrow=c(4,1))
melted_is2acc_Roth_Average2 = melt(is2acc_Roth_Average2)
is2acc_order2 = order( melted_is2acc_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc_Roth_Average2$value[is2acc_order2][1:15] ) 
barplot(HIPPIE[is2acc_order2][1:15])	
barplot(BioGrid[is2acc_order2][1:15])	
barplot(KnownPPI[is2acc_order2][1:15])	
dev.off()

#
# Can be improved more !
#
### AutoActivationCorrection3 + Symmetry
sym_is2acc_Roth_Average2 = SymmetricSum(is2acc_Roth_Average2,0)
roc(KnownPPI,sym_is2acc_Roth_Average2)
sym_is2acc3_Roth_Average2 = SymmetricSum(is2acc3_Roth_Average2,0)
roc(KnownPPI,sym_is2acc3_Roth_Average2)

par(mfrow=c(4,1))	
sym_is2acc3_order = order( sym_is2acc3_Roth_Average2, decreasing = TRUE )
barplot(sym_is2acc3_Roth_Average2[sym_is2acc3_order][1:400])	
barplot(HIPPIE[sym_is2acc3_order][1:400])	
barplot(BioGrid[sym_is2acc3_order][1:400])	
barplot(KnownPPI[sym_is2acc3_order][1:400])	

pdf( "Fig3E.melted_is2acc3_Roth_Average2.pdf" )
par(mfrow=c(4,1))
melted_sym_is2acc3_Roth_Average2 = melt(sym_is2acc3_Roth_Average2)
sym_is2acc3_order2 = order( melted_sym_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_sym_is2acc3_Roth_Average2$value[sym_is2acc3_order2][1:400] ) 
barplot(HIPPIE[sym_is2acc3_order2][1:400])	
barplot(BioGrid[sym_is2acc3_order2][1:400])	
barplot(KnownPPI[sym_is2acc3_order2][1:400])	
write.table(melted_sym_is2acc3_Roth_Average2,file="/Users/jyang/melted_sym_is2acc3_Roth_Average2.txt",sep="\t")
dev.off()

pdf( "Fig3E.melted_sym_is2acc3_Roth_Average2.detail.pdf" )
par(mfrow=c(4,1))
melted_sym_is2acc3_Roth_Average2 = melt(sym_is2acc3_Roth_Average2)
sym_is2acc3_order2 = order( melted_sym_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_sym_is2acc3_Roth_Average2$value[sym_is2acc3_order2][1:15] ) 
barplot(HIPPIE[sym_is2acc3_order2][1:15])	
barplot(BioGrid[sym_is2acc3_order2][1:15])	
barplot(KnownPPI[sym_is2acc3_order2][1:15])	
dev.off()



precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
bm_list = c()
for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ bm_list = c( bm_list, GetInformedness(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }

dev.off()

plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2  )
points( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2  )
plot( 0:100 /10, bm_list, col = "red", pch = 15, cex = 2  )




precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
bm_list = c()
for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ bm_list = c( bm_list, GetInformedness(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }

dev.off()

pdf( "is2acc3_Roth_Average2.f1_mcc.pdf" )
plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
points( 0:100 /10, mcc_list, col = "red", pch = 17, cex = 2  )
dev.off()
plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2  )
plot( 0:100 /10, bm_list, col = "red", pch = 15, cex = 2  )

View( cbind(0:100 /10, sensitivity_list, specificity_list) )



















#========================================================================             
# Firgure3E - Final
# without 1/8 Aba experiment, 1st experiment
# Top 400 rank list (with Average IS of 12 experiments)
#
#========================================================================

is2_Roth_Average2 =  is2_bm2_A + is2_bm2_Q + is2_bm3_A + is2_bm3_Q + is2_bm3_QL + is2_bm4_S4A + is2_bm4_SQ + is2_bm6_SQ + is2_bm7_SQ 
is2_Roth_Average2 = is2_Roth_Average2 / 9
is2acc_Roth_Average2 <<- AutoActivationCorrection2( is2_Roth_Average2 )	
is2acc3_Roth_Average2 <<- AutoActivationCorrection3( is2_Roth_Average2 )	

Fig3D_ACC2( is2_Roth_Average2, KnownPPI )
Fig3D_ACC3( is2_Roth_Average2, KnownPPI, filename="/Users/jyang/melted_is2acc3_Roth_Average2.final.txt"  )	# FINAL?

simple_plotMatrix(is2acc_Roth_Average2,0,4)
simple_plotMatrix(is2acc3_Roth_Average2,0,4)


par(mfrow=c(4,1))	
is2acc_order = order( is2acc_Roth_Average2, decreasing = TRUE )
barplot(is2acc_Roth_Average2[is2acc_order][1:400])	
barplot(HIPPIE[is2acc_order][1:400])	
barplot(BioGrid[is2acc_order][1:400])	
barplot(KnownPPI[is2acc_order][1:400])	

pdf( "Fig3E.melted_is2acc_Roth_Average2.pdf" )
par(mfrow=c(4,1))
melted_is2acc_Roth_Average2 = melt(is2acc_Roth_Average2)
is2acc_order2 = order( melted_is2acc_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc_Roth_Average2$value[is2acc_order2][1:400] ) 
barplot(HIPPIE[is2acc_order2][1:400])	
barplot(BioGrid[is2acc_order2][1:400])	
barplot(KnownPPI[is2acc_order2][1:400])	
write.table(melted_is2acc_Roth_Average2,file="/Users/jyang/melted_is2acc_Roth_Average2.txt",sep="\t")
dev.off()

pdf( "Fig3E.melted_is2acc_Roth_Average2.detail.pdf" )
par(mfrow=c(4,1))
melted_is2acc_Roth_Average2 = melt(is2acc_Roth_Average2)
is2acc_order2 = order( melted_is2acc_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc_Roth_Average2$value[is2acc_order2][1:15] ) 
barplot(HIPPIE[is2acc_order2][1:15])	
barplot(BioGrid[is2acc_order2][1:15])	
barplot(KnownPPI[is2acc_order2][1:15])	
dev.off()

pdf( "Fig3E.melted_is2acc3_Roth_Average2.full.pdf" )
par(mfrow=c(4,1))
melted_is2acc3_Roth_Average2 = melt(is2acc3_Roth_Average2)
is2acc3_order2 = order( melted_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc3_Roth_Average2$value[is2acc3_order2] ) 
barplot(HIPPIE[is2acc3_order2])	
barplot(BioGrid[is2acc3_order2])	
barplot(KnownPPI[is2acc3_order2])	
write.table(melted_is2acc3_Roth_Average2,file="/Users/jyang/melted_is2acc3_Roth_Average2.txt",sep="\t")
dev.off()

melted_is2acc3_Roth_Average2_FINAL_TABLE = cbind(melted_is2acc3_Roth_Average2$value[is2acc3_order2], HIPPIE[is2acc3_order2],BioGrid[is2acc3_order2])
View( melted_is2acc3_Roth_Average2_FINAL_TABLE )
write.table(melted_is2acc3_Roth_Average2_FINAL_TABLE,file="/Users/jyang/melted_is2acc3_Roth_Average2_FINAL_TABLE.txt",sep="\t")

## For boundary
xx = melted_is2acc3_Roth_Average2$value[is2acc3_order2] 
xx[400] = 10
xx[167] = 10
barplot(xx)


pdf( "Fig3E.melted_is2acc3_Roth_Average2.pdf" )
par(mfrow=c(4,1))
melted_is2acc3_Roth_Average2 = melt(is2acc3_Roth_Average2)
is2acc3_order2 = order( melted_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc3_Roth_Average2$value[is2acc3_order2][1:400] ) 
barplot(HIPPIE[is2acc3_order2][1:400])	
barplot(BioGrid[is2acc3_order2][1:400])	
barplot(KnownPPI[is2acc3_order2][1:400])	
write.table(melted_is2acc3_Roth_Average2,file="/Users/jyang/melted_is2acc3_Roth_Average2.txt",sep="\t")
dev.off()

pdf( "Fig3E.melted_is2acc3_Roth_Average2.detail.pdf" )
par(mfrow=c(4,1))
melted_is2acc3_Roth_Average2 = melt(is2acc3_Roth_Average2)
is2acc3_order2 = order( melted_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc3_Roth_Average2$value[is2acc3_order2][1:15] ) 
barplot(HIPPIE[is2acc3_order2][1:15])	
barplot(BioGrid[is2acc3_order2][1:15])	
barplot(KnownPPI[is2acc3_order2][1:15])	
dev.off()

pdf( "Fig3E.melted_is2acc3_Roth_Average2.cutoff_1.2.pdf" )
par(mfrow=c(4,1))
melted_is2acc3_Roth_Average2 = melt(is2acc3_Roth_Average2)
is2acc3_order2 = order( melted_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_is2acc3_Roth_Average2$value[is2acc3_order2][1:167] ) 
barplot(HIPPIE[is2acc3_order2][1:167])	
barplot(BioGrid[is2acc3_order2][1:167])	
barplot(KnownPPI[is2acc3_order2][1:167])	
dev.off()



#
# Can be improved more !
#
### AutoActivationCorrection3 + Symmetry
sym_is2acc_Roth_Average2 = SymmetricSum(is2acc_Roth_Average2,0)
roc(KnownPPI,sym_is2acc_Roth_Average2)
sym_is2acc3_Roth_Average2 = SymmetricSum(is2acc3_Roth_Average2,0)
roc(KnownPPI,sym_is2acc3_Roth_Average2)

par(mfrow=c(4,1))	
sym_is2acc3_order = order( sym_is2acc3_Roth_Average2, decreasing = TRUE )
barplot(sym_is2acc3_Roth_Average2[sym_is2acc3_order][1:400])	
barplot(HIPPIE[sym_is2acc3_order][1:400])	
barplot(BioGrid[sym_is2acc3_order][1:400])	
barplot(KnownPPI[sym_is2acc3_order][1:400])	

pdf( "Fig3E.melted_is2acc3_Roth_Average2.pdf" )
par(mfrow=c(4,1))
melted_sym_is2acc3_Roth_Average2 = melt(sym_is2acc3_Roth_Average2)
sym_is2acc3_order2 = order( melted_sym_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_sym_is2acc3_Roth_Average2$value[sym_is2acc3_order2][1:400] ) 
barplot(HIPPIE[sym_is2acc3_order2][1:400])	
barplot(BioGrid[sym_is2acc3_order2][1:400])	
barplot(KnownPPI[sym_is2acc3_order2][1:400])	
write.table(melted_sym_is2acc3_Roth_Average2,file="/Users/jyang/melted_sym_is2acc3_Roth_Average2.txt",sep="\t")
dev.off()

pdf( "Fig3E.melted_sym_is2acc3_Roth_Average2.detail.pdf" )
par(mfrow=c(4,1))
melted_sym_is2acc3_Roth_Average2 = melt(sym_is2acc3_Roth_Average2)
sym_is2acc3_order2 = order( melted_sym_is2acc3_Roth_Average2$value, decreasing = TRUE )
barplot( melted_sym_is2acc3_Roth_Average2$value[sym_is2acc3_order2][1:15] ) 
barplot(HIPPIE[sym_is2acc3_order2][1:15])	
barplot(BioGrid[sym_is2acc3_order2][1:15])	
barplot(KnownPPI[sym_is2acc3_order2][1:15])	
dev.off()



precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
bm_list = c()
for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }
for( i in 0:100 ){ bm_list = c( bm_list, GetInformedness(RefineMatrix(KnownPPI),RefineMatrix(is2acc_Roth_Average2), i/10) ) }

dev.off()

plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2  )
points( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2  )
plot( 0:100 /10, bm_list, col = "red", pch = 15, cex = 2  )




precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
bm_list = c()
for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }
for( i in 0:100 ){ bm_list = c( bm_list, GetInformedness(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average2), i/10) ) }

dev.off()

pdf( "is2acc3_Roth_Average2.f1_mcc.pdf" )
plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
points( 0:100 /10, mcc_list, col = "red", pch = 17, cex = 2  )
dev.off()
plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2  )
plot( 0:100 /10, bm_list, col = "red", pch = 15, cex = 2  )




#========================================================================             
# Firgure3F - Final
#
# Show: Sampling sensitivity 
#
#========================================================================


# 163
Fig3D_ACC3( is2_Roth_Average2, KnownPPI )
#Fig3D_ACC3( is2acc3_Roth_Average2, KnownPPI )

is2acc3_bm2_A = Fig3D_ACC3( is2_bm2_A, KnownPPI )	# FINAL? 3.1, 75
is2acc3_bm2_Q = Fig3D_ACC3( is2_bm2_Q, KnownPPI )	# FINAL? 3.3, 95
is2acc3_bm3_A = Fig3D_ACC3( is2_bm3_A, KnownPPI )	# FINAL? 2.2, 104
is2acc3_bm3_Q = Fig3D_ACC3( is2_bm3_Q, KnownPPI )	# FINAL? 3.4, 84
is2acc3_bm3_QL = Fig3D_ACC3( is2_bm3_QL, KnownPPI )	# FINAL? 0.8, 236
is2acc3_bm4_S4A = Fig3D_ACC3( is2_bm4_S4A, KnownPPI )	# FINAL? 3, 84
is2acc3_bm4_SQ = Fig3D_ACC3( is2_bm4_SQ, KnownPPI )	# FINAL? 3, 91
is2acc3_bm6_SQ = Fig3D_ACC3( is2_bm6_SQ, KnownPPI )	# FINAL? 3.9, 82
is2acc3_bm7_SQ = Fig3D_ACC3( is2_bm7_SQ, KnownPPI )	# FINAL? 3.7, 88

labels = pos_weak_index[is2_order]


#predictions1 = is2acc3_bm1_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_1 = roc( labels, predictions1  )  # 0.782
#predictions2 = is2acc3_bm1_RA[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_2 = roc( labels, predictions2  )   # 0.8239
predictions3 = is2acc3_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_3 = roc( labels, predictions3  )    # 0.8077
predictions4 = is2acc3_bm2_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_4 =  roc( labels, predictions4 )    # 0.7959
predictions5 = is2acc3_bm3_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_5 = roc( labels, predictions5  )    # 0.8571
predictions6 = is2acc3_bm3_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_6 =  roc( labels, predictions6  )    # 0.8035
predictions7 = is2acc3_bm3_QL[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_7 =  roc( labels, predictions7  )    # 0.7787
predictions8 = is2acc3_bm4_S4A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_8 = roc( labels, predictions8  )   # 0.8472
predictions9 = is2acc3_bm4_S8A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_9 = roc( labels, predictions9  )   # 0.7682
predictions10 = is2acc3_bm4_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_10 = roc( labels, predictions10  )   # 0.7765
predictions11 = is2acc3_bm6_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_11 = roc( labels, predictions11  )  # 0.7993
predictions12 = is2acc3_bm7_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc_12 = roc( labels, predictions12  ) # 0.8061


prefidction_final = is2acc3_Roth_Average2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; 
sum(is2acc3_Roth_Average2>=1.1)    # 174
sum( labels[which(prefidction_final >= 1.1)] )      # 71 / 94 = 0.755 
sum( labels[which(prefidction_final >= 1.2)] )      # 71 / 94 = 0.755 

sum(is2acc3_bm2_A>=1.2)     # 184
sum(is2acc3_bm2_Q>=1.2)     # 213 
sum(is2acc3_bm3_A>=1.2)     # 163
sum(is2acc3_bm3_Q>=1.2)     # 229
sum(is2acc3_bm3_QL>=1.2)    # 205
sum(is2acc3_bm4_S4A>=1.2)   # 287
sum(is2acc3_bm4_SQ>=1.2)    # 232
sum(is2acc3_bm6_SQ>=1.2)    # 198
sum(is2acc3_bm7_SQ>=1.2)    # 198


sum( labels[which(predictions1 >= 1.2)] )      # ?
sum( labels[which(predictions2 >= 1.2)] )      # ? 
sum( labels[which(predictions3 >= 1.2)] )      # 65
sum( labels[which(predictions4 >= 1.2)] )      # 65
sum( labels[which(predictions5 >= 1.2)] )      # 62
sum( labels[which(predictions6 >= 1.2)] )      # 65
sum( labels[which(predictions7 >= 1.2)] )      # 56
sum( labels[which(predictions8 >= 1.2)] )      # 62
sum( labels[which(predictions9 >= 1.2)] )      # ?
sum( labels[which(predictions10 >= 1.2)] )     # 63 
sum( labels[which(predictions11 >= 1.2)] )      # 64
sum( labels[which(predictions12 >= 1.2)] )      # 65 

is2acc3_overlap = (is2acc3_bm2_A>=1.2) + (is2acc3_bm2_Q>=1.2) + (is2acc3_bm3_A>=1.2) + (is2acc3_bm3_Q>=1.2) + (is2acc3_bm3_QL>=1.2) + (is2acc3_bm4_S4A>=1.2) + (is2acc3_bm4_SQ>=1.2) + (is2acc3_bm6_SQ>=1.2) + (is2acc3_bm7_SQ>=1.2)        
simple_plotMatrix(is2acc3_overlap,0,9)        
melted_is2acc3_overlap = melt(is2acc3_overlap)

#	par(mfrow=c(3,1)) 
#	is2_order = order( is2acc3_overlap[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ], decreasing = TRUE )
#	barplot( is2acc3_Roth_Average2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] )
#    barplot( is2acc3_overlap[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
#	barplot( pos_weak_index[is2_order] )        
#    
#	par(mfrow=c(2,1)) 
#	barplot( HIPPIE[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
#	barplot( BioGrid[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] )
#	dev.off()
#
#	par(mfrow=c(2,1)) 
#	is2_order = order( is2acc3_overlap[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ], decreasing = TRUE )
#	barplot( is2acc3_overlap[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
#	barplot( pos_weak_index[is2_order] )     

Fig3c_df = data.frame( 
        avgIS_final = is2acc3_Roth_Average2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ],
        overlap_cnt = is2acc3_overlap[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ],
        spot_test = pos_weak_index[1:192],
        hippie = HIPPIE[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ],
        biogrid = BioGrid[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ]
)
View(Fig3c_df)
Fig3c_df_sorted1 = Fig3c_df[order(Fig3c_df$avgIS_final,decreasing = TRUE),]
Fig3c_df_sorted2 = Fig3c_df_sorted1[order(Fig3c_df_sorted1$overlap_cnt,decreasing = TRUE),]
View(Fig3c_df_sorted2)

pdf( "Fig3c_part1.pdf")
par(mfrow=c(3,1)) 
barplot(Fig3c_df_sorted2$avgIS_final)
barplot(Fig3c_df_sorted2$overlap_cnt)
barplot(Fig3c_df_sorted2$spot_test)
dev.off()

pdf( "Fig3c_part2.pdf")
par(mfrow=c(3,1)) 
barplot(Fig3c_df_sorted2$hippie)
barplot(Fig3c_df_sorted2$biogrid)
dev.off()


View( cbind(pos_weak_index[is2_order], is2acc3_overlap[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) )

#===============================
# For Text - Final
#===============================

precision_list = c()
sensitivity_list = c()
specificity_list = c()
for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(labels,prefidction_final, i/10) ) }
for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(labels,prefidction_final, i/10) ) }
for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(labels,prefidction_final, i/10) ) }

View( cbind(0:100 /10, precision_list, sensitivity_list, specificity_list) )
GetInformation(labels,prefidction_final, 1.2) 


# total positives                  
sum(labels) #[1] 94     


Fig3f_raw_table = cbind( labels, predictions3, predictions4, predictions5, predictions6, predictions7, predictions8, predictions10, predictions11, predictions12 )
View( Fig3f_raw_table )     
write.table(Fig3f_raw_table,file="Fig3f_raw_table.txt",sep="\t")


Fig3f_raw_table = cbind( labels, predictions3, predictions4, predictions5, predictions6, predictions7, predictions8, predictions10, predictions11, predictions12, prefidction_final )
View( Fig3f_raw_table )     
write.table(Fig3f_raw_table,file="Fig3f_raw_table.2.txt",sep="\t")



#samp_sensitivity = data.frame( count=1:3, average = c( 63.0, 83.75, 84.4754901961 )/94*100, stdev = c( 2.74873708375, 2.96624528843, 3.03967026961 ) )
samp_sensitivity = data.frame( count=1:3, average = c( 88.0, 94.1, 95.6 ), stdev = c( 5.49, 3.56, 3.28 ) )

plot(samp_sensitivity[1:2], ylim=c(50,100), cex = 1, col = "dark blue", xlab="Number of Screens", ylab="Sample Sensitivity")
lines(samp_sensitivity[1:2], col = "dark blue")     
                
x = samp_sensitivity[,1]
y = samp_sensitivity[,2]   
std = samp_sensitivity[,3]    

segments(x, y-std,x, y+std, col = "dark blue")
epsilon = 0.02
segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")





#==========================================================
# For Text
#
#==========================================================
> sum( RefineMatrix(KnownPPI)[which( RefineMatrix(is2acc3_Roth_Average2) >=1.2 ) ] )
[1] 63
> sum( RefineMatrix(is2acc3_Roth_Average2) >=1.2 ) )
Error: unexpected ')' in "sum( RefineMatrix(is2acc3_Roth_Average2) >=1.2 ) )"
> sum( RefineMatrix(is2acc3_Roth_Average2) >=1.2 ) 
[1] 163
> 63/163
[1] 0.3865031

# FP rate same as crY2H method
random_vector = as.vector( RefineMatrix(is2acc3_Roth_Average2) )
PPI_knwon_vector = as.vector( RefineMatrix(KnownPPI) )
random_negative_vector = random_vector[which(PPI_knwon_vector==FALSE)]
set.seed(43)
fp_rates = c()
positive_ppi_pair_cnt = 163
for( i in 1:100 ){
    random_negative_Roth = sample( random_negative_vector, positive_ppi_pair_cnt )
    fp_cnt = sum( random_negative_Roth >=1.2 )  # 0.61% 
    fp_rates = c( fp_rates, fp_cnt/positive_ppi_pair_cnt*100.0 )
}
mean(fp_rates)
#[1] 1.97546 %
sd(fp_rates)
#[1] 1.033883 %

set.seed(43)
random_known_Roth = sample( PPI_knwon_vector, 163 )




'''
is2_Roth_Average3 =  is2_bm1_A + is2_bm1_RA + is2_bm2_A + is2_bm2_Q + is2_bm3_A + is2_bm3_Q + is2_bm3_QL + is2_bm4_S4A + is2_bm4_SQ + is2_bm6_SQ + is2_bm7_SQ 
is2_Roth_Average3 = is2_Roth_Average3 / 11
is2acc2_Roth_Average3 <<- AutoActivationCorrection2( is2_Roth_Average3 )	
is2acc3_Roth_Average3 <<- AutoActivationCorrection3( is2_Roth_Average3 )	

precision_list = c()
sensitivity_list = c()
specificity_list = c()
mcc_list = c()
f1_list = c()
acc_list = c()
bm_list = c()
for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }
for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }
for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }
for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }
for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }
for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }
for( i in 0:100 ){ bm_list = c( bm_list, GetInformedness(RefineMatrix(KnownPPI),RefineMatrix(is2acc3_Roth_Average3), i/10) ) }

dev.off()

pdf( "is2acc3_Roth_Average3.f1_mcc.pdf" )
plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
points( 0:100 /10, mcc_list, col = "red", pch = 17, cex = 2  )
dev.off()
plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2  )
plot( 0:100 /10, bm_list, col = "red", pch = 15, cex = 2  )









	## Overlap
    signal_to_noise = 1:7 * 0	
	for( v in 1:7 ){
	    i = v
	    #v = v / 2
	    print(v)
		#overlap = (ris2_3>v)+(ris2_4>v)+(ris2_5>v)+(ris2_6>v)+(ris2_7>v)+(ris2_8>v)+(ris2_9>v)+(ris2_10>v)+(ris2_11>v)+(ris2_12>v)
		overlap = (ris2_1>v)+(ris2_2>v)+(ris2_3>v)+(ris2_4>v)+(ris2_5>v)+(ris2_6>v)+(ris2_7>v)+(ris2_8>v)+(ris2_9>v)+(ris2_10>v)+(ris2_11>v)+(ris2_12>v)
		overlap2 = (is2_1>v)+(is2_2>v)+(is2_3>v)+(is2_4>v)+(is2_5>v)+(is2_6>v)+(is2_7>v)+(is2_8>v)+(is2_9>v)+(is2_10>v)+(is2_11>v)+(is2_12>v)
		x=hist(overlap)
		print(x$breaks)
		print(x$counts)
		barplot(x$counts[2:11])
		print(x$counts[1])
		print(x$counts[2])
		print(sum(x$counts) - x$counts[1]- x$counts[2])
		print( (sum(x$counts) - x$counts[1]- x$counts[2])/x$counts[2] )
        signal_to_noise[i] = (sum(x$counts) - x$counts[1]- x$counts[2])/x$counts[2]
		PPI_PerformanceCheck2(overlap2,KnownPPI)
		for( v2 in 1:7 ){
			print(v2)
			overlap3 = overlap2
			overlap3[which(overlap2<v2)] = 0
			PPI_PerformanceCheck2(overlap3,KnownPPI)
		}
		simple_plotMatrix(overlap,0,8)
	}
	barplot(signal_to_noise)
}
'''



#===============================
# Sub-sampling
#===============================


SubSamplingFig3 <-function()){
    #is2_Roth_Average2 =    is2_bm2_A + 
    #                       is2_bm2_Q + 
    #                       is2_bm3_A + 
    #                       is2_bm3_Q + 
    #                       is2_bm3_QL + 
    #                       is2_bm4_S4A + 
    #                       is2_bm4_SQ + 
    #                       is2_bm6_SQ + 
    #                       is2_bm7_SQ 

    SubSamplingMatrixTest( bm1_W )
    # [1] 308.748
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.4576836 0.7179296 0.8179954 0.8969450 0.9450417 0.9714622 0.9855292 0.9930094 0.9964302 0.9981970 0.9990845

    SubSamplingMatrixTest( bm1_RW )
    # [1] 347.0455
    # [1] 0.7882321 0.9582536 0.9795821 0.9901918 0.9950834 0.9974491 0.9986969 0.9993446 0.9996684 0.9998358 0.9999185

    SubSamplingMatrixTest( bm2_W ) # * work
    # [1] 150.4109
    # [1] 0.3880178 0.7176074 0.8307380 0.8987349 0.9465447 0.9720684 0.9854148 0.9928494 0.9962346 0.9981276 0.9990673

    SubSamplingMatrixTest( bm3_W ) # * work
    # [1] 224.7485
    # [1] 0.4518888 0.7319623 0.8309588 0.9047718 0.9515102 0.9750634 0.9866918 0.9935832 0.9968425 0.9984228 0.9992219

    SubSamplingMatrixTest( bm4_SW ) # * work
    # [1] 238.4045
    # [1] 0.4101534 0.7149865 0.8199269 0.8930060 0.9455548 0.9708755 0.9852450 0.9926619 0.9962001 0.9980608 0.9990158

    SubSamplingMatrixTest( bm5_SW )
    # [1] 435.5077
    # [1] 0.4110799 0.7308692 0.8414320 0.9068528 0.9501777 0.9740359 0.9870007 0.9937292 0.9967789 0.9983970 0.9992112

    SubSamplingMatrixTest( bm6_SW ) # * work
    # [1] 245.0197
    # [1] 0.4367912 0.7484425 0.8386175 0.9062779 0.9502322 0.9742940 0.9866391 0.9936979 0.9966760 0.9983540 0.9991687

    SubSamplingMatrixTest( bm7_SW ) # * work
    # [1] 315.6517
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.4141453 0.7327630 0.8334488 0.9027606 0.9501358 0.9736608 0.9868872 0.9936437 0.9966315 0.9983357 0.9991543



    SubSamplingMatrixTest( bm1_A )
    # [1] 242.0288
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.8639238 0.9725648 0.9858178 0.9912779 0.9956547 0.9975253 0.9988800 0.9995010 0.9997540 0.9998904 0.9999403

    SubSamplingMatrixTest( bm1_RA )
    # [1] 296.9951
    # [1] 0.8654932 0.9728076 0.9859045 0.9918713 0.9955833 0.9979107 0.9989842 0.9995454 0.9997638 0.9998896 0.9999394

    SubSamplingMatrixTest( bm2_A ) # * work
    # [1] 316.2714
    # [1] 0.9487533 0.9856435 0.9925150 0.9962973 0.9981756 0.9989649 0.9993494 0.9997806 0.9999077 0.9999555 0.9999765

    SubSamplingMatrixTest( bm2_Q ) # * work
    # [1] 128.0505
    # [1] 0.9048921 0.9778116 0.9870742 0.9933583 0.9969056 0.9983011 0.9990120 0.9996501 0.9998339 0.9999186 0.9999628

    SubSamplingMatrixTest( bm3_A ) # * work
    # [1] 153.8511
    # [1] 0.8769606 0.9716735 0.9878467 0.9921933 0.9964556 0.9980948 0.9990341 0.9995664 0.9998249 0.9999088 0.9999523

    SubSamplingMatrixTest( bm3_Q ) # * work
    # [1] 199.727
    # [1] 0.9233280 0.9837891 0.9911583 0.9958845 0.9978295 0.9989075 0.9994651 0.9996985 0.9998517 0.9999341 0.9999655

    SubSamplingMatrixTest( bm3_QL ) # * work
    # [1] 117.2554
    # [1] 0.9429685 0.9871510 0.9925873 0.9970698 0.9986315 0.9992843 0.9996309 0.9998338 0.9999308 0.9999657 0.9999782

    SubSamplingMatrixTest( bm4_S4A ) # * work
    # [1] 188.3225
    # [1] 0.8380926 0.9571933 0.9777017 0.9880397 0.9945209 0.9972051 0.9983713 0.9992588 0.9996597 0.9998317 0.9999203

    SubSamplingMatrixTest( bm4_S8A )
    # [1] 179.8999
    # [1] 0.5676598 0.8555189 0.9231547 0.9576047 0.9810423 0.9901957 0.9950077 0.9975963 0.9987920 0.9994207 0.9997086

    SubSamplingMatrixTest( bm4_SQ ) # * work
    # [1] 99.30457
    # [1] 0.8546917 0.9758633 0.9868873 0.9939708 0.9975621 0.9985668 0.9991655 0.9996448 0.9998364 0.9999219 0.9999574

    SubSamplingMatrixTest( bm6_SQ ) # * work
    # [1] 195.3383
    # [1] 0.9230139 0.9801467 0.9886758 0.9943960 0.9971094 0.9985032 0.9991715 0.9996155 0.9998510 0.9999306 0.9999683

    SubSamplingMatrixTest( bm7_SQ ) # * work
    # [1] 176.0602
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.9168854 0.9814225 0.9878279 0.9937448 0.9970801 0.9986275 0.9992585 0.9996947 0.9998618 0.9999345 0.9999661


    ## IS score - comparison
    m2_W_1X = GenerateRandomMatrix(bm2_W, 6084)
    m2_A_1X = GenerateRandomMatrix(bm2_A, 6084)
    m2_Q_1X = GenerateRandomMatrix(bm2_Q, 6084)

    is2_bm2_A_1X <<- InteractionScores(m2_W_1X,NormalMixture2(m2_A_1X),1,9 )	  
	is2_bm2_Q_1X <<- InteractionScores(m2_W_1X,NormalMixture2(m2_Q_1X),1,9 )	  

    CorMatrix(is2_bm2_A, is2_bm2_A_1X)
    # [1] 0.8803816
    CorMatrix(is2_bm2_Q, is2_bm2_Q_1X)
    # [1] 0.8589957

    CorMatrix(AutoActivationCorrection3(is2_bm2_A), AutoActivationCorrection3(is2_bm2_A_1X))
    # [1] 0.8732903
    CorMatrix(AutoActivationCorrection3(is2_bm2_Q), AutoActivationCorrection3(is2_bm2_Q_1X))
    # [1] 0.8476226



    ## IS score - comparison
    m3_W_1X = GenerateRandomMatrix(bm3_W, 6084)
    m3_A_1X = GenerateRandomMatrix(bm3_A, 6084)
    m3_Q_1X = GenerateRandomMatrix(bm3_Q, 6084)

    is2_bm3_A_1X <<- InteractionScores(m3_W_1X,NormalMixture2(m3_A_1X),1,9 )	  
	is2_bm3_Q_1X <<- InteractionScores(m3_W_1X,NormalMixture2(m3_Q_1X),1,9 )	  

    CorMatrix(is2_bm3_A, is2_bm3_A_1X)
    # [1] 0.9060741
    CorMatrix(is2_bm3_Q, is2_bm3_Q_1X)
    # [1] 0.8425248

    CorMatrix(AutoActivationCorrection3(is2_bm3_A), AutoActivationCorrection3(is2_bm3_A_1X))
    # [1] 0.8707725
    CorMatrix(AutoActivationCorrection3(is2_bm3_Q), AutoActivationCorrection3(is2_bm3_Q_1X))
    # [1] 0.8327507


    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm2_A), KnownPPI)
    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm2_A_1X), KnownPPI)

    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm2_Q), KnownPPI)
    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm2_Q_1X), KnownPPI)

    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm3_A), KnownPPI)
    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm3_A_1X), KnownPPI)

    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm3_Q), KnownPPI)
    PPI_PerformanceCheck2(AutoActivationCorrection3(is2_bm3_Q_1X), KnownPPI)



    # how many ppi pairs we could address with 1 MiSeq run
    #    	Raw read count	R1/R2 Mapping	PPI read count	Sqrt
    # Exp 1 (X71)	10277730	7997868	3618354	1902.2
    # Exp 2 (X71)	10378457	8208354	4231921	2057.2
    # Exp 3 (X71)	6321707	5349595	3200376	1789.0
    # Exp 4 (X71)	7227284	5731862	2679138	1636.8
    # Exp 5 (X71)	7444265	5923796	2991575	1729.6

}











#============================
# Fig S7. Comparison liquid gel and solid plate
#============================

solid_R1R2 = c(74.65, 80.29, 74.68, 81.81, 80.58, 80.20, 72.26)  
liquid_R1R2 = c(86.16, 85.39, 80.64, 81.68, 76.70, 81.65, 76.75)  
solid_PPI = c(39.38, 33.41, 35.51, 47.44, 37.74, 42.43, 33.14) 
liquid_PPI = c(56.60, 49.35, 42.03, 39.40, 34.51, 44.73, 34.00)

pdf("liquid_solid_comparison.pdf")              
par(mfrow=c(1,2))
boxplot( liquid_R1R2, solid_R1R2, names = c("Liquid", "Solid"), ylab="R1/R2 Mapping (%)" )  
boxplot( liquid_PPI, solid_PPI, names = c("Liquid", "Solid"), ylab="PPI read (%)" )   
dev.off()
                    

t.test( liquid_R1R2, solid_R1R2 )  
t.test( liquid_PPI, solid_PPI )

wilcox.test( liquid_R1R2, solid_R1R2 )  
wilcox.test( liquid_PPI, solid_PPI ) 



#========================
# sample compeletness
#========================
BasicStats(RefineMatrix(bm2_W+bm3_W)) # 92.0
BasicStats(RefineMatrix(bm2_W+bm4_SW)) # 92.9
BasicStats(RefineMatrix(bm2_W+bm6_SW)) # 90.1
BasicStats(RefineMatrix(bm2_W+bm7_SW)) # 90.5

BasicStats(RefineMatrix(bm3_W+bm4_SW)) # 93.8
BasicStats(RefineMatrix(bm3_W+bm6_SW)) # 91.7
BasicStats(RefineMatrix(bm3_W+bm7_SW)) # 91.8

BasicStats(RefineMatrix(bm4_SW+bm6_SW)) # 92.3
BasicStats(RefineMatrix(bm4_SW+bm7_SW)) # 92.6

BasicStats(RefineMatrix(bm6_SW+bm7_SW)) # 89.4

mean( c(92.0, 92.9, 90.1, 90.5, 93.8, 91.7, 91.8, 92.3, 92.6, 89.4) ) # 91.7




pdf("liquid_solid_comparison.sample_pair_complexity.pdf")              
par(mfrow=c(1,2))
boxplot( c(90.9, 84.6, 85.7), c(87.0,89.50) ,names=c("Liquid-gel","Solid"), ylab="Sample Complexity (%)", ylim=c(80,100) )
boxplot( c(97.9, 95.2, 95.8), c(96.3,97.4) ,names=c("Liquid-gel","Solid"), ylab="Pair Complexity (%)", ylim=c(80,100) )
dev.off()

> mean(c(90.9, 84.6, 85.7))
[1] 87.06667
> mean(c(87.0,89.50))
[1] 88.25
> mean(c(97.9, 95.2, 95.8))
[1] 96.3
> mean(c(96.3,97.4))
[1] 96.85






#=================================
# Toxicity
#=================================

par(mfrow=c(2,1))
bait_toxicity = ( rowSums(RefineMatrix0(bm1_W)==0) ) /76*100
prey_toxicity = ( colSums(RefineMatrix0(bm1_W)==0) ) /76*100

pdf("toxicity.pdf")
par(mfrow=c(2,1))
barplot( ( ( rowSums(RefineMatrix0(bm1_W)==0) ) )/76*100, ylim=c(0,100) )
barplot( ( ( colSums(RefineMatrix0(bm1_W)==0) ) )/76*100, ylim=c(0,100) )
dev.off()

barplot( ( sort( colSums(RefineMatrix0(bm1_W)==0) ) )/76*100 , ylim=c(0,100) )
barplot( ( sort( rowSums(RefineMatrix0(bm1_W)==0) ) )/76*100, ylim=c(0,100) )


plot( bait_toxicity, prey_toxicity )




