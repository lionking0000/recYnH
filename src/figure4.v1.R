#================================================================================================================================================  
#
#	Functions for Fig4
#            
#================================================================================================================================================

Fig4D_NO_ACC <- function(average, KnownPPI_P170, UseSymmetryScore=FALSE, filename = FALSE){
    no_acc_average = average
    if (UseSymmetryScore){
        no_acc_average = SymmetricSum(no_acc_average,0)
    }    
    melted_no_acc_average = melt(no_acc_average[8:173,8:173])
    if (filename != FALSE){
        write.table(melted_no_acc_average,file=filename,sep="\t")
    }
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    for( i in 0:100 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(no_acc_average), i/10) ) }
    f1_list[is.na(f1_list)] = 0.0

    plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - no_acc_average"  )
    plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2, main="f1 - no_acc_average"  )
    plot( 0:100 /10, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - no_acc_average" )
    plot( 0:100 /10, precision_list, col = "green", pch = 15, cex = 2, main="precision - no_acc_average"  )
    
    plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2, main="acc - no_acc_average"  )
    plot( 0:100 /10, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - no_acc_average"  )

    print( roc(KnownPPI_P170, no_acc_average) )
    print( max(f1_list) )
    is_cutoff = (1:100/10)[order(f1_list,decreasing = TRUE)][1]
    print ( is_cutoff )
    print ( sum(no_acc_average>is_cutoff) )
    PlotParDomainInteractions(no_acc_average, is_cutoff )
    View( cbind(0:100 /10, f1_list) )
}

Fig4D_ACC2 <- function(average, KnownPPI_P170, UseSymmetryScore=FALSE, filename = FALSE){
    acc2_average = AutoActivationCorrection2( average )	
    if (UseSymmetryScore){
        acc2_average = SymmetricSum(acc2_average,0)
    }
    melted_acc2_average = melt(acc2_average[8:173,8:173])
    if (filename != FALSE){
        write.table(melted_acc2_average,file=filename,sep="\t")
    }
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    for( i in 0:100 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc2_average), i/10) ) }
    f1_list[is.na(f1_list)] = 0.0
    
    plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - acc2_average"  )
    plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2, main="f1 - acc2_average"  )
    plot( 0:100 /10, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - acc2_average" )
    plot( 0:100 /10, precision_list, col = "green", pch = 15, cex = 2, main="precision - acc2_average"  )
    
    plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2, main="acc - acc2_average"  )
    plot( 0:100 /10, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - acc2_average"  )

    print( roc(KnownPPI_P170, acc2_average) )
    print( max(f1_list) )
    is_cutoff = (1:100/10)[order(f1_list,decreasing = TRUE)][1]
    print ( is_cutoff )
    print ( sum(acc2_average>is_cutoff) )
    PlotParDomainInteractions(acc2_average, is_cutoff )
    View( cbind(0:100 /10, f1_list) )
}


AAC1 <- function(m){
    return( AutoActivationCorrection(m) )
}

AAC2 <- function(m){
    return( AutoActivationCorrection2(m) )
}

AAC3 <- function(m){
    return( AutoActivationCorrection3(m) )
}


Fig4D_ACC3 <- function(average, KnownPPI_P170, UseSymmetryScore=FALSE, filename = FALSE){
    acc3_average = AutoActivationCorrection3( average )	
    if (UseSymmetryScore){
        acc3_average = SymmetricSum(acc3_average,0)
    }   
    melted_acc3_average = melt(acc3_average[8:173,8:173])
    if (filename != FALSE){
        write.table(melted_acc3_average,file=filename,sep="\t")
    }
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    for( i in 0:100 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    for( i in 0:100 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(acc3_average), i/10) ) }
    f1_list[is.na(f1_list)] = 0.0

    pdf( "Fig4D_ACC3.f1_mcc.pdf" )
    plot( 0:100 /10, mcc_list, col = "red", pch = 17, cex = 2  )
    points( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
    dev.off()

    plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2, main="f1 - acc3_average"  )
    points( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - acc3_average"  )
    
    plot( 0:100 /10, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - acc3_average" )
    plot( 0:100 /10, precision_list, col = "green", pch = 15, cex = 2, main="precision - acc3_average"  )
    
    plot( 0:100 /10, acc_list, col = "red", pch = 15, cex = 2, main="acc - acc3_average"  )
    plot( 0:100 /10, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - acc3_average"  )

    plot( 0:100 /10, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - acc3_average"  )
    points( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2, main="f1 - acc3_average"  )


    print( roc(KnownPPI_P170, acc3_average) )
    print( max(f1_list) )
    is_cutoff = (1:100/10)[order(f1_list,decreasing = TRUE)][1]
    print ( is_cutoff )
    print ( sum(acc3_average>is_cutoff) )
    View( melt(acc3_average))
    PlotParDomainInteractions(acc3_average, is_cutoff )
    View( cbind(0:100 /10, f1_list) )
}

Fig4DLoadAverage <- function(){
    #average <<- is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_S8A + is2_bP170_4_SQ + is2_bP170_5_SA8 + is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    #average <<- average / SUM_MASK
    #average[is.na(average)] <<- 0.0
    #average[is.infinite(average)] <<- 0.0
    

    #=============================================
    # average2
    average2 <<- is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_S8A + is2_bP170_4_SQ + is2_bP170_5_SA8 + is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average2 <<- average2 / 10


    #=============================================
    # average2a
    # remove only is2_bP170_5_SA8
    average2a <<- is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_S8A + is2_bP170_4_SQ +  is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average2a <<- average2a / 9               


    #=============================================
    # average3
    #
    # Remove is2_bP170_4_S8A and is2_bP170_5_SA8 
    average3 <<- is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_SQ + is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average3 <<- average3 / 8


    #=============================================
    # average4
    average4 <<- is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_SQ + is2_bP170_5_SQ
    average4 <<- average4 / 4


    #=============================================
    # average5
    aac2_average5 <<- AAC2(is2_bP170_1_Q) + AAC2(is2_bP170_2_Q) + AAC2(is2_bP170_3_Q) + AAC2(is2_bP170_4_S4A) + AAC2(is2_bP170_4_S8A) + AAC2(is2_bP170_4_SQ) + AAC2(is2_bP170_5_SA8) + AAC2(is2_bP170_5_SQ) + AAC2(is2_bP170_6_SA1) + AAC2(is2_bP170_6_SA2)
    aac2_average5 <<- aac2_average5 / 10
    aac3_average5 <<- AAC3(is2_bP170_1_Q) + AAC3(is2_bP170_2_Q) + AAC3(is2_bP170_3_Q) + AAC3(is2_bP170_4_S4A) + AAC3(is2_bP170_4_S8A) + AAC3(is2_bP170_4_SQ) + AAC3(is2_bP170_5_SA8) + AAC3(is2_bP170_5_SQ) + AAC3(is2_bP170_6_SA1) + AAC3(is2_bP170_6_SA2)
    aac3_average5 <<- aac3_average5 / 10
}                       


#================================================================================================================================================ 







#========================================================================             
# Firgure4D - Final
#
# Show top 100
#
#========================================================================
    LoadP170Data()

	is2_bP170_1_Q <<- InteractionScores(bP170_1_W,NormalMixture2( bP170_1_Q ),1,4)     
	is2_bP170_2_Q <<- InteractionScores(bP170_2_W,NormalMixture2( bP170_2_Q ),1,4)        
	is2_bP170_3_Q <<- InteractionScores(bP170_3_W,NormalMixture2( bP170_3_Q ),1,4)  
	is2_bP170_4_S4A <<- InteractionScores(bP170_4_SW,NormalMixture2( bP170_4_S4A ),1,4)           
	is2_bP170_4_S8A <<- InteractionScores(bP170_4_SW,NormalMixture2( bP170_4_S8A ),1,4)  
	is2_bP170_4_SQ <<- InteractionScores(bP170_4_SW,NormalMixture2( bP170_4_SQ ),1,4)                                                             
	is2_bP170_5_SA8 <<- InteractionScores(bP170_5_SW,NormalMixture2( bP170_5_SA8 ),1,4)  
	is2_bP170_5_SQ <<- InteractionScores(bP170_5_SW,NormalMixture2( bP170_5_SQ ),1,4) 
	is2_bP170_6_SA1 <<- InteractionScores(bP170_6_SW1,NormalMixture2( bP170_6_SA1 ),1,4)  
	is2_bP170_6_SA2 <<- InteractionScores(bP170_6_SW2,NormalMixture2( bP170_6_SA2 ),1,4)

    Fig4DLoadAverage()

    Fig4D_ACC2(average,KnownPPI_P170)   # ??
    Fig4D_ACC3(average,KnownPPI_P170)   # ??

    # use all
    Fig4D_ACC2(average2,KnownPPI_P170,FALSE,"/Users/jyang/melted_is2acc2_bP170_Average2.txt") # ok, 0.5264, 1.8, 261
    Fig4D_ACC3(average2,KnownPPI_P170) # no mark-pard3, 0.5282, 1.7, 173
    
    # remove only is2_bP170_5_SA8
    Fig4D_ACC2(average2a,KnownPPI_P170,FALSE,"/Users/jyang/melted_is2acc2_bP170_Average2a.txt")  # ok, 0.5463, 1.7, 310
    Fig4D_ACC3(average2a,KnownPPI_P170)  # no mark-pard3, 0.5394, 1.8, 183

    # Remove is2_bP170_4_S8A and is2_bP170_5_SA8 
    Fig4D_ACC2(average3,KnownPPI_P170)  # ok, 0.5462, 1.8, 209
    Fig4D_ACC3(average3,KnownPPI_P170,FALSE,"/Users/jyang/melted_is2acc3_bP170_Average3.txt")  # ok, 0.5358, 1.6, 235   # FINAL
    is2acc3_bP170_Average3 = AAC3(average3)
    melted_is2acc3_bP170_Average3 = melt( is2acc3_bP170_Average3 )
    xx = sort(melted_is2acc3_bP170_Average3$value, decreasing = TRUE)
    barplot(xx)
    xx[100] = 10 
    xx[236] = 10 
    barplot(xx)           
	
	# 2018/03/09 # of detection
	is2acc3_bP170_Average3_overlap = (is2_bP170_1_Q>=1.6) + (is2_bP170_2_Q>=1.6) + (is2_bP170_3_Q>=1.6) + (is2_bP170_4_S4A>=1.6) + (is2_bP170_4_SQ>=1.6) + (is2_bP170_5_SQ>=1.6) + (is2_bP170_6_SA1>=1.6) + (is2_bP170_6_SA2>=1.6)
    melted_is2acc3_bP170_Average3_overlap = melt( is2acc3_bP170_Average3_overlap )


    # FP rate
    is2acc3_bP170_Average3
    # FP rate same as crY2H method
    random_vector = as.vector( is2acc3_bP170_Average3 )
    PPI_knwon_vector = as.vector( KnownPPI_P170 )
    random_negative_vector = random_vector[which(PPI_knwon_vector==FALSE)]
    set.seed(43)
    fp_rates = c()
    for( i in 1:100 ){
        random_negative_bP170 = sample( random_negative_vector, 235 )
        fp_cnt = sum( random_negative_bP170 >=1.6 )  # 0.61% 
        fp_rates = c( fp_rates, fp_cnt/235*100.0 )
    }
    mean(fp_rates)
    #[1] 0.78 %
    sd(fp_rates)
    #[1] 0.61 %


    # use Q
    Fig4D_ACC2(average4,KnownPPI_P170)  # no mark-pard3 0.5416, 1.5, 418 
    Fig4D_ACC3(average4,KnownPPI_P170)  # ok ,0.5422, 1.2, 371


    # use all
    Fig4D_ACC2(aac2_average5,KnownPPI_P170) # ok, 0.5279, 1.6, 294
    Fig4D_ACC3(aac2_average5,KnownPPI_P170) # no mark-pard3, 0.5273, 1.7, 171

    # use all
    Fig4D_ACC2(aac3_average5,KnownPPI_P170) # no mark-pard3, 0.539, 1.7, 187
    Fig4D_ACC3(aac3_average5,KnownPPI_P170) # no mark-pard3, 0.5311, 1.6, 172

    # use all
    Fig4D_NO_ACC(aac2_average5,KnownPPI_P170) # ok, 0.543, 1.8, 349
    Fig4D_NO_ACC(aac3_average5,KnownPPI_P170) # no mark-pard3, 1.7, 193


    #=================
    # Use symmetry sum value
    Fig4D_ACC2(average,KnownPPI_P170,TRUE)   # no, 0.5125, 2.4, 644
    Fig4D_ACC3(average,KnownPPI_P170,TRUE)   # no, 0.5059, 2.2, 455

    # use all
    Fig4D_ACC2(average2,KnownPPI_P170,TRUE)   # no, 0.5232, 2.6, 285
    Fig4D_ACC3(average2,KnownPPI_P170,TRUE)   # no, 0.5281, 2.2, 251
    
    # remove only is2_bP170_5_SA8
    Fig4D_ACC2(average2a,KnownPPI_P170,TRUE)   # no, 0.5532, 2.9, 251
    Fig4D_ACC3(average2a,KnownPPI_P170,TRUE)   # no, 0.541, 2.3, 251

    # Remove is2_bP170_4_S8A and is2_bP170_5_SA8 
    Fig4D_ACC2(average3,KnownPPI_P170,TRUE)   # no, 0.5526, 2.8, 299
    Fig4D_ACC3(average3,KnownPPI_P170,TRUE)   # ??

    # use Q
    Fig4D_ACC2(average4,KnownPPI_P170,TRUE)   # ??
    Fig4D_ACC3(average4,KnownPPI_P170,TRUE)   # ??

    # use all
    Fig4D_ACC2(aac2_average5,KnownPPI_P170,TRUE)   # ??
    Fig4D_ACC3(aac2_average5,KnownPPI_P170,TRUE)   # ??

    # use all
    Fig4D_ACC2(aac3_average5,KnownPPI_P170,TRUE)   # ??
    Fig4D_ACC3(aac3_average5,KnownPPI_P170,TRUE)   # ??

    # use all
    Fig4D_NO_ACC(aac2_average5,KnownPPI_P170,TRUE)   # ??
    Fig4D_NO_ACC(aac3_average5,KnownPPI_P170,TRUE)   # ??
























#==============================================================================

	gis2_bP170_1_Q <<- InteractionScores(bP170_1_W,GammaMixture2( bP170_1_Q ),1,4)     
	gis2_bP170_2_Q <<- InteractionScores(bP170_2_W,GammaMixture2( bP170_2_Q ),1,4)        
	gis2_bP170_3_Q <<- InteractionScores(bP170_3_W,GammaMixture2( bP170_3_Q ),1,4)  
	gis2_bP170_4_S4A <<- InteractionScores(bP170_4_SW,GammaMixture2( bP170_4_S4A ),1,4)           
	gis2_bP170_4_S8A <<- InteractionScores(bP170_4_SW,GammaMixture2( bP170_4_S8A ),1,4)  
	gis2_bP170_4_SQ <<- InteractionScores(bP170_4_SW,GammaMixture2( bP170_4_SQ ),1,4)                                                             
	gis2_bP170_5_SA8 <<- InteractionScores(bP170_5_SW,GammaMixture2( bP170_5_SA8 ),1,4)  
	gis2_bP170_5_SQ <<- InteractionScores(bP170_5_SW,GammaMixture2( bP170_5_SQ ),1,4) 
	gis2_bP170_6_SA1 <<- InteractionScores(bP170_6_SW1,GammaMixture2( bP170_6_SA1 ),1,4)  
	gis2_bP170_6_SA2 <<- InteractionScores(bP170_6_SW2,GammaMixture2( bP170_6_SA2 ),1,4)

    average4 = gis2_bP170_2_Q + gis2_bP170_3_Q + gis2_bP170_4_SQ + gis2_bP170_5_SQ
    average4 = average4 / 4
    is2acc2_bP170_Average4 <<- AutoActivationCorrection2( average4 )	


    roc(RefineMatrixP170_2(KnownPPI_P170), RefineMatrixP170_2(is2_bP170_2_Q))
    roc(RefineMatrixP170_2(KnownPPI_P170), RefineMatrixP170_2(is2acc2_bP170_Average4))
    roc(KnownPPI_P170[8:173,8:173], is2acc2_bP170_Average4[8:173,8:173])

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }

    plot( 0:20 /2, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - is2acc2_bP170_Average4"  )
    plot( 0:20 /2, f1_list, col = "green", pch = 15, cex = 4, main="f1 - is2acc2_bP170_Average4"  )
    plot( 0:20 /2, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - is2acc2_bP170_Average4" )
    plot( 0:20 /2, precision_list, col = "green", pch = 15, cex = 4, main="precision - is2acc2_bP170_Average4"  )
    
    plot( 0:20 /2, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average4"  )
    plot( 0:20 /2, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average4"  )




















    # 10 exp
	# Average with Sebastian way
	SUM_MASK = matrix(0,nrow=173,ncol=173)       
	SUM_MASK[which(is2_bP170_1_Q>0.0)] = 1
	SUM_MASK[which(is2_bP170_2_Q>0.0)] = SUM_MASK[which(is2_bP170_2_Q>0.0)] + 1
	SUM_MASK[which(is2_bP170_3_Q>0.0)] = SUM_MASK[which(is2_bP170_3_Q>0.0)] + 1
	SUM_MASK[which(is2_bP170_4_S4A>0.0)] = SUM_MASK[which(is2_bP170_4_S4A>0.0)] + 1
    SUM_MASK[which(is2_bP170_4_S8A>0.0)] = SUM_MASK[which(is2_bP170_4_S8A>0.0)] + 1
    SUM_MASK[which(is2_bP170_4_SQ>0.0)] = SUM_MASK[which(is2_bP170_4_SQ>0.0)] + 1
    SUM_MASK[which(is2_bP170_5_SA8>0.0)] = SUM_MASK[which(is2_bP170_5_SA8>0.0)] + 1
    SUM_MASK[which(is2_bP170_5_SQ>0.0)] = SUM_MASK[which(is2_bP170_5_SQ>0.0)] + 1
    SUM_MASK[which(is2_bP170_6_SA1>0.0)] = SUM_MASK[which(is2_bP170_6_SA1>0.0)] + 1
    SUM_MASK[which(is2_bP170_6_SA2>0.0)] = SUM_MASK[which(is2_bP170_6_SA2>0.0)] + 1

    average = is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_S8A + is2_bP170_4_SQ + is2_bP170_5_SA8 + is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average = average / SUM_MASK
    average[is.na(average)] = 0.0

    simple_plotMatrix(average>3,0,4)

    is2acc2_bP170_Average <<- AutoActivationCorrection2( average )	
    simple_plotMatrix(is2acc2_bP170_Average,0,4)

    melted_is2acc2_bP170_Average <<- melt(is2acc2_bP170_Average[8:173,8:173])

    simple_plotMatrix(melted_is2acc2_bP170_Average,0,4)
    melted_is2acc2_bP170_Average <<- melt(melted_is2acc2_bP170_Average[8:173,8:173])
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])

    pdf("Fig4D.is2acc2_bP170_Average.barplot.pdf")
    par(mfrow=c(4,1)) 
    is2_order = order( melted_is2acc2_bP170_Average$value, decreasing = TRUE )
    barplot( melted_is2acc2_bP170_Average$value[is2_order][1:100] ) 
    barplot( melted_HIPPIE_P170$value[is2_order][1:100] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:100] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:100] )   
    dev.off()

    # P-value
    sum(melted_KnownPPI_P170$value)  # [1] 585
    sum(melted_KnownPPI_P170$value[is2_order][1:100] )  #[1] 6
    dim(melted_KnownPPI_P170)  # [1] 27556     3
    # top 100
    # hypergeo 6 27556 585 100 1 ==> P-value = 0.00548974

    # P-value
    sum(melted_KnownPPI_P170$value)  # [1] 585
    sum(melted_KnownPPI_P170$value[is2_order][1:400] )  #[1] 16
    dim(melted_KnownPPI_P170)  # [1] 27556     3
    # top 400
    # hypergeo 16 27556 585 400 1 ==> P-value = 0.00563196

    melted_KnownPPI_P170$value[is2_order][1:30]
    library(taRifx)
    sort(melted_is2acc2_bP170_Average, f= ~ -value, decreasing = TRUE )

    #=============================================
    # average2a
    # remove only is2_bP170_5_SA8

    average2a = is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_S8A + is2_bP170_4_SQ +  is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average2a = average2a / 9
    is2acc2_bP170_Average2a <<- AutoActivationCorrection2( average2a )	
    is2acc2_bP170_Average2b <<- AutoActivationCorrection3( average2a )	
    simple_plotMatrix(is2acc2_bP170_Average2a,0,4)
    melted_is2acc2_bP170_Average2a <<- melt(is2acc2_bP170_Average2a[8:173,8:173])
    melted_is2acc2_bP170_Average2b <<- melt(is2acc2_bP170_Average2b[8:173,8:173])
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])
    write.table( melted_is2acc2_bP170_Average2a, )
    write.table(melted_is2acc2_bP170_Average2a,file="/Users/jyang/melted_is2acc2_bP170_Average2a.txt",sep="\t")
    write.table(melted_is2acc2_bP170_Average2b,file="/Users/jyang/melted_is2acc2_bP170_Average2b.txt",sep="\t")

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2a), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - is2acc2_bP170_Average2a"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 2, main="f1 - is2acc2_bP170_Average2a"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - is2acc2_bP170_Average2a" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 2, main="precision - is2acc2_bP170_Average2a"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average2a"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average2a"  )

    # cut off = 1.65 0.073817762 according to F1-score

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2b), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - is2acc2_bP170_Average2b"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 2, main="f1 - is2acc2_bP170_Average2b"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - is2acc2_bP170_Average2b" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 2, main="precision - is2acc2_bP170_Average2b"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average2b"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average2b"  )

    # cut off =  [35,]  1.70 0.074172185 according to F1-score


    #=============================================
    # average2
    average2 = is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_S8A + is2_bP170_4_SQ + is2_bP170_5_SA8 + is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average2 = average2 / 10
    is2acc2_bP170_Average2 <<- AutoActivationCorrection2( average2 )	
    is2acc3_bP170_Average2 <<- AutoActivationCorrection3( average2 )	
    
    simple_plotMatrix(is2acc2_bP170_Average2,0,4)
    melted_is2acc2_bP170_Average2 <<- melt(is2acc2_bP170_Average2[8:173,8:173])
    melted_is2acc3_bP170_Average2 <<- melt(is2acc3_bP170_Average2[8:173,8:173])
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])

    write.table(melted_is2acc2_bP170_Average2,file="/Users/jyang/melted_is2acc2_bP170_Average2.txt",sep="\t")
    write.table(melted_is2acc3_bP170_Average2,file="/Users/jyang/melted_is2acc3_bP170_Average2.txt",sep="\t")
    

    pdf("Fig4D.is2acc2_bP170_Average2.barplot.pdf")
    par(mfrow=c(4,1)) 
    is2_order = order( melted_is2acc2_bP170_Average2$value, decreasing = TRUE )
    barplot( melted_is2acc2_bP170_Average2$value[is2_order][1:100] ) 
    barplot( melted_HIPPIE_P170$value[is2_order][1:100] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:100] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:100] )   
    dev.off()


    # P-value
    sum(melted_KnownPPI_P170$value)  # [1] 585
    sum(melted_KnownPPI_P170$value[is2_order][1:100] )  #[1] 17
    dim(melted_KnownPPI_P170)  # [1] 27556     3
    # top 100
    # hypergeo 17 27556 585 100 1 ==> P-value = 3.61738e-12



    # P-value
    sum(melted_KnownPPI_P170$value)  # [1] 585
    sum(melted_KnownPPI_P170$value[is2_order][1:400] )  #[1] 17
    dim(melted_KnownPPI_P170)  # [1] 27556     3
    # top 400
    # hypergeo 32 27556 585 400 1 ==> P-value = 3.94161e-11

     melted_KnownPPI_P170$value[is2_order][1:30]
    # [1] FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE
    # [25] FALSE FALSE FALSE  TRUE FALSE FALSE
    sort(melted_is2acc2_bP170_Average2, f= ~ -value, decreasing = TRUE )


    is2_order = order( melted_is2acc2_bP170_Average2$value, decreasing = TRUE )
    barplot( melted_is2acc2_bP170_Average2$value[is2_order][1:400] ) 
    barplot( melted_HIPPIE_P170$value[is2_order][1:400] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:400] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:400] )   

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average2), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - is2acc2_bP170_Average2"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 2, main="f1 - is2acc2_bP170_Average2"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - is2acc2_bP170_Average2" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 2, main="precision - is2acc2_bP170_Average2"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average2"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average2"  )

    # 1.70 0.073985680 is the best according to f1 score


    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average2), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 2, main="mcc - is2acc3_bP170_Average2"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 2, main="f1 - is2acc3_bP170_Average2"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 2, main="sensitivity - is2acc3_bP170_Average2" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 2, main="precision - is2acc3_bP170_Average2"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc3_bP170_Average2"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc3_bP170_Average2"  )

    # 1.70 0.073985680 is the best according to f1 score


    #=============================================
    # average3 (Finale!)
    #
    # Remove is2_bP170_4_S8A and is2_bP170_5_SA8 
    average3 = is2_bP170_1_Q + is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_S4A + is2_bP170_4_SQ + is2_bP170_5_SQ + is2_bP170_6_SA1 + is2_bP170_6_SA2
    average3 = average3 / 8
    is2acc2_bP170_Average3 <<- AutoActivationCorrection2( average3 )	
    is2acc3_bP170_Average3 <<- AutoActivationCorrection3( average3 )	
    simple_plotMatrix(is2acc3_bP170_Average3,0,4)
    melted_is2acc3_bP170_Average3 <<- melt(is2acc3_bP170_Average3[8:173,8:173])
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])
    pdf("Fig4D.is2acc3_bP170_Average3.barplot.pdf")
    par(mfrow=c(4,1)) 
    is2_order = order( melted_is2acc3_bP170_Average3$value, decreasing = TRUE )
    barplot( melted_is2acc3_bP170_Average3$value[is2_order][1:400] ) 
    lines(x=c(0,500),y=c(1.6,1.6),col="red")
    lines(x=c(280,280),y=c(1.6,10),col="red")
    barplot( melted_HIPPIE_P170$value[is2_order][1:400] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:400] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:400] )   
    dev.off()
    write.table(melted_is2acc3_bP170_Average3,file="/Users/jyang/melted_is2acc3_bP170_Average3.txt",sep="\t")

    Fig4D(average3)

    #==========================
    # Test symmetry
    sym_is2acc2_bP170_Average3 = SymmetricSum(is2acc2_bP170_Average3,0)
    roc(HIPPIE_P170,is2acc2_bP170_Average3)
    roc(HIPPIE_P170,sym_is2acc2_bP170_Average3)
    melted_sym_is2acc2_bP170_Average3 <<- melt(sym_is2acc2_bP170_Average3[8:173,8:173])
    pdf("Fig4D.is2acc2_bP170_Average3.barplot.pdf")
    par(mfrow=c(4,1)) 
    is2_order = order( melted_sym_is2acc2_bP170_Average3$value, decreasing = TRUE )
    barplot( melted_sym_is2acc2_bP170_Average3$value[is2_order][1:100] ) 
    barplot( melted_HIPPIE_P170$value[is2_order][1:100] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:100] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:100] )   
    sum(melted_KnownPPI_P170$value[is2_order][1:100]) # 22
    sum(melted_KnownPPI_P170$value[is2_order][1:400]) # 38
    dev.off()
    write.table(melted_sym_is2acc2_bP170_Average3,file="/Users/jyang/melted_sym_is2acc2_bP170_Average3.txt",sep="\t")

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(sym_is2acc2_bP170_Average3), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - sym_is2acc2_bP170_Average3"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 4, main="f1 - sym_is2acc2_bP170_Average3"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - sym_is2acc2_bP170_Average3" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 4, main="precision - sym_is2acc2_bP170_Average3"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - sym_is2acc2_bP170_Average3"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - sym_is2acc2_bP170_Average3"  )
    max(f1_list)
    cbind( 0:200/20,f1_list)

    max(mcc_list)
    cbind( 0:200/20,mcc_list)   
    #==========================


    write.table(melted_is2acc2_bP170_Average3,file="/Users/jyang/melted_is2acc2_bP170_Average3.txt",sep="\t")

    pdf("Fig4D.is2acc2_bP170_Average3.barplot.pdf")
    par(mfrow=c(4,1)) 
    is2_order = order( melted_is2acc2_bP170_Average3$value, decreasing = TRUE )
    barplot( melted_is2acc2_bP170_Average3$value[is2_order][1:100] ) 
    barplot( melted_HIPPIE_P170$value[is2_order][1:100] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:100] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:100] )   
    sum(melted_KnownPPI_P170$value[is2_order][1:100]) # 18
    sum(melted_KnownPPI_P170$value[is2_order][1:400]) # 30
    dev.off()



    roc(KnownPPI_P170, is2acc2_bP170_Average3)
    roc(KnownPPI_P170[8:173,8:173], is2acc2_bP170_Average3[8:173,8:173])
    roc(RefineMatrixP170_2(KnownPPI_P170), RefineMatrixP170_2(is2acc2_bP170_Average3))

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - is2acc2_bP170_Average3"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 4, main="f1 - is2acc2_bP170_Average3"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - is2acc2_bP170_Average3" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 4, main="precision - is2acc2_bP170_Average3"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average3"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average3"  )

    max(f1_list)
    cbind( 0:200/20,f1_list)

    max(mcc_list)
    cbind( 0:200/20,mcc_list)


    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc2_bP170_Average3), i/2) ) }

    plot( 0:20 /2, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - is2acc2_bP170_Average3"  )
    plot( 0:20 /2, f1_list, col = "green", pch = 15, cex = 4, main="f1 - is2acc2_bP170_Average3"  )
    plot( 0:20 /2, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - is2acc2_bP170_Average3" )
    plot( 0:20 /2, precision_list, col = "green", pch = 15, cex = 4, main="precision - is2acc2_bP170_Average3"  )
    
    plot( 0:20 /2, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average3"  )
    plot( 0:20 /2, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average3"  )


    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }
    for( i in 0:20 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(HIPPIE_P170),RefineMatrixP170_2(is2acc3_bP170_Average3), i/2) ) }

    plot( 0:20 /2, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - is2acc3_bP170_Average3"  )
    plot( 0:20 /2, f1_list, col = "green", pch = 15, cex = 4, main="f1 - is2acc3_bP170_Average3"  )
    plot( 0:20 /2, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - is2acc3_bP170_Average3" )
    plot( 0:20 /2, precision_list, col = "green", pch = 15, cex = 4, main="precision - is2acc3_bP170_Average3"  )
    
    plot( 0:20 /2, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc3_bP170_Average3"  )
    plot( 0:20 /2, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc3_bP170_Average3"  )



   #=============================================
    # average4
    average4 = is2_bP170_2_Q + is2_bP170_3_Q + is2_bP170_4_SQ + is2_bP170_5_SQ
    average4 = average4 / 4
    is2acc2_bP170_Average4 <<- AutoActivationCorrection2( average4 )	
    is2acc3_bP170_Average4 <<- AutoActivationCorrection3( average4 )	

    melted_is2acc2_bP170_Average4 = melt(is2acc2_bP170_Average4)
    write.table(melted_is2acc2_bP170_Average3,file="/Users/jyang/melted_is2acc2_bP170_Average4.txt",sep="\t")

    melted_is2acc3_bP170_Average4 = melt(is2acc3_bP170_Average4)
    write.table(melted_is2acc2_bP170_Average3,file="/Users/jyang/melted_is2acc3_bP170_Average4.txt",sep="\t")

    roc(RefineMatrixP170_2(KnownPPI_P170), RefineMatrixP170_2(is2_bP170_2_Q))
    roc(RefineMatrixP170_2(KnownPPI_P170), RefineMatrixP170_2(is2acc2_bP170_Average4))
    roc(KnownPPI_P170[8:173,8:173], is2acc2_bP170_Average4[8:173,8:173])

    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:20 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }
    for( i in 0:20 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc2_bP170_Average4), i/2) ) }

    plot( 0:20 /2, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - is2acc2_bP170_Average4"  )
    plot( 0:20 /2, f1_list, col = "green", pch = 15, cex = 4, main="f1 - is2acc2_bP170_Average4"  )
    plot( 0:20 /2, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - is2acc2_bP170_Average4" )
    plot( 0:20 /2, precision_list, col = "green", pch = 15, cex = 4, main="precision - is2acc2_bP170_Average4"  )
    
    plot( 0:20 /2, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc2_bP170_Average4"  )
    plot( 0:20 /2, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc2_bP170_Average4"  )



    precision_list = c()
    sensitivity_list = c()
    specificity_list = c()
    mcc_list = c()
    f1_list = c()
    acc_list = c()
    fdr_list = c()
    for( i in 0:200 ){ precision_list = c( precision_list, GetPrecision(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }
    for( i in 0:200 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }
    for( i in 0:200 ){ specificity_list = c( specificity_list, GetSpecificity(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }
    for( i in 0:200 ){ mcc_list = c( mcc_list, GetMCC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }
    for( i in 0:200 ){ f1_list = c( f1_list, GetF1Score(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }
    for( i in 0:200 ){ acc_list = c( acc_list, GetACC(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }
    for( i in 0:200 ){ fdr_list = c( fdr_list, GetFDR(RefineMatrixP170_2(KnownPPI_P170),RefineMatrixP170_2(is2acc3_bP170_Average4), i/20) ) }

    plot( 0:200 /20, mcc_list, col = "red", pch = 15, cex = 4, main="mcc - is2acc3_bP170_Average4"  )
    plot( 0:200 /20, f1_list, col = "green", pch = 15, cex = 4, main="f1 - is2acc3_bP170_Average4"  )
    plot( 0:200 /20, sensitivity_list, col = "red", pch = 15, cex = 4, main="sensitivity - is2acc3_bP170_Average4" )
    plot( 0:200 /20, precision_list, col = "green", pch = 15, cex = 4, main="precision - is2acc3_bP170_Average4"  )
    
    plot( 0:200 /20, acc_list, col = "red", pch = 15, cex = 2, main="acc - is2acc3_bP170_Average4"  )
    plot( 0:200 /20, fdr_list, col = "red", pch = 15, cex = 2, main="fdr - is2acc3_bP170_Average4"  )

    max(f1_list)
    cbind( 0:200/20,f1_list)

    max(mcc_list)
    cbind( 0:200/20,mcc_list)


    #=============================================
    # average5
    aac2_average5 = AAC2(is2_bP170_1_Q) + AAC2(is2_bP170_2_Q) + AAC2(is2_bP170_3_Q) + AAC2(is2_bP170_4_S4A) + AAC2(is2_bP170_4_S8A) + AAC2(is2_bP170_4_SQ) + AAC2(is2_bP170_5_SA8) + AAC2(is2_bP170_5_SQ) + AAC2(is2_bP170_6_SA1) + AAC2(is2_bP170_6_SA2)
    aac3_average5 = AAC3(is2_bP170_1_Q) + AAC3(is2_bP170_2_Q) + AAC3(is2_bP170_3_Q) + AAC3(is2_bP170_4_S4A) + AAC3(is2_bP170_4_S8A) + AAC3(is2_bP170_4_SQ) + AAC3(is2_bP170_5_SA8) + AAC3(is2_bP170_5_SQ) + AAC3(is2_bP170_6_SA1) + AAC3(is2_bP170_6_SA2)
    is2acc2_bP170_Average5 = aac2_average5 / 10
    is2acc3_bP170_Average5 = aac3_average5 / 10
    # <<- AutoActivationCorrection2( average2 )	
    simple_plotMatrix(is2acc2_bP170_Average5,0,4)
    melted_is2acc2_bP170_Average5 <<- melt(is2acc2_bP170_Average5[8:173,8:173])
    simple_plotMatrix(is2acc3_bP170_Average5,0,4)
    melted_is2acc3_bP170_Average5 <<- melt(is2acc3_bP170_Average5[8:173,8:173])
    melted_HIPPIE_P170 = melt(HIPPIE_P170[8:173,8:173])
    melted_BioGrid_P170 = melt(BioGrid_P170[8:173,8:173])
    melted_KnownPPI_P170 = melt(KnownPPI_P170[8:173,8:173])

    write.table(melted_is2acc2_bP170_Average5,file="/Users/jyang/melted_is2acc2_bP170_Average5.txt",sep="\t")
    write.table(melted_is2acc2_bP170_Average5,file="/Users/jyang/melted_is2acc3_bP170_Average5.txt",sep="\t")

    pdf("Fig4D.melted_is2acc2_bP170_Average5.barplot.pdf")
    par(mfrow=c(4,1)) 
    is2_order = order( melted_is2acc2_bP170_Average2$value, decreasing = TRUE )
    barplot( melted_is2acc2_bP170_Average2$value[is2_order][1:100] ) 
    barplot( melted_HIPPIE_P170$value[is2_order][1:100] ) 
    barplot( melted_BioGrid_P170$value[is2_order][1:100] ) 
    barplot( melted_KnownPPI_P170$value[is2_order][1:100] )   
    dev.off()

	# Partial domain
    # cdc42_p154	17 
	# cdc42_p167	18
	# cdc42_p168	19    
	# mark2_p155	103
	# mark2_p156	104
	# mark2_p157	105
	# pard3b_p150   121 (full length)
	# pard3b_p158   122
	# pard3b_p159   123
	# pard3b_p160   124
	# pard6g_p152   125 (full length)
	# pard6g_p161   126
	# pard6g_p162   127
	# pard6g_p163   128
	# prkci_p153	135 (full length)
	# prkci_p164    136
	# prkci_p165    137
	# prkci_p166    138
    # ywhae_p169	171 (full length)
	# ywhae_p169	172
	# ywhae_p170	173
    par_domain_index = c(17,18,19,103,104,105,121,122,123,124,125,126,127,128,135,136,137,138,171,172,173)
    is2acc2_par_domain_Average5 = is2acc2_bP170_Average5[par_domain_index,par_domain_index]
    View(is2acc2_par_domain_Average5)
    simple_plotMatrix(is2acc2_par_domain_Average5,0,9)


PlotParDomainInteractions <- function(m, cutoff=1.6){
    #m = AAC3(m)
    par_domain_index = c(17,18,19,103,104,105,121,122,123,124,125,126,127,128,135,136,137,138,171,172,173)
    m2 = m[par_domain_index,par_domain_index]
    View(m2)
    simple_plotMatrix(m2,0,9)
    simple_plotMatrix(m2>cutoff,0,9)
}
    



    is2_bP170_1_Q <<- InteractionScores(bP170_1_W,NormalMixture2( bP170_1_Q ),1,4)     
	is2_bP170_2_Q <<- InteractionScores(bP170_2_W,NormalMixture2( bP170_2_Q ),1,4)        
	is2_bP170_3_Q <<- InteractionScores(bP170_3_W,NormalMixture2( bP170_3_Q ),1,4)  
	is2_bP170_4_S4A <<- InteractionScores(bP170_4_SW,NormalMixture2( bP170_4_S4A ),1,4)           
	is2_bP170_4_S8A <<- InteractionScores(bP170_4_SW,NormalMixture2( bP170_4_S8A ),1,4)  
	is2_bP170_4_SQ <<- InteractionScores(bP170_4_SW,NormalMixture2( bP170_4_SQ ),1,4)                                                             
	is2_bP170_5_SA8 <<- InteractionScores(bP170_5_SW,NormalMixture2( bP170_5_SA8 ),1,4)  
	is2_bP170_5_SQ <<- InteractionScores(bP170_5_SW,NormalMixture2( bP170_5_SQ ),1,4) 
	is2_bP170_6_SA1 <<- InteractionScores(bP170_6_SW1,NormalMixture2( bP170_6_SA1 ),1,4)  
	is2_bP170_6_SA2 <<- InteractionScores(bP170_6_SW2,NormalMixture2( bP170_6_SA2 ),1,4)

    PlotParDomainInteractions( is2_bP170_1_Q )
    PlotParDomainInteractions( is2_bP170_2_Q )
    PlotParDomainInteractions( is2_bP170_3_Q )
    PlotParDomainInteractions( is2_bP170_4_SQ )
    PlotParDomainInteractions( is2_bP170_5_SQ )
    PlotParDomainInteractions( is2_bP170_4_S4A )
    PlotParDomainInteractions( is2_bP170_4_S8A )
    PlotParDomainInteractions( is2_bP170_6_SA1 )
    PlotParDomainInteractions( is2_bP170_6_SA2 )



#===============================
# Sub-sampling
#===============================
SubSamplingFig4 <-function()){
    SubSamplingMatrixTest( bP170_1_W, 173*173 )
    # [1] 103.9542
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.4215782 0.7164512 0.8218303 0.9006979 0.9468755 0.9724838 0.9857249 0.9928469 0.9964846 0.9982378 0.9991207

    SubSamplingMatrixTest( bP170_2_W, 173*173 )
    # [1] 85.75759
    # [1] 0.5661058 0.8222459 0.8957739 0.9450278 0.9713987 0.9853617 0.9925009 0.9962171 0.9981061 0.9990524 0.9995281

    SubSamplingMatrixTest( bP170_3_W, 173*173 )
    # [1] 62.99616
    # [1] 0.6587371 0.8875425 0.9379617 0.9692375 0.9838239 0.9919460 0.9958481 0.9979037 0.9989558 0.9994993 0.9997380

    SubSamplingMatrixTest( bP170_4_SW, 173*173 )
    # [1] 63.83705
    # [1] 0.5532055 0.8304773 0.9028228 0.9493398 0.9734278 0.9862919 0.9929726 0.9964988 0.9982528 0.9991511 0.9995739

    SubSamplingMatrixTest( bP170_5_SW, 173*173 )
    # [1] 58.87701
    # [1] 0.5482939 0.8250627 0.9003366 0.9477419 0.9723864 0.9856722 0.9927143 0.9964064 0.9982430 0.9991230 0.9995612

    SubSamplingMatrixTest( bP170_6_SW1, 173*173 )
    # [1] 49.4079
    # [1] 0.5231145 0.8153941 0.8922641 0.9435595 0.9702516 0.9844334 0.9922482 0.9961774 0.9980994 0.9990567 0.9995217
    
    SubSamplingMatrixTest( bP170_6_SW2, 173*173 )
    # [1] 47.99181
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.5251846 0.8154513 0.8917816 0.9449631 0.9709014 0.9852624 0.9924869 0.9962463 0.9981057 0.9990462 0.9995093



    SubSamplingMatrixTest( bP170_1_Q, 173*173 )
    # [1] 113.7682
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.9704835 0.9942382 0.9972993 0.9984837 0.9991956 0.9996782 0.9998162 0.9998924 0.9999490 0.9999754 0.9999902

    SubSamplingMatrixTest( bP170_2_Q, 173*173 )
    # [1] 81.44743
    # [1] 0.9772102 0.9952467 0.9971063 0.9984451 0.9991678 0.9996391 0.9998181 0.9999055 0.9999478 0.9999743 0.9999873
    
    SubSamplingMatrixTest( bP170_3_Q, 173*173 )
    # [1] 54.34208
    # [1] 0.9884331 0.9971031 0.9985203 0.9990535 0.9996448 0.9998560 0.9999122 0.9999549 0.9999824 0.9999884 0.9999951

    SubSamplingMatrixTest( bP170_4_S4A, 173*173 )
    # [1] 62.07154
    # [1] 0.9274409 0.9826825 0.9916055 0.9958160 0.9978757 0.9991134 0.9995155 0.9997344 0.9998708 0.9999320 0.9999673

    SubSamplingMatrixTest( bP170_4_SQ, 173*173 )
    # [1] 40.03936
    # [1] 0.9678369 0.9930912 0.9967900 0.9981400 0.9990756 0.9996499 0.9998086 0.9999060 0.9999526 0.9999760 0.9999891

    SubSamplingMatrixTest( bP170_5_SQ, 173*173 )
    # [1] 60.36259
    # [1] 0.9586315 0.9925083 0.9964550 0.9983394 0.9990004 0.9995870 0.9997528 0.9998868 0.9999465 0.9999710 0.9999849
 
    SubSamplingMatrixTest( bP170_6_SA1, 173*173 )
    # [1] 50.16375
    # [1] 0.9384942 0.9860755 0.9927377 0.9964859 0.9981367 0.9990847 0.9995304 0.9997598 0.9998805 0.9999446 0.9999696

    SubSamplingMatrixTest( bP170_6_SA2, 173*173 )
    # [1] 45.10371
    # [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
    # [1] 0.9450615 0.9868654 0.9934561 0.9968445 0.9982771 0.9991855 0.9996027 0.9998051 0.9998986 0.9999481 0.9999744





    ## IS score - comparison
    m1_W_1X = GenerateRandomMatrix(bP170_1_W, 173*173)
    m1_Q_1X = GenerateRandomMatrix(bP170_1_Q, 173*173)

    is2_m1_P170_Q_1X <<- InteractionScores(m1_W_1X,NormalMixture2(m1_Q_1X),1,9 )	  
 
    CorMatrix(is2_bP170_1_Q, is2_m1_P170_Q_1X)
    # [1] 0.7730144

    CorMatrix(AutoActivationCorrection3(is2_bP170_1_Q), AutoActivationCorrection3(is2_m1_P170_Q_1X))
    # [1] 0.4550471

    # KnownPPI_P170 = HIPPIE_P170 | BioGrid_P170
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_1_Q), HIPPIE_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m1_P170_Q_1X), HIPPIE_P170)

    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_1_Q), BioGrid_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m1_P170_Q_1X), BioGrid_P170)

    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_1_Q), KnownPPI_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m1_P170_Q_1X), KnownPPI_P170)



    ## IS score - comparison
    m4_W_1X = GenerateRandomMatrix(bP170_4_SW, 173*173)
    m4_Q_1X = GenerateRandomMatrix(bP170_4_SQ, 173*173)

    is2_m4_P170_Q_1X <<- InteractionScores(m4_W_1X,NormalMixture2(m4_Q_1X),1,9 )	  
 
    CorMatrix(is2_bP170_4_SQ, is2_m4_P170_Q_1X)
    # [1] 0.7225343

    CorMatrix(AutoActivationCorrection3(is2_bP170_4_SQ), AutoActivationCorrection3(is2_m4_P170_Q_1X))
    # [1] 0.6886603

    # KnownPPI_P170 = HIPPIE_P170 | BioGrid_P170
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_4_SQ), HIPPIE_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m4_P170_Q_1X), HIPPIE_P170)

    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_4_SQ), BioGrid_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m4_P170_Q_1X), BioGrid_P170)

    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_4_SQ), KnownPPI_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m4_P170_Q_1X), KnownPPI_P170)



    ## IS score - comparison
    m6_W1_1X = GenerateRandomMatrix(bP170_6_SW1, 173*173)
    m6_A1_1X = GenerateRandomMatrix(bP170_6_SA1, 173*173)

    is2_m6_P170_A1_1X <<- InteractionScores(m6_W1_1X,NormalMixture2(m6_A1_1X),1,9 )	  
 
    CorMatrix(is2_bP170_6_SA1, is2_m6_P170_A1_1X)
    # [1] 0.7225343

    CorMatrix(AutoActivationCorrection3(is2_bP170_6_SA1), AutoActivationCorrection3(is2_m6_P170_A1_1X))
    # [1] 0.6886603

    # KnownPPI_P170 = HIPPIE_P170 | BioGrid_P170
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_6_SA1), HIPPIE_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m6_P170_A1_1X), HIPPIE_P170)

    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_6_SA1), BioGrid_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m6_P170_A1_1X), BioGrid_P170)

    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_bP170_6_SA1), KnownPPI_P170)
    PPI_PerformanceCheck_P170(AutoActivationCorrection3(is2_m6_P170_A1_1X), KnownPPI_P170)

}