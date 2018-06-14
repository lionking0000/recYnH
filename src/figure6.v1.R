LoadKRPData()

#========================================================================             
# Firgure6A - Final
#
# Show: 
#
#========================================================================
	KRP_3AT_is2_average = ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi))/3
	KRP_3AT_is2acc3_average = AutoActivationCorrection3(KRP_3AT_is2_average)
	simple_plotMatrix( ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi))/3,0, filename="KRP.3AT.is2_average.pdf" ) 
	simple_plotMatrix( KRP_3AT_is2acc3_average,0, filename="KRP.3AT.is2acc3_average.pdf" ) 



#========================================================================             
# Firgure6A - Final
# with no 3AT
# Show: 
#
#========================================================================
	KRP_3AT_is2_average2 = ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi)+(is2_KRP_6_WH25_rpi))/4
	KRP_3AT_is2acc3_average2 = AutoActivationCorrection3(KRP_3AT_is2_average2)
	simple_plotMatrix( KRP_3AT_is2acc3_average2,0, filename="KRP.3AT.is2acc3_average2.pdf" ) 
	KRP_Performance(KRP_3AT_is2acc3_average2)	# 0.7557

	KRP_3AT_is2_average3 = ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_6_WH25_rpi))/2

	KRP_3AT_is2_average4 = ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_6_WH25_rpi))/3

	# single
	KRP_Performance(AAC3(is2_KRP_5_WH4AT1_rpi), filename="roc.KRP.3AT.is2_KRP_5_WH4AT1_rpi.pdf" ) # 0.7566 (3-AT single exp.)
	KRP_Performance(AAC3(is2_KRP_5_WH4AT2_rpi), filename="roc.KRP.3AT.is2_KRP_5_WH4AT2_rpi.pdf" ) # 0.7603 (3-AT single exp.)
	KRP_Performance(AAC3(is2_KRP_5_WH2AT_rpi), filename="roc.KRP.3AT.is2_KRP_5_WH2AT_rpi.pdf" )  # 0.7595 (3-AT single exp.)

	KRP_Performance(AAC3(is2_KRP_6_WH25_rpi), filename="roc.KRP.3AT.is2_KRP_6_WH25_rpi.pdf" ) 	# 0.7151 ( no 3-AT single exp. )
	
	KRP_Performance(AAC3(is2_KRP_1_WH_rpi)) # no 3AT # AUC 0.6978 (not used)
	KRP_Performance(AAC3(is2_KRP_2_WH_rpi)) # no 3AT # AUC 0.7218 (not used)
	KRP_Performance(AAC3(is2_KRP_3_WH_rpi)) # no 3AT # AUC 0.723 (not used)
	KRP_Performance(AAC3(is2_KRP_4_WH_rpi)) # no 3AT # AUC 0.7299 (not used)
	

	# double
	KRP_Performance(AAC3(KRP_3AT_is2_average3), filename="roc.KRP.3AT.is2acc3_average3.pdf" ) 	#0.7545
	
	# three
	KRP_Performance(AAC3(KRP_3AT_is2_average4), filename="roc.KRP.3AT.is2acc3_average4.pdf" ) 	# 0.7555
	
	# final
	KRP_Performance(KRP_3AT_is2acc3_average2, filename="roc.KRP.3AT.is2acc3_average2.pdf" ) 	# 0.7557 (average IS)



#========================================================================             
# Firgure6C - Final
# with no 3AT
# Show: 
#
#========================================================================
	KRP_Expected_labels = matrix(0,nrow=17,ncol=19)   
	diag( KRP_Expected_labels ) = 1

	precision_list = c()
	sensitivity_list = c()
	specificity_list = c()
	mcc_list = c()
	f1_list = c()
	acc_list = c()
	bm_list = c()
	for( i in 0:100 ){ precision_list = c( precision_list, GetPrecision(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }
	for( i in 0:100 ){ sensitivity_list = c( sensitivity_list, GetSensitivity(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }
	for( i in 0:100 ){ specificity_list = c( specificity_list, GetSpecificity(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }
	for( i in 0:100 ){ mcc_list = c( mcc_list, GetMCC(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }
	for( i in 0:100 ){ f1_list = c( f1_list, GetF1Score(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }
	for( i in 0:100 ){ acc_list = c( acc_list, GetACC(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }
	for( i in 0:100 ){ bm_list = c( bm_list, GetInformedness(KRP_Expected_labels,KRP_3AT_is2acc3_average2, i/10) ) }


	pdf( "KRP_3AT_is2acc3_average2.f1_mcc.pdf" )
	plot( 0:100 /10, f1_list, col = "green", pch = 15, cex = 2  )
	points( 0:100 /10, mcc_list, col = "red", pch = 17, cex = 2  )
	dev.off()

	View( cbind(0:100,f1_list) )
	simple_plotMatrix( KRP_3AT_is2acc3_average2 >= 2.1,0, 9 ) # filename="KRP.3AT.is2acc3_average2.pdf" ) 
	#2.1                          


#========================================================================             
# Firgure6F - Final
# with no 3AT
# Show: 
#
#========================================================================
RPI_usable_reads = c( sum(KRP_5_W_rpi),
					sum(KRP_5_WH4AT1_rpi),
					sum(KRP_5_WH4AT2_rpi),
					sum(KRP_5_WH2AT_rpi),
					sum(KRP_6_W25_rpi),
					sum(KRP_6_WH25_rpi)
					)

#[jyang@ant-login9 Nele]$ zcat 2017-11-15_MiSeq/S1_W_R1.fastq.gz | wc
#9463464 11829330 850478071
#[jyang@ant-login9 Nele]$ zcat 2017-11-15_MiSeq/S2_H.25AT_R1.fastq.gz | wc
#11379592 14224490 1022683900
#[jyang@ant-login9 Nele]$ zcat 2017-11-15_MiSeq/S3_H.25AT_R1.fastq.gz | wc
#9457816 11822270 849976157
#[jyang@ant-login9 Nele]$ zcat 2017-11-15_MiSeq/S4_H.5AT_R1.fastq.gz | wc
#11326232 14157790 1017883374
#[jyang@ant-login9 Nele]$ zcat 2017-11-20_MiSeq/S3_W25_R1.fastq.gz | wc
#10989432 13736790 987604022
#[jyang@ant-login9 Nele]$ zcat 2017-11-20_MiSeq/S4_WH25_R1.fastq.gz | wc
#10639796 13299745 956182808

Total_reads = c(9463464, 11379592, 9457816, 11326232, 10989432, 10639796 ) / 4.0

pdf("RPI_usable_reads.pdf")
boxplot(RPI_usable_reads/Total_reads)
dev.off()

mean(RPI_usable_reads/Total_reads)    #  ==>  [1] 0.5056674







#===============================
# Sub-sampling
#===============================
SubSamplingFig6 <-function()){
    SubSamplingMatrixTest( KRP_1_W_rpi, 17*19 )
	# [1] 12070.25
	# [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
 	# [1] 0.3933986 0.7522880 0.8498570 0.9180408 0.9603807 0.9795825 0.9878518 0.9929355 0.9976588 0.9985614 0.9991803

	SubSamplingMatrixTest( KRP_2_W_rpi, 17*19 )
	# [1] 4720.539
 	# [1] 0.4428156 0.7312261 0.8226212 0.9032787 0.9482229 0.9731037 0.9868459 0.9921469 0.9968752 0.9981743 0.9989734

	SubSamplingMatrixTest( KRP_3_W_rpi, 17*19 )
	# [1] 4058.7
 	# [1] 0.4432820 0.7569450 0.8471961 0.8979367 0.9359333 0.9664537 0.9849387 0.9902299 0.9963023 0.9982325 0.9988720

	SubSamplingMatrixTest( KRP_4_W_rpi, 17*19 )
	# [1] 6436.632
 	# [1] 0.3131525 0.6880250 0.7706062 0.8748581 0.9332977 0.9643718 0.9823829 0.9915383 0.9964222 0.9979962 0.9988760

	SubSamplingMatrixTest( KRP_5_W_rpi, 17*19 )
	# [1] 2413.05
 	# [1] 0.4098168 0.7463464 0.8464924 0.9137678 0.9523318 0.9750346 0.9859821 0.9925068 0.9970814 0.9982773 0.9990220

	SubSamplingMatrixTest( KRP_6_W10_rpi, 17*19 )
	# [1] 2440.146
 	# [1] 0.3467965 0.7841432 0.8446543 0.9118582 0.9476957 0.9718191 0.9846551 0.9928346 0.9971815 0.9982472 0.9991197

	SubSamplingMatrixTest( KRP_6_W25_rpi, 17*19 )
	# [1] 3216.675
 	# [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
 	# [1] 0.5107730 0.8158008 0.8665701 0.9222192 0.9502440 0.9745133 0.9872010 0.9943778 0.9978500 0.9985736 0.9992278




    SubSamplingMatrixTest( KRP_1_WH_rpi, 17*19 )
	# [1] 10631.1
 	# [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
 	# [1] 0.6512913 0.8923796 0.9406485 0.9621347 0.9812740 0.9907864 0.9954439 0.9976598 0.9989130 0.9992318 0.9997059

	SubSamplingMatrixTest( KRP_2_WH_rpi, 17*19 )
	# [1] 6168.514
 	# [1] 0.6899272 0.8959181 0.9356525 0.9629658 0.9799731 0.9916869 0.9946518 0.9971795 0.9984898 0.9991967 0.9996299

	SubSamplingMatrixTest( KRP_3_WH_rpi, 17*19 )
	# [1] 4533.257
 	# [1] 0.6406172 0.8827726 0.9214838 0.9498135 0.9780151 0.9909815 0.9940668 0.9970027 0.9986162 0.9991900 0.9995795

	SubSamplingMatrixTest( KRP_4_WH_rpi, 17*19 )
	# [1] 13450.88
 	# [1] 0.7822825 0.9200100 0.9446335 0.9655971 0.9821164 0.9932423 0.9945093 0.9977242 0.9990771 0.9993848 0.9997334

	SubSamplingMatrixTest( KRP_5_WH4AT1_rpi, 17*19 )
	# [1] 5253.378
 	# [1] 0.7158712 0.9549284 0.9777916 0.9867013 0.9929050 0.9976265 0.9984156 0.9994152 0.9996674 0.9997794 0.9998714

	SubSamplingMatrixTest( KRP_5_WH4AT2_rpi, 17*19 )
	# [1] 4247.48
 	# [1] 0.7600910 0.9495225 0.9807205 0.9836482 0.9917094 0.9974304 0.9984039 0.9993398 0.9996902 0.9997689 0.9998665

	SubSamplingMatrixTest( KRP_5_WH2AT_rpi, 17*19 )
	# [1] 5131.932
 	# [1] 0.7475609 0.9469496 0.9739225 0.9801570 0.9922613 0.9970000 0.9981283 0.9992525 0.9996356 0.9996776 0.9998664

	SubSamplingMatrixTest( KRP_6_WH10_rpi, 17*19 )
	# [1] 3921.613
 	# [1] 0.6723847 0.8891288 0.9331451 0.9534546 0.9742897 0.9885253 0.9934524 0.9970588 0.9988842 0.9992421 0.9996329

	SubSamplingMatrixTest( KRP_6_WH25_rpi, 17*19 )
	# [1] 4647.046
 	# [1] 0.1		0.5	      1         2         4         8         16        32        64        128       256             		
 	# [1] 0.6607082 0.8629340 0.9179954 0.9331774 0.9715179 0.9878764 0.9924377 0.9960942 0.9985954 0.9991568 0.9995721

	##is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi)+(is2_KRP_6_WH25_rpi
	KRP_SubSampling_WH = matrix(0,4,11)
	KRP_SubSampling_WH[1,] = SubSamplingMatrixTest( KRP_5_WH4AT1_rpi, 17*19 )
	KRP_SubSampling_WH[2,] = SubSamplingMatrixTest( KRP_5_WH4AT2_rpi, 17*19 )
	KRP_SubSampling_WH[3,] = SubSamplingMatrixTest( KRP_5_WH2AT_rpi, 17*19 )
	KRP_SubSampling_WH[4,] = SubSamplingMatrixTest( KRP_6_WH25_rpi, 17*19 )
	plot(c(0.1,0.5,1,2,4,8,16,32,64,128,256),colSums(KRP_SubSampling_WH)/4)
	
	x = c(0.1,0.5,1,2,4,8,16,32,64) #,128,256)
	y = (colSums(KRP_SubSampling_WH)/4)[1:9]
	std = (colSds(KRP_SubSampling_WH))[1:9]

	par(mfrow=c(2,3))
	plot(x,y,ylim=c(0.5,1.0))
	segments(x, y-std,x, y+std, col = "dark blue")
	epsilon = 1
	segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
	segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")



    ## IS score - comparison
    KRP_1_W_1X = GenerateRandomMatrix(KRP_1_W_rpi, 17*19)
    KRP_1_WH_1X = GenerateRandomMatrix(KRP_1_WH_rpi, 17*19)

	is_KRP_1_WH_rpi_1X = InteractionScores(KRP_1_W_1X,KRP_1_WH_1X,1,9 )	  
    is2_KRP_1_WH_rpi_1X = InteractionScores(KRP_1_W_1X,NormalMixture2(KRP_1_WH_1X),1,9 )	  
 
    CorMatrix(is2_KRP_1_WH_rpi, is2_KRP_1_WH_rpi_1X)
    # [1] 0.5566174

    CorMatrix(AutoActivationCorrection3(is2_KRP_1_WH_rpi), AutoActivationCorrection3(is2_KRP_1_WH_rpi_1X))
    # [1] 0.5583527

    # KnownPPI_P170 = HIPPIE_P170 | BioGrid_P170
    KRP_Performance(AutoActivationCorrection3(is_KRP_1_WH_rpi_1X))
    KRP_Performance(AutoActivationCorrection3(is2_KRP_1_WH_rpi))
	KRP_Performance(AutoActivationCorrection3(is2_KRP_1_WH_rpi_1X))



    ## IS score - comparison
    KRP_5_W_1X = GenerateRandomMatrix(KRP_5_W_rpi, 17*19*100)
    KRP_5_WH2AT_1X = GenerateRandomMatrix(KRP_5_WH2AT_rpi, 17*19*100)

    is_KRP_5_WH2AT_rpi_1X = InteractionScores(KRP_5_W_1X,KRP_5_WH2AT_1X,1,9 )	  
    is2_KRP_5_WH2AT_rpi_1X = InteractionScores(KRP_5_W_1X,NormalMixture2(KRP_5_WH2AT_1X),1,9 )	  
 
    CorMatrix(is2_KRP_5_WH2AT_rpi, is2_KRP_5_WH2AT_rpi_1X)
    # [1] 0.5566174

    CorMatrix(AutoActivationCorrection3(is2_KRP_5_WH2AT_rpi), AutoActivationCorrection3(is2_KRP_5_WH2AT_rpi_1X))
    # [1] 0.5583527

    # KnownPPI_P170 = HIPPIE_P170 | BioGrid_P170
	KRP_Performance(AutoActivationCorrection3(is_KRP_5_WH2AT_rpi_1X))
	KRP_Performance(AutoActivationCorrection3(is2_KRP_5_WH2AT_rpi))
	KRP_Performance(AutoActivationCorrection3(is2_KRP_5_WH2AT_rpi_1X))

	InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH2AT_rpi),1,9 )	  

}






#===================================
# FP rate same as crY2H method
#===================================
random_vector = as.vector( KRP_3AT_is2acc3_average2 )
RPI_knwon_vector = as.vector( KRP_Expected_labels )
random_negative_vector = random_vector[which(RPI_knwon_vector==FALSE)]
set.seed(43)
fp_rates = c()
fp_cnts = c()
positive_ppi_pair_cnt = 11
for( i in 1:100 ){
    random_negative_Roth = sample( random_negative_vector, positive_ppi_pair_cnt )
    fp_cnt = sum( random_negative_Roth >=2.1 )  # 0.61% 
	fp_cnts = c( fp_cnts, fp_cnt )
    fp_rates = c( fp_rates, fp_cnt/positive_ppi_pair_cnt*100.0 )
}
mean(fp_rates)
#[1] 0.18 %
sd(fp_rates)
#[1] 1.28 %



