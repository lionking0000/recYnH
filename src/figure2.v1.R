######################################################################################################################################  
#
# Version 1.0
# figures.R
#
#
######################################################################################################################################


###################################################################################################################################### 
## Initilize functions and load data
source( "/Volumes/users/lserrano/jyang/work/Mireia/src/init.R" )    

######################################################################################################################################


#=====================================
# Fig2a.v1                            
#=====================================

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/   
# These are color-blind-friendly palettes, one with gray, and one with black.
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")                  
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

usable_reads = data.frame( names=c("R1/R2\nmapping","PPI\nread"), average = c( 79.76, 42.18 ), stdev = c( 5.00, 7.76 ) )

pdf("usable_reads.75x.pdf",width=7, height=7)
par(mfrow=c(1,2))
barCenters <- barplot(height = usable_reads$average, beside = TRUE, las = 1, ylim = c(0, 100), cex.names = 1.5,    font.lab = 2,
                      cex.axis = 1.5,
                      cex.lab = 1.8,
                      col = cbPalette[1:2],
                      main = "Coverage",
                      ylab = "% of total reads",
                      xlab = "Usable reads",   
					  names.arg=usable_reads$names,
                      border = "black", axes = TRUE,
                      legend.text = TRUE,
                      args.legend = list(title = "No. Cylinders", 
                                         x = 100,
                                         y = 80,
                                         bty = "y",
                                         cex = .7))
# Add in line segments and caps for error bars
segments(barCenters, usable_reads$average - usable_reads$stdev, barCenters, usable_reads$average + usable_reads$stdev, lwd = 1.5)
arrows(barCenters, usable_reads$average - usable_reads$stdev, barCenters, usable_reads$average + usable_reads$stdev, lwd = 1.5, angle = 90, code = 3, length = 0.05)
dev.off()


#=====================================
# Fig2b.v1
#=====================================

# 75x
#p<- ggplot(samp_complx, aes(x=count, y=average, fill=cbPalette[1:5])) + geom_bar(stat="identity", color="black", position=position_dodge()) + geom_errorbar(aes(ymin=average-stdev, ymax=average+stdev), width=.2, position=position_dodge(.9))        

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/   
# These are color-blind-friendly palettes, one with gray, and one with black.
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")                  
# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

samp_complx = data.frame( count=1:5, average = c( 86.78, 91.29, 93.23, 94.43, 95.29 ), stdev = c( 3.91, 1.76, 1.09, 0.69, 0 ) )
pair_complx = data.frame( count=1:5, average = c( 94.19, 96.85, 97.94, 98.54, 98.94 ), stdev = c( 2.23, 0.96, 0.57, 0.34, 0.0 ) )

# Construct our barplot using our means matrix
# Specify some custom colors
# Play around with legend positioning    
pdf("coverage.75x.pdf",width=7, height=7)
par(mfrow=c(1,2))
barCenters <- barplot(height = samp_complx$average, beside = TRUE, las = 1, ylim = c(0, 100), cex.names = 1.5,    font.lab = 2,
                      cex.axis = 1.5,
                      cex.lab = 1.8,
                      col = cbPalette[1:5],
                      main = "Coverage",
                      ylab = "Sampling complexity",
                      xlab = "Number of screens",   
					  names.arg=samp_complx$count,
                      border = "black", axes = TRUE,
                      legend.text = TRUE,
                      args.legend = list(title = "No. Cylinders", 
                                         x = 100,
                                         y = 80,
                                         bty = "y",
                                         cex = .7))
# Add in line segments and caps for error bars
segments(barCenters, samp_complx$average - samp_complx$stdev, barCenters, samp_complx$average + samp_complx$stdev, lwd = 1.5)
arrows(barCenters, samp_complx$average - samp_complx$stdev, barCenters, samp_complx$average + samp_complx$stdev, lwd = 1.5, angle = 90, code = 3, length = 0.05)

barCenters <- barplot(height = pair_complx$average, beside = TRUE, las = 1, ylim = c(0, 100), cex.names = 1.5,  
                      cex.axis = 1.5,
                      cex.lab = 1.8,
                      col = cbPalette[1:5],
                      main = "Coverage",
                      ylab = "Pair complexity",
                      xlab = "Number of screens",   
					  names.arg=pair_complx$count,
                      border = "black", axes = TRUE,
                      legend.text = TRUE,
                      args.legend = list(title = "No. Cylinders", 
                                         x = 100,
                                         y = 80,
                                         bty = "y",
                                         cex = .7))
# Add in line segments and caps for error bars
segments(barCenters, pair_complx$average - pair_complx$stdev, barCenters, pair_complx$average + pair_complx$stdev, lwd = 1.5)
arrows(barCenters, pair_complx$average - pair_complx$stdev, barCenters, pair_complx$average + pair_complx$stdev, lwd = 1.5, angle = 90, code = 3, length = 0.05)   
dev.off()        


#========================================================================             
# Firgure2C - Final
#
# BD-Fused proteins against Read Count
#
#========================================================================

pdf("Sample Complexity.pdf")

plot(samp_complx[1:2], ylim=c(80,100), cex = 1, col = "dark blue", xlab="Number of Screens", ylab="Sample Complexity")
lines(samp_complx[1:2], col = "dark blue")     
points(pair_complx[1:2],col = "dark red", cex = 1)
lines(pair_complx[1:2], col = "dark red")  
                
x = samp_complx[,1]
y = samp_complx[,2]   
std = samp_complx[,3]    

segments(x, y-std,x, y+std, col = "dark blue")
epsilon = 0.02
segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")

x = pair_complx[,1]
y = pair_complx[,2]   
std = pair_complx[,3]    

segments(x, y-std,x, y+std,col = "dark red")
epsilon = 0.02
segments(x-epsilon,y-std,x+epsilon,y-std,col = "dark red")
segments(x-epsilon,y+std,x+epsilon,y+std,col = "dark red")

dev.off()


#========================================================================             
# Firgure2B - Final
#
# BD-Fused proteins against Read Count
#
#========================================================================
if (TRUE){
	# Fig 2B
	AA_Index = setdiff( 1:78, c(43,75) )
	AA_from_AD_empty_screening_A2 = Roth75_AA_Test$A2_R1 / Roth75_AA_Test$W_R1
	AA_from_AD_empty_screening_Q = Roth75_AA_Test$Q_R1 / Roth75_AA_Test$W_R1 
	AA_from_AD_empty_screening = ( AA_from_AD_empty_screening_A2 + AA_from_AD_empty_screening_Q ) / 2
	
	AA_from_full_screening = (rowSums(bm1_A)/sum(bm1_A))/(rowSums(bm1_W)/sum(bm1_W))
	names(AA_from_AD_empty_screening) = names(AA_from_full_screening)
	#plot( AA_from_AD_empty_screening_A2[AA_Index], AA_from_full_screening[AA_Index] )
	#plot( AA_from_AD_empty_screening_Q[AA_Index], AA_from_full_screening[AA_Index] )
	plot( AA_from_AD_empty_screening[AA_Index], AA_from_full_screening[AA_Index] )
	#cor( AA_from_AD_empty_screening_A2[AA_Index], AA_from_full_screening[AA_Index] )
	#cor( AA_from_AD_empty_screening_Q[AA_Index], AA_from_full_screening[AA_Index] )	
	cor( AA_from_AD_empty_screening[AA_Index], AA_from_full_screening[AA_Index] )
	cor.test( AA_from_AD_empty_screening[AA_Index], AA_from_full_screening[AA_Index] )
	
	mean( AA_from_AD_empty_screening[AA_Index] )
	mean( AA_from_full_screening[AA_Index]  )
	sort( AA_from_full_screeningAA_Index] )

 	(37041 + 78740 + 80491 + 151955 + 575790) / sum( rowSums(bm1_A) ) # 62.7

	pdf("Fig2B.pdf",width=4, height=4, family="ArialMT")
	Roth75_AA = cbind((Roth75_AA_Test$A2_R1/sum(Roth75_AA_Test$A2_R1))/(Roth75_AA_Test$W_R1/sum(Roth75_AA_Test$W_R1)), (Roth75_AA_Test$Q_R1/sum(Roth75_AA_Test$Q_R1))/(Roth75_AA_Test$W_R1/sum(Roth75_AA_Test$W_R1)))
    x = 1:78
    y = apply( Roth75_AA, 1, mean )       
    z = (rowSums(bm1_A)/sum(bm1_A))/(rowSums(bm1_W)/sum(bm1_W))  # From row sum
    std = apply( Roth75_AA, 1, sd )   
	ste = std / sqrt(2)                                                                              
	epsilon = 1
	if (FALSE){
		# std
		plot(x,y,cex=1,xlab="BD-Fused proteins",ylab="Read Count Ratio")
		segments(x, y-std,x, y+std, col = "dark blue")
		segments(x-epsilon,y-std,x+epsilon,y-std, col = "dark blue")
		segments(x-epsilon,y+std,x+epsilon,y+std, col = "dark blue")  
	}else{
		# ste
		plot(x,y,cex=1,xlab="BD-Fused proteins",ylab="Read Count Ratio",ylim=c(0,90))
		segments(x, y-ste,x, y+ste, col = "dark blue")
		segments(x-epsilon,y-ste,x+epsilon,y-ste, col = "dark blue")
		segments(x-epsilon,y+ste,x+epsilon,y+ste, col = "dark blue")  
	}
	points(x, (rowSums(bm1_A)/sum(bm1_A))/(rowSums(bm1_W)/sum(bm1_W)), cex=1, col="dark red")
	dev.off()      
	
	# Fig 2B - Supplement
	pdf("Fig2B.Supplement.pdf",width=5, height=5, family="ArialMT")        
	par(mfrow=c(2,3))
	plot(disorder$DISORDER.Residues, y, xlab="# of Disorder Residues",ylab="Read Count Ratio")
	plot(disorder$DISORDER.., y, xlab="Disorder Residues Percentage (%)",ylab="Read Count Ratio")  
	plot(disorder$Total.Residues, y, xlab="# of Total Residues",ylab="Read Count Ratio")   
	
	cor.test(disorder$DISORDER.Residues, y) # r = 0.06, p-value = 0.5581          
	cor.test(disorder$DISORDER.., y) # r = 0.07, p-value = 0.543         
	cor.test(disorder$Total.Residues, y) # r = 0.32, p-value = 0.004082    
	
	plot(disorder$DISORDER.Residues, z, xlab="# of Disorder Residues",ylab="Read Count Ratio")
	plot(disorder$DISORDER.., z, xlab="Disorder Residues Percentage (%)",ylab="Read Count Ratio")
	plot(disorder$Total.Residues, z, xlab="# of Total Residues",ylab="Read Count Ratio")

	cor.test(disorder$DISORDER.Residues, z)
	cor.test(disorder$DISORDER.., z)
	cor.test(disorder$Total.Residues, z)
	dev.off()	
}


#========================================================================             
# Firgure2D - Final
#
# BD-Fused proteins against Read Count
#
#========================================================================
convertNullMatrix(bm2_W[5:18,5:18],1, cellwidth = 10, cellheight = 10, filename1 = "bm2_W.pdf", filename2 = "bm2_W.null.pdf")  
convertFreqMatrix(bm2_Q[5:18,5:18],1, cellwidth = 10, cellheight = 10, filename = "bm2_Q.pdf" )  
InteractionScores(bm2_W[5:18,5:18],bm2_Q[5:18,5:18],1,10, cellwidth = 10, cellheight = 10, filename="is_bm2_Q.pdf")
PairInteractionScores(bm2_W[5:18,5:18],bm2_Q[5:18,5:18],1,10, cellwidth = 10, cellheight = 10, filename="pis_bm2_Q.pdf")
simple_plotLog10Matrix(NormalMixture2(bm2_Q[5:18,5:18]),1,10, cellwidth = 10, cellheight = 10, filename="nm2_bm2_Q.pdf")          
simple_plotLog10Matrix(LogNormalMixture(bm2_Q[5:18,5:18],threshold = 0.1),1,10, cellwidth = 10, cellheight = 10, filename="nm3_bm2_Q.pdf")          
InteractionScores(bm2_W[5:18,5:18],NormalMixture2(bm2_Q[5:18,5:18]),1,10, cellwidth = 10, cellheight = 10, filename="is2_bm2_Q.pdf")                   
is2_Fig2D = InteractionScores(bm2_W[5:18,5:18],NormalMixture2(bm2_Q[5:18,5:18]),1,10, cellwidth = 10, cellheight = 10, filename="is2_bm2_Q.pdf")   
is2acc3_Fig2D = AutoActivationCorrection3( is2_Fig2D )
simple_plotMatrix( is2acc3_Fig2D, 0, 9, filename="is2acc3_bm2_Q.pdf")   
NormalMixture2(bm2_Q[5:18,5:18])          
GammaMixture2(bm2_Q[5:18,5:18])        
DrawFig2NoiseFilter()



#========================================================================             
# Firgure2E - Final
#
# Show: Rawdata, Roth, Roth+Filter on 2 experiments
#
#========================================================================
	
	refine_bm2_W = RefineMatrix(bm2_W)
	refine_bm2_A = RefineMatrix(bm2_A)
	
	refine_bm3_W = RefineMatrix(bm3_W)
	refine_bm3_A = RefineMatrix(bm3_A)
	                                          
	is_refine_bm2_A = InteractionScores(refine_bm2_W,refine_bm2_A,1,9 )		     
	is_refine_bm3_A = InteractionScores(refine_bm3_W,refine_bm3_A,1,9 )	 

	is2_refine_bm2_A = InteractionScores(refine_bm2_W,NormalMixture2(refine_bm2_A),1,9 )		     
	is2_refine_bm3_A = InteractionScores(refine_bm3_W,NormalMixture2(refine_bm3_A),1,9 )  # 0.851969
	
	# Fig 2.f
	plot(m1,m2)
	abline( lm( as.vector(m1) ~ as.vector(m2) ) )
		
	pdf("Log2.refine_bm2_bm3.comp.pdf")
	CorMatrix(refine_bm2_A,refine_bm3_A, LogScale=TRUE)
	dev.off()
	    
	pdf("is.refine_bm2_bm3.comp.pdf")
	CorMatrix(is_refine_bm2_A,is_refine_bm3_A)
	dev.off()

	pdf("is2.refine_bm2_bm3.comp.pdf")
	CorMatrix(is2_refine_bm2_A,is2_refine_bm3_A)
	dev.off()	


#========================================================================             
# Firgure2F - Final
#
# Show: Boxplot over all experiments for rawdata, roth, roth + filter
#
#========================================================================

    LoadRefinedRothData()
    LoadRefinedRoth_IS()
    LoadRefinedRoth_IS_AAC()
	LoadRefinedRoth_IS_AAC2()
	LoadRefinedRoth_IS_AAC3()
    LoadRefinedRoth_IS2()
    LoadRefinedRoth_IS2_AAC()
    LoadRefinedRoth_IS2_AAC2()
    LoadRefinedRoth_IS2_AAC3()

    refine_selection_cor_list = c()
    is_cor_list = c()
    isacc_cor_list = c()
    isacc2_cor_list = c()
    isacc3_cor_list = c()
	is2_cor_list = c()
    is2acc_cor_list = c()
    is2acc2_cor_list = c()
	is2acc3_cor_list = c()
	
	index = 1
	for( i in 1:11 ){
		for( j in (i+1):12 ){
			refine_selection_cor_list = c( refine_selection_cor_list, CorMatrix(refine_selection_roth[[i]],refine_selection_roth[[j]],showPlot=FALSE) )
            is_cor_list = c( is_cor_list, CorMatrix(is_roth[[i]],is_roth[[j]],showPlot=FALSE) )
            isacc_cor_list = c( isacc_cor_list, CorMatrix(isaac_roth[[i]],isaac_roth[[j]],showPlot=FALSE) )
            isacc2_cor_list = c( isacc2_cor_list, CorMatrix(isaac2_roth[[i]],isaac2_roth[[j]],showPlot=FALSE) )
            isacc3_cor_list = c( isacc3_cor_list, CorMatrix(isaac3_roth[[i]],isaac3_roth[[j]],showPlot=FALSE) )
            is2_cor_list = c( is2_cor_list, CorMatrix(is2_roth[[i]],is2_roth[[j]],showPlot=FALSE) )
            is2acc_cor_list = c( is2acc_cor_list, CorMatrix(is2aac_roth[[i]],is2aac_roth[[j]],showPlot=FALSE) )
            is2acc2_cor_list = c( is2acc2_cor_list, CorMatrix(is2aac2_roth[[i]],is2aac2_roth[[j]],showPlot=FALSE) )
			is2acc3_cor_list = c( is2acc3_cor_list, CorMatrix(is2aac3_roth[[i]],is2aac3_roth[[j]],showPlot=FALSE) )
			index = index + 1
		}
	}

    # boxplot(refine_selection_cor_list, is_cor_list, isacc_cor_list,  isacc2_cor_list,  isacc3_cor_list, is2_cor_list, is2acc_cor_list, is2acc2_cor_list, is2acc3_cor_list)
    pdf("figure2f.correlation.boxplot.pdf")
    boxplot(refine_selection_cor_list, is_cor_list, is2_cor_list, ylab="Correlation", xlab="Cases", names=c("Raw read count","IS","IS with noise filter"))
    dev.off()


