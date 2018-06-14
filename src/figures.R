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
#Firgure2 - Screen data processing and effect of replicates (X75)   
#========================================================================
# Statistics on how many “useful reads” we have
# Pre-screening round with X75 used to identify auto-activators and toxicity.
# Generation of the interaction score by two normalisation steps
# Comparison -W normalisation or NOT
# Comparison to Roth analysis
# Statistics on sampling AND pair complexity and screen saturation with increasing sample number (it has to be clear what pair complexity is)
# Statistics on assay and sampling sensitivity with increasing sample number
# Correlation of conditions between experiments
# Compare single dropout VS multiple dropbout
# Compare tech rep to full replicate
# Is correlation affected by normalisation?
# How many of ORFeome can we distinguish with certain DNA length and analysis of current available ORFeomes (human; yeast; mouse etc).

# Read count   
if (TRUE){    
   fq_m1_W = 4042954
   fq_m1_A = 3180405
   fq_m1_RW = 4328429
   fq_m1_RA = 3698558
   fq_m2_W = 2323914
   fq_m2_A = 5760045
   fq_m2_Q = 2193771
   fq_m3_W = 2882223
   fq_m3_A = 2479944
   fq_m3_Q = 2863852
   fq_m3_QL = 2152438
   fq_m4_SW = 2562522
   fq_m4_S8A = 2076630
   fq_m4_S4A = 2321645
   fq_m4_SQ = 1437540
   fq_m5_SW = 7271213
   fq_m6_SW = 3783566
   fq_m6_SQ = 3443718
   fq_m7_SW = 4293368
   fq_m7_SQ = 3150897
   
   total_no_selection_reads = c( sum(fq_m1_W), sum(fq_m1_RW), sum(fq_m2_W), sum(fq_m3_W), sum(fq_m4_SW), sum(fq_m5_SW), sum(fq_m6_SW), sum(fq_m7_SW) )   
	#[1] 4042954 4328429 2323914 2882223 2562522 7271213 3783566 4293368
   total_selection_reads = c( sum(fq_m1_A), sum(fq_m1_RA), sum(fq_m2_A), sum(fq_m2_Q), sum(fq_m3_A), sum(fq_m3_Q), sum(fq_m3_QL), sum(fq_m4_S4A), sum(fq_m4_S8A), sum(fq_m4_SQ), sum(fq_m6_SQ), sum(fq_m7_SQ) )
    # [1] 3180405 3698558 5760045 2193771 2479944 2863852 2152438 2321645 2076630 1437540 3443718 3150897  
   total_reads = c( total_no_selection_reads, total_selection_reads )

   usable_no_selection_reads = c( sum(bm1_W), sum(bm1_RW), sum(bm2_W), sum(bm3_W), sum(bm4_SW), sum(bm5_SW), sum(bm6_SW), sum(bm7_SW) )   # 1878423 2111425  915100 1367370 1450453 1490700 1920425
    # [1] 1878423 2111425  915100 1367370 1450453 2649629 1490700 1920425
   usable_selection_reads = c( sum(bm1_A), sum(bm1_RA), sum(bm2_A), sum(bm2_Q), sum(bm3_A), sum(bm3_Q), sum(bm3_QL), sum(bm4_S4A), sum(bm4_S8A), sum(bm4_SQ), sum(bm6_SQ), sum(bm7_SQ) )
   # 1472503 1806918 1924195  779059  936030 1215139  713382 1145754 1094511  604169 1188438 1071150    
   usable_reads = c( usable_no_selection_reads, usable_selection_reads ) 
        
   plot( total_reads, usable_reads )     
   scatter.smooth(x=total_reads, y=usable_reads, main="Usable reads ~ Total reads")    
   reads.lm = lm(usable_reads ~ total_reads)
   plot(reads.lm)    
   summary(reads.lm)   
   
   mean( usable_reads/total_reads )  # 0.4246064

   pdf( "Usable_Reads.percent.pdf") 
   boxplot(usable_reads/total_reads*100.0)
   dev.off()

   ggplot(data.frame(total_reads=total_reads,usable_reads=usable_reads), aes(x=total_reads, y=usable_reads)) + geom_point(color='#2980B9', size = 4) + geom_smooth(method=lm, color='#2C3E50')

   reg.conf.intervals(total_reads, usable_reads)

   plot( usable_no_selection_reads/total_no_selection_reads)
   mean( usable_no_selection_reads/total_no_selection_reads) 

   plot( usable_selection_reads/total_selection_reads)
   mean( usable_selection_reads/total_selection_reads)       

   boxplot( usable_no_selection_reads/total_no_selection_reads, usable_selection_reads/total_selection_reads )


   # chart correlation
   sum(bm1_W), sum(bm1_RW), sum(bm2_W), sum(bm3_W), sum(bm4_SW), sum(bm5_SW), sum(bm6_SW), sum(bm7_SW)
   no_selection_data = data.frame( Exp1=as.vector(log2(bm1_W+1)),Exp2=as.vector(log2(bm1_RW+1)),Exp3=as.vector(log2(bm2_W+1)),Exp4=as.vector(log2(bm3_W+1)),Exp5=as.vector(log2(bm4_SW+1)),Exp6=as.vector(log2(bm5_SW+1)),Exp7=as.vector(log2(bm6_SW+1)),Exp8=as.vector(log2(bm7_SW+1)) )        
   pdf("Roth.no_selection.pdf",width = 600, height = 600); chart.Correlation(no_selection_data,method="pearson"); dev.off()     
   
   # Exp1=as.vector(log2(bm1_A+1)),Exp2=as.vector(log2(bm1_RA+1)) ## containing auto activators
   selection_data = data.frame( Exp3=as.vector(log2(bm2_A+1)),Exp4=as.vector(log2(bm2_Q+1)),Exp5=as.vector(log2(bm3_A+1)),Exp6=as.vector(log2(bm3_Q+1)),Exp7=as.vector(log2(bm3_QL+1)), Exp8=as.vector(log2(bm4_S4A+1)), Exp9=as.vector(log2(bm4_S8A+1)), Exp10=as.vector(log2(bm4_SQ+1)), Exp11=as.vector(log2(bm6_SQ+1)), Exp11=as.vector(log2(bm7_SQ+1)) )        
   pdf("Roth.selection.pdf",width = 600, height = 600); chart.Correlation(selection_data,method="pearson"); dev.off()   
    

   library(mixtools)
   wait = as.vector(bm2_A+1)      
   #wait = as.vector(log2(bm2_A+1))      
   #wait = faithful$wait
   mixmdl = normalmixEM(wait)
   plot(mixmdl,which=2)
   lines(density(wait), lty=2, lwd=2)    

   mixmdl= gammamixEM( as.vector(bm1_RA+1) )
   plot(mixmdl,which=2)
   lines(density(wait), lty=2, lwd=2)      
 


}



### http://kusanagi.hatenablog.jp/entry/2017/01/26/152922

dmixnormal<-function(x,lambda,m1,m2,s1,s2){
     y<-lambda*dnorm(x,m1,s1)+(1-lambda)*dnorm(x,m2,s2)
     y
}

dmixgamma<-function(x,lambda,a1,a2,b1,b2){
     y<-lambda*dgamma(x,a1,b1)+(1-lambda)*dgamma(x,a2,b2)
     y
}

set.seed(0)
dat<-c(rgamma(100,4,1),rgamma(300,10,.2))
hist(dat,breaks=20, main="",col="lightblue")
# Normalization performance     

dat = log2(as.vector(bm1_RA+1))

hist(dat,breaks=20, main="",col="lightblue",freq=F,ylim=c(0,.06),xlim=c(0,100))
x<-seq(0,20,.01)  
lambda= 0.6527965
a1=0.9218977 
a2=9.924739e-02  
b1=9.8971579 
b2=9.501456e+03
lines(x,dmixgamma(x,lambda,a1,a2,b1,b2),lwd=2,col=2)   

lines(x,dmixgamma(x,0.25,3.05,11.54,1/1.23,1/4.37),lwd=2,col=2)

mixmdl= gammamixEM( log2(as.vector(bm1_RA+1)),3 )             

$lambda
[1] 0.8893176 0.1106824

$gamma.pars
         comp.1       comp.2
alpha  0.484698    0.3890987
beta  80.101317 4655.0671008

x<-seq(0,20,.01)  
lambda= 0.8893176
a1=0.484698 
a2= 0.3890987  
b1=80.101317
b2=4655.0671008  
lines(x,dmixgamma(x,lambda,a1,a2,b1,b2),lwd=2,col=2)


NormalMixture2           
if (TRUE){    
	w_performance_summary = matrix(nrow=8,ncol=4)
	w_performance_summary[1,] = PPI_PerformanceCheck2(bm1_W,KnownPPI)
	w_performance_summary[2,] = PPI_PerformanceCheck2(bm1_RW,KnownPPI)
	w_performance_summary[3,] = PPI_PerformanceCheck2(bm2_W,KnownPPI)
	w_performance_summary[4,] = PPI_PerformanceCheck2(bm3_W,KnownPPI)
	w_performance_summary[5,] = PPI_PerformanceCheck2(bm4_SW,KnownPPI)
	w_performance_summary[6,] = PPI_PerformanceCheck2(bm5_SW,KnownPPI)  
	w_performance_summary[7,] = PPI_PerformanceCheck2(bm6_SW,KnownPPI)    
	w_performance_summary[8,] = PPI_PerformanceCheck2(bm7_SW,KnownPPI)
	
	wh_performance_summary = matrix(nrow=12,ncol=4)
	wh_performance_summary[1,] = PPI_PerformanceCheck2(bm1_A,KnownPPI)
	wh_performance_summary[2,] = PPI_PerformanceCheck2(bm1_RA,KnownPPI)
	wh_performance_summary[3,] = PPI_PerformanceCheck2(bm2_A,KnownPPI)
	wh_performance_summary[4,] = PPI_PerformanceCheck2(bm2_Q,KnownPPI)
	wh_performance_summary[5,] = PPI_PerformanceCheck2(bm3_A,KnownPPI)
	wh_performance_summary[6,] = PPI_PerformanceCheck2(bm3_Q,KnownPPI)
	wh_performance_summary[7,] = PPI_PerformanceCheck2(bm3_QL,KnownPPI)
	wh_performance_summary[8,] = PPI_PerformanceCheck2(bm4_S4A,KnownPPI)
	wh_performance_summary[9,] = PPI_PerformanceCheck2(bm4_S8A,KnownPPI)
	wh_performance_summary[10,] = PPI_PerformanceCheck2(bm4_SQ,KnownPPI)    
	wh_performance_summary[11,] = PPI_PerformanceCheck2(bm6_SQ,KnownPPI)    
	wh_performance_summary[12,] = PPI_PerformanceCheck2(bm7_SQ,KnownPPI)    
	
	wh_performance_summary2 = matrix(nrow=12,ncol=4) 
	wh_performance_summary2[1,] = PPI_PerformanceCheck2(NormalMixture2(bm1_A),KnownPPI)
	wh_performance_summary2[2,] = PPI_PerformanceCheck2(NormalMixture2(bm1_RA),KnownPPI)  
	wh_performance_summary2[3,] = PPI_PerformanceCheck2(NormalMixture2(bm2_A),KnownPPI)  
	wh_performance_summary2[4,] = PPI_PerformanceCheck2(NormalMixture2(bm2_Q),KnownPPI)  
	wh_performance_summary2[5,] = PPI_PerformanceCheck2(NormalMixture2(bm3_A),KnownPPI)  
	wh_performance_summary2[6,] = PPI_PerformanceCheck2(NormalMixture2(bm3_Q),KnownPPI)  
	wh_performance_summary2[7,] = PPI_PerformanceCheck2(NormalMixture2(bm3_QL),KnownPPI)  
	wh_performance_summary2[8,] = PPI_PerformanceCheck2(NormalMixture2(bm4_S4A),KnownPPI)  
	wh_performance_summary2[9,] = PPI_PerformanceCheck2(NormalMixture2(bm4_S8A),KnownPPI)  
	wh_performance_summary2[10,] = PPI_PerformanceCheck2(NormalMixture2(bm4_SQ),KnownPPI)  
	wh_performance_summary2[11,] = PPI_PerformanceCheck2(NormalMixture2(bm6_SQ),KnownPPI)  
	wh_performance_summary2[12,] = PPI_PerformanceCheck2(NormalMixture2(bm7_SQ),KnownPPI)  
	   
	is_performance_summary = matrix(nrow=12,ncol=4)
	is_performance_summary[1,] = PPI_PerformanceCheck2(InteractionScores(bm1_W,bm1_A,0.5),KnownPPI)
	is_performance_summary[2,] = PPI_PerformanceCheck2(InteractionScores(bm1_RW,bm1_RA,0.5),KnownPPI)
	is_performance_summary[3,] = PPI_PerformanceCheck2(InteractionScores(bm2_W,bm2_A,0.5),KnownPPI)
	is_performance_summary[4,] = PPI_PerformanceCheck2(InteractionScores(bm2_W,bm2_Q,0.5),KnownPPI)
	is_performance_summary[5,] = PPI_PerformanceCheck2(InteractionScores(bm3_W,bm3_A,0.5),KnownPPI)
	is_performance_summary[6,] = PPI_PerformanceCheck2(InteractionScores(bm3_W,bm3_Q,0.5),KnownPPI)
	is_performance_summary[7,] = PPI_PerformanceCheck2(InteractionScores(bm3_W,bm3_QL,0.5),KnownPPI)
	is_performance_summary[8,] = PPI_PerformanceCheck2(InteractionScores(bm4_SW,bm4_S4A,0.5),KnownPPI)
	is_performance_summary[9,] = PPI_PerformanceCheck2(InteractionScores(bm4_SW,bm4_S8A,0.5),KnownPPI)
	is_performance_summary[10,] = PPI_PerformanceCheck2(InteractionScores(bm4_SW,bm4_SQ,0.5),KnownPPI)    
	is_performance_summary[11,] = PPI_PerformanceCheck2(InteractionScores(bm6_SW,bm6_SQ,0.5),KnownPPI)    
	is_performance_summary[12,] = PPI_PerformanceCheck2(InteractionScores(bm7_SW,bm7_SQ,0.5),KnownPPI)    
	
	is_performance_summary2 = matrix(nrow=12,ncol=4)
	is_performance_summary2[1,] = PPI_PerformanceCheck2(InteractionScores(bm1_W,NormalMixture2(bm1_A),0.5),KnownPPI)
	is_performance_summary2[2,] = PPI_PerformanceCheck2(InteractionScores(bm1_RW,NormalMixture2(bm1_RA),0.5),KnownPPI)
	is_performance_summary2[3,] = PPI_PerformanceCheck2(InteractionScores(bm2_W,NormalMixture2(bm2_A),0.5),KnownPPI)
	is_performance_summary2[4,] = PPI_PerformanceCheck2(InteractionScores(bm2_W,NormalMixture2(bm2_Q),0.5),KnownPPI)
	is_performance_summary2[5,] = PPI_PerformanceCheck2(InteractionScores(bm3_W,NormalMixture2(bm3_A),0.5),KnownPPI)
	is_performance_summary2[6,] = PPI_PerformanceCheck2(InteractionScores(bm3_W,NormalMixture2(bm3_Q),0.5),KnownPPI)
	is_performance_summary2[7,] = PPI_PerformanceCheck2(InteractionScores(bm3_W,NormalMixture2(bm3_QL),0.5),KnownPPI)
	is_performance_summary2[8,] = PPI_PerformanceCheck2(InteractionScores(bm4_SW,NormalMixture2(bm4_S4A),0.5),KnownPPI)
	is_performance_summary2[9,] = PPI_PerformanceCheck2(InteractionScores(bm4_SW,NormalMixture2(bm4_S8A),0.5),KnownPPI)
	is_performance_summary2[10,] = PPI_PerformanceCheck2(InteractionScores(bm4_SW,NormalMixture2(bm4_SQ),0.5),KnownPPI)    
	is_performance_summary2[11,] = PPI_PerformanceCheck2(InteractionScores(bm6_SW,NormalMixture2(bm6_SQ),0.5),KnownPPI)    
	is_performance_summary2[12,] = PPI_PerformanceCheck2(InteractionScores(bm7_SW,NormalMixture2(bm7_SQ),0.5),KnownPPI)	
	
	pis_1 = PairInteractionScores(bm1_W,bm1_A,0.5)
	pis_2 = PairInteractionScores(bm1_RW,bm1_RA,0.5)
	pis_3 = PairInteractionScores(bm2_W,bm2_A,0.5)
	pis_4 = PairInteractionScores(bm2_W,bm2_Q,0.5)
	pis_5 = PairInteractionScores(bm3_W,bm3_A,0.5)
	pis_6 = PairInteractionScores(bm3_W,bm3_Q,0.5)
	pis_7 = PairInteractionScores(bm3_W,bm3_QL,0.5)
	pis_8 = PairInteractionScores(bm4_SW,bm4_S4A,0.5)
	pis_9 = PairInteractionScores(bm4_SW,bm4_S8A,0.5)
	pis_10 = PairInteractionScores(bm4_SW,bm4_SQ,0.5)
	pis_11 = PairInteractionScores(bm6_SW,bm6_SQ,0.5)    
	pis_12 = PairInteractionScores(bm7_SW,bm7_SQ,0.5)
	pis_sum = pis_1+pis_2+pis_3+pis_4+pis_5+pis_6+pis_7+pis_8+pis_9+pis_10+pis_11+pis_12
	
	pis_performance_summary = matrix(nrow=12,ncol=4)
	pis_performance_summary[1,] = PPI_PerformanceCheck2(PairInteractionScores(bm1_W,bm1_A,0.5),KnownPPI)
	pis_performance_summary[2,] = PPI_PerformanceCheck2(PairInteractionScores(bm1_RW,bm1_RA,0.5),KnownPPI)
	pis_performance_summary[3,] = PPI_PerformanceCheck2(PairInteractionScores(bm2_W,bm2_A,0.5),KnownPPI)
	pis_performance_summary[4,] = PPI_PerformanceCheck2(PairInteractionScores(bm2_W,bm2_Q,0.5),KnownPPI)
	pis_performance_summary[5,] = PPI_PerformanceCheck2(PairInteractionScores(bm3_W,bm3_A,0.5),KnownPPI)
	pis_performance_summary[6,] = PPI_PerformanceCheck2(PairInteractionScores(bm3_W,bm3_Q,0.5),KnownPPI)
	pis_performance_summary[7,] = PPI_PerformanceCheck2(PairInteractionScores(bm3_W,bm3_QL,0.5),KnownPPI)
	pis_performance_summary[8,] = PPI_PerformanceCheck2(PairInteractionScores(bm4_SW,bm4_S4A,0.5),KnownPPI)
	pis_performance_summary[9,] = PPI_PerformanceCheck2(PairInteractionScores(bm4_SW,bm4_S8A,0.5),KnownPPI)
	pis_performance_summary[10,] = PPI_PerformanceCheck2(PairInteractionScores(bm4_SW,bm4_SQ,0.5),KnownPPI)    
	pis_performance_summary[11,] = PPI_PerformanceCheck2(PairInteractionScores(bm6_SW,bm6_SQ,0.5),KnownPPI)    
	pis_performance_summary[12,] = PPI_PerformanceCheck2(PairInteractionScores(bm7_SW,bm7_SQ,0.5),KnownPPI)
	
	nis_performance_summary = matrix(nrow=12,ncol=4)
	nis_performance_summary[1,] = PPI_PerformanceCheck2(NewInteractionScores(bm1_W,bm1_A,0.5),KnownPPI)
	nis_performance_summary[2,] = PPI_PerformanceCheck2(NewInteractionScores(bm1_RW,bm1_RA,0.5),KnownPPI)
	nis_performance_summary[3,] = PPI_PerformanceCheck2(NewInteractionScores(bm2_W,bm2_A,0.5),KnownPPI)
	nis_performance_summary[4,] = PPI_PerformanceCheck2(NewInteractionScores(bm2_W,bm2_Q,0.5),KnownPPI)
	nis_performance_summary[5,] = PPI_PerformanceCheck2(NewInteractionScores(bm3_W,bm3_A,0.5),KnownPPI)
	nis_performance_summary[6,] = PPI_PerformanceCheck2(NewInteractionScores(bm3_W,bm3_Q,0.5),KnownPPI)
	nis_performance_summary[7,] = PPI_PerformanceCheck2(NewInteractionScores(bm3_W,bm3_QL,0.5),KnownPPI)
	nis_performance_summary[8,] = PPI_PerformanceCheck2(NewInteractionScores(bm4_SW,bm4_S4A,0.5),KnownPPI)
	nis_performance_summary[9,] = PPI_PerformanceCheck2(NewInteractionScores(bm4_SW,bm4_S8A,0.5),KnownPPI)
	nis_performance_summary[10,] = PPI_PerformanceCheck2(NewInteractionScores(bm4_SW,bm4_SQ,0.5),KnownPPI)    
	nis_performance_summary[11,] = PPI_PerformanceCheck2(NewInteractionScores(bm6_SW,bm6_SQ,0.5),KnownPPI)    
	nis_performance_summary[12,] = PPI_PerformanceCheck2(NewInteractionScores(bm7_SW,bm7_SQ,0.5),KnownPPI)
	
	cis_performance_summary = matrix(nrow=12,ncol=4)
	cis_performance_summary[1,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm1_W,bm1_A,0.5),KnownPPI)
	cis_performance_summary[2,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm1_RW,bm1_RA,0.5),KnownPPI)
	cis_performance_summary[3,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm2_W,bm2_A,0.5),KnownPPI)
	cis_performance_summary[4,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm2_W,bm2_Q,0.5),KnownPPI)
	cis_performance_summary[5,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm3_W,bm3_A,0.5),KnownPPI)
	cis_performance_summary[6,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm3_W,bm3_Q,0.5),KnownPPI)
	cis_performance_summary[7,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm3_W,bm3_QL,0.5),KnownPPI)
	cis_performance_summary[8,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm4_SW,bm4_S4A,0.5),KnownPPI)
	cis_performance_summary[9,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm4_SW,bm4_S8A,0.5),KnownPPI)
	cis_performance_summary[10,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm4_SW,bm4_SQ,0.5),KnownPPI)    
	cis_performance_summary[11,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm6_SW,bm6_SQ,0.5),KnownPPI)    
	cis_performance_summary[12,] = PPI_PerformanceCheck2(ColumnWiseInteractionScores(bm7_SW,bm7_SQ,0.5),KnownPPI)	
                                                                                                                        
	# AUC 
	pdf("Roth.PerformanceTest.AUC.pdf")
    boxplot( w_performance_summary[,1], wh_performance_summary[,1], wh_performance_summary2[,1], is_performance_summary[,1], is_performance_summary2[,1],pis_performance_summary[,1],nis_performance_summary[,1],cis_performance_summary[,1])  
	dev.off()
	mean(is_performance_summary[,1])
    #[1] 0.5209035    0.5112899 (RefineMatrix)
	mean(pis_performance_summary[,1])
	#[1] 0.6311815    0.6101566 (RefineMatrix)  
	mean(nis_performance_summary[,1])
	#[1] 0.6111054    0.6034389 (RefineMatrix)  
	mean(cis_performance_summary[,1])
	#[1] 0.5428511    0.5312352 (RefineMatrix)  

	# T-Test   
	pdf("Roth.PerformanceTest.TTest.pdf")  
	boxplot( w_performance_summary[,2], wh_performance_summary[,2], wh_performance_summary2[,2], is_performance_summary[,2], is_performance_summary2[,2],pis_performance_summary[,2],nis_performance_summary[,2],cis_performance_summary[,2]) 
	dev.off()                                                                                   
}



#=====================================
# Fig2f.v1
#=====================================

rrbm1_A = NormalMixture(refine_bm1_A)                                    
i1_A = InteractionScores(refine_bm1_W, rrbm1_A, 1, 7)

rrbm2_A = NormalMixture(refine_bm2_A)  
i2_A = InteractionScores(refine_bm2_W, rrbm2_A, 1, 7)

rrbm2_Q = NormalMixture(refine_bm2_Q)  
i2_Q = InteractionScores(refine_bm2_W, rrbm2_Q, 1, 7)


rrbm2_A = NormalMixture(refine_bm2_A)  
i2_A = InteractionScores(refine_bm2_W, rrbm2_A, 1, 7)


## GammaMixture Models                                                 
gm1_A = GammaMixture(refine_bm1_A)
gm2_A = GammaMixture(refine_bm2_A)
gm2_Q = GammaMixture(refine_bm2_Q)
gm3_A = GammaMixture(refine_bm3_A)
gm3_Q = GammaMixture(refine_bm3_Q)
gm3_QL = GammaMixture(refine_bm3_QL)   
gm4_S4A = GammaMixture(refine_bm4_S4A)
gm4_S8A = GammaMixture(refine_bm4_S8A)
gm4_SQ = GammaMixture(refine_bm4_SQ)

simple_plotLog2Matrix(gm1_A+gm2_A+gm2_Q+gm3_A+gm3_Q+gm3_QL+gm4_S4A+gm4_S8A+gm4_SQ,1,7)
gmm_sum = gm1_A+gm2_A+gm2_Q+gm3_A+gm3_Q+gm3_QL+gm4_S4A+gm4_S8A+gm4_SQ       
refine_HIPPIE = RefineMatrix(HIPPIE)      
PPI_PerformanceCheck2(gmm_sum,refine_HIPPIE)
 
pdf("gmm_overlap.pdf")  
gmm_overlap = (gm1_A>1)+(gm2_A>1)+(gm2_Q>1)+(gm3_A>1)+(gm3_Q>1)+(gm3_QL>1)+(gm4_S4A>1)+(gm4_S8A>1)+(gm4_SQ>1)
simple_plotMatrix(gmm_overlap,0,7)
dev.off()                           
PPI_PerformanceCheck2(gmm_overlap,refine_HIPPIE)


 
pdf("cva1.refined.pdf")  
cva1_m3=drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/BFG_Y2H/Yachie_Petsalaki_Data_S3/output/CVA_1_+3AT.tsv.refined.txt",0.5,7)   
dev.off()
cva1_Row_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,26,28,29,30,31,32,33,35,36,37,39,41,42,45,46,47,48,50,51,52,53,54,56,57,58,59,60,61,62,64,66,67,70,71,72,73,74 );
cva1_Col_Index = c( 1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,21,24,26,28,29,31,32,33,35,36,37,38,41,42,45,46,47,48,50,51,52,53,56,57,58,60,61,62,63,64,66,67,68,70,71,72,74 );    
Null_cva1 = cva1_m1 * 0 #matrix(0,nrow=78,ncol=78)
cva1_Row_Not_Index = (1:78)[-cva1_Row_Index]
cva1_Col_Not_Index = (1:78)[-cva1_Col_Index]   
Null_cva1[cva1_Row_Not_Index,] = 1     
Null_cva1[,cva1_Col_Not_Index] = 1  
simple_plotMatrix(Null_cva1, 1, 7)


pdf("bm5.SW.pdf")    
bm5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt.new",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
dev.off()       

                                                                                       
#=====================================
# Fig.WorkFlow.ai - Fig 2
#=====================================
convertNullMatrix(bm2_W[5:18,5:18],1, cellwidth = 10, cellheight = 10, filename1 = "bm2_W.pdf", filename2 = "bm2_W.null.pdf")  
convertFreqMatrix(bm2_Q[5:18,5:18],1, cellwidth = 10, cellheight = 10, filename = "bm2_Q.pdf" )  
InteractionScores(bm2_W[5:18,5:18],bm2_Q[5:18,5:18],1,10, cellwidth = 10, cellheight = 10, filename="is_bm2_Q.pdf")
PairInteractionScores(bm2_W[5:18,5:18],bm2_Q[5:18,5:18],1,10, cellwidth = 10, cellheight = 10, filename="pis_bm2_Q.pdf")
simple_plotLog10Matrix(NormalMixture2(bm2_Q[5:18,5:18]),1,10, cellwidth = 10, cellheight = 10, filename="nm2_bm2_Q.pdf")          
simple_plotLog10Matrix(LogNormalMixture(bm2_Q[5:18,5:18],threshold = 0.1),1,10, cellwidth = 10, cellheight = 10, filename="nm3_bm2_Q.pdf")          
InteractionScores(bm2_W[5:18,5:18],NormalMixture2(bm2_Q[5:18,5:18]),1,10, cellwidth = 10, cellheight = 10, filename="is2_bm2_Q.pdf")                   
NormalMixture2(bm2_Q[5:18,5:18])          
GammaMixture2(bm2_Q[5:18,5:18])        
DrawFig2NoiseFilter()


#=====================================
# Fig4.v1
#=====================================

# CorMatrix(nibP170_4_SQ,nibP170_5_SQ)
# 0.835768853281392               

xx_nibP170_4_SQ = nibP170_4_SQ * t(nibP170_4_SQ)
simple_plotLog2Matrix(xx_nibP170_4_SQ,1,4)
simple_plotMatrix(xx_nibP170_4_SQ>5,0,4)

xx_nibP170_5_SQ = nibP170_5_SQ * t(nibP170_5_SQ)
simple_plotLog2Matrix(xx_nibP170_5_SQ,1,4)
simple_plotMatrix(xx_nibP170_5_SQ>5,0,4)

simple_plotMatrix((xx_nibP170_4_SQ>5)+(xx_nibP170_5_SQ>5),0,4)

pdf("163x.core_ppi.pdf")
simple_plotMatrix(((xx_nibP170_4_SQ>5)+(xx_nibP170_5_SQ>5))>1,0,4)
dev.off()                                                  


P170_is_sum = (ibP170_1_Q>3)+(ibP170_2_Q>3)+(ibP170_3_Q>3)+(ibP170_4_SQ>3)+(ibP170_5_SQ>3)   
simple_plotMatrix(P170_is_sum[1:58,1:58],0,7, filename = "P170_1-1.pdf" )  
simple_plotMatrix(P170_is_sum[59:116,1:58],0,7, filename = "P170_2-1.pdf" )  
simple_plotMatrix(P170_is_sum[117:173,1:58],0,7, filename = "P170_3-1.pdf" )        
                                                                          
simple_plotMatrix(P170_is_sum[1:58,59:116],0,7, filename = "P170_1-2.pdf" )  
simple_plotMatrix(P170_is_sum[59:116,59:116],0,7, filename = "P170_2-2.pdf" )  
simple_plotMatrix(P170_is_sum[117:173,59:116],0,7, filename = "P170_3-2.pdf" )

simple_plotMatrix(P170_is_sum[1:58,117:173],0,7, filename = "P170_1-3.pdf" )  
simple_plotMatrix(P170_is_sum[59:116,117:173],0,7, filename = "P170_2-3.pdf" )  
simple_plotMatrix(P170_is_sum[117:173,117:173],0,7, filename = "P170_3-3.pdf" )


# EGFR
if (TRUE){                 
	set.seed(43)
	gm_EGFR1_4A = GammaMixture(bEGFR1_4A)      
	ibEGFR1_4A= InteractionScores(bEGFR1_W,gm_EGFR1_4A,1,4)
	#ibEGFR1_4A= InteractionScores(bEGFR1_W,bEGFR1_4A,1,4)
	#nibEGFR1_4A = NewInteractionScores(bEGFR1_W,gm_EGFR1_4A,1,4)
	
	gm_EGFR1_8A = GammaMixture(bEGFR1_8A)      
	ibEGFR1_8A= InteractionScores(bEGFR1_W,gm_EGFR1_8A,1,4)
	#ibEGFR1_8A= InteractionScores(bEGFR1_W,bEGFR1_8A,1,4)
	#nibEGFR1_8A = NewInteractionScores(bEGFR1_W,gm_EGFR1_8A,1,4)
	
	gm_EGFR1_Q = GammaMixture(bEGFR1_Q)      
	ibEGFR1_Q= InteractionScores(bEGFR1_W,gm_EGFR1_Q,1,4)     
	#ibEGFR1_Q= InteractionScores(bEGFR1_W,bEGFR1_Q,1,4) 
	#nibEGFR1_Q = NewInteractionScores(bEGFR1_W,gm_EGFR1_Q,1,4)   
	
	EGFR_overlap = (nibEGFR1_4A>3)+(nibEGFR1_8A>3)+(nibEGFR1_Q>3)   
	
	EGFR_overlap = (ibEGFR1_4A>3)+(ibEGFR1_8A>3)+(ibEGFR1_Q>3) 
	EGFR_is_sum = (ibEGFR1_4A)+(ibEGFR1_8A)+(ibEGFR1_Q) 

	# remove Auto-activator   
	# (rowSums(EGFR_overlap))[ rowSums(EGFR_overlap) >50 ]
	# 109_PLSCR1   117_PTK2B     124_PXN  177_ZNF259    178_DOK1 202_APBB1IP      22_CRK    38_EPS15      42_FOS     44_GAB1     9_BCAR1 
    #     418          59         218          91          82         244          62          68         356          83         163
    EGFR_row_index = setdiff( 1:182, (1:182)[ rowSums(EGFR_overlap) >50 ] )
    simple_plotMatrix(EGFR_overlap[EGFR_row_index,],0,4)    
                
	(ibEGFR1_4A[141,])[ibEGFR1_4A[141,]>=3] 
	(ibEGFR1_4A[,141])[ibEGFR1_4A[,141]>=3]   
	(ibEGFR1_8A[141,])[ibEGFR1_8A[141,]>=3]   
	(ibEGFR1_8A[,141])[ibEGFR1_8A[,141]>=3]   
	(ibEGFR1_Q[141,])[ibEGFR1_Q[141,]>=3]   
	(ibEGFR1_Q[,141])[ibEGFR1_Q[,141]>=3]   
	
	# KRAS binder
	(EGFR_overlap[141,])[EGFR_overlap[141,]>=1]              
	#> (EGFR_overlap[141,])[EGFR_overlap[141,]>=1]
	#   X128_RALB(o)    X180_RIN2(o)    X182_RGL1(o)    X187_BRAF(o) X198_RAPGEF4  X201_RASSF5(o) 
	#           1            1            1            1            1            1  
	#	X128_RALB  X148_SHOC2   X180_RIN2   X182_RGL1   X187_BRAF X201_RASSF5   X84_MAPK8 
    #      1           1           1           1           1           2           1
   	(EGFR_overlap[,141])[EGFR_overlap[,141]>=1]
 	#109_PLSCR1     124_PXN   159_SPRY2   176_YWHAB    178_DOK1 198_RAPGEF4 202_APBB1IP    38_EPS15      42_FOS     9_BCAR1      TOPBP1 
    #      2           2           1           1           1           1           1           1           2           2           1    
	#109_PLSCR1    124_PXN  176_YWHAB   178_DOK1     42_FOS    9_BCAR1 
	#         3          2          1          1          3          2
	
	simple_plotMatrix( EGFR_overlap, 0, 4  ) 
	EGFR_symmetric = ( EGFR_overlap >= 1 ) * ( t(EGFR_overlap) >= 1 )  
	EGFR_non_symmetric = (-1*(EGFR_symmetric-1))
	
	simple_plotMatrix( EGFR_symmetric, 0, 4 )
	simple_plotMatrix( EGFR_non_symmetric, 0, 4 )
	                   
	simple_plotMatrix( EGFR_overlap*EGFR_symmetric, 0, 4  ) 
	simple_plotMatrix( EGFR_overlap*EGFR_non_symmetric, 0, 4  )
	
	#x = hist(EGFR_overlap*EGFR_symmetric)     
	x = hist((EGFR_overlap*EGFR_symmetric)[EGFR_row_index,])  
	30269     0     0     0    86     0     0     0     0    15     0     0     0     0    24
	(24/(86+15+24))
	y = hist((EGFR_overlap*EGFR_non_symmetric)[EGFR_row_index,])
	29868     0     0     0   451     0     0     0     0    42     0     0     0     0    33
	(33/(451+42+33))	                                     
		
		
	## overlap   
	simple_plotMatrix(EGFR_overlap[1:60,1:60],0,7, filename = "EGFR_1-1.pdf" )
	simple_plotMatrix(EGFR_overlap[61:121,1:60],0,7, filename = "EGFR_2-1.pdf")
	simple_plotMatrix(EGFR_overlap[122:182,1:60],0,7, filename = "EGFR_3-1.pdf")  
	
	simple_plotMatrix(EGFR_overlap[1:60,61:121],0,7, filename = "EGFR_1-2.pdf")
	simple_plotMatrix(EGFR_overlap[61:121,61:121],0,7, filename = "EGFR_2-2.pdf")
	simple_plotMatrix(EGFR_overlap[122:182,61:121],0,7, filename = "EGFR_3-2.pdf")  
	
	simple_plotMatrix(EGFR_overlap[1:60,122:182],0,7, filename = "EGFR_1-3.pdf")
	simple_plotMatrix(EGFR_overlap[61:121,122:182],0,7, filename = "EGFR_2-3.pdf")
	simple_plotMatrix(EGFR_overlap[122:182,122:182],0,7, filename = "EGFR_3-3.pdf")     
	
	# Save file
	melted_EGFR_overlap = melt(EGFR_overlap)
	write.table(melted_EGFR_overlap,file="/Users/jyang/melted_EGFR_overlap.IS.txt",sep="\t")
	
	          
	## is-sum         
	EGFR_is_sum = (ibEGFR1_4A)+(ibEGFR1_8A)+(ibEGFR1_Q)   
	simple_plotMatrix(EGFR_is_sum[1:60,1:60],0,7, filename = "isEGFR_1-1.pdf" )
	simple_plotMatrix(EGFR_is_sum[61:121,1:60],0,7, filename = "isEGFR_2-1.pdf")
	simple_plotMatrix(EGFR_is_sum[122:182,1:60],0,7, filename = "isEGFR_3-1.pdf")  
	
	simple_plotMatrix(EGFR_is_sum[1:60,61:121],0,7, filename = "isEGFR_1-2.pdf")
	simple_plotMatrix(EGFR_is_sum[61:121,61:121],0,7, filename = "isEGFR_2-2.pdf")
	simple_plotMatrix(EGFR_is_sum[122:182,61:121],0,7, filename = "isEGFR_3-2.pdf")  
	
	simple_plotMatrix(EGFR_is_sum[1:60,122:182],0,7, filename = "isEGFR_1-3.pdf")
	simple_plotMatrix(EGFR_is_sum[61:121,122:182],0,7, filename = "isEGFR_2-3.pdf")
	simple_plotMatrix(EGFR_is_sum[122:182,122:182],0,7, filename = "isEGFR_3-3.pdf")  
	
	simple_plotMatrix(((nibEGFR1_4A>3)+(nibEGFR1_8A>3)+(nibEGFR1_Q>3))>1,0,4)    
	
	# NIB
	simple_plotMatrix(((nibEGFR1_4A>3)),0,4)
	simple_plotMatrix(((nibEGFR1_8A>3)),0,4)
	simple_plotMatrix(((nibEGFR1_QD>3)),0,4)

	# IB
	simple_plotMatrix(((ibEGFR1_4A>3)),0,4)
	simple_plotMatrix(((ibEGFR1_8A>3)),0,4)
	simple_plotMatrix(((ibEGFR1_QD>3)),0,4)        
}
   

### P170 to cytoscape

row_names = rownames(P170_is_sum)
col_names = colnames(P170_is_sum)
for (i in 1:173) {
	for (j in 1:173) {
		if( P170_is_sum[i,j] >=2 ){
			print ( sprintf( "%s %s %d", row_names[i], col_names[j], P170_is_sum[i,j] ) )
		}
	}
}
     


# PAR Network          
REF_INDEX = c(92,93,120)
PAR_INDEX = c(17,18,19,103,104,105,121,122,123,124,125,126,127,128,135,136,137,138,171,172,173)
simple_plotMatrix(P170_is_sum[PAR_INDEX,PAR_INDEX],0,7 )  
simple_plotMatrix(P170_is_sum[PAR_INDEX,PAR_INDEX]+pis170_SIG[PAR_INDEX,PAR_INDEX],0,7,filename="P170.PAR.pdf",cellwidth=10,cellheight=10 )  
Matrix2List(P170_is_sum[PAR_INDEX,PAR_INDEX]+pis170_SIG[PAR_INDEX,PAR_INDEX],filename="ParSubNetwork.txt")
Matrix2List(P170_is_sum+pis170_SIG,filename="P170_cutoff.3.txt",cutoff=3)   
Matrix2List(P170_is_sum+pis170_SIG,filename="P170_cutoff.5.txt",cutoff=5)  
Matrix2List(P170_is_sum+pis170_SIG,filename="P170_cutoff.6.txt",cutoff=6)  
Matrix2List(P170_is_sum+pis170_SIG,filename="P170_cutoff.7.txt",cutoff=7)  
Matrix2List(P170_is_sum+pis170_SIG,filename="P170_cutoff.8.txt",cutoff=8)  
Matrix2List(P170_is_sum+pis170_SIG,filename="P170_cutoff.9.txt",cutoff=9)  

ADAPTOR_INDEX = c(14,15,94,159)                  
simple_plotMatrix(P170_is_sum[ADAPTOR_INDEX,ADAPTOR_INDEX],0,7 )  
 
CYTOSKELETON_INDEX = c(8,9,164,165)
simple_plotMatrix(P170_is_sum[CYTOSKELETON_INDEX,CYTOSKELETON_INDEX],0,7 )  

MAP_INDEX = c(67,96,97,98,102)
simple_plotMatrix(P170_is_sum[MAP_INDEX,MAP_INDEX],0,7 )       


MOTOR_INDEX = c(48,49,50,51,52,68,69,70,71,72,73,74,75,76,77,78,80,81,82,83,84,85,86,87,88,89,90,107,108)
simple_plotMatrix(P170_is_sum[MOTOR_INDEX,MOTOR_INDEX],0,7 )

plusTIP_INDEX = c(10,11,12,13,21,22,23,24,25,27,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,53,54,55,58,59,60,62,63,64,65,66,79,91,95,106,109,110,111,112,113,114,115,116,117,118,119,129,130,131,132,133,134,139,140,141,143,144,145,146,147,148,150,151,152,153,154,155,156,157,158,161,162,166,167,168,169,170)
simple_plotMatrix(P170_is_sum[plusTIP_INDEX,plusTIP_INDEX],0,7 ) 

CONTROL_INDEX = c(16,20,26,28,47,61,92,93,99,100,101,120,142,149,163)
simple_plotMatrix(P170_is_sum[CONTROL_INDEX,CONTROL_INDEX],0,7 ) 

DOMAIN_INDEX = c(79,74,107,108,17,18,19,103,104,105,121,122,123,124,125,126,127,128,135,136,137,138,171,172,173)
simple_plotMatrix(P170_is_sum[DOMAIN_INDEX,DOMAIN_INDEX],0,7 )
                                                    

## R75 InsteractionScores
is1A = InteractionScores(bm1_W,bm1_A,1,4)
is2A = InteractionScores(bm2_W,bm2_A,1,4)
is2Q = InteractionScores(bm2_W,bm2_Q,1,4)
is3A = InteractionScores(bm3_W,bm3_A,1,4)   
is3Q = InteractionScores(bm3_W,bm3_Q,1,4) 
is4SA4 = InteractionScores(bm4_SW,bm4_S4A,1,4)   
is4SA8 = InteractionScores(bm4_SW,bm4_S8A,1,4)   
is4SQ = InteractionScores(bm4_SW,bm4_SQ,1,4)   
simple_plotMatrix( is1A + is2A + is2Q + is3A + is3Q + is4SA4 + is4SA8 + is4SQ, 0, 7 )
simple_plotMatrix( (is1A>3) + (is2A>3) + (is2Q>3) + (is3A>3) + (is3Q>3) + (is4SA4>3) + (is4SA8>3) + (is4SQ>3), 0, 7, filename="isRoth75.pdf" ,cellwidth = 10,cellheight = 10)       
simple_plotMatrix( ((is1A>3) + (is2A>3) + (is2Q>3) + (is3A>3) + (is3Q>3) + (is4SA4>3) + (is4SA8>3) + (is4SQ>3))>3, 0, 7 )

## R75 pair-wise normalization
pis1A = PairInteractionScores(bm1_W,bm1_A,1,4)
pis2A = PairInteractionScores(bm2_W,bm2_A,1,4)
pis2Q = PairInteractionScores(bm2_W,bm2_Q,1,4)
pis3A = PairInteractionScores(bm3_W,bm3_A,1,4)   
pis3Q = PairInteractionScores(bm3_W,bm3_Q,1,4) 
pis4SA4 = PairInteractionScores(bm4_SW,bm4_S4A,1,4)   
pis4SA8 = PairInteractionScores(bm4_SW,bm4_S8A,1,4)   
pis4SQ = PairInteractionScores(bm4_SW,bm4_SQ,1,4)   
simple_PlotMatrix( pis1A + pis2A + pis2Q + pis3A + pis3Q + pis4A4 + pis4A8 + pis4SQ, 0, 7 )
simple_plotMatrix( (pis1A>3) + (pis2A>3) + (pis2Q>3) + (pis3A>3) + (pis3Q>3) + (pis4SA4>3) + (pis4SA8>3) + (pis4SQ>3), 0, 7, filename="pisRoth75.pdf",cellwidth = 10,cellheight = 10)       
simple_plotMatrix( ((pis1A>3) + (pis2A>3) + (pis2Q>3) + (pis3A>3) + (pis3Q>3) + (pis4SA4>3) + (pis4SA8>3) + (pis4SQ>3))>3, 0, 7 )



pis_EGFR1_4A = PairInteractionScores(bEGFR1_W,gm_EGFR1_4A,1,4)                
pis_EGFR1_8A = PairInteractionScores(bEGFR1_W,gm_EGFR1_8A,1,4)                
pis_EGFR1_Q = PairInteractionScores(bEGFR1_W,gm_EGFR1_Q,1,4)  
pis_EGFR1_SIG = (pis_EGFR1_4A>3)+(pis_EGFR1_8A>3)+(pis_EGFR1_Q>3)              
                                                                             

simple_plotMatrix(pis_EGFR1_SIG[1:60,1:60],0,7, filename = "pisEGFR_1-1.pdf" )
simple_plotMatrix(pis_EGFR1_SIG[61:121,1:60],0,7, filename = "pisEGFR_2-1.pdf")
simple_plotMatrix(pis_EGFR1_SIG[122:182,1:60],0,7, filename = "pisEGFR_3-1.pdf")  

simple_plotMatrix(pis_EGFR1_SIG[1:60,61:121],0,7, filename = "pisEGFR_1-2.pdf")
simple_plotMatrix(pis_EGFR1_SIG[61:121,61:121],0,7, filename = "pisEGFR_2-2.pdf")
simple_plotMatrix(pis_EGFR1_SIG[122:182,61:121],0,7, filename = "pisEGFR_3-2.pdf")  

simple_plotMatrix(pis_EGFR1_SIG[1:60,122:182],0,7, filename = "pisEGFR_1-3.pdf")
simple_plotMatrix(pis_EGFR1_SIG[61:121,122:182],0,7, filename = "pisEGFR_2-3.pdf")
simple_plotMatrix(pis_EGFR1_SIG[122:182,122:182],0,7, filename = "pisEGFR_3-3.pdf")




## R170 pair-wise normalization  
pis170_1Q = PairInteractionScores(bP170_1_W,bP170_1_Q,1,4) 
pis170_2Q = PairInteractionScores(bP170_2_W,bP170_2_Q,1,4) 
pis170_3Q = PairInteractionScores(bP170_3_W,bP170_3_Q,1,4) 
pis170_4Q = PairInteractionScores(bP170_4_SW,bP170_4_SQ,1,4) 
pis170_5Q = PairInteractionScores(bP170_5_SW,bP170_5_SQ,1,4) 
pis170_SIG = (pis170_1Q>3)+(pis170_2Q>3)+(pis170_3Q>3)+(pis170_4Q>3)+(pis170_5Q>3)                    
simple_plotMatrix(pis170_SIG,0,4, filename = NA )     
                                                        
simple_plotMatrix(pis170_SIG[1:58,1:58],0,7, filename = "pisP170_1-1.pdf" )  
simple_plotMatrix(pis170_SIG[59:116,1:58],0,7, filename = "pisP170_2-1.pdf" )  
simple_plotMatrix(pis170_SIG[117:173,1:58],0,7, filename = "pisP170_3-1.pdf" )        
                                                                          
simple_plotMatrix(pis170_SIG[1:58,59:116],0,7, filename = "pisP170_1-2.pdf" )  
simple_plotMatrix(pis170_SIG[59:116,59:116],0,7, filename = "pisP170_2-2.pdf" )  
simple_plotMatrix(pis170_SIG[117:173,59:116],0,7, filename = "pisP170_3-2.pdf" )

simple_plotMatrix(pis170_SIG[1:58,117:173],0,7, filename = "pisP170_1-3.pdf" )  
simple_plotMatrix(pis170_SIG[59:116,117:173],0,7, filename = "pisP170_2-3.pdf" )  
simple_plotMatrix(pis170_SIG[117:173,117:173],0,7, filename = "pisP170_3-3.pdf" )



simple_plotMatrix(pis170_SIG[PAR_INDEX,PAR_INDEX],0,7 )  
                                                                      



# For no selection condition     
# Exp1 = 2017_02_22 -W condition
# Exp2 = 2017_03_03 -W condition
# Exp3 = 2017_08_22 -W condition
# Exp4 = 2017_06_08 -W / Seaprep condition  
P170_3_data = data.frame( Exp1=as.vector(log2(bP170_3_SW[P166,P166]+1)),Exp2=as.vector(log2(bP170_3_S8A[P166,P166]+1)),Exp3=as.vector(log2(bP170_3_S4A[P166,P166]+1)),Exp4=as.vector(log2(bP170_3_SQ[P166,P166]+1)))        
pdf("P170_3_data.@S1.V1.pdf",width = 600, height = 600); chart.Correlation(P170_4_data,method="pearson"); dev.off()

P170_4_data = data.frame( Exp1=as.vector(log2(bP170_4_SW[P166,P166]+1)),Exp2=as.vector(log2(bP170_4_S8A[P166,P166]+1)),Exp3=as.vector(log2(bP170_4_S4A[P166,P166]+1)),Exp4=as.vector(log2(bP170_4_SQ[P166,P166]+1)))        
pdf("P170_4_data.@S1.V1.pdf",width = 600, height = 600); chart.Correlation(P170_4_data,method="pearson"); dev.off()

P170_5_data = data.frame( Exp1=as.vector(log2(bP170_5_SW+1)),Exp2=as.vector(log2(bP170_5_SA8+1)),Exp4=as.vector(log2(bP170_4_SQ+1)))        
pdf("P170_5_data.@S1.V1.pdf",width = 600, height = 600); chart.Correlation(P170_5_data,method="pearson"); dev.off()                         



# Normalization methods - For Fig2
if (TRUE){
    eg_bm1_W = bm4_SW[5:18,5:18]    
    eg_bm1_A = bm4_SQ[5:18,5:18] 
	nm2 = eg_bm1_A
	rcnt=14
	colSumsMat2 = ( (1:rcnt * 0 + 1) %*% t(colSums(nm2)) )
	colSumsMat3 = nm2 /  colSumsMat2 
	simple_plotLog10Matrix(eg_bm1_A,1) 
	simple_plotMatrix(eg_bm1_A,0)   
	simple_plotMatrix(colSumsMat3,0)        
	
}

    
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
	     
	pdf("192spot.is2_bm2_A.auc.pdf")          
	par(mfrow=c(1,1)) 
	labels = pos_weak_index[is2_order]
	predictions = is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] 
	roc_obj = roc( labels, predictions  )
	if( TRUE ){plot(roc_obj) }   
	dev.off()   
	predictions = is2_bm1_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )    # 0.782
	predictions = is2_bm1_RA[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )   # 0.8239
	predictions = is2_bm2_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )    # 0.8077
	predictions = is2_bm2_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )    # 0.7959
	predictions = is2_bm3_A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )    # 0.8571
	predictions = is2_bm3_Q[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )    # 0.8035
	predictions = is2_bm3_QL[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )   # 0.7787
	predictions = is2_bm4_S4A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )  # 0.8472
	predictions = is2_bm4_S8A[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )  # 0.7682
	predictions = is2_bm4_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )   # 0.7765
	predictions = is2_bm6_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )   # 0.7993
	predictions = is2_bm7_SQ[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ; roc( labels, predictions  )   # 0.8061
	
	boxplot( c(0.782,0.8239,0.8077,0.7959,0.8571,0.8035,0.7787,0.8472,0.7682,0.7765,0.7993,0.8061) )  
	pdf("192spot.is2.auc.box.pdf")      
	boxplot( c(0.782,0.8239,0.8077,0.7959,0.8571,0.8035,0.7787,0.8472,0.7682,0.7765,0.7993,0.8061),ylim=c(0.7,0.9) )
	dev.off()
	
	## Fig 3c
	pdf("192spot.is2_overlap2.pdf")
	v = 3
	is2_overlap2 = (is2_bm1_A>v)+(is2_bm1_RA>v)+(is2_bm2_A>v)+(is2_bm2_Q>v)+(is2_bm3_A>v)+(is2_bm3_Q>v)+(is2_bm3_QL>v)+(is2_bm4_S4A>v)+(is2_bm4_S8A>v)+(is2_bm4_SQ>v)+(is2_bm6_SQ>v)+(is2_bm7_SQ>v)   


	par(mfrow=c(4,1)) 
	is2_order = order( is2_overlap2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ], decreasing = TRUE )
	barplot( is2_overlap2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
	barplot( pos_weak_index[is2_order] )        
	barplot( HIPPIE[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] ) 
	barplot( BioGrid[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] )
	dev.off()
	
	pdf("192spot.is2_overlap2.auc.pdf")          
	par(mfrow=c(1,1)) 
	labels = pos_weak_index[is2_order]
	predictions = is2_overlap2[ ((pos_weak_ALL[,2]-1)*78+pos_weak_ALL[,1])[,1] ][is2_order] 
	roc_obj = roc( labels, predictions  )
	if( TRUE ){plot(roc_obj) }   
	dev.off()
	
	
		
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
 


# rRNA experiments correlation
if (TRUE){
	brRNA_X_WH = data.frame( Exp1 = as.vector(brRNA1_WH), Exp2 = as.vector(brRNA2_WH), Exp3 = as.vector(brRNA3_WH) )
	pdf("brRNA_X_WH.log2.pdf")
	chart.Correlation(log2(brRNA_X_WH+1))           
	dev.off()
	pdf("brRNA_X_WH.pdf")
	chart.Correlation((brRNA_X_WH))           
	dev.off()
	
	brRNA_X_W = data.frame( Exp1 = as.vector(brRNA1_W), Exp2 = as.vector(brRNA2_W), Exp3 = as.vector(brRNA3_W) )    
	pdf("brRNA_X_W.log2.pdf")
	chart.Correlation(log2(brRNA_X_W+1))
	dev.off()
	pdf("brRNA_X_W.pdf")
	chart.Correlation((brRNA_X_W))
	dev.off()
	
}



# KRP experiments correlation
if (TRUE){
	KRP_X_WH_rpi = data.frame( Exp1 = as.vector(KRP_1_WH_rpi), Exp2 = as.vector(KRP_2_WH_rpi), Exp3 = as.vector(KRP_3_WH_rpi), Exp4 = as.vector(KRP_4_WH_rpi) )
	pdf("KRP_X_WH_rpi.log2.pdf")
	chart.Correlation(log2(KRP_X_WH_rpi+1))
	dev.off()
	
	KRP_X_W_rpi = data.frame( Exp1 = as.vector(KRP_1_W_rpi), Exp2 = as.vector(KRP_2_W_rpi), Exp3 = as.vector(KRP_3_W_rpi), Exp4 = as.vector(KRP_4_W_rpi) )
	pdf("KRP_X_W_rpi.log2.pdf")
	chart.Correlation(log2(KRP_X_W_rpi+1))
	dev.off()    
	
	KRP_X_ALL_rpi = data.frame( Exp1 = as.vector(KRP_1_W_rpi), Exp2 = as.vector(KRP_2_W_rpi), Exp3 = as.vector(KRP_3_W_rpi), Exp4 = as.vector(KRP_4_W_rpi), Exp1H = as.vector(KRP_1_WH_rpi), Exp2H = as.vector(KRP_2_WH_rpi), Exp3H = as.vector(KRP_3_WH_rpi), Exp4H = as.vector(KRP_4_WH_rpi) )         
	chart.Correlation(log2(KRP_X_ALL_rpi+1))        
	chart.Correlation((KRP_X_ALL_rpi))    
	
	KRP_X_ALL_rpi = data.frame( Exp1 = as.vector(KRP_1_W_rpi), Exp2 = as.vector(KRP_2_W_rpi), Exp3 = as.vector(KRP_3_W_rpi), Exp4 = as.vector(KRP_4_W_rpi), Exp1H = as.vector(KRP_1_WH_rpi), Exp2H = as.vector(KRP_2_WH_rpi), Exp3H = as.vector(KRP_3_WH_rpi), Exp4H = as.vector(KRP_4_WH_rpi) )         
	chart.Correlation(log2(KRP_X_ALL_rpi+1))        
	chart.Correlation((KRP_X_ALL_rpi))    
	
}                                 




#========================================
# Correlation with disorder vs interaction count
#========================================       
DrawCorrelationWithInteractionCount_P170 <-function(m1,m2){
	##
	# m1 = sw
	# m2 = sq
	library(readxl)
	Final_Tables <- read_excel("Desktop/Final_Tables.xlsx", sheet = "Disopred_P170")
	View(Final_Tables)
    
	#56, 57, 160
    temp_index = setdiff( 1:170, c(1:7,56,57,160) ) 
	#read_ratio =  (rowSums(bP170_5_SQ)/sum(bP170_5_SQ)) /(rowSums(bP170_4_SW)/sum(bP170_4_SW))     
	read_ratio =  (rowSums(m2)/sum(m2)) /(rowSums(m1)/sum(m1))     
	
	#read_ratio[ is.na( read_ratio ) ]  = 0       
	
	plot( read_ratio[temp_index], Final_Tables$`DISORDER %`[temp_index] )
    plot( read_ratio[temp_index], Final_Tables$`DISORDER Residues`[temp_index] )
    plot( read_ratio[temp_index], Final_Tables$`Total Residues`[temp_index] )

	cor( read_ratio[temp_index], Final_Tables$`DISORDER %`[temp_index] )
    cor( read_ratio[temp_index], Final_Tables$`DISORDER Residues`[temp_index] )
	cor( read_ratio[temp_index], Final_Tables$`Total Residues`[temp_index] )                     

	print( cor.test( read_ratio[temp_index], Final_Tables$`DISORDER %`[temp_index], na.rm = TRUE ) )
    print( cor.test( read_ratio[temp_index], Final_Tables$`DISORDER Residues`[temp_index], na.rm = TRUE ) )
    print( cor.test( read_ratio[temp_index], Final_Tables$`Total Residues`[temp_index], na.rm = TRUE ) )                     
} 
                                        
if (FALSE){                                                               
	DrawCorrelationWithInteractionCount_P170( bP170_1_W, bP170_1_Q )    	
	DrawCorrelationWithInteractionCount_P170( bP170_2_W, bP170_2_Q )
	DrawCorrelationWithInteractionCount_P170( bP170_3_W, bP170_3_Q )  
   	DrawCorrelationWithInteractionCount_P170( bP170_4_SW, bP170_4_SQ )               
   	DrawCorrelationWithInteractionCount_P170( bP170_5_SW, bP170_5_SQ )               
}



DetectAutoActivations <- function(is_m,cex=1.0){       
	m3 = t(is_m) - is_m
	boxplot(m3,cex=cex)        
	apply(m3,1,mean)
}

AutoActivator <- function( x ){     
	#x= ((rowSums(m2))/sum(m2))/((rowSums(m1+1))/sum(m1+1))
	
	#summarizing and plotting y     
	
	summary(x)
	hist(x, breaks = 20, col = rgb(0,0,1,0.5))
	boxplot(x, col = rgb(0,0,1,0.5), main = "Boxplot of x")
	shapiro.test(x)
	par(mfrow=c(1,1))
	qqnorm(x, main = "Normal QQ Plot - x")
	qqline(x, col = "red")           
	
	#dixon.test(x,opposite = TRUE)          
	#chisq.out.test(x,variance = var(x),opposite = TRUE)     
	
	           
	y = x[is.na(x)==FALSE]  
	print(y)       
	plot(y,main="y")
	output = Zeros(dim(m1)[1])
	index = 1
	for ( yy in sort(y)){yyy=chisq.out.test(y[y<=yy],variance = var(y[y<=yy]));print(yyy$p.value*length(y[y<=yy]));output[index]=yyy$p.value*length(y[y<=yy]);index = index+1}
	plot( -log10(output ) );	 
	#simple_plotLog2Matrix(m1,1,7)
	#simple_plotLog2Matrix(m2,1,7)
}

AutoActivationCorrection <- function(is_m){
	t_is_m = t(is_m)    
	temp1 = is_m - apply(is_m,1,median)        
	sigma1 = (rowSums( temp1>0 )+1) / (dim(temp1)[2]+1)
	temp2 = t_is_m - apply(t_is_m,1,median)        
	sigma2 = (rowSums( temp2>0 )+1) / (dim(temp2)[2]+1)
	
	#final = ( is_m - apply(is_m,1,median) ) / sqrt(sigma)        
	final = ( is_m - (2* apply(is_m,1,median) ) ) 
	final[which(final<0)] = 0
	return( final )
}   

 
plot(dnbinom(y,size=k1,mu=mu1))  

((prob*	dnbinom(y,size=k1,mu=mu1))+((1-prob)*dnbinom(y,size=k2,mu=mu2)))


bayes[is.na(bayes)] = 10e-20




# https://stats.stackexchange.com/questions/244529/three-components-negative-binomial-mixture-with-r
mixnbinom=function(y,k1,mu1,k2,mu2,prob,eps=1/100000)
{ 
	new.parms=c(k1,mu1,k2,mu2,prob) 
	err=1 
	iter=1 
	maxiter=100 
	hist(y,probability=T,nclass=30,col="lightgrey",main="The EM algorithm") 
	xvals=seq(min(y),max(y),1) 
	lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+(1-prob)*dnbinom(xvals,size=k2,mu=mu2),col="green") 
	while(err>eps){ 
		if(iter<=maxiter){ 
			lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+ 
				(1-prob)*dnbinom(xvals,size=k2,mu=mu2),lty=3) 
		} 
		bayes=(prob*dnbinom(y,size=k1,mu=mu1))/((prob*dnbinom(y,size=k1,mu=mu1))+((1-prob)*dnbinom(y,size=k2,mu=mu2)))        
		bayes[is.na(bayes)] = 10e-100    
		mu1=sum(bayes*y)/sum(bayes) 
		mu2=sum((1-bayes)*y)/sum((1-bayes)) 
		var1=sum(bayes*(y-mu1)^2)/sum(bayes) 
		var2=sum((1-bayes)*(y-mu2)^2)/sum((1-bayes)) 
		k1=abs(mu1/((var1/mu1)-1)) 
		k2=abs(mu2/((var2/mu2)-1)) 
		prob=mean(bayes) 
		old.parms=new.parms 
		new.parms=c(k1,mu1,k2,mu2,prob) 
		err=max(abs((old.parms-new.parms)/new.parms)) 
		iter=iter+1 
	} 
	lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+(1-prob)*dnbinom(xvals,size=k2,mu=mu2),col="red") 
	print(list(k1=k1,mu1=mu1,k2=k2,mu2=mu2,p=prob,iter=iter,err=err))   
	return(list(k1=k1,mu1=mu1,k2=k2,mu2=mu2,p=prob,iter=iter,err=err))
}
k1 = 1
mu1 = 1
k2 = 10
mu2 = 10
prob = 0.9
y = bm2_Q[5:18,5:18]
output = mixnbinom(y,k1,mu1,k2,mu2,prob)
k1 = output$k1
mu1 = output$mu1
k2 = output$k2
mu2 = output$mu2
prob = output$p
hist(y,probability=T,nclass=30,col="lightgrey",main="The EM algorithm",xlim=c(0,1000), breaks=1000) 
xvals=seq(0,1000,1) 
lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1),col="green") 
lines(xvals,(1-prob)*dnbinom(xvals,size=k2,mu=mu2),col="red") 
lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+(1-prob)*dnbinom(xvals,size=k2,mu=mu2),col="black") 


# https://stats.stackexchange.com/questions/244529/three-components-negative-binomial-mixture-with-r
mixnbinom_new=function(y,k1,mu1,mu2,sd2,prob,eps=1/100000)
{ 
	new.parms=c(k1,mu1,k2,mu2,prob) 
	err=1 
	iter=1 
	maxiter=100 
	hist(y,probability=T,nclass=30,col="lightgrey",main="The EM algorithm") 
	xvals=seq(min(y),max(y),1) 
	lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+(1-prob)*dnorm(xvals,mean=mu2,sd=sd2),col="green") 
	while(err>eps){ 
		if(iter<=maxiter){ 
			lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+ 
				(1-prob)*dnbinom(xvals,size=k2,mu=mu2),lty=3) 
		} 
		bayes=(prob*dnbinom(y,size=k1,mu=mu1))/((prob*dnbinom(y,size=k1,mu=mu1))+((1-prob)*dnorm(y,mean=mu2,sd=sd2)))        
		bayes[is.na(bayes)] = 10e-100    
		mu1=sum(bayes*y)/sum(bayes) 
		mu2=sum((1-bayes)*y)/sum((1-bayes)) 
		var1=sum(bayes*(y-mu1)^2)/sum(bayes) 
		sd2=sum((1-bayes)*(y-mu2)^2)/sum((1-bayes)) 
		k1=abs(mu1/((var1/mu1)-1)) 
		k2=abs(mu2/((var2/mu2)-1)) 
		prob=mean(bayes) 
		old.parms=new.parms 
		new.parms=c(k1,mu1,k2,mu2,prob) 
		err=max(abs((old.parms-new.parms)/new.parms)) 
		iter=iter+1 
	} 
	lines(xvals,prob*dnbinom(xvals,size=k1,mu=mu1)+(1-prob)*dnbinom(xvals,size=k2,mu=mu2),col="red") 
	print(list(k1=k1,mu1=mu1,k2=k2,mu2=mu2,p=prob,iter=iter,err=err))   
	return(list(k1=k1,mu1=mu1,k2=k2,mu2=mu2,p=prob,iter=iter,err=err))
}
	
