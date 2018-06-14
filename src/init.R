#===================================================================================================================================== 
#
# Version 1.0
# init.R
#
#
#=====================================================================================================================================  

#=====================================================================================================================================  
## Lord libraries    

library("RColorBrewer")
require(pheatmap)
require(gtools)
#require(gplots)
require(reshape2)
#library(psych)
library(clipr)
#library(swfscMisc)               # for uniform.test    
library(PerformanceAnalytics)
library(mixtools)
library(pROC)    
#library(d3heatmap)
library(outliers)
library(readxl)  
library(mixtools) 
library(matrixStats)

#=====================================================================================================================================  
# Reference Informations                                                                                                 
#                 
# matrix[x,y]
# x --> row (from 1..)
# y --> col (from 1..)
#

#read.clipboard(header = FALSE) ==> copy from clipboard
                                     
# Ray Information
# https://pymolwiki.org/index.php/Ray                                                       


#=====================================================================================================================================  
## Key Parameters
                  
LOG_SCALE = 10

#=====================================================================================================================================  
## Functions
  
Zeros <- function( n ){
	return ( 1:n * 0 )
}                     
                            
NormalizeMatrix <- function(m){
	return ( RefineMatrix(m)/sum(RefineMatrix(m)) )
}

NonZeroMinimum <- function( m ){
   return( min( m[m!=0] ) )
}   
            
Ones <- function( len ){
	output = 1:len / 1:len
	return( output )
}                   

Zeros <- function( len ){
	output = 1:len * 0
	return (output )
}    
            
### 
# Confidence Interval around a Linear Regression Line
# https://rpubs.com/aaronsc32/regression-confidence-prediction-intervals
reg.conf.intervals <- function(x, y) {
  n <- length(y) # Find length of y to use as sample size
  lm.model <- lm(y ~ x) # Fit linear model
  
  # Extract fitted coefficients from model object
  b0 <- lm.model$coefficients[1]
  b1 <- lm.model$coefficients[2]
  
  # Find SSE and MSE
  sse <- sum((y - lm.model$fitted.values)^2)
  mse <- sse / (n - 2)
  
  t.val <- qt(0.975, n - 2) # Calculate critical t-value
  
  # Fit linear model with extracted coefficients
  x_new <- 1:max(x)
  y.fit <- b1 * x_new + b0
  
  # Find the standard error of the regression line
  se <- sqrt(sum((y - y.fit)^2) / (n - 2)) * sqrt(1 / n + (x - mean(x))^2 / sum((x - mean(x))^2))
  
  # Fit a new linear model that extends past the given data points (for plotting)
  x_new2 <- 1:max(x + 100)
  y.fit2 <- b1 * x_new2 + b0
  
  # Warnings of mismatched lengths are suppressed
  slope.upper <- suppressWarnings(y.fit2 + t.val * se)
  slope.lower <- suppressWarnings(y.fit2 - t.val * se)
  
  # Collect the computed confidence bands into a data.frame and name the colums
  bands <- data.frame(cbind(slope.lower, slope.upper))
  colnames(bands) <- c('Lower Confidence Band', 'Upper Confidence Band')
  
  # Plot the fitted linear regression line and the computed confidence bands
  plot(x, y, cex = 1.75, pch = 21, bg = 'gray')
  lines(y.fit2, col = 'black', lwd = 2)
  lines(bands[1], col = 'blue', lty = 2, lwd = 2)
  lines(bands[2], col = 'blue', lty = 2, lwd = 2)
  
  return(bands)
}
            

### Need to improve!!                
# 
# http://www.dummies.com/programming/big-data/data-science/hypothesis-tests-for-data-outliers/
# 
AutoActivator <- function( m1, m2 ){     
	x= ((rowSums(m2))/sum(m2))/((rowSums(m1+1))/sum(m1+1))
	
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
     
background_removeMatrix <- function( data, beta, alpha, fontsize, colorname = "Blues"){
	mixture = normalmixEM(as.vector(data))          
	mu = min(mixture$mu)
	sigma = min(mixture$sigma)  
	print(mixture) 
	print((mu+(sigma*2))*beta)
	simple_plotLog2Matrix( data-((mu+(sigma*2))*beta), alpha, fontsize, colorname=colorname )
}   

plotMoG2 <- function( data, min_v, max_v ){                    
	#mix_data = normalmixEM(data,k=3,epsilon=1e-3,fast=TRUE) 
	#mix_data0 = normalmixEM(data,k=2) 
	
	#data2 = data - min( mix_data0$mu )
	data2 = data
	     
	hdata = hist(data2,100,freq=FALSE, xlim=c(min_v,max_v))
	
	mix_data = normalmixEM(data2,k=2)        
	
	#x1 <- seq(min_v,max_v,length=100)*mix_data$sigma[1] + mix_data$mu[1]
	x1 <- seq(min_v,max_v,length=100) # + mix_data$mu[1]
	hx1 <- dnorm(x1,mix_data$mu[1],mix_data$sigma[1]) * mix_data$lambda[1]  
 
	#x2 <- seq(min_v,max_v,length=100)*mix_data$sigma[2] + mix_data$mu[2]
	x2 <- seq(min_v,max_v,length=100) # + mix_data$mu[2]
	hx2 <- dnorm(x2,mix_data$mu[2],mix_data$sigma[2]) * mix_data$lambda[2]   
	
	lines(x1, hx1, col="red", lwd=2)	
	lines(x2, hx2, col="Blue", lwd=2)  
	
	return (mix_data)  
}                 
  
plotMoG3 <- function( data, min_v, max_v ){                    
	#mix_data = normalmixEM(data,k=3,epsilon=1e-3,fast=TRUE) 
	#mix_data0 = normalmixEM(data,k=3) 
	
	#data2 = data - min( mix_data0$mu )
	data2 = data
	     
	hdata = hist(data2,100,freq=FALSE, xlim=c(min_v,max_v))
	
	mix_data = normalmixEM(data2,k=3)        
	
	#x1 <- seq(min_v,max_v,length=100)*mix_data$sigma[1] + mix_data$mu[1]
	x1 <- seq(min_v,max_v,length=100) # + mix_data$mu[1]
	hx1 <- dnorm(x1,mix_data$mu[1],mix_data$sigma[1]) * mix_data$lambda[1]  
 
	#x2 <- seq(min_v,max_v,length=100)*mix_data$sigma[2] + mix_data$mu[2]
	x2 <- seq(min_v,max_v,length=100) # + mix_data$mu[2]
	hx2 <- dnorm(x2,mix_data$mu[2],mix_data$sigma[2]) * mix_data$lambda[2]   
	
	hx3 <- dnorm(x2,mix_data$mu[3],mix_data$sigma[3]) * mix_data$lambda[3]   
	
	lines(x1, hx1, col="red", lwd=2)	
	lines(x2, hx2, col="Green", lwd=2) 
	lines(x2, hx3, col="Blue", lwd=2)    
	
	return (mix_data)
}

BackgroundSignal <-function( x, min=0, max=5 ){    
	if (FALSE){  
		y = c( as.vector( log2(x[x!=0]) ), as.vector( log2(x[x!=0]) * (-1) ) )     
		#y = c( as.vector( log2(x+1) ), as.vector( log2(x+1) * (-1) ) )     
		hist(y) 
		plotMoG3(y,-15,15)
		return( y )
	}   
	if (TRUE){      
		y = as.vector( log2(x+1) )
		hist(y) 
		plotMoG2(y,min=min,max=max)
		return( y )
	}
}

ToxicityTest <-function(m, column_start, column_end, threshold = 0){
	x1 = c() # to see how far from uniform-distribution (we noticed that toxic column often not uniform)
	x2 = c() # to check the sum of read
	x3 = c() # to check how many 0 reads
	x4 = c() # to check how many rows are under threshold         
	nrows = dim(m)[1]
	ncols = dim(m)[2]
    for (i in column_start:column_end)
    {                                                                              
	    #print(i)                                                                   
	    tryCatch( {y=uniform.test( hist(m[,i]) ); x1[i] = y$statistic}, error=function(error_message){return(NULL)} )
	    x2[i] = sum(m[,i])     
		x3[i] = sum(m[,i]==0) / nrows 
		x4[i] = sum(m[,i]<=threshold) / nrows 
    }                
	names(x1) = colnames(m[,column_start:column_end])  
	names(x2) = colnames(m[,column_start:column_end])  
	names(x3) = colnames(m[,column_start:column_end])  
	names(x4) = colnames(m[,column_start:column_end])  
    x1[which(is.na(x))] = -1     
	print("x1 - Uniform-test statistics")
	print(x1)
	print("x2 - PPI read count")
	print(x2)
	print("x3 - Percentage of zeros") 
	print(x3)
	print("x4 - Percentage of rows under threshold") 
	print(x4)
    plot(x1,x2,xlab="Uniform-test statistics",ylab="PPI read count")     
	plot(x3,x2,xlab="Percentage of zeros",ylab="PPI read count")  
	plot(x4,x2,xlab="Percentage of rows under threshold",ylab="PPI read count")    
	plot(x3,x1,xlab="Percentage of zeros",ylab="Uniform-test statistics")  
	
	return (data.frame(read_count=x2,uniform_test=x1,zero_percentage=x3,threshold_percentage=x4))
}

makeNullMatrix <- function( ppi_output_path, alpha ){
	NullMatrix <- read.delim( ppi_output_path )   
	row.names(NullMatrix) <- NullMatrix$DB.Read.1....AD.Read.2.  
	rcnt = dim(NullMatrix)[2]
	NullMatrix <- NullMatrix[,2:rcnt]      
	NullMatrix_matrix <- data.matrix(NullMatrix) + alpha
	total_reads = sum(NullMatrix_matrix)  
	freq_matrix = NullMatrix_matrix / total_reads   
	rows = rowSums(freq_matrix)
	cols = colSums(freq_matrix)
	reconst_matrix = rows%*%t(cols)
	pheatmap(log10(freq_matrix), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))   
	pheatmap(log10(reconst_matrix), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))   
	return (reconst_matrix)
}      

makeFreqMatrix <- function( ppi_output_path, alpha ){
	FreqMatrix <- read.delim( ppi_output_path )   
	row.names(FreqMatrix) <- FreqMatrix$DB.Read.1....AD.Read.2.  
	rcnt = dim(FreqMatrix)[2]
	FreqMatrix <- FreqMatrix[,2:rcnt]      
	FreqMatrix_matrix <- data.matrix(FreqMatrix)  
	total_reads = sum(FreqMatrix_matrix)  
	freq_matrix = FreqMatrix_matrix / total_reads       
	log2_freq_matrix = log2(freq_matrix)
	log2_freq_matrix[which(log2_freq_matrix==-Inf)] = NA 
	log2_freq_matrix[which(log2_freq_matrix==Inf)] = NA  
	#pheatmap(log10(freq_matrix), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues")) 
	#heatmap.2(log2_freq_matrix, Rowv = F, Colv = F, trace = "none", na.color = "Black",col = colorRampPalette(c("white","blue"))(50))     
	heatmap.2(log2_freq_matrix, Rowv = F, Colv = F, trace = "none", na.color = "Black",col = brewer.pal(9,"Blues") ) 
	#heatmap.2(log2_freq_matrix, Rowv = F, Colv = F, trace = "none", na.color = "Black",col = colorRampPalette(c("white","#08306B"))(50))   
	return (freq_matrix)
}           

getClipboardMatrix <- function(header=TRUE){
	x = read.clipboard(header=header)
	View(x)
	row.names(x) = x[,1]
	rcnt = dim(x)
	rcnt = dim(x)[2]
	x = x[,2:rcnt]
	heatmap.2(as.matrix(x), Rowv = F, Colv = F, trace = "none", na.color = "Black",col = brewer.pal(9,"Blues") ) 
	return (x)
}    
                                                      
write.clipboard <- function(x, sep="\t", row.names=FALSE){
	write_clip(x) 
}
  
#=====================================================
# convertNullMatrix
#
#=====================================================   
convertNullMatrix <- function( matrix, alpha, filename1=NA, filename2=NA, cellwidth = NA, cellheight = NA ){      
	NullMatrix_matrix <- data.matrix(matrix) + alpha
	total_reads = sum(NullMatrix_matrix)  
	freq_matrix = NullMatrix_matrix / total_reads   
	rows = rowSums(freq_matrix)
	cols = colSums(freq_matrix)
	reconst_matrix = rows%*%t(cols)
	pheatmap(log10(freq_matrix), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), filename=filename1, cellwidth = cellwidth, cellheight = cellwidth)   
	pheatmap(log10(reconst_matrix), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), filename=filename2, cellwidth = cellwidth, cellheight = cellwidth)   
	return (reconst_matrix) 
}

#=====================================================
# convertFreqMatrix
#
#=====================================================   
convertFreqMatrix <- function( matrix, alpha, filename=NA, cellwidth = NA, cellheight = NA ){   
	FreqMatrix_matrix <- data.matrix(matrix) + alpha   
	total_reads = sum(FreqMatrix_matrix)  
	freq_matrix = FreqMatrix_matrix / total_reads       
	log10_freq_matrix = log10(freq_matrix)
	log10_freq_matrix[which(log10_freq_matrix==-Inf)] = NA 
	log10_freq_matrix[which(log10_freq_matrix==Inf)] = NA  
	heatmap.2(log10_freq_matrix, Rowv = F, Colv = F, trace = "none", na.color = "Black",col = brewer.pal(9,"Blues"), )  
	pheatmap(log10_freq_matrix, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), filename=filename, cellwidth = cellwidth, cellheight = cellwidth)  
	return (freq_matrix)	
}   
        
AutoActivationCorrection <- function(is_m){    
	temp = is_m - apply(is_m,1,median)        
	simple_plotMatrix( temp, 0 )    
	sigma = (rowSums( temp>0 )+1) / (dim(temp)[2]+1)
	final = ( is_m - apply(is_m,1,median) ) / sigma        
	#final = ( is_m - (2* apply(is_m,1,median) ) ) 
	final[which(final<0)] = 0
	return( final )
}   

AutoActivationCorrection2 <- function(is_m){    
	final = is_m - apply(is_m,1,median)       
	final[which(final<0)] = 0 
	return (final)
}  

AutoActivationCorrection3 <- function(is_m){    
	final = is_m - apply(is_m,1,quantile,3/4)      
	final[which(final<0)] = 0 
	return (final)
	
	#median_values = apply(is_m,1,median)       
	#final = is_m - median_values
	#final[which(final<0)] = 0 
	#for( i in 1:length(median_values) ){
	#	if ( median_values[i] > 0 ){
	#		final[i,] = final[i,] / 2 
	#	}
	#}
	#return (final)
}  


RawInteractionScores <- function(m1,m2, alpha,show=TRUE, filename=NA, cellwidth = NA, cellheight = NA){
	# m1 : no selection matrix
	# m2 : selection matrix  
	                 
	NullMatrix_matrix1 <- data.matrix(m1) # + alpha  
	total_reads1 = sum(NullMatrix_matrix1)  
	freq_matrix1 = NullMatrix_matrix1 / total_reads1   
	rows1 = rowSums(freq_matrix1)
	cols1 = colSums(freq_matrix1)
	reconst_matrix1 = rows1%*%t(cols1)  # null
	
	total_reads2 = sum(m2)
	freq_matrix2 = m2 / total_reads2   
	                  
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

NewInteractionScores <- function(m1,m2, alpha,fontsize=9, sd_type=1, filename=NA, cellwidth = NA, cellheight = NA){  
	# this considers standard deviation
	# m1 : no selection matrix
	# m2 : selection matrix  
	                 
	NullMatrix_matrix1 <- data.matrix(m1) + alpha  
	total_reads1 = sum(NullMatrix_matrix1)  
	freq_matrix1 = NullMatrix_matrix1 / total_reads1   
	#rows1 = rowSums(freq_matrix1) 
	#row.sds <- apply(freq_matrix1, 1, sd)  
	#cols1 = colSums(freq_matrix1)
	#col.sds <- apply(freq_matrix1, 2, sd)
	#reconst_matrix1 = rows1%*%t(cols1)  # null
	#reconst_sd_null_matrix1 = row.sds%*%t(row.sds)  # null   
	#reconst_sd_null_matrix2 = col.sds%*%t(col.sds)  # null                            
	#reconst_sd_null_matrix3 = sqrt( reconst_sd_null_matrix1 + reconst_sd_null_matrix2 )
	                                
	# proper sd calculation : http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm
	rows1 = rowSums(freq_matrix1) 
	row.sds <- apply(freq_matrix1, 1, sd)  
	row.sds <- row.sds ^ 2
	cols1 = colSums(freq_matrix1)
	col.sds <- apply(freq_matrix1, 2, sd)
	col.sds <- col.sds ^ 2
	reconst_matrix1 = rows1%*%t(cols1)  # null
	reconst_sd_null_matrix1 = row.sds%*%t(rep(1,length(cols1)))  # null   
	reconst_sd_null_matrix2 = rep(1,length(rows1))%*%t(col.sds)  # null                            
	reconst_sd_null_matrix3 = sqrt( reconst_sd_null_matrix1 + reconst_sd_null_matrix2 )
	 
	
	total_reads2 = sum(m2)
	freq_matrix2 = m2 / total_reads2   
	#row.sds2 <- apply(freq_matrix2, 1, sd)  
	#col.sds2 <- apply(freq_matrix2, 2, sd)
	#reconst_sd_matrix1 = row.sds2%*%t(row.sds2)  # null   
	#reconst_sd_matrix2 = col.sds2%*%t(col.sds2)  # null                            
	#reconst_sd_matrix4 = sqrt( reconst_sd_matrix1 + reconst_sd_matrix2 )
	                
	if (sd_type==1){
		reconst_sd_matrix_final = reconst_sd_null_matrix3;     
	}
	else{
		reconst_sd_matrix_final = sqrt( reconst_sd_null_matrix1 + reconst_sd_null_matrix2 + reconst_sd_matrix1 + reconst_sd_matrix2 ); 
	}   
	        
	is_matrix = ( freq_matrix2 - reconst_matrix1 ) / reconst_sd_matrix_final 
	
	#log2_is_matrix = log2( is_matrix )
	#log2_is_matrix[ which(log2_is_matrix==-Inf)] = NA
	#log2_is_matrix[ which(log2_is_matrix==Inf)] = NA   
	#log2_is_matrix[ which(log2_is_matrix=="NaN")] = NA
	
	#heatmap.2(log2_is_matrix, Rowv = F, Colv = F, trace = "none", na.color = "Black",col = colorRampPalette(c("grey","white","blue"))(21)) 
	          
	is_matrix[which(is_matrix=="NaN")] = 0
	is_matrix[which(is_matrix==Inf)] = 0  
	is_matrix[which(is_matrix<0)] = 0  
	
	pheatmap(log2(is_matrix+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth)  
	return(log2(is_matrix+alpha))
	
	#return( is_matrix )
}                                  


Epistasis <- function(m1, alpha,fontsize=9, sd_type=1, filename1=NA, filename2=NA,cellwidth = NA, cellheight = NA){  
	# this considers standard deviation
	# m1 : no selection matrix
	                 
	NullMatrix_matrix1 <- data.matrix(m1) #+ alpha  
	total_reads1 = sum(NullMatrix_matrix1)  
	freq_matrix1 = NullMatrix_matrix1 / total_reads1    # Observed Pair Frequency
                                    
	rows1 = rowSums(freq_matrix1) 
	row.sds <- apply(freq_matrix1, 1, sd)  
	row.sds <- row.sds ^ 2
	cols1 = colSums(freq_matrix1)
	col.sds <- apply(freq_matrix1, 2, sd)
	col.sds <- col.sds ^ 2
	reconst_matrix1 = rows1%*%t(cols1)  # null
	reconst_sd_null_matrix1 = row.sds%*%t(rep(1,length(cols1)))  # null   
	reconst_sd_null_matrix2 = rep(1,length(rows1))%*%t(col.sds)  # null                            
	reconst_sd_null_matrix3 = sqrt( reconst_sd_null_matrix1 + reconst_sd_null_matrix2 )
	
	plot( rows1, main="rows" )
	plot( cols1, main="cols" )
	reconst_matrix1 = rows1%*%t(cols1)  # Expected Pair Frequency
                      
	pheatmap(reconst_matrix1, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, fontsize=fontsize, filename=filename1, cellwidth = cellwidth, cellheight = cellwidth)
	
	dev_matrix1 = ( freq_matrix1 - reconst_matrix1 ) / reconst_sd_null_matrix3
	              
	pheatmap(dev_matrix1, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, fontsize=fontsize, filename=filename2, cellwidth = cellwidth, cellheight = cellwidth)  
	return(dev_matrix1)
}
               

                                       


InteractionScores <- function(m1,m2, alpha,fontsize=9, filename=NA,cellwidth = NA, cellheight = NA){              
	raw_is = RawInteractionScores(m1,m2,alpha,FALSE)
	raw_is[which(raw_is=="NaN")] = 0
	raw_is[which(raw_is==Inf)] = 0  
	raw_is[which(raw_is==-Inf)] = 0 
	pheatmap(log2(raw_is+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth)  
	return(log2(raw_is+alpha))
} 


PairInteractionScores <- function(m1,m2, alpha,fontsize=9, filename=NA,cellwidth = NA, cellheight = NA){       
	freq_matrix1 = m1 /sum(m1)    # Observed Pair Frequency
	freq_matrix2 = m2 /sum(m2)             
	raw_is = freq_matrix2 %/% freq_matrix1
	raw_is[which(raw_is=="NaN")] = 0
	raw_is[which(raw_is==Inf)] = 0  
	raw_is[which(raw_is==-Inf)] = 0                                   
	pheatmap(log2(raw_is+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth )  
	return(log2(raw_is+alpha))
} 


                                                                                                              
drawMatrix <- function( ppi_output_path, alpha,fontsize=9,interactive=FALSE,draw=TRUE){
	NullMatrix <- read.delim( ppi_output_path, comment.char = "#" )   
	row.names(NullMatrix) <- NullMatrix$DB.Read.1....AD.Read.2.  
	rcnt = dim(NullMatrix)[2]
	NullMatrix <- NullMatrix[,2:rcnt]      
	NullMatrix_matrix <- data.matrix(NullMatrix)  
	if (sum(abs(NullMatrix_matrix)) != 0){                                 
		if (interactive == TRUE){  
			#d3heatmap(log10(NullMatrix_matrix+alpha))           
			d3heatmap(log10(NullMatrix_matrix+alpha),Rowv=NULL,Colv=NULL) 
		}else{
			if (draw == TRUE ){
				pheatmap(log10(NullMatrix_matrix+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize)  
			}
		}   
	}
	return (NullMatrix_matrix)
}  
     
if (FALSE){
	ppi_output_path = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt.new.2"    
	ppi_output_path = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S1_W.ppi.txt"         
	ppi_output_path = "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S2_Q.ppi.txt"          
	alpha =0.5
	fontsize=9
	interactive = TRUE
	NullMatrix <- read.delim( ppi_output_path, comment.char = "#" )   
	row.names(NullMatrix) <- NullMatrix$DB.Read.1....AD.Read.2.  
	rcnt = dim(NullMatrix)[2]
	NullMatrix <- NullMatrix[,2:rcnt]      
	NullMatrix_matrix <- data.matrix(NullMatrix)        
	#d3heatmap(log10(NullMatrix_matrix+alpha))
	#help("d3heatmap")
	d3heatmap(log10(NullMatrix_matrix+alpha),Rowv=NULL,Colv=NULL)
	#d3heatmap(log10(NullMatrix_matrix+alpha),Rowv=NULL,Colv=NULL,fontsize=8)
}
 
freq_plotLog10Matrix <- function( m, alpha,fontsize=9, cluster_rows=FALSE,cluster_cols = FALSE, filename = NA, colorname = "Blues", cellwidth = NA, cellheight = NA , breaks = NA, step = 9 ){     
	m2 = (m + alpha) / sum( ( m + alpha ) )    
	#color = c( "#FFFFFF", brewer.pal(step-1,colorname) )         
	#color = brewer.pal(step,colorname)                             
	
	colfunc <- colorRampPalette(c("white", "darkblue"))
	color = colfunc(step)
	#heatmap.2((KRP_1_WH_rpi),trace="none",symkey=FALSE, density.info="none",cexRow=0.9,col=colfunc(30),Rowv=TRUE, Colv=TRUE,dendrogram='none')
	               
	print( max(log10(m2)))
	print( min(log10(m2)))
	
	pheatmap(log10(m2), cluster_rows=cluster_rows, cluster_cols = cluster_cols, color = color, fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth , breaks = breaks )   
}    



simple_plotLog10Matrix <- function( m, alpha,fontsize=9, cluster_rows=FALSE,cluster_cols = FALSE, filename = NA, colorname = "Blues", cellwidth = NA, cellheight = NA , breaks = NA ){    
	color = brewer.pal(9,colorname)                      
	pheatmap(log10(m+alpha), cluster_rows=cluster_rows, cluster_cols = cluster_cols, color = color, fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth , breaks = breaks )   
}
     
heatmap2_plotLog2Matrix <- function( m, alpha,fontsize=9, cluster_rows=FALSE,cluster_cols = FALSE, beta = 0, filename = NA, cellwidth = NA, cellheight = NA   ){          
	 m2 = m - beta
	 m2[m2<0] = 0         
	 if ( filename == "" ){              
	 	heatmap.2(log2(m2+alpha), sclae="none",Rowv = F, Colv = F, trace = "none", na.color = "Black", fontsize=fontsize, col = colorRampPalette(c("white","white","blue"))(21), cellwidth = cellwidth, cellheight = cellwidth )      
	 }else{   
		pdf( filename )
		heatmap.2(log2(m2+alpha), sclae="none",Rowv = F, Colv = F, trace = "none", na.color = "Black", fontsize=fontsize, col = colorRampPalette(c("white","white","blue"))(21), cellwidth = cellwidth, cellheight = cellwidth )         
		dev.off()
	}
}       


simple_plotLog2Matrix <- function( m, alpha,fontsize=9, cluster_rows=FALSE,cluster_cols = FALSE, beta = 0, filename = NA, colorname = "Blues", cellwidth = NA, cellheight = NA, breaks = NA   ){    
	 # colorname = "YlOrRd"      
	 m2 = m - beta
	 m2[m2<0] = 0            
	 pheatmap(log2(m2+alpha), cluster_rows=cluster_rows, cluster_cols = cluster_cols, color = brewer.pal(9,colorname), fontsize=fontsize, filename=filename , cellwidth = NA, cellheight = NA, breaks = breaks  )  
}

simple_plotMatrix <- function( m, alpha,fontsize=9, cluster_rows=FALSE,cluster_cols = FALSE, filename = NA, cellwidth = NA, cellheight = NA  ){    
	pheatmap((m+alpha), cluster_rows=cluster_rows, cluster_cols = cluster_cols, color = brewer.pal(9,"Blues"), fontsize=fontsize, filename = filename, cellwidth = cellwidth, cellheight = cellwidth   )  
}

plotLog10Matrix <- function( m, alpha, colorList, breaksList ){                       
	 #breaksList = seq(0, 5, by = 0.1)
	 #colorList = brewer.pal(9,"Blues")   
	 #colorList = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
	 pheatmap(log10(m+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = colorList, breaks = breaksList )   
}  

plotLog2Matrix <- function( m, alpha, color = NA, breaks = NA ){                       
	 #breaksList = seq(0, 12, by = 0.1)
	 #colorList = brewer.pal(9,"Blues")   
	 #colorList = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))
	 pheatmap(log2(m+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = color, breaks = breaks )   
}  

SetMax <- function(mat, limit){
	mat2 = mat #matrix(mat, nrow(mat), ncol(mat))
	#names(mat2) <- row.names(mat)
	#row.names(mat2) <- row.names(mat)
	mat2[which(mat2>limit)] = limit
	return(mat2)
} 

SetMin <- function(mat, limit){
	mat2 = matrix(mat)
	mat2[which(mat2<limit)] = limit
	return(mat2)
}
  
CorMatrix <-function(m1,m2,xlab="",ylab="",main="", method="pearson",lwd=0.5, cex=0.3, LogScale = FALSE, showABLine = TRUE, showPlot = TRUE ) {   
	if (showPlot == FALSE){
		if (LogScale == TRUE){
			return( cor( as.vector(log2(m1+1)), as.vector(log2(m2+1)), method = method ) )		
		}else{
			return( cor( as.vector(m1), as.vector(m2), method = method ) )
		}	
	}
	if (LogScale == TRUE){
		plot(log2(m1+1),log2(m2+1),xlab=xlab,ylab=ylab,main=main,lwd=lwd,cex=cex)  
		if ( showABLine ){
			abline( lm( as.vector(log2(m1+1)) ~ as.vector(log2(m2+1)) ) )
		}
		return( cor( as.vector(log2(m1+1)), as.vector(log2(m2+1)), method = method ) )		
	}else{
		plot(m1,m2,xlab=xlab,ylab=ylab,main=main,lwd=lwd,cex=cex)
		if ( showABLine ){
			abline( lm( as.vector(m1) ~ as.vector(m2) ) )
		}
		return( cor( as.vector(m1), as.vector(m2), method = method ) )
	}
}
                                    
# Added on 2018/03/08 (while revision)
CorMatrixWithoutZeros <-function(m1,m2,xlab="",ylab="",main="", method="pearson",lwd=0.5, cex=0.3, LogScale = FALSE, showABLine = TRUE, showPlot = TRUE, cutoff = 0.0 ) {    
	non_zero_index = intersect( which(m1>cutoff), which(m2>cutoff) )  
	m1 = m1[non_zero_index]
	m2 = m2[non_zero_index]
	if (showPlot == FALSE){
		if (LogScale == TRUE){
			return( cor( as.vector(log2(m1+1)), as.vector(log2(m2+1)), method = method ) )		
		}else{
			return( cor( as.vector(m1), as.vector(m2), method = method ) )
		}	
	}
	if (LogScale == TRUE){
		plot(log2(m1+1),log2(m2+1),xlab=xlab,ylab=ylab,main=main,lwd=lwd,cex=cex)  
		if ( showABLine ){
			abline( lm( as.vector(log2(m1+1)) ~ as.vector(log2(m2+1)) ) )
		}
		return( cor( as.vector(log2(m1+1)), as.vector(log2(m2+1)), method = method ) )		
	}else{   
		plot(m1,m2,xlab=xlab,ylab=ylab,main=main,lwd=lwd,cex=cex)
		if ( showABLine ){
			abline( lm( as.vector(m1) ~ as.vector(m2) ) )
		}
		return( cor( as.vector(m1), as.vector(m2), method = method ) )
	}
}   

SaveHeatMap <- function(mat, outpath) {
	pheatmap(mat, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
	png(filename = outpath, width = 960, height = 960 )
	pheatmap(mat, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
	dev.off()
}

DrawHeatMap <- function(mat) {
	pheatmap(mat, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
	pheatmap(log10(mat), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
}

ReadMatrix <- function(filepath) {
	mat <- read.delim(filepath, header = TRUE)
		View(mat)
		row.names(mat) <- mat[1:78,1]
		mat <- mat[,2:79]
		mat <- data.matrix(mat)
		View(mat)
		return(mat)
		pheatmap(mat, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
		pheatmap(log10(mat), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
		return(mat)
}  

ReadMatrixP170 <- function(filepath) {
	mat <- read.delim(filepath, header = TRUE)
		View(mat)
		row.names(mat) <- mat[1:173,1]
		mat <- mat[,2:174]
		mat <- data.matrix(mat)
		View(mat)
		return(mat)
		pheatmap(mat, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
		pheatmap(log10(mat), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"))
		return(mat)
}

BasicStats <- function(m){                
	BDs = dim(m)[1]      
	BDs_Zeros = length( which( rowSums(m) == 0 ) )
	ADs = dim(m)[2]      
	ADs_Zeros = length( which( colSums(m) == 0 ) )     
	All = BDs * ADs
	Zeros = length(which(m==0))
	print( sprintf("BD cnt = %d, %0.1f %% detected (%d)", BDs, (1.0-(BDs_Zeros/BDs))*100.0, BDs - BDs_Zeros ) )
	print( sprintf("AD cnt = %d, %0.1f %% detected (%d)", ADs, (1.0-(ADs_Zeros/ADs))*100.0, ADs - ADs_Zeros ) )
	print( sprintf("ALL cnt = %d", All) )
	print( sprintf("ZERO cnt = %d", Zeros) )   
	print( sprintf("Detected cnt = %d", All-Zeros ) )
	print( sprintf("%% detected cnt = %0.1f", (All-Zeros)/All*100.0) )
	print( sprintf("sum = %d", sum(m)) )
	print( sprintf("max read = %d", max(m)) ) 
	print( sprintf("min read = %d", min(m)) ) 
	print( sprintf("max AD = %d", max(colSums(m))) ) 
	print( sprintf("min AD = %d", min(colSums(m))) ) 
	print( sprintf("max BD = %d", max(rowSums(m))) ) 
	print( sprintf("max BD = %d", min(rowSums(m))) ) 
	
	par(mfrow=c(2,2))
	plot(sort(rowSums(m)), xlab="Genes (sorted by sum of reads)", ylab="Marginal sum of reads", main="Bindinding Domain")
	plot(sort(colSums(m)), xlab="Genes (sorted by sum of reads)", ylab="Marginal sum of reads", main="Activation Domain")  
	hist(log(m,2), xlab="log2 Pair-reads", ylab="Frequency", main="Pair-reads distribution")   
	barplot(c(All-Zeros,Zeros),xlab="Cases",ylab="Gene pairs (count)",names=c("detected","not detected"),main="Coverage")
}              

BasicStats2 <- function(m1, m2, INT_CUTOFF=0){
	print( "============== m1 ==============")
	BasicStats(m1)     
	print( "============== m2 ==============")      
	BasicStats(m1)           
	
	iCount1 = (m1>INT_CUTOFF)
	iCount2 = (m2>INT_CUTOFF)
	iCount = iCount1 + iCount1

	print( sum(iCount1) )
	print( sum(iCount2) )
	x = hist(iCount)  
	print( x )
	                
	#BDs = dim(m)[1]      
	#BDs_Zeros = length( which( rowSums(m) == 0 ) )
	#ADs = dim(m)[2]      
	#ADs_Zeros = length( which( colSums(m) == 0 ) )     
	#All = BDs * ADs
	#Zeros = length(which(m==0))
	#print( sprintf("BD cnt = %d, %0.1f %% detected (%d)", BDs, (1.0-(BDs_Zeros/BDs))*100.0, BDs - BDs_Zeros ) )
	#print( sprintf("AD cnt = %d, %0.1f %% detected (%d)", ADs, (1.0-(ADs_Zeros/ADs))*100.0, ADs - ADs_Zeros ) )
	#print( sprintf("ALL cnt = %d", All) )
	#print( sprintf("ZERO cnt = %d", Zeros) )   
	#print( sprintf("Detected cnt = %d", All-Zeros ) )
	#print( sprintf("%% detected cnt = %0.1f", (All-Zeros)/All*100.0) )
	#print( sprintf("sum = %d", sum(m)) )
	#print( sprintf("max read = %d", max(m)) ) 
	#print( sprintf("min read = %d", min(m)) ) 
	#print( sprintf("max AD = %d", max(colSums(m))) ) 
	#print( sprintf("min AD = %d", min(colSums(m))) ) 
	#print( sprintf("max BD = %d", max(rowSums(m))) ) 
	#print( sprintf("max BD = %d", min(rowSums(m))) ) 
	
	par(mfrow=c(2,2))
	plot(sort(rowSums(m)), xlab="Genes (sorted by sum of reads)", ylab="Marginal sum of reads", main="Bindinding Domain")
	plot(sort(colSums(m)), xlab="Genes (sorted by sum of reads)", ylab="Marginal sum of reads", main="Activation Domain")  
	hist(log(m,2), xlab="log2 Pair-reads", ylab="Frequency", main="Pair-reads distribution")   
	barplot(c(All-Zeros,Zeros),xlab="Cases",ylab="Gene pairs (count)",names=c("detected","not detected"),main="Coverage")
}


SaveMatrixPlot <- function( m, alpha, filepath, fontsize = 4 ){
	#png("/Users/jyang/Dropbox (CRG)/Collaboration/Mireia/NGS_Results/R75_MGj48/R71_76.S1.png",width=800,height=750)
	#simple_plotLog2Matrix(m1_3,1)
	#dev.off()
	png(filepath,width=800,height=750)    
	simple_plotLog2Matrix(m,alpha,fontsize)
	dev.off()                                   
}            
  
RefineMatrix0 <-function(m){  
	# These are never cloned into pENTR223: 
	# R75_41,WDR7,WD repeat domain 7
	# R75_72,PLEKHG7,pleckstrin homology domain containing, family G (with RhoGef domain) member 7
	         
	row_remove = c(43,75)
	row_index = setdiff(1:78, row_remove)
	
	col_remove = c(43,75)
	col_index = setdiff(1:78, col_remove)
	
	m2 = m[row_index,]    
	m3 = m2[,col_index]
	
	return (m3)
}


RefineMatrix <-function(m){  
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
	         
	row_remove = c(2,19,26,43,51,68,75)
	row_index = setdiff(1:78, row_remove)
	
	col_remove = c(43,75)
	col_index = setdiff(1:78, col_remove)
	
	m2 = m[row_index,]    
	m3 = m2[,col_index]
	
	return (m3)
}

RefineMatrix2 <-function(m){  
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
	# Remove toxic proteins (read count<1000):     
	#             AP2B1  BCL2L2   CTBP1   IKZF1   IKZF5     MME    PCNA PLEKHG7   RBPMS    STX4   VAMP2   VAMP3    WDR7     p53 
    # Read count: 422     145      76     889     314     282     159       0     772     470      96     333       0     540

	         
	row_remove = c(2,19,26,43,51,68,75)
	row_index = setdiff(1:78, row_remove)
	                       
	col_remove = c(2, 4, 19, 29, 30, 34, 41, 43, 53, 65, 73, 74, 75, 78 )
	col_index = setdiff(1:78, col_remove)
	
	m2 = m[row_index,]    
	m3 = m2[,col_index]
	
	return (m3)
} 

RefineMatrixP170_2 <-function(m){  
	# Removed these genes in the BD-vector
	# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
	#           
	# spike-in
	# 1:7, 1:7    
	#         
	# AA remove (1st)
	# fez1_p143		56
	# fez2_p144		57
	# tacc3_p108    160
	row_remove = c(1:7, 15, 56, 57, 160 )
	row_index = setdiff(1:173, row_remove)
	                       
	col_remove = c(1:7)
	col_index = setdiff(1:173, col_remove)
	
	m2 = m[row_index,]    
	m3 = m2[,col_index]
	
	return (m3)
}


RefineMatrixP170 <-function(m){  
	# Removed these genes in the BD-vector
	# Autoactivators, clone only in pAWH, remove from pBKWH since they are auto-activator
	#           
	# spike-in
	# 1:7, 1:7    
	#         
	# AA remove (1st)
	# fez1_p143		56
	# fez2_p144		57
	# tacc3_p108    160
	#
	# AA remove (2nd)
	# nin_p46		113      
	# ywhae_p151    171      
	# rpa2_p68 		150
	# map4_p106		98
	# kifap3_p132	86
	# rbm5_p66      147
	# cask_p142		15
	#
	# Partial domain
	# cdc42_p167	18
	# cdc42_p168	19    
	# mark2_p155	103
	# mark2_p156	104
	# mark2_p157	105
	# pard3b_p150   121
	# pard3b_p158   122
	# pard3b_p159   123
	# pard3b_p160   124
	# pard6g_p152   125
	# pard6g_p161   126
	# pard6g_p162   127
	# pard6g_p163   128
	# prkci_p153	135
	# prkci_p164    136
	# prkci_p165    137
	# prkci_p166    138
	# ywhae_p169	172
	# ywhae_p170	173
	
	row_remove = c(1:7, 15, 56, 57, 86, 98, 113, 147, 150, 160, 171, 18, 19, 103, 104, 105, 121, 122, 123, 124, 125, 126, 127, 128, 135, 136, 137, 138, 172, 173 )
	row_index = setdiff(1:173, row_remove)
	                       
	col_remove = c(1:7, 18, 19, 103, 104, 105, 121, 122, 123, 124, 125, 126, 127, 128, 135, 136, 137, 138, 172, 173 )
	col_index = setdiff(1:173, col_remove)
	
	m2 = m[row_index,]    
	m3 = m2[,col_index]
	
	return (m3)
}


# panel.smooth function is built in.
# panel.cor puts correlation in upper panels, size proportional to correlation
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}


# PPI_PerformanceCheck
PPI_PerformanceCheck <-function( no_selection_file = "", selection_file = "", PPI_Data = HIPPIE)
{
	m1_1 = drawMatrix(no_selection_file,0.5) # no selection ; 		46_WC_R12D12 
	m2_1 = drawMatrix(selection_file,0.5) # selection    ;    	46_AC_R12D12  
	is2_1 = InteractionScores(m1_1,m2_1,1.0)   
	boxplot( is2_1[ which( PPI_Data == 1 ) ], is2_1[ which( PPI_Data == 0 ) ] ) 
	t_test = t.test( is2_1[ which( PPI_Data == 1 ) ], is2_1[ which( PPI_Data == 0 ) ] ) 
	output = sprintf( "%f-%d-%d", -log10(t_test$p.value), sum(m1_1), sum(m2_1) )   
	
	labels = c( Ones(sum(PPI_Data == 1)), Zeros(sum(PPI_Data == 0)) )    
	predictions = c(is2_1[ which( PPI_Data == 1 ) ], is2_1[ which( PPI_Data == 0 ) ]) 
	roc_obj = roc( labels, predictions  )
	if( TRUE ){plot(roc_obj) }
	print ( auc(roc_obj) )
	print(output)
}           
                                                       
## Using Refined matrix
PPI_PerformanceCheck2 <-function( interaction_scores, PPI_Data = HIPPIE, Refined2 = FALSE)
{                                                                                        
	if (Refined2 == TRUE){
		is2_1 = RefineMatrix2( interaction_scores )
		PPI_1 = RefineMatrix2( PPI_Data )
	}else{
	   	is2_1 = RefineMatrix( interaction_scores )
		PPI_1 = RefineMatrix( PPI_Data )
	}  
	boxplot( is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ] ) 
	t_test = t.test( is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ] ) 
	output = sprintf( "%f__%f__%f", -log10(t_test$p.value), mean(is2_1[ which( PPI_1 == 1 ) ]), mean(is2_1[ which( PPI_1 == 0 ) ]) )   
	
	labels = c( Ones(sum(PPI_1 == 1)), Zeros(sum(PPI_1 == 0)) )    
	predictions = c(is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ]) 
	roc_obj = roc( labels, predictions  )
	if( TRUE ){plot(roc_obj) }
	print ( auc(roc_obj) )
	print(output)       
	return( c( auc(roc_obj)[1], -log10(t_test$p.value), mean(is2_1[ which( PPI_1 == 1 ) ]), mean(is2_1[ which( PPI_1 == 0 ) ]) ) )
}


## Using Refined matrix & but check for IS_cut, overlap_cut
PPI_PerformanceCheck3 <-function( interaction_scores, PPI_Data = HIPPIE, Refined2 = FALSE)
{                                                                                        
	if (Refined2 == TRUE){
		is2_1 = RefineMatrix2( interaction_scores )
		PPI_1 = RefineMatrix2( PPI_Data )
	}else{
	   	is2_1 = RefineMatrix( interaction_scores )
		PPI_1 = RefineMatrix( PPI_Data )
	}  
	boxplot( is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ] ) 
	t_test = t.test( is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ] ) 
	output = sprintf( "%f__%f__%f", -log10(t_test$p.value), mean(is2_1[ which( PPI_1 == 1 ) ]), mean(is2_1[ which( PPI_1 == 0 ) ]) )   
	  
	return( c( -1, -log10(t_test$p.value), mean(is2_1[ which( PPI_1 == 1 ) ]), mean(is2_1[ which( PPI_1 == 0 ) ]) ) )
}
        

## Using Refined matrix & but check for IS_cut, overlap_cut
PPI_PerformanceCheck_P170 <-function( interaction_scores, PPI_Data = HIPPIE, Refined2 = FALSE)
{                                                                                        
   	is2_1 = RefineMatrixP170( interaction_scores )
	PPI_1 = RefineMatrixP170( PPI_Data )
	boxplot( is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ] ) 
	t_test = t.test( is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ] ) 
	output = sprintf( "%f__%f__%f", -log10(t_test$p.value), mean(is2_1[ which( PPI_1 == 1 ) ]), mean(is2_1[ which( PPI_1 == 0 ) ]) )   
	  
	labels = c( Ones(sum(PPI_1 == 1)), Zeros(sum(PPI_1 == 0)) )    
	predictions = c(is2_1[ which( PPI_1 == 1 ) ], is2_1[ which( PPI_1 == 0 ) ]) 
	roc_obj = roc( labels, predictions  )
	if( TRUE ){plot(roc_obj) }
	print ( auc(roc_obj) )
	print(output)       
	return( c( auc(roc_obj)[1], -log10(t_test$p.value), mean(is2_1[ which( PPI_1 == 1 ) ]), mean(is2_1[ which( PPI_1 == 0 ) ]), mean(is2_1[ which( PPI_1 == 1 ) ])/ mean(is2_1[ which( PPI_1 == 0 ) ]) ) )
}

    
# Only works for symmetric matrix
SymmetryCut <- function( matrix, cutoff ){
	sym = t(matrix)
	output = ( matrix >= cutoff ) + ( sym >= cutoff )
	return(output)
}

# parameter scan
ParameterScanLogP <-function( file1 = "", file2 = "", len_start=10, len_end=30, len_step=2, count_start=3, count_end=10, count_step=1, PPI_Data = HIPPIE)
{
	for ( len in seq(len_start, len_end-1, by = len_step) ){                              
		for ( count in seq(count_start, count_end-1, by = count_step) ){                              
			
			no_selection_file = sprintf("%s.%d.%d.refine.txt", file1, count, len )
			selection_file = sprintf("%s.%d.%d.refine.txt", file2, count, len )
			m1_1 = drawMatrix(no_selection_file,0.5) # no selection ; 		46_WC_R12D12 
			m2_1 = drawMatrix(selection_file,0.5) # selection    ;    	46_AC_R12D12  
			is2_1 = InteractionScores(m1_1,m2_1,1.0)   
			boxplot( is2_1[ which( PPI_Data == 1 ) ], is2_1[ which( PPI_Data == 0 ) ], names = c("Intereactions","No Interactions"), ylab="Interaction Score (IS)" ) 
			t_test = t.test( is2_1[ which( PPI_Data == 1 ) ], is2_1[ which( PPI_Data == 0 ) ] ) 
			output = sprintf( "%d-%d-%f-%d-%d-%f-%f", len, count, -log10(t_test$p.value), sum(m1_1), sum(m2_1), mean(is2_1[ which( PPI_Data == 1 ) ]), mean(is2_1[ which( PPI_Data == 0 ) ]) )
			print(output)
		}
	}
}  


detailed_ppi <-function(filepath){            
	result = c()                       
	class_names=c("prc_prc","prc_p","p_prc","p_p","prc_v","prc_vrc","p_v","p_vrc","v_prc","vrc_prc","v_p","vrc_p","vrc_vrc","vrc_v","v_vrc","v_v") 
	output = sprintf( "%s.real.txt", filepath ); print( output )   
    result$prc_prc = drawMatrix( output, 0.5)
    output = sprintf( "%s.prc_p.txt", filepath ); print( output )
    result$prc_p = drawMatrix(output,0.5)       
    output = sprintf( "%s.p_prc.txt", filepath ); print( output )
    result$p_prc = drawMatrix(output,0.5)    
    output = sprintf( "%s.p_p.txt", filepath ); print( output )
    result$p_p = drawMatrix(output,0.5)       
    output = sprintf( "%s.prc_v.txt", filepath ); print( output ) 
    result$prc_v = drawMatrix(output,0.5)         
    output = sprintf( "%s.prc_vrc.txt", filepath ); print( output ) 
    result$prc_vrc = drawMatrix(output,0.5)    
    output = sprintf( "%s.p_v.txt", filepath ); print( output )  
    result$p_v = drawMatrix(output,0.5)
    output = sprintf( "%s.p_vrc.txt", filepath ); print( output )  
    result$p_vrc = drawMatrix(output,0.5)  
    output = sprintf( "%s.v_prc.txt", filepath ); print( output )  
    result$v_prc = drawMatrix(output,0.5)
    output = sprintf( "%s.vrc_prc.txt", filepath ); print( output )  
    result$vrc_prc = drawMatrix(output,0.5)
    output = sprintf( "%s.v_p.txt", filepath ); print( output )  
    result$v_p = drawMatrix(output,0.5)    
    output = sprintf( "%s.vrc_p.txt", filepath ) ; print( output ) 
    result$vrc_p = drawMatrix(output,0.5)
    output = sprintf( "%s.vrc_vrc.txt", filepath ); print( output )  
    result$vrc_vrc = drawMatrix(output,0.5)
    output = sprintf( "%s.vrc_v.txt", filepath ) ; print( output ) 
    result$vrc_v = drawMatrix(output,0.5)
    output = sprintf( "%s.v_vrc.txt", filepath ); print( output )  
    result$v_vrc = drawMatrix(output,0.5)
    output = sprintf( "%s.v_v.txt", filepath ) ; print( output ) 
    result$v_v = drawMatrix(output,0.5)
    index = 1
    for (x in result){ sprintf("%s\t%d",class_names[index],sum(x)); index = index + 1 }
    return ( result )
}

## Correlation matrix with p-values. See http://goo.gl/nahmV for documentation of this function
cor.prob <- function (X, dfr = nrow(X) - 2) {
  R <- cor(X, use="pairwise.complete.obs")
  above <- row(R) < col(R)
  r2 <- R[above]^2
  Fstat <- r2 * dfr/(1 - r2)
  R[above] <- 1 - pf(Fstat, 1, dfr)
  R[row(R) == col(R)] <- NA
  R
}

## Use this to dump the cor.prob output to a 4 column matrix
## with row/column indices, correlation, and p-value.
## See StackOverflow question: http://goo.gl/fCUcQ
flattenSquareMatrix <- function(m) {
  if( (class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.") 
  if(!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(i = rownames(m)[row(m)[ut]],
             j = rownames(m)[col(m)[ut]],
             cor=t(m)[ut],
             p=m[ut])
}


RowIndex <-function(m, name){
	return( which( rownames(m) == name ) )
}     

ColIndex <-function(m, name){
	return( which( colnames(m) == name ) )
}     
     

# IRE/IRP      

CheckRPIControls <-function( m1, m2 ){
	ire = RowIndex(m1,"IRE")   # RNA -- row
	irp = ColIndex(m1, "IRP")  # Protein -- col
	print( m2[ ire, irp ] )
	#[1] 376249
	print( m1[ ire, irp ] )
	#[1] 476
	print( m2[ ire, irp ]/m1[ ire, irp ] )
	#[1] 790.4391                         
}       

GetRowColName <- function( m, index ){
	#row_i = ceiling( index / dim(m)[1] )
	#col_i = index %% dim(m)[1]
	col_i = ceiling( index / dim(m)[1] )
	row_i = index %% dim(m)[1]
	
	#print( rownames(m)[row_i] )
	#print( colnames(m)[col_i] )    
	#print( m[index] )         
	rnames = rownames(m)[row_i]
	cnames = colnames(m)[col_i]
	values = m[index]         
	return( data.frame(rnames, row_i, cnames, col_i, values, index) )
	return ( cbind( rnames, row_i, cnames, col_i, values, index ) )
}

## Circular diagram
CircularPlot <-function(m){
	##
	##install packages if not already done so (uncomment)
	##
	# install.packages("circlize")
	# install.packages("dplyr")
	# install.packages("tidyr")

	##
	##read in table and define matrix (m) and reference data.frame (df1)
	##

	df0 <- read.table(system.file("science", "region_custom.txt", package = "migest"), skip=2, stringsAsFactors=FALSE)

	#select and rename labels
	library("dplyr")

	df1 <- df0 %>% select(1:3) %>% rename(order = V1, rgb = V2, region = V3) %>% mutate(region = gsub("_", " ", region))

	#flow matrix
	m <- as.matrix(df0[,-(1:3)]/1e06)
	dimnames(m) <- list(orig = df1$region, dest = df1$region)


	#drop small flows
	m[m<=quantile(m,0.65)]<-0


	#sort regions and create colours
	library("tidyr")
	df1 <- df1 %>% arrange(order) %>% separate(rgb, c("r","g","b")) %>% mutate(col = rgb(r, g, b, max=255), max = rowSums(m)+colSums(m))


	#plot using chordDiagram
	library("circlize")
	circos.clear()
	par(mar = rep(0, 4), cex=0.9)

	circos.par(start.degree = 90, gap.degree = 4)

	chordDiagram(x = m, directional = 1, order = df1$region, 
	              grid.col = df1$col, annotationTrack = "grid", 
	              transparency = 0.25,  annotationTrackHeight = c(0.1, 0.1),
	              diffHeight  = -0.04)


	#add in labels and axis
	circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
	   xlim = get.cell.meta.data("xlim")
	   sector.index = get.cell.meta.data("sector.index")
	   circos.text(mean(xlim), 2.5, sector.index, facing = "bending")
	   circos.axis("top", major.at = seq(0, max(xlim)), minor.ticks=1, labels.away.percentage = 0.2, labels.niceFacing = FALSE )
	 }, bg.border = NA)


	circos.clear()
}
       

## CircosPlot diagram
CircosPlot <-function(m){
    dev.off();                     
	RCircos.Set.Plot.Area();
	RCircos.Chromosome.Ideogram.Plot();
	
	#data(RCircos.Gene.Label.Data);
	#name.col <- 4;
	#side <- "in";
	#track.num <- 1;
	#RCircos.Gene.Connector.Plot(RCircos.Gene.Label.Data,track.num, side);
	#track.num <- 2;
	#RCircos.Gene.Name.Plot(RCircos.Gene.Label.Data,+ name.col,track.num, side);
    
	data(RCircos.Link.Data);
	track.num <- 1;
	#RCircos.Link.Plot(RCircos.Link.Data, track.num, TRUE);
	data(RCircos.Ribbon.Data);           
	#RCircos.Link.Data                   
	RCircos.Ribbon.Plot(ribbon.data=RCircos.Link.Data,track.num=1, by.chromosome=FALSE, twist=FALSE);
	RCircos.Ribbon.Plot(ribbon.data=RCircos.Ribbon.Data,track.num=1, by.chromosome=FALSE, twist=FALSE);
	
	#dev.off();
	
}         

## CircosPlot diagram
CircosParPlot <-function(m){                      
	pdf("CircosParPlot.v1.pdf")
	Par.CytoBandIdeogram <- read.delim( "~/Desktop/Par.txt")      
	Par.Domain.Info <- read.delim( "~/Desktop/ParDomainInfo.txt")
	Par.Domain.Info[,2] = Par.Domain.Info[,2] * 100000 
	Par.Domain.Info[,3] = Par.Domain.Info[,3] * 100000
	Par.Ribbon.Data <- read.delim( "~/Desktop/Par.Ribon.txt")  
	Par.Ribbon.Data[,2] = Par.Ribbon.Data[,2] * 100000 
	Par.Ribbon.Data[,3] = Par.Ribbon.Data[,3] * 100000
	Par.Ribbon.Data[,5] = Par.Ribbon.Data[,5] * 100000 
	Par.Ribbon.Data[,6] = Par.Ribbon.Data[,6] * 100000

	Par.CytoBandIdeogram[,2] = Par.CytoBandIdeogram[,2] * 100000 
	Par.CytoBandIdeogram[,3] = Par.CytoBandIdeogram[,3] * 100000

	chr.exclude <- NULL;
	cyto.info <- Par.CytoBandIdeogram;
	tracks.inside <- 5;
	tracks.outside <- 0;           
	track.height <- 4;
	RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside, track.height);
 
	dev.off();                     
	RCircos.Set.Plot.Area();
	RCircos.Chromosome.Ideogram.Plot();
                                              
    #RCircos.Link.Plot(Par.Ribbon.Data, track.num=3, TRUE);
    side = "in"      
	RCircos.Gene.Connector.Plot(Par.Domain.Info,track.num=1, side);                                   
	RCircos.Gene.Name.Plot(Par.Domain.Info, name.col=4,track.num=2, side);
	RCircos.Ribbon.Plot(ribbon.data=Par.Ribbon.Data,track.num=3, by.chromosome=FALSE, twist=FALSE);                
	dev.off()
}

 


MaskMatrix <-function(m, mask, mask_value = 0, set_value = 0.0){
	x = m
	x[mask==mask_value] = set_value
	return(x)
}      


GenomicMutationSmoothCurves <- function(){
	# adapted     
	library(readr)
	ACC_somatic_mutation_hist <- read_delim("~/Downloads/ACC_somatic_mutation_hist.txt", "\t", escape_double = FALSE, trim_ws = TRUE)  
	
	chr1 = ACC_somatic_mutation_hist[ ACC_somatic_mutation_hist$Chr == 1, ]
	index = sort(chr1$start_pos, index.return = TRUE )$ix
	chr1.DF = data.frame(Pos=chr1$start_pos[index],Mutations=chr1$Num[index])
	Loess.Fit <- loess(Mutations ~ Pos, data = chr1.DF )       
	plot(chr1.DF$Pos, chr1.DF$Mutations)
	Predicted.Values <- predict(Loess.Fit, newdata = data.frame(Pos = seq(1,max(chr1$start_pos+10000000),100000) ) )   
	lines(seq(1,max(chr1$start_pos+10000000),100000),Predicted.Values,col="red",lwd=3)
	
	
	Table.For.Plotting= data.frame(Scores = chr1.DF$Mutations, Index= chr1.DF$Pos )
    
	# original code
	Loess.DF <- data.frame(Score = Table.For.Plotting$Scores, Index = 1:length(Table.For.Plotting$Scores))
	Loess.Fit <- loess(Score ~ Index, data = Loess.DF, span = 0.01)
	j <- order(Loess.DF$Index)
	Predicted.Values <- predict(Loess.Fit, newdata = data.frame(Index = seq(1,length(Table.For.Plotting$Scores),0.1)) )
	lines(seq(1,length(Table.For.Plotting$Scores),0.1),Predicted.Values,col="red",lwd=3)
                                   
}

GammaMixture <-function(data, return.output = FALSE, LogScale = TRUE){      
	# data = refine_bm1_Q
	# se1_W       
	set.seed(47)                                                    
	if (LogScale == TRUE){
		#output = gammamixEM(as.vector(data+1),lambda = NULL, alpha = c(1,10), beta = c(100,100), k = 2, epsilon = 0.01, maxit = 10000, maxrestarts=200, verb = TRUE)       
		output = gammamixEM(as.vector(data+1),lambda = NULL, alpha = c(0.1,10), beta = c(100,100), k = 2, epsilon = 0.001, maxit = 10000, maxrestarts=200, verb = TRUE)  
		print((output))     
	}else{
		output = gammamixEM(as.vector(data),lambda = NULL, alpha = c(0.1,10), beta = c(100,100), k = 2, epsilon = 0.01, maxit = 10000, maxrestarts=200, verb = TRUE) 
		print((output))              
	}
	plot(log10(as.vector(data+1)),output$posterior[,1])
	plot(log10(as.vector(data+1)),output$posterior[,2])	 
	data[output$posterior[,2] < 0.5] = 0         
	simple_plotLog2Matrix(data,1,7)   
	print(output)              
	if ( return.output ){
		return ( output )    
	}   
	else{
		return (data)
	} 
	
	# gm_out_EGFR1_Q = GammaMixture(bEGFR1_Q,return.output = TRUE)                                                                                                                                                      
	# Theta <- list(pdf1 = "gamma", theta1.1 = c(1/100, 1/100, 1/100,1/100), theta2.1 = c(200, 400, 600, 800))
	# gamma1est <- REBMIX(Dataset = as.vector(bEGFR1_Q+1), Preprocessing = "Parzen window", cmax = 8, Criterion = c("AIC", "BIC"), pdf = "gamma")
}


GammaMixture2 <-function(data, return.output = FALSE, LogScale = FALSE){      
	set.seed(47)                                                    
	if (LogScale == TRUE){
		output = gammamixEM(log10(as.vector(data+1)),lambda = c(0.9,0.1), verb = TRUE)  
		print((output))     
	}else{
		output = gammamixEM(as.vector(data+1),lambda = c(0.9,0.1), verb = TRUE) 
		print((output))              
	}
	plot(log10(as.vector(data+1)),output$posterior[,1])
	points(log10(as.vector(data+1)),output$posterior[,2],col="red")	   
	means = output$gamma.pars[1,] *  output$gamma.pars[2,]
	data[output$x < means[1]] = 0  
	data[output$posterior[,2] < 0.5] = 0         
	simple_plotLog2Matrix(data,1,7)   
	print(output)              
	if ( return.output ){
		return ( output )    
	}   
	else{
		return (data)
	} 
	
}



NormalMixture <-function(data, return.output = FALSE){    
	# data = refine_bm1_Q
	# se1_W
	output = normalmixEM(as.vector(data+1))       
	plot(log10(as.vector(data+1)),output$posterior[,1])
	plot(log10(as.vector(data+1)),output$posterior[,2])	    
	
	if ( output$mu[2] > output$mu[1] ){
		data[output$posterior[,2] < 0.9] = 0
	}else{
		data[output$posterior[,1] < 0.9] = 0
	}
	simple_plotLog2Matrix(data,1,7)                 
	if ( return.output ){
		return ( output )    
	}   
	else{
		return (data)
	}
}
  
# Remove 0 values, and gaussian mixture to find signal
NormalMixture2 <-function(data, show.plot = FALSE, threshold = 0.5, non_zero_cut = 0){        
	set.seed(47)   
	# data = refine_bm1_Q
	# se1_W             
	output2 = data
	non_zero_data = data[which(data>non_zero_cut)]
	#output = normalmixEM(log2(as.vector(non_zero_data)))  
	output = normalmixEM(as.vector(non_zero_data))
	#output$x = log2(output$x)
	print( output$mu )
	print( output$sigma )
	if (show.plot){
		plot(output,which=2)
	}
	if (output$mu[2] > output$mu[1]){
		non_zero_data[output$posterior[,2]< threshold] = 0
		output2[which(data>0)] = non_zero_data 
    }else{
		non_zero_data[output$posterior[,1]< threshold] = 0
		output2[which(data>0)] = non_zero_data
    } 
    
	#plot(output,which=2)  
	#print(output)    
	#plot(log10(as.vector(non_zero_data)),output$posterior[,1])
	#plot(log10(as.vector(non_zero_data)),output$posterior[,2])	    
	
	#if ( output$mu[2] > output$mu[1] ){
	#	non_zero_data[output$posterior[,2] < 0.9] = 0
	#}else{
	#	non_zero_data[output$posterior[,1] < 0.9] = 0
	#}
	
	return ( output2 )    
}    

# Remove 0 values, and gaussian mixture to find signal
NormalMixture2_Detail <-function(data, show.plot = FALSE, threshold = 0.5, non_zero_cut = 0){        
	set.seed(47)           
	output2 = data
	non_zero_data = data[which(data>non_zero_cut)]
	output = normalmixEM(as.vector(non_zero_data))
	print( output$mu )
	print( output$sigma )
	if (show.plot){
		plot(output,which=2)
	}
	if (output$mu[2] > output$mu[1]){
		non_zero_data[output$posterior[,2]< threshold] = 0
		output2[which(data>0)] = non_zero_data 
    }else{
		non_zero_data[output$posterior[,1]< threshold] = 0
		output2[which(data>0)] = non_zero_data
    } 
	return ( output )    
}
      
   
FilterMatrix <-function( m, threshold = 100 ){
	m[which(m<threshold)] = 0.0
	return(m)  
}   
                           
chippy1 <- function(x,factor){      
	xx = 10^x          
	o = dnorm(xx, mean=56.65694, sd=56.98714) 
	#temp_index = which(x < 0)
	#x[temp_index] = 1
	#o = dnorm(log10(x), mean=56.65694, sd=56.98714)    
	#o[temp_index] = 0
	
	return (o*factor)             
}
                 
chippy2 <- function(x, factor){      
	xx = 10^x          
	o = dnorm(xx, mean=5477.78462, sd=6154.90499 ) 
	#temp_index = which(x < 0)
	#x[temp_index] = 1
	#o = dnorm(log10(x), mean=56.65694, sd=56.98714)    
	#o[temp_index] = 0
	
	return (o*factor)             
}


# Fig 2 - noise filter
DrawFig2NoiseFilter <- function(){       
	pdf("noise_filter.pdf")
	par(mfrow=c(3,1))
	set.seed(47)     
	data = bm2_Q[5:18,5:18]
	non_zero_data = sort(data[which(data>0)])    
	hist( log10(non_zero_data), breaks=50, probability = TRUE, ylim=c(0,1)  )
	output = normalmixEM(as.vector(non_zero_data))  
	# > output$mu   
	# [1]   56.65694 5477.78462
	# > output$sigma   
	# [1]   56.98714 6154.90499
	curve(chippy1(x,factor=1000), 0, 5, add=TRUE, col="blue") 
	curve(chippy2(x,factor=1000), 0, 5, add=TRUE, col="red")     
	noise_points = cbind(log10(output$x),output$posterior[,1])
	plot(log10(output$x),output$posterior[,1])
	plot(log10(output$x),output$posterior[,2])
	dev.off()	                 
	
	pdf("noise_filter.2.pdf")
	par(mfrow=c(1,1))    
	hist( log10(non_zero_data/sum(non_zero_data)), breaks=30, probability = TRUE, ylim=c(0,1)  )            
	lines(log10(output$x/sum(non_zero_data)),output$posterior[,1],col="blue")
	lines(log10(output$x/sum(non_zero_data)),output$posterior[,2],col="red")
	dev.off()
	  
	# Log10
	# set.seed(47)     
	# data = bm2_Q[5:18,5:18]
	# non_zero_data = data[which(data>0)]    
	# hist( log10(non_zero_data), breaks=30, probability = TRUE, ylim=c(0,1)  )
	# output2 = normalmixEM(log10(as.vector(non_zero_data)))  
	# plot(output2,which=2,breaks=15,xlim=c(0,5))    
}

# Remove 0 values, and log-gaussian mixture to find signal
LogNormalMixture <-function(data, show.plot = FALSE, threshold = 0.5){        
	set.seed(47)           
	
	output2 = data
	non_zero_data = data[which(data>0)]
	output = normalmixEM(log10(as.vector(non_zero_data)))
	print( output$mu )
	print( output$sigma )
	if (show.plot){
		plot(output,which=2)
	}
	if (output$mu[2] > output$mu[1]){
		non_zero_data[output$posterior[,2]< threshold] = 0
		output2[which(data>0)] = non_zero_data 
    }else{
		non_zero_data[output$posterior[,1]< threshold] = 0
		output2[which(data>0)] = non_zero_data
    } 
	
	return ( output2 )    
}

                        
                             
RowWiseCorrelation <-function( A, B ){
	cA <- A - rowMeans(A)
	cB <- B - rowMeans(B)
	sA <- sqrt(rowMeans(cA^2))
	sB <- sqrt(rowMeans(cB^2))

	return ( rowMeans(cA * cB) / (sA * sB) )
}


Matrix2List <-function( data, cutoff=2, filename=NA ){
	if ( is.na(filename) == FALSE ){
		sink( filename )
	}
	row_names = rownames(data)
	col_names = colnames(data)
	for (i in 1:length(row_names)) {
		for (j in 1:length(col_names)) {
			if( data[i,j] >= cutoff ){                
				print ( sprintf( "%s %s %d", row_names[i], col_names[j], data[i,j] ) )
			}
		}
	}                  
	if ( is.na(filename) == FALSE ){
		sink()
	} 
}
                           
 
ColumnWiseInteractionScores <- function(m1,m2,alpha, fontsize=9, rowNoralization=FALSE, filename=NA){      
	rcnt = dim(m1)[1]            
	nm1 = m1 #/ sum(m1)
	nm2 = m2 #/ sum(m2)
	
	# 1. normalized by column sum in selection matrix
	colSumsMat2 = ( (1:rcnt * 0 + 1) %*% t(colSums(nm2)) )
	colNorm_m2 = nm2 /  colSumsMat2	           
	
	# 2. normalized by row sum in no-selection matrix
	colSumsMat1 = ( (1:rcnt * 0 + 1) %*% t(colSums(nm1)) )
	colNorm_m1 = nm1 /  colSumsMat1                          
	colNorm_m1[ which(colNorm_m1=="NaN")] = 0                  
	  
	rowFreq = rowSums(colNorm_m1) / length(colnames(colNorm_m1))
	
    if ( rowNoralization == FALSE ){
		simple_plotMatrix(colNorm_m2,0,fontsize,filename=filename)    
		return( colNorm_m2 )
    }else{
		refined = colNorm_m2/ rowFreq  
		refined[ which(refined == "NaN") ] = 0 
		refined[ is.na(refined) ] = 0 
		refined[ which(refined == Inf) ] = max( refined[ which(refined != Inf) ] ) + 1.0
	   
		simple_plotMatrix(refined,0,fontsize,filename=filename)               
	}
	return (refined)
}


HomoDimerTest <-function(is, removeToxicGenes = FALSE, alpha = 0.1, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE){     
	HomoPositive = c(7,17,29,49,52,53,55,61,70)  
	HomoPositiveWithoutToxic = c(7,17,49,52,55,61,70)

	HomoNegative = c(4,10,13,14,22,24,27,30,31,34,39,41,45,48,58,65,67,69,71,73,74,78)
	HomoNegativeWithoutToxic = c(10,13,14,22,24,27,31,39,45,48,58,65,67,69,71)
	                                                                             
	HomoBioGridPos = c(6,7,8,14,17,19,20,21,26,29,30,33,34,35,37,39,41,42,46,49,51,52,53,55,57,58,59,61,65,69,70,71,72,76,78)
	HomoBioGridNoInfo = c(1,2,3,4,5,9,10,11,12,13,15,16,18,22,23,24,25,27,28,31,32,36,38,40,43,44,45,47,48,50,54,56,60,62,63,64,66,67,68,73,74,75,77)
	      
	ToxicIndex = c(2, 4, 19, 29, 30, 34, 41, 43, 51, 52, 53, 58, 65, 69, 71, 73, 74, 75, 78)     
	
	if (removeToxicGenes == TRUE){
		HomoBioGridPos = setdiff(HomoBioGridPos, ToxicIndex)       # remove toxic genes
		HomoBioGridNoInfo = setdiff(HomoBioGridNoInfo, ToxicIndex)    # remove toxic genes
    }
                          
    ratios = diag(is)
	       
	labels = c( Ones(length(HomoPositiveWithoutToxic)), Zeros(length(HomoNegativeWithoutToxic)) )    
    predictions = c(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic]) 
    roc_obj = roc( labels, predictions  )   

	if (roc_plot1) { plot(roc_obj) }
	print ( auc(roc_obj) )
	#print( wilcox.test( predictions, labels ) )

	if (box_plot1) { boxplot(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic], names=c("HomoPositiveW/OToxic","HomoNegativeW/OToxic")) }    
	
	print ( t.test(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic]) )      
	print ( wilcox.test(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic]) )

	labels = c( Ones(length(HomoBioGridPos)), Zeros(length(HomoBioGridNoInfo)) )    
    predictions = c(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo]) 
    roc_obj = roc( labels, predictions  )
    if (roc_plot2) { plot(roc_obj) }
	print ( auc(roc_obj) )
	
	if ( box_plot2 ) { boxplot(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo], names=c("HomoBioGridPos","HomoBioGridNoInfo")) }    

	print ( t.test(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo]) )      
	print ( wilcox.test(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo]) )
	
	return ( ratios )
}


SetUnionAnalysis <-function(m1, m2){ 
	#
	# m1 and m2 is boolean matrix
	#
	# SetUnionAnalysis( as.vector(is2_E2_WHA1>2), as.vector(is2_E1_Q>2) )
	#
	diff1 = m1 - m2
	diff2 = m2 - m1
	add = m1 + m2
	print( sum(m1) ) # m1 alone
	print( sum(m2) ) # m2 alone
	print( length(which(diff1==1)) ) # m1 only
	print( length(which(diff2==1)) ) # m2 only
	print( length(which(add==2)) ) # common
	print( length(which(add>0)) ) # union
}

PlotHomoDimerTest <- function(sa2_W,sa2_A,sa2_Q,sa3_W,sa3_A,sa3_Q,sa4_SW,sa4_S4A,sa4_SQ){
	par(mfrow=c(3,8))
	HomoDimerNormalization(sa2_W, sa2_A, method="colSums", roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")
	
	par(mfrow=c(3,8))
	HomoDimerNormalization(sa2_W, sa2_A, method="diag", roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")
	
		
	par(mfrow=c(3,8))
	HomoDimerTest(sa2_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerTest(sa2_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerTest(sa3_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerTest(sa3_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerTest(sa4_SW, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerTest(sa4_SW, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")	    
	
	par(mfrow=c(3,8))
	HomoDimerTest(sa2_A, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerTest(sa2_Q, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)   
	mtext("Homo dimer exp-2Q")
	HomoDimerTest(sa3_A, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)       
	mtext("Homo dimer exp-3A")
	HomoDimerTest(sa3_Q, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)   
	mtext("Homo dimer exp-3Q")
	HomoDimerTest(sa4_S4A, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-4A4")
	HomoDimerTest(sa4_SQ, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-4Q")
}                             

HomoDimerNormalization <-function(m1, m2, method="colSums", removeToxicGenes = FALSE, alpha = 0.1, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE){     
	HomoPositive = c(7,17,29,49,52,53,55,61,70)  
	HomoPositiveWithoutToxic = c(7,17,49,52,55,61,70)

	HomoNegative = c(4,10,13,14,22,24,27,30,31,34,39,41,45,48,58,65,67,69,71,73,74,78)
	HomoNegativeWithoutToxic = c(10,13,14,22,24,27,31,39,45,48,58,65,67,69,71)
	                                                                             
	HomoBioGridPos = c(6,7,8,14,17,19,20,21,26,29,30,33,34,35,37,39,41,42,46,49,51,52,53,55,57,58,59,61,65,69,70,71,72,76,78)
	HomoBioGridNoInfo = c(1,2,3,4,5,9,10,11,12,13,15,16,18,22,23,24,25,27,28,31,32,36,38,40,43,44,45,47,48,50,54,56,60,62,63,64,66,67,68,73,74,75,77)
	      
	ToxicIndex = c(2, 4, 19, 29, 30, 34, 41, 43, 51, 52, 53, 58, 65, 69, 71, 73, 74, 75, 78)     
	
	if (removeToxicGenes == TRUE){
		HomoBioGridPos = setdiff(HomoBioGridPos, ToxicIndex)       # remove toxic genes
		HomoBioGridNoInfo = setdiff(HomoBioGridNoInfo, ToxicIndex)    # remove toxic genes
    }

	if (method == "diag"){
    	ratios = diag(m2+alpha)/diag(m1+alpha)
	}else if( method == "colSums"){       
    	ratios = diag(m2+alpha)/colSums(m1+alpha)       
    }
           
	labels = c( Ones(length(HomoPositiveWithoutToxic)), Zeros(length(HomoNegativeWithoutToxic)) )    
    predictions = c(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic]) 
    roc_obj = roc( labels, predictions  )   

	if( roc_plot1 ){ plot(roc_obj) }
	print ( auc(roc_obj) )
	#print( wilcox.test( predictions, labels ) )
                          
	print( ratios[HomoPositiveWithoutToxic] )
	print( ratios[HomoNegativeWithoutToxic] )
	if( box_plot1 ){ boxplot(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic], names=c("HomoPositiveWithoutToxic","HomoNegativeWithoutToxic")) }   
	
	print ( t.test(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic]) )      
	print ( wilcox.test(ratios[HomoPositiveWithoutToxic], ratios[HomoNegativeWithoutToxic]) )

	labels = c( Ones(length(HomoBioGridPos)), Zeros(length(HomoBioGridNoInfo)) )    
    predictions = c(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo]) 
    roc_obj = roc( labels, predictions  )
    if( roc_plot2 ){plot(roc_obj) }
	print ( auc(roc_obj) )
	
	if( box_plot2 ){ boxplot(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo], names=c("HomoBioGridPos","HomoBioGridNoInfo")) }    

	print ( t.test(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo]) )      
	print ( wilcox.test(ratios[HomoBioGridPos], ratios[HomoBioGridNoInfo]) )
	
	return ( ratios )
}   

KRP_Performance <- function(m1, filename = NA ){
	labels = matrix(0,nrow=17,ncol=19)   
	diag( labels ) = 1
    roc_obj = roc( as.vector(labels), as.vector(m1)  )  
  	if ( is.na(filename) ){
		plot(roc_obj)  
	}else{
		pdf(filename)
    	plot(roc_obj)                      
		dev.off()
	}
	return ( roc_obj )
}


# Calculate sample and pair complexity from no selection matrix
SampleComplexity <- function(m1, filename = NA ){              
	print( BasicStats(RefineMatrix(m1)) )
	print( BasicStats(RefineMatrix(m1+t(m1))) )
}


GenerateRandomMatrix <-function(m, count = 6084){
	set.seed(43)

	nrow = dim(m)[1]
	ncol = dim(m)[2]

	#col_index = c()
	#row_index = c()
	col_index = 1:sum(m) * 0
	row_index = 1:sum(m) * 0

	cnt_index = 1
	for (i in 1:nrow){
		for (j in 1:ncol){
			cnt = m[i,j]
			if ( cnt != 0 ){
				#row_index = c( row_index, rep(i,cnt) )
				#col_index = c( col_index, rep(j,cnt) )
				row_index[cnt_index:(cnt_index+cnt-1)] = rep(i,cnt)
				col_index[cnt_index:(cnt_index+cnt-1)] = rep(j,cnt)
				cnt_index = cnt_index + cnt
			}
		}
	}

	random_i_list = sample(1:length(col_index), count, replace = TRUE, prob = NULL)

	random_m = matrix(0,nrow,ncol)

	for (i in random_i_list){
		row = row_index[ i ]
		col = col_index[ i ]
		random_m[ row, col ] = random_m[ row, col ] + 1
	}

	return( random_m )
}


SubSamplingMatrixTest <-function( m, matrix_size = 6084 )
{
    
    m_01X = GenerateRandomMatrix(m, matrix_size/10 )
    m_05X = GenerateRandomMatrix(m, matrix_size/2 )
    m_1X = GenerateRandomMatrix(m, matrix_size )
    m_2X = GenerateRandomMatrix(m, matrix_size*2 )
    m_4X = GenerateRandomMatrix(m, matrix_size*4 )
    m_8X = GenerateRandomMatrix(m, matrix_size*8 )
    m_16X = GenerateRandomMatrix(m, matrix_size*16 )
    m_32X = GenerateRandomMatrix(m, matrix_size*32 )
    m_64X = GenerateRandomMatrix(m, matrix_size*64 )
    m_128X = GenerateRandomMatrix(m, matrix_size*128 )
    m_256X = GenerateRandomMatrix(m, matrix_size*256 )

    print( sum(m)/matrix_size ) # 242.0 x

    r_01X = CorMatrix( m, m_01X )
    r_05X = CorMatrix( m, m_05X )
    r_1X = CorMatrix( m, m_1X )
    r_2X = CorMatrix( m, m_2X )
    r_4X = CorMatrix( m, m_4X )
    r_8X = CorMatrix( m, m_8X )
    r_16X = CorMatrix( m, m_16X )
    r_32X = CorMatrix( m, m_32X )
    r_64X = CorMatrix( m, m_64X )
    r_128X = CorMatrix( m, m_128X )
    r_256X = CorMatrix( m, m_256X )

    plot( c(0.1,0.5,1,2,4,8,16,32,64,128,256), c(r_01X,r_05X,r_1X,r_2X,r_4X,r_8X,r_16X,r_32X,r_64X,r_128X,r_256X),xlab="Coverage X",ylab="Correlation (r)" )
    return (c(r_01X,r_05X,r_1X,r_2X,r_4X,r_8X,r_16X,r_32X,r_64X,r_128X,r_256X))
}
                         
#coefficient of variation
CoV <-function(m){
	return (sd(m)/mean(m))
}     

######################################################################################################################################
#
# Data read
#                        
LOAD_OLD_DATA = FALSE
       
                                    
if (FALSE){
	# PPI database
	#HIPPIE <- ReadMatrix( "~/Dropbox (CRG)/Code/jslib/Mireia/HIPPIE/HIPPIE.R75.txt" )
	HIPPIE <- ReadMatrix( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/HIPPIE/HIPPIE.R75.txt" )
	pheatmap((HIPPIE), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=6)

	#BioGrid <- ReadMatrix( "~/Dropbox (CRG)/Code/jslib/Mireia/BioGrid/BioGrid.R75.txt" )
	BioGrid <- ReadMatrix( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/BioGrid/BioGrid.R75.txt" )
	pheatmap((BioGrid), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=7)   

	KnownPPI = HIPPIE | BioGrid        
	KnownPPI2 = BioGrid | HIPPIE # | BioGrid         
	
	write.table(melt(KnownPPI2),file="KnownPPI.melted.txt",sep="\t")
	
	# PPI database
	#HIPPIE_P170 <- ReadMatrixP170( "~/Dropbox (CRG)/Code/jslib/Mireia/HIPPIE/HIPPIE.P170.txt" )        
	HIPPIE_P170 <- ReadMatrixP170( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/HIPPIE/HIPPIE.P170.txt" )        
	pheatmap((HIPPIE_P170), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=4)

	#BioGrid_P170 <- ReadMatrixP170( "~/Dropbox (CRG)/Code/jslib/Mireia/BioGrid/BioGrid.P170.txt" )       
	BioGrid_P170 <- ReadMatrixP170( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/BioGrid/BioGrid.P170.txt" )       
	pheatmap((BioGrid_P170), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=4)   
	      
	KnownPPI_P170 = HIPPIE_P170 | BioGrid_P170
	pheatmap(KnownPPI_P170, cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=4)   

                                                                                                                 
	# ALL
	#BioGrid_P170_ALL <- ReadMatrixP170( "~/Dropbox (CRG)/Code/jslib/Mireia/BioGrid/BioGrid.P170.ALL.txt" )    
	BioGrid_P170_ALL <- ReadMatrixP170( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/BioGrid/BioGrid.P170.ALL.txt" )    
	pheatmap((BioGrid_P170_ALL), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=4)
}


## OLD data
if (LOAD_OLD_DATA == TRUE){   
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )   
	bm1_W_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  

	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; ACCRD         "sum = 2024754"  
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL_old = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt",0.5) # selection; 53_QL4RD 18715	"sum = 734766"	 

	xx1p =PairInteractionScores(bm1_W_old,bm1_A_old,1,7)
	xx1 = InteractionScores(bm1_W_old,bm1_A_old,1,7)  
	
	xx2p =PairInteractionScores(bm2_W_old,bm2_A_old,1,7)
	xx2 = InteractionScores(bm2_W_old,bm2_A_old,1,7)
	
	xx3p =PairInteractionScores(bm3_W_old,bm3_Q_old,1,7)
	xx3 = InteractionScores(bm3_W_old,bm3_Q_old,1,7)
}



## first
if (FALSE){
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )     
	bm1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt.new.2",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt.new.2",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt.new.2",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt.new.2",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  
	
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt.new.2",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt.new.2",0.5) # selection; ACCRD         "sum = 2024754"
	bm2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt.new.2",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt.new.2",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt.new.2",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt.new.2",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt.new.2",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      
	
	# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
	bm4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt.new.2",0.5) # Roth; -W;                 "sum = 1491017"
	bm4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt.new.2",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	bm4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt.new.2",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
	bm4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt.new.2",0.5) # Roth; -W/-H/-A;      "sum = 618177"
	          
	# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;                                                                        
	bm5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt.new.2",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
	       
	# 2017-10-30_MiSeq; Roth75 - autoactivator test
	S1_BWD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S1_BWD_R1.blastn.cnt.txt")
	S2_BA2D_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S2_BA2D_R1.blastn.cnt.txt")
	S3_BQD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S3_BQD_R1.blastn.cnt.txt")
	S4_AWD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S4_AWD_R2.blastn.cnt.txt")

	# 2017-11-03_MiSeq; Roth75 - technical repeat
	S1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S1_W.ppi.txt.new.2",0.5)
	S2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S2_Q.ppi.txt.new.2",0.5)
	S3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S3_W.ppi.txt.new.2",0.5)
	S4_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S4_Q.ppi.txt.new.2",0.5)	
	
	
	# For homo-dimer                             
	View( cbind(colSums(bm1_W), colSums(bm2_W), diag(bm1_W), diag(bm2_W), diag(bm1_A), diag( bm1_RA), diag(bm2_A), diag(bm2_Q), diag(bm2_Q)) )   
	View( cbind(colSums(bm3_W), colSums(bm4_SW), diag(bm3_W), diag(bm4_SW), diag(bm3_A), diag( bm3_Q), diag(bm4_S8A), diag(bm4_S4A), diag(bm4_SQ)) )   
	
	PlotHomoDimerTest(bm2_W,bm2_A,bm2_Q,bm3_W,bm3_A,bm3_Q,bm4_SW,bm4_S4A,bm4_SQ)   
} 

GetSymmetricMatrix <- function(m, cutoff){
	# m = TRUE/FALSE matrix
	output = matrix(0, nrow=dim(m)[1], ncol=dim(m)[2] )
	output[ m >= cutoff ] = 1
	output[ t(m) >= cutoff ] = output[ t(m) >= cutoff  ] + 1
	return( output )
}                 

SymmetricSum <- function(m, cutoff = 0){
	m2 = m
	m2[ which(m2 < cutoff ) ] = 0
	m2 = m2 + t(m2)
	return (m2)
}                     

## first
if (FALSE){
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )     
	bm1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  
	
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; ACCRD         "sum = 2024754"
	bm2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      
	
	# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
	bm4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt",0.5) # Roth; -W;                 "sum = 1491017"
	bm4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	bm4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
	bm4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt",0.5) # Roth; -W/-H/-A;      "sum = 618177"
	          
	# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;                                                                        
	bm5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
	       
	# 2017-10-30_MiSeq; Roth75 - autoactivator test
	S1_BWD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S1_BWD_R1.blastn.cnt.txt")
	S2_BA2D_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S2_BA2D_R1.blastn.cnt.txt")
	S3_BQD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S3_BQD_R1.blastn.cnt.txt")
	S4_AWD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S4_AWD_R2.blastn.cnt.txt")

	# 2017-11-03_MiSeq; Roth75 - technical repeat
	S1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S1_W.ppi.txt",0.5)
	S2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S2_Q.ppi.txt",0.5)
	S3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S3_W.ppi.txt",0.5)
	S4_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S4_Q.ppi.txt",0.5)	
	
	
	# For homo-dimer                             
	View( cbind(colSums(bm1_W), colSums(bm2_W), diag(bm1_W), diag(bm2_W), diag(bm1_A), diag( bm1_RA), diag(bm2_A), diag(bm2_Q), diag(bm2_Q)) )   
	View( cbind(colSums(bm3_W), colSums(bm4_SW), diag(bm3_W), diag(bm4_SW), diag(bm3_A), diag( bm3_Q), diag(bm4_S8A), diag(bm4_S4A), diag(bm4_SQ)) )   
	
	PlotHomoDimerTest(bm2_W,bm2_A,bm2_Q,bm3_W,bm3_A,bm3_Q,bm4_SW,bm4_S4A,bm4_SQ)   
} 


## Blastn_All_Ref
# Use all ref seq
if (FALSE){
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )     
	bm1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn_All_Ref/17543_S1.ppi.txt",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn_All_Ref/17544_S2.ppi.txt",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn_All_Ref/17545_S3.ppi.txt",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn_All_Ref/17546_S4.ppi.txt",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  
	
	is1_A = InteractionScores(bm1_W,bm1_A,1.0)        
	PPI_PerformanceCheck2(is1_A, HIPPIE)
	
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn_All_Ref/S1.ppi.txt",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn_All_Ref/S2.ppi.txt",0.5) # selection; ACCRD         "sum = 2024754"
	bm2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn_All_Ref/S3.ppi.txt",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn_All_Ref/S1.ppi.txt",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn_All_Ref/S2.ppi.txt",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn_All_Ref/S3.ppi.txt",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn_All_Ref/S7.ppi.txt",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      
	
	# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
	bm4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn_All_Ref/S64_SWD.ppi.txt",0.5) # Roth; -W;                 "sum = 1491017"
	bm4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn_All_Ref/S64_SA4D.ppi.txt",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	bm4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn_All_Ref/S64_SA8D.ppi.txt",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
	bm4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn_All_Ref/S64_SQD.ppi.txt",0.5) # Roth; -W/-H/-A;      "sum = 618177"
	          
	# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;                                                                        
	bm5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn_All_Ref/S52.ppi.txt",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
	       
	# 2017-10-30_MiSeq; Roth75 - autoactivator test
	S1_BWD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn_All_Ref/S1_BWD_R1.blastn.cnt.txt")
	S2_BA2D_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn_All_Ref/S2_BA2D_R1.blastn.cnt.txt")
	S3_BQD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn_All_Ref/S3_BQD_R1.blastn.cnt.txt")
	S4_AWD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn_All_Ref/S4_AWD_R2.blastn.cnt.txt")

	# 2017-11-03_MiSeq; Roth75 - technical repeat
	S1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn_All_Ref/S1_W.ppi.txt",0.5)
	S2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn_All_Ref/S2_Q.ppi.txt",0.5)
	S3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn_All_Ref/S3_W.ppi.txt",0.5)
	S4_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn_All_Ref/S4_Q.ppi.txt",0.5)	
	
	
	# For homo-dimer                             
	View( cbind(colSums(bm1_W), colSums(bm2_W), diag(bm1_W), diag(bm2_W), diag(bm1_A), diag( bm1_RA), diag(bm2_A), diag(bm2_Q), diag(bm2_Q)) )   
	View( cbind(colSums(bm3_W), colSums(bm4_SW), diag(bm3_W), diag(bm4_SW), diag(bm3_A), diag( bm3_Q), diag(bm4_S8A), diag(bm4_S4A), diag(bm4_SQ)) )   
	
	PlotHomoDimerTest(bm2_W,bm2_A,bm2_Q,bm3_W,bm3_A,bm3_Q,bm4_SW,bm4_S4A,bm4_SQ)   
}


#checkDiagonal_new2 <-function(){   
if (FALSE){
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )     
	bm1_W2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt.new.2",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt.new.2",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt.new.2",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt.new.2",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  
	
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt.new.2",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt.new.2",0.5) # selection; ACCRD         "sum = 2024754"
	bm2_Q2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt.new.2",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt.new.2",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt.new.2",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt.new.2",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt.new.2",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      
	
	# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
	bm4_SW2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt.new.2",0.5) # Roth; -W;                 "sum = 1491017"
	bm4_S4A2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt.new.2",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	bm4_S8A2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt.new.2",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
	bm4_SQ2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt.new.2",0.5) # Roth; -W/-H/-A;      "sum = 618177"
	          
	# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;                                                                        
	bm5_SW2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt.new.2",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
	
	sort( diag(bm1_A2) )
	sort( diag(bm2_A2) )
	sort( diag(bm3_A2) )
	sort( diag(bm4_S4A2) )
	
	# For homo-dimer
	View( cbind(colSums(bm3_W2), colSums(bm4_SW2), diag(bm3_W2), diag(bm4_SW2), diag(bm3_A2), diag( bm3_Q2), diag(bm4_S8A2), diag(bm4_S4A2), diag(bm4_SQ2)) )     
}   
                                                                                                                                        

## Sebastian Mauer way
if (FALSE){
	sm4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SWD.ppi.txt",0.5) # Roth; -W;                 "sum = 1491017"
	sm4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SA4D.ppi.txt",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	sm4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SA8D.ppi.txt",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"                                                                                                       
	sm4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SQD.ppi.txt",0.5)

	sm4_SW_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SWD_full.ppi.txt",0.5) # Roth; -W;                 "sum = 1491017"
	sm4_S4A_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SA4D_full.ppi.txt",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	sm4_S8A_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SA8D_full.ppi.txt",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"                                                                                                       
	sm4_SQ_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Sebastian/S64_SQD_full.ppi.txt",0.5)	
	
}   

#========================
# LoadRothData
#========================
LoadRothData <- function(){
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )     
	bm1_W <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt.new",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt.new",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt.new",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt.new",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  
	
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt.new",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt.new",0.5) # selection; ACCRD         "sum = 2024754"
	bm2_Q <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt.new",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt.new",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt.new",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt.new",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt.new",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      
	
	# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
	bm4_SW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt.new",0.5) # Roth; -W;                 "sum = 1491017"
	bm4_S4A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt.new",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	bm4_S8A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt.new",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
	bm4_SQ <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt.new",0.5) # Roth; -W/-H/-A;      "sum = 618177"
		
	# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;                                                                        
	bm5_SW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt.new",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
	                         
	# 2017-11-03_MiSeq; Roth75 - technical repeat
	bm6_SW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S1_W.ppi.txt.new",0.5)
	bm6_SQ <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S2_Q.ppi.txt.new",0.5)
	bm7_SW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S3_W.ppi.txt.new",0.5)
	bm7_SQ <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S4_Q.ppi.txt.new",0.5)

	# 2017-10-30_MiSeq; Roth75 - autoactivator test     (with Empty vectors)
	S1_BWD_R1.cnt <<- read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S1_BWD_R1.blastn.cnt.txt",header=FALSE)
	S2_BA2D_R1.cnt <<- read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S2_BA2D_R1.blastn.cnt.txt",header=FALSE)	# AA test
	S3_BQD_R1.cnt <<- read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S3_BQD_R1.blastn.cnt.txt",header=FALSE)		# AA test
	S4_AWD_R2.cnt <<- read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S4_AWD_R2.blastn.cnt.txt",header=FALSE) 		# toxic protein           

	is2_bm1_A <<- InteractionScores(bm1_W,NormalMixture2(bm1_A),1,9 )	  
	is2_bm1_RA <<- InteractionScores(bm1_W,NormalMixture2(bm1_RA),1,9 )	  
	is2_bm2_A <<- InteractionScores(bm2_W,NormalMixture2(bm2_A),1,9 )	  
	is2_bm2_Q <<- InteractionScores(bm2_W,NormalMixture2(bm2_Q),1,9 )	  
	is2_bm3_A <<- InteractionScores(bm3_W,NormalMixture2(bm3_A),1,9 )	  
	is2_bm3_Q <<- InteractionScores(bm3_W,NormalMixture2(bm3_Q),1,9 )	  
	is2_bm3_QL <<- InteractionScores(bm3_W,NormalMixture2(bm3_QL),1,9 )	  
	is2_bm4_S4A <<- InteractionScores(bm4_SW,NormalMixture2(bm4_S4A),1,9 )	  
	is2_bm4_S8A <<- InteractionScores(bm4_SW,NormalMixture2(bm4_S8A),1,9 )	  
	is2_bm4_SQ <<- InteractionScores(bm4_SW,NormalMixture2(bm4_SQ),1,9 )	  
	is2_bm6_SQ <<- InteractionScores(bm6_SW,NormalMixture2(bm6_SQ),1,9 )	  
	is2_bm7_SQ <<- InteractionScores(bm7_SW,NormalMixture2(bm7_SQ),1,9 )	 

	# PPI database
	#HIPPIE <<- ReadMatrix( "~/Dropbox (CRG)/Code/jslib/Mireia/HIPPIE/HIPPIE.R75.txt" )    
	HIPPIE <<- ReadMatrix( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/HIPPIE/HIPPIE.R75.txt" )    
	
	pheatmap((HIPPIE), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=6)

	#BioGrid <<- ReadMatrix( "~/Dropbox (CRG)/Code/jslib/Mireia/BioGrid/BioGrid.R75.txt" )
	BioGrid <<- ReadMatrix( "/Volumes/Element_Backup/Dropbox_CRG/Code/jslib/Mireia/BioGrid/BioGrid.R75.txt" )
	pheatmap((BioGrid), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=7)   
	      
	KnownPPI <<- HIPPIE | BioGrid        
}

LoadRefinedRothData <- function(){
	## First do LoadRothData
	refine_bm1_W <<- RefineMatrix(bm1_W) 
	refine_bm1_A <<- RefineMatrix(bm1_A)  
	refine_bm1_RW <<- RefineMatrix(bm1_RW) 
	refine_bm1_RA <<- RefineMatrix(bm1_RA) 
	
	refine_bm2_W <<- RefineMatrix(bm2_W)
	refine_bm2_A <<- RefineMatrix(bm2_A)
	refine_bm2_Q <<- RefineMatrix(bm2_Q)
	
	refine_bm3_W <<- RefineMatrix(bm3_W)
	refine_bm3_A <<- RefineMatrix(bm3_A)
	refine_bm3_Q <<- RefineMatrix(bm3_Q)
	refine_bm3_QL <<- RefineMatrix(bm3_QL)
	
	refine_bm4_SW <<- RefineMatrix(bm4_SW)
	refine_bm4_S4A <<- RefineMatrix(bm4_S4A)
	refine_bm4_S8A <<- RefineMatrix(bm4_S8A)
	refine_bm4_SQ <<- RefineMatrix(bm4_SQ)
	
	refine_bm5_SW <<- RefineMatrix(bm5_SW)	  

	refine_bm6_SW <<- RefineMatrix(bm6_SW)	  
	refine_bm6_SQ <<- RefineMatrix(bm6_SQ)

	refine_bm7_SW <<- RefineMatrix(bm7_SW)	  
	refine_bm7_SQ <<- RefineMatrix(bm7_SQ)	 

	refine_selection_roth <<- list( refine_bm1_A, refine_bm1_RA,
		refine_bm2_A, refine_bm2_Q, 
		refine_bm3_A, refine_bm3_Q, refine_bm3_QL,
		refine_bm4_S4A, refine_bm4_S8A, refine_bm4_SQ,
		refine_bm6_SQ, refine_bm7_SQ
	)	  
}

LoadRefinedRoth_IS <- function(){
    is_refine_bm1_A <<- InteractionScores(refine_bm1_W,refine_bm1_A,1,9 )	
    is_refine_bm1_RA <<- InteractionScores(refine_bm1_RW,refine_bm1_RA,1,9 )	
    
	is_refine_bm2_A <<- InteractionScores(refine_bm2_W,refine_bm2_A,1,9 )	
	is_refine_bm2_Q <<- InteractionScores(refine_bm2_W,refine_bm2_Q,1,9 )	

	is_refine_bm3_A <<- InteractionScores(refine_bm3_W,refine_bm3_A,1,9 )	
	is_refine_bm3_Q <<- InteractionScores(refine_bm3_W,refine_bm3_Q,1,9 )
	is_refine_bm3_QL <<- InteractionScores(refine_bm3_W,refine_bm3_QL,1,9 )

	is_refine_bm4_S4A <<- InteractionScores(refine_bm4_SW,refine_bm4_S4A,1,9 )	
	is_refine_bm4_S8A <<- InteractionScores(refine_bm4_SW,refine_bm4_S8A,1,9 )
	is_refine_bm4_SQ <<- InteractionScores(refine_bm4_SW,refine_bm4_SQ,1,9 )

	is_refine_bm6_SQ <<- InteractionScores(refine_bm6_SW,refine_bm6_SQ,1,9 )	
	is_refine_bm7_SQ <<- InteractionScores(refine_bm7_SW,refine_bm7_SQ,1,9 )	

	is_roth <<- list( is_refine_bm1_A, is_refine_bm1_RA,
		is_refine_bm2_A, is_refine_bm2_Q, 
		is_refine_bm3_A, is_refine_bm3_Q, is_refine_bm3_QL,
		is_refine_bm4_S4A, is_refine_bm4_S8A, is_refine_bm4_SQ,
		is_refine_bm6_SQ, is_refine_bm7_SQ
	)	
}

LoadRefinedRoth_IS2 <- function(){
    is2_refine_bm1_A <<- InteractionScores(refine_bm1_W,NormalMixture2(refine_bm1_A),1,9 )	
    is2_refine_bm1_RA <<- InteractionScores(refine_bm1_RW,NormalMixture2(refine_bm1_RA),1,9 )	
    
	is2_refine_bm2_A <<- InteractionScores(refine_bm2_W,NormalMixture2(refine_bm2_A),1,9 )	
	is2_refine_bm2_Q <<- InteractionScores(refine_bm2_W,NormalMixture2(refine_bm2_Q),1,9 )	

	is2_refine_bm3_A <<- InteractionScores(refine_bm3_W,NormalMixture2(refine_bm3_A),1,9 )	
	is2_refine_bm3_Q <<- InteractionScores(refine_bm3_W,NormalMixture2(refine_bm3_Q),1,9 )
	is2_refine_bm3_QL <<- InteractionScores(refine_bm3_W,NormalMixture2(refine_bm3_QL),1,9 )

	is2_refine_bm4_S4A <<- InteractionScores(refine_bm4_SW,NormalMixture2(refine_bm4_S4A),1,9 )	
	is2_refine_bm4_S8A <<- InteractionScores(refine_bm4_SW,NormalMixture2(refine_bm4_S8A),1,9 )
	is2_refine_bm4_SQ <<- InteractionScores(refine_bm4_SW,NormalMixture2(refine_bm4_SQ),1,9 )

	is2_refine_bm6_SQ <<- InteractionScores(refine_bm6_SW,NormalMixture2(refine_bm6_SQ),1,9 )	
	is2_refine_bm7_SQ <<- InteractionScores(refine_bm7_SW,NormalMixture2(refine_bm7_SQ),1,9 )	

	is2_roth <<- list( is2_refine_bm1_A, is2_refine_bm1_RA,
		is2_refine_bm2_A, is2_refine_bm2_Q, 
		is2_refine_bm3_A, is2_refine_bm3_Q, is2_refine_bm3_QL,
		is2_refine_bm4_S4A, is2_refine_bm4_S8A, is2_refine_bm4_SQ,
		is2_refine_bm6_SQ, is2_refine_bm7_SQ
	)
	
}

LoadRefinedRoth_IS_AAC <- function(){
	# Auto-activation corrected
    isaac_refine_bm1_A <<- AutoActivationCorrection( is_refine_bm1_A )	
    isaac_refine_bm1_RA <<- AutoActivationCorrection( is_refine_bm1_RA )	
    
	isaac_refine_bm2_A <<- AutoActivationCorrection( is_refine_bm2_A )	
	isaac_refine_bm2_Q <<- AutoActivationCorrection( is_refine_bm2_Q )	

	isaac_refine_bm3_A <<- AutoActivationCorrection( is_refine_bm3_A )	
	isaac_refine_bm3_Q <<- AutoActivationCorrection( is_refine_bm3_Q )
	isaac_refine_bm3_QL <<- AutoActivationCorrection( is_refine_bm3_QL )

	isaac_refine_bm4_S4A <<- AutoActivationCorrection( is_refine_bm4_S4A )	
	isaac_refine_bm4_S8A <<- AutoActivationCorrection( is_refine_bm4_S8A )
	isaac_refine_bm4_SQ <<- AutoActivationCorrection( is_refine_bm4_SQ )

	isaac_refine_bm6_SQ <<- AutoActivationCorrection( is_refine_bm6_SQ )	
	isaac_refine_bm7_SQ <<- AutoActivationCorrection( is_refine_bm7_SQ )	

	isaac_roth <<- list( isaac_refine_bm1_A, isaac_refine_bm1_RA,
		isaac_refine_bm2_A, isaac_refine_bm2_Q, 
		isaac_refine_bm3_A, isaac_refine_bm3_Q, isaac_refine_bm3_QL,
		isaac_refine_bm4_S4A, isaac_refine_bm4_S8A, isaac_refine_bm4_SQ,
		isaac_refine_bm6_SQ, isaac_refine_bm7_SQ
	)
}

LoadRefinedRoth_IS_AAC2 <- function(){
	# Auto-activation corrected
    isaac2_refine_bm1_A <<- AutoActivationCorrection2( is_refine_bm1_A )	
    isaac2_refine_bm1_RA <<- AutoActivationCorrection2( is_refine_bm1_RA )	
    
	isaac2_refine_bm2_A <<- AutoActivationCorrection2( is_refine_bm2_A )	
	isaac2_refine_bm2_Q <<- AutoActivationCorrection2( is_refine_bm2_Q )	

	isaac2_refine_bm3_A <<- AutoActivationCorrection2( is_refine_bm3_A )	
	isaac2_refine_bm3_Q <<- AutoActivationCorrection2( is_refine_bm3_Q )
	isaac2_refine_bm3_QL <<- AutoActivationCorrection2( is_refine_bm3_QL )

	isaac2_refine_bm4_S4A <<- AutoActivationCorrection2( is_refine_bm4_S4A )	
	isaac2_refine_bm4_S8A <<- AutoActivationCorrection2( is_refine_bm4_S8A )
	isaac2_refine_bm4_SQ <<- AutoActivationCorrection2( is_refine_bm4_SQ )

	isaac2_refine_bm6_SQ <<- AutoActivationCorrection2( is_refine_bm6_SQ )	
	isaac2_refine_bm7_SQ <<- AutoActivationCorrection2( is_refine_bm7_SQ )	

	isaac2_roth <<- list( isaac_refine_bm1_A, isaac_refine_bm1_RA,
		isaac2_refine_bm2_A, isaac2_refine_bm2_Q, 
		isaac2_refine_bm3_A, isaac2_refine_bm3_Q, isaac2_refine_bm3_QL,
		isaac2_refine_bm4_S4A, isaac2_refine_bm4_S8A, isaac2_refine_bm4_SQ,
		isaac2_refine_bm6_SQ, isaac2_refine_bm7_SQ
	)
}

LoadRefinedRoth_IS_AAC3 <- function(){
	# Auto-activation corrected
    isaac3_refine_bm1_A <<- AutoActivationCorrection3( is_refine_bm1_A )	
    isaac3_refine_bm1_RA <<- AutoActivationCorrection3( is_refine_bm1_RA )	
    
	isaac3_refine_bm2_A <<- AutoActivationCorrection3( is_refine_bm2_A )	
	isaac3_refine_bm2_Q <<- AutoActivationCorrection3( is_refine_bm2_Q )	

	isaac3_refine_bm3_A <<- AutoActivationCorrection3( is_refine_bm3_A )	
	isaac3_refine_bm3_Q <<- AutoActivationCorrection3( is_refine_bm3_Q )
	isaac3_refine_bm3_QL <<- AutoActivationCorrection3( is_refine_bm3_QL )

	isaac3_refine_bm4_S4A <<- AutoActivationCorrection3( is_refine_bm4_S4A )	
	isaac3_refine_bm4_S8A <<- AutoActivationCorrection3( is_refine_bm4_S8A )
	isaac3_refine_bm4_SQ <<- AutoActivationCorrection3( is_refine_bm4_SQ )

	isaac3_refine_bm6_SQ <<- AutoActivationCorrection3( is_refine_bm6_SQ )	
	isaac3_refine_bm7_SQ <<- AutoActivationCorrection3( is_refine_bm7_SQ )	

	isaac3_roth <<- list( isaac3_refine_bm1_A, isaac3_refine_bm1_RA,
		isaac3_refine_bm2_A, isaac3_refine_bm2_Q, 
		isaac3_refine_bm3_A, isaac3_refine_bm3_Q, isaac3_refine_bm3_QL,
		isaac3_refine_bm4_S4A, isaac3_refine_bm4_S8A, isaac3_refine_bm4_SQ,
		isaac3_refine_bm6_SQ, isaac3_refine_bm7_SQ
	)
}



LoadRefinedRoth_IS2_AAC <- function(){
	# Auto-activation corrected
    is2aac_refine_bm1_A <<- AutoActivationCorrection( is2_refine_bm1_A )	
    is2aac_refine_bm1_RA <<- AutoActivationCorrection( is2_refine_bm1_RA )	
    
	is2aac_refine_bm2_A <<- AutoActivationCorrection( is2_refine_bm2_A )	
	is2aac_refine_bm2_Q <<- AutoActivationCorrection( is2_refine_bm2_Q )	

	is2aac_refine_bm3_A <<- AutoActivationCorrection( is2_refine_bm3_A )	
	is2aac_refine_bm3_Q <<- AutoActivationCorrection( is2_refine_bm3_Q )
	is2aac_refine_bm3_QL <<- AutoActivationCorrection( is2_refine_bm3_QL )

	is2aac_refine_bm4_S4A <<- AutoActivationCorrection( is2_refine_bm4_S4A )	
	is2aac_refine_bm4_S8A <<- AutoActivationCorrection( is2_refine_bm4_S8A )
	is2aac_refine_bm4_SQ <<- AutoActivationCorrection( is2_refine_bm4_SQ )

	is2aac_refine_bm6_SQ <<- AutoActivationCorrection( is2_refine_bm6_SQ )	
	is2aac_refine_bm7_SQ <<- AutoActivationCorrection( is2_refine_bm7_SQ )	

	is2aac_roth <<- list( is2aac_refine_bm1_A, is2aac_refine_bm1_RA,
		is2aac_refine_bm2_A, is2aac_refine_bm2_Q, 
		is2aac_refine_bm3_A, is2aac_refine_bm3_Q, is2aac_refine_bm3_QL,
		is2aac_refine_bm4_S4A, is2aac_refine_bm4_S8A, is2aac_refine_bm4_SQ,
		is2aac_refine_bm6_SQ, is2aac_refine_bm7_SQ
	)
}


LoadRefinedRoth_IS2_AAC2 <- function(){
	# Auto-activation corrected
    is2aac2_refine_bm1_A <<- AutoActivationCorrection2( is2_refine_bm1_A )	
    is2aac2_refine_bm1_RA <<- AutoActivationCorrection2( is2_refine_bm1_RA )	
    
	is2aac2_refine_bm2_A <<- AutoActivationCorrection2( is2_refine_bm2_A )	
	is2aac2_refine_bm2_Q <<- AutoActivationCorrection2( is2_refine_bm2_Q )	

	is2aac2_refine_bm3_A <<- AutoActivationCorrection2( is2_refine_bm3_A )	
	is2aac2_refine_bm3_Q <<- AutoActivationCorrection2( is2_refine_bm3_Q )
	is2aac2_refine_bm3_QL <<- AutoActivationCorrection2( is2_refine_bm3_QL )

	is2aac2_refine_bm4_S4A <<- AutoActivationCorrection2( is2_refine_bm4_S4A )	
	is2aac2_refine_bm4_S8A <<- AutoActivationCorrection2( is2_refine_bm4_S8A )
	is2aac2_refine_bm4_SQ <<- AutoActivationCorrection2( is2_refine_bm4_SQ )

	is2aac2_refine_bm6_SQ <<- AutoActivationCorrection2( is2_refine_bm6_SQ )	
	is2aac2_refine_bm7_SQ <<- AutoActivationCorrection2( is2_refine_bm7_SQ )	

	is2aac2_roth <<- list( is2aac2_refine_bm1_A, is2aac2_refine_bm1_RA,
		is2aac2_refine_bm2_A, is2aac2_refine_bm2_Q, 
		is2aac2_refine_bm3_A, is2aac2_refine_bm3_Q, is2aac2_refine_bm3_QL,
		is2aac2_refine_bm4_S4A, is2aac2_refine_bm4_S8A, is2aac2_refine_bm4_SQ,
		is2aac2_refine_bm6_SQ, is2aac2_refine_bm7_SQ
	)
}



LoadRefinedRoth_IS2_AAC3 <- function(){
	# Auto-activation corrected
    is2aac3_refine_bm1_A <<- AutoActivationCorrection3( is2_refine_bm1_A )	
    is2aac3_refine_bm1_RA <<- AutoActivationCorrection3( is2_refine_bm1_RA )	
    
	is2aac3_refine_bm2_A <<- AutoActivationCorrection3( is2_refine_bm2_A )	
	is2aac3_refine_bm2_Q <<- AutoActivationCorrection3( is2_refine_bm2_Q )	

	is2aac3_refine_bm3_A <<- AutoActivationCorrection3( is2_refine_bm3_A )	
	is2aac3_refine_bm3_Q <<- AutoActivationCorrection3( is2_refine_bm3_Q )
	is2aac3_refine_bm3_QL <<- AutoActivationCorrection3( is2_refine_bm3_QL )

	is2aac3_refine_bm4_S4A <<- AutoActivationCorrection3( is2_refine_bm4_S4A )	
	is2aac3_refine_bm4_S8A <<- AutoActivationCorrection3( is2_refine_bm4_S8A )
	is2aac3_refine_bm4_SQ <<- AutoActivationCorrection3( is2_refine_bm4_SQ )

	is2aac3_refine_bm6_SQ <<- AutoActivationCorrection3( is2_refine_bm6_SQ )	
	is2aac3_refine_bm7_SQ <<- AutoActivationCorrection3( is2_refine_bm7_SQ )	

	is2aac3_roth <<- list( is2aac3_refine_bm1_A, is2aac3_refine_bm1_RA,
		is2aac3_refine_bm2_A, is2aac3_refine_bm2_Q, 
		is2aac3_refine_bm3_A, is2aac3_refine_bm3_Q, is2aac3_refine_bm3_QL,
		is2aac3_refine_bm4_S4A, is2aac3_refine_bm4_S8A, is2aac3_refine_bm4_SQ,
		is2aac3_refine_bm6_SQ, is2aac3_refine_bm7_SQ
	)
}


## Final used
if (FALSE){    
	#PPI_PerformanceCheck("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt.new.2", "/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt.new.2", HIPPIE)
	
	LoadRothData()
                       
    Roth75_AA_Test = data.frame( Index = 1:78, W_R1 = S1_BWD_R1.cnt[,2], A2_R1 = S2_BA2D_R1.cnt[,2], Q_R1 = S3_BQD_R1.cnt[,2], A2_R1_Ratio = S2_BA2D_R1.cnt[,2]/S1_BWD_R1.cnt[,2], Q_R1_Ratio = S3_BQD_R1.cnt[,2]/S1_BWD_R1.cnt[,2], A_R2 = S4_AWD_R2.cnt[,2] )
    row.names( Roth75_AA_Test ) = S1_BWD_R1.cnt[,1]      

	ps( "Roth75_AA_Test.1.ps" ); plot(Roth75_AA_Test$A2_R1/Roth75_AA_Test$W_R1); dev.off();
    ps( "Roth75_AA_Test.2.ps" ); plot(Roth75_AA_Test$Q_R1/Roth75_AA_Test$W_R1); dev.off();
	ps( "Roth75_AA_Test.3.ps" ); plot(Roth75_AA_Test$A_R2); dev.off();     
		 
	S1_BWD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S1_BWD_R2.blastn.cnt.txt")
	S2_BA2D_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S2_BA2D_R2.blastn.cnt.txt")
	S3_BQD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S3_BQD_R2.blastn.cnt.txt")
	S4_AWD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Blastn/S4_AWD_R1.blastn.cnt.txt")	
    
	sS1_BWD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S1_BWD_R1.blastn.cnt.txt")
	sS2_BA2D_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S2_BA2D_R1.blastn.cnt.txt")
	sS3_BQD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S3_BQD_R1.blastn.cnt.txt")
	sS4_AWD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S4_AWD_R2.blastn.cnt.txt")

	sS1_BWD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S1_BWD_R2.blastn.cnt.txt")
	sS2_BA2D_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S2_BA2D_R2.blastn.cnt.txt")
	sS3_BQD_R2.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S3_BQD_R2.blastn.cnt.txt")
	sS4_AWD_R1.cnt = read.delim("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-30_MiSeq/Sebastian/S4_AWD_R1.blastn.cnt.txt")
		        
	se_S1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S1_W.ppi.txt",0.5)
	se_S2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S2_Q.ppi.txt",0.5)  
	se_S3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S3_W.ppi.txt",0.5)  
	se_S4_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S4_Q.ppi.txt",0.5)  
	

	se_S1_W_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S1_W_full.ppi.txt",0.5)
	se_S2_Q_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S2_Q_full.ppi.txt",0.5)  
	se_S3_W_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S3_W_full.ppi.txt",0.5)  
	se_S4_Q_full = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Sebastian/S4_Q_full.ppi.txt",0.5)

	is2_bm1_A = InteractionScores(bm1_W,NormalMixture2(bm1_A),1,9 )	  
	is2_bm1_RA = InteractionScores(bm1_W,NormalMixture2(bm1_RA),1,9 )	  
	is2_bm2_A = InteractionScores(bm2_W,NormalMixture2(bm2_A),1,9 )	  
	is2_bm2_Q = InteractionScores(bm2_W,NormalMixture2(bm2_Q),1,9 )	  
	is2_bm3_A = InteractionScores(bm3_W,NormalMixture2(bm3_A),1,9 )	  
	is2_bm3_Q = InteractionScores(bm3_W,NormalMixture2(bm3_Q),1,9 )	  
	is2_bm3_QL = InteractionScores(bm3_W,NormalMixture2(bm3_QL),1,9 )	  
	is2_bm4_S4A = InteractionScores(bm4_SW,NormalMixture2(bm4_S4A),1,9 )	  
	is2_bm4_S8A = InteractionScores(bm4_SW,NormalMixture2(bm4_S8A),1,9 )	  
	is2_bm4_SQ = InteractionScores(bm4_SW,NormalMixture2(bm4_SQ),1,9 )	  
	is2_bm6_SQ = InteractionScores(bm6_SW,NormalMixture2(bm6_SQ),1,9 )	  
	is2_bm7_SQ = InteractionScores(bm7_SW,NormalMixture2(bm7_SQ),1,9 )	  
	            
	
	##================ Fig 2 Last (AUC cutoff)                     
	# single drop out
	is2_bm2_A
	is2_bm3_A
	is2_bm4_S4A 
	is2_bm4_S8A
	
	# double drop out   
	is2_bm2_Q
	is2_bm3_Q  
	is2_bm3_QL
	is2_bm4_SQ
	is2_bm6_SQ  
	is2_bm7_SQ
	
	Cutoff_Evaluation1 <- function( IS_cut = 1, overlap_cut = 1 ){ 
		#PPI_PerformanceCheck2(is2_bm2_A+is2_bm3_A+is2_bm4_S4A+is2_bm4_S8A,KnownPPI)            
	    p = PPI_PerformanceCheck2((is2_bm2_A>IS_cut)+(is2_bm3_A>IS_cut)+(is2_bm4_S4A>IS_cut)+(is2_bm4_S8A>IS_cut),KnownPPI)            
		return(p)
	}   
	
	Cutoff_Evaluation2 <- function( IS_cut = 1, overlap_cut = 1 ){ 
		#PPI_PerformanceCheck2(is2_bm2_A+is2_bm3_A+is2_bm4_S4A+is2_bm4_S8A,KnownPPI)            
	    p = PPI_PerformanceCheck2((is2_bm2_Q>IS_cut)+(is2_bm3_Q>IS_cut)+(is2_bm3_QL>IS_cut)+(is2_bm4_SQ>IS_cut)+(is2_bm6_SQ>IS_cut)+(is2_bm7_SQ>IS_cut),KnownPPI)            
		return(p)
	}  

	Cutoff_Evaluation3 <- function( IS_cut = 1, overlap_cut = 1 ){ 
		p = PPI_PerformanceCheck3( ((is2_bm2_A>IS_cut)+(is2_bm3_A>IS_cut)+(is2_bm4_S4A>IS_cut)+(is2_bm4_S8A>IS_cut))>=overlap_cut,KnownPPI)            
		return(p)
	}   
	
	Cutoff_Evaluation4 <- function( IS_cut = 1, overlap_cut = 1 ){          
	    p = PPI_PerformanceCheck3( ((is2_bm2_Q>IS_cut)+(is2_bm3_Q>IS_cut)+(is2_bm3_QL>IS_cut)+(is2_bm4_SQ>IS_cut)+(is2_bm6_SQ>IS_cut)+(is2_bm7_SQ>IS_cut))>=overlap_cut,KnownPPI)            
		return(p)
	}
	                                                                                                                          
	performance_test1 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for single drop out
	performance_test2 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for double drop out
	for (i in 1:50){       
		IS_cutoff = i / 10
		performance_test1[i,1] = IS_cutoff   
		performance_test2[i,1] = IS_cutoff                                                                                     
		p1 = Cutoff_Evaluation1(i/10)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
		performance_test1[i,2:5] = p1
		p2 = Cutoff_Evaluation2(i/10)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
		performance_test2[i,2:5] = p2
	}
	
	par(mfrow=c(2,2))
	plot(performance_test1[,1],performance_test1[,2])
	plot(performance_test1[,1],performance_test1[,3])
	plot(performance_test1[,1],performance_test1[,4])
	plot(performance_test1[,1],performance_test1[,5])
	par(mfrow=c(2,2))
	plot(performance_test2[,1],performance_test2[,2])
	plot(performance_test2[,1],performance_test2[,3])
	plot(performance_test2[,1],performance_test2[,4])
	plot(performance_test2[,1],performance_test2[,5])


	performance_test1 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for single drop out
	performance_test2 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for double drop out
	for (i in 1:5){       
		IS_cutoff = i 
		for (j in 1:5){     
			Overlap_cutoff = j
			performance_test1[i+j*5,1] = IS_cutoff
			performance_test1[i+j*5,2] = Overlap_cutoff   
			performance_test2[i+j*5,1] = IS_cutoff                                                                                     
			performance_test2[i+j*5,2] = Overlap_cutoff                                                                                     
			p1 = Cutoff_Evaluation3(i,j)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
			performance_test1[i+j*5,3:5] = p1[2:4]
			p2 = Cutoff_Evaluation4(i,j)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
			performance_test2[i+j*5,3:5] = p2[2:4]
		}
	}
	
	par(mfrow=c(2,2))
	plot(performance_test1[,1],performance_test1[,2])
	plot(performance_test1[,1],performance_test1[,3])
	plot(performance_test1[,1],performance_test1[,4])
	plot(performance_test1[,1],performance_test1[,5])
	par(mfrow=c(2,2))
	plot(performance_test2[,1],performance_test2[,2])
	plot(performance_test2[,1],performance_test2[,3])
	plot(performance_test2[,1],performance_test2[,4])
	plot(performance_test2[,1],performance_test2[,5])
	
	single_dropout = matrix(performance_test1[6:30,3],nrow=5,ncol=5)
	double_dropout = matrix(performance_test2[6:30,3],nrow=5,ncol=5)
	rownames( single_dropout ) = 1:5
	colnames( single_dropout ) = 1:5 
	
	rownames( double_dropout ) = 1:5
	colnames( double_dropout ) = 1:5
		
	simple_plotMatrix(single_dropout,0,filename="single_dropout.cutoff.pdf")   
	simple_plotMatrix(double_dropout,0,filename="double_dropout.cutoff.pdf")      
	##================ Fig 2 Last (AUC cutoff)          
	 	
	
	refine_bm1_W = RefineMatrix(bm1_W) 
	refine_bm1_A = RefineMatrix(bm1_A)  
	refine_bm1_RW = RefineMatrix(bm1_RW) 
	refine_bm1_RA = RefineMatrix(bm1_RA) 
	
	refine_bm2_W = RefineMatrix(bm2_W)
	refine_bm2_A = RefineMatrix(bm2_A)
	refine_bm2_Q = RefineMatrix(bm2_Q)
	
	refine_bm3_W = RefineMatrix(bm3_W)
	refine_bm3_A = RefineMatrix(bm3_A)
	refine_bm3_Q = RefineMatrix(bm3_Q)
	refine_bm3_QL = RefineMatrix(bm3_QL)
	
	refine_bm4_SW = RefineMatrix(bm4_SW)
	refine_bm4_S4A = RefineMatrix(bm4_S4A)
	refine_bm4_S8A = RefineMatrix(bm4_S8A)
	refine_bm4_SQ = RefineMatrix(bm4_SQ)
	
	refine_bm5_SW = RefineMatrix(bm5_SW)	   
	                                          
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
	               
	
	# For homo-dimer
	View( cbind(colSums(bm3_W), colSums(bm4_SW), diag(bm3_W), diag(bm4_SW), diag(bm3_A), diag( bm3_Q), diag(bm4_S8A), diag(bm4_S4A), diag(bm4_SQ)) )
	
	# PPI database
	HIPPIE <- ReadMatrix( "~/Dropbox (CRG)/Code/jslib/Mireia/HIPPIE/HIPPIE.R75.txt" )
	pheatmap((HIPPIE), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=6)

	BioGrid <- ReadMatrix( "~/Dropbox (CRG)/Code/jslib/Mireia/BioGrid/BioGrid.R75.txt" )
	pheatmap((BioGrid), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=7)   
	      
	KnownPPI = HIPPIE | BioGrid
	
	#heatmap2_plotLog2Matrix(jb3_A,1,filename="jb3_A.heatmap.pdf")        
	#heatmap2_plotLog2Matrix(jb3_W,1,filename="jb3_W.heatmap.pdf") 
	refine_KnownPPI = RefineMatrix(KnownPPI)
	simple_plotMatrix(refine_KnownPPI,0,7)      
	
	
	barplot(S1_BWD_R1.cnt[,2])
	barplot(S2_BA2D_R1.cnt[,2])
	barplot(S3_BQD_R1.cnt[,2])
	barplot(S4_AWD_R2.cnt[,2])       
	# R75_11,AP2B1,adaptor-related protein complex 2, beta 1 subunit (BA2D,BQD)
	# R75_30,CTBP1,C-terminal binding protein 1                      (BA2D,)
	# R75_52,RBBP8,retinoblastoma binding protein 8                  (BQD)
	# R75_58,FHL3,four and a half LIM domains 3                      (BA2D,)
	# R75_70,TEX11,testis expressed 11                               (BA2D,BQD)
	
	# PIAS1                                                          (BA2D,)
	# CNOT7															 (BA2D,) 
	
	bqd = (S3_BQD_R1.cnt[,2]/sum(S3_BQD_R1.cnt[,2]))/(S1_BWD_R1.cnt[,2]/sum(S1_BWD_R1.cnt[,2])) 
	bad = (S2_BA2D_R1.cnt[,2]/sum(S2_BA2D_R1.cnt[,2]))/(S1_BWD_R1.cnt[,2]/sum(S1_BWD_R1.cnt[,2])) 
	names(bad) = S3_BQD_R1.cnt[,1]
	names(bqd) = S3_BQD_R1.cnt[,1]        
	                     
	bad[ is.na(bad) ] = 0                  
	bqd[ is.na(bqd) ] = 0
	simple_plotLog2Matrix(rbind(bad,bqd),1,filename="AutoActivity.Screen.pdf")
	
	sbqd = (sS3_BQD_R1.cnt[,2]/sum(sS3_BQD_R1.cnt[,2]))/(sS1_BWD_R1.cnt[,2]/sum(sS1_BWD_R1.cnt[,2])) 
	sbad = (sS2_BA2D_R1.cnt[,2]/sum(sS2_BA2D_R1.cnt[,2]))/(sS1_BWD_R1.cnt[,2]/sum(sS1_BWD_R1.cnt[,2])) 
	names(sbad) = sS3_BQD_R1.cnt[,1]
	names(sbqd) = sS3_BQD_R1.cnt[,1]       
	
	
	        
	#===========================
	### For Figure 3
	#===========================
	#is1_A = InteractionScores(bm1_W,bm1_A,1.0)        
	#PPI_PerformanceCheck2(is1_A, HIPPIE)      
	
	#refined_HIPPIE = RefineMatrix(HIPPIE)
	#is2_refine_bm2_A = InteractionScores(refine_bm2_W,NormalMixture2(refine_bm2_A),1,9 )		  
	is2_bm2_A = InteractionScores(bm2_W,NormalMixture2(bm2_A),1,9 )		  
	PPI_PerformanceCheck2(is2_bm2_A, HIPPIE)         
	index = order(is2_bm2_A,decreasing = TRUE)
	HIPPIE[index]    
	
	pdf( "is2_bm2_A.1-400.pdf" )
	par(mfrow=c(3,1))
	barplot( is2_bm2_A[index][1:400] )  
	barplot( HIPPIE[index][1:400] )
	barplot( BioGrid[index][1:400] )
	dev.off()
   
    pdf( "is2_bm2_A.1-30.pdf" )    
    par(mfrow=c(3,1))   
	barplot( is2_bm2_A[index][1:30] )  
	barplot( HIPPIE[index][1:30] )
	barplot( BioGrid[index][1:30] )
	dev.off()
	
	GetPairInfo(is2_bm2_A,index[1:30])
	
	sum(is2_bm2_A >5.46)
	simple_plotMatrix(is2_bm2_A>5.96,0,9)
	sum(is2_bm2_A>5.96)
		
}         
 

GetPairInfo <-function(m, index){    
	# row - col
	for(i in index){
		output = sprintf( "%d - %f - %s - %s",i, m[i], rownames(m)[(i-1)%%78+1],rownames(m)[((i-1)/78)+1])
		print(output)
	}
}

## 75X library - Sebastian - local - X75                                
if (FALSE){                                               
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )             
	
	sa2_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_51-WCCRD_S1_sorted_flag113tidyformat.xlsx.txt",0.5) 
	sa2_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_51-A1CCRD_S2_sorted_flag113tidyformat.xlsx.txt",0.5) 
	sa2_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18691_S3_sorted_flag113tidyformat.xlsx.txt",0.5)
		 
	sa3_W = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18689_S1_sorted_flag113tidyformat.xlsx.txt",0.5) 
	sa3_A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18690_S2_sorted_flag113tidyformat.xlsx.txt",0.5) 
	sa3_Q = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_18691_S3_sorted_flag113tidyformat.xlsx.txt",0.5)

	sa4_SW = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21387-64-SWD_S1_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
	sa4_S8A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21389-64-SA8D_S3_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
	sa4_S4A = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21388-64-SA4D_S2_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
	sa4_SQ = drawMatrix("/Users/destine/Dropbox_CRG/MatrixPPi/data/sebastian/local/X75/aligned\ vs\ final\ 150bps/Ref_X75_seq_21390-64-SQD_S4_sorted_flag113tidyformat.xlsx.txt", 0.5, 7 )
	                                                                                                                                                   
	View( cbind(colSums(sa2_W), colSums(sa4_SW), diag(sa2_W), diag(sa4_SW), diag(sa2_A), diag( sa2_Q), diag(sa4_S8A), diag(sa4_S4A), diag(sa4_SQ)) )     
	
	View( cbind(colSums(sa3_W), colSums(sa4_SW), diag(sa3_W), diag(sa4_SW), diag(sa3_A), diag( sa3_Q), diag(sa4_S8A), diag(sa4_S4A), diag(sa4_SQ)) )
	View( cbind(rowSums(sa3_W), rowSums(sa4_SW) ) )
	
	ib3SE_A = InteractionScores(sa3_W,sa3_A,1,4)   
	ib3SE_Q = InteractionScores(sa3_W,sa3_Q,1,4)     
	
	ib4SE_8A = InteractionScores(sa4_SW,sa4_S8A,1,4)          
	ib4SE_4A = InteractionScores(sa4_SW,sa4_S4A,1,4)          
	ib4SE_Q = InteractionScores(sa4_SW,sa4_SQ,1,4)         
	
	PlotHomoDimerTest(sa2_W,sa2_A,sa2_Q,sa3_W,sa3_A,sa3_Q,sa4_SW,sa4_S4A,sa4_SQ) 
	                                                   
	par(mfrow=c(3,4))
	HomoDimerNormalization(sa2_W, sa2_A, roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-4Q")     

	par(mfrow=c(3,4))
	HomoDimerNormalization(sa2_W, sa2_A, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")	   
	
	par(mfrow=c(3,4))
	HomoDimerNormalization(sa2_W, sa2_A, method="diag", roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, method="diag",roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, method="diag",roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, method="diag",roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, method="diag",roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, method="diag",roc_plot1=TRUE, roc_plot2=FALSE, box_plot1=TRUE, box_plot2=FALSE)      
	mtext("Homo dimer exp-4Q")	 
                          
	par(mfrow=c(3,8))
	HomoDimerNormalization(sa2_W, sa2_A, method="colSums", roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, method="colSums",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")
	
	par(mfrow=c(3,8))
	HomoDimerNormalization(sa2_W, sa2_A, method="diag", roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerNormalization(sa2_W, sa2_Q, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerNormalization(sa3_W, sa3_A, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerNormalization(sa3_W, sa3_Q, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerNormalization(sa4_SW, sa4_S4A, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerNormalization(sa4_SW, sa4_SQ, method="diag",roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")
	
		
	par(mfrow=c(3,8))
	HomoDimerTest(sa2_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerTest(sa2_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerTest(sa3_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerTest(sa3_W, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerTest(sa4_SW, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerTest(sa4_SW, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")	    
	
	par(mfrow=c(3,8))
	HomoDimerTest(sa2_A, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerTest(sa2_Q, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)   
	mtext("Homo dimer exp-2Q")
	HomoDimerTest(sa3_A, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)       
	mtext("Homo dimer exp-3A")
	HomoDimerTest(sa3_Q, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)   
	mtext("Homo dimer exp-3Q")
	HomoDimerTest(sa4_S4A, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-4A4")
	HomoDimerTest(sa4_SQ, roc_plot1=TRUE, roc_plot2=TRUE, box_plot1=TRUE, box_plot2=TRUE)    
	mtext("Homo dimer exp-4Q")
	
	par(mfrow=c(3,4))
	HomoDimerTest(sa2_W, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerTest(sa2_W, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerTest(sa3_W, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerTest(sa3_W, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerTest(sa4_SW, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerTest(sa4_SW, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")	    
	
	par(mfrow=c(3,4))
	HomoDimerTest(sa2_A, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)    
	mtext("Homo dimer exp-2A")  
	HomoDimerTest(sa2_Q, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-2Q")
	HomoDimerTest(sa3_A, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3A")
	HomoDimerTest(sa3_Q, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-3Q")
	HomoDimerTest(sa4_S4A, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4A4")
	HomoDimerTest(sa4_SQ, roc_plot1=FALSE, roc_plot2=TRUE, box_plot1=FALSE, box_plot2=TRUE)      
	mtext("Homo dimer exp-4Q")
}

LoadP170Data <- function(){
	bP170_1_W <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S53.ppi.txt.new",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3250858
	bP170_1_Q <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S54.ppi.txt.new",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3553164
    #gm_bP170_1_Q <<- GammaMixture(bP170_1_Q)      
	#ibP170_1_Q <<- InteractionScores(bP170_1_W,gm_bP170_1_Q,1,4)
	#nibP170_1_Q <<- NewInteractionScores(bP170_1_W,gm_bP170_1_Q,1,4)     
	                                              
	bP170_2_W <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_W.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "sum = 2629509"
	bP170_2_Q <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_Q.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in        "sum = 2507307"
	#gm_bP170_2_Q <<- GammaMixture(bP170_2_Q)                     
	#ibP170_2_Q <<- InteractionScores(bP170_2_W,bP170_2_Q,1,4) 
	#ibP170_2_Q <<- InteractionScores(bP170_2_W,gm_bP170_2_Q,1,4)
	#nibP170_2_Q <<- NewInteractionScores(bP170_2_W,gm_bP170_2_Q,1,4)  
	   
	#bP170_3_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
	bP170_3_W <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
	bP170_3_Q <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
	#gm_bP170_3_Q <<- GammaMixture(bP170_3_Q)           
	#ibP170_3_Q <<- InteractionScores(bP170_3_W,bP170_3_Q,1,4) 
	#ibP170_3_Q <<- InteractionScores(bP170_3_W,gm_bP170_3_Q,1,4)
	#nibP170_3_Q <<- NewInteractionScores(bP170_3_W,gm_bP170_3_Q,1,4)       
	
	bP170_4_SW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
	bP170_4_S4A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SA4D.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/4; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
	bP170_4_S8A <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SA8D.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/8; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
	bP170_4_SQ <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/-H/-A; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
	#gm_bP170_4_SQ <<- GammaMixture(bP170_4_SQ)      
	#ibP170_4_SQ <<- InteractionScores(bP170_4_SW,gm_bP170_4_SQ,1,4)
	#nibP170_4_SQ <<- NewInteractionScores(bP170_4_SW,gm_bP170_4_SQ,1,4)  
	#ibP170_4_S4A <<- InteractionScores(bP170_4_SW,bP170_4_S4A,1,4) 

	bP170_5_SW <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W;             # ==> 1836282 out of ? (48.4%)
	bP170_5_SA8 <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SAD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/8;  # ==>  out of  (43.7%)
	bP170_5_SQ <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/-H/-A; R75_01 # ==>  out of  (47.4%)
	#gm_bP170_5_SQ <<- GammaMixture(bP170_5_SQ)          
	#is_P170_5_SQ <<- InteractionScores(bP170_5_SW,bP170_5_SQ,1,4)            
	#ibP170_5_SQ <<- InteractionScores(bP170_5_SW,gm_bP170_5_SQ,1,4)
	#nibP170_5_SQ <<- NewInteractionScores(bP170_5_SW,gm_bP170_5_SQ,1,4)  
	#pisP170_5_SQ  <<- PairInteractionScores(bP170_5_SW,bP170_5_SQ,1,4)      
	
	#bP170_6_SW1_OLD <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S1_WD.ppi.txt",0.5,4)
	#bP170_6_SW2_OLD <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S2_WD.ppi.txt",0.5,4)      
	#bP170_6_SA1_OLD <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S3_AD.ppi.txt",0.5,4)    
	#bP170_6_SA2_OLD <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S4_AD.ppi.txt",0.5,4)    
	#ibP170_6_SA1_OLD <<- InteractionScores(bP170_6_SW1_OLD,bP170_6_SA1_OLD,1,4)      
	#ibP170_6_SA2_OLD <<- InteractionScores(bP170_6_SW2_OLD,bP170_6_SA2_OLD,1,4)      
	
	bP170_6_SW1 <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S1_WD.ppi.txt.new",0.5,4)
	bP170_6_SW2 <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S2_WD.ppi.txt.new",0.5,4)      
	bP170_6_SA1 <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S3_AD.ppi.txt.new",0.5,4)    
	bP170_6_SA2 <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S4_AD.ppi.txt.new",0.5,4)    
	#ibP170_6_SA1 <<- InteractionScores(bP170_6_SW1,bP170_6_SA1,1,4)      
	#ibP170_6_SA2 <<- InteractionScores(bP170_6_SW2,bP170_6_SA2,1,4)	 
	
	#pisP170_6_SA1 <<- PairInteractionScores(bP170_6_SW1,bP170_6_SA1,1,4)      
	#pisP170_6_SA2 <<- PairInteractionScores(bP170_6_SW2,bP170_6_SA2,1,4)	 

    # ============================
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
}

## 163X library  - Final Used                              
if (FALSE){ 
	LoadP170Data()

	#===============================
	# Save
	
	melted_is2_bP170_1_Q = melt( is2_bP170_1_Q )  
	melted_is2_bP170_2_Q = melt( is2_bP170_2_Q )  
	melted_is2_bP170_3_Q = melt( is2_bP170_3_Q )  
	melted_is2_bP170_4_S4A= melt( is2_bP170_4_S4A )  
	melted_is2_bP170_4_S8A = melt( is2_bP170_4_S8A )  
	melted_is2_bP170_4_SQ = melt( is2_bP170_4_SQ )  
	melted_is2_bP170_5_SA8 = melt( is2_bP170_5_SA8 )  
	melted_is2_bP170_5_SQ = melt( is2_bP170_5_SQ )  
	melted_is2_bP170_6_SA1 = melt( is2_bP170_6_SA1 )  
	melted_is2_bP170_6_SA2 = melt( is2_bP170_6_SA2 )
	
	write.table(melted_is2_bP170_1_Q,file="/Users/jyang/melted_is2_bP170_1_Q.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_2_Q,file="/Users/jyang/melted_is2_bP170_2_Q.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_3_Q,file="/Users/jyang/melted_is2_bP170_3_Q.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_4_S4A,file="/Users/jyang/melted_is2_bP170_4_S4A.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_4_S8A,file="/Users/jyang/melted_is2_bP170_4_S8A.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_4_SQ,file="/Users/jyang/melted_is2_bP170_4_SQ.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_5_SA8,file="/Users/jyang/melted_is2_bP170_5_SA8.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_5_SQ,file="/Users/jyang/melted_is2_bP170_5_SQ.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_6_SA1,file="/Users/jyang/melted_is2_bP170_6_SA1.IS2.txt",sep="\t")	 
	write.table(melted_is2_bP170_6_SA2,file="/Users/jyang/melted_is2_bP170_6_SA2.IS2.txt",sep="\t")	 
	
	melted_is2_bP170_ALL = data.frame( P1 = melted_is2_bP170_1_Q$Var1, P2 = melted_is2_bP170_1_Q$Var2, Q1 = melted_is2_bP170_1_Q$value, Q2 = melted_is2_bP170_2_Q$value, Q3 = melted_is2_bP170_3_Q$value, A44 = melted_is2_bP170_4_S4A$value, A48 = melted_is2_bP170_4_S8A$value, Q4 = melted_is2_bP170_4_SQ$value, A58 = melted_is2_bP170_5_SA8$value, Q5 = melted_is2_bP170_5_SQ$value, A61 = melted_is2_bP170_6_SA1$value, A62 = melted_is2_bP170_6_SA2$value )
	write.table(melted_is2_bP170_ALL,file="/Users/jyang/melted_is2_bP170_ALL.IS2.txt",sep="\t")	 
	
	
	
	is2_bP170_overlap = (is2_bP170_2_Q > 3) + (is2_bP170_3_Q > 3) + (is2_bP170_4_S4A > 3) + (is2_bP170_4_S8A > 3) + (is2_bP170_4_SQ > 3) + (is2_bP170_5_SA8 > 3) + (is2_bP170_5_SQ > 3) + (is2_bP170_6_SA1 > 3) + (is2_bP170_6_SA2 > 3)    
	simple_plotMatrix(is2_bP170_overlap,0,4)
	
	is2_data = data.frame( is2_bP170_2_Q=as.vector(is2_bP170_2_Q), is2_bP170_3_Q=as.vector(is2_bP170_3_Q),is2_bP170_4_S4A=as.vector(is2_bP170_4_S4A), is2_bP170_4_S8A=as.vector(is2_bP170_4_S8A), is2_bP170_4_SQ=as.vector(is2_bP170_4_SQ), is2_bP170_5_SA8=as.vector(is2_bP170_5_SA8), is2_bP170_5_SQ=as.vector(is2_bP170_5_SQ), is2_bP170_6_SA1=as.vector(is2_bP170_6_SA1), is2_bP170_6_SA2=as.vector(is2_bP170_6_SA2) )	
	pdf("bP170.is2.correlation.pdf",width = 600, height = 600); chart.Correlation(is2_data,method="pearson"); dev.off()   
	
	# Save file
	melted_is2_bP170_overlap = melt(is2_bP170_overlap)
	write.table(melted_is2_bP170_overlap,file="/Users/jyang/melted_is2_bP170_overlap.IS2.txt",sep="\t")	  
	
	
	## Auto activators
	cbind( 1:173, rowSums(is2_bP170_overlap) )
	
	simple_plotLog10Matrix(RefineMatrixP170(bP170_5_SQ),0.5,4)
	                  
	PPI_PerformanceCheck_P170( bP170_5_SQ, HIPPIE_P170 )   
	PPI_PerformanceCheck_P170( ibP170_5_SQ, HIPPIE_P170 )   
	PPI_PerformanceCheck_P170( is_P170_5_SQ, HIPPIE_P170 )   
	PPI_PerformanceCheck_P170( nibP170_5_SQ, HIPPIE_P170 )   
	PPI_PerformanceCheck_P170( pisP170_5_SQ, HIPPIE_P170 )       
	PPI_PerformanceCheck_P170( is2_bP170_5_SQ, HIPPIE_P170 )

	PPI_PerformanceCheck_P170( is2_bP170_5_SQ, BioGrid_P170 )

	
	write.table(melt(HIPPIE_P170),file="/Users/jyang/melted_HIPPIE_P170.txt",sep="\t")             
	write.table(melt(RefineMatrixP170(is2_bP170_2_Q)),file="/Users/jyang/melted_is2_bP170_2_Q.txt",sep="\t")	   
	                  
	
	P170_PPI_DATA = HIPPIE_P170
	#P170_PPI_DATA = BioGrid_P170_ALL
	
	## Fig 4 - Performance Check
	Cutoff_Evaluation1_P170 <- function( IS_cut = 1, overlap_cut = 1, P170_PPI_DATA = BioGrid_P170_ALL ){                             
	    p = PPI_PerformanceCheck_P170((is2_bP170_4_S4A>=IS_cut)+(is2_bP170_4_S8A>=IS_cut)+(is2_bP170_5_SA8>=IS_cut)+(is2_bP170_6_SA1>=IS_cut)+(is2_bP170_6_SA2>=IS_cut),P170_PPI_DATA)            
		#p = PPI_PerformanceCheck_P170((is2_bP170_4_S4A>=IS_cut)+(is2_bP170_6_SA1>=IS_cut)+(is2_bP170_6_SA2>=IS_cut),P170_PPI_DATA)            
		return(p)
	}   
	
	Cutoff_Evaluation2_P170 <- function( IS_cut = 1, overlap_cut = 1, P170_PPI_DATA = BioGrid_P170_ALL  ){        
	    p = PPI_PerformanceCheck_P170((is2_bP170_2_Q>=IS_cut)+(is2_bP170_3_Q>=IS_cut)+(is2_bP170_4_SQ>=IS_cut)+(is2_bP170_5_SQ>=IS_cut), P170_PPI_DATA)   
		#p = PPI_PerformanceCheck_P170( (is2_bP170_4_SQ>=IS_cut)+(is2_bP170_5_SQ>=IS_cut), P170_PPI_DATA)            
		return(p)
	}         

	Cutoff_Evaluation3_P170 <- function( IS_cut = 1, overlap_cut = 1, P170_PPI_DATA = BioGrid_P170_ALL  ){                             
	    p = PPI_PerformanceCheck_P170( (is2_bP170_4_S4A>=IS_cut)+(is2_bP170_4_S8A>=IS_cut)+(is2_bP170_5_SA8>=IS_cut)+(is2_bP170_6_SA1>=IS_cut)+(is2_bP170_6_SA2>=IS_cut)+(is2_bP170_2_Q>=IS_cut)+(is2_bP170_3_Q>=IS_cut)+(is2_bP170_4_SQ>=IS_cut)+(is2_bP170_5_SQ>=IS_cut), P170_PPI_DATA)            
		return(p)
	}   
	
	#Cutoff_Evaluation4_P170 <- function( IS_cut = 1, overlap_cut = 1, P170_PPI_DATA = BioGrid_P170_ALL  ){        
	#    p = PPI_PerformanceCheck_P170( ( (is2_bP170_2_Q>IS_cut)+(is2_bP170_3_Q>IS_cut)+(is2_bP170_4_SQ>IS_cut)+(is2_bP170_5_SQ>IS_cut) ) >= overlap_cnt, P170_PPI_DATA)   
	#	#p = PPI_PerformanceCheck_P170( (is2_bP170_4_SQ>IS_cut)+(is2_bP170_5_SQ>IS_cut), P170_PPI_DATA)            
	#	return(p)
	#}         
	                                                                                                                          
	performance_test1 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for single drop out
	for (i in 1:50){       
		IS_cutoff = i / 10
		performance_test1[i,1] = IS_cutoff   
		p1 = Cutoff_Evaluation1_P170(i/10, 1, P170_PPI_DATA)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
		performance_test1[i,2:5] = p1[1:4]
	}

	performance_test2 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for double drop out
	for (i in 1:50){       
		IS_cutoff = i / 10
		performance_test2[i,1] = IS_cutoff                                                                                     
		p2 = Cutoff_Evaluation2_P170(i/10, 1, P170_PPI_DATA)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
		performance_test2[i,2:5] = p2[1:4]
	}

	performance_test3 = cbind(1:50,1:50,1:50,1:50,1:50) * 0 # for double drop out
	for (i in 1:50){       
		IS_cutoff = i / 10
		performance_test3[i,1] = IS_cutoff                                                                                     
		p3 = Cutoff_Evaluation3_P170(i/10, 1, P170_PPI_DATA)   # auc, -log10 p_value, Positive PPI mean, Negative PPI mean                  
		performance_test3[i,2:5] = p3[1:4]
	}
	
	par(mfrow=c(2,2))
	plot(performance_test1[,1],performance_test1[,2])
	plot(performance_test1[,1],performance_test1[,3])
	plot(performance_test1[,1],performance_test1[,4])
	plot(performance_test1[,1],performance_test1[,5])
	par(mfrow=c(2,2))
	plot(performance_test2[,1],performance_test2[,2])
	plot(performance_test2[,1],performance_test2[,3])
	plot(performance_test2[,1],performance_test2[,4])
	plot(performance_test2[,1],performance_test2[,5])
	par(mfrow=c(2,2))
	plot(performance_test3[,1],performance_test3[,2])
	plot(performance_test3[,1],performance_test3[,3])
	plot(performance_test3[,1],performance_test3[,4])
	plot(performance_test3[,1],performance_test3[,5])        
	
	IS_cut = 2.0
	overall_P170 = (is2_bP170_4_S4A>=IS_cut)+(is2_bP170_4_S8A>=IS_cut)+(is2_bP170_5_SA8>=IS_cut)+(is2_bP170_6_SA1>=IS_cut)+(is2_bP170_6_SA2>=IS_cut)+(is2_bP170_2_Q>=IS_cut)+(is2_bP170_3_Q>=IS_cut)+(is2_bP170_4_SQ>=IS_cut)+(is2_bP170_5_SQ>=IS_cut)
	simple_plotMatrix(overall_P170>=4,0,4)      
	            
	## for gene selection
	IS_cut = 1.6
	overall_P170_over1_6 = (is2_bP170_4_S4A>=IS_cut)+(is2_bP170_4_S8A>=IS_cut)+(is2_bP170_5_SA8>=IS_cut)+(is2_bP170_6_SA1>=IS_cut)+(is2_bP170_6_SA2>=IS_cut)+(is2_bP170_2_Q>=IS_cut)+(is2_bP170_3_Q>=IS_cut)+(is2_bP170_4_SQ>=IS_cut)+(is2_bP170_5_SQ>=IS_cut)
	simple_plotMatrix(overall_P170_over1_6>=4,0,4)    
	melted_overall_P170_over1_6 = melt( overall_P170_over1_6 )
	write.table(melted_overall_P170_over1_6,file="/Users/jyang/melted_overall_P170_over1_6.txt",sep="\t")   
	
	                                                  
	overall_P170_melt = melt( overall_P170 )                  
	BioGrid_P170_melt = melt(BioGrid_P170)       
	HIPPIE_P170_melt = melt(HIPPIE_P170)
	overall_P170_df = data.frame(P1 = overall_P170_melt[,1],P2 = overall_P170_melt[,2], overlap = overall_P170_melt[,3], BioGridAll = BioGrid_P170_melt[,3], HIPPIE = HIPPIE_P170_melt[,3] )
	view(overall_P170_df)
	BioGrid_P170_df = data.frame(melt(BioGrid_P170))
} 


## 163X library                                
if (FALSE){ 
	bP170_1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S53.ppi.txt.new.2",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3250858
	bP170_1_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-12_MiSeq/Blastn/S54.ppi.txt.new.2",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3553164
    gm_bP170_1_Q = GammaMixture(bP170_1_Q)      
	ibP170_1_Q = InteractionScores(bP170_1_W,gm_bP170_1_Q,1,4)
	nibP170_1_Q = NewInteractionScores(bP170_1_W,gm_bP170_1_Q,1,4)       
	                                              
	bP170_2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_W.ppi.txt.new.2",0.5,4) # BD P170-4 x AD P170-4; -W; R75_01 ~ R75_07 spike in                      "sum = 2629509"
	bP170_2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-07-03_MiSeq/Blastn/60_Q.ppi.txt.new.2",0.5,4) # BD P170-4 x AD P170-4; -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in        "sum = 2507307"
	gm_bP170_2_Q = GammaMixture(bP170_2_Q)      
	ibP170_2_Q = InteractionScores(bP170_2_W,gm_bP170_2_Q,1,4)
	nibP170_2_Q = NewInteractionScores(bP170_2_W,gm_bP170_2_Q,1,4)       
	   
	#bP170_3_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
	bP170_3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
	bP170_3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-15_MiSeq/Blastn/S61_PQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Plate); -W/-H/-A/Aba 1/2; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
	gm_bP170_3_Q = GammaMixture(bP170_3_Q)      
	ibP170_3_Q = InteractionScores(bP170_3_W,gm_bP170_3_Q,1,4)
	nibP170_3_Q = NewInteractionScores(bP170_3_W,gm_bP170_3_Q,1,4)       
	
	bP170_4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W; R75_01 ~ R75_07 spike in                "# ==> 2909222 out of 6012622 (48.4%)
	bP170_4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SA4D.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/4; R75_01 ~ R75_07 spike in      		  "# ==> 1956131 out of 4473498 (43.7%)
	bP170_4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SA8D.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/8; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
	bP170_4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-28_MiSeq/Blastn/S68_SQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/-H/-A; R75_01 ~ R75_07 spike in  "# ==> 1684777 out of 3553019 (47.4%)
	gm_bP170_4_SQ = GammaMixture(bP170_4_SQ)      
	ibP170_4_SQ = InteractionScores(bP170_4_SW,gm_bP170_4_SQ,1,4)
	nibP170_4_SQ = NewInteractionScores(bP170_4_SW,gm_bP170_4_SQ,1,4)  
	ibP170_4_S4A = InteractionScores(bP170_4_SW,bP170_4_S4A,1,4)                                                                       

	bP170_5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SWD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W;             # ==> 1836282 out of ? (48.4%)
	bP170_5_SA8 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SAD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/Aba 1/8;  # ==>  out of  (43.7%)
	bP170_5_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-30_MiSeq/Blastn/S65_SQD.ppi.txt.new",0.5,4) # BD P170-4 x AD P170-4 (Liquid); -W/-H/-A; R75_01 # ==>  out of  (47.4%)
	gm_bP170_5_SQ = GammaMixture(bP170_5_SQ)      
	ibP170_5_SQ = InteractionScores(bP170_5_SW,gm_bP170_5_SQ,1,4)
	nibP170_5_SQ = NewInteractionScores(bP170_5_SW,gm_bP170_5_SQ,1,4)   
	
	bP170_6_SW1_OLD = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S1_WD.ppi.txt",0.5,4)
	bP170_6_SW2_OLD = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S2_WD.ppi.txt",0.5,4)      
	bP170_6_SA1_OLD = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S3_AD.ppi.txt",0.5,4)    
	bP170_6_SA2_OLD = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S4_AD.ppi.txt",0.5,4)    
	ibP170_6_SA1_OLD = InteractionScores(bP170_6_SW1_OLD,bP170_6_SA1_OLD,1,4)      
	ibP170_6_SA2_OLD = InteractionScores(bP170_6_SW2_OLD,bP170_6_SA2_OLD,1,4)      
	
	bP170_6_SW1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S1_WD.ppi.txt.new.2",0.5,4)
	bP170_6_SW2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S2_WD.ppi.txt.new.2",0.5,4)      
	bP170_6_SA1 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S3_AD.ppi.txt.new.2",0.5,4)    
	bP170_6_SA2 = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-16_MiSeq/Blastn/S4_AD.ppi.txt.new.2",0.5,4)    
	ibP170_6_SA1 = InteractionScores(bP170_6_SW1,bP170_6_SA1,1,4)      
	ibP170_6_SA2 = InteractionScores(bP170_6_SW2,bP170_6_SA2,1,4)	 
	
	pisP170_6_SA1 = PairInteractionScores(bP170_6_SW1,bP170_6_SA1,1,4)      
	pisP170_6_SA2 = PairInteractionScores(bP170_6_SW2,bP170_6_SA2,1,4)	 
	
}
  

## A463-MGj69.RBP-MAP library                                
if (FALSE){     
	# 452*452
	bA463_1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-17_MiSeq/S1_W.ppi.txt",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 3106500
	bA463_1_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-17_MiSeq/S2_Q.ppi.txt",0.5,4)    # BD P170-4 x AD P170-4 (Plate); ==> 2158636  
	
	smA463_1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-17_MiSeq/Sebastian/S1_W.ppi.txt",0.5,4)    # BD 461 x AD 461 (Plate); ==> 2899475
	smA463_1_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-17_MiSeq/Sebastian/S2_Q.ppi.txt",0.5,4)    # BD 461 x AD 461 (Plate); ==> 1964393
	
    gm_bA463_1_Q = GammaMixture(bA463_1_Q) 
    nm_bA463_1_Q = NormalMixture2(bA463_1_Q)   
    is_bA463_1_Q = InteractionScores(bA463_1_W,bA463_1_Q,1,4)  
    is_gm_bA463_1_Q = InteractionScores(bA463_1_W,gm_bA463_1_Q,1,4)  
	is2_bA463_1_Q = InteractionScores(bA463_1_W,nm_bA463_1_Q,1,4) 
	simple_plotMatrix(ibA463_1_Q>3,0,4)
}
         

## EGFR library                                
if (FALSE){     
	bEGFR1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S1_WD.ppi.txt.new.2", 0.5, 4 )
	bEGFR1_4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S2_WD.ppi.txt.new.2", 0.5, 4 )     
	bEGFR1_8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S3_A8D.ppi.txt.new.2", 0.5, 4 )     
	bEGFR1_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-10-09_MiSeq/Blastn/S4_QD.ppi.txt.new.2", 0.5, 4 )     
	
	is_EGFR1_4A = InteractionScores(bEGFR1_W, bEGFR1_4A,1,4)
	is_EGFR1_8A = InteractionScores(bEGFR1_W, bEGFR1_8A,1,4)
	is_EGFR1_Q = InteractionScores(bEGFR1_W, bEGFR1_Q,1,4)  
	
	is2_EGFR1_4A = InteractionScores(bEGFR1_W, NormalMixture2(bEGFR1_4A),1,4)
	is2_EGFR1_8A = InteractionScores(bEGFR1_W, NormalMixture2(bEGFR1_8A),1,4)
	is2_EGFR1_Q = InteractionScores(bEGFR1_W, NormalMixture2(bEGFR1_Q),1,4)  
	
	
	CorMatrix(is_EGFR1_Q,is_EGFR1_4A)   
	cis_EGFR1_4A = ColumnWiseInteractionScores(bEGFR1_W, bEGFR1_4A,1,4,rowNoralization = TRUE)
	cis_EGFR1_8A = ColumnWiseInteractionScores(bEGFR1_W, bEGFR1_8A,1,4,rowNoralization = TRUE)
	cis_EGFR1_Q = ColumnWiseInteractionScores(bEGFR1_W, bEGFR1_Q,1,4,rowNoralization = TRUE)     
	
	Rowwise = cor(t(bEGFR1_W))[ upper.tri( cor(t(bEGFR1_W)) ) ]
	Columnwise = cor((bEGFR1_W))[ upper.tri( cor((bEGFR1_W)) ) ]  
	pdf("EGFR.W.col_row_correlation.pdf")
	boxplot(Rowwise,Columnwise, names=c("All Row-wise pairs","All Column-wise pairs"), ylab="Pearson Correlation (r)")
	dev.off()   
	
	withoutAA = setdiff(1:182,c(6,  15,  23,  57,  62,  72,  73,  97, 111, 115, 123, 124, 127, 128, 170))
	Rowwise = cor(t(bEGFR1_4A[withoutAA,]))[ upper.tri( cor(t(bEGFR1_4A[withoutAA,])) ) ]                                                                                                                   
	Columnwise = cor((bEGFR1_4A[withoutAA,]))[ upper.tri( cor((bEGFR1_4A[withoutAA,])) ) ]  
	pdf("EGFR.4A.col_row_correlation.pdf")
	boxplot(Rowwise,Columnwise, names=c("All Row-wise pairs","All Column-wise pairs"), ylab="Pearson Correlation (r)")    
	dev.off()     
	
	withoutAA = setdiff(1:182,c(6,  15,  23,  57,  62,  72,  73,  97, 111, 115, 123, 124, 127, 128, 170))
	Rowwise = cor(t(bEGFR1_8A[withoutAA,]))[ upper.tri( cor(t(bEGFR1_8A[withoutAA,])) ) ]                                                                                                                   
	Columnwise = cor((bEGFR1_8A[withoutAA,]))[ upper.tri( cor((bEGFR1_8A[withoutAA,])) ) ]  
	pdf("EGFR.8A.col_row_correlation.pdf")
	boxplot(Rowwise,Columnwise, names=c("All Row-wise pairs","All Column-wise pairs"), ylab="Pearson Correlation (r)")    
	dev.off()
	
	withoutAA = setdiff(1:182,c(6,  15,  23,  57,  62,  72,  73,  97, 111, 115, 123, 124, 127, 128, 170))
	Rowwise = cor(t(bEGFR1_Q[withoutAA,]))[ upper.tri( cor(t(bEGFR1_Q[withoutAA,])) ) ]                                                                                                                   
	Columnwise = cor((bEGFR1_Q[withoutAA,]))[ upper.tri( cor((bEGFR1_Q[withoutAA,])) ) ]  
	pdf("EGFR.Q.col_row_correlation.pdf")
	boxplot(Rowwise,Columnwise, names=c("All Row-wise pairs","All Column-wise pairs"), ylab="Pearson Correlation (r)")    
	dev.off()     
	
	# KRAS
	is2_EGFR1_4A[141,which(is2_EGFR1_4A[141,]>3)]
	is2_EGFR1_4A[which(is2_EGFR1_4A[,141]>3),141]
	
	is2_EGFR1_8A[141,which(is2_EGFR1_8A[141,]>3)]
	is2_EGFR1_8A[which(is2_EGFR1_8A[,141]>3),141]
	
	is2_EGFR1_Q[141,which(is2_EGFR1_Q[141,]>3)]
	is2_EGFR1_Q[which(is2_EGFR1_Q[,141]>3),141]
	
}

## Final
LoadKRPData <- function(){
	#KRP_ROW_INDEX = c(18,1:16)
	#KRP_COL_INDEX = c(19,23,26,30,20,32,27,29,34,21,22,24,28,25,31,35,33,17,36)     

	NA_domain_genes = c("Irp", "Ffh", "Protein N", "RplK", "Rpl7ae", "RpsH", "SelB" )
	KH_domain_genes = c( "Nova2", "Nova1", "GTPase Era" ) 
	dsRBM_domain_genes = c( "Adar2", "Adar2", "Rnt1" )
	RRM_domain_genes = c( "Ncl", "Nab3", "Celf1", "Elavl1" )
	sort( NA_domain_genes )
	sort( KH_domain_genes )
	sort( dsRBM_domain_genes )
	sort( RRM_domain_genes )

	# domain collected
	#KRP_ROW_INDEX <<- c(18,2,4,6,3,7,12,15,8,1,9,10,11,14,16,13,5)
	#KRP_COL_INDEX <<- c(19,26,20,27,30,29,28,35,34,23,21,22,24,31,33,25,32,17,36)     
	#gene_names = c("Irp","Ffh","Protein N","Rpl7ae","RplK","RpsH","selB","Adar2","Adar2","RNT1","GTPase Era","Nova1","Nova2","CELF1","Elavl","NAB3","Ncl","GFP","MBP")

	#KRP_ROW_INDEX <<- c(18,7,2,3,4,12,6,1,8,9,10,11,5,13,14,15,16)   
	#KRP_COL_INDEX <<- c(19,29,26,30,20,28,27,23,34,21,22,24,32,25,31,35,33,17,36)
	#gene_names = c("Irp", "RpsH", "Ffh", "RplK", "Protein N", "selB", "Rpl7ae", "Nova2", "Nova1", "Adar2", "Adar2", "RNT1", "Ncl", "NAB3", "CELF1", "GTPase Era", "Elavl", "GFP", "MBP" )


	KRP_ROW_INDEX <<- c(18,7,2,3,4,12,6,8,1,15,9,10,11,5,13,14,16)   
	KRP_COL_INDEX <<- c(19,29,26,30,20,28,27,34,23,35,21,22,24,32,25,31,33,17,36)
	gene_names = c("Irp","RpsH","Ffh","RplK","Protein N","selB","Rpl7ae","Nova1","Nova2","GTPase Era","Adar2","Adar2","RNT1","Ncl","NAB3","CELF1","Elavl","GFP","MBP")

	

	KRP_1_W_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-19_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4,draw=FALSE)  # 3898811    Samples were switched.
	KRP_1_WH_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-19_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4,draw=FALSE)  # 3433923    Samples were switched.
	KRP_1_W_rpi <<- KRP_1_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_1_WH_rpi <<- KRP_1_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	colnames(KRP_1_WH_rpi) = gene_names

	KRP_2_W_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S5_W.exact.rpi.txt",0.5,4,draw=FALSE)   # 1524737
	KRP_2_WH_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S6_WH.exact.rpi.txt",0.5,4,draw=FALSE) # 1992468
	KRP_2_W_rpi <<- KRP_2_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_2_WH_rpi <<- KRP_2_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]	

	KRP_3_W_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S7_WD.exact.rpi.txt",0.5,4,draw=FALSE)   # 1310970
	KRP_3_WH_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S8_HD.exact.rpi.txt",0.5,4,draw=FALSE) # 1464242
	KRP_3_W_rpi <<- KRP_3_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_3_WH_rpi <<- KRP_3_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	
	KRP_4_W_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-02_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4,draw=FALSE)   # 2079100
	KRP_4_WH_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-02_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4,draw=FALSE) # 4344649
	KRP_4_W_rpi <<- KRP_4_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_4_WH_rpi <<- KRP_4_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]

	KRP_5_W_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4,draw=FALSE)   # 779456
	KRP_5_WH4AT1_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S2_H.25AT.exact.rpi.txt",0.5,4,draw=FALSE) # 1696841     
	KRP_5_WH4AT2_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S3_H.25AT.exact.rpi.txt",0.5,4,draw=FALSE) # 1371988      
	KRP_5_WH2AT_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S4_H.5AT.exact.rpi.txt",0.5,4,draw=FALSE) # 1657655      
	
	KRP_5_W_rpi <<- KRP_5_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_5_WH4AT1_rpi <<- KRP_5_WH4AT1_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_5_WH4AT2_rpi <<- KRP_5_WH4AT2_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_5_WH2AT_rpi <<- KRP_5_WH2AT_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]

	KRP_6_W10_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S1_W10.exact.rpi.txt",0.5,4,draw=TRUE)   # 788210
	KRP_6_WH10_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S2_WH10.exact.rpi.txt",0.5,4,draw=TRUE) # 1266682     
	KRP_6_W25_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S3_W25.exact.rpi.txt",0.5,4,draw=TRUE) # 1038997      
	KRP_6_WH25_rpi <<- drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S4_WH25.exact.rpi.txt",0.5,4,draw=TRUE) # 1500996      
	
	KRP_6_W10_rpi <<- KRP_6_W10_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_6_WH10_rpi <<- KRP_6_WH10_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_6_W25_rpi <<- KRP_6_W25_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_6_WH25_rpi <<- KRP_6_WH25_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]	  


	is2_KRP_1_WH_rpi <<- InteractionScores(KRP_1_W_rpi,NormalMixture2(KRP_1_WH_rpi,threshold = 0.5),1,9 )		
	is2_KRP_2_WH_rpi <<- InteractionScores(KRP_2_W_rpi,NormalMixture2(KRP_2_WH_rpi,threshold = 0.5),1,9 )		
	is2_KRP_3_WH_rpi <<- InteractionScores(KRP_3_W_rpi,NormalMixture2(KRP_3_WH_rpi,threshold = 0.5),1,9 )		
	is2_KRP_4_WH_rpi <<- InteractionScores(KRP_4_W_rpi,NormalMixture2(KRP_4_WH_rpi,threshold = 0.5),1,9 )
	is2_KRP_5_WH4AT1_rpi <<- InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH4AT1_rpi,threshold = 0.5),1,9 )		# Area under the curve: 
	is2_KRP_5_WH4AT2_rpi <<- InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH4AT2_rpi,threshold = 0.5),1,9 )		# Area under the curve:
	is2_KRP_5_WH2AT_rpi <<- InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH2AT_rpi,threshold = 0.5),1,9 )		# Area under the curve:
	is2_KRP_6_WH10_rpi <<- InteractionScores(KRP_6_W10_rpi,NormalMixture2(KRP_6_WH10_rpi,threshold = 0.5),1,9 )		# Area under the curve:
	is2_KRP_6_WH25_rpi <<- InteractionScores(KRP_6_W25_rpi,NormalMixture2(KRP_6_WH25_rpi,threshold = 0.5),1,9 )		# Area under the curve:
	
}
   
TestGammaMixture2 <- function(){     
	## For KRP
	GammaMixture2(KRP_1_WH_rpi,LogScale = FALSE)    
	GammaMixture2(KRP_2_WH_rpi,LogScale = FALSE)    
	GammaMixture2(KRP_3_WH_rpi,LogScale = FALSE)    
	GammaMixture2(KRP_4_WH_rpi,LogScale = FALSE)    
	
	GammaMixture2(KRP_5_WH4AT1_rpi,LogScale = FALSE)    
	GammaMixture2(KRP_5_WH4AT2_rpi,LogScale = FALSE) 
	GammaMixture2(KRP_5_WH2AT_rpi,LogScale = FALSE)       
    NormalMixture2(KRP_5_WH2AT_rpi)

	GammaMixture2(KRP_6_W10_rpi,LogScale = FALSE)    
	#GammaMixture2(KRP_6_WH10_rpi,LogScale = FALSE) 
	#GammaMixture2(KRP_6_W25_rpi,LogScale = FALSE)     
	GammaMixture2(KRP_6_WH25_rpi,LogScale = FALSE)    
	
	InteractionScores(KRP_5_W_rpi,GammaMixture2(KRP_5_WH4AT1_rpi,LogScale = FALSE),1)
	InteractionScores(KRP_5_W_rpi,GammaMixture2(KRP_5_WH4AT2_rpi,LogScale = FALSE),1)
	InteractionScores(KRP_5_W_rpi,GammaMixture2(KRP_5_WH2AT_rpi,LogScale = FALSE),1)    
	
	
	KRP_Performance(InteractionScores(KRP_5_W_rpi,GammaMixture2(KRP_5_WH2AT_rpi,LogScale = FALSE),1))
	
    ## For P170
	gm2_is_P170_1 = InteractionScores(bP170_1_W,GammaMixture2(bP170_1_Q,LogScale = FALSE),1,4)
	gm2_is_P170_2 = InteractionScores(bP170_2_W,GammaMixture2(bP170_2_Q,LogScale = FALSE),1,4)
	gm2_is_P170_3 = InteractionScores(bP170_3_W,GammaMixture2(bP170_3_Q,LogScale = FALSE),1,4)
	gm2_is_P170_4 = InteractionScores(bP170_4_SW,GammaMixture2(bP170_4_SQ,LogScale = FALSE),1,4)
	gm2_is_P170_5 = InteractionScores(bP170_5_SW,GammaMixture2(bP170_5_SQ,LogScale = FALSE),1,4)
	gm2_is_P170_6A1 = InteractionScores(bP170_6_SW1,GammaMixture2(bP170_6_SA1,LogScale = FALSE),1,4)  
	gm2_is_P170_6A2 = InteractionScores(bP170_6_SW2,GammaMixture2(bP170_6_SA2,LogScale = FALSE),1,4)

	CorMatrix(gm2_is_P170_6A1,gm2_is_P170_6A2)
	CorMatrix(bP170_6_SA1,bP170_6_SA2)    
	CorMatrix(bP170_6_SW1,bP170_6_SW2)    

	## For R75
	gm2_is_bm1_A = InteractionScores(bm1_W,GammaMixture2(bm1_A,LogScale = FALSE),1,4)        
	gm2_is_bm2_A = InteractionScores(bm2_W,GammaMixture2(bm2_A,LogScale = FALSE),1,4)        
	gm2_is_bm2_Q = InteractionScores(bm2_W,GammaMixture2(bm2_Q,LogScale = FALSE),1,4)        
	gm2_is_bm3_A = InteractionScores(bm3_W,GammaMixture2(bm3_A,LogScale = FALSE),1,4)        
	gm2_is_bm3_Q = InteractionScores(bm3_W,GammaMixture2(bm3_Q,LogScale = FALSE),1,4)        
	gm2_is_bm3_QL = InteractionScores(bm3_W,GammaMixture2(bm3_QL,LogScale = FALSE),1,4) 
	gm2_is_bm4_S4A = InteractionScores(bm4_SW,GammaMixture2(bm4_S4A,LogScale = FALSE),1,4)   
	gm2_is_bm4_S8AL = InteractionScores(bm4_SW,GammaMixture2(bm4_S8A,LogScale = FALSE),1,4)   
	gm2_is_bm4_SQ = InteractionScores(bm4_SW,GammaMixture2(bm4_SQ,LogScale = FALSE),1,4)   			       
	gm2_is_bm6_SQ = InteractionScores(bm6_SW,GammaMixture2(bm6_SQ,LogScale = FALSE),1,4)      
	gm2_is_bm7_SQ = InteractionScores(bm7_SW,GammaMixture2(bm7_SQ,LogScale = FALSE),1,4)      		  
}                                              
                                                   

if (FALSE){
	# 2016-12-22_MiSeq; Roth75-exp1; MGj46 ( - R75_41,WDR7 / - R75_72,PLEKHG7 )     
	bm1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17543_S1.ppi.txt.new",0.5) # no selection ; 		46_WC_R12D12    "sum = 1968532"
	bm1_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17544_S2.ppi.txt.new",0.5) # selection    ;    	46_AC_R12D12    "sum = 1541502" 
	bm1_RW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17545_S3.ppi.txt.new",0.5) # no selection, RCA 	46_WRC_R12D12   "sum = 2190445" 
	bm1_RA = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2016-12-22_MiSeq/Blastn/17546_S4.ppi.txt.new",0.5) # selection, RCA    	46_ARC_R12D12   "sum = 1877953"  
	
	# 2017-02-22_MiSeq; Roth75-exp2; R75_MGj51 ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )   
	bm2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S1.ppi.txt.new",0.5) # no selection; R75_MGj51  "sum = 960245"
	bm2_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S2.ppi.txt.new",0.5) # selection; ACCRD         "sum = 2024754"
	bm2_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-02-22_MiSeq/Blastn/S3.ppi.txt.new",0.5) # selection; QCCRD         "sum = 813912"     # check 6503 - BCL2L2 - BIK
	
	# 2017-03-03_MiSeq; Roth75-exp3; [53]  ( - R75_41,WDR7 / - R75_72,PLEKHG7 / - AA )  
	bm3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S1.ppi.txt.new",0.5) # no selection; 53_WCCRD	18689	24,8%	"sum = "
	bm3_A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S2.ppi.txt.new",0.5) # selection; 53_ACCRD	18690	24,8%	"sum = "
	bm3_Q = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S3.ppi.txt.new",0.5) # selection; 53_QCCRD	18691	24,8%	"sum = 1255002"         ## NOT USE (LOW QUALITY)
	bm3_QL = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-03-03_MiSeq/Blastn/S7.ppi.txt.new",0.5) # selection; 53_QL4RD 18715	"sum = 734766"      
	
	# 2017-08-22_MiSeq; Roth75-exp4 with Seaprep; [61] 
	bm4_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SWD.ppi.txt.new",0.5) # Roth; -W;                 "sum = 1491017"
	bm4_S4A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA4D.ppi.txt.new",0.5) # Roth; -W/Aba 1/4;   "sum = 1177404"
	bm4_S8A = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SA8D.ppi.txt.new",0.5) # Roth; -W/Aba 1/8;     "sum = 1125083"
	bm4_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-08-22_MiSeq/Blastn/S64_SQD.ppi.txt.new",0.5) # Roth; -W/-H/-A;      "sum = 618177"
		
	# 2017-06-08_MiSeq; Roth75-exp5 with Seaprep only -W;                                                                        
	bm5_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-06-08_MiSeq/Blastn/S52.ppi.txt.new",0.5) # Roth; -W; Seaprep; R75_01 ~ R75_07 spike in                  "sum = 2801344"   
	                         
	# 2017-11-03_MiSeq; Roth75 - technical repeat
	bm6_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S1_W.ppi.txt.new",0.5)
	bm6_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S2_Q.ppi.txt.new",0.5)
	bm7_SW = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S3_W.ppi.txt.new",0.5)
	bm7_SQ = drawMatrix("/Volumes/users/lserrano/jyang/work/Mireia/src/output/2017-11-03_MiSeq/Blastn/S4_Q.ppi.txt.new",0.5)
}
	
	


if (FALSE){     

		           
    #===========================================================================
	### RNA - PROTEIN INTERACTION                                               
	#===========================================================================         
	# KRP	single gibson cloning	SeaPrep	23082	P2516	23/10/17	13_5_W      # OK
	# KRP	single gibson cloning 	SeaPrep	23083	P2516	23/10/17	13_5_WH     # OK
	# KRP	single gibson cloning	SeaPrep	23237		30/10/2017	14_5_W          # OK
	# KRP	single gibson cloning	SeaPrep	23238		30/10/2017	14_5_WH         # OK
	# KRP	single gibson cloning	SeaPrep	23239		30/10/2017	73_WD           # OK
	# KRP	single gibson cloning	SeaPrep	23240		30/10/2017	73_HD           # OK
	# KRP	single gibson cloning	SeaPrep	23344	P2545	3/11/2017	16_5_W      # OK
	# KRP	single gibson cloning	SeaPrep	23345	P2545	3/11/2017	16_5_WH     # OK
	# KRP	single gibson cloning	SeaPrep				
	# KRP	single gibson cloning	SeaPrep				
	# KRP	single gibson cloning	SeaPrep				
	# KRP	single gibson cloning	SeaPrep

	KRP_ROW_INDEX = c(18,1:16)
	KRP_COL_INDEX = c(19,23,26,30,20,32,27,29,34,21,22,24,28,25,31,35,33,17,36)     
	
	KRP_1_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-19_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4,draw=FALSE)  # 3898811    Samples were switched.
	KRP_1_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-19_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4,draw=FALSE)  # 3433923    Samples were switched.
	KRP_1_W_rpi = KRP_1_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_1_WH_rpi = KRP_1_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]     
	simple_plotLog10Matrix(KRP_1_W_rpi,1,9, filename = "KRP_1_W_rpi.pdf")  
    simple_plotLog10Matrix(KRP_1_W_rpi,1,9, filename = "KRP_1_W_rpi.pdf", breaks=seq(0,6,by=6/9))
    simple_plotLog10Matrix(KRP_1_WH_rpi,1,9, filename = "KRP_1_WH_rpi.pdf",breaks=seq(0,6,by=6/9)) 
                                              
	#freq_plotLog10Matrix(KRP_1_WH_rpi,1,9,step=5, breaks=seq(-6,-0.8,by=5.2/5)) 
	              
	freq_plotLog10Matrix(KRP_1_W_rpi,1,9,filename = "KRP_1_W_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100)) 
    freq_plotLog10Matrix(KRP_1_WH_rpi,1,9,filename = "KRP_1_WH_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100))

    # KRP_2 and KRP_3 are technical duplicates

	KRP_2_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S5_W.exact.rpi.txt",0.5,4,draw=FALSE)   # 1524737
	KRP_2_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S6_WH.exact.rpi.txt",0.5,4,draw=FALSE) # 1992468
	KRP_2_W_rpi = KRP_2_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_2_WH_rpi = KRP_2_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
    simple_plotLog2Matrix(KRP_2_W_rpi,1,9, filename = "KRP_2_W_rpi.pdf", breaks=seq(0,6,by=6/9))           
    simple_plotLog2Matrix(KRP_2_WH_rpi,1,9, filename = "KRP_2_WH_rpi.pdf", breaks=seq(0,6,by=6/9))          

	freq_plotLog10Matrix(KRP_2_W_rpi,1,9,filename = "KRP_2_W_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100)) 
    freq_plotLog10Matrix(KRP_2_WH_rpi,1,9,filename = "KRP_2_WH_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100))
	
	KRP_3_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S7_WD.exact.rpi.txt",0.5,4,draw=FALSE)   # 1310970
	KRP_3_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-30_MiSeq/Blastn/S8_HD.exact.rpi.txt",0.5,4,draw=FALSE) # 1464242
	KRP_3_W_rpi = KRP_3_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_3_WH_rpi = KRP_3_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
    simple_plotLog2Matrix(KRP_3_W_rpi,1,9, filename = "KRP_3_W_rpi.pdf", breaks=seq(0,6,by=6/9))           
    simple_plotLog2Matrix(KRP_3_WH_rpi,1,9, filename = "KRP_3_WH_rpi.pdf", breaks=seq(0,6,by=6/9))           

	freq_plotLog10Matrix(KRP_3_W_rpi,1,9,filename = "KRP_3_W_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100)) 
    freq_plotLog10Matrix(KRP_3_WH_rpi,1,9,filename = "KRP_3_WH_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100))

	
	KRP_4_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-02_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4,draw=FALSE)   # 2079100
	KRP_4_WH_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-02_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4,draw=FALSE) # 4344649
	KRP_4_W_rpi = KRP_4_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_4_WH_rpi = KRP_4_WH_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]					 
    simple_plotLog2Matrix(KRP_4_W_rpi,1,9, filename = "KRP_4_W_rpi.pdf")       
    simple_plotLog2Matrix(KRP_4_WH_rpi,1,9, filename = "KRP_4_WH_rpi.pdf")          

	freq_plotLog10Matrix(KRP_4_W_rpi,1,9,filename = "KRP_4_W_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100)) 
    freq_plotLog10Matrix(KRP_4_WH_rpi,1,9,filename = "KRP_4_WH_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100))

                                    
	# 3AT test
	KRP_5_W_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4,draw=FALSE)   # 779456
	KRP_5_WH4AT1_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S2_H.25AT.exact.rpi.txt",0.5,4,draw=FALSE) # 1696841     
	KRP_5_WH4AT2_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S3_H.25AT.exact.rpi.txt",0.5,4,draw=FALSE) # 1371988      
	KRP_5_WH2AT_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-15_MiSeq/Blastn/S4_H.5AT.exact.rpi.txt",0.5,4,draw=FALSE) # 1657655      
	
	KRP_5_W_rpi = KRP_5_W_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_5_WH4AT1_rpi = KRP_5_WH4AT1_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_5_WH4AT2_rpi = KRP_5_WH4AT2_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_5_WH2AT_rpi = KRP_5_WH2AT_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]         
	
    simple_plotLog2Matrix(KRP_5_W_rpi,1,9, filename = "KRP_5_W_rpi.pdf")       
    simple_plotLog2Matrix(KRP_5_WH4AT1_rpi,1,9, filename = "KRP_5_WH4AT1_rpi.pdf")      
	simple_plotLog2Matrix(KRP_5_WH4AT2_rpi,1,9, filename = "KRP_5_WH4AT2_rpi.pdf")      
	simple_plotLog2Matrix(KRP_5_WH2AT_rpi,1,9, filename = "KRP_5_WH2AT_rpi.pdf")          

	freq_plotLog10Matrix(KRP_5_W_rpi,1,9,filename = "KRP_5_W_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100)) 
    freq_plotLog10Matrix(KRP_5_WH4AT1_rpi,1,9,filename = "KRP_5_WH4AT1_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100))

	
	# DNA amount test
	KRP_6_W10_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S1_W10.exact.rpi.txt",0.5,4,draw=TRUE)   # 788210
	KRP_6_WH10_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S2_WH10.exact.rpi.txt",0.5,4,draw=TRUE) # 1266682     
	KRP_6_W25_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S3_W25.exact.rpi.txt",0.5,4,draw=TRUE) # 1038997      
	KRP_6_WH25_rpi = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-11-20_MiSeq/Blastn/S4_WH25.exact.rpi.txt",0.5,4,draw=TRUE) # 1500996      
	
	KRP_6_W10_rpi = KRP_6_W10_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
	KRP_6_WH10_rpi = KRP_6_WH10_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_6_W25_rpi = KRP_6_W25_rpi[KRP_ROW_INDEX,KRP_COL_INDEX] 
	KRP_6_WH25_rpi = KRP_6_WH25_rpi[KRP_ROW_INDEX,KRP_COL_INDEX]
		
    simple_plotLog2Matrix(KRP_6_W10_rpi,1,9, filename = "KRP_6_W10_rpi.pdf")       
    simple_plotLog2Matrix(KRP_6_WH10_rpi,1,9, filename = "KRP_6_WH10_rpi.pdf")      
	simple_plotLog2Matrix(KRP_6_W25_rpi,1,9, filename = "KRP_6_W25_rpi.pdf")      
	simple_plotLog2Matrix(KRP_6_WH25_rpi,1,9, filename = "KRP_6_WH25_rpi.pdf")          

	#freq_plotLog10Matrix(KRP_5_W_rpi,1,9,filename = "KRP_5_W_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100)) 
    #freq_plotLog10Matrix(KRP_5_WH4AT1_rpi,1,9,filename = "KRP_5_WH4AT1_rpi.pdf",step=100, breaks=seq(-7,-0.95,by=6.05/100))	   				 



	cis_KRP_1_WH_rpi = ColumnWiseInteractionScores(KRP_1_W_rpi,KRP_1_WH_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_1_WH_rpi.pdf" )   # Area under the curve: 0.7753
	cis_KRP_2_WH_rpi = ColumnWiseInteractionScores(KRP_2_W_rpi,KRP_2_WH_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_2_WH_rpi.pdf")   # Area under the curve: 0.9123
	cis_KRP_3_WH_rpi = ColumnWiseInteractionScores(KRP_3_W_rpi,KRP_3_WH_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_3_WH_rpi.pdf")   # Area under the curve: 0.8376
	cis_KRP_4_WH_rpi = ColumnWiseInteractionScores(KRP_4_W_rpi,KRP_4_WH_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_4_WH_rpi.pdf")   # Area under the curve: 0.8722     
	cis_KRP_5_WH4AT1_rpi = ColumnWiseInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT1_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_5_WH4AT1_rpi.pdf")   # Area under the curve: 0.8178     
	cis_KRP_5_WH4AT2_rpi = ColumnWiseInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT2_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_5_WH4AT2_rpi.pdf")   # Area under the curve: 0.8062    
	cis_KRP_5_WH2AT_rpi = ColumnWiseInteractionScores(KRP_5_W_rpi,KRP_5_WH2AT_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_5_WH2AT_rpi.pdf")   # Area under the curve: 0.8952   
    cis_KRP_6_WH10_rpi = ColumnWiseInteractionScores(KRP_6_W10_rpi,KRP_6_WH10_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_6_WH10_rpi.pdf")		# Area under the curve:
	cis_KRP_6_WH25_rpi = ColumnWiseInteractionScores(KRP_6_W25_rpi,KRP_6_WH25_rpi,1,9,rowNoralization = FALSE, filename = "cis_KRP_6_WH25_rpi.pdf")		# Area under the curve:	  
	
	mean( c( 0.7753,	0.9123,	0.8376,	0.8722 ) )	# 0.84935	             
	
	cis_KRP_1_WH_rpi_row = ColumnWiseInteractionScores(KRP_1_W_rpi,KRP_1_WH_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_1_WH_rpi_row.pdf")  # Area under the curve: 0.7764
	cis_KRP_2_WH_rpi_row = ColumnWiseInteractionScores(KRP_2_W_rpi,KRP_2_WH_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_2_WH_rpi_row.pdf")  # Area under the curve: 0.9373
	cis_KRP_3_WH_rpi_row = ColumnWiseInteractionScores(KRP_3_W_rpi,KRP_3_WH_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_3_WH_rpi_row.pdf")  # Area under the curve: 0.8431
	cis_KRP_4_WH_rpi_row = ColumnWiseInteractionScores(KRP_4_W_rpi,KRP_4_WH_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_4_WH_rpi_row.pdf")  # Area under the curve: 0.8714   
	cis_KRP_5_WH4AT1_rpi_row = ColumnWiseInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT1_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_5_WH4AT1_rpi_row.pdf")   # Area under the curve:      
	cis_KRP_5_WH4AT2_rpi_row = ColumnWiseInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT2_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_5_WH4AT2_rpi_row.pdf")   # Area under the curve:     
	cis_KRP_5_WH2AT_rpi_row = ColumnWiseInteractionScores(KRP_5_W_rpi,KRP_5_WH2AT_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_5_WH2AT_rpi_row.pdf")   # Area under the curve: 
    cis_KRP_6_WH10_rpi_row = ColumnWiseInteractionScores(KRP_6_W10_rpi,KRP_6_WH10_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_6_WH10_rpi_row.pdf")		# Area under the curve:
	cis_KRP_6_WH25_rpi_row = ColumnWiseInteractionScores(KRP_6_W25_rpi,KRP_6_WH25_rpi,1,9,rowNoralization = TRUE, filename = "cis_KRP_6_WH25_rpi_row.pdf")		# Area under the curve:
	mean( c(0.7764,	0.9373,	0.8431,	0.8714) )	# 0.85705
	
	simple_plotMatrix(cis_KRP_1_WH_rpi>0.3,0)  
	simple_plotMatrix(cis_KRP_2_WH_rpi>0.3,0)  
	simple_plotMatrix(cis_KRP_3_WH_rpi>0.3,0)  
	simple_plotMatrix(cis_KRP_4_WH_rpi>0.3,0)   
	simple_plotMatrix(cis_KRP_5_WH4AT1_rpi>0.3,0)   
	simple_plotMatrix(cis_KRP_5_WH4AT2_rpi>0.3,0)   
	simple_plotMatrix(cis_KRP_5_WH2AT_rpi>0.3,0)   
	
	
	
	# Pair-wise
	simple_plotLog2Matrix(KRP_4_WH_rpi/(KRP_4_W_rpi+1),1,9)       
	simple_plotLog2Matrix(KRP_4_WH_rpi/(KRP_4_W_rpi+1),1,9)       
	simple_plotLog2Matrix(KRP_4_WH_rpi/(KRP_4_W_rpi+1),1,9)       
	simple_plotLog2Matrix(KRP_4_WH_rpi/(KRP_4_W_rpi+1),1,9)                     
	
	pis_KRP_1_WH_rpi = PairInteractionScores(KRP_1_W_rpi,KRP_1_WH_rpi,1,9, filename = "pis_KRP_1_WH_rpi.pdf")	   
	pis_KRP_2_WH_rpi = PairInteractionScores(KRP_2_W_rpi,KRP_2_WH_rpi,1,9, filename = "pis_KRP_2_WH_rpi.pdf")		# Area under the curve: 0.7657
	pis_KRP_3_WH_rpi = PairInteractionScores(KRP_3_W_rpi,KRP_3_WH_rpi,1,9, filename = "pis_KRP_3_WH_rpi.pdf")		# Area under the curve: 0.728
	pis_KRP_4_WH_rpi = PairInteractionScores(KRP_4_W_rpi,KRP_4_WH_rpi,1,9, filename = "pis_KRP_4_WH_rpi.pdf")		# Area under the curve: 0.7418
	pis_KRP_5_WH4AT1_rpi = PairInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT1_rpi,1,9, filename = "pis_KRP_5_WH4AT1_rpi.pdf")		# Area under the curve: 
	pis_KRP_5_WH4AT2_rpi = PairInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT2_rpi,1,9, filename = "pis_KRP_5_WH4AT2_rpi.pdf")		# Area under the curve:
	pis_KRP_5_WH2AT_rpi = PairInteractionScores(KRP_5_W_rpi,KRP_5_WH2AT_rpi,1,9, filename = "pis_KRP_5_WH2AT_rpi.pdf")		# Area under the curve:
    pis_KRP_6_WH10_rpi = PairInteractionScores(KRP_6_W10_rpi,KRP_6_WH10_rpi,1,9, filename = "pis_KRP_6_WH10_rpi.pdf")		# Area under the curve:
	pis_KRP_6_WH25_rpi = PairInteractionScores(KRP_6_W25_rpi,KRP_6_WH25_rpi,1,9, filename = "pis_KRP_6_WH25_rpi.pdf")		# Area under the curve:
	   	
	nis_KRP_1_WH_rpi = NewInteractionScores(KRP_1_W_rpi,KRP_1_WH_rpi,1,9, filename = "nis_KRP_1_WH_rpi.pdf")		# Area under the curve: 0.7345
	nis_KRP_2_WH_rpi = NewInteractionScores(KRP_2_W_rpi,KRP_2_WH_rpi,1,9, filename = "nis_KRP_2_WH_rpi.pdf")		# Area under the curve: 0.7657
	nis_KRP_3_WH_rpi = NewInteractionScores(KRP_3_W_rpi,KRP_3_WH_rpi,1,9, filename = "nis_KRP_3_WH_rpi.pdf")		# Area under the curve: 0.728
	nis_KRP_4_WH_rpi = NewInteractionScores(KRP_4_W_rpi,KRP_4_WH_rpi,1,9, filename = "nis_KRP_4_WH_rpi.pdf")		# Area under the curve: 0.7418    
	nis_KRP_5_WH4AT1_rpi = NewInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT1_rpi,1,9, filename = "nis_KRP_5_WH4AT1_rpi.pdf")		# Area under the curve: 
	nis_KRP_5_WH4AT2_rpi = NewInteractionScores(KRP_5_W_rpi,KRP_5_WH4AT2_rpi,1,9, filename = "nis_KRP_5_WH4AT2_rpi.pdf")		# Area under the curve:
	nis_KRP_5_WH2AT_rpi = NewInteractionScores(KRP_5_W_rpi,KRP_5_WH2AT_rpi,1,9, filename = "nis_KRP_5_WH2AT_rpi.pdf")		# Area under the curve:  
    nis_KRP_6_WH10_rpi = NewInteractionScores(KRP_6_W10_rpi,KRP_6_WH10_rpi,1,9, filename = "nis_KRP_6_WH10_rpi.pdf")		# Area under the curve:
	nis_KRP_6_WH25_rpi = NewInteractionScores(KRP_6_W25_rpi,KRP_6_WH25_rpi,1,9, filename = "nis_KRP_6_WH25_rpi.pdf")		# Area under the curve:   

	is2_KRP_1_WH_rpi = InteractionScores(KRP_1_W_rpi,NormalMixture2(KRP_1_WH_rpi,threshold = 0.5),1,9, filename = "is2_KRP_1_WH_rpi.pdf")		
	is2_KRP_2_WH_rpi = InteractionScores(KRP_2_W_rpi,NormalMixture2(KRP_2_WH_rpi,threshold = 0.5),1,9, filename = "is2_KRP_2_WH_rpi.pdf")		
	is2_KRP_3_WH_rpi = InteractionScores(KRP_3_W_rpi,NormalMixture2(KRP_3_WH_rpi,threshold = 0.5),1,9, filename = "is2_KRP_3_WH_rpi.pdf")		
	is2_KRP_4_WH_rpi = InteractionScores(KRP_4_W_rpi,NormalMixture2(KRP_4_WH_rpi,threshold = 0.5),1,9, filename = "is2_KRP_4_WH_rpi.pdf")
	is2_KRP_5_WH4AT1_rpi = InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH4AT1_rpi,threshold = 0.5),1,9, filename = "is2_KRP_5_WH4AT1_rpi.pdf")		# Area under the curve: 
	is2_KRP_5_WH4AT2_rpi = InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH4AT2_rpi,threshold = 0.5),1,9, filename = "is2_KRP_5_WH4AT2_rpi.pdf")		# Area under the curve:
	is2_KRP_5_WH2AT_rpi = InteractionScores(KRP_5_W_rpi,NormalMixture2(KRP_5_WH2AT_rpi,threshold = 0.5),1,9, filename = "is2_KRP_5_WH2AT_rpi.pdf")		# Area under the curve:
	is2_KRP_6_WH10_rpi = InteractionScores(KRP_6_W10_rpi,NormalMixture2(KRP_6_WH10_rpi,threshold = 0.5),1,9, filename = "is2_KRP_6_WH10_rpi.pdf")		# Area under the curve:
	is2_KRP_6_WH25_rpi = InteractionScores(KRP_6_W25_rpi,NormalMixture2(KRP_6_WH25_rpi,threshold = 0.5),1,9, filename = "is2_KRP_6_WH25_rpi.pdf")		# Area under the curve:
	
		
	is_KRP_1_WH_rpi = InteractionScores(KRP_1_W_rpi,KRP_1_WH_rpi,1,9, filename = "is_KRP_1_WH_rpi.pdf")		
	is_KRP_2_WH_rpi = InteractionScores(KRP_2_W_rpi,KRP_2_WH_rpi,1,9, filename = "is_KRP_2_WH_rpi.pdf")		
	is_KRP_3_WH_rpi = InteractionScores(KRP_3_W_rpi,KRP_3_WH_rpi,1,9, filename = "is_KRP_3_WH_rpi.pdf")		
	is_KRP_4_WH_rpi = InteractionScores(KRP_4_W_rpi,KRP_4_WH_rpi,1,9, filename = "is_KRP_4_WH_rpi.pdf")
	is_KRP_5_WH4AT1_rpi = InteractionScores(KRP_5_W_rpi,KRP_5_WH4AT1_rpi,1,9, filename = "is_KRP_5_WH4AT1_rpi.pdf")		# Area under the curve: 
	is_KRP_5_WH4AT2_rpi = InteractionScores(KRP_5_W_rpi,KRP_5_WH4AT2_rpi,1,9, filename = "is_KRP_5_WH4AT2_rpi.pdf")		# Area under the curve:
	is_KRP_5_WH2AT_rpi = InteractionScores(KRP_5_W_rpi,KRP_5_WH2AT_rpi,1,9, filename = "is_KRP_5_WH2AT_rpi.pdf")		# Area under the curve:			
    is_KRP_6_WH10_rpi = InteractionScores(KRP_6_W10_rpi,KRP_6_WH10_rpi,1,9, filename = "is_KRP_6_WH10_rpi.pdf")		# Area under the curve:
	is_KRP_6_WH25_rpi = InteractionScores(KRP_6_W25_rpi,KRP_6_WH25_rpi,1,9, filename = "is_KRP_6_WH25_rpi.pdf")		# Area under the curve:

	   	
	# PERFORMANCE CHECK
	# KRP_Performance( NormalMixture2(KRP_1_WH_rpi,threshold = 0.5) )
	
	KRP_Performance( cis_KRP_1_WH_rpi, filename="cis_KRP_1_WH_rpi.auc.pdf" )        # Area under the curve: 0.7753
	KRP_Performance( cis_KRP_2_WH_rpi, filename="cis_KRP_2_WH_rpi.auc.pdf" )         # Area under the curve: 0.9123
	KRP_Performance( cis_KRP_3_WH_rpi, filename="cis_KRP_3_WH_rpi.auc.pdf" )         # Area under the curve: 0.8376
	KRP_Performance( cis_KRP_4_WH_rpi, filename="cis_KRP_4_WH_rpi.auc.pdf" ) 		   # Area under the curve: 0.8722    
    KRP_Performance( cis_KRP_5_WH4AT1_rpi, filename="cis_KRP_5_WH4AT1_rpi.auc.pdf" ) 		   # Area under the curve: 0.8178
	KRP_Performance( cis_KRP_5_WH4AT2_rpi, filename="cis_KRP_5_WH4AT2_rpi.auc.pdf" ) 		   # Area under the curve: 0.8062
	KRP_Performance( cis_KRP_5_WH2AT_rpi, filename="cis_KRP_5_WH2AT_rpi.auc.pdf" ) 		   # Area under the curve: 0.8952
	KRP_Performance( cis_KRP_6_WH10_rpi, filename="cis_KRP_6_WH10_rpi.auc.pdf" )        # Area under the curve:0.8231
	KRP_Performance( cis_KRP_6_WH25_rpi, filename="cis_KRP_6_WH25_rpi.auc.pdf" )         # Area under the curve: 0.8696
	
	KRP_Performance( cis_KRP_1_WH_rpi_row, filename="cis_KRP_1_WH_rpi_row.auc.pdf" ) 		# Area under the curve: 0.7764   
	KRP_Performance( cis_KRP_2_WH_rpi_row, filename="cis_KRP_2_WH_rpi_row.auc.pdf" ) 		# Area under the curve: 0.9373   
	KRP_Performance( cis_KRP_3_WH_rpi_row, filename="cis_KRP_3_WH_rpi_row.auc.pdf" ) 		# Area under the curve: 0.8431   
	KRP_Performance( cis_KRP_4_WH_rpi_row, filename="cis_KRP_4_WH_rpi_row.auc.pdf" ) 		# Area under the curve: 0.8714
    KRP_Performance( cis_KRP_5_WH4AT1_rpi_row, filename="cis_KRP_5_WH4AT1_rpi_row.auc.pdf" ) 		   # Area under the curve: 0.8232
	KRP_Performance( cis_KRP_5_WH4AT2_rpi_row, filename="cis_KRP_5_WH4AT2_rpi_row.auc.pdf" ) 		   # Area under the curve: 0.8047
	KRP_Performance( cis_KRP_5_WH2AT_rpi_row, filename="cis_KRP_5_WH2AT_rpi_row.auc.pdf" ) 		   # Area under the curve: 0.8998
	KRP_Performance( cis_KRP_6_WH10_rpi_row, filename="cis_KRP_6_WH10_rpi_row.auc.pdf" )        # Area under the curve: 0.8295
	KRP_Performance( cis_KRP_6_WH25_rpi_row, filename="cis_KRP_6_WH25_rpi_row.auc.pdf" )         # Area under the curve: 0.794
			                                         
	KRP_Performance( pis_KRP_1_WH_rpi, filename="pis_KRP_1_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7345  
	KRP_Performance( pis_KRP_2_WH_rpi, filename="pis_KRP_2_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7657   
	KRP_Performance( pis_KRP_3_WH_rpi, filename="pis_KRP_3_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.728   
	KRP_Performance( pis_KRP_4_WH_rpi, filename="pis_KRP_4_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7418    
    KRP_Performance( pis_KRP_5_WH4AT1_rpi, filename="pis_KRP_5_WH4AT1_rpi.auc.pdf" ) 		   # Area under the curve: 0.6895
	KRP_Performance( pis_KRP_5_WH4AT2_rpi, filename="pis_KRP_5_WH4AT2_rpi.auc.pdf" ) 		   # Area under the curve: 0.6944
	KRP_Performance( pis_KRP_5_WH2AT_rpi, filename="pis_KRP_5_WH2AT_rpi.auc.pdf" ) 		   # Area under the curve: 0.7563
	KRP_Performance( pis_KRP_6_WH10_rpi, filename="pis_KRP_6_WH10_rpi.auc.pdf" )        # Area under the curve:0.6387
	KRP_Performance( pis_KRP_6_WH25_rpi, filename="pis_KRP_6_WH25_rpi.auc.pdf" )         # Area under the curve: 0.7061
			
	KRP_Performance( nis_KRP_1_WH_rpi, filename="nis_KRP_1_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7545
	KRP_Performance( nis_KRP_2_WH_rpi, filename="nis_KRP_2_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7584
	KRP_Performance( nis_KRP_3_WH_rpi, filename="nis_KRP_3_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7412
	KRP_Performance( nis_KRP_4_WH_rpi, filename="nis_KRP_4_WH_rpi.auc.pdf" ) 		# Area under the curve: 0.7457
    KRP_Performance( nis_KRP_5_WH4AT1_rpi, filename="nis_KRP_5_WH4AT1_rpi.auc.pdf" ) 		   # Area under the curve:  0.7366
	KRP_Performance( nis_KRP_5_WH4AT2_rpi, filename="nis_KRP_5_WH4AT2_rpi.auc.pdf" ) 		   # Area under the curve: 0.7355
	KRP_Performance( nis_KRP_5_WH2AT_rpi, filename="nis_KRP_5_WH2AT_rpi.auc.pdf" ) 		   # Area under the curve: 0.7908
	KRP_Performance( nis_KRP_6_WH10_rpi, filename="nis_KRP_6_WH10_rpi.auc.pdf" )        # Area under the curve: 0.739
	KRP_Performance( nis_KRP_6_WH25_rpi, filename="nis_KRP_6_WH25_rpi.auc.pdf" )         # Area under the curve: 0.70.7591
		
	KRP_Performance( is_KRP_1_WH_rpi, filename="is_KRP_1_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7624
	KRP_Performance( is_KRP_2_WH_rpi, filename="is_KRP_2_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7814
	KRP_Performance( is_KRP_3_WH_rpi, filename="is_KRP_3_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7334
	KRP_Performance( is_KRP_4_WH_rpi, filename="is_KRP_4_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7459 
    KRP_Performance( is_KRP_5_WH4AT1_rpi, filename="is_KRP_5_WH4AT1_rpi.auc.pdf" ) 		   # Area under the curve: 0.7256
	KRP_Performance( is_KRP_5_WH4AT2_rpi, filename="is_KRP_5_WH4AT2_rpi.auc.pdf" ) 		   # Area under the curve: 0.6959
	KRP_Performance( is_KRP_5_WH2AT_rpi, filename="is_KRP_5_WH2AT_rpi.auc.pdf" ) 		   # Area under the curve: 0.8033
	KRP_Performance( is_KRP_6_WH10_rpi, filename="is_KRP_6_WH10_rpi.auc.pdf" )        # Area under the curve: 0.7341
	KRP_Performance( is_KRP_6_WH25_rpi, filename="is_KRP_6_WH25_rpi.auc.pdf" )         # Area under the curve: 0.7389
	
	KRP_Performance( is2_KRP_1_WH_rpi, filename="is2_KRP_1_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7624
	KRP_Performance( is2_KRP_2_WH_rpi, filename="is2_KRP_2_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7814
	KRP_Performance( is2_KRP_3_WH_rpi, filename="is2_KRP_3_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7334
	KRP_Performance( is2_KRP_4_WH_rpi, filename="is2_KRP_4_WH_rpi.auc.pdf" ) 			# Area under the curve: 0.7459 
    KRP_Performance( is2_KRP_5_WH4AT1_rpi, filename="is2_KRP_5_WH4AT1_rpi.auc.pdf" ) 		   # Area under the curve: 0.7256
	KRP_Performance( is2_KRP_5_WH4AT2_rpi, filename="is2_KRP_5_WH4AT2_rpi.auc.pdf" ) 		   # Area under the curve: 0.7603
	KRP_Performance( is2_KRP_5_WH2AT_rpi, filename="is2_KRP_5_WH2AT_rpi.auc.pdf" ) 		   # Area under the curve: 0.7595
	KRP_Performance( is2_KRP_6_WH10_rpi, filename="is2_KRP_6_WH10_rpi.auc.pdf" )        # Area under the curve: 0.717
	KRP_Performance( is2_KRP_6_WH25_rpi, filename="is2_KRP_6_WH25_rpi.auc.pdf" )         # Area under the curve: 0.7168
			     
    pdf( "KRP.AUC_BOX.pdf")
    boxplot( c(0.7624,0.7814,0.7334,0.7459),c(0.7345,0.7657,0.728,0.7418),c(0.7545,0.7584,0.7412,0.7457),c(0.7764,0.9373,0.8431,0.8714), names=c("Roth\'s","Element-wise","mRoth\'s","Column-wise"),ylab="AUC",xlab="Normalization Methods")  
	dev.off()       
	t.test(c(0.7624,0.7814,0.7334,0.7459,0.7345,0.7657,0.728,0.7418,0.7545,0.7584,0.7412,0.7457),c(0.7764,0.9373,0.8431,0.8714))  # p-value = 0.04697
	
	# Figire - Chart Correlation To see the effect
	pdf( "KRP_1_W_rpi.chart_cor.pdf")
	chart.Correlation(KRP_1_W_rpi) 
	dev.off()
	
	pdf( "KRP_1_W_rpi.t_chart_cor.pdf")
	chart.Correlation(t(KRP_1_W_rpi)) 
	dev.off()

	pdf( "KRP_1_WH_rpi.chart_cor.pdf")
	chart.Correlation(KRP_1_WH_rpi) 
	dev.off()
	
	pdf( "KRP_1_WH_rpi.t_chart_cor.pdf")
	chart.Correlation(t(KRP_1_WH_rpi)) 
	dev.off()      
	
	
	# Figure - Column correlation
	                                   
	# 136 = (17*17-17)/2
	RNAwise = c( (sum( cor(t(KRP_1_W_rpi))>0.95 ) - 17 )/2, (sum( cor(t(KRP_2_W_rpi))>0.95 ) - 17 )/2, (sum( cor(t(KRP_3_W_rpi))>0.95 ) - 17 )/2, (sum( cor(t(KRP_4_W_rpi))>0.95 ) - 17 )/2 )  
	# 171 = (19*19-19)/2
	#Proteinwise = c( (sum( cor((KRP_1_W_rpi))>0.95 ) - 19 )/2, (sum( cor((KRP_2_W_rpi))>0.95 ) - 19 )/2, (sum( cor((KRP_3_W_rpi))>0.95 ) - 19 )/2, (sum( cor((KRP_4_W_rpi))>0.95 ) - 19 )/2 ) 
	#pdf("KRP.col_row_correlation.pdf")
	#boxplot(RNAwise/136*100.0,Proteinwise/171*100.0, names=c("All RNA pairs","All Protein pairs"), ylab="Percentage of cases r > 0.95")
	#dev.off()   
	                                                     
	RNAwise = cor(t(KRP_2_W_rpi))[ upper.tri( cor(t(KRP_2_W_rpi)) ) ]
	Proteinwise = cor((KRP_2_W_rpi))[ upper.tri( cor((KRP_2_W_rpi)) ) ]  
	pdf("KRP.col_row_correlation.pdf")
	boxplot(RNAwise,Proteinwise, names=c("All RNA pairs","All Protein pairs"), ylab="Pearson Correlation (r)")
	dev.off()  
	t.test(RNAwise, Proteinwise) # p-value < 2.2e-16
	
	# 136 = (17*17-17)/2
	RNAwise2 = c( (sum( cor(t(KRP_1_WH_rpi))>0.95 ) - 17 )/2, (sum( cor(t(KRP_2_WH_rpi))>0.95 ) - 17 )/2, (sum( cor(t(KRP_3_WH_rpi))>0.95 ) - 17 )/2, (sum( cor(t(KRP_4_WH_rpi))>0.95 ) - 17 )/2 )  
	# 171 = (19*19-19)/2
	#Proteinwise2 = c( (sum( cor((KRP_1_WH_rpi))>0.95 ) - 19 )/2, (sum( cor((KRP_2_WH_rpi))>0.95 ) - 19 )/2, (sum( cor((KRP_3_WH_rpi))>0.95 ) - 19 )/2, (sum( cor((KRP_4_WH_rpi))>0.95 ) - 19 )/2 ) 
	#pdf("KRP.col_row_correlation.pdf")
	#boxplot(RNAwise/136*100.0,Proteinwise/171*100.0, names=c("All RNA pairs","All Protein pairs"), ylab="Percentage of cases r > 0.95")
	#dev.off()   
	                                                     
	RNAwise = cor(t(KRP_2_WH_rpi))[ upper.tri( cor(t(KRP_2_WH_rpi)) ) ]
	Proteinwise = cor((KRP_2_WH_rpi))[ upper.tri( cor((KRP_2_WH_rpi)) ) ]  
	pdf("KRP.WH.col_row_correlation.pdf")
	boxplot(RNAwise,Proteinwise, names=c("All RNA pairs","All Protein pairs"), ylab="Pearson Correlation (r)")
	dev.off()   
	        
	
	# comparison with DNA amount and 3AT
	DIAGONAL_MASK = matrix(0,nrow=17,ncol=19)
	diag(DIAGONAL_MASK) = 1
	
	boxplot(is2_KRP_6_WH25_rpi[which(DIAGONAL_MASK==0)],is2_KRP_6_WH10_rpi[which(DIAGONAL_MASK==0)],is2_KRP_5_WH2AT_rpi[which(DIAGONAL_MASK==0)])   
	t.test(is2_KRP_6_WH25_rpi[which(DIAGONAL_MASK==0)],is2_KRP_6_WH10_rpi[which(DIAGONAL_MASK==0)]) # p-value = 0.9613
	t.test(is2_KRP_5_WH2AT_rpi[which(DIAGONAL_MASK==0)],is2_KRP_6_WH10_rpi[which(DIAGONAL_MASK==0)]) # p-value = 1.225e-05        
	boxplot(log2(is2_KRP_6_WH25_rpi[which(DIAGONAL_MASK==0)]+1),log2(is2_KRP_6_WH10_rpi[which(DIAGONAL_MASK==0)]+1),log2(is2_KRP_5_WH2AT_rpi[which(DIAGONAL_MASK==0)]+1))       
                                                                                                                          
	snr1 = mean(is2_KRP_1_WH_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_1_WH_rpi[which(DIAGONAL_MASK==0)])	# 30.54821  
	snr2 = mean(is2_KRP_2_WH_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_2_WH_rpi[which(DIAGONAL_MASK==0)])	# 9.06232   
	snr3 = mean(is2_KRP_3_WH_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_3_WH_rpi[which(DIAGONAL_MASK==0)])	# 9.07446
	snr4 = mean(is2_KRP_4_WH_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_4_WH_rpi[which(DIAGONAL_MASK==0)])	# 10.20112   	
	snr5_1 = mean(is2_KRP_5_WH4AT1_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_5_WH4AT1_rpi[which(DIAGONAL_MASK==0)])	# 52.1762
	snr5_2 = mean(is2_KRP_5_WH4AT2_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_5_WH4AT2_rpi[which(DIAGONAL_MASK==0)])	# 75.82921
	snr5_3 = mean(is2_KRP_5_WH2AT_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_5_WH2AT_rpi[which(DIAGONAL_MASK==0)])   # 84.68874 
    snr6_1 =  mean(is2_KRP_6_WH10_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_6_WH10_rpi[which(DIAGONAL_MASK==0)])   # 8.211385     	 
	snr6_2 = mean(is2_KRP_6_WH25_rpi[which(DIAGONAL_MASK==1)]) / mean(is2_KRP_6_WH25_rpi[which(DIAGONAL_MASK==0)])   # 8.948108     
	pdf("KRP.SNR.pdf")
	barplot(c(snr6_1,snr6_2,snr5_1,snr5_2,snr5_3),xlab="Conditions", ylab="Signal To Noise Ratio (SNR)", names=c("25ng","100ng","1/4AT","1/4AT","1/2AT"))
	dev.off()        
	
	simple_plotMatrix((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi),0, filename="KRP.3AT.is2_sum.pdf")    
	simple_plotMatrix((is2_KRP_5_WH4AT1_rpi>3)+(is2_KRP_5_WH4AT2_rpi>3)+(is2_KRP_5_WH2AT_rpi>3),0, filename="KRP.3AT.is2_over3.pdf")  
	simple_plotMatrix((is2_KRP_5_WH4AT1_rpi>2)+(is2_KRP_5_WH4AT2_rpi>2)+(is2_KRP_5_WH2AT_rpi>2),0, filename="KRP.3AT.is2_over2.pdf")    
	
	# Flow chart

	simple_plotLog10Matrix( NormalMixture2(KRP_5_WH2AT_rpi/1657614.0,threshold = 0.5), 1/1657614.0, filename="KRP_5_WH2AT_rpi.filtered.log10Frq.pdf" )
	simple_plotLog10Matrix( KRP_5_WH2AT_rpi/1657614.0, 1/1657614.0, filename="KRP_5_WH2AT_rpi.log10Frq.pdf" )
	simple_plotLog10Matrix( KRP_5_W_rpi/779415.0, 1/779415.0, filename="KRP_5_W_rpi.log10Frq.pdf" )    
	   
	# Average with Sebastian way
	SUM_MASK = matrix(0,nrow=17,ncol=19)       
	SUM_MASK[which(is2_KRP_5_WH4AT1_rpi>0.0)] = 1
	SUM_MASK[which(is2_KRP_5_WH4AT2_rpi>0.0)] = SUM_MASK[which(is2_KRP_5_WH4AT2_rpi>0.0)] + 1
	SUM_MASK[which(is2_KRP_5_WH2AT_rpi>0.0)] = SUM_MASK[which(is2_KRP_5_WH2AT_rpi>0.0)] + 1
	
	simple_plotMatrix( ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi))/SUM_MASK,0, filename="KRP.3AT.is2_sebastian_sum.pdf" ) 

	KRP_3AT_is2_average = ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi))/3
	KRP_3AT_is2acc3_average = AutoActivationCorrection3(KRP_3AT_is2_average)
	simple_plotMatrix( ((is2_KRP_5_WH4AT1_rpi)+(is2_KRP_5_WH4AT2_rpi)+(is2_KRP_5_WH2AT_rpi))/3,0, filename="KRP.3AT.is2_average.pdf" ) 
	simple_plotMatrix( KRP_3AT_is2acc3_average,0, filename="KRP.3AT.is2acc3_average.pdf" ) 
}  
    
# rRNA data 
if (FALSE){  
	RBP_INDEX = c(1,59,85,19,27,47,68,15,34,9,18,37,46,48,54,57,39,38,87,3,4,5,6,7,8,10,11,12,13,14,16,17,20,21,22,23,24,25,26,28,29,30,31,32,33,35,36,40,41,42,43,44,45,49,50,51,52,53,55,56,58,60,61,62,63,64,65,66,67,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,86,88,89,90,91,92,93,94,95,96,97,101,2,98,99,100)
	RNA_INDEX = c(111,222,284,295,306,317,328,339,350,11,22,33,44,55,66,77,88,99,110,122,133,144,155,166,177,188,199,210,221,233,244,255,266,277,279,280,281,282,283,285,286,287,288,289,290,291,292,293,294,296,297,298,299,300,301,302,303,304,305,307,308,309,310,311,312,313,314,315,316,318,319,320,321,322,323,324,325,326,327,329,330,331,332,333,334,335,336,337,338,340,341,342,343,344,345,346,347,348,349,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,23,24,25,26,27,28,29,30,31,32,34,35,36,37,38,39,40,41,42,43,45,46,47,48,49,50,51,52,53,54,56,57,58,59,60,61,62,63,64,65,67,68,69,70,71,72,73,74,75,76,78,79,80,81,82,83,84,85,86,87,89,90,91,92,93,94,95,96,97,98,100,101,102,103,104,105,106,107,108,109,112,113,114,115,116,117,118,119,120,121,123,124,125,126,127,128,129,130,131,132,134,135,136,137,138,139,140,141,142,143,145,146,147,148,149,150,151,152,153,154,156,157,158,159,160,161,162,163,164,165,167,168,169,170,171,172,173,174,175,176,178,179,180,181,182,183,184,185,186,187,189,190,191,192,193,194,195,196,197,198,200,201,202,203,204,205,206,207,208,209,211,212,213,214,215,216,217,218,219,220,223,224,225,226,227,228,229,230,231,232,234,235,236,237,238,239,240,241,242,243,245,246,247,248,249,250,251,252,253,254,256,257,258,259,260,261,262,263,264,265,267,268,269,270,271,273,278,276,272,274,275)                                       
    
	# phusi   
	brRNA1_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-08-24_MiSeq/Blastn/S1_W.rpi.txt",0.5,4)
	brRNA1_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-08-24_MiSeq/Blastn/S2_WH.rpi.txt",0.5,4)
	
	brRNA2_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-09-07_MiSeq/Blastn/N6_W.rpi.txt",0.5,4)
	brRNA2_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-09-07_MiSeq/Blastn/N6_WH.rpi.txt",0.5,4)
	
	brRNA3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S1_W.rpi.txt",0.5,4)
	brRNA3_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S2_WH.rpi.txt",0.5,4)
	#brRNA3_W = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S1_W.exact.rpi.txt",0.5,4)
	#brRNA3_WH = drawMatrix("/Volumes/users/lserrano/jyang/work/Nele/src/output/2017-10-13_MiSeq/Blastn/S2_WH.exact.rpi.txt",0.5,4)    
	                                                                       
    simple_plotLog10Matrix(brRNA1_W[1:350,351:451][RNA_INDEX,RBP_INDEX],1,2,filename="rRNA.exp1.W.log10.pdf")
    simple_plotLog10Matrix(brRNA1_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1,2,filename="rRNA.exp1.WH.log10.pdf")    
    simple_plotLog10Matrix(brRNA2_W[1:350,351:451][RNA_INDEX,RBP_INDEX],1,2,filename="rRNA.exp2.W.log10.pdf")    
    simple_plotLog10Matrix(brRNA2_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1,2,filename="rRNA.exp2.WH.log10.pdf")    
    simple_plotLog10Matrix(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX],1,2,filename="rRNA.exp3.W.log10.pdf")    
    simple_plotLog10Matrix(brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1,2,filename="rRNA.exp3.WH.log10.pdf")    
	
	is_rRNA1 = InteractionScores(brRNA1_W[1:350,351:451],brRNA1_WH[1:350,351:451],1,4)             
	is_rRNA2 = InteractionScores(brRNA2_W[1:350,351:451],brRNA2_WH[1:350,351:451],1,4)             
	is_rRNA3 = InteractionScores(brRNA3_W[1:350,351:451],brRNA3_WH[1:350,351:451],1,4)             
	
	is_rRNA1_reorder = InteractionScores(brRNA1_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA1_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="is_rRNA.exp1.pdf")   
	is_rRNA2_reorder = InteractionScores(brRNA2_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA2_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="is_rRNA.exp2.pdf")   
	is_rRNA3_reorder = InteractionScores(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="is_rRNA.exp3.pdf")   
	
	pis_rRNA1_reorder = PairInteractionScores(brRNA1_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA1_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="pis_rRNA.exp1.pdf")   
	pis_rRNA2_reorder = PairInteractionScores(brRNA2_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA2_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="pis_rRNA.exp2.pdf")   
	pis_rRNA3_reorder = PairInteractionScores(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="pis_rRNA.exp3.pdf")
	
	nis_rRNA1_reorder = NewInteractionScores(brRNA1_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA1_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="nis_rRNA.exp1.pdf")   
	nis_rRNA2_reorder = NewInteractionScores(brRNA2_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA2_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="nis_rRNA.exp2.pdf")   
	nis_rRNA3_reorder = NewInteractionScores(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="nis_rRNA.exp3.pdf")
	
	cis_rRNA1_reorder = ColumnWiseInteractionScores(brRNA1_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA1_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="cis_rRNA.exp1.pdf")   
	cis_rRNA2_reorder = ColumnWiseInteractionScores(brRNA2_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA2_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="cis_rRNA.exp2.pdf")   
	cis_rRNA3_reorder = ColumnWiseInteractionScores(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX],brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX],1.0,2,filename="cis_rRNA.exp3.pdf")   
	
	is2_rRNA1_reorder = InteractionScores(brRNA1_W[1:350,351:451][RNA_INDEX,RBP_INDEX],NormalMixture2(brRNA1_WH[1:350,351:451][RNA_INDEX,RBP_INDEX]),1.0,2,filename="is2_rRNA.exp1.pdf")   
	is2_rRNA2_reorder = InteractionScores(brRNA2_W[1:350,351:451][RNA_INDEX,RBP_INDEX],NormalMixture2(brRNA2_WH[1:350,351:451][RNA_INDEX,RBP_INDEX]),1.0,2,filename="is2_rRNA.exp2.pdf")   
	is2_rRNA3_reorder = InteractionScores(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX],NormalMixture2(brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX]),1.0,2,filename="is2_rRNA.exp3.pdf")	
	
	m1 = brRNA3_W[1:350,351:451]
	m2 = brRNA3_WH[1:350,351:451]
	cis_rRNA3 = ColumnWiseInteractionScores(brRNA3_W[1:350,351:451],brRNA3_WH[1:350,351:451],1,2)                       
	cis_rRNA3 = ColumnWiseInteractionScores(brRNA3_W[1:350,351:451][RNA_INDEX,RBP_INDEX] ,brRNA3_WH[1:350,351:451][RNA_INDEX,RBP_INDEX] ,1,2, filename="cis.2017-10-13_MiSeq.rRNA.pdf") 
	simple_plotMatrix(cis_rRNA3>50,0,2)        
	
	simple_plotMatrix( (is2_rRNA1_reorder>3)+(is2_rRNA2_reorder>3)+(is2_rRNA3_reorder>3),0,4)   

	simple_plotMatrix( ((is2_rRNA1_reorder>3)+(is2_rRNA2_reorder>3)+(is2_rRNA3_reorder>3))>=3,0,4,filename="is2_rRNA.3.over3.pdf")     
	simple_plotMatrix( ((is2_rRNA1_reorder>3)+(is2_rRNA2_reorder>3)+(is2_rRNA3_reorder>3))>=2,0,4,filename="is2_rRNA.3.over2.pdf") 
	simple_plotMatrix( ((is2_rRNA1_reorder>3)+(is2_rRNA2_reorder>3)+(is2_rRNA3_reorder>3))>=1,0,4,filename="is2_rRNA.3.over1.pdf")
		 
	simple_plotMatrix( ((is2_rRNA1_reorder>2)+(is2_rRNA2_reorder>2)+(is2_rRNA3_reorder>2))>=3,0,4,filename="is2_rRNA.2.over3.pdf")     
	simple_plotMatrix( ((is2_rRNA1_reorder>2)+(is2_rRNA2_reorder>2)+(is2_rRNA3_reorder>2))>=2,0,4,filename="is2_rRNA.2.over2.pdf") 
	simple_plotMatrix( ((is2_rRNA1_reorder>2)+(is2_rRNA2_reorder>2)+(is2_rRNA3_reorder>2))>=1,0,4,filename="is2_rRNA.2.over1.pdf") 
	             
	simple_plotMatrix( ((is2_rRNA1_reorder>1.5)+(is2_rRNA2_reorder>1.5)+(is2_rRNA3_reorder>1.5))>=3,0,4,filename="is2_rRNA.1.5.over3.pdf")    
    simple_plotMatrix( ((is2_rRNA1_reorder>1.5)+(is2_rRNA2_reorder>1.5)+(is2_rRNA3_reorder>1.5))>=2,0,4,filename="is2_rRNA.1.5.over2.pdf")  
    simple_plotMatrix( ((is2_rRNA1_reorder>1.5)+(is2_rRNA2_reorder>1.5)+(is2_rRNA3_reorder>1.5))>=1,0,4,filename="is2_rRNA.1.5.over1.pdf")     
}



if (FALSE){
	# EB1-3 pulldown proteome

	EB1_1 = c(0.56,	2.7)
	EB1_2 = c(0.51,0.54,0.7,2.03,0.09,0.23,0.07,1.76)
	EB1_3 = c(0.66,6.68,0.09,0.34,0.12,0.18,0.96,0.29,0.31,0.07)

	EB2_1 = c(1.47,	1.34)
	EB2_2 = c(0.45,1.02,1.71,1.35,0.23,0.6,0.47,0.62,1.59,2.81,1.97)
	EB2_3 = c(0.5,0.36,0.68,1.54,0.5)
	
	EB3_1 = c(0.79,1,3.57)
	EB3_2 = c(0.78,1.43,1.04,2.92,0.86,0.14,0.03,1.51,2.47,0.38,0.71)
	EB3_3 = c(8.51,0.25,0.28,0.11,0.25,0.15,0.4,0.06)
	    
	L1 = c(EB1_1, EB2_1, EB3_1)
	L2 = c(EB1_2, EB2_2, EB3_2)
	L3 = c(EB1_3, EB2_3, EB3_3)
	
	boxplot(EB1_1,EB1_2,EB1_3)
	boxplot(EB2_1,EB2_2,EB2_3)
	boxplot(EB3_1,EB3_2,EB3_3)
	
	GRP = c(Ones(length(L1)),Ones(length(L2))*2,Ones(length(L3))*3) 
	LAYER = c(L1,L2,L3)         
	data = data.frame(GRP=GRP,LAYER=LAYER)        
	
	pdf( "EB1_3_pulldown.EMPAI.pdf")
	boxplot(LAYER ~ GRP, data = data, lwd = 2, ylab = 'EMPAI Score')
	stripchart(LAYER ~ GRP, data = data, vertical = TRUE, method = "jitter", add = TRUE, pch = 20, col = 'blue')   
	dev.off()
	
}

