#===================================================================================================================================== 
#
# Version 1.0
# init.R
#
#
#=====================================================================================================================================  

#=====================================================================================================================================  
## Load libraries    

library("RColorBrewer")
library(pheatmap)
library(mixtools)

# Remove 0 values, and gaussian mixture to find signal
NormalMixture2 <-function(data, show.plot = FALSE, threshold = 0.5, non_zero_cut = 0){        
	set.seed(47)   
       
	output2 = data
	non_zero_data = data[which(data>non_zero_cut)]
	output = normalmixEM(as.vector(non_zero_data))

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


RawInteractionScores <- function(m1,m2, alpha,show=FALSE, filename=NA, cellwidth = NA, cellheight = NA){
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


InteractionScores <- function(m1,m2, alpha,fontsize=9, filename=NA,cellwidth = NA, cellheight = NA){              
	raw_is = RawInteractionScores(m1,m2,alpha,FALSE)
	raw_is[which(raw_is=="NaN")] = 0
	raw_is[which(raw_is==Inf)] = 0  
	raw_is[which(raw_is==-Inf)] = 0 
	#pheatmap(log2(raw_is+alpha), cluster_rows=FALSE, cluster_cols = FALSE, color = brewer.pal(9,"Blues"), fontsize=fontsize, filename=filename, cellwidth = cellwidth, cellheight = cellwidth)  
	return(log2(raw_is+alpha))
} 

drawMatrix <- function( ppi_output_path, alpha,fontsize=9,interactive=FALSE,draw=FALSE){
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
