######################################################################################################################################  
#
# Version 1.0
# visualization.R
#
#
######################################################################################################################################


source( "./src/init.R" )


#######################
# Read Argument
#######################
#print("[recYnH] Visualization : Read Argument")
argv=commandArgs(T)
#print( argv[1] )	# args.program # the experiments type ('Y2H'|'Y3H')
#print( argv[2] )	# args.matrix1 # the interaction matrix of selection condition
#print( argv[3] )	# args.matrix2 # the interaction matrix of non-selection condition
#print( argv[4] )	# args.output # the output folder name
#print( argv[5] )	# args.name # output name

IS_output_path <- sprintf( "%s/%s.is.txt", argv[4], argv[5] )
NIS_output_path <- sprintf( "%s/%s.nis.txt", argv[4], argv[5] )

PNG_is_output_path <- sprintf( "%s/%s.is.png", argv[4], argv[5] )
PNG_nis_output_path <- sprintf( "%s/%s.nis.png", argv[4], argv[5] )

#m1 = drawMatrix("../share/output/test/recYnH.ppi.txt",0.5) 
m1 = drawMatrix( argv[2], 0.5 ) 
#m2 = drawMatrix("../share/output/test/recYnH.ppi.txt",0.5) 
m2 = drawMatrix( argv[3], 0.5) 
is2 = InteractionScores(m1,m2,1.0)  
#nis2 = NewInteractionScores(m1,m2,1.0)  
nis2 = InteractionScores(m1,NormalMixture2(m2),1.0) 

#CorMatrix(is2,nis2)
#SaveMatrixPlot(is2, 1.0, PNG_is_output_path, 4 ) 
#SaveMatrixPlot(nis2, 1.0, PNG_nis_output_path, 4 ) 

write.table(is2,file=IS_output_path,sep="\t")
write.table(nis2,file=NIS_output_path,sep="\t")

quit()
