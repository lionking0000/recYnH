import os
VERBOSE = False

def run_cmd( cmd ):
    if VERBOSE: print cmd
    os.system( cmd )

def read_R_generated_matrix( file ):
    rownames = []
    colnames = []
    m = []
    f = open( file )
    for line in f.xreadlines():
        if line[0] == "#": continue
        rownames = line[:-1].split("\t")
        break
    for line in f.xreadlines():
        fields = line[:-1].split("\t")
        colnames.append( fields[0] )
        m.append( [ float(x) for x in fields[1:] ] )
    f.close()
    return m, rownames, colnames

def save_matrix( output_file_path, m, rownames, colnames ):
    # save originalIS
    fout = open( output_file_path, "w" )
    fout.write( "# This file contains 'Interaction Scores' generated by rec-YnH score function\n" )  
    output = "DB(Read 1) \ AD(Read 2)"	
    for name in rownames:
        output += "\t"
        output += name
    fout.write( output + "\n" )
    i = 0
    for x in m:
        j = 0
        output = colnames[i]
        for y in x:
            output += "\t"
            output += "%f" % y
            j += 1
        i += 1
        output += "\n"
        fout.write( output )
    fout.close()

def run( args ):
    #print args.program # the experiments type ('Y2H'|'Y3H')
    #print args.matrix1 # the interaction matrix of selection condition
    #print args.matrix2 # the interaction matrix of non-selection condition
    #print args.output # the output folder name
    #print args.name # output name

    if args.matrix2 == None:
        args.matrix2 = args.matrix1
    
    [ dirname1, m1 ] = os.path.split( args.matrix1 )
    [ dirname2, m2 ] = os.path.split( args.matrix2 )

    if args.output == None:
        args.output = dirname1
    
    # make output folder
    if os.path.exists( args.output ) == False:
        os.makedirs( args.output )
    
    temp_dir = os.path.join( args.output, "tmp" )
    if os.path.exists( temp_dir ) == False:
        os.makedirs( temp_dir )

    print "[ Normalizing matrix ]", args.matrix1, args.matrix2
    #cmd = "Rscript /usr/local/bin/visualization.R %s %s %s %s %s" % ( args.program, args.matrix1, args.matrix2, args.output, args.name )
    cmd = "Rscript ./src/visualization.R %s %s %s %s %s" % ( args.program, args.matrix1, args.matrix2, args.output, args.name )
    run_cmd( cmd )


    NIS_output_path = "%s/%s.nis.txt" % ( args.output, args.name )
    print "[ Saving Interaction Score matrix ]", NIS_output_path
    #IS_output_path = "%s/%s.is.txt" % ( args.output, args.name )
    
    #if os.path.exists( IS_output_path ):
    #    m, rownames, colnames = read_R_generated_matrix( IS_output_path )
    #    save_matrix( IS_output_path, m, rownames, colnames )

    if os.path.exists( NIS_output_path ):
        m, rownames, colnames = read_R_generated_matrix( NIS_output_path )
        save_matrix( NIS_output_path, m, rownames, colnames )
