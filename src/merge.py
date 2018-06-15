import os
VERBOSE = False

def run_cmd( cmd ):
    if VERBOSE: print cmd
    os.system( cmd )

def run( args ):
    print args.program # the experiments type ('Y2H'|'Y3H')
    print args.matrix1 # the interaction matrix of selection condition
    print args.matrix2 # the interaction matrix of non-selection condition
    print args.output # the output folder name
    print args.name # output name

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

        
    cmd = "Rscript /usr/local/bin/visualization.R %s %s %s %s %s" % ( args.program, args.matrix1, args.matrix2, args.output, args.name )
    run_cmd( cmd )
