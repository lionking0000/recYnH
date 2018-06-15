#!/opt/conda/bin/python2.7
import cPickle

def save( p, file ):
    '''
        save( p, file )
    '''
    fout = open( file, "w" )
    cPickle.dump( p, fout )
    fout.close()

def load( file ):
    '''
        load( file )
    '''
    fin = open( file )
    p = cPickle.load( fin )
    fin.close()
    return p