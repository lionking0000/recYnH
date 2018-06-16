#!/opt/conda/bin/python2.7
import re

def read_fasta( filepath ):
    f = open(filepath,"r")
    data = f.readlines()
    f.close()

    sequence_dic = {}
    id = ""
    for line in data:
        line = line.strip()
        if line[0] == '>':
            id = line[1:]
            sequence_dic[ id ] = ""
        else:
            sequence_dic[id] += line
    return sequence_dic

def read_fasta_file( filepath ):
    return read_fasta( filepath )

def RevComplement( seq ):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def read( filepath, pattern = None, id_list = None ):
    out = {}
    f = open( filepath )

    if pattern == None:
        id = ""
        seq = ""
        for line in f.xreadlines():
            if line[0] == ">":
                if id != "":
                    if id_list != None:
                        id_list.append( id )
                    out[ id ] = seq
                    seq = ""
                id = line[1:-1]
            else:
                seq += line.strip()
        if id != "":
            if id_list != None:
                id_list.append( id )
            out[ id ] = seq
    else:
        id = ""
        seq = ""
        bFind = False
        for line in f.xreadlines():
            if line[0] == ">":
                found = re.search( pattern, line[1:-1] )
                if found != None:
                    if id != "":
                        if id_list != None:
                            id_list.append( id )
                        out[ id ] = seq
                        seq = ""
                    id = line[1:-1]
                    bFind = True
                else:
                    bFind = False
                continue

            if bFind == False:
                continue    
            else:
                seq += line.strip()
        if id != "":
            if id_list != None:
                id_list.append( id )
            out[ id ] = seq

    f.close()

    return out

def read_all( filepath ):
    out = {}
    f = open( filepath )

    id = ""
    seq = ""
    index = 0
    for line in f.xreadlines():
        if line[0] == ">":
            if id != "":
                out[ id ] = seq
                seq = ""
            id = line[1:-1] + "\t%d" % index
            index += 1
        else:
            seq += line.strip()
    if id != "":
        out[ id ] = seq
    return out