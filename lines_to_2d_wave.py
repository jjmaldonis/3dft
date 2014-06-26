import sys
import numpy as np


def main():
    """ takes a 3d matrix (of equal dimensions) written out on a single line
        (from fortran) and writes it to a file of format Igor Text Wave """
    infile = sys.argv[1]
    outfile = sys.argv[2]

    #with open(infile) as f:
    #    content = f.readlines()
    npix = 0
    with open(infile) as f:
        for line in f:
            if(len(line.strip()) > 0):
                npix += 1
    print("Number of pixels in each direction (assumed to be equal): {0}".format(npix))

    of = open(outfile,'w')
    of.write('IGOR\nWAVES/D/N=({0},{0}) {1}\nBEGIN\n'.format(npix,outfile[:outfile.index('.')]))
    i = 0
    with open(infile) as f:
        print("Opened infile for reading...")
        for line in f:
            line = line.strip()
            if(len(line) > 0):
                line = line.split()
                if(len(line) != npix):
                    raise Exception("Hey! Increase the spacing between the numbers in the fortran output! {0}".format(len(line)))
                of.write(' '.join(line)+'\n')
                i += 1
        if(i != npix):
            raise Exception("Number of lines was incorrect! {0}".format(i))
    of.write('END\n')
    of.write('X SetScale/P x 0,1,"", {0}; SetScale/P y 0,1,"", {0}'.format(outfile[:outfile.index('.')]))
    of.close()

if __name__ == "__main__":
    main()
