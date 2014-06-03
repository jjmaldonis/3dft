import sys
import numpy as np


def main():
    """ takes a 3d matrix (of equal dimensions) written out on a single line
        (from fortran) and writes it to a file of format Igor Text Wave """
    infile = sys.argv[1]
    outfile = sys.argv[2]

    #with open(infile) as f:
    #    content = f.readlines()
    with open(infile) as f:
        for i, l in enumerate(f):
            pass
    i += 1

    npix = int(round(i**(1.0/3.0)))
    print("Number of pixels in each direction (assumed to be equal): {0}".format(npix))

    of = open(outfile,'w')
    of.write('IGOR\nWAVES/N=({0},{0},{0})\t {1}\nBEGIN\n'.format(npix,outfile[:outfile.index('.')]))
    i = 0
    subwave = ['                ']*(npix**2)
    with open(infile) as f:
        print("Opened infile for reading...")
        for line in f:
            subwave[ i%npix**2 ] = line.strip()
            if( i > 0 and i%npix**2 == 0):
                of.write('\t{0}\n'.format("\t".join(subwave)))
                print(i)
            i += 1
    of.write('\t{0}\n'.format("\t".join(subwave)))
    print(i)
    of.write('END\n')
    of.write('X SetScale/P x 0,1,"", {0}; SetScale/P y 0,1,"", {0}; SetScale/P z 0,1,"", {0}; SetScale d 0,0,"", {0}'.format(outfile[:outfile.index('.')]))
    of.close()

    #content = [float(x.strip()) for x in content]

    #wave = []
    #subwave = []
    #n = -1
    #for i,x in enumerate(content):
    #    if( i>0 and i%npix == 0.0 ):
    #        n += 1
    #        wave.append(subwave)
    #        subwave = []
    #    subwave.append(x)
    #wave.append(subwave)

    #wave = [ [] for x in xrange(npix**2) ]
    #print(len(wave))
    #l = len(wave)
    #n = -1
    #for i,x in enumerate(content):
    #    if(i%npix**2 == 0): n += 1
    #    #print([i,l,n,(n*npix) + (i%npix)])
    #    wave[(n*npix) + (i%npix)].append(x)

    #of = open(outfile,'w')
    #of.write('IGOR\nWAVES/N=({0},{0},{0})\t {1}\nBEGIN\n'.format(npix,outfile[:outfile.index('.')]))

    #wave = [ [str(x) for x in subwave] for subwave in wave]

    #for subwave in wave:
    #    of.write('\t{0}\n'.format("\t".join(subwave)))

    #of.write('END\n')
    #of.write('X SetScale/P x 0,1,"", {0}; SetScale/P y 0,1,"", {0}; SetScale/P z 0,1,"", {0}; SetScale d 0,0,"", {0}'.format(outfile[:outfile.index('.')]))
    #of.close()


if __name__ == "__main__":
    main()
