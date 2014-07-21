import sys


def lines_to_3d_wave(infile,outfile):
    """ takes a 3d matrix (of equal dimensions) written out on a single line
        (from fortran) and writes it to a file of format Igor Text Wave """
    #with open(infile) as f:
    #    content = f.readlines()
    npix = 0
    with open(infile) as f:
        for line in f:
            if(len(line.strip()) > 0):
                npix += 1
    print("Number of pixels in each direction (assumed to be equal): {0}".format(npix))

    of = open(outfile,'w')
    of.write('IGOR\nWAVES/N=({0},{0},{0})\t {1}\nBEGIN\n'.format(npix,outfile[:outfile.index('.')]))
    #i = 0
    #subwave = ['                ']*(npix**2)
    with open(infile) as f:
        print("Opened infile {0} for reading...".format(infile))
        for i,line in enumerate(f):
            #subwave[ i%npix**2 ] = line.strip()
            #if( i > 0 and i%npix**2 == 0):
            #    of.write('\t{0}\n'.format("\t".join(subwave)))
            #    print(i)
            #i += 1
            line = line.strip().split()
            if(len(line) != npix**2):
                print("Wrong number of elements on line {0}! {1} {2}".format(i,len(line),npix**2))
            of.write(' '.join(line)+'\n')
            #of.write(' '.join(line.strip().split())+'\n')
    #of.write('\t{0}\n'.format("\t".join(subwave)))
    of.write('END\n')
    of.write('X SetScale/P x 0,1,"", {0}; SetScale/P y 0,1,"", {0}; SetScale/P z 0,1,"", {0}; SetScale d 0,0,"", {0}'.format(outfile[:outfile.index('.')]))
    of.close()

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    lines_to_3d_wave(infile,outfile)

if __name__ == "__main__":
    main()
