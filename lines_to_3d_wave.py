import sys


def lines_to_3d_wave(infile,outfile):
    """ takes a 3d matrix (of equal dimensions) written out on a single line
        (from fortran) and writes it to a file of format Igor Text Wave """

    of = open(outfile,'w')
    with open(infile) as f:
        print("Opened infile {0} for reading...".format(infile))
        for i,line in enumerate(f):
            line = line.strip().split()
            if(i == 0):
                npixx, npixy, npixz = tuple([int(x) for x in line])
                print("Number of pixels in each direction: {0} {1} {2}".format(npixx, npixy, npixz))
                of.write('IGOR\nWAVES/N=({0},{1},{2})\t {3}\nBEGIN\n'.format(npixx,npixy,npixz,outfile[:outfile.index('.')]))
            else:
                if(len(line) != npixx*npixy):
                    print("Wrong number of elements on line {0}! You have {1} and there should be {2}.".format(i,len(line),npixx*npixy))
                of.write(' '.join(line)+'\n')
    of.write('END\n')
    of.write('X SetScale/P x 0,1,"", {0}; SetScale/P y 0,1,"", {0}; SetScale/P z 0,1,"", {0}; SetScale d 0,0,"", {0}'.format(outfile[:outfile.index('.')]))
    of.close()

def main():
    infile = sys.argv[1]
    outfile = sys.argv[2]
    lines_to_3d_wave(infile,outfile)

if __name__ == "__main__":
    main()
