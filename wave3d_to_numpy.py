import sys
import re
import numpy as np


def wave3d_to_numpy(igortext_filename):
    k = 0
    i = 0
    for numlines,line in enumerate(open(igortext_filename,'rU').readlines()):
        if(numlines > 2):
            for j,x in enumerate(line.split()):
                #print(k,i,j,x,numlines)
                data[k][i][j] = float(x)
            if(i == dims[1]-1):
                k += 1
                i = -1
                if(k == dims[0]):
                    break
            i += 1
        elif(numlines == 1):
            line = line.strip()
            dims = re.findall(r'([0-9]*[,)])', line)
            dims = tuple(int(x[:-1]) for x in dims)
            data = np.zeros(dims,dtype=float)
    return data


def main():
    print(wave3d_to_numpy(sys.argv[1]))

if __name__ == '__main__':
    main()

