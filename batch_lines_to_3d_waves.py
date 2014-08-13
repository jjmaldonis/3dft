import sys
from lines_to_3d_wave import lines_to_3d_wave

def main():
    base = sys.argv[1]
    jobid = sys.argv[2]
    base = [base[0:base.index('spot')+4],base[base.index('_'+jobid):]]
    for i in range(0,20):
        #print(str(i).join(base),str(i).join(base)[:-4]+'.txt')
        lines_to_3d_wave(str(i).join(base),str(i).join(base)[:-4]+'.txt')


if __name__ == '__main__':
    main()
