import sys
import subprocess
from lines_to_3d_wave import lines_to_3d_wave
sys.path.insert(0, '/home/maldonis/model_analysis/scripts')
from ift_atom_selection import atom_selection

def main():
    if(len(sys.argv) != 6):
        raise Exception("Usage:  modelfile, prefix, jobid, start, stop")
    modelfile = sys.argv[1]
    prefix = sys.argv[2]
    jobid = sys.argv[3]
    start = int(sys.argv[4])
    stop = int(sys.argv[5])
    types = ['mgrid','ft_onespot']

    for i in range(start,stop+1):
        for t in types:
            file = prefix + str(i) + '_' + t + '_' + jobid + '.gfx'
            lines_to_3d_wave(file, prefix + str(i) + '_' + t + '_' + jobid + '.txt')
            if(t == 'mgrid'):
                p = subprocess.Popen(['/home/maldonis/3dft/stdev',file], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                for nextline in iter(p.stdout.readline, ""):
                    sys.stdout.write(nextline)
                    sys.stdout.flush()
                poutput = p.stdout.read()
                perr = p.stderr.read()
                preturncode = p.wait()
                if(preturncode != 0):
                    print("stdev exit status: "+str(preturncode))
                    raise Exception("stdev failed on file {0}!".format(file)+perr)
                atom_selection(modelfile, 'stdev.gfx', prefix + str(i) + '_' + jobid)
                print('')


if __name__ == '__main__':
    main()
