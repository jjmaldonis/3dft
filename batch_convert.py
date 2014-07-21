import sys
import subprocess
from lines_to_3d_wave import lines_to_3d_wave
sys.path.insert(0, '/home/maldonis/model_analysis/scripts')
from ift_atom_selection import atom_selection
from 3d_rot import calc_rot_array_from_hkl, rot
from model import Model

def main():
    if(len(sys.argv) != 3):
        raise Exception("Usage: paramfile, jobid")
    paramfile = sys.argv[1]
    jobid = sys.argv[2]
    types = ['mgrid','ft_onespot']
    with open(paramfile) as f:
        params = f.readlines()
    params = [x.strip() for x in params]
    modelfile = params[0]
    stop = int(params[1])

    for i in range(0,stop):
        j = i*4 +2
        for t in types:
            prefix = params[j] + t + '_' + jobid
            file = prefix + '.gfx'
            lines_to_3d_wave(file, prefix + '.txt')
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
                new_model_file_base_name = params[j] + jobid
                atom_selection(modelfile, 'stdev.gfx', new_model_file_base_name)
                xc = int(params[j+1].split()[2])
                yc = int(params[j+2].split()[2])
                zc = int(params[j+3].split()[2])
                m = Model(new_model_file_base_name+'.xyz')
                rot_arr = calc_rot_array_from_hkl(xc,yc,zc)
                rot(m,rot_arr)
                m.write_real_xyz(new_model_file_base_name+'.rotated.real.xyz')
                print('')


if __name__ == '__main__':
    main()
