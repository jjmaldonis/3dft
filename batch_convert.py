import sys
import subprocess
from lines_to_3d_wave import lines_to_3d_wave
sys.path.insert(0, '/home/maldonis/model_analysis/scripts')
from ift_atom_selection import atom_selection
from rot_3d import calc_rot_array_from_hkl, rot
from model import Model
from math import sqrt

def main():
    if(len(sys.argv) != 4):
        raise Exception("Usage: paramfile, jobid, npix for center positions of spots")
    paramfile = sys.argv[1]
    jobid = sys.argv[2]
    npix = int(sys.argv[3])
    #types = ['mgrid','ft_onespot1','ft_onespot2']
    types = ['mgrid']
    with open(paramfile) as f:
        params = f.readlines()
    params = [x.strip() for x in params]
    modelfile = params[0]
    stop = int(params[1])

    for i in range(0,stop):
        j = i*7 +2
        for t in types:
            prefix = params[j] + t + '_' + jobid
            file = prefix + '.gfx'
            lines_to_3d_wave(file, prefix + '.txt')
            if(t == 'mgrid'):
                p = subprocess.Popen(['/home/maldonis/3dft/stdev',file,jobid], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
                worked = atom_selection(modelfile, 'stdev_'+jobid+'.gfx', 256, new_model_file_base_name)
                x0,y0,z0 = tuple([float(x) for x in params[j+1].split()])
                sx,sy,sz = tuple([float(x) for x in params[j+2].split()])
                cxy,cxz,cyz = tuple([float(x) for x in params[j+3].split()])
                xc = int(params[j+4].split()[2])
                yc = int(params[j+5].split()[2])
                zc = int(params[j+6].split()[2])
                gvec = float(params[j+6].split()[3])
                print("g-vector length = {0}".format(gvec))
                print("The more correct g-vector length == {0}".format(sqrt((npix/4.0-x0+1)**2 + (npix/4.0-y0+1)**2 + (npix/4.0-z0)**2+1)*3.0/(npix/2)))
                if(worked == 0):
                    m = Model(new_model_file_base_name+'.xyz')
                    rot_arr = calc_rot_array_from_hkl(npix/2-xc,npix/2-yc,npix/2-zc)
                    rot(m,rot_arr)
                    m.write_real_xyz(new_model_file_base_name+'.rotated.real.xyz')
                print('')


if __name__ == '__main__':
    main()
