import sys
import subprocess
from lines_to_3d_wave import lines_to_3d_wave
import shlex

def main():
    if(len(sys.argv) != 4):
        raise Exception("Usage: paramfile, jobid, npix")
    paramfile = sys.argv[1]
    jobid = sys.argv[2]
    npix = int(sys.argv[3])
    with open(paramfile) as f:
        params = f.readlines()
    params = [x.strip() for x in params]
    stop = int(params[1])

    for i in range(0,stop):
        j = i*7 +2
        prefix = params[j] + jobid
        modelfile = 'submodels/' + prefix + '.xyz'
        print('qsub ft_slurm.sh {0} {1} {2}'.format(modelfile,prefix+'_512_',npix))
        #args = shlex.split('qsub ft_slurm.sh {0} {1} {2}'.format(modelfile,prefix+'_512_',npix))
        print(' '.join(['qsub ../ft_slurm.sh',modelfile,prefix+'_512_',str(npix)]))
        p = subprocess.Popen(['qsub','../ft_slurm.sh',modelfile,prefix+'_512_',str(npix)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #print(' '.join(['../3dft',modelfile,prefix+'_512_',str(npix)]))
        #p = subprocess.Popen(['../3dft',modelfile,prefix+'_512_',str(npix)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        #p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        for nextline in iter(p.stdout.readline, ""):
            sys.stdout.write(nextline)
            sys.stdout.flush()
        poutput = p.stdout.read()
        perr = p.stderr.read()
        preturncode = p.wait()
        if(preturncode != 0):
            print("qsub ft_slurm.sh exit status: "+str(preturncode))
            raise Exception("qsub ft_slurm.sh failed: {0}!".format(perr))
        #print('')


if __name__ == '__main__':
    main()
