import sys, os, time
import subprocess
import shlex


def run_subproc(args):
    """ Run subprocess given args as a command line string"""
    print(args)
    #return 5000
    args = shlex.split(args)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pcomm = p.communicate()
    poutput = pcomm[0] #p.stdout.read()
    perr = pcomm[1] #p.stderr.read()
    print(poutput)
    print(perr)
    preturncode = p.wait()
    if(preturncode != 0):
        raise Exception("{0} failed!".format(args[0]))
    return pcomm

def main():
    modelfile = sys.argv[1]
    outbase = sys.argv[2]
    try:
        numpix = sys.argv[3]
    except:
        numpix = 512

    pcomm = run_subproc("sbatch /home/maldonis/3dft/submits/3dft.sh {0} {1} {2}".format(modelfile, outbase, numpix))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

    pcomm = run_subproc("sbatch --dependency=afterok:{2} /home/maldonis/3dft/submits/3dft_analysis.sh {0} {1}ft.gfx".format(modelfile, outbase, jobid))
    #pcomm = run_subproc("sbatch /home/maldonis/3dft/submits/3dft_analysis.sh {0} {1}ft.gfx".format(modelfile, outbase))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

    pcomm = run_subproc("sbatch --dependency=afterok:{0} /home/maldonis/3dft/submits/ift.sh paramfile.txt".format(jobid))
    #pcomm = run_subproc("sbatch /home/maldonis/3dft/submits/ift.sh paramfile.txt")
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

    pcomm = run_subproc("sbatch --dependency=afterok:{0} /home/maldonis/3dft/submits/batch_convert.sh paramfile.txt {0} {1}".format(jobid, numpix))
    #pcomm = run_subproc("sbatch /home/maldonis/3dft/submits/batch_convert.sh paramfile.txt {0} {1}".format(jobid, numpix))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

if __name__ == '__main__':
    main()
