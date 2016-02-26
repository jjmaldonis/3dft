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

    head, _ = os.path.split(outbase)

    pcomm = run_subproc("sbatch -D {head} /home/maldonis/3dft/submits/3dft.sh {modelfile} {outbase} {numpix}".format(
            head=head,
            modelfile=modelfile,
            outbase=outbase,
            numpix=numpix
    ))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

    pcomm = run_subproc("sbatch -D {head} --dependency=afterok:{jobid} /home/maldonis/3dft/submits/3dft_analysis.sh {modelfile} {outbase}ft.gfx {head}".format(
            head=head,
            modelfile=modelfile,
            outbase=outbase,
            jobid=jobid
    ))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

    pcomm = run_subproc("sbatch -D {head} --dependency=afterok:{jobid} /home/maldonis/3dft/submits/ift.sh {paramfile}".format(
            head=head,
            jobid=jobid,
            paramfile=os.path.join(head, 'paramfile.txt')
    ))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])

    pcomm = run_subproc("sbatch -D {head} --dependency=afterok:{jobid} /home/maldonis/3dft/submits/batch_convert.sh {paramfile} {jobid} {numpix}".format(
            head=head,
            jobid=jobid,
            numpix=numpix,
            paramfile=os.path.join(head, 'paramfile.txt')
    ))
    jobid = int(pcomm[0].strip().split('\n').pop().split()[-1])


if __name__ == '__main__':
    main()
