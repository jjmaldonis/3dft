import os
import sys
import shlex
import subprocess


def run_subproc(args):
    """ Run subprocess given args as a command line string"""
    print(args)
    #return 5000
    args = shlex.split(args)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pcomm = p.communicate()
    poutput = pcomm[0] #p.stdout.read()
    perr = pcomm[1] #p.stderr.read()
    print(poutput.rstrip())
    print(perr.rstrip())
    preturncode = p.wait()
    if(preturncode != 0):
        raise Exception("{0} failed!".format(args[0]))
    return pcomm


def main():
    dir = sys.argv[1]
    for root, subFolders, files in os.walk(dir):
        for f in files:
            if '.xyz' not in f:
                continue
            outbase = f.split('.')
            outbase = [b for b in outbase if 'model_final' in b][0]
            outbase = os.path.abspath(os.path.join(root, outbase))
            f = os.path.abspath(os.path.join(root, f))
            command = 'python /home/maldonis/3dft/submits/submit.py {} {}'.format(f, outbase)
            run_subproc(command)

        
if __name__ == '__main__':
    main()
