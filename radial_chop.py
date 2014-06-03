import sys
import numpy as np

def radial_chop(mat3d,xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax):
    xbinsize = float(len(mat3d))/(xmax-xmin)
    ybinsize = float(len(mat3d))/(ymax-ymin)
    zbinsize = float(len(mat3d))/(zmax-zmin)

    xx = np.arange(xmin,xmax,xbinsize,dtype=float)
    yy = np.arange(xmin,ymax,ybinsize,dtype=float)
    zz = np.arange(zmin,zmax,zbinsize,dtype=float)

    for i,yzplane in enumerate(mat3d):
        for j,zlist in enumerate(yzplane):
            for k,z in enumerate(zlist):
                r = sqrt(xx[i]**2 + yy[j]**2 + zz[k]**2)
                if( r => rmin and r =< rmax):
                    print('{0} {1} {2} {3}'.format(xx[i],yy[j],zz[k],z))


def main():
    with open(sys.argv[1]) as f:
        content = f.readlines()
    content = [ [ [float(zl[i]) for i in xrange(len(zl))] for zl in yzl] for yzl in content]

    radial_chop(content,-1.5,1.5,-1.5,1.5,-1.5,1.5,0.34118,0.42353)
