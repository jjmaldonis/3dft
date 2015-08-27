"""

This is the python version of the Igor code that runs the spot finding and fitting on the results of 3dft.

Necessary functions to write:
IsoAverage3D
RadiusChop
makeSphere
#makeSphereColorWave
#VolumeStdDev
#StdDev3d
#VolumeSum

FindSpots
gauss_fit_particle
gauss_fit_particles
#regenerateSpotParameters
create_output_wave
create_output_particle
fitSpotTo3DGauss
SaveParamfile
#create_windowed_spots
#create_unwindowed_spots
MDsort
FindFirstMin
gauss3d
#ShowXYZLineProfiles
#ShowXYZLineProfile
#ImageLineProfile3D
#try_automagic_spot_fix
AppendGvecsToVk
linearInterp
cut_ft_smaller_for_isosurface
window_edge_checking


Essentialfunctions:
IsoAverage3D
RadiusChop
gauss_fit_particle
gauss_fit_particles
create_output_wave
create_output_particle
fitSpotTo3DGauss
SaveParamfile
gauss3d

"""

import sys,os
import re
import numpy as np
import matplotlib.pyplot as  plt


def wave3d_to_numpy(igortext_filename):
    k = 0
    i = 0
    for numlines,line in enumerate(open(igortext_filename,'rU').readlines()):
        if(numlines > 2):
            for j,x in enumerate(line.split()):
                #print(k,i,j,x,numlines)
                data[i][j][k] = float(x)
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



# I want this Scale class to include a naming scheme like:
# DimSize, DimDelta, DimOffset, etc.
#class Scale(object):
#    """ An Igor scale-like object. """
#    def __init__(self,scale):
        

class Wave:
    """ An Igor wave-like object.

    Attributes:
        self.dimensions: The number of dimensions of the wave
        self.data: A numpy array containing the data
        self.scale: A tuple of 3-tuples where each 3-tuple contains the start and end scale range for the dimension's
            scale, followed by the step value, for each i-th dimension. The number of 2-tuples is equal to the number of
            dimensions.
            Default scale for each dimension is (0.,1.,...).

    Methods:
        get_scale:
        
    """

    def __init__(self,data,scale=None):
        """ Creates an Igor wave-like object.

        Args:
            data: The raw wave data; any dimension is accepted
            scale: (optional) The scale for each dimension in tuple format (start_value, end_value, step) for each 
                dimension, also in tuple format. If scale is not given, it defaults to (0.,1.,...) for each dimension.

        Returns:
            An Igor wave-like object
        """
        if(type(data) == type([])): data = np.array(data)
        self.dimensions = len(data.shape)
        self.data = data
        self.shape = self.data.shape
        # default to (0.,1.,) if scale is not given
        if(scale == None):
            scale = [(0.,1.,1./(self.data.shape[i]-1)) for i in range(self.dimension)]
        self.scale = scale
        self.full_scale = [self.set_full_scale(d) for d in range(self.dimensions)]

    def __getitem__(self,tup):
        if(isinstance(tup,int) or isinstance(tup,'numpy.int64')):
            return self.data[tup]
        elif(len(tup) == 2):
            x,y = tup
            return self.data[x][y]
        elif(len(tup) == 3):
            x,y,z = tup
            return self.data[x][y][z]
        else:
            raise Exception("__getitem__ for Wave object has not been implemented for waves of this dimension, you must access .data manually")
    def __setitem__(self,tup,val):
        if(isinstance(tup,int) or isinstance(tup,'numpy.int64')):
            self.data[tup] = val
        elif(len(tup) == 2):
            x,y = tup
            self.data[x][y] = val
        elif(len(tup) == 3):
            x,y,z = tup
            self.data[x][y][z] = val
        else:
            raise Exception("__setitem__ for Wave object has not been implemented for waves of this dimension, you must access .data manually")
        

    def get_scale(self,coord,pixelsGiven=True):
        """ Returns the scale of the coordinate in each dimension as a numpy array. """
        if(pixelsGiven):
            scale = np.array([ self.scale[i][0] + self.scale[i][2]*coord[i] for i in range(len(coord)) ])
            #scale = np.array( self.scale[i][0] + self.scale[i][2]*coord[i] for i in range(len(coord)) )
        else:
            # TODO
            scale = None
        return scale
    
    def set_scale(self,scale):
        self.scale = scale

    def set_full_scale(self,dim):
        """ Returns a 1D array of the entire scale of the data for the given dimension from start to finish by step.
            Args:
                dim: dimension to calculate the scale for. Value should be on of (1,2,3,...,ndim)
        """
        return np.array([self.scale[dim][0] + self.scale[dim][2]*i for i in range(self.shape[dim])])

    def x2pnt(self,dim,x):
        scale = self.full_scale[dim]
        idx = (np.abs(scale-x)).argmin()
        return idx


# I should probably move this to be a function of the wave class
# Then I can dynamically write it for multiple dimensions I think
def IsoAverage3D(data,strip_width):
    """ Calculates a 1D array of the average intensity at each radius position.

    Args:
        data: An Igor wave-like object of dimension=3
        strip_width: The number of pixels to average over

    Returns:
        A 1D array of the average intensity at each radius position.

    """

    dim = min( min( data.shape[0],data.shape[1]), data.shape[2] ) / (2*strip_width) + 1
    iso_av = np.zeros( dim, dtype=float)
    npix = np.zeros( dim, dtype=float)
    end_scale = min( min( abs(data.scale[0][0]), abs(data.scale[1][0]) ), abs(data.scale[2][0]) )
    scale = ( (0, end_scale, float(end_scale)/(npix.shape[0]-1)) ,)
    iso_av = Wave(iso_av, scale)
    npix = Wave(npix, scale)

    for i in range(data.shape[0]):
        print(i)
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                r = np.linalg.norm(data.get_scale(tuple([i,j,k])))
                #print(data.scale[0][2])
                #print(data.get_scale(tuple([i,j,k])))
                #print(i,j,k,data.get_scale(tuple([i,j,k])),r,npix.x2pnt(0,r))
                # If r is within the bounds of the data wave = bounds of the iso wave
                if(r < iso_av.scale[0][1]+iso_av.scale[0][2]):
                    #print(r,npix.x2pnt(0,r))
                    npix[npix.x2pnt(0,r)] += 1
                    iso_av[iso_av.x2pnt(0,r)] += data[i][j][k]
                    #print(data[i][j][k])
    print(iso_av.scale)
    print(iso_av.data)
    print(npix.data)
    scale = iso_av.scale
    iso_av = np.divide(iso_av.data,npix.data)
    iso_av = Wave(iso_av,scale)
    print(list(iso_av.data))
    return iso_av



def RadiusChop(data,rmin,rmax):
    """ Chops the data wave by setting everything to 0 except the points radially between rmin and rmax.
    Args:
        data: An Igor wave-like object
        rmin: The starting radius r
        rmax: The ending radius r
    Returns:
        A new Igor wave-like object.
    """
    radius_chop = np.zeros( data.shape, dtype=float)
    for i in range(data.shape[0]):
        print(i)
        for j in range(data.shape[1]):
            for k in range(data.shape[2]):
                r = np.linalg.norm(data.get_scale(tuple([i,j,k])))
                if( (r<=rmax) and (r>=rmin) ):
                    radius_chop[i][j][k] = data[i][j][k]


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
  
    output:
        the smoothed signal
        
    example:
  
    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
  
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
  
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y
  
  

def plot(data):
    plt.plot(data.data)
    plt.show()
  
  

def FindSpots(data, numstdevs):
    """ """
    scale = (-1.5,1.5,3./(data.shape[0]-1))
    data.set_scale(tuple([scale,scale,scale]))
    iso_av = IsoAverage3D(data,2)
    scale = iso_av.scale
    iso_av = smooth(iso_av.data,3)
    iso_av = Wave(iso_av,scale)
    #plot(iso_av)
    min_index  = np.r_[True, iso_av.data[1:] < iso_av.data[:-1]] & np.r_[iso_av.data[:-1] < iso_av.data[1:], True]
    print(min_index)
    try:
        min_index = [i for i,x in enumerate(min_index) if x][0]
    except:
        min_index = 0
    print(min_index)







def main():
    wave = wave3d_to_numpy(sys.argv[1])
    scales = (-1.5,1.5,3./64.)
    wave = Wave(wave, (scales, scales, scales) )
    #IsoAverage3D(wave,2)
    FindSpots(wave,2)
    

if __name__ == '__main__':
    main()
