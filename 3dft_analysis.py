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
from scipy import ndimage
from numba import jit
import matplotlib.pyplot as  plt
import scipy.signal as signal


def wave3d_to_numpy(igortext_filename):
    """ Takes in an .itx filename and loads the data into a numpy array, which is returned """
    k = 0
    i = 0
    for numlines,line in enumerate(open(igortext_filename,'rU').readlines()):
        if(numlines > 2):
            for j,x in enumerate(line.split()):
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
            data = np.zeros(dims,dtype=float) # Create the array of the correct size
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
            scale, followed by the step value, for each i-th dimension. The number of 3-tuples is equal to the number of
            dimensions.
            Default scale for each dimension is (0., npix, 1.)
        self.shape: The numpy shape of the data

    Methods:
        get_scale:
        
    """

    def __init__(self, data, scale=None, full_scale=None):
        """ Creates an Igor wave-like object.

        Args:
            data: The raw wave data; any dimension is accepted
            scale: (optional) The scale for each dimension in tuple format (start_value, end_value, step) for each 
                dimension, also in tuple format. If scale is not given, it defaults to (0., npix, 1.) for each dimension.

        Returns:
            An Igor wave-like object
        """
        if isinstance(data, list): data = np.array(data)
        self.dimensions = len(data.shape)
        self.data = data
        self.shape = self.data.shape
        if scale is None: # Default to (0., npix, 1.) if scale is not given
            scale = [(0., (self.data.shape[i]-1), 1.) for i in range(self.dimensions)]
        self.scale = scale
        self._full_scale = full_scale # This will only be used to define a custom (non-linear) scale. In this case the user must also set 'scale' to 'Custom'

    # If I get rid of the axis keyword and reduce the functionality then I would make this a property
    def calculate_full_scale(self, axis=None):
        """ Returns a 1D array of the entire scale of the data for the given dimension from start to finish by step.
            Args:
                dim: dimension to calculate the scale for. Value should be on of (1,2,3,...,ndim)
        """
        if self.scale == 'Custom': return self._full_scale
        if axis is not None or self.dimensions == 1:
            if self.dimensions == 1:
                dim = 0
            else:
                dim = axis
            return np.array([self.scale[dim][0] + self.scale[dim][2]*i for i in range(self.shape[dim])])
        else:
            #return [np.array([self.scale[dim][0] + self.scale[dim][2]*i for i in range(self.shape[dim])]) for dim in range(self.dimensions)]
            return np.array( [self.calculate_full_scale(dim) for dim in range(self.dimensions)] )

    @property
    def full_scale(self):
        if self.scale != 'Custom':
            return self.calculate_full_scale()
        else:
            return self._full_scale

    def __getitem__(self, tup):
        return self.data[tup]

    def __setitem__(self, tup, val):
        self.data[tup] = val

    # I have no idea what the point of this is but I use it in RadiusChop
    def get_scale(self, coord, pixelsGiven=True):
        """ Returns the scale of the coordinate in each dimension as a numpy array. """
        if(pixelsGiven):
            scale = np.array([ self.scale[i][0] + self.scale[i][2]*coord[i] for i in range(len(coord)) ])
            #scale = np.array( self.scale[i][0] + self.scale[i][2]*coord[i] for i in range(len(coord)) )
        else:
            # TODO
            scale = None
        return scale
    
    def x2index(self, x, dim=0):
        """ Converts a point in the domain to its nearest index. """
        scale = self.calculate_full_scale()[dim]
        idx = (np.abs(scale-x)).argmin()
        return idx

    def x2val(self, x, dim=0):
        """ Converts a value in the domain to a value in the range. No interpolation is done, it just finds the nearest point. """
        return self[self.x2index(x, dim)]

    def __call__(self, x, dim=0):
        """ Shorthand for x2val(x, dim) """
        return self.x2val(x, dim)

    def val2index(self, val, dim=0):
        """ Converts a value in the range to a value in the domain. No interpolation is done, it just finds the nearest piont. """
        idx = (np.abs(self.data-val)).argmin()
        return idx

    def val2x(self, val, dim=0):
        """ Converts a value in the range to its nearest index. """
        idx = self.val2index(val, dim)
        return self.full_scale[idx] # This might be a bit buggy for higher dimensional arrays (ie > 1)

    def save(self, outfile, ext="itx"):
        """ writes this wave to a file of format Igor Text Wave """
        # Need to reverse the order of the data for Igor
        swapped = self.data
        for i in list(reversed(range(self.dimensions)))[1:]:
            swapped = np.swapaxes(swapped, i-1, i)
        # Now we save the array with the changed viewpoint
        with open(outfile, 'w') as of:
            shape = self.shape
            if len(shape) == 1: shape = '({0})'.format(shape[0]) # Need to remove the extra , at the end of the tuple
            of.write('IGOR\nWAVES/N={0}\t {1}\nBEGIN\n'.format(self.shape, outfile[:outfile.index('.')]))
            self._save_data(swapped, of)
            #for arr in swapped:
            #    for row in arr:
            #        row = [str(x) for x in row]
            #        of.write('\t'.join(row)+'\n')
            of.write('END\n')
            # We can even set the scale correctly
            if self.dimensions < 4:
                axis_names = ['x', 'y', 'z', 'd']
                s = 'X ' + ' '.join( ['SetScale/P {0} 0,{1},"", {2};'.format(axis_names[i], self.scale[i][2], outfile[:outfile.index('.')]) for i in range(len(self.scale))] )
                #of.write('X SetScale/P x 0,{0},"", {3}; SetScale/P y 0,{1},"", {3}; SetScale/P z 0,{2},"", {3};'.format(self.scale[0][2], self.scale[1][2], self.scale[2][2], outfile[:outfile.index('.')]))
                of.write(s)

    def _save_data(self, data, of):
        if type(data[0]) in [int, float, np.float64, np.int64, np.int32, np.float32]:
            of.write('\t'.join([str(x) for x in data])+'\n')
        else:
            for row in data:
                self._save_data(row, of)


# I should probably move this to be a function of the wave class
# Then I can dynamically write it for multiple dimensions I think
#@jit
def IsoAverage3D(wave):
    """ Calculates a 1D array of the average intensity at each radius position.

    Args:
        wave: An Igor wave-like object of dimension=3

    Returns:
        A 1D array of the average intensity at each radius position.

    """
    # Convert to spherical coordinates for fast summation
    sph = to_spherical(wave)
    iso_av = np.array( [np.sum(sph[:,i,:])/(sph.shape[1]+sph.shape[2]) for i in range(sph.shape[0])] ) # Not sure why I need [:,i,:] rather than [i,:,:]
    #for x in iso_av:
    #    print(x)
    end_scale = min( min( abs(wave.scale[0][0]), abs(wave.scale[1][0]) ), abs(wave.scale[2][0]) )
    scale = ( (0, end_scale, float(end_scale)/(iso_av.shape[0]-1)) ,)
    iso_av = Wave(iso_av,scale)
    return iso_av


#@jit
def condense_wave(wave, width):
    """ Condenses a 1D wave to be shorter by averaging its neighboring elements.
        Setting the scale doesn't work great right now. """
    # Apply width
    new = np.array( [sum(wave[i:i+width])/(width) for i in range(0, wave.shape[0], width)] )
    end_scale = min( [abs(wave.scale[i][1]) for i in range(wave.dimensions)] )
    scale = ( (0, end_scale, float(end_scale)/(new.shape[0]-1)) ,)
    return Wave(new, scale)


#@jit
def to_spherical(wave):
    scale = np.array( wave.full_scale )
    rmax = max( abs(np.amin(scale[0,:])), abs(np.amax(scale[0,:])) ) # rmax should really take into account the diagonal distance, so this essentially cuts off the corners
    sph_scale = Wave(np.zeros(scale.shape), scale)
    R = np.linspace(0.0, rmax, sph_scale.shape[1])
    THETA = np.linspace(0, np.pi, sph_scale.shape[1])
    PHI = np.linspace(-np.pi, np.pi, sph_scale.shape[1])
    sgrid = np.meshgrid(R, THETA, PHI)
    #sgrid = np.array( np.meshgrid(R, THETA, PHI) )
    cgrid = np.zeros_like(sgrid)
    cgrid[0] = sgrid[0]*np.sin(sgrid[2])*np.cos(sgrid[1])
    cgrid[1] = sgrid[0]*np.sin(sgrid[2])*np.sin(sgrid[1])
    cgrid[2] = sgrid[0]*np.cos(sgrid[2])
    npix = scale.shape[-1]-2
    m = npix/(np.amax(scale)-np.amin(scale))
    b = -m*np.amin(scale)
    c2grid = m*cgrid + b
    sph = ndimage.map_coordinates(wave.data, c2grid, cval=0.0)
    sph = Wave(sph, scale=sph_scale)
    #print(np.amin(wave.data))
    #print(np.amax(wave.data))
    #print(np.amin(c2grid))
    #print(np.amax(c2grid))
    #print(np.amin(sph.data))
    #print(np.amax(sph.data))
    return sph


class _WithinRadius(object):
    """ This is a helper class to make RadiusChop run much faster """
    def __init__(self, wave, rmin, rmax):
        self.wave = wave
        self.full_scale = wave.full_scale
        self.rmin2 = rmin**2
        self.rmax2 = rmax**2

    @np.vectorize
    def __call__(self, i,j,k):
        return self.wave[i,j,k] if self.rmin2 < self.full_scale[0][i]**2+self.full_scale[1][j]**2+self.full_scale[2][k]**2 < self.rmax2 else 0.0
        

def RadiusChop(wave, rmin, rmax):
    """ Chops the data wave by setting everything to 0 except the points radially between rmin and rmax.
    Args:
        data: A Wave object
        rmin: The starting radius r
        rmax: The ending radius r
    Returns:
        A new Igor wave-like object.
    """
    within_radius = _WithinRadius(wave, rmin, rmax)
    data = np.fromfunction(lambda i,j,k: within_radius(klass,i,j,k), shape=wave.shape, dtype=np.int32)
    return Wave(data, wave.scale)


def smooth(x, window_len=11, window='hanning'):
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


def FindSpots(wave, numstdevs):
    """ """
    scale = (-1.5,1.5,3./(wave.shape[0]-1))
    wave.scale = (scale,scale,scale,)

    # Get first dip in iso_av for RadiusChop
    iso_av = IsoAverage3D(wave)
    peaks = FindPeaks(iso_av)
    highest = max(peaks)
    highest_index = iso_av.val2index(highest)
    artifact_cutoff = iso_av.val2x( min(iso_av[0:highest_index]) )
    print("Data below R = {0} is bad due to artifacts from the 3D FT, cutting that out...".format(artifact_cutoff))
    chopped = RadiusChop(wave, artifact_cutoff, 1.5)
    chopped.save('output.itx')


def FindPeaks(wave):
    maxs = signal.argrelextrema(wave.data, np.greater)[0]
    positions = [wave.full_scale[x] for x in maxs]
    intensities = [wave[x] for x in maxs]
    return Wave(intensities, scale='Custom', full_scale=positions)


def FindValleys(wave):
    maxs = signal.argrelextrema(wave.data, np.less)[0]
    positions = [wave.full_scale[x] for x in maxs]
    intensities = [wave[x] for x in maxs]
    return Wave(intensities, scale='Custom', full_scale=positions)



def main():
    wave = wave3d_to_numpy(sys.argv[1])
    scales = (-1.5,1.5,3./64.)
    wave = Wave(wave, (scales, scales, scales) )
    iso_av = IsoAverage3D(wave)
    iso_av_stripped = condense_wave(iso_av, 4)
    #print("Full scale")
    #print([list(x) for x in iso_av.full_scale])
    #print("Data")
    #print(list(iso_av.data))
    FindSpots(wave, 15)

if __name__ == '__main__':
    main()
