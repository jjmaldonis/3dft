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

import sys, os, copy
import re
import numpy as np
from scipy import ndimage
from numba import jit
import matplotlib.pyplot as  plt
import scipy.signal as signal
import scipy.stats as stats
import scipy.ndimage as ndimage
import scipy.optimize as optimize
from collections import namedtuple


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


def wave3d_to_numpy(igortext_filename):
    """ Takes in an .itx filename and loads the data into a numpy array, which is returned """
    k = 0
    i = 0
    for numlines,line in enumerate(open(igortext_filename,'rU').readlines()):
        if(numlines > 2):
            print(numlines)
            for j,x in enumerate(line.split()):
                j = j % data.shape[1]
                i = int(i/data.shape[1])
                data[i][j][k] = float(x)
            if(i == dims[1]-1):
                k += 1
                i = -1
                if(k == dims[0]):
                    break
            i += 1
            print(i,j,k)
        elif(numlines == 1):
            line = line.strip()
            dims = re.findall(r'([0-9]*[,)])', line)
            dims = tuple(int(x[:-1]) for x in dims)
            data = np.zeros(dims,dtype=float) # Create the array of the correct size
    return data



class Scale1D(object):
    def __init__(self, start, stop=None, step=None, func=None, custom=False):
        self.func = func
        if func is None:
            self.func = lambda x: x
        if custom is False:
            self.start = custom[0]
            self.end = custom[-1]
            self.step = None
        else:
            self.start = self.__start = start
            self.stop = self.__stop = stop
            self.step = self.__step = step

    @property
    def values(self):
        if custom is not False:
            return custom
        else:
            return np.array([self.func(x) for x in drange(self.start, self.stop, self.step)])


class Scale(object):
    """ An Igor scale-like object.

        parameters:

            tup:     For 1D data: (start, stop, step, func=None)
                     For multidimensional data:
                        ( (xstart, xstop, xstep, func=None), (ystart, ystop, ystep, func=None), ...)
                     The optional 'func' inside each scale tuple will be evaluated on each intermediate
                     value in range(start, stop, step) to allow for a non-linear scale. The default
                     None will simply evaulate f(x)=x on each of these range values.
                     Default is None, but this should always be overwritten unless 'custom' is provided.

            custom:  If your data does not have a functionally evaluatable domain, you can define a
                     custom scale object that defines each axis in full.
                     For example, custom=( [0, 2, 4, 8, 16], [0, 3, 9, 27] ) defines a custom log
                     scale for a dataset of shape 5x4.

            Note:    You must define one and only one of 'tup' or 'custom'.
        
        usage:
        
            self.start:  If the data is 1D, this is the value of 'start' in the tuple.
                         Otherwise this is a tuple of the start values given for each dimension.
            self.stop:   Similar to self.start
            self.step:   Similar to self.start
            self.scale:  Returns the fully enumerated scale in each dimension.
            self[...]:   You can get the domain of your data easily by simply passing in the indices
                         at which you want the domain.
    """
    def __init__(self, tup=None, custom=False):
        if custom is False and tup is None:
            raise Exception("You must define a scale somehow!")
        elif custom is not False and tup is not None:
            raise Exception("You cannot define both a custom scale and a functionally-evaluatable one.")
        self.custom = custom # False or the custom full scale you want to return

        if tup is not None and isinstance(tup[0], tuple): # Multidimensional
            self._scale = tuple([Scale1D(*scale)] for scale in tup)
            self.dimensions = len(self._scale)
            self.__meshgrid = np.meshgrid(*self.scale)
            self.start = tuple( [scale.start for scale in self._scale] )
            self.stop = tuple( [scale.stop for scale in self._scale] )
            self.step = tuple( [scale.step for scale in self._scale] )
        else: # 1D
            self._scale = Scale1D(*tup)
            self.dimensions = 1
            self.__meshgrid = self._scale
            self.start = self._scale.start
            self.stop = self._scale.stop
            self.step = self._scale.step

    @property
    def scale(self):
        if self.custom != False: return self.custom
        if axis is not None or self.dimensions == 1:
            return self._scale.values
        else:
            return np.array( [scale.values for scale in self._scale] )

    def __getitem__(self, coord):
        return self.__meshgrid[coord]



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

    def val2index(self, val):
        """ Converts a value in the range to a value in the domain. No interpolation is done, it just finds the nearest piont. """
        idx = (np.abs(self.data-val)).argmin()
        if self.dimensions > 1: idx = np.unravel_index(idx, self.data.shape)
        return idx

    def val2x(self, val):
        """ Converts a value in the range to its nearest index. """
        idx = self.val2index(val)
        return self.full_scale[idx] # This might be a bit buggy for higher dimensional arrays (ie > 1)

    def save(self, outfile):
        """ writes this wave to a file of one of these formats: Igor itx wave; numpy array"""
        filename, extension = os.path.splitext(outfile)
        if extension == '.itx':
            self._save_itx(outfile)
        elif extension == '.npy':
            np.save(outfile, self.data)
        else:
            raise Exception("Cannot save array, filetype not supported. Provide a supported file extension.")
    
    def _save_itx(self, outfile):
        """ Save as Igor wave """
        # Need to reverse the order of the data for Igor
        swapped = self.data
        for i in list(reversed(range(self.dimensions)))[1:]:
            swapped = np.swapaxes(swapped, i-1, i)
        # Now we save the array with the changed viewpoint
        with open(outfile, 'w') as of:
            shape = self.shape
            if len(shape) == 1:
                shape = '({0})'.format(shape[0]) # Need to remove the extra , at the end of the tuple
            of.write('IGOR\nWAVES/N={0}\t {1}\nBEGIN\n'.format(shape, outfile[:outfile.index('.')]))
            self._save_data(swapped, of)
            of.write('END\n')
            # We can even set the scale correctly
            if self.dimensions < 4:
                axis_names = ['x', 'y', 'z', 'd']
                s = 'X ' + ' '.join( ['SetScale/P {0} 0,{1},"", {2};'.format(axis_names[i], self.scale[i][2], outfile[:outfile.index('.')]) for i in range(len(self.scale))] )
                of.write(s)

    def _save_data(self, data, of):
        """ Goes through the data and saves it regardless of dimension. This is a recursive function so it needs to be its own function. """
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
    print("  Transforming to spherical coordinates...")
    sph = to_spherical(wave)
    print("  Doing radial summation for IsoAv...")
    iso_av = np.array( [np.sum(sph[:,i,:])/(sph.shape[1]+sph.shape[2]) for i in range(sph.shape[0])] ) # Not sure why I need [:,i,:] rather than [i,:,:]
    end_scale = min( min( abs(wave.scale[0][0]), abs(wave.scale[1][0]) ), abs(wave.scale[2][0]) )
    scale = ( (0, end_scale, float(end_scale)/(iso_av.shape[0]-1)) ,)
    iso_av = Wave(iso_av,scale)
    print("  Finished IsoAv!")
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


@jit
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
    data = np.fromfunction(lambda i,j,k: within_radius(within_radius,i,j,k), shape=wave.shape, dtype=np.int32)
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


def FindSpots(wave, numstdevs, iso_av=None, verbose=True):
    """ """
    scale = (-1.5,1.5,3./(wave.shape[0]-1))
    wave.scale = (scale,scale,scale,)

    # Get first dip in iso_av for RadiusChop
    if iso_av is None: # Optionally included so that I can pre-calculate it and then reuse it later if I want to versus it being discarded after use in this function. This also allows me to use condense_wave on it first (or some other operation) if I want to for some reason.
        iso_av = IsoAverage3D(wave)
    peaks = FindPeaks(iso_av)
    highest = max(peaks)
    highest_index = iso_av.val2index(highest)
    artifact_cutoff = iso_av.val2x( min(iso_av[0:highest_index]) )
    if verbose: print("Data below R = {0} is bad due to artifacts from the 3D FT, cutting that out for further calculations...".format(artifact_cutoff))
    chopped = RadiusChop(wave, artifact_cutoff, 1.5)

    if verbose: print("Thresholding/masking chopped image...")
    stdev = np.std(chopped.data)
    mean = np.mean(chopped.data)
    if verbose: print("  Mean: {0}".format(mean))
    above = mean + numstdevs*stdev
    threshold = (chopped.data > above).astype(np.int) # This is a 0/1 mask used to find spots in the next section

    if verbose: print("Finding separated spots... (Spots with small volumes will be excluded. Diagonals count as connected.)") # TODO
    labeled_array, num_features = ndimage.measurements.label(threshold, structure=[
        [[1,1,1], [1,1,1], [1,1,1]],
        [[1,1,1], [1,1,1], [1,1,1]],
        [[1,1,1], [1,1,1], [1,1,1]]
        ]) # Allow diagonals as connections -- TODO: maybe not?
    slices = ndimage.measurements.find_objects(labeled_array)
    spots = []
    # Create the spots!
    for slice in slices:
        volume = np.sum(threshold[slice])
        if volume >= 9:
            spot = {
                'pixels': threshold[slice], #TODO Print out a larger slice
                'data': wave.data[slice],
                'slices': slice,
                'volume': volume,
                'summed_intensity': np.sum(wave.data[slice]),
                'max_intensity': np.amax(wave.data[slice]),
                }
            spot['index_center'] = np.unravel_index(np.argmax(wave.data[slice]), spot['data'].shape)
            spot['index_center'] = (spot['index_center'][0] + slice[0].start, spot['index_center'][1] + slice[1].start, spot['index_center'][2] + slice[2].start)
            spot['domain_center'] = (wave.full_scale[0][spot['index_center'][0]], wave.full_scale[1][spot['index_center'][1]], wave.full_scale[2][spot['index_center'][2]]) # TODO this is ugly, I should implement something inside the Wave object to make it much prettier
            spot['|g|'] = np.linalg.norm( spot['domain_center'] )
            spots.append(spot)
    spots.sort(key=lambda s: s['max_intensity'], reverse=True)
    if verbose:
        for i,spot in enumerate(spots):
            print("Spot #{0}".format(i))
            print("  |g| = {0}".format(spot['|g|']))
            print("  max intensity = {0}".format(spot['max_intensity']))
            print("  position = {0}".format(spot['domain_center']))

    if verbose: print("Calculating average background intensity...")
    for slice in slices:
        chopped[slice] = 0.0
    mean = np.mean(chopped.data)
    if verbose: print("  Mean: {0}".format(mean))

    # Identify spot-pairs -- TODO for now I will just remove every-other one like I was doing before...
    # assert len(spots) % 2 == 0
    #spot_pairs = 
    spots = spots[0::2]
    return spots


def expand_spot(wave, spot, extra):
    assert extra % 2 == 0
    new = copy.copy(spot)
    slicex = slice(spot['slices'][0].start-extra/2, spot['slices'][0].stop+extra/2)
    slicey = slice(spot['slices'][1].start-extra/2, spot['slices'][1].stop+extra/2)
    slicez = slice(spot['slices'][2].start-extra/2, spot['slices'][2].stop+extra/2)
    new['slices'] = (slicex, slicey, slicez)
    new['data'] = wave[new['slices']]
    new['index_center'] = (new['index_center'][0]+extra/2, new['index_center'][1]+extra/2, new['index_center'][2]+extra/2)
    new['domain_center'] = (wave.full_scale[0][new['index_center'][0]], wave.full_scale[1][new['index_center'][1]], wave.full_scale[2][new['index_center'][2]]) # TODO this is ugly, I should implement something inside the Wave object to make it much prettier
    new['pixels'] = None
    new['summed_intensity'] = np.sum(new['data'])
    return new


def fit_spot(wave, spot, function='gauss', hold_sigmas=False, hold_centers=False):
    """ There might be a bug somewhere because the best fit is always when extra=0. I would expect that to be the best but it shouldn't always be the case. """
    best = {'spot':None, 'popt':None, 'perr':None, 'residual':1e6, 'extra':0}
    for extra in range(0, 10, 2):
        new = expand_spot(wave, spot, extra)
        try:
            tup = fit_spot_gauss(wave, new)
        except Error as error:
            print(error)
            continue
        if tup is not None:
            popt, perr, res = tup
            print(extra, res)
            if res < best['residual']:
                best['spot'] = new
                best['popt'] = popt
                best['perr'] = perr
                best['residual'] = res
                best['extra'] = extra
        else:
            continue
    return best


def fit_spot_gauss(wave, spot, hold_sigmas=False, hold_centers=False):
    # Normalize the image to have a max intensity of 1, but don't overwrite the original data
    image = np.copy(spot['data'])
    image /= spot['max_intensity']

    # Do the fit
    if np.sum(image.size) > 9:
        # Set widths
        xw, yw, zw = image.shape
        # The spot we fit is going to be the small one so change the center indicies to reflect that
        x0, y0, z0 = [spot['index_center'][i] - spot['slices'][i].start for i in range(len(spot['slices']))]

        # Set initial guesses for parameters (in order of parameters specified in the fitting function)
        params = [xw/2.0, yw/2.0, zw/2.0, 0.1, 0.1, 0.1, x0/2.0, y0/2.0, z0/2.0]
        # Setup the data in a format that can be read and vectorized properly for curve_fit
        data4D = np.array( [ [i,j,k,image[i,j,k]] for i in range(image.shape[0]) for j in range(image.shape[1]) for k in range(image.shape[2]) ] )
        try:
            popt, pcov = optimize.curve_fit(VectorizedGaussian3D, data4D[:,:3], data4D[:,3], p0=params, maxfev=4000)
            perr = np.sqrt(np.diag(pcov))
            res = residuals(image, Gaussian3D, popt)
            return popt, perr, res
        except RuntimeError:
            pass
    return None


def residuals(data, func, parameters):
    return sum([abs(func(i, *parameters)-x) for i,x in np.ndenumerate(data)])/data.size


def Gaussian3D(coord, sx,sy,sz, cxy,cxz,cyz, x0,y0,z0):
    x,y,z = coord
    return np.exp(  -1/(2* (-1+cxy**2+cxz**2+cyz**2-2*cxy*cxz*cyz)) * ( (cyz**2-1)*(x-x0)**2/sx**2 + (cxz**2-1)*(y-y0)**2/sy**2 + (cxy**2-1)*(z-z0)**2/sz**2 + (2*cxy*(x-x0)*(y-y0)-2*cxz*cyz*(x-x0)*(y-y0))/(sx*sy) + (2*cxz*(x-x0)*(z-z0)-2*cxy*cyz*(x-x0)*(z-z0))/(sx*sz) + (2*cyz*(y-y0)*(z-z0)-2*cxy*cxz*(y-y0)*(z-z0))/(sy*sz) )  )


def VectorizedGaussian3D(coord, sx,sy,sz, cxy,cxz,cyz, x0,y0,z0):
    coord = coord.T
    x,y,z = coord
    val = np.exp(  -1/(2* (-1+cxy**2+cxz**2+cyz**2-2*cxy*cxz*cyz)) * ( (cyz**2-1)*(x-x0)**2/sx**2 + (cxz**2-1)*(y-y0)**2/sy**2 + (cxy**2-1)*(z-z0)**2/sz**2 + (2*cxy*(x-x0)*(y-y0)-2*cxz*cyz*(x-x0)*(y-y0))/(sx*sy) + (2*cxz*(x-x0)*(z-z0)-2*cxy*cyz*(x-x0)*(z-z0))/(sx*sz) + (2*cyz*(y-y0)*(z-z0)-2*cxy*cxz*(y-y0)*(z-z0))/(sy*sz) )  )
    return val.flatten()


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



def create_output(modelfilename, spots, filename, num_spots=None):
    if num_spots is None: num_spots = len(spots)
    of = open(filename, 'w')
    of.write('{0}\n'.format(modelfilename))
    of.write('{0}\n'.format(num_spots))
    for i,spot in enumerate(spots):
        of.write('spot{0}_\n'.format(i))
        of.write('{0} {1} {2} {3}\n'.format(spot['fit'].x0,  spot['fit'].y0,  spot['fit'].z0,  spot['|g|']))
        of.write('{0} {1} {2} {3}\n'.format(spot['fit'].sx,  spot['fit'].sy,  spot['fit'].sz,  spot['|g|']))
        of.write('{0} {1} {2} {3}\n'.format(spot['fit'].cxy, spot['fit'].cxz, spot['fit'].cyz, spot['|g|']))
    of.close()


def load_data(infile):
    """ read a file into a numpy array. these formats are supported: Igor wave (.itx); numpy array (.npy)"""
    filename, extension = os.path.splitext(infile)
    if extension == '.itx':
        return wave3d_to_numpy(infile)
    elif extension == '.npy':
        return np.load(infile)
    else:
        raise Exception("Cannot load data, filetype not supported.")



def main():
    modelfile = sys.argv[1]
    print("Loading data...")
    wave = load_data(sys.argv[2])
    scales = (-1.5,1.5,3./64.)
    wave = Wave(wave, (scales, scales, scales) )
    print("Loaded data!")

    print("Calculating annular average...")
    iso_av = IsoAverage3D(wave)
    iso_av.save('iso_av.itx')
    iso_av_stripped = condense_wave(iso_av, 4)

    spots = FindSpots(wave, 10, iso_av=iso_av, verbose=True)

    Fit = namedtuple('Fit', ['x0', 'y0', 'z0', 'sx', 'sy', 'sz', 'cxy', 'cxz', 'cyz', 'residual', 'perr'])
    for i,spot in enumerate(spots):
        Wave(spot['data']).save('spot{0}.npy'.format(i))
        Wave(spot['data']).save('spot{0}.itx'.format(i))
        print("Starting new spot:  {0}".format(i))
        print("  center = {0}".format(spot['index_center']))
        print("  max intensity = {0}".format(spot['max_intensity']))
        print("  |g| = {0}".format(spot['|g|']))
        print("  volume of spot = {0}".format(spot['volume']))
        print("  shape of spot = {0}".format(spot['data'].shape))
        best_fit = fit_spot(wave, spot)
        popt = best_file['popt']
        fit = Fit(*popt, residual=best_fit['residual'], perr=best_fit['perr'])
        spot['fit'] = fit
        if best_fit is not None:
            print(best_fit['extra'])
            print(best_fit['popt'])
            print(best_fit['perr'])
            print(best_fit['residual'])
        else:
            print("  Fit failed :(")
        if i > 20: break

    create_output(modelfile, spots, 'paramfile.txt')


if __name__ == '__main__':
    main()
