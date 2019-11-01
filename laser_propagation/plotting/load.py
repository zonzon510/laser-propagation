"""Module for loading the output data from a simulation"""
from os.path import join, abspath, dirname, isfile
import configparser
import glob
import numpy
from scipy.constants import epsilon_0, c, pi
from scipy.integrate import trapz


class Results:
    """Provides convenient access to all of the output data files"""
    def __init__(self, filename):
        self.filename = filename
        self.config = load_config(filename)
        self.folder = self.config['folder']

        self.distances = load(join(self.folder, self.config['results/distance']))
        self.time = load(join(self.folder, self.config['results/time']))
        self.radius = load(join(self.folder, self.config['results/radius']))
        self.omega = load(join(self.folder, self.config['results/omega']))
        self.wavelength = load(join(self.folder, self.config['results/wavelength']))
        self.kperp = load(join(self.folder, self.config['results/kperp']))


    def electric_field(self, requested_distance):
        """Returns the closest distance and electric field"""
        i = self._find_nearest(requested_distance)
        filename = self._make_filename(self.config['results/temporal_field'], i)
        E = load_binary(filename, len(self.radius), len(self.time))
        z = self.distances[i]
        return z, E


    def electron_density(self, requested_distance):
        """Returns the closest distance and electron density"""
        i = self._find_nearest(requested_distance)
        filename = self._make_filename(self.config['results/electron_density'], i)
        Rho = load_binary(filename, len(self.radius), len(self.time), dtype=float)
        z = self.distances[i]
        return z, Rho


    def spectral_field(self, requested_distance):
        """Return the closest distance and spectral field"""
        i = self._find_nearest(requested_distance)
        filename = self._make_filename(self.config['results/spectral_field'], i)
        A = load_binary(filename, len(self.kperp), len(self.omega))
        z = self.distances[i]
        return z, A


    def spectrum(self, requested_distance):
        """Returns the closest distance and spectrum"""
        i = self._find_nearest(requested_distance)
        filename = self._make_filename(self.config['results/spectral_field'], i)
        A = load_binary(filename, len(self.kperp), len(self.omega))
        z = self.distances[i]


        I = 1/2*epsilon_0*c*abs(A)**2
        k = self.kperp

        # To integrate in kperp space, we use cylindrical rings where
        # each ring is centered at 1/2 (k[i] + k[i-1])
        N = len(k)
        weights = numpy.zeros(N)

        if N == 1:
            # if only one point in k, the weight is a cylinder with radius = 3/2 k
            weights[0] = pi * (1.5 * k[0])**2
        else:
            weights[0] = pi * ((k[1] + k[0]) / 2)**2
            for i in range(1, N-1):
                weights[i] = pi * ((k[i+1] + k[i]) / 2)**2 - weights[i-1]

            # k-point beyond last one is found by extrapolation
            weights[-1] = pi * ((3*k[-1] - k[-2]) / 2)**2 - weights[-2]

        spec = weights.dot(I)
        return z, spec
    
    
    def energy(self):
        """Returns the distances and energy values along propagation"""
        z, energies = load(join(self.folder, self.config['results/energy']))
        return z, energies


    def max_intensity(self):
        """Returns the distances and maximum intensity values along propagation"""
        z, maxI = load(join(self.folder, self.config['results/max_intensity']))
        return z, maxI

    
    def max_electron_density(self):
        """Returns the distances and maximum electron density values along propagation"""
        z, max_den = load(join(self.folder, self.config['results/max_density']))
        return z, max_den
    
    
    def _make_filename(self, name, i):
        return join(self.folder, '{}{:03d}.dat'.format(name, i))


    def _find_nearest(self, requested_distance):
        i = abs(requested_distance - self.distances).argmin()
        return i


    
def convert(string):
    """Attempts to convert string into an int or float. If possible
returns the converted value, otherwise returns the original string.

    """
    try:
        value = int(string)
        return value
    except ValueError:
        pass

    try:
        value = float(string)
        return value
    except ValueError:
        return string

    
def load_config(filename):
    """Load the parameters from a configuration file that has the syntax
[section]
key = value

and return them as a dictionary in the form
{'section/key': value}

"""
    if not isfile(filename):
        raise OSError('Config file not found: {}'.format(filename))
    cfg = configparser.ConfigParser()
    cfg.read(filename)
    d = {}
    for section in cfg.sections():
        for key in cfg[section]:
            k = '/'.join((section, key))
            v = cfg[section][key]
            d[k] = convert(v)

    # add some extra parameters
    d['folder'] = dirname(abspath(filename))
    return d
    

def load(filename):
    """Reads plain text multi-column file and returns each column as its own array"""
    return numpy.loadtxt(filename, unpack=True, ndmin=1)


def load_complex(filename):
    """Loads data from a plain text file where columns alternate between
real and imaginary

    """
    a = numpy.loadtxt(filename, dtype=float, ndmin=1)
    a.dtype = complex
    return a


def load_binary(filename, nrows, ncols, dtype=complex):
    """Load data from a binary file"""
    temporal = numpy.fromfile(filename, dtype=dtype)
    temporal.shape = (nrows, ncols)
    return temporal



        
