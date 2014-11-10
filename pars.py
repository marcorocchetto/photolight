'''
ETAS Parameter File
'''

# Parameters
params = {}
params['pssexlimit'] = 15     # minimum number of stars to proceed with extraction in each frame
params['edgelimit'] = 30      # chip edge limits in pixels. Stars in the edges of the chip will be excluded
params['masternstack'] = 5    # number of frames to stack to create master frame. Might be less if frames are misaligned
params['masternstack_min'] = 2 # minimum number of frames to stack to create master frame
params['maxstars'] = 50       # max number of stars to keep in the master catalogue
params['minstars'] = 2        # minimum number of stars to keep in the master catalogue
params['framemaxoffset'] = 3  # maximum offset allowed for frames in arcminutes.
params['detect_tresh'] = 15       # initial detection threshold for sextractor. Default is 15
params['onlinecat'] = 'GSC-2.3'    # online catalogue to use for plate solution. Default is 'GSC-2.3'
params['aperture_step'] = 0.25    # minimum step between different apertures
params['aperture_range'] = [2, 15]    # factor for range of apertures: avg_fwhm +/- aperture_range*avg_fwhm

# Iraf combine params
iraf_imcombine = {}
iraf_imcombine['scale'] = 'mode'
iraf_imcombine['weight'] = 'mode'
iraf_imcombine['combine'] = 'median'
iraf_imcombine['Stdout'] = 1

# Iraf photometry params
irafphot = {}
irafphot['salgori'] = 'median'
irafphot['annulus'] = '25'
irafphot['dannulu'] = '15'
irafphot['calgori'] = 'centroid'
irafphot['cbox'] = '10'
irafphot['maxshif'] = '3'

# Switch on/off executions of different modules
photometry = {}
photometry['compute_fluxes'] = True # optimization of ensemble + relative fluxes, for all stars
photometry['scint_noise_corr'] = True  # rescale the scintillation noise, based on rms
photometry['sigscalefactor'] = True  # rescale error bars, based on reduced chi^2 distribution of comparison stars


# Observatories
class observatory:
    name = ''
    fullname = ''
    elevation = 0
    longitud = 0
    latitud = 0
    timezone = 0

observatories = {}
# here you can add new observatories

#observatories['ulo'] = observatory()
#observatories['ulo'].name = 'ULO'
#observatories['ulo'].fullname = 'University of London Observatory'
#observatories['ulo'].elevation = 82.0
#observatories['ulo'].latitud = 51.6133
#observatories['ulo'].longitud = 0.24
#observatories['ulo'].timezone = 0

## Telescopes
class telescope:
    name = ''
    fullname = ''
    aperture = 0
    rdnoise= 0        # readout noise
    egain =  0        # electron gain
    aperture = 0


    pmat = {} # approximate matrix for plate solution
    pmat['CD1_1'] = ''
    pmat['CD1_2'] = ''
    pmat['CD2_1'] = ''
    pmat['CD2_2'] = ''

telescopes = {}
# here you can add new telescopes

#telescopes['c14w'] = telescope()
#telescopes['c14w'].name = 'c14w'
#telescopes['c14w'].fullname = 'C14 West'
#telescopes['c14w'].rdnoise = 10
#telescopes['c14w'].egain = 2.53
#telescopes['c14w'].aperture = 35
#telescopes['c14w'].pmat['CD1_1'] = ' 2.37780000000E-004'
#telescopes['c14w'].pmat['CD1_2'] = ' 6.73650000000E-006'
#telescopes['c14w'].pmat['CD2_1'] = '-6.73040000000E-006'
#telescopes['c14w'].pmat['CD2_2'] =  '2.37999000000E-004'

# Filters
class filter:
    name = ''
    midwav = 0

filters = {}
#here you can add new filters

#filters['rc'] = filter()
#filters['rc'].name = 'Rc'
#filters['rc'].midwav = 650


# This file has been optimiesed for use with Mac OS X

import os, inspect, logging, sys

# initilise logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# specify platform
# Currently available: mac
platform = 'mac'

# diractory names
idir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))  # installation dir
wdir = '/robo/'                                        # working directory
oecloc = idir + '/extern/oec/systems/'                 # Open Exoplanet Catalogue Database location
aconfdir = idir + '/astromatic_config/'                # Astromatic configuration files (don't change!)
vartoolsdir = idir + '/extern/%s/vartools/' % platform # vartools ephemfile, leapsecfile, planetdatafile file location

# external executable
jktld = idir + '/extern/%s/jktld/jktld' % platform                 # jktld executable
vartools = idir + '/extern/%s/vartools/vartools' % platform        # vartools executable
sex = idir + '/extern/%s/sex' % platform                           # sexextractor executable
scamp = idir + '/extern/%s/scamp' % platform                       # scamp executable
missfits = idir + '/extern/%s/missfits' % platform                 # missfits executable

# Misc TO BE IMPROVED!
extlist=['.FIT','.FITS','.FTS','.fts','.fits','.fit']   # recognised FITS file extensions
platekeys = ['CTYPE1', 'CRVAL1', 'CRPIX1', 'CTYPE2', 'CRVAL2', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'] # plate solution
stkeys = ['DATE-OBS', 'EXPTIME', 'OBJCTRA', 'OBJCTDEC', 'EPOCH'] # other keys
extrakeys = ['EAIRMASS']  # add here other keys
