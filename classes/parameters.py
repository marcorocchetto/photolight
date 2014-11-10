'''

Exolight 0.1 - Default parameters

'''

import os
import inspect
import logging
import sys

# initilise logging
logging.basicConfig(level=logging.INFO)

class Parameters():

    # specify platform. Available: mac
    platform = 'mac'

    # directory names
    idir = '/Users/marco/Dropbox/workspace/packages-dev/photolight'  # installation dir
    wdir = '/robo/'                                        # working directory
    oecdir = idir + '/extern/oec/systems/'                 # Open Exoplanet Catalogue Database location
    aconfdir = idir + '/astromatic_config/'                # Astromatic configuration files (don't change!)
    vartoolsdir = idir + '/extern/%s/vartools/' % platform # vartools ephemfile, leapsecfile, planetdatafile file location

    # external executable
    jktld = idir + '/extern/%s/jktld/jktld' % platform                 # jktld executable
    vartools = idir + '/extern/%s/vartools/vartools' % platform        # vartools executable
    sex = idir + '/extern/%s/sex' % platform                           # sexextractor executable
    scamp = idir + '/extern/%s/scamp' % platform                       # scamp executable
    missfits = idir + '/extern/%s/missfits' % platform                 # missfits executable

    # Fits keys. @todo: Needs to be moved/improved!
    extlist=['.FIT','.FITS','.FTS','.fts','.fits','.fit']   # recognised FITS file extensions
    platekeys = ['CTYPE1', 'CRVAL1', 'CRPIX1', 'CTYPE2', 'CRVAL2', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2'] # plate solution
    stkeys = ['DATE-OBS', 'EXPTIME', 'OBJCTRA', 'OBJCTDEC', 'EPOCH'] # other keys
    extrakeys = ['EAIRMASS']  # add here other keys

    # Parameters specific to the class Dataset
    dataset = {}
    dataset['pssexlimit'] = 15     # minimum number of stars to proceed with extraction in each frame
    dataset['edgelimit'] = 30      # chip edge limits in pixels. Stars in the edges of the chip will be excluded
    dataset['masternstack'] = 5    # number of frames to stack to create master frame. Might be less if frames are misaligned
    dataset['masternstack_min'] = 2 # minimum number of frames to stack to create master frame
    dataset['maxstars'] = 50       # max number of stars to keep in the master catalogue
    dataset['minstars'] = 2        # minimum number of stars to keep in the master catalogue
    dataset['framemaxoffset'] = 3  # maximum offset allowed for frames in arcminutes.
    dataset['detect_tresh'] = 15       # initial detection threshold for sextractor. Default is 15
    dataset['onlinecat'] = 'GSC-2.3'    # online catalogue to use for plate solution. Default is 'GSC-2.3'

    # Iraf combine params
    iraf_imcombine = {}
    iraf_imcombine['scale'] = 'mode'
    iraf_imcombine['weight'] = 'mode'
    iraf_imcombine['combine'] = 'median'
    iraf_imcombine['Stdout'] = 1

    # Iraf photometry parameters
    irafphot = {}
    irafphot['salgori'] = 'median'
    irafphot['annulus'] = '25'
    irafphot['dannulu'] = '15'
    irafphot['calgori'] = 'centroid'
    irafphot['cbox'] = '10'
    irafphot['maxshif'] = '3'
    irafphot['aperture_range'] = [3,12] # apertures range [min, max] (in px)
    irafphot['aperture_step'] = 0.25       # step between apertures (in px)
    #irafphot['apsizes'] = [1,2,3]      # Use only specified apertures. Must be a list (in px)

    # Switch on/off executions of different modules
    photometry = {}
    photometry['compute_fluxes'] = True    # optimization of ensemble + relative fluxes, for all stars
    photometry['scint_noise_corr'] = True  # rescale the scintillation noise, based on rms
    photometry['sigscalefactor'] = True    # rescale error bars, based on reduced chi^2 distribution of comparison stars
