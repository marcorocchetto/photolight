
import os
import sys
import logging
import numpy as np
from StringIO import StringIO
from ConfigParser import SafeConfigParser


class Parameters():

    def __init__(self, parfile):

        #config file parser
        self.parser = SafeConfigParser()
        self.parser.read(parfile)

        self.default_parser = SafeConfigParser()
        self.default_parser.read('default.par')

        self.verbose = self.getpar('General', 'verbose', 'bool')
        self.wdir = os.path.abspath(os.path.expanduser(self.getpar('General', 'wdir')))

        logging.basicConfig(filename=os.path.abspath(os.path.expanduser(self.getpar('General', 'log'))), level=logging.INFO)

        # define a Handler which writes INFO messages or higher to the sys.stderr
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

        try:
            self.odir = self.getpar('General', 'odir')
        except:
            self.odir = None

        # executables
        self.vartools = self.getpar('Executables', 'vartools')
        self.sex = self.getpar('Executables', 'sex')

        # some directories
        self.idir = os.path.dirname(os.path.realpath(__file__))
        self.vartoolsdir = self.idir + '/../vartools'
        self.aconfdir = self.idir + '/../astromatic_config'

        # Target
        try:
            self.target = {'name': self.parser.get('Target', 'name'),
                           'ra': self.parser.get('Target', 'ra'),
                           'dec': self.parser.get('Target', 'dec')}
        except:
            logging.error('Target parameters could not be found. Include name, ra, dec.')
            sys.exit()

        # Observatory
        try:
            self.observatory = {'name':      self.parser.get('Observatory', 'name'),
                                'elevation': self.parser.getfloat('Observatory', 'elevation'),
                                'longitude': self.parser.getfloat('Observatory', 'longitude'),
                                'latitude':  self.parser.getfloat('Observatory', 'latitude')}
        except:
            logging.error('Observatory parameters could not be found. Include name, elevation, longitude, latitude.')
            sys.exit()

        # Telescope
        try:
            self.telescope = {'name':     self.parser.get('Telescope', 'name'),
                              'egain':    self.parser.getfloat('Telescope', 'egain'),
                              'rdnoise':  self.parser.getfloat('Telescope', 'rdnoise'),
                              'aperture': self.parser.getfloat('Telescope', 'aperture')}
        except:
            logging.error('Observatory parameters could not be found. Include name, egain, rdnoise, aperture.')
            sys.exit()

        # Filter
        try:
            self.filter = {'name':   self.parser.get('Filter', 'name'),
                           'midwav': self.parser.getfloat('Filter', 'midwav')}
        except:
            logging.error('Filter parameters could not be found. Include name, midwav.')
            sys.exit()

        # Dataset parameters
        self.dataset = {'run_all_sex':            self.getpar('Dataset', 'run_all_sex', 'bool'),
                        'edge_limit':             self.getpar('Dataset', 'edge_limit', 'float'),
                        'master_nstack':          self.getpar('Dataset', 'master_nstack', 'int'),
                        'master_max_stars':       self.getpar('Dataset', 'master_max_stars', 'int'),
                        'frame_max_offset':       self.getpar('Dataset', 'frame_max_offset', 'float'),
                        'detect_tresh':           self.getpar('Dataset', 'detect_tresh', 'int'),
                        'ignore_unsolved_frames': self.getpar('Dataset', 'ignore_unsolved_frames', 'bool'),
                        'min_frames':             self.getpar('Dataset', 'min_frames', 'int')}

        self.header_egain     = self.getpar('HeadersKeys','egain', 'list-str')
        self.header_dateobs   = self.getpar('HeadersKeys','dateobs', 'list-str')
        self.header_imagetyp  = self.getpar('HeadersKeys','imagetyp', 'list-str')
        self.header_filter    = self.getpar('HeadersKeys','filter', 'list-str')
        self.header_ccdtemp    = self.getpar('HeadersKeys','ccdtemp', 'list-str')
        self.header_ccdtemp_units    = self.getpar('HeadersKeys','ccdtemp_units')
        self.header_object   = self.getpar('HeadersKeys','object', 'list-str')
        self.header_objctra   = self.getpar('HeadersKeys','objctra', 'list-str')
        self.header_objctdec  = self.getpar('HeadersKeys','objctdec', 'list-str')
        self.header_naxis1  = self.getpar('HeadersKeys','naxis1', 'list-str')
        self.header_naxis2  = self.getpar('HeadersKeys', 'naxis2', 'list-str')


        self.platesolve = { 'platesolve': self.getpar('PlateSolve','platesolve', 'bool'),
                            'overwrite': self.getpar('PlateSolve','overwrite', 'bool'),
                            'solvefield_path': self.getpar('PlateSolve','solvefield_path'),
                            'data_dir': self.getpar('PlateSolve','data_dir'),
                            'use_coordinates': self.getpar('PlateSolve','use_coordinates', 'bool'),
                            'search_radius': self.getpar('PlateSolve','search_radius', 'float'),
                            'flags': self.getpar('PlateSolve','flags')}

        # Iraf combine params
        self.iraf_imcombine = {'scale':   self.getpar('IRAF combine', 'scale'),
                               'weight':  self.getpar('IRAF combine', 'weight'),
                               'combine': self.getpar('IRAF combine', 'combine')}

        self.iraf_phot = {'irafcall': self.getpar('IRAF photometry', 'irafcall', 'bool'),
                          'salgori':  self.getpar('IRAF photometry', 'salgori'),
                          'annulus':  self.getpar('IRAF photometry', 'annulus'),
                          'dannulu':  self.getpar('IRAF photometry', 'dannulu'),
                          'calgori':  self.getpar('IRAF photometry', 'calgori'),
                          'cbox':     self.getpar('IRAF photometry', 'cbox'),
                          'maxshif':  self.getpar('IRAF photometry', 'maxshif')}

        self.photometry = {'apertures':  self.getpar('Photometry', 'apertures', 'list-float'),
                           'multithread': self.getpar('Photometry', 'multithread', 'bool'),
                           'tool': self.getpar('Photometry', 'tool')}

        # Switch on/off executions of different modules
        self.modules_run = {'compute_fluxes': self.getpar('General', 'compute_fluxes', 'bool'),       # optimization of ensemble + relative fluxes, for all stars
                            'bad_frames_exclusion': self.getpar('General', 'bad_frames_exclusion', 'bool'),   # rescale the scintillation noise, based on rms
                            'scint_noise_corr': self.getpar('General', 'scint_noise_corr', 'bool'),   # rescale the scintillation noise, based on rms
                            'rescale_errorbars': self.getpar('General', 'rescale_errorbars', 'bool')} # rescale error bars, based on red chi^2 distribution of comp. stars



    def getpar(self, sec, par, type=None):

        # get parameter from user defined parser. If parameter is not found there, load the default parameter
        # the default parameter file parser is self.default_parser, defined in init

        try:

            if type == None:
                try:
                    return self.parser.get(sec, par)
                except:
                    return self.default_parser.get(sec, par)
            elif type == 'float':
                try:
                    return self.parser.getfloat(sec, par)
                except Exception,e:
                    return self.default_parser.getfloat(sec, par)

            elif type == 'bool':
                try:
                    return self.parser.getboolean(sec, par)
                except:
                    return self.default_parser.getboolean(sec, par)
            elif type == 'int':
                try:
                    return self.parser.getint(sec, par)
                except:
                    return self.default_parser.getint(sec, par)
            elif type == 'list-str':
                try:
                    l = self.parser.get(sec,par).split(',')
                    return [str(m).strip() for m in l]
                except:
                    l = self.default_parser.get(sec,par).split(',')
                    return [str(m).strip() for m in l]
            elif type == 'list-float':
                try:
                    l = self.parser.get(sec,par).split(',')
                    return [float(m) for m in l]
                except:
                    l = self.default_parser.get(sec,par).split(',')
                    return [float(m) for m in l]
            elif type == 'list-int':
                try:
                    l = self.parser.get(sec,par).split(',')
                    return [int(m) for m in l]
                except:
                    l = self.default_parser.get(sec,par).split(',')
                    return [int(m) for m in l]
            else:
                logging.error('Cannot set parameter %s in section %s. Parameter type %s not recognized. Set to None' (par, sec, type))
                return None
        except:
            logging.error('Cannot set parameter %s in section %s. Set to None' % (par, sec))
            return None