from ..library.functions import *


import os
import logging
import sys
import imp
import oecpy

from parameters import *

# start logging
etaslog = logging.getLogger('etaslog')

class Target(object):
    def __init__(self, name, pars=None, ra='', dec=''):

        # ra and dec are used to insert planets that are not in the catalogue
        etaslog.info('Initialize class Target')

        if not pars:
            self.pars = Parameters()
        else:
            self.pars = pars

        if ra != '' and dec != '':
            self.ra_d = sex_to_deg(ra, 'ra')
            self.dec_d = sex_to_deg(dec, 'dec')
            self.ra = deg_to_sex(ra=self.ra_d)
            self.dec = deg_to_sex(dec=self.dec_d)
            self.name = name
            etaslog.info('Set target name to %s' % (self.name))
            etaslog.info('Set target coordinates to %s %s' % (self.ra, self.dec))
            self.mname = self.name.replace(' ', '_')

class Planet(Target):


    def __init__(self, name, ra='', dec=''):

        Target.__init__(self, ra, dec)

        etaslog.info('Initialize subclass Planet')

        # load OEC database of exo-planets
        etaslog.info('Loading parameters from Open Exoplanet Catalogue')
        exocat = oecpy.OECDatabase(self.pars.oecdir)
        self.oecplanet = exocat.searchPlanet(name)

        if not self.oecplanet:
            etaslog.info('The planet was not found in the catalogue')
            etaslog.info('Use ra and dec from Target class (input arguments)')
            if not (self.ra or self.dec):
                etaslog.error('RA and/or Dec have not been specified for this object. As the planet %s is not'
                              'recognized in the Open Exoplanet Catalogue you need to provide RA and Dec in the '
                              'arguments of your class instance. '
                              'Use e.g.: planet = photo.Planet("PlanetName", ra="10:12:45.3", dec="-45:17:50")')
        else:
            self.oecplanet.ra_d = sex_to_deg(self.oecplanet.ra, 'ra')
            self.oecplanet.dec_d = sex_to_deg(self.oecplanet.dec, 'dec')
            self.name = self.oecplanet.name
            self.ra_d = self.oecplanet.ra_d
            self.dec_d = self.oecplanet.dec_d
            self.dec_d = self.oecplanet.dec_d
            self.ra = self.oecplanet.ra
            self.dec = self.oecplanet.dec

        self.mname = self.name.replace(' ', '_')

        directory = os.path.join(self.pars.wdir, self.mname)
        if not os.path.isdir(directory):
            try:
                etaslog.info('Creating directory for the planet %s: %s' % (self.name, directory))
                os.mkdir(directory)
            except Exception, e:
                etaslog.error('Failed to create the folder %s. Cannot continue' % directory, exc_info=True)
                sys.exit()
        else:
             etaslog.info('The directory %s already exists' % directory)

        etaslog.info('Planet instance correctly initialised')