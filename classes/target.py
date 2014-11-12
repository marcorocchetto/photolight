from functions import *

import logging
import sys

class Target(object):
    def __init__(self, pars=None):

        if not pars:
            logging.error('You need to specify a parameter file!')
            sys.exit()

        self.pars = pars

        # ra and dec are used to insert planets that are not in the catalogue
        logging.info('Initialize class Target')

        name = pars.target['name']
        ra = pars.target['ra']
        dec = pars.target['dec']
        if not name or not ra or not dec:
            logging.error('You need to specify the name and coordinates of your object!')
            sys.exit()

        self.ra_d = sex_to_deg(ra, 'ra')
        self.dec_d = sex_to_deg(dec, 'dec')
        self.ra = deg_to_sex(ra=self.ra_d)
        self.dec = deg_to_sex(dec=self.dec_d)
        self.name = name
        logging.info('Set target name to %s' % self.name)
        logging.info('Set target coordinates to %s %s' % (self.ra, self.dec))
        self.mname = self.name.replace(' ', '_')
