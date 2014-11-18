import sys
import optparse
import logging
import pickle
import os
import numpy as np

curdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(curdir, 'classes'))
sys.path.append(os.path.join(curdir, 'library'))

import parameters, target, dataset, photometry
from parameters import Parameters
from dataset import Dataset
from photometry import Photometry

import functions, img_scale, list_frames
from functions import *
from img_scale import *
from list_frames import *

def main():

    parser = optparse.OptionParser()
    parser.add_option('-p', '--parfile',
                      dest="param_filename",
    )
    parser.add_option('-i', '--in',
                      dest="photometry_instance",
    )
    parser.add_option('-a', '--aperture',
                      dest="aperture_size",
    )
    parser.add_option('-o', '--out',
                      dest="output_filename",
    )
    parser.add_option('-d', '--data',
                      dest="print_data",
    )

    options, remainder = parser.parse_args()

    param_filename = os.path.abspath(os.path.expanduser(options.param_filename))
    photometry_instance = os.path.abspath(os.path.expanduser(options.photometry_instance))
    output_filename = os.path.abspath(os.path.expanduser(options.output_filename))


    pars = Parameters(param_filename)

    try:
        photometry = pickle.load(open(photometry_instance))
        logging.info('Photometry instance correctly initialised')
    except:
        logging.error('Cannot load photometry instance. Check argument -i')
        sys.exit()

    logging.info('Target is %s (RA %s Dec %s), star n. %i in master catalogue' %
                 (photometry.dataset.object.name, photometry.dataset.object.ra, photometry.dataset.object.dec,
                  photometry.dataset.targetid))

    if options.print_data == 'mag':
        bjd = np.asarray(photometry.dataset.bjd)
        mag = np.nan_to_num(photometry.aperture[float(options.aperture_size)].mag)
        #np.savetxt(output_filename, np.column_stack((bjd,mag)), fmt='%.7e')

        f = open(output_filename, 'wb')
        for n in range(photometry.nobs):
            line = '%.10e' % bjd[n]
            for i in range(3):
                line += ' %.7e' % mag[n, i]
            f.write(line + '\n')
        f.close()

    if options.print_data == 'mag_err':
        bjd = np.asarray(photometry.dataset.bjd)
        mag = np.nan_to_num(photometry.aperture[float(options.aperture_size)].mag[:,:3])
        sigmag = np.nan_to_num(photometry.aperture[float(options.aperture_size)].sigmag[:,:3])

        f = open(output_filename, 'wb')
        for n in range(photometry.nobs):
            line = '%.10e' % bjd[n]
            for i in range(3):
                line += ' %.7e %.7e' % (mag[n, i], sigmag[n, i])
            f.write(line + '\n')
        f.close()


if __name__ == '__main__':
    main()