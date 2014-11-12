#! /usr/bin/python -W ignore

import sys
import optparse
import logging
import pickle
import os

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
    parser.add_option('-o', '--out',
                      dest="output_filename",
    )
    parser.add_option('-i', '--in',
                      dest="dataset_instance",
    )

    options, remainder = parser.parse_args()

    pars = Parameters(options.param_filename)

    try:
        dataset = pickle.load(open(options.dataset_instance))
        logging.info('Dataset instance correctly initialised')
    except:
        logging.error('Cannot load dataset instance. Check argument -i')
        sys.exit()

    logging.info('Target is %s (RA %s Dec %s), star n. %i in master catalogue' %
                 (dataset.object.name, dataset.object.ra, dataset.object.dec, dataset.targetid))

    photometry = Photometry(dataset, pars)

    pickle.dump(photometry, open(options.output_filename, 'wb'))
    logging.info('Photometry instance saved in %s' % options.output_filename)


if __name__ == '__main__':
    main()