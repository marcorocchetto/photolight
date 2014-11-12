#! /usr/bin/python -W ignore

import sys
import optparse
import logging
import pickle
import os

curdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(curdir, 'classes'))
sys.path.append(os.path.join(curdir, 'library'))

import parameters, target, dataset
from parameters import Parameters
from target import Target
from dataset import Dataset

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
    parser.add_option('-s', '--source',
                      dest="source_directory",
    )

    options, remainder = parser.parse_args()

    pars = Parameters(options.param_filename)

    target = Target(pars)
    dataset = Dataset(target, odir=options.source_directory)
    dataset.create_master_frame()

    pickle.dump(dataset, open(options.output_filename, 'wb'))
    logging.info('Dataset instance saved in %s' % options.output_filename)

    logging.info('Master frame fits, png and catalogue files are saved in %s' %
                 os.path.join(dataset.pars.wdir, dataset.masterdir))


if __name__ == '__main__':
    main()