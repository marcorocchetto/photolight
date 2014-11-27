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
    parser.add_option('-s', '--source',
                      dest="source_directory",
    )

    options, remainder = parser.parse_args()

    param_filename = os.path.abspath(os.path.expanduser(options.param_filename))
    output_filename = os.path.abspath(os.path.expanduser(options.output_filename))
    source_directory = os.path.abspath(os.path.expanduser(options.source_directory))

    pars = Parameters(param_filename)

    target = Target(pars)
    dataset = Dataset(target, odir=source_directory)
    dataset.create_master_frame()
    photometry = Photometry(dataset, pars)

    pickle.dump(photometry, open(output_filename, 'wb'))
    logging.info('Photometry instance saved in %s' % output_filename)

if __name__ == '__main__':
    main()