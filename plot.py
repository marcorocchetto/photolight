
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
from plotting import Plotting

import functions, img_scale, list_frames
from functions import *
from img_scale import *
from list_frames import *

def main():

    parser = optparse.OptionParser()
    parser.add_option('-i', '--in',
                      dest="phot_class",
    )
    parser.add_option('-a', '--aperture',
                      dest="aperture",
    )
    parser.add_option('-o', '--output',
                      dest="out_plot",
    )

    options, remainder = parser.parse_args()

    print 'load ', options.phot_class

    photometry = pickle.load(open(options.phot_class))

    plot = Plotting(photometry, options.out_plot)
    print 'Plot noise fire'

    plot.plot_noise(float(options.aperture))

if __name__ == '__main__':
    main()