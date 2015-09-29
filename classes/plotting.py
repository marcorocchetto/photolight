import os
import shutil
import math
import sys
import numpy as np
import logging
from robuststatistics import *
from scipy import stats
from functions import *
import matplotlib.pylab as plt


np.set_printoptions(threshold=np.nan)

class Plotting(object):

    def __init__(self, photometry, plot_out='', pars=None):

        logging.info('Initialize class Plotting')

        self.p = photometry
        self.d = photometry.dataset
        self.t = photometry.target
        self.plot_out = plot_out

    def plot_noise(self, apsize):

        logging.info('Plotting noise budget')

        p = self.p.aperture[apsize]

        plt.scatter(p.medmag, p.rms)
        plt.plot(p.mags_model, p.sigmastar_model, label='Star noise')
        plt.plot(p.mags_model, p.sigmasky_model, label='Sky noise')
        plt.plot(p.mags_model, p.sigmaread_model, label='Readout noise')
        plt.plot(p.mags_model, p.sigscint_model, label='Scint. noise')
        plt.yscale('log')
        plt.plot(p.mags_model, p.sigtot_model, label='Tot. noise')
        plt.legend(loc=4)

        plt.savefig(os.path.join(self.plot_out, 'noise_budget.pdf'))
        plt.show()

    def plot_all_flux(self, apsize, xval='bjd', maxstars=None, savefolder=None):

        if xval == 'bjd':
            x = self.d.bjd
            xlabel = 'BJD'
        elif xval == 'hjd':
            x = self.d.bjd
            xlabel = 'HJD'
        elif xval == 'airmass':
            x = self.d.bjd
            xlabel = 'Airmass'

        if not maxstars or max > self.d.nstar:
            maxstars = self.d.nstar

        for nstar in range(maxstars):
            if nstar == self.d.targetid:
                starlabel = self.t.name
            else:
                starlabel = 'n. %i' % nstar

            if nstar in self.p.aperture[apsize].star_exclude:
                starlabel += ' (EXCLUDED)'

            plt.figure()
            plt.title('%s - Star %s - Aperture %.1f px' % (self.t.name, starlabel, apsize))
            plt.plot(x, self.p.aperture[apsize].eflux[nstar], 'o')
            plt.xlabel(xlabel)
            plt.ylabel('Relative flux')
            if savefolder:
                plt.savefig(os.path.join(savefolder, '%s_eflux_star%i.pdf' % (xval, nstar)))
            else:
                plt.show()


    def plot_star_flux(self, nstar=None, apsize=None, xval='bjd'):

        if not nstar:
            nstar = self.d.targetid
        if not apsize:
            apsize = self.p.apsizes[0]

        if xval == 'bjd':
            x = self.d.bjd
            xlabel = 'BJD'
        elif xval == 'hjd':
            x = self.d.bjd
            xlabel = 'HJD'
        elif xval == 'airmass':
            x = self.d.bjd
            xlabel = 'Airmass'

        plt.figure()
        plt.title('%s - Star n %i - Aperture %.1f px' % (self.t.name, nstar, apsize))
        plt.plot(x, self.p.aperture[apsize].eflux[nstar], 'o')
        plt.xlabel(xlabel)
        plt.ylabel('Relative flux')
        plt.show()
