import os
import shutil
import math
import glob
import alipy
import subprocess
import sys
import numpy as np
import logging
import matplotlib.pyplot as plt

import imp

from robuststatistics import *

from multiprocessing import Process, Queue
from scipy import stats
from pyraf import iraf

from ..library.functions import *
from ..library.list_frames import *

np.set_printoptions(threshold=np.nan)

etaslog = logging.getLogger('etaslog')

class Photometry(object):

    def __init__(self, dataset, apsizes=[], irafcall=True, multithread=True):

        etaslog.info('Initialize class Photometry')

        # inherit parameters class from dataset
        self.pars = dataset.pars

        self.dataset = dataset
        self.targetid = dataset.targetid
        self.object = dataset.object
        self.nobs = dataset.nobs
        self.nstar = dataset.nstar

        # initialize variables used to store photometry with different apertures
        self.xc = {}
        self.yc = {}
        self.skystd = {}
        self.psky = {}
        self.fluxsum = {}
        self.parea = {}
        self.flux = {}

        # set iraf photometry working directory
        self.irafphotdir = os.path.join(self.dataset.telescopedir, 'irafphot')
        directory = os.path.join(self.pars.wdir, self.irafphotdir)
        if not os.path.isdir(directory):
            try:
                etaslog.info('Creating directory %s' % directory)
                os.mkdir(directory)
            except Exception, e:
                etaslog.error('Failed to create the folder %s. Cannot continue'
                              % directory, exc_info=True)
                sys.exit()
        else:
             etaslog.info('The directory %s already exists' % directory)

        ######################################
        # Initialize and run IRAF photometry #
        ######################################

        if self.dataset.frames:

            # create lists for iraf photometry
            f2 = open(os.path.join(self.pars.wdir, self.irafphotdir, 'outputlist'), 'w')
            f3 = open(os.path.join(self.pars.wdir, self.irafphotdir, 'coordlist'), 'w')
            for n in self.dataset.frames:
                basename = os.path.splitext(os.path.basename(self.dataset.frames[n]['path']))[0]
                shutil.copy(os.path.join(self.pars.wdir, self.dataset.masterdir, 'master.final'),
                            os.path.join(self.pars.wdir, self.irafphotdir, basename + '.coord'))
                f2.write(os.path.join(self.pars.wdir, self.irafphotdir, basename + '.phot') + '\n')
                f3.write(os.path.join(self.pars.wdir, self.irafphotdir, basename + '.coord') + '\n')
            f2.close()
            f3.close()

            if isinstance(apsizes, list):
                self.apsizes = apsizes
            else:
                # determine aperture sizes from default values
                if 'apsizes' in self.pars.irafphot:
                    if self.pars.irafphot['apsizes'] is list:
                        apsizes = self.pars.irafphot['apsizes']
                else:
                    apsizes = np.arange(self.pars.irafphot['aperture_range'][0],
                                        self.pars.irafphot['aperture_range'][1],
                                        self.pars.irafphot['aperture_step'])


            if self.apsizes is not list or len(self.apsizes) == 0:
                etaslog.error('The aperture sizes you specified are not valid')

            apsizestr = ', '.join([str(i) for i in apsizes])
            etaslog.info('Photometry will be performed using the following apertures: %s' % apsizestr)


            # run iraf photometry
            self.run_iraf(apsizes, irafcall)

            # multithreading: if enabled, each aperture is executed as a separate thread.
            # the aperture tasks are run using an instance of the class Aperture, which is then saved
            # into the current photometry class as self.aperture[apsize]

            if not multithread: # do not use multithreads

                self.aperture = {}
                for apsize in apsizes:
                    # initialise instance of Aperture
                    aperture = Aperture(self, apsize)
                    if self.pars.photometry['compute_fluxes']:
                        aperture.compute_fluxes()
                    del aperture.p
                    self.aperture[apsize] = aperture

                    # print out logs
                    for log in aperture.aplog:
                        etaslog.info(log)
                    del aperture.aplog
            else:

                # run each aperture on a separate core
                gotQueues = dict()

                def empty_queues(jobs):
                    for i, job in enumerate(jobs):
                        if not queues[i].empty():
                            if i in gotQueues:
                                gotQueues[i] += queues[i].get()
                            else:
                                gotQueues[i] = queues[i].get()

                queues = [Queue() for apsize in apsizes]
                jobs = [MultiThread(apsize, self, queues[i]) for i, apsize in enumerate(apsizes)]
                for job in jobs:
                    job.start()

                import time
                while any([jj.is_alive() for jj in jobs]):
                    time.sleep(0.01)  # 0.01 Wait a while before next update. Slow down updates for really long runs.
                    empty_queues(jobs)

                # store apertures into current photometry class
                self.aperture = {}
                for i, apsize in enumerate(apsizes):
                    self.aperture[apsize] = gotQueues[i]

                    # print out logs
                    for log in self.aperture[apsize].aplog:
                        etaslog.info(log)
                    del self.aperture[apsize].aplog

    def run_iraf(self, apsizes, irafcall=True):

        apsizestr = ', '.join([str(i) for i in apsizes])

        if irafcall:

            etaslog.info('Perform IRAF photometry')

            for n in self.dataset.frames:
                basename = os.path.splitext(os.path.basename(self.dataset.frames[n]['path']))[0]
                # remove phot files
                if os.path.isfile(os.path.join(self.pars.wdir, self.irafphotdir, basename + '.phot')):
                    os.remove(os.path.join(self.pars.wdir, self.irafphotdir, basename + '.phot'))
                if os.path.isfile(os.path.join(self.pars.wdir, self.irafphotdir, basename + '.phot2')):
                    os.remove(os.path.join(self.pars.wdir, self.irafphotdir, basename + '.phot2'))

            # initialize IRAF
            iraf.noao(Stdout=1)
            iraf.digiphot(Stdout=1)
            iraf.daophot(Stdout=1, verify='no')
            iraf.tables(Stdout=1)
            iraf.ttools(Stdout=1)
            iraf.astutil(Stdout=1)

            # set some parameters for iraf photometry
            iraf.fitskypars.setParam('salgori', self.pars.irafphot['salgori'])
            iraf.fitskypars.setParam('annulus', self.pars.irafphot['annulus'])
            iraf.fitskypars.setParam('dannulu', self.pars.irafphot['dannulu'])
            iraf.centerpars.setParam('calgori', self.pars.irafphot['calgori'])
            iraf.centerpars.setParam('cbox',    self.pars.irafphot['cbox'])
            iraf.centerpars.setParam('maxshif', self.pars.irafphot['maxshif'])
            iraf.photpars.setParam('apertur', apsizestr)

            o = iraf.phot(image='@'+os.path.join(self.pars.wdir, self.dataset.framesdir, 'imagelist'), wcsin='world',
                          coords='@'+os.path.join(self.pars.wdir, self.irafphotdir, 'coordlist'),
                          output='@'+os.path.join(self.pars.wdir, self.irafphotdir, 'outputlist'),
                          interactive='no', verify='no', Stdout=1)

            for n in self.dataset.frames:
                # iraf breaks the lines with \ within photometry for a single image (and *\ between different apertures)
                # a single * denotes the end of photometry for a given file.
                # remove these line breaks, so to obtain one line per file with multiple aperture photometry

                basename = os.path.splitext(os.path.basename(self.dataset.frames[n]['path']))[0]
                filename = os.path.join(self.pars.wdir, self.irafphotdir, basename) + '.phot2'
                if os.path.isfile(filename):
                    os.remove(filename)
                o = open(filename, 'a')  # open for append
                for line in open(os.path.join(self.pars.wdir, self.irafphotdir, basename) + '.phot'):
                    if not line[0:1] == '#':
                        if line[-3:-1] == "*\\":
                            writeline = line[:-4].replace('\n', '')
                        elif line[-2:-1] == "\\":
                            writeline = line[:-4].replace('\n', '')
                        else:
                            writeline = line[:-3] + '\n'
                        o.write(writeline)
                o.close()

        # column number in IRAF .phot files for photometry of first aperture
        xc_col = 6
        yc_col = 7
        stdev_col = 15
        nsky_col = 17
        fluxsum_col = 26
        parea_col = 27
        flux_col = 28

        etaslog.info('Save IRAF output to arrays')

        # save IRAF output to arrays.
        for apsize in apsizes:

            xc = []  # x postion of star
            yc = []  # y position of star
            stdev = []  # sky annulus standard deviation
            nsky = []  # number of pixels in annulus
            fluxsum = []  # counts from annulus + aperture
            parea = []  # aperture area
            flux = []  # aperture counts

            self.xc[apsize] = np.empty((self.nobs, self.nstar))
            self.yc[apsize] = np.empty((self.nobs, self.nstar))
            self.skystd[apsize] = np.empty((self.nobs, self.nstar))
            self.psky[apsize] = np.empty((self.nobs, self.nstar))
            self.fluxsum[apsize] = np.empty((self.nobs, self.nstar))
            self.parea[apsize] = np.empty((self.nobs, self.nstar))
            self.flux[apsize] = np.empty((self.nobs, self.nstar))

            for n in self.dataset.frames:

                basename = os.path.splitext(os.path.basename(self.dataset.frames[n]['path']))[0]
                filename = os.path.join(self.pars.wdir, self.irafphotdir, basename) + '.phot2'

                frame = np.loadtxt(filename,
                                   usecols=(xc_col, yc_col, stdev_col, nsky_col, fluxsum_col, parea_col, flux_col),
                                   converters = {xc_col: converter, yc_col: converter, stdev_col: converter,
                                                 nsky_col: converter, fluxsum_col: converter, parea_col: converter,
                                                 flux_col: converter},
                                   unpack=True)

                self.xc[apsize][n] = frame[0]
                self.yc[apsize][n] = frame[1]
                self.skystd[apsize][n] = frame[2]
                self.psky[apsize][n] = frame[3]
                self.fluxsum[apsize][n] = frame[4]
                self.parea[apsize][n] = frame[5]
                self.flux[apsize][n] = frame[6]

            fluxsum_col += 8
            parea_col += 8
            flux_col += 8


class Aperture(object):

    def __init__(self, photometry, apsize):

        self.p = photometry
        self.apsize = apsize
        self.siggcd = []
        self.mag = []
        self.sigmag = []
        self.wt = []
        self.medmag = []
        self.medsigccd = []
        self.medwt = []
        self.med = []
        self.ensmag = []
        self.sigens = []
        self.targetmed = []
        self.emagdiff = []
        self.esigdiff = []
        self.esigflux = []
        self.sigccd = []
        self.signoise = []
        self.eflux = []
        self.star_exclude = []
        self.sigsky = []
        self.sigscalefactor = []
        self.scintnoisecorr = []

        self.rms = []
        self.mags_model = []
        self.sigmastar_model = []
        self.sigmasky_model = []
        self.sigmaread_model = []
        self.sigscint_model = []
        self.sigtot = []

        self.pars = self.p.pars
        self.targetid = self.p.targetid
        self.nobs = self.p.nobs
        self.nstar = self.p.nstar
        self.airmass = self.p.dataset.airmass
        self.aplog = [] # log messages are stored here
        self.sigscint_predicted = self.p.dataset.sigscint
        self.rdnoise = self.p.dataset.params['telescope']['rdnoise']

        # calculate poisson + background noise for this aperture
        nfactor = self.p.parea[apsize]*(1.+(self.p.parea[apsize]/self.p.psky[apsize]))
        starflux = self.p.dataset.params['telescope']['egain']*self.p.flux[apsize]
        self.sigsky = self.p.dataset.params['telescope']['egain']*self.p.skystd[apsize]
        sigstar = np.sqrt(starflux+(nfactor*self.sigsky**2))
        self.sigccd = 1.085736205*sigstar/starflux

        # instrumental magnitudes
        self.mag = 25.0 - 2.5 * np.log10(starflux)

        # perform preliminary outlier rejection (5 sigma from median magnitude)
        #self.outliers = np.empty((self.nobs, self.nstar))

        for i in range(self.nstar):
            ratio = np.abs((self.mag[:,i]-np.nanmedian(self.mag[:,i]))/np.nanstd(self.mag[:,i]))
            self.outliers[:, i] = np.less_equal(ratio, 5)

        # median subtracted magnitudes & errors
        self.medmag = stats.nanmedian(self.mag, 0)
        self.medsigccd = stats.nanmedian(self.sigccd, 0)
        self.med = self.mag - np.array([self.medmag, ]*self.nobs)

    def compute_fluxes(self):


        # default values for error bars scaling factors
        self.sigscalefactor = 1
        self.scintnoisecorr = 1

        log = '[%.2f px] Calculate magnitude ensemble including all stars' % self.apsize
        self.aplog.append(log)
        print log
        self.magnitude_ensemble()

        if self.pars.photometry['scint_noise_corr']:

            log = '[%.2f px] Rescale scintillation noise from RMS distribution' % self.apsize
            print log
            self.aplog.append(log)
            self.scint_noise_corr()

        log = '[%.2f px] Optimize ensemble, after scintillation noise correction' % self.apsize
        print log
        self.aplog.append(log)
        rlm, star_exclude = self.optimize_ensemble(star_exclude=[])

        redchi2 = [rlm[i].polyfit_redchi2 for i in rlm]
        self.sigscalefactor = math.sqrt(np.average(redchi2))
        log = '[%.2f px] Average redchi2 is %.2f. Error bars scaling factor is %.2f' \
              % (self.apsize, np.average(redchi2), self.sigscalefactor)
        print log
        self.aplog.append(log)

        self.magnitude_ensemble(star_exclude=star_exclude)

        log = '[%.2f px] Optimize ensemble, after chi2 adjustment' % self.apsize
        print log
        self.aplog.append(log)
        rlm, self.star_exclude = self.optimize_ensemble(star_exclude=[])

        self.magnitude_ensemble(star_exclude=self.star_exclude)

        redchi2 = [rlm[i].polyfit_redchi2 for i in rlm]
        log = '[%.2f px] Average redchi2 is now %.2f' % (self.apsize, np.average(redchi2))
        print log
        self.aplog.append(log)

        #etaslog.info('[%.2f px] Create outlier rejection map' % apsize
        #self.outlier_rejection(apsize, rlm)

    def magnitude_ensemble(self, refid='planet', star_exclude=[]):

        # refid is the star that will not be included in the ensemble. By default this is the exoplanet
        refid = self.targetid if refid == 'planet' else refid

        # scintillation noise, converted to array (scintillation noise is equal for all star, in each frame)
        self.sigscint_corr = np.array([self.sigscint_predicted]*self.nstar).transpose()*self.scintnoisecorr

        # magnitudes errors, with scaling correction

        self.sigmag = np.sqrt((self.sigccd*self.sigscalefactor)**2+self.sigscint_corr**2)

        # weight
        self.wt = self.outliers * (1/(self.sigmag**2))
        self.medwt = stats.nanmedian(self.wt, 0)

        # create a 2d array with magnitude ensemble for each star in the field
        # it creates nstar ensemble magnitudes,
        # one for each star each ensemble excludes the star being considered and the target (if refid exists)
        # stars can be exluded specifying the input  exclude=[]

        # declare arrays for aperture apsize
        self.ensmag = np.empty((self.nstar, self.nobs))
        self.sigens = np.empty((self.nstar, self.nobs))
        self.targetmed = np.empty((self.nstar, self.nobs))
        self.emagdiff = np.empty((self.nstar, self.nobs))
        self.esigdiff = np.empty((self.nstar, self.nobs))
        self.esigflux = np.empty((self.nstar, self.nobs))
        self.signoise = np.empty((self.nstar, self.nobs))
        self.eflux = np.empty((self.nstar, self.nobs))

        #save a first copy of the weight array
        wt1 = np.copy(self.wt)

        # exclude stars in 'exclude' array (set weight to 0)
        if len(star_exclude) > 0:
            for i in star_exclude:
                wt1[:, i] = 0


        # loop over all stars
        for i in range(self.nstar):


            # save another copy of the weight array and exclude the target star refid (set weight to 0)
            wt2 = np.copy(wt1)
            if isinstance(refid, int):
                wt2[:, refid] = 0

            # exclude the i-th star being considered in the cycle (set weight to 0)
            wt2[:, i] = 0

            # create ensemble magnitudes
            sumens = np.nansum(wt2*self.med, 1)
            sumwt = np.nansum(wt2, 1)
            self.ensmag[i] = sumens / sumwt
            self.sigens[i] = np.sqrt(1/sumwt)
            self.targetmed[i] = self.med[:, i]
            self.emagdiff[i] = self.ensmag[i] - self.targetmed[i]
            sig1 = self.sigmag[:, i]
            sig2 = self.sigens[i]
            self.esigdiff[i] = np.sqrt((sig1*sig1)+(sig2*sig2))
            self.eflux[i] = 10.0**(0.4*self.emagdiff[i])
            self.esigflux[i] = self.esigdiff[i]*self.eflux[i]/1.085736205
            self.signoise[i] = self.eflux[i]/self.esigflux[i]


    def scint_noise_corr(self, refid='planet'):

        # this function computes a correction to the theoretical signal to noise

        refid = self.targetid if refid == 'planet' else refid

        # optimize ensemble, excluding very bad stars
        # rlm[i] is the robust linear model for each star i; star_exclude is a list with the excluded stars
        rlm, star_exclude = self.optimize_ensemble(refid, customlimits=[0.3, 4])

        # rlm[i].polyfit_rms returns the rms of the robust linear fit to each star i
        rms = []
        sigmaccd = []
        medmag = []
        sigsky = []
        sigskymed = stats.nanmedian(self.sigsky, 0)

        for n in rlm:
            rms.append(rlm[n].polyfit_rms)
            sigmaccd.append(self.medsigccd[n])
            medmag.append(self.medmag[n])
            sigsky.append(sigskymed[n])
        medmag = np.asarray(medmag)
        sigsky = np.asarray(sigsky)
        rms = np.asarray(rms)

        # take the median of the squared difference between the RMS of the linear fit for each star
        # and the calculated ccd noise. Divide it by the median of the scintillation noise amongst frames
        # to get the scaling ratio.

        parea = math.pi*self.apsize**2
        psky = math.pi*(float(self.pars.irafphot['annulus'])+float(self.pars.irafphot['dannulu']))**2-\
               math.pi*float(self.pars.irafphot['annulus'])**2
        sigmastar = np.sqrt(10**(-0.4*(medmag-25)))/10**(-0.4*(medmag-25))
        sigmasky = np.sqrt(parea*(1+parea/psky))*np.median(sigsky)/10**-(0.4*(medmag-25))
        sigmaread = math.sqrt(parea)*self.rdnoise/10**-(0.4*(medmag-25))

        partialnoise = sigmastar**2+sigmasky**2+sigmaread**2

        from scipy.optimize import leastsq
        correction = 1
        out = leastsq(self.fit_scint_noise, correction, (partialnoise, rms, medmag))
        correction = out[0][0]

        self.scintnoisecorr = correction
        log = '[%.2f px] Scintillation noise correction factor is %.4f' % (self.apsize, self.scintnoisecorr)
        print log
        self.aplog.append(log)


        # recalculate magnitude ensemble
        self.magnitude_ensemble()

        mags = np.arange(np.min(medmag-1), np.max(medmag+1), 0.01)

        sigmastar = np.sqrt(10**(-0.4*(mags-25)))/10**(-0.4*(mags-25))
        sigmasky = (np.sqrt(parea*(1+parea/psky))*np.median(sigsky))/10**-(0.4*(mags-25))
        sigmaread = math.sqrt(parea)*self.rdnoise/10**-(0.4*(mags-25))
        sigscint = stats.nanmedian(stats.nanmedian(self.sigscint_corr))
        totnoise = np.sqrt(sigmastar**2+sigmasky**2+sigmaread**2+sigscint**2)

        # these can be used to plot noise model versus observed noise.
        self.medmag = medmag
        self.rms = rms
        self.mags_model = mags
        self.sigmastar_model = sigmastar
        self.sigmasky_model = sigmasky
        self.sigmaread_model = sigmaread
        self.sigscint_model = [sigscint]*len(mags)
        self.sigtot_model = np.sqrt(sigmastar**2+sigmasky**2+sigmaread**2+sigscint**2)

        #plt.scatter(medmag, rms)
        #plt.plot(mags, sigmastar, label='Star noise')
        ##plt.plot(mags, sigmasky, label='Sky noise')
        #plt.plot(mags, sigmaread, label='Readout noise')
        #plt.plot([9, 16], [sigscint, sigscint], label='Scint. noise')
        #plt.yscale('log')
        #plt.plot(mags, np.sqrt(sigmastar**2+sigmasky**2+sigmaread**2+sigscint**2), label='Tot. noise')
        #plt.legend(loc=4)
        #plt.show()

    def fit_scint_noise(self, correction, partialnoise, rms, mag):
        sigscint = correction*self.sigscint_predicted
        medsigscint = stats.nanmedian(sigscint, 0)
        totnoise = np.sqrt(partialnoise+medsigscint**2)
        return (rms - totnoise)/rms**2

    def optimize_ensemble(self, refid='planet', chi2limit=0.995, customlimits=[], star_exclude=[]):

        # performs a chi-squared test on each light curve o check if the distributions are consistent with
        # straight lines. T-rex calculates the normalized chi squared values for each frame and performs an
        # optimisation routine excluding all the bad comparison stars.

        refid = self.targetid if refid == 'planet' else refid

        rlm = {}
        maxredchi2 = 0
        maxredchi2id = 0
        minredchi2 = 1000
        minredchi2id = 0

        # loop all stars exept those in star_exclude and the current star and run instance of RobustStatistics
        # for each one
        for i in range(self.nstar):

            if not i in star_exclude and not i == refid:

                x = np.copy(self.p.dataset.airmass)
                y = np.copy(self.eflux[i])
                sig = np.copy(self.esigflux[i])
                rlm[i] = RobustStatistic(x, y, sig, chi2limit=chi2limit, customlimits=customlimits)

                # get maximum and minimum redchi2
                redchi2 = rlm[i].polyfit_redchi2
                if redchi2 > maxredchi2:
                    maxredchi2 = rlm[i].polyfit_redchi2
                    maxredchi2id = i
                if redchi2 < minredchi2:
                    minredchi2 = rlm[i].polyfit_redchi2
                    minredchi2id = i

        if not rlm[maxredchi2id].chi2_within_limits():
            star_exclude.append(maxredchi2id)
            self.magnitude_ensemble( star_exclude=star_exclude)
            log = '[%.2f px] Rejecting star n %i (redchi2 = %.3f, bottomlim = %.2f, uplim = %.2f)' % \
                  (self.apsize, maxredchi2id, maxredchi2,
                  rlm[maxredchi2id].polyfit_chi2limits[0],
                  rlm[maxredchi2id].polyfit_chi2limits[1])
            print log
            self.aplog.append(log)


            return self.optimize_ensemble(refid, chi2limit, customlimits, star_exclude)

        elif not rlm[minredchi2id].chi2_within_limits():
            star_exclude.append(minredchi2id)
            self.magnitude_ensemble(star_exclude=star_exclude)
            log = '[%.2f px] Rejecting star n %i (redchi2 = %.3f, bottomlim = %.2f, uplim = %.2f)' %\
                  (self.apsize, minredchi2id, minredchi2,
                   rlm[minredchi2id].polyfit_chi2limits[0],
                   rlm[minredchi2id].polyfit_chi2limits[1])
            print log
            self.aplog.append(log)
            return self.optimize_ensemble(refid, chi2limit, customlimits, star_exclude)
        else:
            log = '[%.2f px] Optimize ensemble output - Excluded stars: %s' % (self.apsize, star_exclude)
            self.aplog.append(log)
            print log
            return rlm, star_exclude

    def outlier_rejection(self, rlm=False, refid='planet'):

        # work in progress. Experiment?

        if not rlm:
            refid = self.targetid if refid == 'planet' else refid
            rlm = {}
            # loop all stars exept those in star_exclude and the current star and run instance of RobustStatistics
            # for each one
            for i in range(self.nstar):
                if not i in self.star_exclude and not i == refid:
                    x = np.copy(self.airmass)
                    y = np.copy(self.eflux[i])
                    sig = np.copy(self.esigflux[i])
                    rlm[i] = RobustStatistic(x, y, sig)

        residuals = np.empty((self.nstar, self.nobs))
        errors = np.empty((self.nstar, self.nobs))
        weights = np.empty((self.nstar, self.nobs))

        for i in range(self.nstar-1):
            weight = self.medwt[i]
            for j in range(self.nobs-1):
                if not i in self.star_exclude and not i == refid:
                    calculated = rlm[i].afit*self.airmass[j]+rlm[i].bfit
                    observed = self.eflux[i][j]
                    residuals[i][j] = np.abs(observed-calculated)
                    errors[i][j] = self.esigflux[i][j]
                    weights[i][j] = residuals[i][j]*weight
        ratios_frame = np.nanmean(np.abs(weights), 0)
        self.rlm = rlm
        self.rejectionweights = ratios_frame
        for i in range(self.nobs-1):
            print ratios_frame[i]


class MultiThread(Process):

    def __init__(self, apsize, photometry, queue):
        etaslog.info('[%.2f px] Start thread for %.2f px aperture' % (apsize, apsize))
        Process.__init__(self)
        self.pars = photometry.pars
        self.photometry = photometry
        self.apsize = apsize
        self.aperture = Aperture(photometry, apsize)
        self.queue = queue
        return

    def run(self):
        if self.pars.photometry['compute_fluxes']:
            self.aperture.compute_fluxes()
        del self.aperture.p
        self.queue.put(self.aperture)
        return