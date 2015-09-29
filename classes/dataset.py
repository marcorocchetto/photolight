# import modules
import os
import shutil
import math
import alipy
import subprocess
import sys
import numpy as np
import logging
from scipy import stats
from pyraf import iraf
from functions import *
from list_frames import *

class Dataset(object):

    def __init__(self, target, odir=None, params=None):

        import datetime

        logging.info('Initialize class Lightcurve')

        # inherit parameters class from target
        if not params:
            self.params = target.params
        else:
            self.params = params

        # store the target instance locally
        self.target = target

        # source files directory
        if odir:
            self.odir = odir
        else:
            if self.params.odir:
                self.odir = self.params.odir
            else:
                logging.error('Source directory not specified.')
                sys.exit()
        logging.info('Source directory of fits: %s' % self.odir)

        # initialise some instance variables
        self.targetid = None
        self.fwhm = None
        self.avgfwhm = None
        self.stdfwhm = None
        self.ra = None
        self.dec = None
        self.nstar = None
        self.bjd = None
        self.hjd = None
        self.sigscint = None
        self.medsigscint = None



        # get list of frames from source directory
        self.frames = list_frames(self.odir)
        logging.info('There are %i valid FITS frames in the source directory' % len(self.frames))

        # Check number of plate solved frames. Exclude non plate solved frames if ignore_unsolved_frames = True
        npltsolved = 0
        for n in self.frames:
            if self.frames[n]['plate_solved']:
                npltsolved += 1
            else:
                logging.info('Frame %s is not plate solved.' % self.frames[n]['path'])
        logging.info('There are %i plate solved frames' % npltsolved)

        if npltsolved < len(self.frames) and not self.params.dataset['ignore_unsolved_frames']:
            logging.error('Some frames do not have a plate solution. Fix these '
                          'frames or set ignore_unsolved_frames = Off')
            sys.exit()

        # elif npltsolved < len(self.frames):
        #     logging.error('Excluding non plate solved frames')
        #     for n in self.frames:
        #         if not self.frames[n]['plate_solved']:
        #             del self.frames[n] # @todo are we sure we're doing it correctly????
        #

        # # Proceed only if there are sufficient frames MOVE LATER TO CREATE MASTER
        # if npltsolved < self.params.dataset['min_frames']:
        #     logging.error('There are not enough plate solved frames to continue (min_frames = %i)'
        #                   % self.params.dataset['min_frames'])
        #     sys.exit()

        # Check if all frames have same size (in pixel), then save image size
        img_x = -1; img_y = -1
        for n in self.frames:
            outputshape = alipy.align.shape(self.frames[n]['path'], verbose=False)
            if img_x < 0:
                img_x = outputshape[0]
            if img_y < 0:
                img_y = outputshape[1]
            if img_x != outputshape[0] or img_y != outputshape[1]:
                    logging.error('Frames do not have the same size.'
                                  'Frame %s is %.1f x %.1f,'
                                  'Frame %s is %.1f x %.1f,'
                                  % (self.frames[n-1]['path'], img_x, img_y, self.frames[n]['path'], img_x, img_y))
                    sys.exit()
            img_x = outputshape[0]
            img_y = outputshape[1]
        self.img_x = img_x
        self.img_y = img_y
        logging.info('Frames pixel size: %.1f x %.1f' % (self.img_x, self.img_y))

        # set: nobs; airmass and exptime arrays
        self.nobs = len(self.frames)
        self.airmass = np.empty(self.nobs)
        self.exptime = np.empty(self.nobs)
        self.path = {}
        for n in self.frames:
            self.airmass[n] = self.frames[n]['AIRMASS']
            self.exptime[n] = self.frames[n]['EXPTIME']
            self.path[n] = self.frames[n]['path']

        # get mean calendar and julian date of dataset
        ts = 0
        for n in self.frames:
            ts += self.frames[n]['timestamp']
        tsmid = ts/len(self.frames)
        self.cdate = datetime.datetime.fromtimestamp(int(tsmid)).strftime('%Y-%m-%d')
        self.jdate = int(float((tsmid/86400.0)+2440587.5)) # floor approx
        logging.info('Calendar date is %s, Julian date is %s' % (self.cdate, self.jdate))

        # set directory names
        self.datedir = os.path.join(self.target.mname, str(self.cdate))
        self.telescopedir = os.path.join(self.datedir, self.params.telescope['name'])
        self.framesdir = os.path.join(self.telescopedir, 'frames')
        self.masterdir = os.path.join(self.telescopedir, 'master')

        # master frame path
        self.masterfilepath = os.path.join(self.masterdir, 'master.fits')
        self.masterpngpath = os.path.join(self.params.wdir, self.masterdir, 'master.png')

        # create directories for the julian date and telescope
        create_dirs = [os.path.join(self.params.wdir, self.target.mname),
                       os.path.join(self.params.wdir, self.datedir),
                       os.path.join(self.params.wdir, self.telescopedir),
                       os.path.join(self.params.wdir, self.framesdir),
                       os.path.join(self.params.wdir, self.masterdir)]
        for directory in create_dirs:
            if os.path.isdir(directory):
                logging.info('The directory %s already exists. Remove it' % directory)
                shutil.rmtree(directory)
            try:
                logging.info('Creating directory %s' % directory)
                os.mkdir(directory)
            except Exception, e:
                logging.error('Failed to create the directory %s. Cannot continue' % directory, exc_info=True)
                sys.exit()
        logging.info('Dataset instance correctly initialised')

    def create_master_frame(self):

        # Create Master Frame
        logging.info('Creating master frame: %s ' % os.path.join(self.params.wdir, self.masterfilepath))
        if os.path.isfile(os.path.join(self.params.wdir, self.masterfilepath)):
            logging.info('Removing previous version of master frame')
            os.remove(os.path.join(self.params.wdir, self.masterfilepath))

        logging.info('Stack frames using plate solution of first %i frames' % self.params.dataset['master_nstack'])

        # get average RA and DEC of frames
        crval1 = []
        crval2 = []
        for n in self.frames:
            crval1.append(self.frames[n]['CRVAL1'])
            crval2.append(self.frames[n]['CRVAL2'])
        crval1avg = np.average(crval1)
        crval2avg = np.average(crval2)
        logging.info('Average RA is %s' % deg_to_sex(ra=crval1avg))
        logging.info('Average Dec is %s' % deg_to_sex(dec=crval2avg))

        logging.info('Selecting well aligned plate solved frames. Max offset is %.2f arcmin' %
                     self.params.dataset['frame_max_offset'])

        # prepare a file with a list of well aligned plate solved frames. Used by iraf imcombine
        stacklist = os.path.join(self.params.wdir, self.masterdir, 'stacklist')
        tmpf = open(stacklist, 'w')

        # loop the  frames and select only those that are not offset by more than self.params.framemaxoffset
        ra_shift = crval1avg  # average ra
        dec_shift = crval2avg  # average dec
        nstack = 0
        for n in self.frames:
            ra = float(self.frames[n]['CRVAL1'])  # frame ra
            dec = float(self.frames[n]['CRVAL2'])  # frame dec
            offset = math.sqrt(math.fabs(ra-ra_shift)**2+math.fabs(dec-dec_shift)**2)*60
            if offset < self.params.dataset['frame_max_offset']:
                    if nstack < self.params.dataset['master_nstack']:
                        tmpf.write("%s\n" % self.frames[n]['path'])  # write path to file
                        nstack += 1
                    else:
                        break  # close loop when we reach the minimum number of frames
            else:
                logging.warning('Frame %s is offset is %.1f, more than limit (%.1f)' %
                                (self.frames[n]['path'], offset, self.params.dataset['frame_max_offset']))
        tmpf.close()

        # nstack is the number of well aligned plate solved frames that can be stacked
        if nstack == self.params.dataset['master_nstack']:

            logging.info('Crate master frame using iraf.imcombine')
            if os.path.isfile(os.path.join(self.params.wdir, self.masterfilepath)):
                # remove the master file in case it exists (iraf cannot overwrite!)
                logging.warning('The master frame already exists (%s). Remove it.'
                                % os.path.join(self.params.wdir, self.masterfilepath))
                os.remove(os.path.join(self.params.wdir, self.masterfilepath))

            iraf.images()
            iraf.imcombine(input='@%s' % stacklist,
                           output=os.path.join(self.params.wdir, self.masterfilepath),
                           scale=self.params.iraf_imcombine['scale'],
                           weight=self.params.iraf_imcombine['weight'],
                           combine=self.params.iraf_imcombine['combine'],
                           offsets='world',
                           Stdout=0)
        if os.path.isfile(os.path.join(self.params.wdir, self.masterfilepath)):
            logging.info('Master frame correctly created: %s' % os.path.join(self.params.wdir, self.masterfilepath))
        else:
            logging.error('Master frame could not be created')
            sys.exit()


        # Generate master catalogue of stars

        # extract star catalogue from sextractor output
        logging.info('Run sextractor on master frame to obtain coordinates of stars in WCS.')
        logging.info('Detection treshold is %i' % self.params.dataset['detect_tresh'])

        # master.cat contains the WCS coordinates of the stars to perform photometry on
        sexrun = '%s ' % self.params.sex
        sexrun += '%s ' % os.path.join(self.params.wdir, self.masterfilepath)
        sexrun += '-c %s ' % os.path.join(self.params.aconfdir, 'master.sex')
        sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.params.aconfdir, 'master.param')
        sexrun += '-FILTER_NAME %s ' % os.path.join(self.params.aconfdir, 'default.conv')
        sexrun += '-CATALOG_TYPE ASCII '
        sexrun += '-DETECT_THRESH %i ' % self.params.dataset['detect_tresh']
        sexrun += '-CATALOG_NAME %s ' % os.path.join(self.params.wdir, self.masterdir, 'master.cat')


        p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]

        # load catalogue of stars in array using numpy
        mastercat_orig = np.loadtxt(os.path.join(self.params.wdir, self.masterdir,  'master.cat'))
        mastercat_sorted = mastercat_orig[mastercat_orig[:,4].argsort()]  # sort by ISO FLUX (4th column)

        # remove edge star
        logging.info('Remove edge stars. Edge is %i pixels, frame is %ix%i'
                     % (self.params.dataset['edge_limit'], self.img_x, self.img_y))

        mastercat_edge = []
        j = 0
        for line in mastercat_sorted:
            if (line[2] < (self.img_x - self.params.dataset['edge_limit'])) and (line[2] > self.params.dataset['edge_limit']) and \
               (line[3] < (self.img_y - self.params.dataset['edge_limit'])) and (line[3] > self.params.dataset['edge_limit']):  # y axis
                mastercat_edge.append([line[0], line[1], line[2], line[3], line[4], line[5]])
                j += 1
            else:
                logging.info('Exclude edge star at x: %.1f y: %.1f ' % (line[2], line[3]))

        # if the stars in the catalogue are more than the maximum allowed (self.params.dataset['master_max_stars']), take only the last
        # self.params.dataset['master_max_stars'] stars in the list (the list is already sorted by flux in ascending order)
        if j > self.params.dataset['master_max_stars']:
            logging.info("%i stars found (max %i). Exclude faintest stars. " % (j, self.params.dataset['master_max_stars']))
            mastercat_final = []
            for n in range(self.params.dataset['master_max_stars']):
                mastercat_final.append(mastercat_edge[j-n-1])
        else:
            mastercat_final = mastercat_edge

        # check for the presence of the target star in the final master catalogue
        logging.info('Looking for target star...')

        k = 0
        for line in mastercat_final:
            xshift = float(self.target.dec_d - line[1])
            yshift = float(self.target.ra_d - line[0])
            shift = math.sqrt(xshift**2 + yshift**2)
            if shift < 0.005:
                self.targetid = k
                break
            k += 1

        if not self.targetid:

            # the planet was not found in the final list, check if it's been excluded
            logging.info('The target was not found in the final catalogue. Check for fainter stars...')
            for line in mastercat_edge:
                if math.fabs(float(self.target.dec_d - line[1])) < 0.005 and \
                   math.fabs(float(self.target.ra_d - line[0])) < 0.005:
                    # the target is in the original master catalogue of stars. It was too faint to be included
                    # in the catalogue. Add it now and exclude faintest star.
                    del mastercat_final[len(mastercat_final)-1]
                    mastercat_final.append(line)
                    self.targetid = len(mastercat_final)-1
            try:
                self.targetid
            except:
                logging.info('The target was not found! It is not possible to continue with the photometry')
                sys.exit()

        if not self.targetid:
            logging.error('The target was not found in the final catalogue. ')
            sys.exit()

        logging.info('Target detected! %s is star number %i' % (self.target.name, self.targetid))

        # rewrite an iraf.tvmark version of the catalogue
        mastercat_tvmark = []
        for line in mastercat_final:
            mastercat_tvmark.append([line[2], line[3]])

        # convert catalogues to numpy arrays
        mastercat_final = np.array(mastercat_final)
        mastercat_tvmark = np.array(mastercat_tvmark)

        # write output to files
        filename = os.path.join(self.params.wdir, self.masterdir, 'master.final')
        logging.info('Output file with catalogue stars: %s' % filename)
        np.savetxt(filename, mastercat_final)
        filename = os.path.join(self.params.wdir, self.masterdir, 'master.tvmark')
        logging.info('Output file with catalogue stars: %s' % filename)
        np.savetxt(filename, mastercat_tvmark)

        # generate png of master file with catalogue stars
        logging.info('Generate PNG of master file with catalogue stars:  %s' % (self.masterpngpath))
        coords = []
        n = 0
        for line in mastercat_final:
            coord = {}
            coord['x'] = line[2]
            coord['y'] = line[3]
            coord['id'] = str(n)
            coords.append(coord)
            n += 1
        fits_to_png(os.path.join(self.params.wdir, self.masterdir, 'master.fits'),
                    os.path.join(self.params.wdir, self.masterdir, 'master.png'), coords)

        # get coordinates of stars
        logging.info('Get coordinate of stars from the master catalogue of stars')
        coord = np.loadtxt(os.path.join(self.params.wdir, self.masterdir, 'master.final'), unpack=True)
        self.ra = coord[0]
        self.dec = coord[1]
        self.nstar = len(coord[0])


        # Create SExtractor catalogue for all frames and
        if self.params.dataset['run_all_sex']:

            logging.info('Create final SExtractor catalogue for all frames')

            for frame in self.frames:
                basename = os.path.splitext(os.path.basename(self.frames[frame]['path']))[0]
                cataloguename = os.path.join(self.params.wdir, self.framesdir, basename) + '.cat'

                sexrun = '%s ' % self.params.sex
                sexrun += '%s ' % self.frames[frame]['path']
                sexrun += '-c %s ' % os.path.join(self.params.aconfdir, 'master.sex')
                sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.params.aconfdir, 'default.param')
                sexrun += '-FILTER_NAME %s ' % os.path.join(self.params.aconfdir, 'default.conv')
                sexrun += '-CATALOG_TYPE ASCII '
                sexrun += '-DETECT_THRESH 60 '
                sexrun += '-CATALOG_NAME %s ' % cataloguename
                p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]

            #  get FWHM of all stars in all frames
            mastercat = np.loadtxt(os.path.join(self.params.wdir, self.masterdir, 'master.final'))
            self.fwhm = np.empty((self.nobs, self.nstar))
            i = 0
            for frame in self.frames:
                basename = os.path.splitext(os.path.basename(self.frames[frame]['path']))[0]
                cataloguename = os.path.join(self.params.wdir, self.framesdir, basename) + '.cat'
                catalogue = np.loadtxt(cataloguename)
                k = 0
                for line in mastercat:
                    ra = line[0]
                    dec = line[1]
                    ck = False
                    for line1 in catalogue:
                        ra1 = line1[0]
                        dec1 = line1[1]
                        shift = math.sqrt(math.fabs(ra-ra1)**2 + math.fabs(dec-dec1)**2)
                        if shift < 0.005:
                            self.fwhm[i][k] = line1[5]
                            ck = True
                            k += 1
                            break
                    if not ck:
                        self.fwhm[i][k] = np.NaN
                        k += 1
                i += 1

            avg1 = stats.stats.nanmean(self.fwhm, 0)
            avg2 = stats.stats.nanmean(avg1)
            std = stats.stats.nanstd(avg1)
            self.avgfwhm = avg2
            self.stdfwhm = std
            logging.info('Average FWHM is %.1f +/- %.1f' % (self.avgfwhm, self.stdfwhm))

        # Calculate LMST and AIRMASS values for all frames using IRAF
        '''
        logging.info('Updating airmass values in fits headers')

        iraf.noao()
        iraf.astutil()
        iraf.observatory.setParam('obsid', 'obspars')
        iraf.observatory.setParam('command', 'set')
        iraf.observatory.setParam('observatory', 'obspars')
        iraf.observatory.setParam('name', self.params.dataset.observatory['name'])
        iraf.observatory.setParam('longitu', str(self.params.dataset.observatory['longitude']))
        iraf.observatory.setParam('altitud', str(self.params.dataset.observatory['elevation']))
        iraf.observatory.setParam('latitud', str(self.params.dataset.observatory['latitude']))
        iraf.observatory.setParam('override', 'obspars')

        f1 = open(os.path.join(self.params.wdir, self.framesdir, 'imagelist'), 'w')
        for n in self.frames:
            f1.write(self.frames[n]['path'] + '\n')
        f1.close()

        logging.info('Calculating LMST')
        f = open(os.path.join(self.params.wdir, self.framesdir, 'lmst-calc'), 'w')
        f.write("observat = 'obspars'\n")
        f.write("lmst = mst (@'DATE-OBS', @'TIME-OBS', obsdb (observat, 'longitude'))\n")
        f.write("quit")
        f.close()
        iraf.asthedit(images='@'+os.path.join(self.params.wdir, self.framesdir, 'imagelist'),
                      commands=os.path.join(self.params.wdir, self.framesdir, 'lmst-calc'),
                      Stdout=1, verbose='no')

        logging.info('Calculating airmass (an extra header EAIRMASS is added to the fits files)')
        iraf.setairmass(images='@'+os.path.join(self.params.wdir, self.framesdir, 'imagelist'), observa='obspars', \
                        intype='beginning', outtype='effective', st='lmst', \
                        ra='OBJCTRA', dec='OBJCTDEC', equinox='EPOCH', ut='TIME-OBS', \
                        date='DATE-OBS', exposur='EXPTIME', airmass='EAIRMASS', utmiddl='utmiddle', \
                        scale='750.', update='yes', Stdout=1) '''

        # Calculate BJD and HJD for each frame using vartools
        logging.info('Calculate BJD and HJD for each frame')

        bjd_filename = os.path.join(self.params.wdir, self.framesdir, 'bjdtbd')
        hjd_filename = os.path.join(self.params.wdir, self.framesdir, 'hjd')
        utc_filename = os.path.join(self.params.wdir, self.framesdir, 'utc')

        # load UTC for each frame from header
        import datetime
        f = open(utc_filename, 'w')
        for n in self.frames:
            ts = float(self.frames[n]['timestamp']) + float(self.frames[n]['EXPTIME'])/2
            utc = datetime.datetime.fromtimestamp(int(ts)).strftime('%Y-%m-%dT%H:%M:%S')
            f.write('%s %i %i \n' % (utc, 0, 0))
        f.close()

        # ra and dec of target, used to calculate BJD
        ra = self.target.ra_d
        dec = self.target.dec_d

        # run vartools -converttime, write the BJD values to bjd_filename
        runvartools = self.params.vartools + " -i " + utc_filename
        runvartools += " -quiet -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 "
        runvartools += "-converttime input jd inputsys-utc output bjd outputsys-tdb radec fix %f %f " % (ra, dec)
        runvartools += "ephemfile %s/de421.bsp " % self.params.vartoolsdir
        runvartools += "leapsecfile %s/naif0010.tls " % self.params.vartoolsdir
        runvartools += "planetdatafile %s/pck00010.tpc " % self.params.vartoolsdir
        runvartools += "coords fix %f %f %f "  % (self.params.observatory['latitude'],
                                                  self.params.observatory['longitude'],
                                                  self.params.observatory['elevation'])
        runvartools += "-o %s " % bjd_filename

        p = subprocess.Popen(runvartools, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]
        self.bjd = np.loadtxt(bjd_filename, unpack=True)[0]

        # run vartools -converttime, write the HJD values to hjd_filename
        runvartools = self.params.vartools + " -i " + utc_filename
        runvartools += " -quiet -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 "
        runvartools += "-converttime input jd inputsys-utc output hjd outputsys-utc radec fix %f %f " % (ra, dec)
        runvartools += "ephemfile %s/de421.bsp " % self.params.vartoolsdir
        runvartools += "leapsecfile %s/naif0010.tls " % self.params.vartoolsdir
        runvartools += "planetdatafile %s/pck00010.tpc " % self.params.vartoolsdir
        runvartools += "coords fix %f %f %f " % (self.params.observatory['latitude'],
                                                 self.params.observatory['longitude'],
                                                 self.params.observatory['elevation'])
        runvartools += "-o %s " % hjd_filename
        p = subprocess.Popen(runvartools, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]
        self.hjd = np.loadtxt(hjd_filename, unpack=True)[0]

        # Calculate theoretical scintillation noise theoretical (will be added in quadrature to the photometry.ccdnoise)
        logging.info('Calculate theoretical scintillation noise')
        scintwav = (self.params.filter['midwav']/550.0)**(-7.0/12.0)
        scintalti = math.exp(-self.params.observatory['elevation']/8000.0)
        scintap = self.params.telescope['aperture']**(2.0/3.0)
        scintfact = 0.09*scintalti*scintwav/scintap
        self.sigscint = scintfact*(self.airmass**1.75)/np.sqrt(2.*self.exptime)
        self.medsigscint = stats.nanmedian(self.sigscint, 0)

        logging.info('Master frame and master catalogue of stars correctly created.')
