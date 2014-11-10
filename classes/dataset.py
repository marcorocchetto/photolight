
import os
import shutil
import math
import glob
import alipy
import subprocess
import sys
import numpy as np
import logging
import imp

from scipy import stats
from pyraf import iraf

from ..library.functions import *
from ..library.list_frames import *

from target import Target

etaslog = logging.getLogger('etaslog')

class Dataset(object):

    def __init__(self, object, **kwargs):

        import datetime

        etaslog.info('Initialize class Lightcurve')

        # inherit parameters class from object
        self.pars = object.pars

        # store the object instance (planet, star etc) locally
        self.object = object

        # source files directory
        self.odir = kwargs['odir']
        etaslog.info('Source directory of fits: %s' % self.odir)

        # get list of frames from source directory
        oframes = list_frames(self.odir, out='keys', keys=self.pars.platekeys+self.pars.stkeys+['PLTSOLVED'])[0]

        etaslog.info('There are %i valid FITS frames in the source directory' % len(oframes))

        # get mean calendar and julian date of dataset
        ts = 0
        for n in oframes:
            ts += oframes[n]['timestamp']
        tsmid = ts/len(oframes)

        self.cdate = datetime.datetime.fromtimestamp(int(tsmid)).strftime('%Y-%m-%d')
        self.jdate = int(float((tsmid/86400.0)+2440587.5)) # floor approx
        etaslog.info('Calendar date is %s, Julian date is %s' % (self.cdate, self.jdate))

        # set various parameters:
        self.params = {}

        # set observatory, telescope and filter parameters
        reqobjects = ['observatory', 'telescope', 'filter']
        reqparams = {'observatory': ['elevation', 'longitude', 'latitude'],
                     'telescope': ['name', 'egain', 'rdnoise', 'aperture'],
                     'filter': ['midwav']}
        for obj in reqobjects:
            if obj in kwargs:
                for param in reqparams[obj]:
                    if not param in kwargs[obj]:
                        etaslog.error('Parameter "%s" for "%s" missing. Check input parameters of Dataset' %
                                      (param, obj))
                        sys.exit()
                self.params[obj] = kwargs[obj]

        # set other parameters (self.pars.dataset[], look at parameters class)
        for param in self.pars.dataset:
            if param in kwargs:
                if isinstance(kwargs[param], type(self.pars.dataset[param])):
                    self.params[param] = kwargs[param]
                    etaslog.info('Parameter "%s" set to %i' % (param, self.params[param]))
                else:
                    etaslog.warning('Parameter "%s" is not valid. Use default: %s' % (param, str(self.pars.dataset[param])))
            else:
                etaslog.info('Use default parameter for "%s": %s' % (param, str(self.pars.dataset[param])))
                self.params[param] = self.pars.dataset[param]

        # create directories fo the julian date and telescope
        self.datedir = os.path.join(self.object.mname, str(self.cdate))
        self.telescopedir = os.path.join(self.datedir, self.params['telescope']['name'])
        self.framesdir = os.path.join(self.telescopedir, 'frames')
        self.masterdir = os.path.join(self.telescopedir, 'master')
        create_dirs = [os.path.join(self.pars.wdir, self.object.mname),
                       os.path.join(self.pars.wdir, self.datedir),
                       os.path.join(self.pars.wdir, self.telescopedir),
                       os.path.join(self.pars.wdir, self.framesdir),
                       os.path.join(self.pars.wdir, self.masterdir)]
        for directory in create_dirs:
            if os.path.isdir(directory):
                etaslog.info('The directory %s already exists. Remove it' % directory)
                shutil.rmtree(directory)
            try:
                etaslog.info('Creating directory %s' % directory)
                os.mkdir(directory)
            except Exception, e:
                etaslog.error('Failed to create the directory %s. Cannot continue' % directory, exc_info=True)
                sys.exit()

        # Check if all frames have same size (in pixel), then save image size
        img_x = -1; img_y = -1
        for n in oframes:
            outputshape = alipy.align.shape(oframes[n]['path'], verbose=False)
            if img_x < 0:
                img_x = outputshape[0]
            if img_y < 0:
                img_y = outputshape[1]

            if img_x != outputshape[0] or img_y != outputshape[1]:
                    etaslog.error('Frame do not have the same size.'
                                  'Frame %s is %.1f x %.1f,'
                                  'Frame %s is %.1f x %.1f,'
                                  % (oframes[n-1]['path'], img_x, img_y, oframes[n]['path'], img_x, img_y))
                    sys.exit()
            img_x = outputshape[0]
            img_y = outputshape[1]
        self.img_x = img_x
        self.img_y = img_y

        etaslog.info('Frames pixel size: %.1f x %.1f' % (self.img_x, self.img_y))

        ########################################################
        # Copy or link source files to local working directory #
        ########################################################

        etaslog.info('Copy or link source files to local working directory: %s'
                     % os.path.join(self.pars.wdir, self.framesdir))

        # default is to link source files to local working directory. Otherwise copy.
        if not 'copy' in kwargs:
            copy = False
        else:
            copy = kwargs['copy']

        if copy:
            etaslog.info('Copying %i files to working directory ' % len(oframes))
        else:

            import stat
            st = os.stat(os.path.join(self.pars.wdir, self.framesdir))
            if bool(st.st_mode & stat.S_IRGRP):
                etaslog.info('Linking %i files to working directory' % len(oframes))
                etaslog.warning('Be carful linking files... Source files will be modified (plate solution appended to '
                                'unsolved frames, fits headers checked and sanitised, EAIRMASS added to fits '
                                'headers)')
            else:
                etaslog.error('You asked to link the source files to the working directory but you do not have'
                              'write permission on the source files. Change the permissions or copy the files')
                sys.exit()

        for n in oframes:
            filen = os.path.join(self.pars.wdir, self.framesdir, 'd%s.fits' % str(n).zfill(5))
            if not os.path.exists(filen):
                etaslog.debug('From %s to %s' % (oframes[n]['path'], filen))
                if copy:
                    shutil.copy(oframes[n]['path'], filen)
                else:
                    os.symlink(oframes[n]['path'], filen)

        ###################################
        # Check and sanitise fits headers #
        ###################################

        import pyfits
        failed_frames = []
        for n in oframes:
            keys = oframes[n].keys()
            for header in self.pars.stkeys: # check that standard header keys are present
                if not header in keys:
                    # adjust some headers
                    if header in self.pars.stkeys:
                        hdulist = pyfits.open(oframes[n]['path'], mode='update')
                        prihdr = hdulist[0].header

                        if header == 'EPOCH':
                            # assume J2000
                            etaslog.info('EPOCH key missing for frame %s, assume J2000'
                                         % os.path.basename(oframes[n]['path']))
                            prihdr.update('EPOCH', '2000.0')
                        elif header == 'OBJCTRA' or header == 'OBJCTDEC':
                            etaslog.info('OBJCTRA and/or OBJCTDEC keys missing for frame %s, assume coordinates given'
                                         'by the object'
                                         % os.path.basename(oframes[n]['path']))
                            prihdr.update('OBJCTRA', str(self.object.ra))
                            prihdr.update('OBJCTDEC', str(self.object.dec))
                            hdulist.close()

                        hdulist.close()

                    else:
                        etaslog.warning('%s key missing for frame %s. Exclude this frame.'
                                        % (header, os.path.basename(oframes[n]['path'])))
                        failed_frames.append(n)

        for n in failed_frames:
            shutil.move(oframes[n]['path'], oframes[n]['path'] + '.rem')
            del oframes[n]

        #######################
        # Create Master Frame #
        #######################

        self.masterfilepath = os.path.join(self.masterdir, 'master.fits')
        etaslog.info('Creating master frame: %s ' % os.path.join(self.pars.wdir, self.masterfilepath))

        if os.path.isfile(os.path.join(self.pars.wdir, self.masterfilepath)):
            etaslog.info('Removing previous version of master frame')
            os.remove(os.path.join(self.pars.wdir, self.masterfilepath))


        # Get the total number of plate solved frames
        npltsolved = 0
        for n in oframes:
            if oframes[n]['PLTSOLVED']:
                npltsolved += 1

        if npltsolved >= self.params['masternstack']:
            # masternstack is a parameter that determines the number of frames to stack to create a master frame
            # and npltsolved is the total number of plate solved frames.
            # here we check if there are at least masternstack well aligned plate solved frames that can be stacked

            etaslog.info('Stack frames using plate solution of first %i frames' % self.params['masternstack'])

            # get average RA and DEC of frames
            crval1 = []
            crval2 = []
            for n in oframes:
                if oframes[n]['PLTSOLVED']:
                    crval1.append(oframes[n]['CRVAL1'])
                    crval2.append(oframes[n]['CRVAL2'])

            crval1avg = np.average(crval1)
            crval2avg = np.average(crval2)
            etaslog.info('Average RA is %s' % deg_to_sex(ra=crval1avg))
            etaslog.info('Average Dec is %s' % deg_to_sex(dec=crval2avg))

            etaslog.info('Selecting well aligned plate solved frames. Max offset is %.2f arcmin' %
                         self.params['framemaxoffset'])

            # prepare a file with a list of well aligned plate solved frames. Used by iraf imcombine
            stacklist = os.path.join(self.pars.wdir, self.masterdir, 'stacklist')
            tmpf = open(stacklist, 'w')

            # loop the platemsolved frames and select only those that are not offset by more than self.pars.framemaxoffset
            ra_shift = crval1avg # average ra
            dec_shift = crval2avg # average dec
            nstack = 0
            for n in oframes:
                if oframes[n]['PLTSOLVED']:  # only plate solved frames
                    ra = float(oframes[n]['CRVAL1'])  # frame ra
                    dec = float(oframes[n]['CRVAL2'])  # frame dec
                    # check offset
                    offset = math.sqrt(math.fabs(ra-ra_shift)**2+math.fabs(dec-dec_shift)**2)*60
                    if offset < self.params['framemaxoffset']:
                            if nstack < self.params['masternstack']:
                                tmpf.write("%s\n" % oframes[n]['path'])  # write path to file
                                nstack += 1
                            else:
                                break  # close loop if we reached the minimum number of frames
                    else:
                        etaslog.warning('Frame %s is offset is %.1f, more than limit (%.1f)'
                                     % (oframes[n]['path'], offset, self.params['framemaxoffset']))
            tmpf.close()

            # nstack is the number of well aligned plate solved frames that can be stacked
            # if nstack < masternstack then we need to align the frames using alipy.align and then stack them using iraf
            # otherwise stack frames with iraf using platesolution (no alignment required)

            if nstack == self.params['masternstack']:

                etaslog.info('There are enough plate solved and aligned frames. Crate master frame using iraf.imcombine')

                if os.path.isfile(os.path.join(self.pars.wdir, self.masterfilepath)):
                    # remove the master file in case it exists (iraf cannot overwrite!)
                    etaslog.warning('The master frame already exists (%s). Remove it.'
                                    % os.path.join(self.pars.wdir, self.masterfilepath))
                    os.remove(os.path.join(self.pars.wdir, self.masterfilepath))

                iraf.images()
                iraf.imcombine(input='@%s' % stacklist,
                               output=os.path.join(self.pars.wdir, self.masterfilepath),
                               scale=self.pars.iraf_imcombine['scale'],
                               weight=self.pars.iraf_imcombine['weight'],
                               combine=self.pars.iraf_imcombine['combine'],
                               offsets='world',
                               Stdout=self.pars.iraf_imcombine['Stdout'])
                if os.path.isfile(os.path.join(self.pars.wdir, self.masterfilepath)):
                    etaslog.info('Master frame correctly created: %s' % os.path.join(self.pars.wdir, self.masterfilepath))
                else:
                    etaslog.warning('Master frame could not be created. Attempt to align the frames and stack them')

        if not os.path.isfile(os.path.join(self.pars.wdir, self.masterfilepath)):

            if npltsolved < self.params['masternstack']:

                etaslog.warning('There are not enough plate solved and well aligned frames to create a master frame'
                                'using iraf.imcombine and a plate solution.')

            etaslog.info('Select %i frames and align them using alipy ' % self.params['masternstack'])

            alignmasterdir = os.path.join(self.pars.wdir, self.masterdir, 'align')

            if os.path.isdir(alignmasterdir):
                # remove align directory if it is already present,
                shutil.rmtree(alignmasterdir)

            try:
                etaslog.info('Creating directory %s' % alignmasterdir)
                os.mkdir(alignmasterdir)
            except Exception, e:
                etaslog.error('Failed to create the folder %s. Cannot continue' % alignmasterdir, exc_info=True)
                sys.exit()

            # copy frames to align to a separate directory
            j = 0
            for n in oframes:
                if j < self.params['masternstack']:
                    filen = os.path.join(self.pars.wdir, self.masterdir, 'align', 'd%s.fits' % str(j).zfill(5))
                    copyfrom = oframes[n]['path']
                    copyto = os.path.join(self.pars.wdir, self.masterdir, 'align', filen)
                    shutil.copy(copyfrom, copyto)
                    etaslog.debug('Copy %s to %s' % (copyfrom, copyto))
                    j += 1
                else:
                    break

            # align frames with with alipy
            etaslog.info('Run alipy to align %i frames' % self.params['masternstack'])
            images_to_align = sorted(glob.glob(os.path.join(self.pars.wdir, self.masterdir, 'align', '*.fits')))
            ref_image = os.path.join(self.pars.wdir, self.masterdir, 'align', 'd00001.fits')
            outputshape = alipy.align.shape(ref_image, verbose=False)
            identifications = alipy.ident.run(ref_image, images_to_align, visu=False, verbose=False)
            for id in identifications:
                if id.ok == True:
                    aligned_dir = os.path.join(self.pars.wdir, self.masterdir, 'align', 'alipy_aligned')
                    # Variant 1, using only scipy and the simple affine transorm :
                    alipy.align.affineremap(id.ukn.filepath, id.trans, shape=outputshape, makepng=True,
                                            outdir=aligned_dir, verbose=False)
                    # Variant 2, using geomap/gregister, correcting also for distortions :
                    # alipy.align.irafalign(id.ukn.filepath, id.uknmatchstars, id.refmatchstars, shape=outputshape,
                    # makepng=False)

            images_to_stack = glob.glob(os.path.join(self.pars.wdir, self.masterdir, 'align', 'alipy_aligned', '*.fits'))
            if len(images_to_stack) < self.params['masternstack_min']:
                etaslog.error('Only %i frames could be aligned (min is %i). Check that the second frame of the series'
                              'is well aligned, as the reference frame for alignment is always the second one.'
                              % (len(images_to_stack), self.params['masternstack_min']))
                sys.exit()

            # now combine the frames using imcombine
            stacklist = os.path.join(self.pars.wdir, self.masterdir, 'stacklist')
            tmpf = open(stacklist, 'w')
            for path in images_to_stack:
                tmpf.write("%s\n" % path)
            tmpf.close()
            etaslog.info('Combine frames with iraf.imcombine')
            iraf.images()
            iraf.imcombine(input='@' + stacklist,
                           output=os.path.join(self.pars.wdir, self.masterfilepath),
                           scale=self.pars.iraf_imcombine['scale'],
                           weight=self.pars.iraf_imcombine['weight'],
                           combine=self.pars.iraf_imcombine['combine'],
                           Stdout=self.pars.iraf_imcombine['Stdout'])

            etaslog.info('Delete intermediate files')
            remove_dir = os.path.join(self.pars.wdir, self.masterdir, 'align')
            shutil.rmtree(remove_dir)

            # run sextractor on master frame in order to find plate solution with SCAMP
            etaslog.info('Run sextractor on the master frame')
            sexrun = '%s ' % self.pars.sex
            sexrun += '%s ' % os.path.join(self.pars.wdir, self.masterfilepath)
            sexrun += '-c %s ' % os.path.join(self.pars.aconfdir, 'master.sex')
            sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.pars.aconfdir, 'scamp.param')
            sexrun += '-FILTER_NAME %s ' % os.path.join(self.pars.aconfdir, 'default.conv')
            sexrun += '-CATALOG_TYPE FITS_LDAC '
            sexrun += '-DETECT_THRESH 10 '
            sexrun += '-CATALOG_NAME %s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.ps')
            p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output = p.communicate()[1]
            f = output.find('Objects: detected')
            sextracted = int(output[f+40:f+44])
            etaslog.info('%i stars extracted from master frame' % sextracted)

            if sextracted < self.params['minstars']:
                etaslog.erro('There are not enough identified stars in the master frame (min is %i). '
                             'Cannot continue...' % self.params['minstars'])
                sys.exit()

            etaslog.info('Plate solve master frame using SCAMP')

            # get CRVAL1 and CRVAL2 for plate solution
            etaslog.info('Use object coordinates from object class to estimate the catalogue FOV '
                         'for plate solution')
            crval1 = self.object.ra_d
            crval2 = self.object.dec_d

            # create .ahead file with approximate plate solution
            etaslog.info('Create .ahead file with an approximate plate solution')
            crpix1 = self.img_x / 2.0
            crpix2 = self.img_y / 2.0
            aheadf = os.path.join(self.pars.wdir, self.masterdir, 'master.ahead')
            f = open(aheadf, 'w')
            f.write('EQUINOX = 2000.0 \n')
            f.write('EPOCH   = 2000.0 \n')
            f.write('PA      = 0 \n')
            f.write("CTYPE1  = 'RA---TAN'          /\n")
            f.write('CRVAL1  = ' + str(crval1) + '\n')
            f.write('CRPIX1  = ' + str(crpix1) + '\n')
            f.write('CROTA1  = 0 \n')
            f.write("CTYPE2  = 'DEC--TAN'          /\n")
            f.write('CRVAL2  = ' + str(crval2) + ' \n')
            f.write('CRPIX2  = ' + str(crpix2) + '\n')
            f.write('CD1_1   = ' + str(self.params['telescope']['CD1_1']) + '\n')
            f.write('CD1_2   = ' + str(self.params['telescope']['CD1_2']) + ' \n')
            f.write('CD2_1   = ' + str(self.params['telescope']['CD2_1']) + ' \n')
            f.write('CD2_2   = ' + str(self.params['telescope']['CD2_2']) + '')
            f.close()

            # run Scamp to find plate solution
            if not internet_on():
                etaslog.error('Cannot connect to the Internet: SCAMP needs Internet access to download the star'
                              'catalogues. Check your connection and try again')
                sys.exit()

            etaslog.info('Run scamp using online catalogue %s to get plate solution of master frame'
                         % self.params['onlinecat'])

            scamprun = '%s ' % self.pars.scamp
            scamprun += '%s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.ps')
            scamprun += '-c %s ' % os.path.join(self.pars.aconfdir, 'default.scamp')
            scamprun += '-ASTREF_CATALOG %s ' % self.params['onlinecat']
            print scamprun

            p = subprocess.Popen(scamprun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            output = p.communicate()[1]

            # @todo: check output --> need to be sure the frame is plate solved!

            etaslog.info('Run missfits to append plate solution to master frame')
            missrun = self.pars.missfits + " " + os.path.join(self.pars.wdir, self.masterfilepath) + "  \
                      -c " + os.path.join(self.pars.aconfdir, 'default.missfits')
            p = subprocess.Popen(missrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]

            etaslog.warning('Be careful, it is not certain that the plate solution was successful! If something goes'
                            'wrong later, one possible reason is that the master frame is not correctly plate solved')

        ######################################
        # Generate master catalogue of stars #
        ######################################

        # extract star catalogue from sextractor output

        etaslog.info('Run sextractor on master frame to obtain coordinates of stars in WCS.')
        etaslog.info('Detection treshold is %i' % self.params['detect_tresh'])

        # run sextractor to obtain coordinates of stars in WCS
        # master.cat contains the WCS coordinates of the stars to perform photometry on

        sexrun = '%s ' % self.pars.sex
        sexrun += '%s ' % os.path.join(self.pars.wdir, self.masterfilepath)
        sexrun += '-c %s ' % os.path.join(self.pars.aconfdir, 'master.sex')
        sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.pars.aconfdir, 'master.param')
        sexrun += '-FILTER_NAME %s ' % os.path.join(self.pars.aconfdir, 'default.conv')
        sexrun += '-CATALOG_TYPE ASCII '
        sexrun += '-DETECT_THRESH %i ' % self.params['detect_tresh']
        sexrun += '-CATALOG_NAME %s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.cat')
        p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]

        # load catalogue of stars in array using numpy
        mastercat_orig = np.loadtxt(os.path.join(self.pars.wdir, self.masterdir,  'master.cat'))
        mastercat_sorted = mastercat_orig[mastercat_orig[:,4].argsort()]  # sort by ISO FLUX (4th column)

        # remove edge star
        etaslog.info('Remove edge stars. Edge is %i pixels, frame is %ix%i'
                     % (self.params['edgelimit'], self.img_x, self.img_y))

        mastercat_edge = []
        j = 0
        for line in mastercat_sorted:
            if (line[2] < (self.img_x - self.params['edgelimit'])) and (line[2] > self.params['edgelimit']) and \
               (line[3] < (self.img_y - self.params['edgelimit'])) and (line[3] > self.params['edgelimit']):  # y axis
                mastercat_edge.append([line[0], line[1], line[2], line[3], line[4], line[5]])
                j += 1
            else:
                etaslog.info('Exclude edge star at x: %.1f y: %.1f ' % (line[2], line[3]))

        # if the stars in the catalogue are more than the maximum allowed (self.params['maxstars']), take only the last
        # self.params['maxstars'] stars in the list (the list is already sorted by flux in ascending order)
        if j > self.params['maxstars']:
            etaslog.info("%i stars found (max %i). Exclude faintest stars. " % (j, self.params['maxstars']))
            mastercat_final = []
            for n in range(self.params['maxstars']):
                mastercat_final.append(mastercat_edge[j-n-1])
        else:
            mastercat_final = mastercat_edge

        # check for the presence of the target star in the final master catalogue
        if isinstance(self.object, Target):
            etaslog.info('Looking for target star...')
            k=0
            for line in mastercat_final:
                xshift = float(self.object.dec_d - line[1])
                yshift = float(self.object.ra_d - line[0])
                shift = math.sqrt(xshift**2+yshift**2)
                if shift < 0.005:
                    self.targetid = k
                    break
                k += 1
            try:
                self.targetid
            except:
                # the planet was not found in the final list, check if it's been excluded
                etaslog.info('The target was not found in the final catalogue. Check for fainter stars...')
                for line in mastercat_edge:
                    if math.fabs(float(self.object.dec_d - line[1])) < 0.005 and \
                       math.fabs(float(self.object.ra_d - line[0])) < 0.005:
                        # the target is in the original master catalogue of stars. It was too faint to be included
                        # in the catalogue. Add it now and exclude faintest star.
                        del mastercat_final[len(mastercat_final)-1]
                        mastercat_final.append(line)
                        self.targetid = len(mastercat_final)-1
                try:
                    self.targetid
                except:
                    etaslog.info('The target was not found! It is not possible to continue with the photometry')
                    sys.exit()
            etaslog.info('Target detected! %s is star number %i' % (self.object.name, self.targetid))
        else:
            etaslog.info('The current object class does not have a specific target')

        # rewrite an iraf.tvmark version of the catalogue
        mastercat_tvmark = []
        for line in mastercat_final:
            mastercat_tvmark.append([line[2], line[3]])

        # convert catalogues to numpy arrays
        mastercat_final = np.array(mastercat_final)
        mastercat_tvmark = np.array(mastercat_tvmark)

        # write output to files
        filename = os.path.join(self.pars.wdir, self.masterdir, 'master.final')
        etaslog.info('Output file with catalogue stars: %s' % filename)
        np.savetxt(filename, mastercat_final)
        filename = os.path.join(self.pars.wdir, self.masterdir, 'master.tvmark')
        etaslog.info('Output file with catalogue stars: %s' % filename)
        np.savetxt(filename, mastercat_tvmark)

        # generate png of master file with catalogue stars
        self.masterpngpath = os.path.join(self.pars.wdir, self.masterdir, 'master.png')
        etaslog.info('Generate PNG of master file with catalogue stars:  %s' % (self.masterpngpath))
        coords = []
        n = 0
        for line in mastercat_final:
            coord = {}
            coord['x'] = line[2]
            coord['y'] = line[3]
            coord['id'] = str(n)
            coords.append(coord)
            n += 1
        fits_to_png(os.path.join(self.pars.wdir, self.masterdir, 'master.fits'),
                    os.path.join(self.pars.wdir, self.masterdir, 'master.png'), coords)

        ######################
        # Plate solve frames #
        ######################

        etaslog.info('Check that all frames are plate solved. If not, plate solve them, using the stars found in the '
                     'master frame as a reference catalogue.')

        frames = list_frames(os.path.join(self.pars.wdir, self.framesdir))[0]
        platesolve_frames = []
        for n in frames:
            if not frames[n]['PLTSOLVED']:
                platesolve_frames.append(frames[n])

        run_scamp = []
        run_missfits = []
        for frame in platesolve_frames:

            etaslog.info('Frame %s is not plate solved. Create SExtractor catalog of stars. ' % frame['path'])
            frame_basename = os.path.splitext(os.path.basename(frame['path']))[0]
            cataloguename = os.path.join(self.pars.wdir, self.framesdir, frame_basename) + '.scampcat'
            sexrun = '%s ' % self.pars.sex
            sexrun += '%s ' % frame['path']
            sexrun += '-c %s ' % os.path.join(self.pars.aconfdir, 'platesolve.sex')
            sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.pars.aconfdir, 'scamp.param')
            sexrun += '-FILTER_NAME %s ' % os.path.join(self.pars.aconfdir, 'default.conv')
            sexrun += '-CATALOG_TYPE FITS_LDAC '
            sexrun += '-DETECT_THRESH 10 '
            #sexrun += '-DETECT_THRESH %i ' % self.params['detect_tresh']
            sexrun += '-CATALOG_NAME %s ' % cataloguename
            p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

            # determine number of detected objects from SExtractor output and exclude frames with low detections
            output = p.communicate()[1]
            f = output.find('Objects: detected')
            sextracted = int(output[f+40:f+44])
            if sextracted < self.params['pssexlimit']:

                etaslog.warning('Cannot find enough stars in this frame. Exclude it: %s'
                                % os.path.basename(frame['path']))
                shutil.move(frame['path'], frame['path'] + '.rem.fits')
            else:
                etaslog.info('%i stars detected, a plate solution will be attempted on %s. '
                                % (sextracted, os.path.basename(frame['path'])))

                # get crval1/2 and crpix1.2
                outputshape = alipy.align.shape(frame['path'], verbose=False)
                crval1 = sex_to_deg(frame['OBJCTRA'], 'ra')
                crval2 = sex_to_deg(frame['OBJCTDEC'], 'dec')
                img_x = outputshape[0]  # image x size (pixels)
                img_y = outputshape[1]  # image y size (pixels)
                crpix1 = img_x / 2.0
                crpix2 = img_y / 2.0

                # create .ahead file (use with scamp)
                aheadf = os.path.join(self.pars.wdir, self.framesdir, frame_basename) + '.ahead'
                f = open(aheadf, 'w')
                f.write('EQUINOX = 2000.0 \n')  # assume J2000
                f.write('EPOCH   = 2000.0 \n')
                f.write('PA      = 0 \n')
                f.write("CTYPE1  = 'RA---TAN'          /\n")
                f.write('CRVAL1  = ' + str(crval1) + '\n')
                f.write('CRPIX1  = ' + str(crpix1) + '\n')
                f.write('CROTA1  = 0 \n')
                f.write("CTYPE2  = 'DEC--TAN'          /\n")
                f.write('CRVAL2  = ' + str(crval2) + ' \n')
                f.write('CRPIX2  = ' + str(crpix2) + '\n')
                f.write('CD1_1   = ' + str(self.params['telescope']['CD1_1']) + '\n')
                f.write('CD1_2   = ' + str(self.params['telescope']['CD1_2']) + ' \n')
                f.write('CD2_1   = ' + str(self.params['telescope']['CD2_1']) + ' \n')
                f.write('CD2_2   = ' + str(self.params['telescope']['CD2_2']) + '')
                f.close()
                run_scamp.append(cataloguename)
                run_missfits.append(frame['path'])

        if not len(run_scamp) > 0:
            etaslog.info('There are no frames to plate solve')
        else:
            # create reference catalogue for plate solution from master frame (to use as an input catalogue for scamp)
            etaslog.info('Create reference catalogue for plate solution from master frame. This will be used as an '
                         'input catalogue by scamp')

            sexrun = '%s ' % self.pars.sex
            sexrun += '%s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.fits')
            sexrun += '-c %s ' % os.path.join(self.pars.aconfdir, 'master.sex')
            sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.pars.aconfdir, 'scamp.param')
            sexrun += '-FILTER_NAME %s ' % os.path.join(self.pars.aconfdir, 'default.conv')
            sexrun += '-CATALOG_TYPE FITS_LDAC '
            sexrun += '-DETECT_THRESH 60 '
            sexrun += '-CATALOG_NAME %s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.platesolve')
            p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]

            # run scamp
            for file in run_scamp:     # run scamp
                refcatpath = os.path.join(self.pars.wdir, self.masterdir, 'master.platesolve')
                refcat = glob.glob(refcatpath)
                if refcat:
                    etaslog.info('Run scamp on frame %s' % os.path.basename(file))

                    scamprun = '%s ' % self.pars.scamp
                    scamprun += '%s ' % file
                    scamprun += '-c %s ' % os.path.join(self.pars.aconfdir, 'default.scamp')
                    scamprun += '-ASTREF_CATALOG FILE %s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.platesolve')
                    scamprun += '-ASTREFCAT_NAME %s ' % os.path.join(self.pars.wdir, self.masterdir, 'master.platesolve')
                    p = subprocess.Popen(scamprun, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, shell=True).communicate()[0]
                else:
                    etaslog.error('Cannot find reference catalogue: %s' % refcatpath)
                    sys.exit()

            #run missfits
            for file in run_missfits:     # run missfits
                etaslog.info('Run missfits on frame %s' % os.path.basename(file))

                missrun = '%s ' % self.pars.missfits
                missrun += '%s ' % file
                missrun += '-c %s ' % os.path.join(self.pars.aconfdir, 'default.missfits')
                p = subprocess.Popen(missrun, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE, shell=True).communicate()[0]

        # Save the final list of frames
        etaslog.info('Obtain the final list of plate solved frames')
        self.frames = list_frames(os.path.join(self.pars.wdir, self.framesdir))[0]
        self.nobs = len(self.frames)

        # get coordinates of stars
        etaslog.info('Get coordinate of stars from the master catalogue of stars')
        coord = np.loadtxt(os.path.join(self.pars.wdir, self.masterdir, 'master.final'), unpack=True)
        self.ra = coord[0]
        self.dec = coord[1]
        self.nstar = len(coord[0])

        ############################################################
        # Create SExtractor catalogue for all frames and get FWHM  #
        ############################################################

        etaslog.info('Create final SExtractor catalogue for all frames')

        for frame in self.frames:
            basename = os.path.splitext(os.path.basename(self.frames[frame]['path']))[0]
            cataloguename = os.path.join(self.pars.wdir, self.framesdir, basename) + '.cat'

            sexrun = '%s ' % self.pars.sex
            sexrun += '%s ' % self.frames[frame]['path']
            sexrun += '-c %s ' % os.path.join(self.pars.aconfdir, 'master.sex')
            sexrun += '-PARAMETERS_NAME %s ' % os.path.join(self.pars.aconfdir, 'default.param')
            sexrun += '-FILTER_NAME %s ' % os.path.join(self.pars.aconfdir, 'default.conv')
            sexrun += '-CATALOG_TYPE ASCII '
            sexrun += '-DETECT_THRESH 60 '
            sexrun += '-CATALOG_NAME %s ' % cataloguename

            p = subprocess.Popen(sexrun, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]

        etaslog.info('Create SExtractor catalogue for all frames')

        mastercat = np.loadtxt(os.path.join(self.pars.wdir, self.masterdir, 'master.final'))

        self.fwhm = np.empty((self.nobs, self.nstar))
        i = 0
        for frame in self.frames:
            basename = os.path.splitext(os.path.basename(self.frames[frame]['path']))[0]
            cataloguename = os.path.join(self.pars.wdir, self.framesdir, basename) + '.cat'
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
                if ck == False:
                    self.fwhm[i][k] = np.NaN
                    k += 1
            i += 1

        avg1 = stats.stats.nanmean(self.fwhm, 0)
        avg2 = stats.stats.nanmean(avg1)
        std = stats.stats.nanstd(avg1)
        self.avgfwhm = avg2
        self.stdfwhm = std

        etaslog.info('Average FWHM is %.1f +/- %.1f' % (self.avgfwhm, self.stdfwhm))

        ###############################################################
        # Calculate LMST and AIRMASS values for all frames using IRAF #
        ###############################################################

        etaslog.info('Updating airmass values in fits headers')

        iraf.noao()
        iraf.astutil()
        iraf.observatory.setParam('obsid', 'obspars')
        iraf.observatory.setParam('command', 'set')
        iraf.observatory.setParam('observatory', 'obspars')
        iraf.observatory.setParam('name', self.params['observatory']['name'])
        iraf.observatory.setParam('longitu', str(self.params['observatory']['longitude']))
        iraf.observatory.setParam('altitud', str(self.params['observatory']['elevation']))
        iraf.observatory.setParam('latitud', str(self.params['observatory']['latitude']))
        iraf.observatory.setParam('override', 'obspars')

        f1 = open(os.path.join(self.pars.wdir, self.framesdir, 'imagelist'), 'w')
        for n in self.frames:
            f1.write(self.frames[n]['path'] + '\n')
        f1.close()

        etaslog.info('Calculating LMST')
        f = open(os.path.join(self.pars.wdir, self.framesdir, 'lmst-calc'), 'w')
        f.write("observat = 'obspars'\n")
        f.write("lmst = mst (@'DATE-OBS', @'TIME-OBS', obsdb (observat, 'longitude'))\n")
        f.write("quit")
        f.close()
        iraf.asthedit(images='@'+os.path.join(self.pars.wdir, self.framesdir, 'imagelist'),
                      commands=os.path.join(self.pars.wdir, self.framesdir, 'lmst-calc'),
                      Stdout=1, verbose='no')

        etaslog.info('Calculating airmass (an extra header EAIRMASS is added to the fits files)')
        iraf.setairmass(images='@'+os.path.join(self.pars.wdir, self.framesdir, 'imagelist'), observa='obspars', \
                        intype='beginning', outtype='effective', st='lmst', \
                        ra='OBJCTRA', dec='OBJCTDEC', equinox='EPOCH', ut='TIME-OBS', \
                        date='DATE-OBS', exposur='EXPTIME', airmass='EAIRMASS', utmiddl='utmiddle', \
                        scale='750.', update='yes', Stdout=1)

        # update self.frames with new AIRMASS key
        self.frames = list_frames(os.path.join(self.pars.wdir, self.framesdir))[0]
        self.airmass = np.empty((self.nobs))
        self.exptime = np.empty((self.nobs))
        self.path = {}

        for n in self.frames:
            self.airmass[n] = self.frames[n]['EAIRMASS']
            self.exptime[n] = self.frames[n]['EXPTIME']
            self.path[n] = self.frames[n]['path']

        ########################################
        # Calculate BJD and HJD for each frame #
        #######################################

        etaslog.info('Calculate BJD and HJD for each frame')

        bjd_filename = os.path.join(self.pars.wdir, self.framesdir, 'bjdtbd')
        hjd_filename = os.path.join(self.pars.wdir, self.framesdir, 'hjd')
        utc_filename = os.path.join(self.pars.wdir, self.framesdir, 'utc')

        # load UTC for each frame from header
        import datetime
        f = open(utc_filename, 'w')
        for n in self.frames:
            ts = float(self.frames[n]['timestamp']) + float(self.frames[n]['EXPTIME'])/2
            utc = datetime.datetime.fromtimestamp(int(ts)).strftime('%Y-%m-%dT%H:%M:%S')
            f.write('%s %i %i \n' % (utc, 0, 0))
        f.close()

        # ra and dec of target, used to calculate BJD
        ra = self.object.ra_d
        dec = self.object.dec_d

        # run vartools -converttime, write the BJD values to bjd_filename
        runvartools = self.pars.vartools + " -i " + utc_filename
        runvartools += " -quiet -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 "
        runvartools += "-converttime input jd inputsys-utc output bjd outputsys-tdb radec fix %f %f " % (ra, dec)
        runvartools += "ephemfile %s/de421.bsp " % self.pars.vartoolsdir
        runvartools += "leapsecfile %s/naif0010.tls " % self.pars.vartoolsdir
        runvartools += "planetdatafile %s/pck00010.tpc " % self.pars.vartoolsdir
        runvartools += "coords fix %f %f %f "  % (self.params['observatory']['latitude'],
                                                  self.params['observatory']['longitude'],
                                                  self.params['observatory']['elevation'])
        runvartools += "-o %s " % bjd_filename

        #p = subprocess.Popen(runvartools, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]
        os.system(runvartools)

        self.bjd = np.loadtxt(bjd_filename, unpack=True)[0]

        # run vartools -converttime, write the HJD values to hjd_filename
        runvartools = self.pars.vartools + " -i " + utc_filename
        runvartools += " -quiet -readformat 0 inpututc '%Y-%M-%DT%h:%m:%s' 1 2 3 "
        runvartools += "-converttime input jd inputsys-utc output hjd outputsys-utc radec fix %f %f " % (ra, dec)
        runvartools += "ephemfile %s/de421.bsp " % self.pars.vartoolsdir
        runvartools += "leapsecfile %s/naif0010.tls " % self.pars.vartoolsdir
        runvartools += "planetdatafile %s/pck00010.tpc " % self.pars.vartoolsdir
        runvartools += "coords fix %f %f %f " % (self.params['observatory']['latitude'],
                                                 self.params['observatory']['longitude'],
                                                 self.params['observatory']['elevation'])
        runvartools += "-o %s " % hjd_filename
        p = subprocess.Popen(runvartools, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True).communicate()[0]
        self.hjd = np.loadtxt(hjd_filename, unpack=True)[0]

        #############################################
        # Calculate scintillation noise theoretical #
        #############################################

        # this is the expected scintillation noise for each frame, which will be added in quadrature to the ccdnoise
        # this error might be underestimated. The function scint_noise()  in the class Photometry, will attempt
        # to derive a correcting factor.

        etaslog.info('Calculate theoretical scintillation noise')

        # scintillation noise
        scintwav = (self.params['filter']['midwav']/550.0)**(-7.0/12.0)
        scintalti = math.exp(-self.params['observatory']['elevation']/8000.0)
        scintap = self.params['telescope']['aperture']**(2.0/3.0)
        scintfact = 0.09*scintalti*scintwav/scintap

        self.sigscint = scintfact*(self.airmass**1.75)/np.sqrt(2.*self.exptime)
        self.medsigscint = stats.nanmedian(self.sigscint, 0)

        etaslog.info('Dataset instance correctly initialised')