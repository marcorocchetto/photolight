
import numpy as np
import logging
import os
import shutil

from functions import *
from list_frames import *


class Platesolve(object):

    def __init__(self, dataset, overwrite=True, params=None):

        logging.info('Initialising platesolve object')

        # set some objects and variables
        if params is not None:
            self.params = params
        else:
            self.params = dataset.params

        print dataset.params

        self.dataset = dataset
        self.overwrite = overwrite

    def is_platesolved(self, filepath):

        # check if a file is already plate solved
        header = get_header(filepath)

        keys = ['CTYPE1', 'CRVAL1', 'CRPIX1', 'CDELT1', 'CROTA1', 'CTYPE2', 'CRVAL2', 'CRPIX2', 'CDELT2', 'CROTA2']
        for key in keys:
            if not key in header.keys():
                return False
        return True

    def astrometry(self):

        # plate solve frames using astrometry.net routines
        psoln = 0
        for n in self.dataset.frames:

            filepath = self.dataset.frames[n]['path']

            print self.params

            if self.overwrite or (not self.overwrite and not self.is_platesolved(filepath)):
                logging.info('Plate solve file %s using astrometry.net', os.path.basename(filepath))
                run = self.params.platesolve['solvefield_path']
                if self.params.platesolve['use_coordinates']:
                    header = get_header(filepath)
                    print self.params.header_objctra

                    print self.params

                    ra = get_headerval_from_keywords(header, self.params.header_objctra).replace(' ', ':')
                    dec = get_headerval_from_keywords(header, self.params.header_objctdec).replace(' ', ':')
                    run += ' --ra %s --dec %s --radius %f ' % (ra, dec, self.params.platesolve['search_radius'])
                if self.params.platesolve['flags']:
                    run += ' %s ' %  self.params.platesolve['flags']
                run_file = '%s %s' % (run, filepath)
                os.system(run_file)

                # this file exists if the field is solved
                is_solved = os.path.isfile('%s.solved' % os.path.splitext(filepath)[0])
                # file name of the solved field. Only if is solved exists.
                new_filename = '%s.new' % os.path.splitext(filepath)[0]

                # if frame was not platesolved, try with sextractor source extraction
                if not is_solved:
                    run += ' --use-sextractor '
                    run_file = '%s %s' % (run, filepath)
                    os.system(run_file)
                    # this file exists if the field is solved
                    is_solved = os.path.isfile('%s.solved' % os.path.splitext(filepath)[0])
                    if not is_solved:
                        logging.error('Could not plate solve  the frame %s' % os.path.basename(filepath))
                        continue

                # remove aux files
                exts = ['-objs.png', '.axy', '-indx.xyls', '-objs.png', '.corr', '.rdls', '.match', '.wcs', '.solved']
                for ext in exts:
                    rmfile = '%s%s' % (os.path.splitext(filepath)[0], ext)
                    if os.path.isfile(rmfile):
                        os.remove(rmfile)

                # save the new file with WCS
                os.remove(filepath)
                shutil.move('%s.new' % os.path.splitext(filepath)[0], filepath)

                logging.info('Frame %s correctly plate solved.' % os.path.basename(filepath))

                psoln += 1

        logging.info('%i frames have been plate solved using astrometry.net' % psoln)