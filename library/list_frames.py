
def list_frames(dir, exclude=[]):

    import os
    import time
    import pyfits


    platekeys = ['CTYPE1', 'CRVAL1', 'CRPIX1', 'CTYPE2', 'CRVAL2', 'CRPIX2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    stkeys = ['DATE-OBS', 'EXPTIME', 'AIRMASS', 'CRVAL1', 'CRVAL2']

    valid_frames = {}
    k = 0

    for path, subdirs, files in os.walk(dir):   # Careful: this also consider subdirectories!
        files.sort() # sort by filename
        for name in files:
            if not name in exclude:
                frame = os.path.join(path, name)
                valid_frames[k] = {}
                valid_frames[k]['path'] = frame
                header = pyfits.getheader(frame)    # Read header

                # store all headers specified from input
                for key in stkeys:
                    if key in header:
                        valid_frames[k][key] = header[key]
                    else:
                        print 'There are missing header keys' # @todo convert to log
                        return False # return if one frame does not have a

                # DATE-OBS is converted into a UNIX timestamp and a key named "timestamp" is added
                if 'DATE-OBS' in header:
                    obsdate = header["DATE-OBS"]
                    dtime = obsdate.split(".")
                    try:
                        #convert DATEOBS in unix timestamp
                        timestamp = time.mktime(time.strptime(dtime[0], "%Y-%m-%dT%H:%M:%S"))
                    except:
                        #convert DATEOBS in unix timestamp
                        obsdate = header["DATE"]
                        dtime = obsdate.split(".")
                        try:
                            timestamp = time.mktime(time.strptime(dtime[0], "%Y-%m-%dT%H:%M:%S"))
                        except:
                            print 'Error reading DATE-OBS' # @todo convert to log
                            return False
                    valid_frames[k]['timestamp'] = timestamp
                    valid_frames[k]['DATE-OBS'] = obsdate

                # add an extra key called "plate_solved" to the output dictionary for this frame.
                # value is True if the frame is plate solved.
                for key in platekeys:
                    # it checks if the headers specified in platekeys are present
                    if not key in header:
                        valid_frames[k]['plate_solved'] = False  # the frame is not platesolved
                        break

                if not 'plate_solved' in valid_frames[k]:
                    valid_frames[k]['plate_solved'] = True  # the frame is platesolved

                k += 1

    return valid_frames

