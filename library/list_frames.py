from .. import pars

def list_frames(dir, out='full', exclude=[], keys=[]):

    """
    Create a list of frame headers:
    each element of the list is a dictionary with keys corresponding to the FITS header required.
    each dictionary also contains the path of the fits file

    Input:

        dir     Source directory, where FITS files are stored
        out     If set to 'full' all default headers are stored. Note that the default headers are defined in the
                pars.py file in the variables platekeys, stkeys and extrakeys.
                If set to 'keys' then specify headers in the keys input variable as a list. Note that the header
                must be part of the list of headers defined in pars.py
        exclude Exclude specific files
        keys    List of headers if out='keys'

        Extrakeys are also added to each dictionary:

            'path', contains the path of the file

            If out='full' or out='keys' and PLTSOLVED is in keys, an extra key is addted to the dictionary
            called PLTSOLVED. True if the frame is plate solved, False if it is not

            If out='full' or out='keys' and PLTSOLVED is in keys, an extra key is addted to the dictionary
            called timestamp. It includes the UNIX timestamp of the DATE-OBS header.

    Output: tuple(validFrames,failedFrames)

    """
    import os
    import time
    import pyfits

    # these are the headers that can be read. All other headers will not be read!
    allkeys = pars.platekeys + pars.stkeys + pars.extrakeys

    failedFrames = []
    validFrames = {}
    k = 0

    for path, subdirs, files in os.walk(dir):   # Careful: this also consider subdirectories!
        files.sort() # sort by filename
        for name in files:
            if not name in exclude and os.path.splitext(name)[1] in pars.extlist:
                frame = os.path.join(path, name)
                validFrames[k] = {}
                validFrames[k]['path'] = frame

                if out == 'full' or out == 'keys':

                    header = pyfits.getheader(frame)    # Read header

                    # store all headers specified from input
                    for key in allkeys:
                        if out == 'full' or key in keys:
                            if key in header:
                                validFrames[k][key] = header[key]
                            else:
                                # skip header
                                continue

                    # DATE-OBS is converted into a UNIX timestamp and a key named "timestamp" is added
                    # to the output dictionary for this frame
                    if out == 'full' or 'DATE-OBS' in keys:
                        if 'DATE-OBS' in header:
                            obsdate = header["DATE-OBS"]
                            dtime = obsdate.split(".")
                            try:
                                timestamp = time.mktime(time.strptime(dtime[0], "%Y-%m-%dT%H:%M:%S")) #convert DATEOBS in unix timestamp
                            except:
                                obsdate = header["DATE"]
                                dtime = obsdate.split(".")
                                try: timestamp = time.mktime(time.strptime(dtime[0], "%Y-%m-%dT%H:%M:%S")) #convert DATEOBS in unix timestamp
                                except:
                                    print 'Error reading DATE-OBS'
                                    continue
                            validFrames[k]['timestamp'] = timestamp
                            validFrames[k]['DATE-OBS'] = obsdate

                    # add an extra key called "PLTSOLVED" to the output dictionary for this frame.
                    # value is True if the frame is platesolved.
                    if out == 'full' or 'PLTSOLVED' in keys:
                        for key in pars.platekeys:
                            # it checks if the headers specified in pars.platekeys are present
                            if not key in header:
                                validFrames[k]['PLTSOLVED'] = False  # the frame is not platesolved
                                break
                        if not 'PLTSOLVED' in validFrames[k]:
                            validFrames[k]['PLTSOLVED'] = True  # the frame is platesolved

                k += 1
    return validFrames, failedFrames
