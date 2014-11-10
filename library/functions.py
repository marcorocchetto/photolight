'''
ETAS Functions
Last update: 14 Aug 2014
'''


def internet_on():

    # check internet connection

    import urllib2

    try:
        response=urllib2.urlopen('http://74.125.228.100', timeout=1)
        return True
    except urllib2.URLError as err:
        pass

    return False

def converter(x):

    #convert INDEF values with -9999. Used by loadtxt

    import numpy as np

    try:
        if x == 'INDEF' or float(x) == 0. or float(x) < 0.:
            return np.NaN
        else:
            return float(x)
    except:
        return np.NaN


def fits_to_png(inpath, outpath, coords=''):

    import numpy as np
    import img_scale
    import pylab
    import pyfits

    fn = inpath
    sig_fract = 5.0
    percent_fract = 0.01
    hdulist = pyfits.open(fn)
    img_header = hdulist[0].header
    img_data_raw = hdulist[0].data
    hdulist.close()
    width=img_data_raw.shape[0]
    height=img_data_raw.shape[1]
    img_data_raw = np.array(img_data_raw, dtype=float)
    sky, num_iter = img_scale.sky_mean_sig_clip(img_data_raw, sig_fract, percent_fract, max_iter=10)
    img_data = img_data_raw - sky
    min_val = 10.0
    
    new_img = img_scale.log(img_data, scale_min = min_val)
    pylab.imshow(new_img, interpolation='nearest', origin='lower',cmap=pylab.cm.gray)

    if coords != '':
        for coord in coords:
            circle=pylab.Circle((coord['x'],coord['y']),radius=10)
            circle.set_edgecolor( 'red' )
            circle.set_facecolor( 'none' )  # "none" not None
            circle.set_alpha( 0 )
            pylab.gca().add_artist(circle)
            pylab.text(coord['x']+15, coord['y']-20, coord['id'], color='red')
    pylab.axis('off')
    pylab.savefig(outpath)
    pylab.clf()
      
# calculate calendar & julian dates
def calc_date(dir):

    import list_frames
    import datetime

    print 'Calculate mean calendar and Julian dates of files in %s ' % dir

    ts = 0

    frames = list_frames(dir, out='keys',keys=['DATE-OBS'])[0]
    for n in frames:
        ts = ts + frames[n]['timestamp']
        
    tsmid = ts/len(frames)
    cdate = datetime.datetime.fromtimestamp(int(tsmid)).strftime('%Y-%m-%d')
    jdate = float((tsmid/86400.0)+2440587.5)
    print 'Calendar date is %s, Julian date is %s' % (cdate, jdate)
    return cdate, jdate

def sex_to_deg(value, coord):

    # @todo: validate coord!

    import string
    #convert OBJCTRA and OBJCDEC from sexagesimal to decimal degrees
    #assume RA: hh mm ss.sss
    #assume Dec: +/-deg arcmin arcsec.ss
    
    if coord.upper() == 'RA':
        ra=string.split(value, ":")
        if len(ra) < 2: ra = string.split(value, " ")
        hh=float(ra[0])*15
        mm=(float(ra[1])/60)*15
        ss=(float(ra[2])/3600)*15
        return hh+mm+ss
    elif coord.upper() == 'DEC':
        dec=string.split(value, ":")
        if len(dec) < 2: dec = string.split(value, " ")
        hh=abs(float(dec[0]))
        mm=float(dec[1])/60
        ss=float(dec[2])/3600
        if float(dec[0]) < 0:
            return (hh+mm+ss)*(-1)
        else:
            return (hh+mm+ss)

def deg_to_sex(ra='', dec='', round=False):

    RA, DEC = '', ''

    if dec:
        deg = int(dec)
        decM = abs(int((dec-deg)*60))
        if round:
            decS = int((abs((dec-deg)*60)-decM)*60)
        else:
            decS = (abs((dec-deg)*60)-decM)*60
            DEC = '{0:+} {1} {2}'.format(deg, decM, decS)

    if ra:
        raH = int(ra/15)
        raM = int(((ra/15)-raH)*60)
        if round:
            raS = int(((((ra/15)-raH)*60)-raM)*60)
        else:
            raS = ((((ra/15)-raH)*60)-raM)*60
            RA = '{0} {1} {2}'.format(raH, raM, raS)

    if ra and dec:
        return RA, DEC
    else:
        return RA or DEC


def jd_to_date(date):
    import datetime
    ts = float((date - 2440587.5)*86400.0)
    cdate = datetime.datetime.fromtimestamp(int(ts)).strftime('%Y %b %d')            
    return cdate


def round_base(x, base=5):
    return int(base * round(float(x)/base))
