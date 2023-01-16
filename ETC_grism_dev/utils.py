"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *

ETC_GRISM_HOME = "/Users/gael/Documents/PostDoc/Halifax/Work/CASTOR/forecastor/castor_etc/ETC_grism/ETC_grism_dev/"

def mag_to_fnu(mag):
    """
    Converts AB mag to micro Jansky flux densities
    """

    return 10**(23-(mag+48.6)/2.5) * 1e6


def fnu_to_mag(fnu):
    """
    Converts micro Jansky flux densities to AB mag
    """

    return (23 - np.log10(fnu / 1e6)) * 2.5 - 48.6


def fnu_to_flam(fnu, efflamb):
    """
    Converts micro Jansky flux densities to f_lambda flux densities
    """

    #convert mJy to f_lambda
    kste = 10**-29*2.9979*10**18
    f_lambda = fnu * kste / efflamb**2
#    #f_lambda errors
#    deriv_efflamb = - 2. / efflambs**3
#    efflamb_err = 0 #(snr in flambda same as in f_mJy if assume efflamb_err = 0)
#    f_lambda_err = np.sqrt((f_mJy_err * kste / efflambs**2)**2 + (f_mJy * kste * deriv_efflamb * efflamb_err)**2)

    return f_lambda


def flam_to_fnu(f_lambda, efflamb):
    """
    Converts f_lambda flux densities to micro Jansky flux densities
    """

    #convert mJy to f_lambda
    kste = 10**-29*2.9979*10**18
    fnu = f_lambda * efflamb**2 / kste

    return fnu


def mag_to_flam(mag, efflamb):
    """
    Converts AB mag to f_lambda flux densities
    """

    fnu = mag_to_fnu(mag)
    f_lambda = fnu_to_flam(fnu, efflamb)

    return f_lambda


def flam_to_mag(f_lambda, efflamb):
    """
    Converts f_lambda flux densities to AB mag
    """

    fnu = flam_to_fnu(f_lambda, efflamb)
    mag = fnu_to_mag(fnu)

    return mag


def gaussian(x, m, sigma):
    """ Return the normalized Gaussian with standard deviation sigma. """
    c = np.sqrt(2 * np.pi)
    return np.exp(-0.5 * ((x-m) / sigma)**2) / sigma / c

