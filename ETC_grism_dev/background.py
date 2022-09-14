"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *
from . import utils


# ------------------------------- BACKGROUND NOISE VALUES ------------------------------ #

# The flux of the geocoronal emission line [O II] 2471A.
# See <https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds
GEOCORONAL_FLUX_HIGH = 3.0e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_FLUX_AVG = 1.5e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_FLUX_LOW = 7.5e-17  # erg/cm^2/s/arcsec^2
GEOCORONAL_WAVELENGTH = 2471  # angstroms
GEOCORONAL_LINEWIDTH = 0.023  # angstroms

# Average sky background in each passband (Earthshine + zodiacal light)
#SKY_BACKGROUND = {"uv": 26.08, "u": 23.74, "g": 22.60}  # AB mag/arcsec^2
sky_background_mag = {"uv": 27.7831, "u": 24.2583, "g": 22.6063}  # AB mag/arcsec^2
sky_background_flam = {"uv": 1.7186261475789375e-19, "u": 1.835935409219377e-18, "g": 4.452212822832222e-18} #flam/arcsec^2
sky_background_countrate = {}

# (Teledyne e2v CMOS detector)
# Dark current is 0.01 electrons/s/pixel at -50°C and halves for every reduction of 5-6°C.
# CASTOR operates at 180 K, implying:
# dark current = 0.5^((223.15K - 180K) / 6) * 0.01 ~ 1e-4 electrons/s/pixel (negligible)
dark_countrate = 1e-4  # electrons/s/pixel

bias_count = 100 #electrons

read_noise_count = 2 #e-/pix

add_bkgrd_noise_count = 0 #add noise per reads e-/pix

# -------------------------------------------------------------------------------------- #


def tot_unif_noise(self, exposure_time, filter="uv", Nreads=1, Nbin=1):

    #Nreads: total number of read-outs
    #Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used.
    #add_bkgrd_noise: additionnal noise, eg, HST: the background added using the post-flash option in e− pixel-1

    self.total_unif_noise = self.sky_background_countrate[filter] * exposure_time + self.dark_countrate * exposure_time + self.read_noise_count**2 * Nreads / Nbin + self.add_bkgrd_noise_count * Nreads


def recompute_sky_background(self, default_sky=True, user_zodi=None, user_earthshine=None, user_geo=None):
    '''
    Function to create the grism sky background.
    Assuming that the sky background is dispersed over all the pixels in the grism imaging,
    the grism sky_background is simply uniform and the same as the imaging sky background.

    if default_sky is False, user_zodi, user_earthshine, and user_geo are paths to files containing
    the corresponding background in the apropriate format

    returns
        grism_sky: 2D array of the grism sky backfground (zodi+earthshine+geocoronal)
                    same size as self.grism_box
        
    '''

    if default_sky == False:
        path_to_zodi = user_zodi
        path_to_earthshine = user_earthshine
    else:
        path_to_zodi = utils.ETC_GRISM_HOME+'/files/sky_background/zodi.fits'
        path_to_earthshine = utils.ETC_GRISM_HOME+'/files/sky_background/earthshine.fits'

    spectrum_zodi = fits.open(path_to_zodi)
    spectrum_earthshine = fits.open(path_to_earthshine)

    wave_zodi = spectrum_zodi[1].data['WAVELENGTH']
    flux_zodi = spectrum_zodi[1].data['FLUX']
    wave_es = spectrum_earthshine[1].data['WAVELENGTH']
    flux_es = spectrum_earthshine[1].data['FLUX']

    #get total
    if np.prod(wave_zodi == wave_es) == 1:
        wave = wave_zodi
        flux = flux_zodi+flux_es
    else:
        raise("Zodiacl and Earthshine not on same wavelength grid.")

    self.sky_background_flam = {}
    self.sky_background_mag = {}
    self.sky_background_countrate = {}


    for filter in ["uv", "u", "g"]:
        filter_transmission_file = ascii.read(utils.ETC_GRISM_HOME+'files/passbands/passband_castor.'+filter)
        wave_transmission = filter_transmission_file['col1'] * 1e4 #in angtroms
        filter_transmission = filter_transmission_file['col2']

        #interpolate transmission onto finer grid
        dlambda = 1  #in angstroms
        grid_wl = np.arange(np.min(wave_transmission), np.max(wave_transmission)+dlambda, dlambda)
        grid_ftrans = np.interp(grid_wl, wave_transmission, filter_transmission)

        #resample to wavelength grid of filter transmission curve
        flux_resamp = spectres.spectres(grid_wl, wave, flux, spec_errs=None, verbose=False)

        #flux as observed through filter
        flux_numerator_filter = np.sum(flux_resamp * grid_ftrans * grid_wl * dlambda)
        flux_norm_filter = np.sum(grid_ftrans * grid_wl * dlambda)
        flux_filter = flux_numerator_filter/flux_norm_filter

        #convert to magnitudes
        #compute effective wavelengths to convert f_lambda to mag
        efflamb_numerator = np.sum(grid_ftrans * dlambda)
        efflamb_norm = np.sum(grid_ftrans * dlambda / grid_wl**2)
        efflamb = np.sqrt(efflamb_numerator / efflamb_norm)
        mag_filter = utils.flam_to_mag(flux_filter, efflamb)

        #populate sky_background dictionnary
        self.sky_background_flam[filter] = flux_filter  #flam per arcsec**2
        self.sky_background_mag[filter] = mag_filter    #mag per arcsec**2

        #compute background in count rate (electrons per seconds)
        pix_scale = 0.1 #arcsec per pixel

        #photflam is the inverse sensitivity (units: erg cm−2 Å−1 electron−1), as a function of wavelength.
        #This is the scaling factor to transform an instrumental flux in units of electrons per second to a physical flux density, as a function of wavelength.
        #FOR JWST NIRISS, PHOTFLAM = 1.9429001E-20 (ground calibration?) for pivot wavelength 15369.176 A.
        photflam = 2e-19 #arbitrary, Needs to be updated with real value, depending on filter.
        self.sky_background_countrate[filter] = flux_filter * photflam * pix_scale**2

    return 0


def dark(self, dark_current):
    #e/s/pix
    self.dark_countrate = dark_current

def bias(self, bias):
    # electron
    self.bias_count = bias

def read_noise(self, read_noise):
    # electron/pixel
    self.read_noise_count = read_noise

def add_bkgrd_noise(self, add_noise):
    #additionnal noise per reads in e/pix. eg, HST: the background added using the post-flash option in e− pixel-1
    self.add_bkgrd_noise_count = add_noise

GAIN = 2.0  # electron/ADU
