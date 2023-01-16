"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *
from . import utils


# ------------------------------- BACKGROUND NOISE VALUES ------------------------------ #

# The flux of the geocoronal emission line [O II] 2471A.
# See https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds
GEOCORONAL_FLUX_HIGH = 3.0e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_FLUX_AVG = 1.5e-15  # erg/cm^2/s/arcsec^2
GEOCORONAL_FLUX_LOW = 1.5e-17  # erg/cm^2/s/arcsec^2
GEOCORONAL_WAVELENGTH = 2471  # angstroms
GEOCORONAL_LINEWIDTH = 0.023  # angstroms

# Average sky background in each passband without geocoronal (Earthshine + zodiacal light)
#SKY_BACKGROUND = {"uv": 26.08, "u": 23.74, "g": 22.60}  # AB mag/arcsec^2
#sky_background_mag = {"uv": 27.7831, "u": 24.2583, "g": 22.6063}  # AB mag/arcsec^2
#sky_background_flam = {"uv": 1.7186261475789375e-19, "u": 1.835935409219377e-18, "g": 4.452212822832222e-18} #flam/arcsec^2
sky_background_mag = {"uv": 27.7831, "u": 24.2583}  # AB mag/arcsec^2
sky_background_flam = {"uv": 1.7186261475789375e-19, "u": 1.835935409219377e-18} #flam/arcsec^2
sky_background_countrate = {"uv": 0.0005811746388609881, "u": 0.01418100467547261} #e-/s /pix
#Imaging ETC gives sky_background countrates of uv: 0.0004768260787744109, and u: 0.01325195645039117. About the same as the method implemented here.

# Average sky background in each passband with geocoronal (Earthshine + zodiacal light + average geocoronal)
# Here, the geocoronal spectrum includes the Lya, OI, OI] and [OII] lines.
sky_background_mag = {'uv': 25.141716564724156, 'u': 24.25825403868489} # AB mag/arcsec^2
sky_background_flam = {'uv': 1.9576700516350905e-18, 'u': 1.835936959124605e-18} #flam/arcsec^2
sky_background_countrate = {'uv': 0.004671128308536505, 'u': 0.014181008896939849} #e-/s /pix
# Average geocoronal is negligible in the U band, but quite significant in the UV, including the contribution from very short wavelength lines (Lya, OI, OI]).

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
    #total_unif_noise is per pixel

    self.total_unif_noise = self.sky_background_countrate[filter] * exposure_time + self.dark_countrate * exposure_time + self.read_noise_count**2 * Nreads / Nbin + self.add_bkgrd_noise_count * Nreads


def recompute_sky_background(self, default_sky=True, geo_intensity='avg', user_zodi=None, user_earthshine=None, user_geo=None):
    '''
    Function to create the grism sky background.
    Assuming that the sky background is dispersed over all the pixels in the grism imaging,
    the grism sky_background is simply uniform and the same as the imaging sky background.

    'geo' can take value of 'high', 'avg', or 'low', and corresponds to the geocoronal line intensities; see 'generate_geo_spectra()'.

    if default_sky is False, user_zodi, user_earthshine, and user_geo are paths to files containing
    the corresponding background in the apropriate format

    returns
        grism_sky: 2D array of the grism sky backfground (zodi+earthshine+geocoronal)
                    same size as self.grism_box
        
    '''

    if default_sky == False:
        path_to_zodi = user_zodi
        path_to_earthshine = user_earthshine
        path_to_geo = user_geo
    else:
        path_to_zodi = utils.ETC_GRISM_HOME+'/files/sky_background/zodi.fits'
        path_to_earthshine = utils.ETC_GRISM_HOME+'/files/sky_background/earthshine.fits'
        path_to_geo = utils.ETC_GRISM_HOME+'/files/sky_background/geo_'+geo_intensity+'.fits'

    spectrum_zodi = fits.open(path_to_zodi)
    spectrum_earthshine = fits.open(path_to_earthshine)
    spectrum_geo = fits.open(path_to_geo)

    wave_zodi = spectrum_zodi[1].data['WAVELENGTH']
    flux_zodi = spectrum_zodi[1].data['FLUX']
    wave_es = spectrum_earthshine[1].data['WAVELENGTH']
    flux_es = spectrum_earthshine[1].data['FLUX']
    wave_geo = spectrum_geo[1].data['WAVELENGTH']
    flux_geo = spectrum_geo[1].data['FLUX']

    #get total zodi + earthshine
    if np.prod(wave_zodi == wave_es) == 1:
        wave = wave_zodi
        flux = flux_zodi+flux_es
    else:
        raise("Zodiacal and Earthshine not on same wavelength grid.")

    self.sky_background_flam = {}
    self.sky_background_mag = {}
    self.sky_background_countrate = {}

    for filter in ["uv", "u"]:
        filter_transmission_file = ascii.read(utils.ETC_GRISM_HOME+'files/passbands/passband_castor.'+filter)
        wave_transmission = filter_transmission_file['col1'] * 1e4 #in angtroms
        filter_transmission = filter_transmission_file['col2']

        #interpolate transmission onto finer grid
        dlambda = 1  #in angstroms
        grid_wl = np.arange(np.min(wave_transmission), np.max(wave_transmission)+dlambda, dlambda)
        grid_ftrans = np.interp(grid_wl, wave_transmission, filter_transmission)

        #resample sky background to wavelength grid of filter transmission curve
        #If spectrum doesn't fully overlap with wavelength_array, values of 0 are assumed.
        flux_resamp = spectres.spectres(grid_wl, wave, flux, spec_errs=None, fill=0, verbose=False)

        #resample geocoronal spectrum to same wavelength grid
        geo_resamp = spectres.spectres(grid_wl, wave_geo, flux_geo, spec_errs=None, fill=0, verbose=False)
        #add geocoronal spectrum to total sky spectrum
        flux_resamp += geo_resamp

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

        #mirror area
        mirror_diameter = 100 #cm
        mirror_area = np.pi * 0.25 * mirror_diameter**2 #cm2

#        #photflam is the inverse sensitivity (units: erg cm−2 Å−1 electron−1), as a function of wavelength.
#        #This is the scaling factor to transform an instrumental flux in units of electrons per second to a physical flux density, as a function of wavelength.
#        #FOR JWST NIRISS, PHOTFLAM = 1.9429001E-20 (ground calibration?) for pivot wavelength 15369.176 A.
#        photflam = 2e-19 #arbitrary, Needs to be updated with real value, depending on filter.
#        self.sky_background_countrate[filter] = flux_filter / photflam * pix_scale**2 #e-/s

        #NEW METHOD, WITHOUT photflam
        #DIRECTLY USES THE TOTAL SKY BACKGROUND SPECTRUM AND THE FULL GRISM THROUGHTPUT TO CONVERT FROM erg/cm2/s/A to e-/s
        flam_to_photonlam = 5.0341165675427094e7 #factor to convert flux density (erg/cm2/s/A) to photon density (photons/cm2/s/A), with photon_density = flam_to_photonlam * lambda * flux_density
        sky_photon_density = flam_to_photonlam * grid_wl * flux_resamp #in photons/cm2/s/A/arcsec**2
        #get grism throughtput
        grism_throughtput = ascii.read(utils.ETC_GRISM_HOME+'files/grism_files/castor_grism_efficiency_.'+filter+'.txt')
        wavelength_key = grism_throughtput.keys()[0]
        thput_key = grism_throughtput.keys()[1]
        grid_gthput = np.interp(grid_wl, grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key])
        #get sky background electron density
        sky_electron_density = sky_photon_density * grid_gthput #"On detector" sky background in e-/cm2/s/A/arcsec**2
        #get sky spectrum count rate
        sky_countrate = sky_electron_density * mirror_area * pix_scale**2 * dlambda #e-/s
        #each pixel sees the full spectrum when the sky is fully dispersed on the grism detector
        sky_countrate_total = np.nansum(sky_countrate) #sky countrate in e-/s, seen in each pixel (should be the same as the imaging)
        self.sky_background_countrate[filter] = sky_countrate_total

    return 0


def generate_geo_spectra(self, overwrite=False):
    '''
    Function to create the geocoronal emission line spectra for the three line intensities
    as defined in https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds
    Lines are assumed gaussian, but are unresolved at the dlambda=1 angstrom (line widths << 1 A).
    '''

    #line wavelengths in A, in order: Lya, OI, OI], [OII]
    lwaves = [1216, 1302, 1356, 2471]
    #line widths in A
    lwidth = [0.04, 0.013, 0.013, 0.023]
    #line intensity
    high = [6.1e-13, 5.7e-14, 5e-15, 3e-15]
    avg = [3.05e-13, 2.85e-14, 2.50e-15, 1.50e-15]
    low= [6.1e-14, 3.8e-16, 3e-17, 1.5e-17]
#    lya = [6.1e-13, 3.05e-13, 6.1e-14]
#    oi = [5.7e-14, 2.85e-14, 3.8e-16]
#    oib = [5e-15, 2.50e-15, 3e-17]
#    oii = [3e-15, 1.50e-15, 1.5e-17]
    #wavelength grid in A. Use dlambda=1 angstrom, so no need to divide/multiply by dlambda
    ww = np.arange(1000,11001,1)

    #generate emission lines
    lines = [utils.gaussian(ww, l_wave, l_width) for l_wave, l_width in zip(lwaves, lwidth)]
    
    #normalize
    lines_high = [line / np.sum(line) * line_intensity for line, line_intensity in zip(lines, high)]
    lines_avg = [line / np.sum(line) * line_intensity for line, line_intensity in zip(lines, avg)]
    lines_low = [line / np.sum(line) * line_intensity for line, line_intensity in zip(lines, low)]

    #add lines
    lines_high_spec = np.sum(lines_high, axis=0)
    lines_avg_spec = np.sum(lines_avg, axis=0)
    lines_low_spec = np.sum(lines_low, axis=0)

    #make astropy table
    tb_high = Table([ww, lines_high_spec], names=['WAVELENGTH', 'FLUX'], units=['Angstrom', 'flam'], dtype=['float64', 'float64'])
    tb_avg = Table([ww, lines_avg_spec], names=['WAVELENGTH', 'FLUX'], units=['Angstrom', 'flam'], dtype=['float64', 'float64'])
    tb_low = Table([ww, lines_low_spec], names=['WAVELENGTH', 'FLUX'], units=['Angstrom', 'flam'], dtype=['float64', 'float64'])
    
    #write to file
    tb_high.write(utils.ETC_GRISM_HOME+'/files/sky_background/geo_high.fits', overwrite=overwrite)
    tb_avg.write(utils.ETC_GRISM_HOME+'/files/sky_background/geo_avg.fits', overwrite=overwrite)
    tb_low.write(utils.ETC_GRISM_HOME+'/files/sky_background/geo_low.fits', overwrite=overwrite)

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
