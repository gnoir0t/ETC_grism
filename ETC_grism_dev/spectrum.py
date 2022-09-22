"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *
from . import utils

spectrum = (0, 0)
fsps_sp = 0

def spectrum_from_file(path_to_file, wave_key, flux_key, check=True):
    '''
    Grab spectrum from file.
    ASSUMES SPECTRUM IS ALREADY REDSHIFTED.
    Wavelength must be in angstroms.
    Flux must be in ergs/cm2/s/A (normalized or not).
    '''

    spec = ascii.read(path_to_file)
    wave = spec[wave_key] #must be in angstroms
    flux = spec[flux_key] #must be in ergs/cm2/s/A (normalized or not)

    self.spectrum = (wave, flux)

    if check:
        plt.plot(wave, flux, '-k', label='Spectrum')
        plt.xscale('log')
        plt.xlabel('Wavelength (angstroms)')
        plt.ylabel('Flux (ergs/cm2/s/A/Norm)')
        plt.legend()
        plt.show()

    return 0


def spectrum_from_fsps(self, use_fsps_params_dict=True, fsps_params_dict={"tau": 1, "logzsol": 0, "dust2": 0.2}, fsps_age = 5, redshift=0, check=True):
    '''
    Generate spectrum from fsps
    '''
    
    self.fsps_sp = sp

    #use fsps_params_dict
    for param in fsps_params_dict:
        self.fsps_sp.params[param] = fsps_params_dict[param]

    #get
    wave, flux_Lumlamb = self.fsps_sp.get_spectrum(tage=fsps_age, peraa=True) #1Mo of mass formed

    sol_lum_to_ergs_sec = 3.826*10**33
    flux = flux_Lumlamb * sol_lum_to_ergs_sec

    #Get redshifted spectrum
    wave *= (1+redshift)

    if redshift > 0:
        log_mass_form_norm = 0
        scale = 10**log_mass_form_norm

        dist_lum = cosmo.luminosity_distance(redshift).value * 3.08567758128e24 #in cm
        redshift_norm = ((4*np.pi*dist_lum**2)*(1+redshift))
        flux = flux / redshift_norm * scale

    self.spectrum = (wave, flux)

    if check:
        plt.plot(wave, flux, '-k', label='Full Spectrum')
        plt.xscale('log')
        plt.xlabel('Wavelength (angstroms)')
        plt.ylabel('Mass-Normalized Flux (ergs/cm2/s/A/Mo)')
        plt.legend()
        plt.show()

    return 0


def normalize_spectrum(self, magnitude=23, filter_channel="u", check=True):
    '''
    Normalize given spectrum to given magnitude in given CASTOR filter
    '''

    filter_transmission_file = ascii.read(utils.ETC_GRISM_HOME+'files/passbands/passband_castor.'+filter_channel)
    wave_transmission = filter_transmission_file['col1'] * 1e4 #in angtroms
    filter_transmission = filter_transmission_file['col2']

    #interpolate transmission onto finer grid
    dlambda = 1  #in angstroms
    grid_wl = np.arange(np.min(wave_transmission), np.max(wave_transmission)+dlambda, dlambda)
    grid_ftrans = np.interp(grid_wl, wave_transmission, filter_transmission)

    wave, flux = self.spectrum

    #resample to wavelength grid of filter transmission curve
    flux_resamp = spectres.spectres(grid_wl, wave, flux, spec_errs=None, verbose=False)

    #flux as observed through filter
    flux_numerator_filter = np.sum(flux_resamp * grid_ftrans * grid_wl * dlambda)
    flux_norm_filter = np.sum(grid_ftrans * grid_wl * dlambda)
    flux_filter = flux_numerator_filter/flux_norm_filter

    #compute effective wavelengths to convert mag to f_lambda
    efflamb_numerator = np.sum(grid_ftrans * dlambda)
    efflamb_norm = np.sum(grid_ftrans * dlambda / grid_wl**2)
    efflamb = np.sqrt(efflamb_numerator / efflamb_norm)

    #convert target mag to flux
    target_flux = utils.mag_to_flam(magnitude, efflamb)

    #scaling
    scale = target_flux / flux_filter
    print("TARGET_FLAM, NON_NORM_FLAM, SCALING: ", target_flux, flux_filter, scale)

    #normalize and populate self.spectrumwith new wavelength grid and resampled flux
    flux_resamp *= scale
    self.spectrum = (grid_wl, flux_resamp)

    if check:
        wave, flux = self.spectrum
        #flux as observed through filter
        flux_numerator_filter = np.sum(flux * grid_ftrans * wave * dlambda)
        flux_norm_filter = np.sum(grid_ftrans * wave * dlambda)
        flux_filter = flux_numerator_filter/flux_norm_filter

        print("NORM_FLAM, NORM_FNU, NORM_MAG: ", flux_filter, utils.flam_to_fnu(flux_filter, efflamb), utils.flam_to_mag(flux_filter, efflamb))
        print("TARGET_FLAM, TARGET_FNU, TARGET_MAG: ", target_flux, utils.mag_to_fnu(magnitude), magnitude)

        plt.plot(wave, flux, 'k-', label='Normalized Spectrum')
        plt.plot(efflamb, target_flux, 'or', label='Target Photometry')
        plt.plot(grid_wl, grid_ftrans*plt.gca().get_ylim()[1], '--b', label='Filter Curve')
        plt.xlabel('Wavelength (angstroms)')
        plt.ylabel('Flux (ergs/cm2/s/A)')
        plt.legend()
        plt.show()

    return 0




