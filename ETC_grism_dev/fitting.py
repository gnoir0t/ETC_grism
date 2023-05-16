"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *
from . import utils
from . import spectrum
from . import disperser

from scipy import signal
import specutils
from specutils.spectra.spectrum1d import Spectrum1D
from specutils.analysis import correlation
from astropy.nddata import StdDevUncertainty
from astropy import constants as const


#Function to cross-correlate a source spectrum and a template and return a redhift estimate.
def cross_correlate(template_, spec_, apodization_window=0, check=True):

    '''
    inputs:
        template_: specutils Spectrum1D template

        spec_: specutils Spectrum1D observed spectrum, with both flux and uncertainties.

        apodization_window: apodization window for the cross-correlation (see specutils.analysis.correlation.template_correlate documentation)

        check: produces plot of the redshift fit (redshift corresponding to best cross-correlation between template and spectrum) if set to True

    returns:
        z_m, z_fit: the two redshift estimates.
    '''

    corr, lag = specutils.analysis.correlation.template_correlate(spec_, template_, apodization_window=apodization_window)

    # Redshift based on maximum
    index_peak = np.argmax(corr)
    v = lag[index_peak]
    z_m = v / const.c.to('km/s')

    # Redshift based on parabolic fit to mazimum
    n = 8 # points to the left or right of correlation maximum
    peak_lags = lag[index_peak-n : index_peak+n+1].value
    peak_vals = corr[index_peak-n : index_peak+n+1].value
    p = np.polyfit(peak_lags, peak_vals, deg=2)
    roots = np.roots(p)

    v_fit = np.mean(roots) * u.km/u.s # maximum lies at mid point between roots
    z_fit = v_fit / const.c.to('km/s')


    if check==True:
        plt.figure()

        plt.plot(
            template_.spectral_axis * (1. + z_m),
            template_.flux / np.max(template_.flux),
            color='red', label='Template', lw=5)

        plt.plot(spec_.spectral_axis, spec_.flux / np.max(spec_.flux), label='Observed', color='black', zorder=1)
        plt.legend(loc=2)
        plt.xlabel(spec_.spectral_axis.unit)


    return z_m, z_fit


#Function to create a specutils Spectrum1D object
def make_Spectrum1D(spectral_axis_data, spectral_axis_units, flux_data, flux_units, flux_error=None):
    '''
    inputs:
        spectral_axis_data: wavelength array of input spectrum

        spectral_axis_units: astropy-compatible units of spectral_axis_data (eg, u.AA)

        flux_data:  flux array of input spectrum

        flux_units: astropy-compatible units of flux_data (eg, u.Jy)

        flux_error: flux error array of input spectrum (if set to None, Spectrum1D object will be created without error array, eg for a template to fit to an observed spectrum)

    returns:
        spec1D: specutils Spectrum1D object.
    '''
    if flux_error is not None:
        spec1D = Spectrum1D(spectral_axis=spectral_axis_data*spectral_axis_units,
                        flux=flux_data*flux_units,
                        uncertainty=StdDevUncertainty(flux_error),
                        unit=flux_units)
    else:
        spec1D = Spectrum1D(spectral_axis=spectral_axis_data*spectral_axis_units,
                        flux=flux_data*flux_units,
                        unit=flux_units)

    return spec1D



#Function to extract 1D grism data and noise based on 2D source image
def extract_1d(source_image=None, grism_box=None, is_noise=False, extraction_size=None, check=True):

    '''
    inputs:
        2D source image

        2D grism data (signal or noise, if noise is_noise must be set to True)

        extraction_size, custom extraction size in pixels to extract 1D gris data from 2D if a source image is not provided.

        check: produces plots of the 2D and 1D extracted data if set to True

    returns:
        1D extracted grism data

        1D x-axis array.
    '''
    
    if source_image is None:
        extraction_size = extraction_size
        half_extraction_size = int(extraction_size / 2.)
    
    elif source_image is not None:
        half_source_size = int((source_image.shape[0]-1) / 2)
        half_extraction_size = half_source_size
        #Or: fit source profile
        #eg, from petrofit.segmentation import make_catalog
        #threshold = image_stddev * 2
        #kernel_size = 3
        #fwhm = 3
        #npixels = 2**2
        #cat, segm, segm_deblend = make_catalog(
        #    image.data,
        #    threshold=threshold,
        #    deblend=True,
        #    kernel_size=kernel_size,
        #    fwhm=fwhm,
        #    npixels=npixels,
        #    contrast=0.00,
        #    plot=True, vmax=vmax, vmin=vmin
        #)
        #extraction_size = cat[src]['r_eff']

    #extract
    box_center = int((grism_box.shape[0]-1) / 2)
    if is_noise == False:
        grism_1d = np.sum(grism_box[box_center-half_extraction_size:box_center+half_extraction_size+1,:],
                      axis=0)
    elif is_noise == True:
        grism_1d = np.sqrt(np.sum(grism_box[box_center-half_extraction_size:box_center+half_extraction_size+1,:]**2,
                             axis=0))
    grism_1d_x = np.arange(0,len(grism_1d), 1)

    if check == True:
        #The noiseless data
        plt.figure()
        plt.imshow(grism_box, aspect="auto", interpolation="none")
        plt.colorbar(label='Spectrum counts (e-)')
        plt.xlabel('Pixels (Dispersion direction)')
        plt.ylabel('Pixels (Spatial direction)')
        
        #1D spectrum
        plt.figure()
        plt.plot(grism_1d_x, grism_1d, '-k')
        plt.ylabel('Spectrum counts (e-)')
        plt.xlabel('Pixels (Dispersion direction)')


    return grism_1d, grism_1d_x



#Function to forward model a template to observed space ("on detector" CASTOR grism data)
#Done to properly fit template to observed spectrum in cross-correlation step.
def forward_model(spectrum_2_model, source_image, source_seg, normalize=False, source_mag_norm=23, filter_channel_norm="uv", grism_channel_disperse="uv", exposure_time=1, check=True):

    '''
    inputs:
        spectrum_2_model: spectrum.spectrum object of the template/spectrum to forward model

        source_image: source image to use for spectrum forward modelling

        source_seg: source segmentation to use for spectrum forward modelling (only pixels in seg region will be dispersed).

        normalize: whether or not to normalize spectrum to given mag in specific CASTOR band. This shouldn't matter if, for instance, forward-modelled spectrum/template is re-normalized once created.

        source_mag_norm: mag to normalize spectrum to if normalize is set to True

        filter_channel_norm: CASTOR channel to normalize spectrum to if normalize is set to True

        grism_channel_disperse: CASTOR channel to disperse (forward model) spectrum to

        exposure_time: Exposure time to use for forward modelling. Note the forward modelling produces a noiseless spectrum, so this shoudn't affect the results and will only scale the spectrum by 'exposure_time'.

    returns:
        disperser.integrated_grism_box_count: the forward-modelled (noiseless) spectrum
    '''

    #normalize?
    source_mag = source_mag_norm
    filter_channel = filter_channel_norm
    spectrum.spectrum = spectrum_2_model
    if normalize == True:
        spectrum.normalize_spectrum(spectrum, magnitude=source_mag, filter_channel=filter_channel, check=check)

    #disperse
    source_spectrum = spectrum.spectrum
    grism_channel = grism_channel_disperse
    disperser.disperse(disperser, source_image=source_image, source_disperse_region=source_seg,
                   source_spectrum=source_spectrum, grism_channel=grism_channel, check=check)

    #Create first the noiseless integrated grism spectrum.
    exposure_time = exposure_time #in seconds
    disperser.expose(disperser, exposure_time=exposure_time)

    return disperser.integrated_grism_box_count



#generate noisy spectrum based on noise data
def noisy_realisation(grism_data_1d, grism_noise_1d, check=True):

    '''
    inputs:
        grism_data_1d: noiseless 1D signal

        grism_noise_1d: corresponding 1D noise

        check: show noisy spectrum realization if set to True

    returns:
        data_noisy: the noisy spectrum realization
    '''

    noise_real = np.random.normal(0, grism_noise_1d)
    data_noisy = grism_data_1d + noise_real
    
    if check == True:
        #1D spectrum
        plt.figure()
        plt.plot(np.arange(0,len(data_noisy),1), data_noisy, '-k')
        plt.ylabel('Spectrum counts (e-)')
        plt.xlabel('Pixels (Dispersion direction)')
    
    return data_noisy


#generate X realisations of noisy spectrum based on noise data
def noisy_realisation_samples(grism_data_1d, grism_noise_1d, nsamples=500):
    '''
    inputs:
        grism_data_1d: noiseless 1D signal

        grism_noise_1d: corresponding 1D noise

        nsamples: number of realizations to generate

    returns:
        data_noisy_samples: the noisy spectrum realizations
    '''
    data_noisy_samples = []
    
    #loop over noisy_realisation()
    for n in range(nsamples):
        data_noisy = noisy_realisation(grism_data_1d, grism_noise_1d, check=False)
        data_noisy_samples.append(data_noisy)

    return data_noisy_samples


def fit_z_range(obs_spec, obs_spec_err, template_to_fit, source_image, source_seg, extraction_size, z_range=[0.2, 3.5], z_step=0.05, grism_channel_disperse=["uv"], funits='flam', check=True):
    '''
    inputs:
        obs_spec: (wavelength, flux). Wavelength in AA. Flux in ergs/cm2/s/AA, unless funits='rate' then flux in e-/sec.

        obs_spec_err: (wavelength, flux_err). Wavelength in AA (same as obs_spec). Flux err in ergs/cm2/s/AA, unless funits='rate' then flux in e-/sec.

        obs_spec_exp_time: observing time of the observed specrtum obs_spec.

        template_to_fit: (wavelength, flux). 1D rest-frame template. Wavelength in AA. Flux in ergs/cm2/s/AA, unless funits='rate' then flux in e-/sec.

        source_image: source image to use for forward modelling of template. Ideally, should be a model of the source associated to the observed spectrum obs_spec.

        source_seg: source segmentation image to use for forward modelling of template. Ideally, should be the segmentation of the source associated to the observed spectrum obs_spec.

        z_range: redshift range over which to fit template to obs_spec; e.g. z_range=[0.1, 2.2]

        z_step: redshift step of z_range.

        grism_channel_disperse: list of channels to forward model template to; e.g. ["uv", "u"].

        funits: If funits='rate' then fluxes are expected in e-/sec, if not, then ergs/cm2/s/AA are expected.

    returns:

    '''

    #create redshift array
    z_step = z_step
    z_array = np.arange(np.min(z_range), np.max(z_range)+z_step, z_step)

    #forward model template
    spectrum_2_model = template_to_fit
    source_image = source_image
    source_seg = source_seg

    ##LOOP OVER Z_ARRAY
    chi2_z_range = []
    chi2_reduced_z_range = []
    log_likelihood_z_range = []
    for z in z_array:
        template_flam = []

        #redshift the template to appropriate z
        template_wave = spectrum_2_model[0] * (1+z)
        scale = 1
        dist_lum = cosmo.luminosity_distance(z).value * 3.08567758128e24 #in cm
        redshift_norm = ((4*np.pi*dist_lum**2)*(1+z))
        template_flux = spectrum_2_model[1] / redshift_norm * scale
        spectrum_2_model_at_z = (template_wave, template_flux)

        for disp_channel in grism_channel_disperse:
            grism_box = forward_model(spectrum_2_model_at_z, source_image, source_seg, normalize=False, grism_channel_disperse=disp_channel, exposure_time=1, check=check)
            
            #extract 1D
            template_grism_1d, template_grism_1d_x = extract_1d(source_image=None, grism_box=grism_box, is_noise=False, extraction_size=extraction_size, check=check)

            #sensitivity function
            sens = Table.read(utils.ETC_GRISM_HOME+'files/grism_files/grism_sensitivity_'+disp_channel+'.fits')
            sens_wave = sens['WAVELENGTH']
            sens_curve = sens['SENSITIVITY']

            #convert template from counts (e-/s) to flam (erg/cm2/s/A) using sensitivity curve
            if funits=='rate':
                template_grism_flam_1d = template_grism_1d[:len(sens_wave)]
            else:
                template_grism_flam_1d = template_grism_1d[:len(sens_wave)] / sens_curve

            template_flam.append((sens_wave, template_grism_flam_1d))

        #stitch grisms together if need be
        if len(template_flam) > 1:
            waves = [w for w,f in template_flam]
            fluxes = [f for w,f in template_flam]
            template_wave_full = np.concatenate(waves)
            template_flam_full = np.concatenate(fluxes)
            template_flam = (template_wave_full, template_flam_full)
        else:
            template_flam = template_flam[0]
        
        template = template_flam

        #match wavelength array of obs and template spectra
        if np.prod(obs_spec[0] == template[0]) != 1:
            template = spectres.spectres(obs_spec[0], template[0], template[1], fill=np.nan)
            template = (obs_spec[0], template)

        #normalize
        scale_template = np.nansum(template[1])
        template = (obs_spec[0], template_grism_flam_1d/scale_template)
        scale_obs = np.nansum(obs_spec[1])

        #do chi2
        obs = obs_spec[1] / scale_obs
        obs_err = obs_spec_err[1] / scale_obs
        chi2, chi2_reduced = get_chi2(obs, obs_err, template[1], nb_params=1)
        chi2_z_range.append(chi2)
        chi2_reduced_z_range.append(chi2_reduced)

        #get likelihood
        log_likelihood_at_z = get_log_likelihood(chi2)
        log_likelihood_z_range.append(log_likelihood_at_z)

        #save templates
        if z==np.min(z_array):
            best_log_likelihood_at_z = log_likelihood_at_z
            best_template = template
            scaling_template = scale_template
        if log_likelihood_at_z > best_log_likelihood_at_z:
            best_template = template
            best_log_likelihood_at_z = log_likelihood_at_z
            scaling_template = scale_template

        ##
    chi2_z_range = np.array(chi2_z_range)
    chi2_reduced_z_range = np.array(chi2_reduced_z_range)
    log_likelihood_z_range = np.array(log_likelihood_z_range)

    return z_array, chi2_z_range, chi2_reduced_z_range, log_likelihood_z_range, best_template, obs, obs_err, scaling_template, scale_obs



def get_chi2(obs, obs_err, template, nb_params=1):
    chi2 = np.sum( (obs - template)**2 / obs_err**2, axis=0)

    chi2_reduced = chi2 / (len(obs) - nb_params)

    return chi2, chi2_reduced

def get_log_likelihood(chi2):

    return -chi2/2
