"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *
from . import utils
from . import background


def expose(self, exposure_time=1000):
    '''
    Function to simulate a noiseless grism observation for a given integration time.

    exposure_time
        Exposure time in seconds.
        
    '''

    #assumes self.grism_box is in ergs/cm2/s/A
    #self.integrated_grism_box is then in units of ergs/cm2/A
    self.integrated_grism_box = self.grism_box * exposure_time

    #integrated_grism_box_count is in units of electron
    #where inv_sens is the inverse sensitivity (units: erg cm−2 Å−1 electron−1), as a function of wavelength.
    #This is the scaling factor to transform an instrumental flux in units of electrons per second to a physical flux density, as a function of wavelength.
    #For imaging, the inverse_sensitivity is PHOTFLAM. FOR JWST NIRISS, PHOTFLAM = 1.9429001E-20 (ground calibration?) for pivot wavelength 15369.176 A.
    inv_sens = 2e-19 #arbitrary, Needs to be updated with real array.
    self.integrated_grism_box_count = self.integrated_grism_box / inv_sens

    return 0


def total_noise(self, exposure_time=1000, filter="uv", Nreads=1, Nbin=1):
    '''
    Function to generate the total noise of a grism observation for a given integration time.

    exposure_time
        Exposure time in seconds.

    Nreads: total number of read-outs (int).

    Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used (int).

    '''

    self.grism_noise_total = np.zeros_like(self.grism_box)

    #compute the uniform background
    background.tot_unif_noise(background, exposure_time=exposure_time, filter=filter, Nreads=Nreads, Nbin=Nbin)

    #compute total background
    self.grism_noise_total = np.sqrt( self.integrated_grism_box_count
                                        + background.total_unif_noise
                                        )

    return 0


def disperse(self, source_image=None, source_disperse_region=None, source_spectrum=None, grism_channel="u", check=True):
    '''
    Grism disperser function

    source_image
        2D array of fluxes. Relative of absolute.
        Pixel scale needs to be CASTOR pixel scale (no oversampling).

    source_disperse_region
        Boolean 2D array, same size as source_image.
        Pixels with False will be masked.
        Pixels with True will be dispersed.

    source_spectrum
        Source spectrum, in flux densities.
        Source spectrum should be normalized (through spectrum.normalize_spectrum or equivalent), and have same wavelength grid as filter_transmission.

    grism_channel
        Which grism? "uv" or "u"
        
    '''

    if source_image is None:
        source_image = np.copy(source.direct_image)
    else:
        source_image = np.copy(source_image)
    if source_disperse_region is None:
        source_disperse_region = np.copy(source.disperse_region)
    else:
        source_disperse_region = np.copy(source_disperse_region)
    if source_spectrum is None:
        source_spectrum = np.copy(spectrum.spectrum)
    else:
        source_spectrum = np.copy(source_spectrum)

    #get filter transmission
    filter_transmission = ascii.read(utils.ETC_GRISM_HOME+'files/passbands/passband_castor.'+grism_channel)
    #get grism_dispersion
    #grism_dispersion = ascii.read(utils.ETC_GRISM_HOME+'files/grism_files/grism_dispersion_'+grism_channel+'_extrap.txt')
    grism_dispersion = ascii.read(utils.ETC_GRISM_HOME+'files/grism_files/grism_dispersion_'+grism_channel+'.txt')
    #get grism psf profile
    grism_psf_profile = ascii.read(utils.ETC_GRISM_HOME+'files/grism_files/grism_approx_profile_uv.txt')

    #Normalize image to get relative flux in each pixel (to scale spectum)
    norm = np.nansum(source_image[source_disperse_region])
    source_image_norm = source_image / norm

    #Get pixels of source that need to be dispered (i.e., pixels in non-masked region)
    nb_rows_img, nb_columns_img = source_image.shape
    r_i, c_i = np.where(source_disperse_region)
    indices = [(r,c) for r,c in zip(r_i, c_i)]

    #Create grism box that will be populated as we disperse the source sprectum onto it.
    margin_grism_profile = 10./100. #in percent
    margin_grism_dispersion = 5./100. #in percent
    #spatial profile extent
    nb_rows_grism = nb_rows_img + int(np.nanmax(grism_psf_profile['col2']) + np.abs(np.nanmin(grism_psf_profile['col2'])))
    nb_rows_grism = nb_rows_grism + int(nb_rows_grism*margin_grism_profile)*2
    #dispersion direction extent
    nb_columns_grism = int(np.nanmax(grism_dispersion['col1']))
    nb_columns_grism = nb_columns_grism + int(nb_columns_grism*margin_grism_dispersion)*2
    #create box
    self.grism_box = np.zeros((nb_rows_grism, nb_columns_grism))

    #grism offsets wrt direct imaging
    y_offset = int((nb_rows_grism - nb_rows_img) / 2) #in pixels
    #arbitrary offset in dispersion direction
    x_offset = 0  #in pixels

    #Multiply spectrum with filter response curve and grism efficiency curve
    wave, flux = source_spectrum
    grid_ftrans = np.interp(wave, filter_transmission['col1']*1e4, filter_transmission['col2'])
    flux *= grid_ftrans

    #Use dummy grism_efficiency for now (flat 25% transmission)
    grism_efficiency = 0.25
    flux *= grism_efficiency

    #Resample spectrum to pixel dispersion
    #Assume pixel 1 is 3000 Angstrom for u channel and 1500 Angstrom for uv channel
    if grism_channel=="u":
        wave_zp = 3000
    if grism_channel=="uv":
        wave_zp = 1500
    wavelength_array = np.array([wave_zp+np.sum(grism_dispersion['col2'][:i]*10) for i in range(len(grism_dispersion['col2']))])
    flux_resamp = spectres.spectres(wavelength_array, wave, flux, fill=np.nan)

    plt.plot(wave,flux,'-k')
    plt.plot(wavelength_array,flux_resamp+1e-19,'--r')
    plt.show()

    #Make psf profile irradiance spectrum array.
    psf_profile_pix_oversampled = grism_psf_profile['col2']
    psf_profile_pix = np.arange(np.min(psf_profile_pix_oversampled), np.max(psf_profile_pix_oversampled)+1, 1)
    psf_profile_oversampled = grism_psf_profile['col4']
    psf_profile = np.interp(psf_profile_pix, psf_profile_pix_oversampled, psf_profile_oversampled)
    #normalize to conserve total flux smeared due to psf profile
    norm = np.nansum(psf_profile)
    psf_profile_norm = psf_profile / norm
    #Smear spectrum using psf_profile in spatial direction
    spectrum_spatial = np.ones((len(psf_profile_norm),len(wavelength_array))) * flux_resamp * psf_profile_norm[:, None]

    #Loop over the pixel and disperse them (scale and add all 2D irradiance spectrum arrays).
    for indice in indices:

        #get spatial position (pixel indice) in the grism box of the source indice.
        y_indice = indice[0] + y_offset
        #get position of pixel1 in grism box
        x_indice = indice[1] + x_offset
    
        #Scale 2D array according to relative flux of each pixel in direct imaging.
        scale = source_image_norm[indice]
        spectrum_spatial_scaled = spectrum_spatial * scale

        #Add to grism box.
        y_len_spectrum_spatial = len(spectrum_spatial)
        y_start = int(y_indice - ((y_len_spectrum_spatial-1)/ 2))
        y_stop = y_start + y_len_spectrum_spatial
        x_len_spectrum_spatial = len(spectrum_spatial[0])
        x_stop = x_indice+x_len_spectrum_spatial
        #print(y_len_spectrum_spatial, y_stop, y_len_spectrum_spatial, y_indice)
        #print(x_len_spectrum_spatial, x_stop, x_len_spectrum_spatial, x_indice)
        #return grism_box, indice, spectrum_spatial

        self.grism_box[y_start:y_stop, x_indice:x_stop] += spectrum_spatial_scaled

#        plt.figure()
#        plt.imshow(self.grism_box, aspect="auto", interpolation="none")
#        plt.show()

    plt.figure()
    plt.imshow(self.grism_box, aspect="auto", interpolation="none")
    plt.show()

    return 0




#just a test to check if package works locally with the imports.
def test():

    ff = np.array([1,2,3,4,5])
    print(ff)

    spectrum.fsps_sp.params['tau'] = 0.2

    wave, spec = spectrum.fsps_sp.get_spectrum(tage=10, peraa=True)

    plt.figure()
    plt.plot(wave,spec, 'k-')
    plt.xscale('log')
    plt.show()
