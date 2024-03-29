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

    #self.grism_box is in e-/s
    #self.integrated_grism_box_count is then in units of e-
    self.integrated_grism_box_count = self.grism_box * exposure_time

    #NEW METHOD GIVES 'self.grism_box * exposure_time' IN e- (INSTEAD OF ergs/cm2/A BEFORE).
    #NO LONGER NEEDED TO USE THE INVERSE SENSITIVITY FUNCTION.
    #integrated_grism_box_count is in electrons

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
    #get full grism throughtput
    grism_throughtput = ascii.read(utils.ETC_GRISM_HOME+'files/grism_files/castor_grism_efficiency_.'+grism_channel+'.txt')
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
    margin_grism_profile = 0./100. #in percent
    margin_grism_dispersion = 0./100. #in percent
    #spatial profile extent
    nb_rows_grism = nb_rows_img + int(np.nanmax(grism_psf_profile['col2']) + np.abs(np.nanmin(grism_psf_profile['col2'])))
    nb_rows_grism = nb_rows_grism + int(nb_rows_grism*margin_grism_profile)*2
    #dispersion direction extent
    nb_columns_grism = int(np.nanmax(grism_dispersion['col1']))
    nb_columns_grism = nb_columns_grism + nb_columns_img + int((nb_columns_grism+nb_columns_img)*margin_grism_dispersion)*2
    #create box
    self.grism_box = np.zeros((nb_rows_grism, nb_columns_grism))

    #grism offsets wrt direct imaging
    y_offset = int((nb_rows_grism - nb_rows_img) / 2) #in pixels
    #arbitrary offset in dispersion direction
    x_offset = 0  #in pixels

        #NOW USING FULL GRISM THROUGHTPUT
        #    #Multiply spectrum with filter response curve and grism efficiency curve
        #    wave, flux = source_spectrum
        #    grid_ftrans = np.interp(wave, filter_transmission['col1']*1e4, filter_transmission['col2'])
        #    flux *= grid_ftrans
        #
        #    #Use dummy grism_efficiency for now (flat 25% transmission)
        #    grism_efficiency = 0.25
        #    flux *= grism_efficiency
    wave, flux = np.copy(source_spectrum)
    #Converts ergs/cm2/s/A to photons/cm2/s/A
    flam_to_photonlam = 5.0341165675427094e7
    spectrum_photon_density = flam_to_photonlam * wave * flux #in photons/cm2/s/A
    #Multiply spectrum (in photons/cm2/s/A) with full grism throughtput to get e-/cm2/s/A
    wavelength_key = grism_throughtput.keys()[0]
    thput_key = grism_throughtput.keys()[1]
    grid_gthput = np.interp(wave, grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key])
    spectrum_electron_density = spectrum_photon_density * grid_gthput   #"On detector" spectrum in e-/cm2/s/A
    flux = spectrum_electron_density

    #Resample spectrum to pixel dispersion
    #Pixel 1 corresponds to 3000 Angstrom for u channel and 1500 Angstrom for uv channel
    if grism_channel=="u":
        wave_zp = 3000
    if grism_channel=="uv":
        wave_zp = 1500
    self.wavelength_array = np.array([wave_zp+np.sum(grism_dispersion['col2'][:i]*10) for i in range(len(grism_dispersion['col2']))])
    #If spectrum doesn't fully overlap with wavelength_array, values of 0 are assumed.
    flux_resamp = spectres.spectres(wavelength_array, wave, flux, fill=0)

    if check:
        #show 1D spectrum as seen "on detector"
        plt.figure()
        plt.plot(source_spectrum[0],source_spectrum[1],'-k', label='Emitted Spectrum')
        plt.xlim(min(wavelength_array)*0.9, max(wavelength_array)*1.1)
        norm_y = 1.1*max(source_spectrum[1][(source_spectrum[0]>plt.gca().get_xlim()[0]) & (source_spectrum[0]<plt.gca().get_xlim()[1])])
        plt.ylim(0, norm_y)
        plt.plot(grism_throughtput[wavelength_key]*10, 10**grism_throughtput[thput_key] * norm_y, '-', color='grey', label='End-to-end grism throughtput', zorder=-1)
        plt.xlabel('Wavelength (angstroms)')
        plt.ylabel('Flux (ergs/cm2/s/A)')
        plt.legend()
        plt.figure()
        plt.plot(wavelength_array,flux_resamp,'-r', ds='steps-mid', label='"On Detector" Spectrum')
        plt.xlim(min(wavelength_array)*0.9, max(wavelength_array)*1.1)
        plt.xlabel('Wavelength (angstroms)')
        plt.ylabel('"On detector" Electron Flux (e-/cm2/s/A)')
        plt.show()

    #convert flux from electron flux densities (e-/cm2/s/A) to electron flux of e-/cm2/s (multiply by dispersion).
    #This is done here to be able to convert flux densities to counts (electrons) during the exposure.
    #If done at the exposure stage, then it complicates things as different wavelengths with different dispersions then overlap on the same pixels, which makes the conversion impracticle at this stage.
    flux_resamp *= grism_dispersion['col2']*10 #in e-/cm2/s

    #Make psf profile irradiance spectrum array. Currently assumes no PSF variation across the different grisms.
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
        #print("y:", y_len_spectrum_spatial, y_stop, y_len_spectrum_spatial, y_indice)
        #print("x:", x_len_spectrum_spatial, x_stop, x_len_spectrum_spatial, x_indice)
        #print("spectrum_spatial_scaled:", len(spectrum_spatial_scaled))
        #print("y_start:y_stop, x_indice:x_stop", y_stop-y_start, x_stop-x_indice)
        #return grism_box, indice, spectrum_spatial

        self.grism_box[y_start:y_stop, x_indice:x_stop] += spectrum_spatial_scaled

#        if check:
#            plt.figure()
#            plt.imshow(self.grism_box, aspect="auto", interpolation="none")
#            plt.show()

    #mirror area
    mirror_diameter = 100 #cm
    mirror_area = np.pi * 0.25 * mirror_diameter**2 #cm2

    #Multiply by mirror area to get a count rate in e-/s
    self.grism_box *= mirror_area

    if check:
        plt.figure()
        plt.imshow(self.grism_box, aspect="auto", interpolation="none", origin="lower")
        plt.colorbar(label='"On detector" Count Rate (e-/s)')
        plt.xlabel('Pixels (Dispersion direction)')
        plt.ylabel('Pixels (Spatial direction)')
        plt.show()

    return 0


#function to disperse a full scene with multiple objects
def observe_scene(self, scene_direct, scene_seg, spectra, grism_channel, exposure_time, image_fov=False, check=True):
    
    #how many segmentation regions?
    seg_regions = np.unique(scene_seg[scene_seg>0])

    #check if sources and spectra match
    if len(seg_regions) != len(spectra):
        print("ERROR: NUMBER OF SOURCES (SEG REGIONS) AND SPECTRA DO NOT MATCH."
                +" DOUBLE CHECK EACH SOURCE AS AN ASSOCIATED SPECTRUM."
                +" ALTERNATIVELY, THERE MIGHT BE OVERLAPPING SEGMENTATION REGIONS IN THE IMAGE (IF A SINGLE DIRECT IMAGE IS USED),"
                +" OR IN SOME OF THE LAYERS/DIMENSIONS (IF MULTIPLE DIMENSIONS ARE USED)."
                +" MAKE SURE NO SEG REGIONS ARE OVERLAPPING IN EACH IMAGE.")
        return "FAIL"

    #prepare grism scene
    fov_grism = scene_direct.shape[1]
    #for fov_grism: double fov_direct to account for sources as large as entire fov (worst case scenario)
    #and add 500 pixels to account for full grism dispersion.
    #can be trimmed back to original image fov after.
    fov_grism = 2 * scene_direct.shape[1] + 500
    self.grism_scene = np.zeros([fov_grism, fov_grism])
    print("FOV:", fov_grism)
    
    #loop over each segmentation region
    k=0
    for seg_region in seg_regions:
        #select region
        region = scene_seg==seg_region
        region_indexes = np.where(region)
        
        #make box around seg region
        y_min = np.min(region_indexes[0])
        y_max = np.max(region_indexes[0])
        x_min = np.min(region_indexes[1])
        x_max = np.max(region_indexes[1])
        box_seg = region[y_min:y_max+1, x_min:x_max+1]
        box_direct = scene_direct[y_min:y_max+1, x_min:x_max+1]
        
        if check==True:
            plt.figure()
            plt.imshow(np.log10(box_direct), origin='lower', interpolation='none', vmin=-1, vmax=1)
            plt.show()
            plt.figure()
            plt.imshow(box_seg, origin='lower', interpolation='none')
            plt.show()
    
        #disperse
        disperse(self, source_image=box_direct, source_disperse_region=box_seg,
                   source_spectrum=spectra[k], grism_channel=grism_channel, check=check)

        #observe
        expose(self, exposure_time=exposure_time)
        grism_2d_box = self.integrated_grism_box_count
        
        #add trace to grism scene
        x_start = int(np.median([x_min, x_max]))
        x_end = x_start + grism_2d_box.shape[1]
        y_gbox_start = int(grism_2d_box.shape[0] / 2) - int(box_seg.shape[0] / 2)
        y_gbox_end = y_gbox_start + int(box_seg.shape[0])
        self.grism_scene[y_min:y_max+1, x_start:x_end] += grism_2d_box[y_gbox_start:y_gbox_end]
        
        k+=1

    #trim fov back to original image fov?
    if image_fov==True:
        self.grism_scene = self.grism_scene[:scene_direct.shape[1], :scene_direct.shape[0]]

    return 0


#function to get associated noise map of a full scene with multiple objects
def scene_noise(self, exposure_time=1000, filter="uv", Nreads=1, Nbin=1):
    '''
    Function to generate the total noise of a grism scene for a given integration time.

    exposure_time
        Exposure time in seconds.

    Nreads: total number of read-outs (int).

    Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used (int).

    '''

    self.grism_scene_noise = np.zeros_like(self.grism_scene)

    #compute the uniform background
    background.tot_unif_noise(background, exposure_time=exposure_time, filter=filter, Nreads=Nreads, Nbin=Nbin)

    #compute total background
    self.grism_scene_noise = np.sqrt( self.grism_scene
                                        + background.total_unif_noise
                                        )

    return 0


#function to disperse a full scene with multiple objects overlapping in the direct imaging
#It takes as input multi-dimensionnal "scene_direct", "scene_seg", and "spectra", each with non-overlapping sources/segmentations.
#Other inputs same as "observe_scene()".
def observe_scene_multi_dim(self, scene_direct, scene_seg, spectra, grism_channel, exposure_time, image_fov=False, check=True):

    grism_scene_md = []

    #loop over the dimensions
    for scene_direct_i, scene_seg_i, spectra_i in zip(scene_direct, scene_seg, spectra):
        #create scene for each dimension
        observe_scene(self, scene_direct_i, scene_seg_i, spectra_i, grism_channel, exposure_time, image_fov=image_fov, check=check)
        grism_scene_md.append(self.grism_scene)
    
    self.grism_scene_multi_dim = np.sum(grism_scene_md, axis=0)
    
    return 0


#function to to get associated noise map of a full scene with multiple dimensions.
def scene_noise_md(self, exposure_time=1000, filter="uv", Nreads=1, Nbin=1):
    '''
    Function to generate the total noise of a grism scene with multiple dimensions for a given integration time.

    exposure_time
        Exposure time in seconds.

    Nreads: total number of read-outs (int).

    Nbin: the number of detector pixels binned to one read-out pixel when on-chip binning is used (int).

    '''

    self.grism_scene_noise_multi_dim = np.zeros_like(self.grism_scene_multi_dim)

    #compute the uniform background
    background.tot_unif_noise(background, exposure_time=exposure_time, filter=filter, Nreads=Nreads, Nbin=Nbin)

    #compute total background
    self.grism_scene_noise_multi_dim = np.sqrt( self.grism_scene_multi_dim
                                        + background.total_unif_noise
                                        )

    return 0

