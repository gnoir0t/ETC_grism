Note that although the FITS headers for `earthshine.fits` and `zodi.fits` say the spectrum
is in units of erg/s/cm^2/A, the actual units are in erg/s/cm^2/A/arcsec^2. It does not
make sense if the sky background values are not in units of surface brightness. What has
likely happened is that these sky background values are given for a specific aperture
(i.e., an aperture of 1 sq. arcsec), hence the header units are in erg/s/cm^2/A.

P.S. The data are from:
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/Zodi.fits>
and
<https://github.com/spacetelescope/pysynphot/blob/master/pysynphot/data/generic/earthshine.fits>.

The geocoronal spectra ('geo_high.fits', 'geo_avg.fits', and 'geo_low.fits') have been generated using the 'generate_geo_spectra()' function in 'background.py'. It uses the line intensities as defined in <https://hst-docs.stsci.edu/stisihb/chapter-6-exposure-time-calculations/6-6-tabular-sky-backgrounds> and assume they are gaussians. Note however that they are not resolved at the sampling of the generated spectra (dlambda=1 angstrom). Units are erg/cm**2/s/A/arcsec**2 (equivalent to erg/cm**2/s/arcsec**2) given the dlambda=1 angstrom sampling.