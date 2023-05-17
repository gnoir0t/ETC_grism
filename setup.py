"""
Grism ETC
"""
from setuptools import setup

#long_description = __doc__.strip()  # Remove leading and trailing newlines

setup(
    name="ETC_grism_dev",
    version="0.4.0",  # see semantic versioning (<https://semver.org/spec/v2.0.0.html>)
    description="CASTOR GRISM Exposure Time Calculator (ETC)",
#    long_description=long_description,
#    url="",
    author="Gaël Noirot",
    author_email="gael.noirot@smu.ca",
    # maintainer="Gaël Noirot",
    # maintainer_email="gael.noirot@smu.ca",
    packages=[
        "ETC_grism_dev",
        "ETC_grism_dev.files",
        "ETC_grism_dev.files.grism_files",
        "ETC_grism_dev.files.passbands",
        "ETC_grism_dev.files.sky_background",
        "ETC_grism_dev.files.spectra"
    ],
    package_data={
        "ETC_grism_dev.files.grism_files": ["*_profile_uv.txt", "*_dispersion_uv.txt", "*_dispersion_u.txt", "*_efficiency_uv.txt", "*_efficiency_u.txt"],
        "ETC_grism_dev.files.passbands": ["*.uv", "*.u", "*.g"],
        "ETC_grism_dev.files.sky_background": ["earthshine.fits", "zodi.fits", "geo_low.fits", "geo_avg.fits", "geo_high.fits"],
        "ETC_grism_dev.files.spectra": ["fsps_quiescent.txt", "fsps_star_forming.txt"]
    },
    install_requires=["numpy", "matplotlib", "astropy", "spectres", "h5py", "fsps", "specutils"],
    license="GPLv3",
    python_requires=">=3.9",
    platforms=["MacOS"],  # only tested on MacOS.
)
