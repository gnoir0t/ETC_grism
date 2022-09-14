"""
Grism ETC
"""

import time
import os
import sys
import pkg_resources
from tqdm import tqdm

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from astropy.io import ascii, fits
from astropy.table import Table, vstack, hstack
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

import spectres
import h5py

import fsps
sp = fsps.StellarPopulation(compute_vega_mags=False, zcontinuous=1, sfh=4, imf_type=2, dust_type=2, add_neb_emission=True, add_igm_absorption=True)
#fsps IGM absorption is Madau 1995.
