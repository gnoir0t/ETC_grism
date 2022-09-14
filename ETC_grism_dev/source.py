"""
Grism ETC
"""
#from .__init__ import *
from .all_imports import *

direct_image = 0
disperse_region = 0

def direct_image_from_file(self):
    '''
    Grab direct image from file.
    '''
    
    self.direct_image = 5
    self.disperse_region = 5

    return 0


def disperse_region_from_file():
    '''
    Grab region to disperse from file.
    '''

    return 0


def direct_image_castor_etc():
    '''
    Grab direct image using castor etc to generate source.
    '''

    return 0


def disperse_region_castor_etc():
    '''
    Grab region to disperse using castor etc to generate source aperture.
    '''

    return 0
