''' module for handling SEQUOIA grid geometry

classes: Grid
methods: 
uses: 
author: FPS                                                                                            date: September 2019 (documented)
changes:
python 3
'''

import numpy as np

class Grid():
    ''' grid object defines the geometry of the SEQUOIA array and provides methods to compute offsets
    '''
    def __init__(self, spacing=27.9, rotation=46.0):
        ''' constructor for grid object 
        input:
        spacing - spacing of beams in grid in arcsec
        rotation - rotation of array wrt telescope in degrees 
        '''
        self.nhorns = 16
        self.spacing = spacing
        self.rotation = rotation/180.*np.pi
        # these are the offsets in the grid ordered by pixel
        self.RIGHT = np.array([-1.5, -1.5, -1.5, -1.5, -.5, -.5, -.5, -.5, .5, .5, .5, .5, 1.5, 1.5, 1.5, 1.5])
        self.UP = np.array([1.5, .5, -.5, -1.5, 1.5, .5, -.5, -1.5, 1.5, .5, -.5, -1.5, 1.5, .5, -.5, -1.5])

    def azel(self, elev, tracking_beam):
        ''' computes offsets of beams in az-el system with respect to position being tracked
        input:
        elev - elevation in radians 
        tracking_beam - tracking_beam is the pixel number being tracked (set -1 for center of array)
        return:
        azmap - array with azimuth offset positions for each beam
        elmap - array with elevation offset positions for each beam
        '''
        azmap = self.spacing*(self.RIGHT*np.cos(self.rotation-elev)-self.UP*np.sin(self.rotation-elev))
        elmap = self.spacing*(self.RIGHT*np.sin(self.rotation-elev)+self.UP*np.cos(self.rotation-elev))
        if tracking_beam >= 0: # tracking a pixel
            track_az = azmap[tracking_beam]
            track_el = elmap[tracking_beam]
            azmap = azmap - track_az
            elmap = elmap - track_el
            
        return(azmap,elmap)

    def radec(self, elev, parang, tracking_beam):
        ''' computes offsets of beams in ra-dec system with respect to position being tracked
        input:
        elev - elevation in radians 
        parang - paralactic angle in radians
        tracking_beam - tracking_beam is the pixel number being tracked (set -1 for center of array)
        return:
        ramap - array with right ascension offset positions for each beam
        decmap - array with declination offset positions for each beam
        '''
        azmap, elmap = self.azel(elev,tracking_beam)
        ramap = -azmap*np.cos(parang)+elmap*np.sin(parang)
        decmap =+azmap*np.sin(parang)+elmap*np.cos(parang)

        return(ramap,decmap)

