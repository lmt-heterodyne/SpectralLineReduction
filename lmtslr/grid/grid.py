"""
module for handling SEQUOIA grid geometry

classes: Grid
methods: 
uses: 
author: FPS 
date: September 2019 (documented)
changes:
python 3
PJT:     handle AE/RD/LL coordinates
"""

import numpy as np

class Grid():
    """
    Class to define the geometry of the SEQUOIA array and provide 
    methods to compute offsets.
    """
    def __init__(self, rx="Sequoia", spacing=27.9, rotation=46.0):
        """
        Constructor for Grid class.
        Args:
            spacing (float): spacing of beams in grid in arcsec
            rotation (float): rotation of array relative to telescope 
                in degrees 
        """
        if rx == "Msip1mm":
            self.nhorns = 6
            self.spacing = 0
            self.rotation = 0
            # these are the offsets in the grid ordered by pixel
            self.RIGHT = np.array([0, 0, 0, 0, 0, 0])
            self.UP = np.array([0, 0, 0, 0, 0, 0])
        else:
            self.nhorns = 32
            self.spacing = spacing
            self.rotation = rotation/180.*np.pi
            # these are the offsets in the grid ordered by pixel
            self.RIGHT = np.array(2*[-1.5, -1.5, -1.5, -1.5, -.5, -.5, -.5, -.5,
                                    .5, .5, .5, .5, 1.5, 1.5, 1.5, 1.5])
            self.UP = np.array(2*[1.5, .5, -.5, -1.5, 1.5, .5, -.5, -1.5, 1.5, 
                                .5, -.5, -1.5, 1.5, .5, -.5, -1.5])

    def azel(self, elev, tracking_beam):
        """
        Returns the offsets of beams in az-el system with respect to 
        position being tracked.
        Args:
            elev (float): elevation in radians 
            tracking_beam (int): pixel number being tracked (-1 for 
                center of array)
        Returns:
            azmap (array): array with azimuth offset positions for each
                beam
            elmap (array): array with elevation offset positions for 
                each beam
        """
        azmap = self.spacing * (self.RIGHT * np.cos(self.rotation - elev) 
                                - self.UP * np.sin(self.rotation - elev))
        elmap = self.spacing * (self.RIGHT * np.sin(self.rotation - elev) 
                                + self.UP * np.cos(self.rotation - elev))
        if tracking_beam >= 0: # tracking a pixel
            track_az = azmap[tracking_beam]
            track_el = elmap[tracking_beam]
            azmap = azmap - track_az
            elmap = elmap - track_el
            
        return(azmap, elmap)

    def radec(self, elev, parang, tracking_beam):
        """
        Returns the offsets of beams in ra-dec system with respect to 
        position being tracked.
        Args:
            elev (float): elevation in radians 
            parang (float): paralactic angle in radians
            tracking_beam (int): tracking_beam is the pixel number 
                being tracked (set -1 for center of array)
        Returns:
            ramap (array): array with right ascension offset positions 
                           for each beam
            decmap (array): array with declination offset positions 
                            for each beam
        """
        azmap, elmap = self.azel(elev, tracking_beam)
        ramap = - azmap * np.cos(parang) + elmap * np.sin(parang)
        decmap = + azmap * np.sin(parang) + elmap * np.cos(parang)

        return(ramap, decmap)

    def latlon(self, elev, parang, galang, tracking_beam):
        """
        Handle galactic Lat Long
        galang (float): angle between B and Dec measures clockwise. 
                        B aligned with +RA = -90, B aligned with -RA = +90.
        Returns:
            lonmap (array): array with gal.longitude offset positions
                            for each beam
            latmap (array): array with gal.latitude offset positions
                            for each beam
        """
        azmap, elmap = self.azel(elev, tracking_beam)
        totang = parang + galang    # PJT @todo  + or -
        lonmap = - azmap * np.cos(totang) + elmap * np.sin(totang)
        latmap = + azmap * np.sin(totang) + elmap * np.cos(totang)

        return(lonmap,latmap)


