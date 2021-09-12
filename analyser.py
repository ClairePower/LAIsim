"""
    Copyright 2021 Claire Power 
    Released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
    
    analyser
    ------
    analyser.py is part of LAIsim.
    
    analyser contains functions related to basic simulation of analysing samples from the hair module.
"""

import numpy as np
from scipy.signal import find_peaks

class Analyser:
    segments = []
    def __init__(self, segments):
        self.segments = segments
    def find_troughs(self):
        return find_peaks(-self.segments)
    @staticmethod
    def homogenize(segments):
        """
        homogenize takes a list of hair segments, takes the mean of each segment, normalizes to a max of 1 and returns a list
        
        Parameters
        ----------
        segements : List of float
            A segment of hair with concentration varying along it
            
        Returns
        -------
        segements : List of float
            Normalized mean values of each segment of hair
        """
        segs = np.array(list(map(np.mean, segments)))
        max_response = max(segs)
        return segs/max_response
    
