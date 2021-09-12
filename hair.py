"""
    Copyright 2021 Claire Power 
    Released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
    
    hair
    ------
    hair.py  is part of LAIsim.
    
    hair is glue logic to get LAIsim generated curves into appropriate formate for the analyser.
"""

import numpy as np

class Hair:
    """
    Hair strand from concentration curve
    """
    strand = []
    def __init__(self, concs):
        self.strand = concs
    def segment_into_n(self, n):
        return np.array_split(self.strand, n)
    def segment_into_len(self, n):
        return np.array_split(self.strand, int(len(self.strand)/n))
