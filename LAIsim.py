"""
    Copyright 2021 Claire Power 
    Released under CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
    
    LAIsim
    ------
    
    Long acting injectable simulator is used to generate plasma concentration curves
    using known plasma levels or pharmacokinetic paramters.
    
    It generates single dose curves and stacks them at appropriate intervals to 
    simulate multiple doses.
"""

import numpy as np
from scipy.optimize import curve_fit

class LAIsim:
    """
    Flip-flop pharmocokinetic simulator for estiamting plasma levels following single and repeated depot administration.
    """
    name = ''
    curve = None
    plasma = None
    popt = None
    pcov = None
    def __init__(self, name=''):
        self.name = name
    def tune(self, tdata, cdata, duration, p0=[150, 0.05, 0.01]):
        """
        tune tunes the flip-flop model using time and plasma level data.
        
        Parameters
        ----------
        tdata : List of float
            `tdata` time data points in plasma levels.
        cdata : float
            `cdata` concentration data points at times in `tdata`.
        duration : float
            `duration` is the number of days to similate.
        p0 : List of float
            `p0` is passed to Scipy.Optimize.curve_fit()
        """
        self.popt, self.pcov = curve_fit(self.conc_tune, tdata, cdata, p0)
        self.curve = np.array(
            list(
                map(
                    lambda t: self.conc(t, *self.popt),
                    np.arange(0, duration, 1)
                )
            )
        )
    def pk_tune(self, tmax, cmax, thalf, duration=None, p0=[150, 0.05, 0.01]):
        """
        pk_tune tunes the flip-flop model using pharmocokinetic paramaters of time to peak plasma level, peak plasma level and the half life.
        
        Parameters
        ----------
        tmax : float
            `tmax` is the time to the peak plasma concentration in days after a LAI is given.
        cmax : float
            `cmax` is the peak plasma concentration after a LAI is given.
        thalf : float
            `thalf` is the elimination half-life in days.
        duration : float
            `duration` is the number of days to similate. Default is 10 `thalf` after `tmax`.
        p0 : List of float
            `p0` is passed to Scipy.Optimize.curve_fit()
        """
        tdata = [0, tmax, tmax + thalf]#, tmax + 20*thalf]
        cdata = [0, cmax, cmax/2]#, cmax/1024]
        if duration == None:
            duration = tmax + 10 * thalf
        self.tune(tdata, cdata, duration, p0)
    def conc(self, t, d=1, m=0.7, n=0.5):
        """
        Calculates flip-flop kinematics given parameters.
        
        Parameters
        ----------
        t : int
            `t` is time in days.
        d : float
        m : float
        n : float
        
        Returns
        -------
        float
            Concentration at time t.
        """
        return (np.exp(-m*t) - np.exp(-n*t)) * (m*d)/(n-m)
    def conc_tune(self, t, d, m, n):
        """
        Internal version of conc used for scipy.Optimize.curve_fit().
        """
        #def integ(t): return (np.exp(-m*t) - np.exp(-n*t)) * (m*d)/(n-m)
        return [((np.exp(-m*t_i) - np.exp(-n*t_i)) * (m*d)/(n-m)) for t_i in t]
    def simulate_n(self, num_doses, dose_interval):
        """
        Parameters
        ----------
        num_doses : float
            `num_doses` is the number of doses given
       dose_interval : float
           `dose_interval` is the number of days between doses
        """
        nD = 0
        _curve = plasma_level = self.curve
        while nD < num_doses:
            _curve = np.pad(_curve, (dose_interval, 0), 'constant')
            plasma_level = np.pad(plasma_level, (0, dose_interval), 'constant')
            plasma_level = np.add(_curve, plasma_level)
            nD += 1
        self.plasma = plasma_level
        return plasma_level
    def save(self, path=''):
        """
        Parameters
        ----------
        path : str
            `path` is the prefix used when saving data to `path`_{single,repeated}.csv files.
        """
        if path == '':
            if self.name == '':
                raise NameError('If no Name set, then path must be provided')
            path = self.name
        if self.curve.all() == None:
            return
        days = np.arange(0, len(self.curve), 1)
        np.savez(f"{path}_single", days=days, plasma_level=self.curve)
        if self.plasma.all() == None:
            return
        days = np.arange(0, len(self.plasma), 1)
        np.savez(f"{path}_repeated", days=days, plasma_level=self.plasma)
    def parameters(self):
        d,m,n = self.popt
        return d,m,n
