""" Measure the redshift from the spectrum of a galaxy

use the entire Balmer series, since you know it's there
can also use O III 
"""

import numpy as np
from scipy.optimize import curve_fit


def gaussian(x, amp, cen, wid):
    """ Fit a Gaussian """
    return amp * np.exp(-(x-cen)**2 / wid)


def fit_redshift(specfile, l0, z0):
    """
    Given a galaxy spectrum, get the best-fit redshift

    Parameters
    ----------
    specfile: the ascii file with the spectrum
    l0: the rest wavelength of the line in question
    z0: the initial guess for the redshift

    Returns
    -------
    z: best-fit redshift
    zerr: uncertainty on best-fit redshift
    """
    window = 50 # Angstroms

    dat = np.loadtxt(specfile)
    wl_all = dat[:,0]
    flux_all = dat[:,1]

    # Extract the region of interest
    lm = l0 * (z0 + 1)

    choose = np.logical_and(
            wl_all > (lm - window),
            wl_all < (lm + window))
    wl = wl_all[choose]
    flux = flux_all[choose]

    # Fit the Gaussian to this region
    best_vals, covar = curve_fit(gaissna, wl, flux, p0=init_vals)


if __name__=="__main__":
    specfile = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/ZTF18abukavn_20181109_Keck1_v1"
