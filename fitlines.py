""" Measure the redshift from the spectrum of a galaxy

use the entire Balmer series, since you know it's there
can also use O III 
"""

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def gaussian(x, amp, cen, wid):
    """ Fit a Gaussian """
    return amp * np.exp(-(x-cen)**2 / wid)


def fit_redshift(specfile, l0, z0, window):
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

    # Normalize
    flux_norm = flux / np.median(flux)  

    # Fit the Gaussian to this region
    init_vals = [4, lm, 3]
    best_vals, covar = curve_fit(
            gaussian, wl, flux_norm, p0=init_vals)

    # Plot the fit
    plt.plot(wl, flux_norm, c='k')
    xfit = np.linspace(min(wl), max(wl), 1000)
    yfit = gaussian(xfit, best_vals[0], best_vals[1], best_vals[2])
    plt.plot(xfit, yfit, c='r')
    plt.show()

    # Convert the best-fit into an actual redshift
    center = best_vals[1]
    ecenter = covar[1][1]
    z = (center-l0)/l0
    ez = ecenter/l0

    # Return values
    return z, ez


if __name__=="__main__":
    specfile = "/Users/annaho/Dropbox/Projects/Research/ZTF18abukavn/data/spec/ZTF18abukavn/ZTF18abukavn_20181109_Keck1_v1.ascii"
    # Initial guess
    z0 = 0.0322
    # Window size, in angstroms
    window = 20 

    # Fit for centroids of as many lines as you can tolerate
    balmer = np.array([6564.61, 4862.68, 4341.68, 4102.89, 3970.072])
    oiii = np.array([4363, 4932.6, 4960.295, 5008.24]) # O III
    # Strong lines
    lines = np.hstack((balmer[0], balmer[1], oiii[-1]))
    zall = []
    ezall = []

    # Solve
    for line in lines:
        z, ez = fit_redshift(specfile, line, z0, window)
        zall.append(z)
        ezall.append(ez)
    zall = np.array(zall)
    ezall = np.array(ezall)

    # Use the STD of the best fits as the uncertainty
    w = 1/ezall**2
    zmean = np.average(zall, weights=w)
    ezmean = np.std(zall)

    # Print the best-fit redshift, and uncertainty
    print("%s +/- %s" %(np.round(zmean,7), np.round(ezmean, 7)))
