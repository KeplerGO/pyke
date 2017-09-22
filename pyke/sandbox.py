# This functions are not well-tested enough

def make_aperture_mask(self, snr_threshold=5, margin=4):
    """Returns an aperture photometry mask.

    Parameters
    ----------
    snr_threshold : float
        Background detection threshold.
    """

    # Find the pixels that are above the threshold in the median flux image
    median = np.nanmedian(self.flux, axis=0)
    mad = median_absolute_deviation(median[np.isfinite(median)])
    # 1.4826 turns MAD into STDEV for a Gaussian
    mad_cut = 1.4826 * mad * snr_threshold

    region = np.where(median > mad_cut, 1, 0)
    # Label all contiguous regions above the threshold
    labels = scipy.ndimage.label(region)[0]
    # Central pixel coordinate
    centralpix = [1 + median.shape[0] // 2, 1 + median.shape[1] // 2]

    # find brightest pix within margin of central pix
    central_img = median[centralpix[0] - margin: centralpix[0] + margin,
                         centralpix[1] - margin: centralpix[1] + margin]
    # unravel_index converts indices into a tuple of coordinate arrays
    brightestpix = np.unravel_index(central_img.argmax(), central_img.shape)
    bpixy, bpixx = brightestpix

    # Which label corresponds to the brightest pixel?
    regnum = labels[centralpix[0] - margin + bpixy, centralpix[1] - margin + bpixx]

    return labels == regnum

def centroids(self):
    """Returns the centroids for every cadence under a given aperture
    mask.

    Attributes
    ----------

    Returns
    -------
    xc, yc: ndarrays
        centroid positions for every cadence
    """

    xc = np.zeros(self.n_cadences)
    yc = np.zeros(self.n_cadences)
    y, x = np.mgrid[:self.shape[0], :self.shape[1]]

    x = x[self.aperture_mask]
    y = y[self.aperture_mask]

    if self._aperture_flux is None:
        self._aperture_flux = self._get_aperture_flux()

    for i in range(self.n_cadences):
        flux_i = self.flux[i][self.aperture_mask]
        xc[i] = np.nansum(flux_i * x) / self._aperture_flux[i]
        yc[i] = np.nansum(flux_i * y) / self._aperture_flux[i]

    return xc, yc
