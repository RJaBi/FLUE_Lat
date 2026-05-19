from FLUE.compiled import flue_c as fc


def stoutSmearLinks(U, rho, n_sweeps):
    """
    Apply stout-link smearing to a rank-7 gauge field array.

    Parameters
    ----------
    U : numpy.ndarray
        Gauge field with shape (NT, NX, NY, NZ, 4, 3, 3).
    rho : float
        Smearing parameter.
    n_sweeps : int
        Number of stout smearing sweeps.

    Returns
    -------
    numpy.ndarray
        Smeared gauge field of the same shape as `U`.
    """
    usmeared = fc.stoutsmearlinks_c(U, rho, n_sweeps)
    return usmeared
