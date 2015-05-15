import numpy
import orbit
import ma02
import scipy.optimize
import pylab

def lightcurve(JD, P, p, Ttr = 0, Ecc = 0, a = 10, \
               incl = numpy.pi/2, \
               omega = 0, limbd = [0, 0], sec = 0):
    """Calculate transit / eclipse light curve according to Mandel &
    Agol (2002), with linear, quadratic or 4-parameter non-linear limb-darkening.
    TSB: modified the code to calculate the model on a finer grid than the
    original data is sampled on. The model is then averaged over the finer data.
    This accounts for non-linear changes in the brightnesses during ingress
    and egress. """
    #try binning model
    JD = numpy.array(JD)
    avediff = numpy.median(JD[1:] - JD[0:-1])
    #add in 4 extra data points per obs
    JDnew = numpy.zeros(len(JD)*5)
    flagorig = numpy.zeros(len(JD)*5)
    for i in numpy.arange(len(JD)):
        #i*5 + 2
        JDnew[(i*5)] = JD[i] - (2.*avediff*0.25)
        JDnew[(i*5)+1] = JD[i] - (avediff*0.25)
        JDnew[(i*5)+2] = JD[i] 
        JDnew[(i*5)+3] = JD[i] + (avediff*0.25)
        JDnew[(i*5)+4] = JD[i] + (2.*avediff*0.25)

        flagorig[(i*5)+2] = 1.
    JDorig = numpy.copy(JD)
    JD = JDnew

    T0 = orbit.getT0(P, Ttr, omega, Ecc)
    x,y,z = orbit.skycoord(JD, P, T0, Ecc, a, incl, 0, omega)
    r = numpy.sqrt(x**2+y**2)
    res = numpy.ones(len(JD))
    f_uni = 1.0 - ma02.occultuniform(r, p)
    if numpy.size(limbd) == 0:
        f_tr = f_uni
    elif numpy.size(limbd) == 4:
        f_tr = ma02.occultnonlin(r, p, limbd)
    else:
        f_tr = ma02.occultquad(r, p, limbd)
    f_sec = 1.0 - (1.0 - f_uni) * sec / p**2

    #now sum up model
    f_trnew = numpy.zeros(len(JDorig))
    f_secnew = numpy.zeros(len(JDorig))
    znew = numpy.zeros(len(JDorig))
    for i in numpy.arange(len(JDorig)):
        f_trnew[i] = (f_tr[(i*5)] + f_tr[(i*5)+1] + f_tr[(i*5)+2] + 
            f_tr[(i*5)+3] + f_tr[(i*5)+4]) / 5.
        f_secnew[i] = (f_sec[(i*5)] + f_sec[(i*5)+1] + f_sec[(i*5)+2] + 
            f_sec[(i*5)+3] + f_sec[(i*5)+4]) / 5.
        znew[i] = z[(i*5)+2]
    f_tr = f_trnew
    f_sec = f_secnew
    z = znew
    return numpy.select([z > 0, z <= 0], [f_tr, f_sec])

def transit_errfunc_ptia(par_tofit, par_fixed, JD, flux, flux_err):
    p, Ttr, incl, a = par_tofit
    P, Ecc, omega, limbd, sec = par_fixed
    flux_model = lightcurve(JD, P, p, Ttr, Ecc, a, incl, omega, limbd, sec)
    chi2 = numpy.sum((flux - flux_model)**2 / flux_err**2)
    print 'chi2: ', chi2
    return chi2

def transit_fit_example():
    '''Simulate some data containing a transit plus white Gaussian
    noise and fit for the planet / star radius ratio (p = Rp / Rs),
    epoch Ttr, inclination, and ratio of semi-major axis to stellar
    radius (a/Rs), while holding the other parameters fixed at their
    true values.'''
    JD = numpy.arange(100) * 0.01 - 0.5 + pylab.normal(0, 0.005, 100)
    P = 3.0
    Ttr_true = 0.0; Ttr_guess = 0.02
    p_true = 0.1; p_guess = 0.08
    incl_true = 0.99 * numpy.pi / 2.; incl_guess = numpy.pi / 2.
    a_true = 5.0; a_guess = 4.5
    cn = [0, 0.45, 0, 0.05]
    Ecc = 0; omega = 0; sec = 0
    truth = lightcurve(JD, P, p_true, Ttr_true, Ecc, a_true, incl_true, \
                       omega, cn, sec) 
    flux = truth + pylab.normal(0, 0.002, numpy.size(JD))
    flux_err = numpy.ones(numpy.size(JD)) * 0.002
    initial_guess = lightcurve(JD, P, p_guess, Ttr_guess, \
                               Ecc, a_guess, incl_guess, omega, cn, sec)
    p_fit, Ttr_fit, incl_fit, a_fit = \
        scipy.optimize.fmin(transit_errfunc_ptia, \
                            [p_guess, Ttr_guess, incl_guess, a_guess], \
            args = ([P, Ecc, omega, cn, sec], JD, flux, flux_err))
    fit = lightcurve(JD, P, p_fit, Ttr_fit, \
                               Ecc, a_fit, incl_fit, omega, cn, sec)
    pylab.clf()
    pylab.errorbar(JD, flux, yerr = flux_err, fmt = 'k.', label = 'data')
    pylab.plot(JD, truth, 'k-', label = 'true')
    pylab.plot(JD, initial_guess, 'b--', label = 'initial')
    pylab.plot(JD, fit, 'r-', label = 'fit')
    pylab.legend(loc = 'lower right')
    print p_true, p_guess, p_fit
    print Ttr_true, Ttr_guess, Ttr_fit
    print incl_true, incl_guess, incl_fit
    print a_true, a_guess, a_fit
    return
