# Module to calculate relative positions and radial velocities for
# binary system. Written by Suzanne Aigrain.

import numpy

f = numpy.MachAr()
machep = f.eps

def eccan(Ecc, M, Tol = 1.0e-8, Nmax = 50):
    """Calculate eccentric anomaly using Newton-Raphson process."""
    if M < Tol: return M
    x = Ecc * numpy.sin(M) / (1 - Ecc * numpy.cos(M))
    Eo = M + x * (1-x*x/2.)
    Diff = 1
    Flag = 0
    i = 0
    while (Diff > Tol):
        En = Eo + (M  + Ecc * numpy.sin(Eo) - Eo) / (1 - Ecc * numpy.cos(Eo))
        Diff = abs((En - Eo) / Eo)
        Eo = En
        i += 1
        if i >= Nmax:
            if Flag ==1:
                print Ecc, M
                print 'Eccan did not converge'
                return M
            Flag = 1
            i = 0
            Eo = M
            Diff = 1
    return En

def phase(JD, P, T0 = 0.0):
    """Phase-fold array of dates; result in range [0:1]."""
    Phase = ((JD-T0) % P) / P
# ensure > 0
    Phase = numpy.select([Phase >= 0., numpy.ones(len(Phase))], \
                             [Phase, Phase + 1])
    return Phase

def truean(JD, P, T0 = 0, Ecc = 0):
    """Calculate true anomaly for array of dates."""
    Phase = phase(JD, P, T0) # phases
    M = 2 * numpy.pi * Phase # mean anomaly
    if Ecc <= machep:
        return M
    eccanV = numpy.vectorize(eccan)
    E = eccanV(Ecc, M) % (2 * numpy.pi)  # eccentric anomaly
    cosE = numpy.cos(E)
    cosNu = (cosE - Ecc) / (1 - Ecc * cosE)
    Nu = numpy.arccos(cosNu) # true anomaly
    Nu = numpy.select([E <= numpy.pi, numpy.ones(len(Nu))], \
                          [Nu, 2 * numpy.pi - Nu]) # E>pi cases
    return Nu

def truedist(Nu, a, Ecc = 0):
    """True distance from true anomaly, orbital distance & eccentricity."""
    return a * (1 - Ecc**2) / (1 + Ecc * numpy.cos(Nu))

def orbitcoord(JD, P, T0 = 0, Ecc = 0, a = 1):
    """Coordinates in orbital plane. X is towards observer."""
    Nu = truean(JD, P, T0, Ecc)
    r = truedist(Nu, a, Ecc)
    X = r * numpy.cos(Nu)
    Y = r * numpy.sin(Nu)
    return X, Y, Nu

def skycoord(JD, P, T0 = 0, Ecc = 0, a = 1, \
                 incl = numpy.pi/2, Omega = 0, omega = 0):
    """Coordinates in plane of sky. y is North."""
    X, Y, Nu = orbitcoord(JD, P, T0, Ecc, a)
    cosi = numpy.cos(incl)
    sini = numpy.sin(incl)
    cosO = numpy.cos(Omega)
    sinO = numpy.sin(Omega)
    coso = numpy.cos(omega)
    sino = numpy.sin(omega)
    cosxX = - cosi * sinO * sino + cosO * coso
    cosxY = - cosi * sinO * coso - cosO * sino
    cosyX = cosi * cosO * sino + sinO * coso
    cosyY = cosi * cosO * coso - sinO * sino
    x = X * cosxX + Y * cosxY
    y = X * cosyX + Y * cosyY
    z = numpy.sqrt(X**2+Y**2) * sini * numpy.sin(omega+Nu)
    return x, y, z

def radvel(JD, P, K, T0 = 0, V0 = 0, Ecc = 0, omega = 0):
    """Radial velocity (user-defined semi-amplitude)"""
    Nu = truean(JD, P, T0, Ecc)
    Vr = V0 + K * (numpy.cos(omega + Nu) + Ecc * numpy.cos(omega))
    if (K < 0): Vr[:] = -999
    return Vr

def getT0(P, Ttr, omega = 0, Ecc = 0):
    """Compute time of periastron passage from time of transit centre
    and other orbital elements."""
    nu = numpy.pi/2 - omega
    cosnu = numpy.cos(nu)
    cosE = (Ecc + cosnu) / (1 + Ecc * cosnu)
    E = numpy.arccos(cosE)
    if (nu > -numpy.pi) and (nu < 0.0): E = -E  #mcquillan:1/4/10
    M = E - Ecc * numpy.sin(E)
    T0 = Ttr - M * P / (2 * numpy.pi)
    return T0
