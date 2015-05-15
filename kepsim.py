#!/usr/bin/env python

import pylab
from pylab import *
from copy import copy
import sys, math, copy, numpy, kepmsg

# ----------------------------------------
# binary transit/eclipse/tides/boosting phtometry model

def transitModel(anorm,M1,M2,R1,R2,period,inclination,bjd0,eccn,omega,depth,albedo,
                 c1,c2,c3,c4,gamma,contamination,npt,time,exptime,dtype,eclipses,
                 dopboost,ellipsoidal):

# takes account of:
# 1. non-linear limb darkening
# 2. doppler boosting
# 3. ellipsoidal variations
#
# INPUT ARGUMENTS
#            M1: stellar mass (Msun)
#            M2: mass of planet/companion (Msun)
#            R1: stellar radius (Rsun)
#            R2: radius of planet/companion (Rsun)
#        period: period of orbit (days)
#   inclination: orbital inclination (deg; 90 = edge-on)
#          bjd0: center of transit time (Epoch, BJD-2454900)
#         ecosw: e = eccentricity w = orbit angle of periastron
#         esinw: e = eccentricity w = orbit angle of periastron
#         depth: occultation depth (unitless)
#             c: non-linear limb-darkening (4 unitless coefficients)
#         gamma: gamma velocity (m/s)
# contamination: fractional contamination from background sources
#           npt: the number of points for the model to return
#          time: time (center) of each model point (days)
# exptime[npt]: the integration time for each model point (days)
#         dtype: 0 = photometry, 1 = RV
#
# OUTPUT VARIABLES
#        tmodel: the model for each time point
#
# OTHER PARAMETERS
#        nintg: used to find the average flux over the integration time
#            K: amplitude of RV
#         voff: radial velocity offset (m/s)
#       vmodel: model velocities (m/s)
#         vrot: rotational velocity (m/s)
#        Eanom: Eccentric anomaly
#        Manom: Mean anomaly over nintg subsampling of exposure times
# Tanom[nintg]: True anomaly
#       kepler: solves Kepler equations
#  trueanomaly: calcaulates the true anomaly
#  arad[nintg]: distance between star and planet
#         eccn: eccentricity
#          ted: parameter for thermal eclipse (mmag)
#         Psec: period of orbit (sec)
#          Per: period of orbit (days)
#        asemi: semi-major axis (m)
#         incl: orbital inclination (radians)
#   phi[nintg]: orbital phases of exposure sub-sampling (radians)
#     t[nintg]: times of exposure sub-sampling (days)
#    x2[nintg]: position of planet/companion during exposure sub-sampling (days)
#    y2[nintg]: position of planet/companion during exposure sub-sampling (days)

# PHYSICAL CONSTANTS
#           G: gravitational constant m3 kg-1 s-2
#          Cs: speed of light (m/s)
#        Msun: solar mass (kg)
#        Rsun: solar radius (m)
#         fDB: Doppler boosting factor
#          fT: ellipsoidal mass factor
#         tpi: 2 x pi
#        Pid2: pi / 2
#        c[4]: four limb darkening coefficients

# startup parameters

    nintg = 11
    vrot = 7.0e4
    G = 6.67384e-11
    Cs = 2.99792458e8
    Msun = 1.98892e30
    Rsun = 6.955e8
    fDB = 1.896
    fT = 3.37
    tPi = 2.0 * math.pi
    Pid2 = math.pi / 2.0
    
# fractional contamination from background sources

    dilute = copy.copy(contamination)

# body masses

    M1 = M1 * Msun
    M2 = M2 * Msun

# orbital period

    Per = copy.copy(period)
    Psec = copy.copy(period) * 8.64e4

# semi-major axis

    asemi = (Psec * Psec * G * (M1 + M2) / (4.0 * math.pi * math.pi))**(1.0 / 3.0)
   
# body radii

    R1 = abs(R1 * Rsun)
    R2 = abs(R2 * Rsun)

# impact parameter

    bmin = asemi / R1 * math.cos(inclination * math.pi / 180.0)

# eccentric orbit parameters
    
    ecosw = eccn * math.cos(math.pi * omega / 180)
    esinw = eccn * math.sin(math.pi * omega / 180)

# orbital eccentricity and omega

    eccn = math.sqrt(ecosw * ecosw + esinw * esinw)
    if eccn > 1.0: 
        eccn = 0.99
    if eccn == 0.0:
        eccn = 1.0e-10
        w = 1.0e-10
    else:
        w = math.atan(esinw / ecosw)
        if ecosw > 0.0 and  esinw < 0.0:
            w = tPi + w
        elif ecosw < 0.0 and esinw >= 0.0:
            w = math.pi + w
        elif ecosw < 0.0 and esinw < 0.0:
            w = math.pi + w
        if w == 0.0:
            w = 1.0e-10

# starting guess for eccentric anomaly and mean anomaly

           
    Eanom = copy.copy(w)
    Manom = copy.copy(w)

# parameter for thermal eclipse

    ted = depth * 1.0e-6

# non-linear limb darkening

    c = numpy.array([c1,c2,c3,c4],dtype='float32')

# inclination angle

    inclmin = 180.0 * math.tan((R1 + R2) / asemi) / math.pi
    inclmin=90.0-inclmin
    incl = copy.copy(inclination)
    if inclmin >= 0.0 and inclmin <= 90.0:
        if incl > 90.0:
            incl = 180.0 - incl
    incl = math.pi * (90.0 - incl) / 180.0

# observer-star-planet angle. Find phase at centre of transit

    epoch = copy.copy(bjd0)
    Eanom = kepler(Manom,Eanom,eccn)
    phi0 = trueanomaly(eccn,Eanom)

# RV initialization
 
    K = 2.0 * math.pi * G * M2**3 * (math.sin(incl + Pid2))**3 / \
        (Psec * (1.0 - eccn * eccn)**(3.0 / 2.0) * (M1 + M2) * (M1 + M2))
    K = K**(1.0 / 3.0)
    voff = copy.copy(gamma)

# normalization constant for light curve

    norm = math.pi

# initialize center of primary star to 0,0 coordinates

    x1 = 0.0
    y1 = 0.0

# subsampling of individual data points

    dnintg = float(nintg)
    dnintgm1 = 2.0 * dnintg - 2.0

# integration width initialization

    xintold = 0.0

# initialization of projected planet star distance

    y2pold = 0.0

# calculate model for each timestamp

    t = numpy.zeros((nintg),dtype='float64')
    phi = numpy.zeros((nintg),dtype='float64')
    x2 = numpy.zeros((nintg),dtype='float64')
    y2 = numpy.zeros((nintg),dtype='float64')
    Tanom = numpy.zeros((nintg),dtype='float64')
    arad = numpy.zeros((nintg),dtype='float64')
    tmodel = numpy.zeros((npt),dtype='float64')
    for i in range(npt):

# array of sub-sampled times spanning one exposure. times are centered on time[i]

        for j in range(nintg):
            t[j] = time[i] + exptime[i] * (2.0 * float(j) - dnintg - 1.0) / dnintgm1 - epoch
            phi[j]= t[j] / Per - math.floor(t[j] / Per)
            phi[j] = phi[j] * tPi
            Manom = phi[j] + w
            if Manom > tPi:
                Manom = Manom - tPi
            if Manom < 0.0:
                Manom = Manom + tPi
            Eanom = kepler(Manom,Eanom,eccn)
            Tanom[j] = trueanomaly(eccn,Eanom)
            if phi[j] > math.pi:
                phi[j] = phi[j] - tPi
            arad[j] = distance(asemi,eccn,Tanom[j])
            x2[j] = arad[j] * math.sin(Tanom[j] - phi0)
            y2[j] = arad[j] * math.cos(Tanom[j] - phi0) * math.sin(incl)

# photometric data case

        if dtype[i] == 0:

# initialize stellar surface area

            zarea = 0.0
            tflux = 0.0
            for j in range(nintg):

# doppler boosting
             
                if dopboost:
                    Kc = -K * (math.cos(Pid2 + Tanom[j] - phi0) + eccn * math.cos(w))
                    tflux = tflux + fDB * Kc / Cs

# ellipsoidal variations

                if ellipsoidal:
                    tflux = tflux + tides(M1,M2,R1,asemi,incl,Tanom[j],eccn,phi0)

# flux from star + planet/companion

                Ag = albedo * R1 * R1 / (arad[j] * arad[j])
                zarea = zarea + albedomod(t[j],Per,Ag,R1,R2,Tanom[j]-phi0)

# rescale flux and surface area 

            tflux = tflux / dnintg
            zarea = zarea / dnintg

# orbital phase (unitless) and flux change from planet transiting the star

            phase = Tanom[nintg / 2 + 1] - phi0
            if phase > math.pi:
                phase = phase - tPi
            if phase < -math.pi:
                phase = phase + tPi
            if abs(phase) < Pid2:
                (managol,b0,mulimb0,mulimbf,dist) = mandelagol(nintg,R1,R2,x1,x2,y1,y2,c)
                darea = (math.pi * managol + zarea + tflux) / norm
                vrotf = 0.0
            else:
                darea = 0.0
                for j in range(nintg):
                    darea = darea + eclmod2(R1,R2,x1,x2[j],y1,y2[j],zarea,norm,ted)
                    if j == 1:
                        xintold2 = copy.copy(xintold)
                        y2pold2 = copy.copy(y2pold)

# average area and add on delta term

                darea = darea / dnintg + tflux / norm
                xintold = copy.copy(xintold2)
                y2pold = copy.copy(y2pold2)

# Convert relative fluxes to magnitude to match observations

            tmodel[i] = darea * 1.0 + (1.0 - darea) * dilute

# RV data case

        elif dtype[i] == 1:

            tmodel[i] = 0.0
            for j in range(nintg):
                tmodel[i] = tmodel[i] + K * (math.cos(Pid2 + Tanom[j] - phi0) + eccn * math.cos(w))
            tmodel[i] = tmodel[i] / dnintg + voff

    return tmodel * anorm

# ----------------------------------------
# snapshot binary separation during eccentric orbit 

def distance(asep,eccn,Tanom):

    return asep * (1.0 - eccn * eccn) / (1.0 + eccn * math.cos(Tanom))

# ----------------------------------------
# true eccentric anomaly

def trueanomaly(eccn,Eanom):
      
      work1 = math.sqrt((1.0 + eccn) / (1.0 - eccn))
      work2 = math.tan(Eanom / 2.0)
      anomaly = 2.0 * math.atan(work1 * work2)
      
      return anomaly

# ----------------------------------------
# the Kepler equation

def kepler(Manom,Eanom,eccn):

    itmax=100
    thres = 1.0e-6
    Eold = copy.copy(Eanom)
    Eanom = Manom + eccn * math.sin(Eanom)
    diff = abs(1.0 - Eanom / Eold)
    Eold = copy.copy(Eanom)
    i = 0
    while diff > thres and i < itmax:
        Eanom = Manom + eccn * math.sin(Eanom)
        diff = abs(1.0 - Eanom / Eold)
        Eold = copy.copy(Eanom)
        i += 1

    return Eanom
      
# ----------------------------------------
def invkepler(Eanom,Manom,eccn):

    itmax=100
    thres = 1.0e-6
    Mold = Manom
    Manom = Eanom - eccn * math.sin(Manom)
    diff = abs(1.0 - Manom / Mold)
    Mold = Manom
    i = 0
    while diff > thres and i < itmax:
        Manom = Eanom - eccn * math.sin(Manom)
        diff = abs(1.0 - Manom / Mold)
        Mold = Manom
        i += 1
        
    return Manom

# ----------------------------------------
# tides on the surface of a binary star component

def tides(M1,M2,R1,asemi,incl,tanom,eccn,phi0):

# Eric Pfahl, Phil Arras and Bill Paxton, 2008, ApJ, 679, 783
#     M1 - stellar mass
#     M2 - companion mass
#     R1 - stellar radius
#     asemi - semi-major axis
#     incl - orbital inclination
#     tanom - true anomoly
#     eccn - orbital eccentricity
#     w - azimuth angle of perihelion (deg)

    eps = M2 / M1 * (R1 / asemi)**3.0
    d = asemi * (1.0 - eccn * eccn) / (1.0 + eccn * math.cos(tanom))
    Inc = math.pi - (incl - math.pi / 2.0)
    cos2incl = math.cos(Inc) * math.cos(Inc)
    sin2incl = math.sin(Inc) * math.sin(Inc)

    ra = numpy.zeros((4),dtype='float64')
    ad = numpy.zeros((4),dtype='float64')
    f = numpy.zeros((4),dtype='float64')
    P = numpy.zeros((4),dtype='float64')
    lambd = numpy.zeros((4),dtype='float64')
    dJ = numpy.zeros((4),dtype='float64')

    for l in range(2,4):
        ra[l] = (R1 / asemi)**(l - 2)
        ad[l] = (asemi / d)**(l + 1)
        lambd[l] = float(l) + 2.0
    f[2] = -1.3e1 * (1.0 + lambd[2] / 4.0) / 1.0e1
    f[3] = -5.0 * (1.0 + lambd[3] / 1.0e1) / 8.0
    P[2] = 2.5e-1 * (-(3.0 * cos2incl - 1.0) + 3.0 * sin2incl * math.cos(2.0 * (tanom - phi0)))
    P[3] = 1.25e-1 * math.sin(Inc) * \
        (-3.0 * (5.0 * cos2incl - 1.0) * math.cos(tanom - phi0) + 5.0 * \
              sin2incl * math.cos(3.0 * (tanom - phi0)))

    tide = 0.0
    for l in range(2,4):
        dJ[l] = ra[l] * ad[l] * f[l] * P[l]
        tide = tide + dJ[l]
    tide = tide * eps * math.pi

    return tide

# ----------------------------------------
# planet/companion albedo

def albedomod(t,Per,ag,R1,R2,phi):

    phi = phi + math.pi
    if phi > 2.0 * math.pi:
        phi = phi - 2.0 * math.pi

    alpha = abs(phi)
    alpha = alpha - 2.0 * math.pi * float(int(alpha / (2.0 * math.pi)))
    if alpha > math.pi:
        alpha = abs(alpha - 2.0 * math.pi)
    phase = (math.sin(alpha) + (math.pi - alpha) * math.cos(alpha)) / math.pi  # Lambertian Sphere

    return ag * math.pi * R2 * R2 / (R1 * R1) * phase

# ----------------------------------------
# Mandel & Agol transit model

def mandelagol(nintg,R1,R2,x1,x2,y1,y2,c):

# adapted from Mandel and Agol, 2002, ApJ 580, L171
#
# INPUT
#           nintg: number of subsampled data points within one exposure
#              R1: radius of primary star (m)
#              R2: radius of planet/companion star (m)
#              x1: position of primary star
#              y1: position of primary star
#       x2[nintg]: position of planet/companion star
#       y2[nintg]: position of planet/companion star
#            c[4]: limb darkening coefficients
#
# OUTPUT
#         managol: flux relative to unobscured source 
#  mulimb0[nintg]: flux relative to unobscured source from small body model
#  mulimbf[nintg]: flux relative to unobscured source from large body model
#     dist[nintg]: binary separation
#       b0[nintg]: impact parameters (positive number normalized to stellar radius)
#
# OTHER
#              rl: ratio of planet/companion radius to stellar radius

# startup
    
    c1 = c[0]
    c2 = c[1]
    c3 = c[2]
    c4 = c[3]
    rl = R2 / R1
    sflag = 0

# binary separation

    dist = numpy.zeros((nintg),dtype='float64')
    b0 = numpy.zeros((nintg),dtype='float64')
    for i in range(nintg):
        dist[i] = math.sqrt((x2[i] - x1) * (x2[i] - x1) + (y2[i] - y1) * (y2[i] - y1)) / (R1 + R2)
        b0[i] = (R1 + R2) * dist[i] / R1

# small body occultation model

    if sflag < nintg:
        mulimb0 = occultsmall(rl,c1,c2,c3,c4,nintg,b0)

# large body occultation model: dummy

    if sflag < nintg:
        mulimbf = numpy.zeros((nintg),dtype='float32')

    managol = 0.0
    for i in range(nintg):
        managol = managol + mulimb0[i]
    managol = managol / nintg

    return managol, mulimb0, mulimbf, dist, b0

# ----------------------------------------
# small planet occultation model

def occultsmall(p,c1,c2,c3,c4,nz,z):

# INPUT
#      p: ratio of planet radius to stellar radius
#  c1-c4: non-linear limb-darkening coefficients
#     nz: number of subsampled data points within one exposure
#  z[nz]: impact parameters (positive number normalized to stellar radius)
#
# OUTPUT:
# mu[nz]: flux relative to unobscured source for each z

    mu = numpy.zeros((nz),dtype='float64')
    norm = math.pi * (1.0 - c1 / 5.0 - c2 / 3.0 - 3.0 * c3 / 7.0 - c4 / 2.0)
    i1 = 1.0 - c1 - c2 - c3 - c4
    for i in range(nz):
        mu[i] = 1.0
        if z[i] > 1.0 - p and z[i] < 1.0 + p:
            x = 1.0 - (z[i] - p)**2
            tmp = (1.0 - c1 * (1.0 - 0.8 * x**0.25) \
                       - c2 * (1.0 - 2.0 / 3.0 * x**0.5) \
                       - c3 * (1.0 - 4.0 / 7.0 * x**0.75) \
                       - c4 * (1.0 - 0.5 * x))
            mu[i] = 1.0 - tmp * (p**2 * math.acos((z[i] - 1.0) / p) \
                                     - (z[i] - 1.0) * math.sqrt(p**2 - (z[i] - 1.0)**2)) / norm
        if z[i] < 1.0 - p and z[i] != 0.0:
            mu[i] = 1.0 - math.pi * p**2 * iofr(c1,c2,c3,c4,z[i],p) / norm
        if z[i] == 0.0:
            mu[i] = 1.0 - math.pi * p**2 / norm

    return mu

# ----------------------------------------
# small planet occultation model

def iofr(c1,c2,c3,c4,r,p):

    sig1 = math.sqrt(math.sqrt(1.0 - (r - p)**2))
    sig2 = math.sqrt(math.sqrt(1.0 - (r + p)**2))
    return 1.0 - c1 * (1.0 + (sig2**5 - sig1**5) / 5.0 / p / r) \
        - c2 * (1.0 + (sig2**6 - sig1**6) / 6.0 / p / r) \
        - c3 * (1.0 + (sig2**7 - sig1**7) / 7.0 / p / r) \
        - c4 * (p**2 + r**2)

# ----------------------------------------
# small planet occultation model

def eclmod2(R1,R2,x1,x2,y1,y2,zarea,norm,ted):

# put everything in upper (positive) quadrant to make life simple
 
     ax1=abs(x1)
     ax2=abs(x2)
     ay1=abs(y1)
     ay2=abs(y2)
     
# distance between projected center of planet and star
  
     dist = math.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1))
      
# if dist < R1 + R2 then !we have an eclipse 

     sflag = radsolve(R1,R2,x1,x2,y1,y2)
     if sflag == 0:

# get x,y which is the intercept of the stellar radius with a diameter
# of the stellar radius pointed towards the stellar radius origin

         sinth = abs(ax2 - ax1) / dist
         x = ax1 + abs(R1 * sinth)
         y = ay1 + math.sqrt(R1 * R1 - (x - ax1) * (x - ax1))
         d1 = math.sqrt((x - ax1) * (x - ax1) + (y - ay1) * (y - ay1))
         d2 = math.sqrt((x - ax2) * (x - ax2) + (y - ay2) * (y - ay2))
         if d1 <= dist:
             ratio = (R2 - d2) / (2.0 * R2)
         else:
             ratio = (R2 + d2) / (2.0 * R2)
         eclmod = (math.pi + zarea) / norm - ted * ratio
     else:
         y2p = math.sqrt(x2 * x2 + y2 * y2)
         if y2p / R1 + R2 / R1 < 1.0:
             eclmod = (math.pi + zarea) / norm - ted #  inside transit
         else:
             eclmod = (math.pi + zarea) / norm       # outside transit
      
     return eclmod

# ----------------------------------------
# small planet occultation model

def radsolve(R1,R2,x1,x2,y1,y2):
      
    a = x1 - x2
    b = y1 - y2
    temp1 = b * b * ((R1 - R2) * (R1 - R2) - a * a - b * b) * (-(R1 + R2) * (R1 + R2) + a * a + b * b)
    if temp1 > 0.0:
        sflag = 0
        temp2 = 2.0 * (a * a + b * b)
        temp3 = -R1 * R1 * a + R2 * R2 * a + (x1 + x2)*(a * a + b * b)
        x01 = (temp3 + math.sqrt(temp1)) / temp2
        x02 = (temp3 - math.sqrt(temp1)) / temp2
    else:
        sflag = 1

    return sflag
