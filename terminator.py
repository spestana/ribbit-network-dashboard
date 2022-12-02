
# From https://github.com/joergdietrich/Leaflet.Terminator
# re-written from javascript to python

import numpy as np
import pandas as pd
import json

def julian(datetime):
    ''' Calculate the present UTC Julian Date. Function is valid after
        the beginning of the UNIX epoch 1970-01-01 and ignores leap
        seconds.'''
    return (datetime / 86400000) + 2440587.5
	
def GMST(julianDay):
    ''' Calculate Greenwich Mean Sidereal Time according to 
        http://aa.usno.navy.mil/faq/docs/GAST.php '''
    d = julianDay - 2451545.0
    # Low precision equation is good enough for our purposes.
    return (18.697374558 + 24.06570982441908 * d) % 24
	
def sunEclipticPosition(julianDay):
    ''' Compute the position of the Sun in ecliptic coordinates at
        julianDay.  Following
        http://en.wikipedia.org/wiki/Position_of_the_Sun '''
    # Days since start of J2000.0
    n = julianDay - 2451545.0
    # mean longitude of the Sun
    L = 280.460 + 0.9856474 * n
    L %= 360
    # mean anomaly of the Sun
    g = 357.528 + 0.9856003 * n
    g %= 360
    # ecliptic longitude of Sun
    l = L + 1.915 * np.sin(np.radians(g)) + 0.02 * np.sin(np.radians(2 * g))
    # distance from Sun in AU
    R = 1.00014 - 0.01671 * np.cos(np.radians(g)) - 0.0014 * np.cos(np.radians(2 * g))
    return l, R
	
def eclipticObliquity(julianDay):
    ''' Following the short term expression in
        http://en.wikipedia.org/wiki/Axial_tilt#Obliquity_of_the_ecliptic_.28Earth.27s_axial_tilt.29 '''
    n = julianDay - 2451545.0;
    # Julian centuries since J2000.0
    T = n / 36525
    epsilon = 23.43929111 - \
            T * (46.836769 / 3600 - \
                     T * (0.0001831 / 3600 + \
                          T * (0.00200340 / 3600 - \
                               T * (0.576e-6 / 3600 - \
                                    T * 4.34e-8 / 3600))))
    return epsilon
	
def sunEquatorialPosition(l, epsilon):
    '''Compute the Sun's equatorial position from its ecliptic
       position. Inputs are expected in degrees. Outputs are in
       degrees as well. '''
    alpha = np.degrees( np.arctan(np.cos(np.radians(epsilon)) * np.tan(np.radians(l))) )
    delta = np.degrees( np.arcsin(np.sin(np.radians(epsilon)) * np.sin(np.radians(l))) )
    
    lQuadrant = np.floor(epsilon / 90) * 90
    raQuadrant = np.floor(alpha / 90) * 90
    alpha = alpha + (lQuadrant - raQuadrant)
    
    return alpha, delta
	
def hourAngle(lng, alpha, gst):
    ''' Compute the hour angle of the sun for a longitude on
        Earth. Return the hour angle in degrees. '''
    lst = gst + lng / 15
    return lst * 15 - alpha
	
def latitude(ha, delta):
    ''' For a given hour angle and sun position, compute the
        latitude of the terminator in degrees. '''
    lat = np.degrees( np.arctan( -1*np.cos(np.radians(ha)) / np.tan(np.radians(delta)) ) )
    return lat
	
def get_terminator(datetime, resolution):
    julianDay = julian(datetime)
    gst = GMST(julianDay)
    latLng = []
    
    [l, R] = sunEclipticPosition(julianDay)
    epsilon = eclipticObliquity(julianDay)
    alpha, delta = sunEquatorialPosition(l, epsilon)
    
    latLng = np.zeros((len(np.arange(0,720,resolution)),2))
    c = 0
    for i in np.arange(0,720,resolution):
        lng = 360 - i / resolution
        ha = hourAngle(lng, alpha, gst)
        latLng[c] = [-latitude(ha, delta), lng]
        c+=1

    if delta < 0:
       latLng[0] = [90, 360];
       latLng[-1] = [90, -360];
    else:
       latLng[0] = [-90, 360];
       latLng[-1] = [-90, -360];

    terminator = {"type":"FeatureCollection","features":[{"type":"Feature","properties":{},"geometry":{"coordinates":[np.flip(latLng, 1).tolist()],"type":"Polygon"}}]}
    return terminator