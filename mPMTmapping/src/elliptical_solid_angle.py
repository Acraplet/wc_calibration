"""
Created on Fri Dec  8 16:58:00 2017

Copyright (c) 2017-19-22, Sam Lambrick
All rights reserverd.
Subject to the MIT licence.

Contains the nessacery functions to perform analytical solid angle
integrals for ellipses. Also contains a function for integrating and
calculating the solid angle.

The derivation of the functions was performed by Abbas et al. 2015, for an
understanding of the equations used it is recommened to look there. The
variables follow the convention used there.
doi:10.1016/j.nima.2014.10.061

The derivation of the simpler general function (using vector potentials) was
done by John T. Conway 2010. His formual should be used in general, however the
geometric derivation of Abbas et al. can be combined more readily with other
functions.
doi:10.1016/j.nima.2009.11.075

The theta integral from Abbas et al. can be performed analytically. This has
been used here along with the double integral form. The arguments p and q are
distances and should be positive.

Variables used here, using the same naming convention as used by Abbas:
    p - the distance between the centre of the ellipse and the point of
        interest in the plane of the ellipse along the major axis of the
        ellipse
    q - the distance between the centre of the ellipse and the point of
        interest in the plane of the ellipse along the minor axis of the
        ellipse
    a - the semi-major axis of the ellipse
    b - the semi-minor axis of the ellipse
    h - the vertical perpendicular distance between the point of interest and
        the plane of the ellpise
    rho - sqrt(p^2 + q^2), the distance from the centre of the ellipse to the
          point of interest projected into the plane of the ellipse

Solid angle functions:
    calc_omega        - performs a single integration to calculate the
                        solid angle using integral forms from John T. Conway
    solid_angle_calc  - performs a single integration to calculate the solid
                        angle using the integral forms from Abbas et al. with
                        the theta integral performed
    solid_angle_calc2 - performs a double inetegration, using the form directly
                        from Abbas et al., to calculate the solid angles

Calling syntax:
    solid_angle_calc...(p, q, a, b, h)

None of these functions return the error estimate on the integrals, however
it is easy to modify them to do so.

This module makes use of numeric python (numpy) and scientific python (scipy).
"""

import numpy as np
import scipy.integrate as integrate
import sys


def sec(theta):
    """Calculates the sectant, 1/cos(theta)."""
    return(1/np.cos(theta))


## Geometric functions for Abbas ---------------------------------------------#


def r_1(phi, a, b):
    """Calculates the distance from the point of interest to the edge of the
    ellipse, in the plane of the ellipse. For the case p = 0."""

    denominator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = a*b/np.sqrt(denominator)
    return(result)


def r_2(p, phi, a, b):
    """Calculated the distance from the point of interest to the edge of the
    ellipse, in the plane of the ellipse. For the case p < a."""

    numerator1 = -p*b*b*np.cos(phi)
    numerator2 = a*b*np.sqrt((b*np.cos(phi))**2 + (a*a - p*p)*(np.sin(phi))**2)
    denominator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = (numerator1 + numerator2)/denominator
    return(result)


def r_3(p, phi, a, b):
    """Calculates the distance from the point of interest to the edge of the
    ellipse, in the plane of the ellipse. For the case of p = a."""

    numerator = 2*a*b*b*np.cos(phi)
    denominator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = numerator/denominator
    return(result)


# This is the negative of what it is in the paper
def r_4(p, phi, a, b):
    """Calculates the distance from the point of interest to the edge of the
    ellipse projected in the plane of the ellipse. For the case p > a."""

    numerator1 = a*a*p*(np.sin(phi))**2 * sec(phi)
    numerator2 = a*b*np.sqrt((b*np.cos(phi))**2 + (a*a-p*p)*(np.sin(phi))**2)
    demoninator = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    result = ((numerator1 + numerator2)/demoninator) - p*sec(phi)
    return(-result)


def r_5(p, phi, a, b):
    """Calculates the distance between the two interceps of the ellipse on the
    line of integration."""

    numerator = 2*a*b*np.sqrt((b*np.cos(phi))**2 + (a*a-p*p)*(np.sin(phi))**2)
    denominator = a**2 * (np.sin(phi))**2 + b**2 * (np.cos(phi))**2
    result = numerator/denominator
    return(result)


def r_6(p, q, phi, a, b):
    """Calculates the distance between the point of interest to the edge of the
    ellipse for the case rho < a, rho < b."""

    tmp = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    numerator1 = p*b*b*np.cos(phi) + q*a*a*np.sin(phi)
    numerator2 = a*b*np.sqrt(tmp - (p*np.tan(phi) - q)**2 * (np.cos(phi))**2)
    denominator = tmp
    result = (-numerator1 + numerator2)/denominator
    return(result)


def r_8(p, q, phi, a, b):
    """Calculates the distance between the point of interest and the edge of
    the ellipse (projected into the plane of the ellipse) along the line of
    integration."""

    tmp = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    numerator1 = a*a*np.sin(phi)*(p*np.tan(phi) - q)
    numerator2 = a*b*np.sqrt(tmp - (p*np.sin(phi) - q*np.cos(phi))**2)
    denominator = tmp
    result = (numerator1 + numerator2)/denominator - p*sec(phi)
    return(-result)


def r_9(p, q, phi, a, b):
    """Calculated the distance from the edge of the ellipse to the other edge
    along the line of integration for the case rho > a, rho > b, for a specific
    value of phi."""

    tmp = (a*np.sin(phi))**2 + (b*np.cos(phi))**2
    numerator = 2*a*b*np.sqrt(tmp - (p*np.sin(phi) - q*np.cos(phi))**2)
    denominator = tmp
    result = numerator/denominator
    return(result)


def theta_1(h, phi, a, b):
    """Upper limit of the theta integral for the case p = 0."""
    result = np.arctan(r_1(phi, a, b)/h)
    return(result)


def theta_2(p, h, phi, a, b):
    """Upper limit of the theta integral for the case p < a."""
    result = np.arctan(r_2(p, phi, a, b)/h)
    return(result)


def theta_3(p, h, phi, a, b):
    """Upper limit of the theta integral for the case p = a."""
    result = np.arctan(r_3(p, phi, a, b)/h)
    return(result)


def theta_4(p, h, phi, a, b):
    """Lower limit of the theta integral for the case p > a."""
    result = np.arctan(r_4(p, phi, a, b)/h)
    return(result)


def theta_5(p, h, phi, a, b):
    """Upper limit of the theta integral for the case p > a."""
    result = np.arctan((r_4(p, phi, a, b) + r_5(p, phi, a, b))/h)
    return(result)


def theta_6(p, q, h, phi, a, b):
    """Upper limit of integration for the case rho < a, rho < b."""
    result = np.arctan(r_6(p, q, phi, a, b)/h)
    return(result)


def theta_8(p, q, h, phi, a, b):
    """Lower limit of integration for the case rho > a, rho > b."""
    result = np.arctan(r_8(p, q, phi, a, b)/h)
    return(result)


def theta_9(p, q, h, phi, a, b):
    """Upper limit of integration for the case rho > a, rho > b."""
    result = np.arctan((r_8(p, q, phi, a, b) + r_9(p, q, phi, a, b))/h)
    return(result)


## Integrands for Abbas ------------------------------------------------------#

def integrand_1(phi, h, a, b):
    """Solid angle integrand for the case p = 0."""
    result = np.cos(theta_1(h, phi, a, b))
    return(result)


def integrand_2(phi, p, h, a, b):
    """Solid angle integrand for the case p < a."""
    result = np.cos(theta_2(p, h, phi, a, b))
    return(result)


def integrand_3(phi, p, h, a, b):
    """Solid angle integrand for the case p = a."""
    result = 2*np.cos(theta_3(p, h, phi, a, b))
    return(result)


def integrand_4(phi, p, h, a, b):
    """Solid angle integrand for the case p > a."""
    result = 2*(np.cos(theta_4(p, h, phi, a, b)) - np.cos(theta_5(p, h, phi,
                                                                  a, b)))
    return(result)


def integrand_5(phi, p, q, h, a, b):
    """Solid angle integrand for the case rho < a, rho < b."""
    result = np.cos(theta_6(p, q, h, phi, a, b))
    return(result)


def integrand_6(phi, p, q, h, a, b):
    """Solid angle integrand for the case rho > a, rho > b."""
    result = np.cos(theta_8(p, q, h, phi, a, b)) - np.cos(theta_9(p, q, h, phi,
                                                                  a, b))
    return(result)


def phi_max_calc(p, a, b):
    """Upper limit of the phi integral for the case p > a."""
    y = np.sqrt(b**2/(p**2 - a**2))
    result = np.arctan(y)
    return(result)


def phi_prime_max_calc(p, q, a, b):
    """Lowwer limit of the phi integral for the case rho > a, rho > b."""
    arg_numerator = p*q - np.sqrt(-a*a*b*b + b*b*p*p + a*a*q*q)
    result = np.arctan(arg_numerator/(p*p - a*a))
    return(result)


def phi_primeprime_max_calc(p, q, a, b):
    """Upper limit of the phi integral for the case rho > a, rho > b."""
    arg_numerator = p*q + np.sqrt(-a*a*b*b + b*b*p*p + a*a*q*q)
    result = np.arctan(arg_numerator/(p*p - a*a))
    return(result)

## Abbas solid angle calculations --------------------------------------------#

def solid_angle_calc(p, q, h, a, b):
    """Calculates the solid angle subtended by an ellipse at a point on axis of
    the semi-axis of the ellipse. This function is not vectorized as it uses
    scipy.integrate.quad.

    Inputs:
     p - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the ellipse along the major axis.
         Should not be negative.
     q - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the elllipse along the minor
         axis. Should not be negative
     h - The perpendicular distance between the point of interest and the plane
         of the ellipse. Should not be negative.
     a - The semi-axis of the ellipse along which the point of interest lies
         when it is projected into the plane of the ellipse.
     b - The other semi-axis of the ellipse.

    Output:
     omega - The solid angle subtended by the ellipse from the point of
             interest
    """

    if p == 0 and q != 0:
        # No need to preform the more complicated integral if the point of
        # interest lies along either axis of the ellipse
        p = q
        q = 0
        tmp = a
        a = b
        b = tmp

    if q == 0:
        # The cases of the point of interest lying along an axis of the ellipse
        if p == 0:
            omega = np.pi - integrate.quad(integrand_1, 0, np.pi,
                                           args=(h, a, b,))[0]
            omega = 2*omega
        elif p < a:
            omega = np.pi - integrate.quad(integrand_2, 0, np.pi,
                                           args=(p, h, a, b,))[0]
            omega = 2*omega
        elif p == a:
            omega = np.pi - integrate.quad(integrand_3, 0, np.pi/2,
                                           args=(p, h, a, b,))[0]
        else:
            phi_max = phi_max_calc(p, a, b)
            omega = integrate.quad(integrand_4, 0, phi_max,
                                   args=(p, h, a, b,))[0]
    else:
        # The case of the point of interest not lying along an axis of the
        # ellispe
        if (p**2/a**2 + q**2/b**2) < 1:
            # Point inside the ellipse
            omega = 2*np.pi - integrate.quad(integrand_5, 0, 2*np.pi,
                                             args=(p, q, h, a, b,))[0]
        elif p > a or q > b:
            # Point outiside of the bounding box
            if q > b and p <= a:
                tmp = q
                q = p
                p = tmp
                tmp = a
                a = b
                b = tmp
            phi_min = phi_prime_max_calc(p, q, a, b)
            phi_max = phi_primeprime_max_calc(p, q, a, b)
            omega = integrate.quad(integrand_6, phi_min, phi_max,
                                   args=(p, q, h, a, b,))[0]
        else:
            # Point outside of the ellipse but inside the bounding box
            # TODO: derive expression for this
            omega = solid_angle_calc_general(p, q, h, a, b)

    return(omega)


def solid_angle_calc2(p, q, h, a, b):
    """Calculates the solid angle subtended by an ellipse at a point on axis of
    the semi-axis of the ellipse. Uses a double integration, recommended to use
    solid_angle_calc instead. This function is not vectorized as it uses
    scipy.integrate.quad.

    Inputs:
     p - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the ellipse along the major axis
         of the ellipse. Should not be negative.
     q - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the ellipse, along the minor axis
         of the ellipse. Should not be negative
     h - The perpendicular distance between the point of interest and the plane
         of the ellipse. Should not be negative.
     a - The semi-major-axis of the ellipse.
     b - The semi-minor-axis of the ellipse.

    Output:
     omega - The solid angle subtended by the ellipse from the point of
             interest
    """

    if p == 0 and q != 0:
        # No need to preform the more complicated integral if the point of
        # interest lies along either axis of the ellipse
        p = q
        q = 0
        tmp = a
        a = b
        b = tmp

    if q == 0:
        # The cases of the point of interest lying along an axis of the ellipse
        if p == 0:
            th_1 = lambda y: theta_1(h, y, a, b)
            integrand = lambda y,z: np.sin(y)
            omega = 2*integrate.dblquad(integrand, 0, np.pi, lambda x: 0,
                                        th_1)[0]
        elif p < a:
            th_2 = lambda y: theta_2(p, h, y, a, b)
            integrand = lambda y,z: np.sin(y)
            omega = 2*integrate.dblquad(integrand, 0, np.pi, lambda x: 0,
                                        th_2)[0]
        elif p == a:
            th_3 = lambda y: theta_3(p, h, y, a, b)
            integrand = lambda y,z: 2*np.sin(y)
            omega = integrate.dblquad(integrand, 0, np.pi/2, lambda x: 0,
                                      th_3)[0]
        else:
            phi_max = phi_max_calc(p, a, b)
            th_4 = lambda y: theta_4(p, h, y, a, b)
            th_5 = lambda y: theta_5(p, h, y, a, b)
            integrand = lambda y,z: 2*np.sin(y)
            omega = integrate.dblquad(integrand, 0, phi_max, th_4, th_5)[0]
    else:
        # The case of the point of interest not lying along an axis of the
        # ellispe
        if (p**2/a**2 + q**2/b**2) < 1:
            # Point inside the ellipse
            th_6 = lambda y: theta_6(p, q, h, y, a, b)
            integrand = lambda y,z: np.sin(y)
            omega = integrate.dblquad(integrand, 0, 2*np.pi, lambda x: 0,
                                      th_6)[0]
        elif p > a or q > b:
             # Point outiside of the bounding box
            if q > b and p <= a:
                tmp = q
                q = p
                p = tmp
                tmp = a
                a = b
                b = tmp
            th_8 = lambda y: theta_8(p, q, h, y, a, b)
            th_9 = lambda y: theta_9(p, q, h, y, a, b)
            phi_min = phi_prime_max_calc(p, q, a, b)
            phi_max = phi_primeprime_max_calc(p, q, a, b)
            integrand = lambda y,z: np.sin(y)
            omega = integrate.dblquad(integrand, phi_min, phi_max, th_8,
                                      th_9)[0]
        else:
            # Point outside of the ellipse but inside the bounding box
            # TODO: derive and implement this
            omega = solid_angle_calc_general(p, q, h, a, b)

    return(omega)


## Cosine modified integrals -------------------------------------------------#


# Two functions are included, one that performs a double numerical integral
# and the other that has had the theta integral performed analytically. The
# double integral is included for completness, and as a demonstartion of how the
# inclusion of distributions can be achieved.


def cosine_intensity(p, h, a, b):
    """Calculates the intensity of a cosine distribution that enters an
    elliptical appeture."""

    if p < a:
        th_2 = lambda y: theta_2(p, h, y, a, b)
        integrand = lambda y: np.cos(2*th_2(y))
        II = np.pi/2 - 0.25*integrate.quad(integrand, 0, 2*np.pi)[0]
    elif p == a:
        th_3 = lambda y: theta_3(p, h, y, a, b)
        integrand = lambda y: np.cos(2*th_3(y))
        II = np.pi/4 - 0.5*integrate.quad(integrand, 0, np.pi/2)[0]
    else:
        phi_max = phi_max_calc(p, a, b)
        th_4 = lambda y: theta_4(p, h, y, a, b)
        th_5 = lambda y: theta_5(p, h, y, a, b)
        integrand = lambda y: np.cos(2*th_4(y)) - np.cos(2*th_5(y))
        II = 0.5*integrate.quad(integrand, 0, phi_max)[0]

    return(II)


def cosine_intensity2(p, h, a, b):
    """Calculates the intensity of a cosine distribution that enters an
    elliptical appeture. Uses a double integral."""

    if p == 0:
        th_1 = lambda y : theta_1(h, y, a, b)
        integrand = lambda y,z : np.sin(y)*np.cos(y)
        II = 2*integrate.dblquad(integrand, 0, np.pi, lambda x: 0, th_1)[0]
    elif p < a:
        th_2 = lambda y : theta_2(p, h, y, a, b)
        integrand = lambda y,z : np.sin(y)*np.cos(y)
        II = 2*integrate.dblquad(integrand, 0, np.pi, lambda x: 0, th_2)[0]
    elif p == a:
        th_3 = lambda y: theta_3(p, h, y, a, b)
        integrand = lambda y,z : 2*np.sin(y)*np.cos(y)
        II = integrate.dblquad(integrand, 0, np.pi/2, lambda x: 0, th_3)[0]
    else:
        phi_max = phi_max_calc(p, a, b)
        th_4 = lambda y : theta_4(p, h, y, a, b)
        th_5 = lambda y : theta_5(p, h, y, a, b)
        integrand = lambda y,z : 2*np.sin(y)*np.cos(y)
        II = integrate.dblquad(integrand, 0, phi_max, th_4, th_5)[0]

    return(II)


## Conway method functions ---------------------------------------------------#

def integrand_general(phi, p, q, h, a, b):
    """Integrand for the general case from 'John T. Conway' 'Analytic solution
    for the solid angle subtended at a point source radiation vector potential'
    which claims to be general, unlike Abbas et.al. who do not consider the
    case inside the bounding box but outisde the ellipse."""

    tmp = p*p + q*q + h*h + 2*a*p*np.cos(phi) + 2*b*q*np.sin(phi) + \
        (a*np.cos(phi))**2 + (b*np.sin(phi))**2
    term1 = 1 - h/(np.sqrt(tmp))
    tmp = tmp - h*h
    term2 = (a*b + p*b*np.cos(phi) + q*a*np.sin(phi))/(tmp)
    result = term1*term2
    return(result)


def calc_omega(p, q, h, a, b):
    """Calculates the solid angle using the formula derived by John T. Conway
    2010. Is faster than the alternative integral in 'solid_angle_calc'.

    Inputs:
     p - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the ellipse along the major axis.
         Should not be negative.
     q - The distance between the point of interest and the centre of the
         ellipse, projected into the plane of the elllipse along the minor
         axis. Should not be negative
     h - The perpendicular distance between the point of interest and the plane
         of the ellipse. Should not be negative.
     a - The semi-axis of the ellipse along which the point of interest lies
         when it is projected into the plane of the ellipse.
     b - The other semi-axis of the ellipse.

    Output:
     omega - The solid angle subtended by the ellipse from the point of
             interest
    """

    omega = integrate.quad(integrand_general, 0, 2*np.pi,
                           args=(p, q, h, a, b,))[0]
    return(omega)


def projectSourceOntoPMT(sourcePosition, PMTPosition, PMTdirection):
    """Source position, PMT position and direction are 3-vectors
    This is projecting the source position onto the axis of the mPMT
    First define the plane Ax+By+Cz = D of the PMT using its centre and the pointing direction,
    then calculate the projection onto the plane using the PMT direction
    finally, calculate p and q, projection of that point onto two axes that just
    have to be perpendicular to one another
    """
    #need to normalise the PMT direction, otherwise the reconstruction is wrong
    PMTdirection = PMTdirection/np.linalg.norm(PMTdirection)
    A, B, C = PMTdirection[0], PMTdirection[1], PMTdirection[2]
    D = A * PMTPosition[0] + B * PMTPosition[1] + C * PMTPosition[2]

    #the disance that we need to project the source along the negative of the PMT direction
    #for to reach the plane of the PMT
    numerator = (D - A * sourcePosition[0] - B * sourcePosition[1] - C * sourcePosition[2])
    # negative sign because we go the other direction w.r.t the PMT axis
    denominator = - (A * A + B * B + C * C)
    #the perpendicular distance to the plane of the PMT
    R = numerator/denominator
    if R <=10**-5:
        raise ValueError('The source-PMT distance is negative or zero (R=%.2f), is the source inside the PMT?'%R)
    projectionSource = np.array([ - R * PMTdirection[0], - R * PMTdirection[1], - R * PMTdirection[2] ]) + sourcePosition

    #project the point along a single axis and decide the other one is 0 (this is a circle - same dims)
    p = np.linalg.norm(PMTPosition - projectionSource)
    if p <0:
        raise ValueError('p should not be negative, current value: %.2f'%p)
    q = 0
    return R, p, q



### take a circular disk right below the point of interest
p = 0
q = 0
h = 1 #in cm
a = 4 #in cm
b = a #this is a disk, same distance
solid_angle  = calc_omega(p, q, h, a, b)
# print("Solid angle for this disk of radius %.2f at a distance %.2f from the point is %.2f"%(a, h, solid_angle))


#we call the source and PMT position and the PMT direction, from the cpp input
#TODO: cleanup this with some flag inputs instead and a helper function
sourcePosition = np.array([float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])])
PMTPosition = np.array([float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])])
PMTdirection = np.array([float(sys.argv[7]), float(sys.argv[8]), float(sys.argv[9])])

# sourcePosition = np.array([0, 10,0])
# PMTPosition = np.array([0,0,0])
# PMTdirection = np.array([0,1,0])

PMTdirection = PMTdirection/np.linalg.norm(PMTdirection)

R, p, q = projectSourceOntoPMT(sourcePosition, PMTPosition, PMTdirection)

# print("Distance between the source and the projected point: %.2f" %(R))
# print("Distance of the projected point along the axis: %.2f" %(p))
solid_angle  = calc_omega(p, q, R, a, b)

print(solid_angle)
# solid_angle2  = solid_angle_calc2(p, q, R, a, b)
#
# print("Solid angle for this disk of radius %.2f, poiting towards %s at a distance %.2f from the point is %.2f"%(a, PMTdirection, R, solid_angle))
# print("Solid angle for this disk of radius %.2f, poiting towards %s at a distance %.2f from the point is %.2f"%(a, PMTdirection, R, solid_angle2))
