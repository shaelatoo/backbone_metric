#!/usr/bin/env python3
#
#
# ==========================================================================
# 
# angular_separation.py
# 
# ==========================================================================
# 06/16/2021    sij    created 
#
#-----------------------------------------------------------------------------
""" calculates the angular separation between two points on a sphere """

import math
import numpy as np

def angular_separation(phi0, theta0, phi1, theta1, degrees=False, cossep=False):
    """ calculates distance between two points on a sphere, in radians """


    if degrees:    # convert to radians, if needed
        first_term = math.sin(math.radians(theta0))*math.sin(math.radians(theta1))
        second_term = math.cos(math.radians(theta0)) * math.cos(math.radians(theta1)) * \
                  math.cos(math.radians(phi0-phi1))
    else:
        first_term = math.sin(theta0)*math.sin(theta1)
        second_term = math.cos(theta0) * math.cos(theta1) * math.cos(phi0-phi1)

    try:
        dist = np.arccos(first_term + second_term)
    except:
        print(first_term + second_term)

    if cossep:     # return the cosine of the separation instead
        return first_term + second_term

    return dist
