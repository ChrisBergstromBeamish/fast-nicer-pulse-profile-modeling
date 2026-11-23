# -*- coding: utf-8 -*-
"""Constants (cgs and keV). GM of the Sun is an observed value; h and c are exact by definition of the units."""
H_KEV_S = 4.135667696e-18          # Planck's constant, keV s
C_CM_S = 2.99792458e10             # speed of light, cm/s
GM_SUN = 1.32712440018e26          # cm^3/s^2
PC_CM = 3.0856775814913673e18      # parsec, cm
BB_CONST = 2.0 / (H_KEV_S**3 * C_CM_S**2)   # photons / (cm^2 s sr keV^3): black-body photon intensity is BB_CONST E^2/(exp(E/kT)-1)
