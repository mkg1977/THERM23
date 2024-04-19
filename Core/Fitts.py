#!/usr/bin/env python
# coding: utf-8
# **********************************************
# Python modules to work with 
# **********************************************

import os
import sys
import glob
import numpy as np
import pandas as pd
import timeit
import subprocess
import matplotlib.pyplot as plt
import cython
import math
# import constants
# from exception import InvalidThermoModelError
# ..............................................
# from module import ...
from datetime import date
from statistics import mean
from sympy.core.symbol import symbols
from sympy.solvers.solveset import nonlinsolve
from scipy import optimize
from scipy.optimize import curve_fit
from scipy import interpolate


class WilhoitModel:
    """
    A thermodynamics model based on the Wilhoit equation for heat capacity,
    """

    def __init__(self, cp0=0.0, cpInf=0.0, a0=0.0, a1=0.0, a2=0.0, a3=0.0, H0=0.0, S0=0.0, B=500.0):
        self.cp0 = cp0
        self.cpInf = cpInf
        self.B = B
        self.a0 = a0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.H0 = H0
        self.S0 = S0

    def __repr__(self):
        """
        Return a string representation that can be used to reconstruct the
        object.
        """
        return 'cp0=%g, cpInf=%g, a0=%g, a1=%g, a2=%g, a3=%g, H0=%g, S0=%g, B=%g' % (
            self.cp0, self.cpInf, self.a0, self.a1, self.a2, self.a3, self.H0, self.S0, self.B)

    def getHeatCapacities(self, Tlist):
        return np.array([self.getHeatCapacity(T) for T in Tlist], np.float64)

    def getEnthalpies(self, Tlist):
        return np.array([self.getEnthalpy(T) for T in Tlist], np.float64)

    def getEntropies(self, Tlist):
        return np.array([self.getEntropy(T) for T in Tlist], np.float64)

    def getHeatCapacity(self, T):
        """
        Return the Cps for a polynomial grade 3 for extrapolation. (This is used by default)
        """
        cython.declare(y=cython.double)
        y = T / (T + self.B)
        return self.cp0 + (self.cpInf - self.cp0)*y**3*(1 + (y-1) * (self.a0 + y * (self.a1 + y * (self.a2 + y * self.a3))) )
        #return self.cp0 + (self.cpInf - self.cp0)* y * y * (1 + (y-1) * (self.a0 + y * (self.a1 + y * (self.a2 + y * self.a3))) )

    def getHeatCapacity2(self, T):
        """
        Return the Cps for a polynomial grade 4 for extrapolation. (this one is for special cases only)
        """
        cython.declare(y=cython.double)
        y = T / (T + self.B)
        #return self.cp0 + (self.cpInf - self.cp0)*y**3*(1 + (y-1) * (self.a0 + y * (self.a1 + y * (self.a2 + y * self.a3))) )
        #return self.cp0 + (self.cpInf - self.cp0)* y * y * (1 + (y-1) * (self.a0 + y * (self.a1 + y * (self.a2 + y * self.a3))) )
        return self.cp0 + (self.cpInf - self.cp0)*y**4 * (1 + (y-1) * (self.a0 + y * (self.a1 + y * (self.a2 + y * self.a3))) )


    def getEnthalpy(self, T):
        """
        Return the enthalpy in J/mol at the specified temperature `T` in
        K.
        """
        cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double,
                       a2=cython.double, a3=cython.double)
        cython.declare(y=cython.double, y2=cython.double, logBplust=cython.double)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        y2 = y * y
        logBplust = math.log(B + T)
        return self.H0/1000 + cp0 * T - (cpInf - cp0) * T * (y2 * ((3 * a0 + a1 + a2 + a3) / 6.0 +(4 * a1 + a2 + a3) * y / 12.0 + (5 * a2 + a3) * y2 / 20. + a3 * y2 * y / 5.0) + (2 + a0 + a1 + a2 + a3) * (y / 2.0 - 1 + (1 / y - 1) * logBplust))

    def getEntropy(self, T):
        """
        Return the entropy in J/mol*K at the specified temperature `T` in
        K.
        """
        cython.declare(cp0=cython.double, cpInf=cython.double, B=cython.double, a0=cython.double, a1=cython.double,
                       a2=cython.double, a3=cython.double)
        cython.declare(y=cython.double, logt=cython.double, logy=cython.double)
        cp0, cpInf, B, a0, a1, a2, a3 = self.cp0, self.cpInf, self.B, self.a0, self.a1, self.a2, self.a3
        y = T / (T + B)
        logt = math.log(T)
        logy = math.log(y)
        return self.S0 + cpInf * logt - (cpInf - cp0) * (logy + y * (1 + y * (a0 / 2 + y * (a1 / 3 + y * (a2 / 4 + y * a3 / 5)))))

    def residuo(self, B, Tlist, Cplist, linear, nFreq, nRotors, H298, S298):
        # The residual corresponding to the fitToData() method
        # Parameters are the same as for that method
        cython.declare(Cp_fit=np.ndarray)
        self.fitToDataForConstantB(Tlist, Cplist, linear, nFreq, nRotors, B, H298, S298)
        Cp_fit = self.getHeatCapacities(Tlist)
        # Objective function is linear least-squares
        return np.sum((Cp_fit - Cplist) * (Cp_fit - Cplist))

    def fitToData(self, Tlist, Cplist, linear, nFreq, nRotors, H298, S298, B0):
        """
        Fit a Wilhoit model to the data points provided, allowing the
        """
        self.B = B0
        import scipy.optimize
        #scipy.optimize.fminbound(self.residuo, 300, 2000, args=(Tlist, Cplist, linear, nFreq, nRotors, H298, S298))
        scipy.optimize.fminbound(self.residuo, 500, 3500, args=(Tlist, Cplist, linear, nFreq, nRotors, H298, S298))
        return self

    def fitToDataForConstantB(self, Tlist, Cplist, linear, nFreq, nRotors, B, H298, S298):
        """
        Fit a Wilhoit model to the data points provided using a specified value
        """

        cython.declare(y=np.ndarray, A=np.ndarray, b2=np.ndarray, x=np.ndarray)
        R = 1.9872 ## cal;
        # Set the Cp(T) limits as T -> and T -> infinity
        if  linear in ["yes", "Y", "y", "YES", "YEs", "Yes", "yeS", "yES"]:
            self.cp0   = 3.5 * R
            #self.cpInf = (self.cp0 + (nFreq + 1.5 * nRotors) * 2)#(4*R))
            self.cpInf = (self.cp0 + (nFreq + 1.5 * nRotors) * (2*R))
        elif linear in ["no", "N", "n", "NO", "No", "nO"]:
            self.cp0   = 4.0 * R
            #self.cpInf = (self.cp0 + (nFreq + 0.5 * nRotors) * 2)#(4*R))
            self.cpInf = (self.cp0 + (nFreq + 0.5 * nRotors) * (2*R))
        else:
            self.cp0   = 4.0 * R
            #self.cpInf = (self.cp0 + (nFreq + 0.5 * nRotors) * 2)#(2*R))
            self.cpInf = (self.cp0 + (nFreq + 0.5 * nRotors) * (2*R))
        #self.cpInf = 25#(3*8 - 2)*R
        #print(self.cpInf)
        #self.cpInf2 = (self.cp0 + (nFreq + 0.5 * nRotors) * R) #* 0.9)
        # What remains is to fit the polynomial coefficients (a0, a1, a2, a3)
        # This can be done directly - no iteration required
        y = Tlist / (Tlist + B)
        A = np.zeros((len(Cplist), 4), np.float64)
        for j in range(4):
            A[:, j] = (y * y * y - y * y) * y ** j
        Cplist = np.array([x for x in Cplist], np.float64)
        b2 = ((Cplist - self.cp0) / (self.cpInf - self.cp0) - y * y)
        x, residues, rank, s = np.linalg.lstsq(A, b2, rcond=None)
        self.B  = float(B)
        self.a0 = float(x[0])
        self.a1 = float(x[1])
        self.a2 = float(x[2])
        self.a3 = float(x[3])
        self.H0 = 0.0
        self.S0 = 0.0
        self.H0 = H298 - self.getEnthalpy(298.15)
        self.S0 = S298 - self.getEntropy(298.15)

        return self

