#!/usr/bin/env python
# coding: utf-8
# **********************************************
# Python modules to work with 
# **********************************************
import numpy as np


class Thermo:
    """
    This class calculates the thermochemistry properties
    requested by user in input file such as Hf, S and Cp from 300 - 1500 K.
    Entropy correction:  -R*ln(σ) where σ = (σext*σint)/nopt
    where σext is the external symmetry number, σint is the internal rotational symmetry number,
    and nopt is the number of optical isomers. Therefore, R*ln(σ) was added to S298 in eq 6 before
    fitting the group additivity values.
    """

    def __init__(self, dfobj, speciesname, formula, numberofgroups, groupid, quantity, symmetrynumber, numberofrotors):
        self.dfObj = dfobj
        self.SpeciesName = speciesname
        self.Formula = formula
        self.NumberOfGroups = numberofgroups
        self.GroupNames = groupid
        self.Quantity = quantity
        self.SymmetryNumber = symmetrynumber
        self.NumberOfRotors = numberofrotors

    def thermo_props(self):
        hfs            = []
        GroupEnthalpos = []
        ValEnthalpos   = []
        s              = []
        ValEntropos    = []
        cp300          = []
        Valcp300       = []
        cp400          = []
        Valcp400       = []
        cp500          = []
        Valcp500       = []
        cp600          = []
        Valcp600       = []
        cp800          = []
        Valcp800       = []
        cp1000         = []
        Valcp1000      = []
        cp1500         = []
        Valcp1500      = []
        for u in range(len(self.GroupNames)):
            for k in range(len(self.dfObj)):
                if self.GroupNames[u] == self.dfObj.iloc[k, 0]:
                    qty       = float(self.Quantity[u])
                    henthalpy = float(self.dfObj.iloc[k, 1])
                    GroupEnthalpos.append(self.GroupNames[u])
                    ValEnthalpos.append(henthalpy)
                    Sentropy  = float(self.dfObj.iloc[k, 2])
                    ValEntropos.append(Sentropy)
                    CP300     = float(self.dfObj.iloc[k, 3])
                    Valcp300.append(CP300)
                    CP400     = float(self.dfObj.iloc[k, 4])
                    Valcp400.append(CP400)
                    CP500     = float(self.dfObj.iloc[k, 5])
                    Valcp500.append(CP500)
                    CP600     = float(self.dfObj.iloc[k, 6])
                    Valcp600.append(CP600)
                    CP800     = float(self.dfObj.iloc[k, 7])
                    Valcp800.append(CP800)
                    CP1000    = float(self.dfObj.iloc[k, 8])
                    Valcp1000.append(CP1000)
                    CP1500    = float(self.dfObj.iloc[k, 9])
                    Valcp1500.append(CP1500)
                    hfs.append(str(henthalpy * qty))
                    s.append(str(Sentropy * qty))
                    cp300.append(str(CP300 * qty))
                    cp400.append(str(CP400 * qty))
                    cp500.append(str(CP500 * qty))
                    cp600.append(str(CP600 * qty))
                    cp800.append(str(CP800 * qty))
                    cp1000.append(str(CP1000 * qty))
                    cp1500.append(str(CP1500 * qty))
                else:
                    pass
        enthalpy = round(sum([float(i) for i in hfs]), 2)
        entropy  = round(sum([float(i) for i in s]) - 1.9872 * np.log(self.SymmetryNumber), 2)
        # Radicals:
        # Entropy  = round(sum([float(i) for i in S])-(1.98720425864083)*np.log(self.SymmetryNumber),2);
        cpp300   = round(sum([float(i) for i in cp300]), 2)
        cpp400   = round(sum([float(i) for i in cp400]), 2)
        cpp500   = round(sum([float(i) for i in cp500]), 2)
        cpp600   = round(sum([float(i) for i in cp600]), 2)
        cpp800   = round(sum([float(i) for i in cp800]), 2)
        cpp1000  = round(sum([float(i) for i in cp1000]), 2)
        cpp1500  = round(sum([float(i) for i in cp1500]), 2)
        if  cpp400 < cpp300 and cpp500 < cpp600 or cpp400 == 0.0 and cpp500 < cpp600:
            cpp400 = round((cpp300 + cpp500)/2, 2)
        else:
            pass
        if  cpp500 < cpp400 and cpp600 < cpp800 or cpp500 == 0.0 and cpp600 < cpp800:
            cpp500 = round((cpp400 + cpp600)/2, 2)
        else:
            pass
        if  cpp600 < cpp500 and cpp800 < cpp1000 or cpp600 == 0.0 and cpp800 < cpp1000:
            cpp600 = round((cpp500 + cpp800)/2, 2)
        else:
            pass
        if  cpp800 < cpp600 and cpp1000 < cpp1500 or cpp800 == 0.0 and cpp1000 < cpp1500:
            cpp800 = round((cpp600 + cpp1000)/2, 2)
        else:
            pass
        if  cpp1000 < cpp800 or cpp1000 == 0.0:
            cpp1000 = round((cpp800 + cpp1500)/2, 2)
        else:
            pass
        if  cpp1500 < cpp1000 or cpp1500 == 0.0:
            cpp1500 = round((cpp800 + cpp1000)*0.6, 2)
        else:
            pass        
        ZIPTest = zip(GroupEnthalpos, ValEnthalpos, ValEntropos, Valcp300, Valcp400, Valcp500, Valcp600, Valcp800,
                      Valcp1000, Valcp1500)
        return enthalpy, entropy, ZIPTest, cpp300, cpp400, cpp500, cpp600, cpp800, cpp1000, cpp1500

    def thermo_props_rads(self):
        hfs            = []
        GroupEnthalpos = []
        ValEnthalpos   = []
        s              = []
        ValEntropos    = []
        cp300          = []
        Valcp300       = []
        cp400          = []
        Valcp400       = []
        cp500          = []
        Valcp500       = []
        cp600          = []
        Valcp600       = []
        cp800          = []
        Valcp800       = []
        cp1000         = []
        Valcp1000      = []
        cp1500         = []
        Valcp1500      = []
        for u in range(len(self.GroupNames)):
            for k in range(len(self.dfObj)):
                if self.GroupNames[u] == self.dfObj.iloc[k, 0]:
                    qty       = float(self.Quantity[u])
                    henthalpy = float(self.dfObj.iloc[k, 1])
                    GroupEnthalpos.append(self.GroupNames[u])
                    ValEnthalpos.append(henthalpy)
                    Sentropy  = float(self.dfObj.iloc[k, 2])
                    ValEntropos.append(Sentropy)
                    CP300     = float(self.dfObj.iloc[k, 3])
                    Valcp300.append(CP300)
                    CP400     = float(self.dfObj.iloc[k, 4])
                    Valcp400.append(CP400)
                    CP500     = float(self.dfObj.iloc[k, 5])
                    Valcp500.append(CP500)
                    CP600     = float(self.dfObj.iloc[k, 6])
                    Valcp600.append(CP600)
                    CP800     = float(self.dfObj.iloc[k, 7])
                    Valcp800.append(CP800)
                    CP1000    = float(self.dfObj.iloc[k, 8])
                    Valcp1000.append(CP1000)
                    CP1500    = float(self.dfObj.iloc[k, 9])
                    Valcp1500.append(CP1500)
                    hfs.append(str(henthalpy * qty))
                    s.append(str(Sentropy * qty))
                    cp300.append(str(CP300 * qty))
                    cp400.append(str(CP400 * qty))
                    cp500.append(str(CP500 * qty))
                    cp600.append(str(CP600 * qty))
                    cp800.append(str(CP800 * qty))
                    cp1000.append(str(CP1000 * qty))
                    cp1500.append(str(CP1500 * qty))
                else:
                    pass
        enthalpy = round(sum([float(i) for i in hfs]) - 52.103, 2)
        # Radicals:
        # Entropy  = round(sum([float(i) for i in S])-(1.98720425864083)*np.log(self.SymmetryNumber),2);
        entropy  = round(sum([float(i) for i in s]) - 1.9872 * np.log(self.SymmetryNumber), 2)
        cpp300   = round(sum([float(i) for i in cp300]), 2)
        cpp400   = round(sum([float(i) for i in cp400]), 2)
        cpp500   = round(sum([float(i) for i in cp500]), 2)
        cpp600   = round(sum([float(i) for i in cp600]), 2)
        cpp800   = round(sum([float(i) for i in cp800]), 2)
        cpp1000  = round(sum([float(i) for i in cp1000]), 2)
        cpp1500  = round(sum([float(i) for i in cp1500]), 2)
        if  cpp400 < cpp300 and cpp500 < cpp600 or cpp400 == 0.0 and cpp500 < cpp600:
            cpp400 = round((cpp300 + cpp500)/2, 2)
        else:
            pass
        if  cpp500 < cpp400 and cpp600 < cpp800 or cpp500 == 0.0 and cpp600 < cpp800:
            cpp500 = round((cpp400 + cpp600)/2, 2)
        else:
            pass
        if  cpp600 < cpp500 and cpp800 < cpp1000 or cpp600 == 0.0 and cpp800 < cpp1000:
            cpp600 = round((cpp500 + cpp800)/2, 2)
        else:
            pass
        if  cpp800 < cpp600 and cpp1000 < cpp1500 or cpp800 == 0.0 and cpp1000 < cpp1500:
            cpp800 = round((cpp600 + cpp1000)/2, 2)
        else:
            pass
        if  cpp1000 < cpp800 or cpp1000 == 0.0:
            cpp1000 = round((cpp800 + cpp1500)/2, 2)
        else:
            pass
        if  cpp1500 < cpp1000 or cpp1500 == 0.0:
            cpp1500 = round((cpp800 + cpp1000)*0.6, 2)
        else:
            pass                 
        ZIPTest = zip(GroupEnthalpos,ValEnthalpos, ValEntropos, Valcp300, Valcp400, Valcp500, Valcp600, Valcp800,
                      Valcp1000, Valcp1500)
        return enthalpy, entropy, ZIPTest, cpp300, cpp400, cpp500, cpp600, cpp800, cpp1000, cpp1500
