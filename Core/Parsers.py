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
from sympy import *
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from Core.Inputs import Inputs
from datetime import date


class Parser:
    """
        This class collect all GAVs definitions from *.grp files to further
        calculate the thermochemistry properties. Benson group additivity 'GAVs'.
    """

    def __init__(self, groupfiles, h=0.0, s=0.0, cp300=0.0, cp400=0.0, cp500=0.0, cp600=0.0, cp800=0.0, cp1000=0.0, cp1500=0.0, temps=0.0):
        self.GroupFiles = groupfiles
        self.H      = h
        self.S      = s
        self.Cp300  = cp300
        self.Cp400  = cp400
        self.Cp500  = cp500
        self.Cp600  = cp600
        self.Cp800  = cp800
        self.Cp1000 = cp1000
        self.Cp1500 = cp1500
        self.Temps  = temps

    #@property
    def data2(self):
        dfobj = pd.DataFrame(columns=['GAVs', 'Hf', 'S', 'Cp: 300', '400', '500', '600', '800', '1000', '1500'])
        PythonVersion = sys.version
        #print(PythonVersion)
        for p in range(len(self.GroupFiles)):
            if PythonVersion > "3.7":
                datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", encoding='cp1252', engine='python', comment="!", skiprows=[0])
                #datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", encoding='cp1252', comment="!")
            else:
                datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", engine="python", comment="!", skiprows=[0])
            for m in range(len(datafile)):
                #if "," in datafile.iloc[m, 0]:
                #    dfobj = dfobj.append({'GAVs': datafile.iloc[m, 0].split(",")[0],
                #                          'Hf': datafile.iloc[m, 0].split()[1],
                #                          'S': datafile.iloc[m, 0].split()[2],
                #                          'Cp: 300': datafile.iloc[m, 0].split()[3],
                #                          '400': datafile.iloc[m, 0].split()[4],
                #                          '500': datafile.iloc[m, 0].split()[5],
                #                          '600': datafile.iloc[m, 0].split()[6],
                #                          '800': datafile.iloc[m, 0].split()[7],
                #                          '1000': datafile.iloc[m, 0].split()[8],
                #                          '1500': datafile.iloc[m, 0].split()[9]}, ignore_index=True)
                #else:
                #pass
                ###dfobj = dfobj.append({'GAVs': datafile.iloc[m, 0].split()[0],
                ###                      'Hf': datafile.iloc[m, 0].split()[1],
                ###                      'S': datafile.iloc[m, 0].split()[2],
                ###                      'Cp: 300': datafile.iloc[m, 0].split()[3],
                ###                      '400': datafile.iloc[m, 0].split()[4],
                ###                      '500': datafile.iloc[m, 0].split()[5],
                ###                      '600': datafile.iloc[m, 0].split()[6],
                ###                      '800': datafile.iloc[m, 0].split()[7],
                ###                      '1000': datafile.iloc[m, 0].split()[8],
                ###                      '1500': datafile.iloc[m, 0].split()[9]}, ignore_index=True)                    
                Data1Frame = pd.DataFrame({'GAVs': datafile.iloc[m, 0].split()[0],
                                      'Hf': datafile.iloc[m, 0].split()[1],
                                      'S': datafile.iloc[m, 0].split()[2],
                                      'Cp: 300': datafile.iloc[m, 0].split()[3],
                                      '400': datafile.iloc[m, 0].split()[4],
                                      '500': datafile.iloc[m, 0].split()[5],
                                      '600': datafile.iloc[m, 0].split()[6],
                                      '800': datafile.iloc[m, 0].split()[7],
                                      '1000': datafile.iloc[m, 0].split()[8],
                                      '1500': datafile.iloc[m, 0].split()[9]}, index=[0])
                dfobj = pd.concat([dfobj, Data1Frame], ignore_index=True)

        return dfobj
    #@property
    def data3(self):
        PythonVersion = sys.version
        #print(PythonVersion)
        for p in range(len(self.GroupFiles)):
            if PythonVersion > "3.7":
                datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", encoding='cp1252', engine='python', comment="!", skiprows=[0])
                #datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", encoding='cp1252', comment="!")
            else:
                datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", engine="python", comment="!", skiprows=[0])                   
        DF = datafile.replace(r'\t','', regex=True) 
        return DF
    #@property
    def data(self):
        dfobj = pd.DataFrame(columns=['GAVs', 'Hf', 'S', 'Cp: 300', '400', '500', '600', '800', '1000', '1500'])
        PythonVersion = sys.version
        #print(PythonVersion)
        for p in range(len(self.GroupFiles)):
            if PythonVersion > "3.7":
                datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", encoding='cp1252', engine='python', comment="!")
                #datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", encoding='cp1252', comment="!")
            else:
                datafile = pd.read_csv(self.GroupFiles[p], sep="s+\t", engine="python", comment="!")
            for m in range(len(datafile)):
                if "," in datafile.iloc[m, 0]:
                    ###dfobj = dfobj.append({'GAVs': datafile.iloc[m, 0].split(",")[0],
                    ###                      'Hf': datafile.iloc[m, 0].split()[1],
                    ###                      'S': datafile.iloc[m, 0].split()[2],
                    ###                      'Cp: 300': datafile.iloc[m, 0].split()[3],
                    ###                      '400': datafile.iloc[m, 0].split()[4],
                    ###                      '500': datafile.iloc[m, 0].split()[5],
                    ###                      '600': datafile.iloc[m, 0].split()[6],
                    ###                      '800': datafile.iloc[m, 0].split()[7],
                    ###                      '1000': datafile.iloc[m, 0].split()[8],
                    ###                      '1500': datafile.iloc[m, 0].split()[9]}, ignore_index=True)
                    Data1Frame = pd.DataFrame({'GAVs': datafile.iloc[m, 0].split(",")[0],
                                          'Hf': datafile.iloc[m, 0].split()[1],
                                          'S': datafile.iloc[m, 0].split()[2],
                                          'Cp: 300': datafile.iloc[m, 0].split()[3],
                                          '400': datafile.iloc[m, 0].split()[4],
                                          '500': datafile.iloc[m, 0].split()[5],
                                          '600': datafile.iloc[m, 0].split()[6],
                                          '800': datafile.iloc[m, 0].split()[7],
                                          '1000': datafile.iloc[m, 0].split()[8],
                                          '1500': datafile.iloc[m, 0].split()[9]}, index=[0])
                    dfobj = pd.concat([dfobj, Data1Frame], ignore_index=True)

                else:
                    pass #Data1Frame = pd.DataFrame({'GAVs': datafile.iloc[m, 0].split()[0],
                    #                      'Hf': datafile.iloc[m, 0].split()[1],
                    #                      'S': datafile.iloc[m, 0].split()[2],
                    #                      'Cp: 300': datafile.iloc[m, 0].split()[3],
                    #                      '400': datafile.iloc[m, 0].split()[4],
                    #                      '500': datafile.iloc[m, 0].split()[5],
                    #                      '600': datafile.iloc[m, 0].split()[6],
                    #                      '800': datafile.iloc[m, 0].split()[7],
                    #                      '1000': datafile.iloc[m, 0].split()[8],
                    #                      '1500': datafile.iloc[m, 0].split()[9]}, index=[0])
                    #dfobj = pd.concat([dfobj, Data1Frame], ignore_index=True)

        return dfobj

    def fit_termo(self):
        """
        This class calculates the 1st set of coefficients base on the thermochemistry properties
        requested by user in the input file such as Hf(298), S(298) and Cp from 300 - 800 K.
        Ro = 1.987×10^−3 kcal⋅K^−1⋅mol^−1
        Cp(T)/Ro = a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4; solving for 300 - 800 K

        Way to fit the best curve to data:
        popt, pcov  = curve_fit(funcCP, Tlimit, NewDataCP)
        b1 = (popt[0])
        b2 = (popt[1])
        b3 = (popt[2])
        b4 = (popt[3])
        b5 = (popt[4])

        """
        ###CPs1 = []
        ###CPs2 = []
        ###CPs2Fit = []
        #### BEGING definition

        r_kcal_k_mol = 0.0019872
        r_cal_k_mol = r_kcal_k_mol * 1000
        # rj_k_mol = 8.3145
        """ 
        Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
        """
        a1, a2, a3, a4, a5 = symbols('a1, a2, a3, a4, a5', real=True)
        # ...........................................................
        CPSLIST = []
        if self.Cp300 != 0:
           CPSLIST.append(self.Cp300/r_cal_k_mol)
        else:
           CPSLIST.append(00000)
        if self.Cp400 != 0:
           CPSLIST.append(self.Cp400/r_cal_k_mol)
        else:
           CPSLIST.append(11111)
        if self.Cp500 != 0:
           CPSLIST.append(self.Cp500/r_cal_k_mol)
        else:
           CPSLIST.append(22222)
        if self.Cp600 != 0:
           CPSLIST.append(self.Cp600/r_cal_k_mol)
        else:
           CPSLIST.append(33333)
        if self.Cp800 != 0:
           CPSLIST.append(self.Cp800/r_cal_k_mol)
        else:
           CPSLIST.append(44444)
        if self.Cp1000 != 0:
           CPSLIST.append(self.Cp1000/r_cal_k_mol)
        else:
           CPSLIST.append(55555)
        if self.Cp1500 != 0:
           CPSLIST.append(self.Cp1500/r_cal_k_mol)
        else:
           CPSLIST.append(66666)
        # deleting items not suitable for fitting from CPLIST and Temps list
        if 00000 in CPSLIST:
           CPSLIST.pop(0)
           self.Temps.pop(0)
        else:
           pass
        if 11111 in CPSLIST:
           CPSLIST.pop(1)
           self.Temps.pop(1)
        else:
           pass
        if 22222 in CPSLIST:
           CPSLIST.pop(2)
           self.Temps.pop(2)
        else:
           pass
        if 33333 in CPSLIST:
           CPSLIST.pop(3)
           self.Temps.pop(3)
        else:
           pass
        if 44444 in CPSLIST:
           CPSLIST.pop(4)
           self.Temps.pop(4)
        else:
           pass
        if 55555 in CPSLIST:
           CPSLIST.pop(5)
           self.Temps.pop(5)
        else:
           pass
        if 66666 in CPSLIST:
           CPSLIST.pop(6)
           self.Temps.pop(6)
        else:
           pass
        # Now the fitting can be done, at least the 1st fitting *****  FIRST FIT   ******
        #print("\t Lista de Cps: ",   CPSLIST)
        #print("\t Lista de Temps: ", self.Temps)
        def funcCP(T, a1, a2, a3, a4, a5):
            R = 1.98#72  # in cal units
            # Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
            # CP/R here produce by code
            return (a1 + a2 * (T) + a3 * (T) ** 2 + a4 * (T) ** 3 + a5 * (T) ** 4)
        popt, pcov  = curve_fit(funcCP, self.Temps, CPSLIST)
        f1 = (popt[0])
        f2 = (popt[1])
        f3 = (popt[2])
        f4 = (popt[3])
        f5 = (popt[4])
        Tempos = np.zeros(1701)  # np.zeros(471)
        Ti = 300
        DeltaT = 1
        counter = 0
        for u in range(len(Tempos)):
            Tempos[counter] = (Ti + DeltaT * counter)
            counter += 1
        TLT  = Tempos
        CPLT = []
        for h in range(len(TLT)):
            CPLT.append(funcCP(TLT[h], f1, f2, f3, f4, f5))
        #print(TLT)        
        # 2nd fitting *****  SECOND FIT   ******
        popt, pcov  = curve_fit(funcCP, TLT, CPLT)
        a1 = (popt[0])
        a2 = (popt[1])
        a3 = (popt[2])
        a4 = (popt[3])
        a5 = (popt[4])      
        """
        # H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
        H(T)/Ro = a1*T + (a2/2)*T**2 + (a3/3)*T**3 + (a4/4)*T**4 + (a5/5)*T**5 + a6; solving with T = 298 K 
        """
        a6 = (self.H / r_kcal_k_mol) - a1 * 298 - (a2 / 2) * 298 ** 2 - (a3 / 3) * 298 ** 3 - (a4 / 4) * (
            298) ** 4 - (a5 / 5) * 298 ** 5
        """
        S(T)/Ro = a1*ln(T) + a2*T + (a3/2)*T**2 + (a4/3)*T**3 + (a5/4)*T**4 + a7;
        Natural logarithm calculated by 'np.log()'
        """
        a7 = (self.S / r_cal_k_mol) - a1 * np.log(298) - a2 * 298 - (a3 / 2) * 298 ** 2 - (a4 / 3) * 298 ** 3 - (
                a5 / 4) * 298 ** 4
        return (a1, a2, a3, a4, a5, a6, a7)
    def guess_Cps(self):
        """
        This definition calculates the 1st set of coefficients based on the thermochemistry properties;
        However, this guesser can only complete the gaps of a series of T = [300, 400, 500, 600, 800, 1000, 1500]
        It requires at least 4 entries in order to calculate any other CP(T), if the file has gaps wihtout numbers,
        it should be filled with 0.00 values.
        """
        def func(T, a1, a2, a3, a4):
            return a1 + a2 * (T) + a3 * (T) ** 2 + a4 * (T) ** 3
        def Zeroed(CPT, number):
            if CPT == 0:
                ZEROS.append(number)
                #print("\t Its Zero: " + number)
            else:
                NOTZEROS.append(number)

        # ============================
        ZEROS    = []
        NOTZEROS = []
        # ============================
        Zeroed(self.Cp300,   "Cp300")
        Zeroed(self.Cp400,   "Cp400")
        Zeroed(self.Cp500,   "Cp500")
        Zeroed(self.Cp600,   "Cp600")
        Zeroed(self.Cp800,   "Cp800")
        Zeroed(self.Cp1000, "Cp1000")
        Zeroed(self.Cp1500, "Cp1500")
        # ============================
        #print("\t Zeroes:\t", ZEROS)
        #print("\t Notzeroes:\t", NOTZEROS)
        CPNewData = []
        Temps     = []
        for k in NOTZEROS:             
            if   k == "Cp300":             
                   CPNewData.append(self.Cp300)
                   Temps.append(300)
            elif k == "Cp400":             
                   CPNewData.append(self.Cp400)
                   Temps.append(400)
            elif k == "Cp500":             
                   CPNewData.append(self.Cp500)
                   Temps.append(500)
            elif k == "Cp600":             
                   CPNewData.append(self.Cp600)
                   Temps.append(600)
            elif k == "Cp800":             
                   CPNewData.append(self.Cp800)
                   Temps.append(800)
            elif k == "Cp1000":             
                   CPNewData.append(self.Cp1000)
                   Temps.append(1000)
            elif k == "Cp1500":             
                   CPNewData.append(self.Cp1500)
                   Temps.append(1500)

        #CPNewData = [self.Cp300, self.Cp500, self.Cp1000, self.Cp1500]
        #print(CPNewData)
        #print(Temps)
        #Temps = [300, 500, 1000, 1500]
        popt, pcov = curve_fit(func, Temps, CPNewData)
        # print("\t ",popt)
        b1 = (popt[0])
        b2 = (popt[1])
        b3 = (popt[2])
        b4 = (popt[3])
        #print(ZEROS)
        return b1, b2, b3, b4

class CoeffReader:
    """
        This class calculates the thermochemistry properties Hf(T), S(T) and CP(T)
        from 1st and 2nd set of coefficients base on the *.therm or *.dat format files.
    """

    def __init__(self, InputFile, MaxT):
        self.InputFile = InputFile
        self.MaxT      = MaxT

    def reading_coeffs(self):
        """
        This definition plots every specie in *.dat files in individual graphs
        providing the Low Temperature (LT) NASA polynomials and the High Temperature
        (HT) NASA polynomials along with the breaking point and the Absolute Percentage
        Relative Error (APRE) between these two set of polynomial at 1000 K (join point). 
        """
        #print("\t Under construction")
        print("\t Plotting: ")
        print("\t")
        print("\t " + str(self.InputFile))
        print("\t")
        CWD = os.getcwd()
        PandasFile = pd.read_table(CWD + "\\OutputsDir\\" + self.InputFile, names=['col'], encoding='cp1252',engine='python', comment="!")
        #PandasFile = pd.read_csv(CWD + "\\OutputsDir\\" + self.InputFile, names=['col'])
        #PandasFile = pd.read_table(CWD + "\\OutputsDir\\" + self.InputFile)#, header=0)
        #PanFile = pd.DataFrame(PandasFile)
        SpeciesThermos = []
        for k in range(len(PandasFile)):
            FirstData = str(PandasFile.iloc[k,0]).split()[0]
            SecondData = str(PandasFile.iloc[k,0]).split()
            #print(FirstData)
            #print(SecondData)
            if   "C" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "H" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "HE" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "AR" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "H" in FirstData and "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "H" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "O" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData and "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData and "O" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "0G" in SecondData:
                 SpeciesThermos.append(FirstData)
            elif "G" in SecondData:
                 SpeciesThermos.append(FirstData)
            elif "1 " in SecondData:
                 SpeciesThermos.append(FirstData)
            else:
                 #SpeciesThermos.append(FirstData)
                 pass #print("\t Species not found!")
        
        print("\t The next is a list of species found in the file")
        SpeciesThermos = list(set(SpeciesThermos))
        if   "THERMO" in SpeciesThermos: 
             SpeciesThermos.remove("THERMO")
        else:
             pass
        if   "END" in SpeciesThermos: 
             SpeciesThermos.remove("END")
        else:
             pass
        #print(SpeciesThermos)
        for j in SpeciesThermos:
            print("\t "+j, end=' ')
        print("\t ")
        #HowMany2Plot = input("\t How many species do you want to plot? (just hit Enter if you want to plot all of them)\n\t ") or 0
        HowMany2Plot = input("\t How many species do you want to plot? (Type 'all' to plot them all)\n\t ") or 0
        Names2Plot   = []
        try:
            HowMany2Plot = int(HowMany2Plot)
            if  HowMany2Plot == 0:
                print("\t 0 species provided, no plots generated!.")
                print("\t Please type a number or type 'all' to plot them all.")
                print("\t If you press 'Enter' again code will end.")
                TypeO = input("\t ") or "end"
                if   TypeO == "end":
                      pass
                elif TypeO in [0]:
                      print("\t 0 species provided, no plots generated!.")
                      print("\t Ending code.")
                      pass
                elif TypeO in ["ALL", "all", "ALl", "All", "aLL", "alL"]:
                      TypeO = len(SpeciesThermos)
                      for t in range(TypeO):
                          Names2Plot.append(SpeciesThermos[t])
                else:
                      for t in range(TypeO):
                          Names2Plot.append(input("\t Species #" + str(t+1) +" name:"))
                 #pass
            else:
                 for t in range(HowMany2Plot):
                    Names2Plot.append(input("\t Species #" + str(t+1) +" name:"))
        except:
            if HowMany2Plot in ["ALL", "all", "ALl", "All", "aLL", "alL"]:
                 HowMany2Plot = len(SpeciesThermos)
                 for t in range(HowMany2Plot):
                    Names2Plot.append(SpeciesThermos[t])
            else:
                 print("\t Wrong answer, try again... ")
                 pass           
        #print("\t By default this code save figures in a poor quality")
        #print("\t Press 'Enter' to continue with the default settings\n\t or provide the number of dpi (dots per inch) desired")
        #print("\t For example: 800")
        #print("\t This is high quality graph, but code will take longer in save them.")
        #print("\t d.p.i. : dots per inch")
        # ---------------------------------
        #print("PandasFile: ", PandasFile)
        #print("--"*50)
        #print("df: ", df)
        #print("--"*50)
        #print("HowMany2Plot: ", HowMany2Plot)
        #print("Names2Plot: ", Names2Plot)
        #print("Names2Plot len: ", len(Names2Plot))
        df        = pd.DataFrame(PandasFile)
        try:
            df["col"] = df["col"].str.replace("(","_", regex=True)
            df["col"] = df["col"].str.replace(")","_", regex=True)
        except:
            pass
        #print(df)
        Xpecie      = []
        XtempBP     = []
        Abs_ErrorCP = []
        Abs_ErrorH  = []
        Abs_ErrorS  = []
        print("\t Total number of species in Thermo: ", HowMany2Plot)#, len(HowMany2Plot))
        #for y in range(len(Names2Plot)):
        #    print("#",y+1)
        # ---------------------------------
        for u in range(len(Names2Plot)):
            print("--"*50)
            print("Names2Plot: ", Names2Plot)
            print("Names2Plot len: ", len(Names2Plot))
            print("\t" + "--"*10)
            print("\t Plotting NASA polynomials for: " + Names2Plot[u])
            print("\t" + ".."*10)
            # DropingStuff[DropingStuff.col.str.contains("IDT") == False]
            #print(PandasFile[PandasFile.col.str.contains(Names2Plot[u]) == True])
            # ---------------------------------------
            VarSpec = (Names2Plot[u] + " ")
            try:
                SpecieIndex = int((str(df[df.col.str.contains(VarSpec) == True]).split()[1]))
            except:
                VarSpec     = str(VarSpec.replace("(","_")).replace(")","_")
                SpecieIndex = int((str(df[df.col.str.contains(VarSpec) == True]).split()[1]))
            # ---------------------------------------
            #SpecieIndex = int((str(PandasFile[PandasFile.col.str.contains(Names2Plot[u]+" ") == True]).split()[1]))
            #print("\t Loop #" + str(u))
            #print("\t Index specie #" + str(SpecieIndex))
            print("\t Ploynomials extracted:")
            BPRowLr      = []
            BPRow        = (PandasFile.iloc[SpecieIndex,:]).tolist()
            BPRow2       = list(BPRow[0])
            for w in range(len(BPRow2)):
                if w >= 66 and w < 77:
                   BPRowLr.append(BPRow2[w])
            #print(BPRow)
            #print(BPRow2)
            #print(BPRowLr)
            FirstRow     = (PandasFile.iloc[SpecieIndex+1,:]).tolist()
            SecondRow    = (PandasFile.iloc[SpecieIndex+2,:]).tolist()
            ThirdRow     = (PandasFile.iloc[SpecieIndex+3,:]).tolist()
            PreFirstRow  = list(FirstRow[0])
            PreSecondRow = list(SecondRow[0])
            PreThirdRow  = list(ThirdRow[0])
            # Extracting polynomial coefficients from b1 to b5
            b1 = []; b2 = []; b3 = []; b4 = []; b5 = []; 
            for p in range(len(PreFirstRow)):
                if   p < 15:
                     b1.append(PreFirstRow[p])
                elif p >= 15 and p < 30:
                     b2.append(PreFirstRow[p])
                elif p >= 30 and p < 45:
                     b3.append(PreFirstRow[p])
                elif p >= 45 and p < 60:
                     b4.append(PreFirstRow[p])
                elif p >= 60 and p < 75:
                     b5.append(PreFirstRow[p])
            # Extracting polynomial coefficients from b6 to a3
            b6 = []; b7 = []; a1 = []; a2 = []; a3 = []; 
            for p in range(len(PreSecondRow)):
                if   p < 15:
                     b6.append(PreSecondRow[p])
                elif p >= 15 and p < 30:
                     b7.append(PreSecondRow[p])
                elif p >= 30 and p < 45:
                     a1.append(PreSecondRow[p])
                elif p >= 45 and p < 60:
                     a2.append(PreSecondRow[p])
                elif p >= 60 and p < 75:
                     a3.append(PreSecondRow[p])
            # Extracting polynomial coefficients from a4 to a7
            a4 = []; a5 = []; a6 = []; a7 = [];
            for p in range(len(PreThirdRow)):
                if   p < 15:
                     a4.append(PreThirdRow[p])
                elif p >= 15 and p < 30:
                     a5.append(PreThirdRow[p])
                elif p >= 30 and p < 45:
                     a6.append(PreThirdRow[p])
                elif p >= 45 and p < 60:
                     a7.append(PreThirdRow[p])
            b1      = "".join(b1); b2 = "".join(b2); b3 = "".join(b3); b4 = "".join(b4); b5 = "".join(b5); b6 = "".join(b6); b7 = "".join(b7); 
            a1      = "".join(a1); a2 = "".join(a2); a3 = "".join(a3); a4 = "".join(a4); a5 = "".join(a5); a6 = "".join(a6); a7 = "".join(a7); 
            BPRowL  = round(float(("".join(BPRowLr)).strip()))
            BPRowLZ = float(("".join(BPRowLr)).strip()); 
            if BPRowL == 0 and BPRowLZ == 0.0:
               BPRowLr      = []
               BPRow        = (PandasFile.iloc[SpecieIndex,:]).tolist()
               BPRow2       = list(BPRow[0])
               for w in range(len(BPRow2)):
                   if w >= 65 and w < 76:
                      BPRowLr.append(BPRow2[w])
               BPRowL  = round(float(("".join(BPRowLr)).strip()))
               BPRowLZ = float(("".join(BPRowLr)).strip()); 
            else:
               pass

            print("\t Breaking Point (BP) at : " + str(BPRowLZ))
            print("\t a1 : " + str(a1))
            print("\t a2 : " + str(a2))
            print("\t a3 : " + str(a3))
            print("\t a4 : " + str(a4))
            print("\t a5 : " + str(a5))
            print("\t a6 : " + str(a6))
            print("\t a7 : " + str(a7))
            print("\t a8 : " + str(b1))
            print("\t a9 : " + str(b2))
            print("\t a10: " + str(b3))
            print("\t a11: " + str(b4))
            print("\t a12: " + str(b5))
            print("\t a13: " + str(b6))
            print("\t a14: " + str(b7))
             
            # Converting data from NASA to thermochemistry: Cps/R, H and S to plot
            def funcCP(T, m1, m2, m3, m4, m5):
                R = 1.9872  # in cal units
                # Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
                # CP/R here produce by code
                return ((m1) + (m2 * (T)) + (m3 * (T) ** 2) + (m4 * (T) ** 3) + (m5 * (T) ** 4))*R

            def funcH(T, m1, m2, m3, m4, m5, m6):
                T2 = T * T
                T4 = T2 * T2
                R = 0.0019872  # in kcal units
                # H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
                # H here produce by code
                return (m1 + m2 * T/2 + m3 * T2/3 + m4 * T2 * T/4 + m5 * T4/5 + m6 / T) * R * T

            def funcS(T, m1, m2, m3, m4, m5, m7):
                import math
                R = 1.9872  # in cal units
                T2 = T * T
                T4 = T2 * T2
                # S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
                # S here produce by code
                return (m1 * math.log(T) + m2 * T + m3 * T2/2 + m4 * T2 * T/3 + m5 * T4/4 + m7) * R 
            CPs1 = []; Hs1 = []; Ss1 = []; CPs2 = []; Hs2 = []; Ss2 = [];
            #Temperatures2Plot1 = np.arange(300, BPRowL+1, 1).tolist()
            #BPRowL  = 1000
            #BPRowLZ = 1000
            #print(BPRowL)
            Temperatures2Plot1 = np.arange(300, BPRowL+1, 1).tolist()
            Temperatures2Plot2 = np.arange(BPRowL, self.MaxT + 100, 1).tolist()
            CWD2 = os.getcwd()
            OutputPlotPath = CWD2 + "\\OutputsDir\\"
            if not os.path.exists(OutputPlotPath + "\\Plots\\Single\\" + self.InputFile):
                   os.makedirs(OutputPlotPath + "\\Plots\\Single\\" + self.InputFile)
            for d1 in (Temperatures2Plot1):
                CPs1.append(round(funcCP(d1, float(a1), float(a2), float(a3), float(a4), float(a5)),5))
                Hs1.append(round(funcH(d1, float(a1), float(a2), float(a3), float(a4), float(a5), float(a6)),5))
                Ss1.append(round(funcS(d1, float(a1), float(a2), float(a3), float(a4), float(a5), float(a7)),5))
            for d2 in (Temperatures2Plot2):
                CPs2.append(round(funcCP(d2, float(b1), float(b2), float(b3), float(b4), float(b5)),5))
                Hs2.append(round(funcH(d2, float(b1), float(b2), float(b3), float(b4), float(b5), float(b6)),5))
                Ss2.append(round(funcS(d2, float(b1), float(b2), float(b3), float(b4), float(b5), float(b7)),5))
            # printing on screen the Breaking point difference
            IdxNumLT     = Temperatures2Plot1.index(BPRowL)
            IdxNumHT     = Temperatures2Plot2.index(BPRowL)
            print("\t Breaking point...")
            print("\t CP at "+str(float(BPRowLZ))+" 1st set NASA:\t%.3f" % (round(CPs1[IdxNumLT],3)))
            print("\t CP at "+str(BPRowLZ)+" 2nd set NASA:\t%.3f" % (round(CPs2[IdxNumHT],3)))
            AVECPs = abs((CPs2[IdxNumHT] - CPs1[IdxNumLT])/CPs1[IdxNumLT])*100
            print("\t AbsError CPs: %.3f " % (round(AVECPs,3)),"%")
            print("\t H at "+str(BPRowLZ)+" 1st set NASA:\t%.3f" % (round(Hs1[IdxNumLT],3)))
            print("\t H at "+str(BPRowLZ)+" 2nd set NASA:\t%.3f" % (round(Hs2[IdxNumHT],3)))
            AVEHs = abs((Hs2[IdxNumHT] - Hs1[IdxNumLT])/Hs1[IdxNumLT])*100
            print("\t AbsError Hs:  %.3f " % (round(AVEHs,3)),"%")
            print("\t S at "+str(BPRowLZ)+" 1st set NASA:\t%.3f" % (round(Ss1[IdxNumLT],3)))
            print("\t S at "+str(BPRowLZ)+" 2nd set NASA:\t%.3f" % (round(Ss2[IdxNumHT],3)))
            AVESs = abs((Ss2[IdxNumHT] - Ss1[IdxNumLT])/Ss1[IdxNumLT])*100
            print("\t AbsError Ss:  %.3f " % (round(AVESs,3)),"%")
            print("\t")
            print("\t" + "--"*10)
            print("\t Max. Cp1 value:")
            print("\t", max(CPs1))
            print("\t Max. Cp2 value:")
            print("\t", max(CPs2))
            #print("\t AbsError CPs: " + str(round(AVECPs,2) + "%"))
            #BreakPointLT = 
            #
            #plt.plot(Temperatures2Plot1,CPs1, '--', color="black",label=" HT",lw=3)
            # Plotter
            for y in range(2):
                Title = Names2Plot[u]
                title    = "Specie name: " + Title + "\nfile name:" + self.InputFile
                plt.title(title, fontsize=14, weight= 'bold')
                plt.rcParams["font.weight"] = "bold"
                plt.tick_params(axis= 'both', direction='out', length=4, width=2, labelsize=10)
                #
                plt.rcParams["axes.labelweight"] = "bold"
                plt.rcParams["axes.linewidth"] = "3"
                #
                plt.plot(Temperatures2Plot1,CPs1, '-', color="black",label=" LT",lw=3)
                #print(Temperatures2Plot2)
                #print(CPs2)
                plt.plot(Temperatures2Plot2,CPs2, '--', color="black",label=" HT",lw=3)
                IndexNumber1000K = Temperatures2Plot1.index(BPRowL)
                plt.scatter(BPRowL, CPs1[IndexNumber1000K], label=" BP at "+str(BPRowLZ)+" K", s=150, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.ylabel('Heat capacity (cal/mol*K)', fontsize=14, weight= 'bold'); 
                plt.xlabel('Temperature / K', fontsize=14, weight= 'bold'); 
                #plt.ylabel('IDT / $\mu$s', fontsize=16, weight= 'bold');
                plt.legend();
                #plt.tight_layout()
                #plt.show()
                try:
                    Plot = plt.savefig(OutputPlotPath + "\\Plots\\Single\\" + self.InputFile + "\\" + str(Title) + "_" + str(self.InputFile) + "_CP.png")
                    plt.close()
                    if y == 0:
                       print("\t " + str(Title) + "_" + str(self.InputFile) + "_CP.png saved.")
                    else:
                       pass
                except:
                    print("\t Something went wrong, error saving Cp plot...")
                #print(Hs1, Hs2)
                Title = Names2Plot[u]
                title    = "Specie name: " + Title + "\nfile name:" + self.InputFile
                plt.title(title, fontsize=14, weight= 'bold')
                plt.plot(Temperatures2Plot1,Hs1, "-", color="black",label=" LT",lw=3)
                plt.plot(Temperatures2Plot2,Hs2, "--", color="black",label=" HT",lw=3)
                IndexNumber1000K = Temperatures2Plot1.index(BPRowL)
                plt.scatter(BPRowL, Hs1[IndexNumber1000K], label=" BP at "+str(BPRowLZ)+" K", s=150, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.ylabel('Enthalpy (kcal/mol)', fontsize=14, weight= 'bold'); 
                plt.xlabel('Temperature / K', fontsize=14, weight= 'bold'); 
                #plt.ylabel('IDT / $\mu$s', fontsize=16, weight= 'bold');
                plt.legend();
                #plt.tight_layout()
                #plt.show()
                try:
                    Plot = plt.savefig(OutputPlotPath + "\\Plots\\Single\\" + self.InputFile + "\\" + str(Title) + "_" + str(self.InputFile) + "_H.png")
                    plt.close()
                    if y == 0:
                       print("\t " + str(Title) + "_" + str(self.InputFile) + "_H.png  saved.")
                    else:
                       pass
                except:
                    print("\t Something went wrong, error saving H plot...")
                #print(Ss1, Ss2)
                Title = Names2Plot[u]
                title    = "Specie name: " + Title + "\nfile name:" + self.InputFile
                plt.title(title, fontsize=14, weight= 'bold')
                plt.plot(Temperatures2Plot1,Ss1, "-", color="black",label=" LT",lw=3)
                plt.plot(Temperatures2Plot2,Ss2, "--", color="black",label=" HT",lw=3)
                IndexNumber1000K = Temperatures2Plot1.index(BPRowL)
                plt.scatter(BPRowL, Ss1[IndexNumber1000K],color="black",label=" BP at "+str(BPRowLZ)+" K", s=150, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.ylabel('Entropy (cal/mol*K)', fontsize=14, weight= 'bold'); 
                plt.xlabel('Temperature / K', fontsize=14, weight= 'bold'); 
                #plt.ylabel('IDT / $\mu$s', fontsize=16, weight= 'bold');
                plt.legend();
                #plt.tight_layout()
                #plt.show()
                try:
                    Plot = plt.savefig(OutputPlotPath + "\\Plots\\Single\\" + self.InputFile + "\\" + str(Title) + "_" + str(self.InputFile) + "_S.png")
                    plt.close()
                    if y == 0:
                       print("\t " + str(Title) + "_" + str(self.InputFile) + "_S.png  saved.")
                    else:
                       pass
                except:
                    print("\t Something went wrong, error saving S plot...")
            #
            PathPlots = OutputPlotPath + "\\Plots\\Single\\" + self.InputFile + "\\"
            #print(PathPlots)
            GraphsChecker = glob.glob(PathPlots + "\\*.png")
            if  GraphsChecker == 0:
                print("\t Not graphs were found on directory folder...")
                print("\t " + GraphsChecker)
            else:
                pass
            # 
            if AVECPs > 0.1 or AVEHs  > 0.1 or AVESs  > 0.1:
               print("Names2Plot[u]: ",Names2Plot[u])
               Xpecie.append(Names2Plot[u]) 
               XtempBP.append(BPRowLZ) 
               Abs_ErrorCP.append(AVECPs) 
               Abs_ErrorH.append(AVEHs) 
               Abs_ErrorS.append(AVESs) 
            else:
                pass
            #
        print("Xpecie: ", Xpecie)
        #
        if  len(Xpecie) == 0:
            pass
        else:
            CWD2 = os.getcwd()
            OutputPlotPath = CWD2 + "\\OutputsDir\\Plots\\Single\\" 
            DFdata = pd.DataFrame({"SpecieName":Xpecie,"BreakingPoint/K":Xtemp,"AbsErrorCP":Abs_ErrorCP,"AbsErrorH":Abs_ErrorH,"AbsErrorS":Abs_ErrorS})
            #DFdata = pd.DataFrame()
            DFdata.to_csv(PathPlots + "ReportDiscontinuities.log")
        #
        return self

    def MultiPlotter(self, ListOfFiles):
        """
        This definition plots every specie in *.dat files in multi-trend graphs
        providing the Low Temperature (LT) NASA polynomials and the High Temperature
        (HT) NASA polynomials along with the breaking point and the Absolute Percentage
        Relative Error (APRE) between these two set of polynomial at 1000 K (join point). 
        
        The 1st set/collection of data (NASA polynomials extraceted) are visualise as
        Solid Lines (SL) while the 2nd set (comparison) will be plotted in the same color
        code but with different type of line, Dotted Line (DL).
        
        Color code is:
        1st NASA polynomial set data color is "Black".
        2nd NASA polynomial set data color is "Blue".
        Breaking point scatter data color is "Red".
        """
        # date and time
        today = date.today()
        # dd_mm_YY
        dia = today.strftime("%d_%m_%Y")
        print("\t Plotting: ")
        print("\t")
        print("\t " + str(self.InputFile))
        print("\t")
        CWD = os.getcwd()
        PandasFile1 = pd.read_table(CWD + "\\OutputsDir\\" + ListOfFiles[0], names=['col'])
        PandasFile2 = pd.read_table(CWD + "\\OutputsDir\\" + ListOfFiles[1], names=['col'])
        SpeciesThermos = []
        #print("\t God Dammn!")

        for k in range(len(PandasFile1)):
            FirstData  = str(PandasFile1.iloc[k,0]).split()[0]
            SecondData = str(PandasFile1.iloc[k,0]).split()
            #print(FirstData)
            #print(SecondData)
            if   "C" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "H" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "HE" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "AR" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "H" in FirstData and "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "H" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "O" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData and "O" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "C" in FirstData and "H" in FirstData and "O" in FirstData and "N" in FirstData:
                 SpeciesThermos.append(FirstData)
            elif "0G" in SecondData:
                 SpeciesThermos.append(FirstData)
            elif "G" in SecondData:
                 SpeciesThermos.append(FirstData)
            elif "1 " in SecondData:
                 SpeciesThermos.append(FirstData)
            else:
                 #SpeciesThermos.append(FirstData)
                 pass #print("\t Species not found!")
        
        print("\t The next is a list of species found in the file")
        SpeciesThermos = list(set(SpeciesThermos))
        if   "THERMO" in SpeciesThermos: 
             SpeciesThermos.remove("THERMO")
        else:
             pass
        if   "END" in SpeciesThermos: 
             SpeciesThermos.remove("END")
        else:
             pass
        #print(SpeciesThermos)
        for j in SpeciesThermos:
            print("\t "+j, end=' ')
        print("\t ")
        HowMany2Plot = input("\t How many species do you want to plot? (Type 'all' to plot them all)\n\t ") or 0
        Names2Plot = []
        try:
            HowMany2Plot = int(HowMany2Plot)
            if   HowMany2Plot == 0:
                 print("\t 0 specie provided, no plot(s) generated!.")
                 pass
            else:
                 for t in range(HowMany2Plot):
                    while True:
                        try:
                            FileName1  = input("\t Species #" + str(t+1) +" name:")
                            if FileName1 in SpeciesThermos:
                               Names2Plot.append(FileName1)
                               break
                            else:
                                print("\t> " + str(FileName1) + " specie not in the list!")
                                continue
                        except:
                            print("\t> " + str(FileName1) + " specie not in the list!")
                            continue
                    
        except:
            if   HowMany2Plot in ["ALL", "all", "All", "ALl", "alL", "aLL"]:
                 HowMany2Plot = len(SpeciesThermos)
                 for t in range(HowMany2Plot):
                    Names2Plot.append(SpeciesThermos[t])
            else:
                 print("\t> " + str(FileName1) + " is not a valid answer, type a integer next time.")
                 pass                 
        # Extraction of NASA polynomial values from dat files
        for u in range(len(Names2Plot)):
            print("\t" + "--"*10)
            print("\t Plotting NASA polynomials for: " + Names2Plot[u])
            print("\t" + ".."*10)
            SpecieIndex = int((str(PandasFile1[PandasFile1.col.str.contains(Names2Plot[u] + " ") == True]).split()[1]))
            print("\t Polynomials extracted:")
            #
            BPRowLr      = []
            BPRow        = (PandasFile1.iloc[SpecieIndex,:]).tolist()
            BPRow2       = list(BPRow[0])
            for w in range(len(BPRow2)):
                if w >= 66 and w < 77:
                   BPRowLr.append(BPRow2[w])
            #print(BPRow)
            #print(BPRow2)
            #print(BPRowLr)
            FirstRow = (PandasFile1.iloc[SpecieIndex+1,:]).tolist()
            SecondRow = (PandasFile1.iloc[SpecieIndex+2,:]).tolist()
            ThirdRow = (PandasFile1.iloc[SpecieIndex+3,:]).tolist()
            PreFirstRow  = list(FirstRow[0])
            PreSecondRow = list(SecondRow[0])
            PreThirdRow  = list(ThirdRow[0])
            # Extracting polynomial coefficients from b1 to b5
            b1 = []; b2 = []; b3 = []; b4 = []; b5 = []; 
            for p in range(len(PreFirstRow)):
                if   p < 15:
                     b1.append(PreFirstRow[p])
                elif p >= 15 and p < 30:
                     b2.append(PreFirstRow[p])
                elif p >= 30 and p < 45:
                     b3.append(PreFirstRow[p])
                elif p >= 45 and p < 60:
                     b4.append(PreFirstRow[p])
                elif p >= 60 and p < 75:
                     b5.append(PreFirstRow[p])
            # Extracting polynomial coefficients from b6 to a3
            b6 = []; b7 = []; a1 = []; a2 = []; a3 = []; 
            for p in range(len(PreSecondRow)):
                if   p < 15:
                     b6.append(PreSecondRow[p])
                elif p >= 15 and p < 30:
                     b7.append(PreSecondRow[p])
                elif p >= 30 and p < 45:
                     a1.append(PreSecondRow[p])
                elif p >= 45 and p < 60:
                     a2.append(PreSecondRow[p])
                elif p >= 60 and p < 75:
                     a3.append(PreSecondRow[p])
            # Extracting polynomial coefficients from a4 to a7
            a4 = []; a5 = []; a6 = []; a7 = [];
            for p in range(len(PreThirdRow)):
                if   p < 15:
                     a4.append(PreThirdRow[p])
                elif p >= 15 and p < 30:
                     a5.append(PreThirdRow[p])
                elif p >= 30 and p < 45:
                     a6.append(PreThirdRow[p])
                elif p >= 45 and p < 60:
                     a7.append(PreThirdRow[p])
            b1 = "".join(b1); b2 = "".join(b2); b3 = "".join(b3); b4 = "".join(b4); b5 = "".join(b5); b6 = "".join(b6); b7 = "".join(b7); 
            a1 = "".join(a1); a2 = "".join(a2); a3 = "".join(a3); a4 = "".join(a4); a5 = "".join(a5); a6 = "".join(a6); a7 = "".join(a7); 
            #
            BPRowL  = round(float(("".join(BPRowLr)).strip()))
            BPRowLZ = float(("".join(BPRowLr)).strip()); 
            if BPRowL == 0 and BPRowLZ == 0.0:
               BPRowLr      = []
               BPRow        = (PandasFile.iloc[SpecieIndex,:]).tolist()
               BPRow2       = list(BPRow[0])
               for w in range(len(BPRow2)):
                   if w >= 65 and w < 76:
                      BPRowLr.append(BPRow2[w])
               BPRowL  = round(float(("".join(BPRowLr)).strip()))
               BPRowLZ = float(("".join(BPRowLr)).strip()); 
            else:
               pass
            print("\t 1st file NASA polynomials; ")
            print("\t Breaking Point (BP): " + str(BPRowLZ))
            print("\t a1 : " + str(a1))
            print("\t a2 : " + str(a2))
            print("\t a3 : " + str(a3))
            print("\t a4 : " + str(a4))
            print("\t a5 : " + str(a5))
            print("\t a6 : " + str(a6))
            print("\t a7 : " + str(a7))
            print("\t a8 : " + str(b1))
            print("\t a9 : " + str(b2))
            print("\t a10: " + str(b3))
            print("\t a11: " + str(b4))
            print("\t a12: " + str(b5))
            print("\t a13: " + str(b6))
            print("\t a14: " + str(b7))
            # 
            # Second file       
            # Extraction of NASA polynomial values from dat files
            try:
                print("\t" + "--"*10)
                print("\t Plotting NASA polynomials for: " + Names2Plot[u])
                print("\t" + ".."*10)
                SpecieIndex2 = int((str(PandasFile2[PandasFile2.col.str.contains(Names2Plot[u] + " ") == True]).split()[1]))
                print("\t Plynomials extracted:")
                #
                BPRowLr2      = []
                BPRow2        = (PandasFile2.iloc[SpecieIndex2,:]).tolist()
                BPRow22       = list(BPRow2[0])
                for w in range(len(BPRow22)):
                    if w >= 66 and w < 77:
                       BPRowLr2.append(BPRow22[w])
                #print(BPRow)
                #print(BPRow2)
                #print(BPRowLr)
                FirstRow2 = (PandasFile2.iloc[SpecieIndex2+1,:]).tolist()
                SecondRow2 = (PandasFile2.iloc[SpecieIndex2+2,:]).tolist()
                ThirdRow2 = (PandasFile2.iloc[SpecieIndex2+3,:]).tolist()
                PreFirstRow2  = list(FirstRow2[0])
                PreSecondRow2 = list(SecondRow2[0])
                PreThirdRow2  = list(ThirdRow2[0])
                # Extracting polynomial coefficients from b1 to b5
                b12 = []; b22 = []; b32 = []; b42 = []; b52 = []; 
                for p in range(len(PreFirstRow2)):
                    if   p < 15:
                        b12.append(PreFirstRow2[p])
                    elif p >= 15 and p < 30:
                        b22.append(PreFirstRow2[p])
                    elif p >= 30 and p < 45:
                        b32.append(PreFirstRow2[p])
                    elif p >= 45 and p < 60:
                        b42.append(PreFirstRow2[p])
                    elif p >= 60 and p < 75:
                        b52.append(PreFirstRow2[p])
                # Extracting polynomial coefficients from b6 to a3
                b62 = []; b72 = []; a12 = []; a22 = []; a32 = []; 
                for p in range(len(PreSecondRow2)):
                    if   p < 15:
                        b62.append(PreSecondRow2[p])
                    elif p >= 15 and p < 30:
                        b72.append(PreSecondRow2[p])
                    elif p >= 30 and p < 45:
                        a12.append(PreSecondRow2[p])
                    elif p >= 45 and p < 60:
                        a22.append(PreSecondRow2[p])
                    elif p >= 60 and p < 75:
                        a32.append(PreSecondRow2[p])
                # Extracting polynomial coefficients from a4 to a7
                a42 = []; a52 = []; a62 = []; a72 = [];
                for p in range(len(PreThirdRow2)):
                    if   p < 15:
                        a42.append(PreThirdRow2[p])
                    elif p >= 15 and p < 30:
                        a52.append(PreThirdRow2[p])
                    elif p >= 30 and p < 45:
                        a62.append(PreThirdRow2[p])
                    elif p >= 45 and p < 60:
                        a72.append(PreThirdRow2[p])
                b11 = "".join(b12); b21 = "".join(b22); b31 = "".join(b32); b41 = "".join(b42); b51 = "".join(b52); b61 = "".join(b62); b71 = "".join(b72); 
                a11 = "".join(a12); a21 = "".join(a22); a31 = "".join(a32); a41 = "".join(a42); a51 = "".join(a52); a61 = "".join(a62); a71 = "".join(a72); 
                #
                BPRowL2  = round(float(("".join(BPRowLr2)).strip()))
                BPRowLZ2 = float(("".join(BPRowLr2)).strip()); 
                if BPRowL2 == 0 and BPRowLZ2 == 0.0:
                   BPRowLr2      = []
                   BPRow2        = (PandasFile2.iloc[SpecieIndex2,:]).tolist()
                   BPRow22       = list(BPRow2[0])
                   for w in range(len(BPRow2)):
                       if w >= 65 and w < 76:
                          BPRowLr2.append(BPRow22[w])
                   BPRowL2  = round(float(("".join(BPRowLr2)).strip()))
                   BPRowLZ2 = float(("".join(BPRowLr2)).strip()); 
                else:
                   pass
                print("\t 2nd file NASA polynomials; ")
                print("\t Breaking Point (BP): " + str(BPRowLZ2))
                print("\t b1: " + str(b11))
                print("\t b2: " + str(b21))
                print("\t b3: " + str(b31))
                print("\t b4: " + str(b41))
                print("\t b5: " + str(b51))
                print("\t b6: " + str(b61))
                print("\t b7: " + str(b71))
                print("\t a1: " + str(a11))
                print("\t a2: " + str(a21))
                print("\t a3: " + str(a31))
                print("\t a4: " + str(a41))
                print("\t a5: " + str(a51))
                print("\t a6: " + str(a61))
                print("\t a7: " + str(a71))
            except:
                print("\t Specie: " + str(Names2Plot[u]) + " not found in 2nd file, skiping..." )
                print("\t Please be sure the files you provided contain the same species to compare vs" )
                pass
            # Plotting results
            # Converting data from NASA to thermochemistry: Cps/R, H and S to plot
            def funcCP(T, m1, m2, m3, m4, m5):
                R = 1.9872  # in cal units
                # Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
                # CP/R here produce by code
                return ((m1) + (m2 * (T)) + (m3 * (T) ** 2) + (m4 * (T) ** 3) + (m5 * (T) ** 4))*R

            def funcH(T, m1, m2, m3, m4, m5, m6):
                T2 = T * T
                T4 = T2 * T2
                R = 0.0019872  # in kcal units
                # H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
                # H here produce by code
                return (m1 + m2 * T/2 + m3 * T2/3 + m4 * T2 * T/4 + m5 * T4/5 + m6 / T) * R * T

            def funcS(T, m1, m2, m3, m4, m5, m7):
                import math
                R = 1.9#872  # in cal units
                T2 = T * T
                T4 = T2 * T2
                # S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
                # S here produce by code
                return (m1 * math.log(T) + m2 * T + m3 * T2/2 + m4 * T2 * T/3 + m5 * T4/4 + m7) * R 

            CPs1 = []; Hs1 = []; Ss1 = []; CPs2 = []; Hs2 = []; Ss2 = [];
            CPs11 = []; Hs11 = []; Ss11 = []; CPs22 = []; Hs22 = []; Ss22 = [];
            #Temperatures2Plot1 = np.arange(300, 1100, 100).tolist()
            Temperatures2Plot1  = np.arange(300, BPRowL+1, 1).tolist()
            Temperatures2Plot12 = np.arange(300, BPRowL2+1, 1).tolist()
            #print(Temperatures2Plot1)
            #Temperatures2Plot2 = np.arange(1000, 3100, 100).tolist()
            Temperatures2Plot2  = np.arange(BPRowL, 5100, 1).tolist()
            Temperatures2Plot22 = np.arange(BPRowL2, 5100, 1).tolist()
            #print(Temperatures2Plot2)
            CWD2 = os.getcwd()
            OutputPlotPath = CWD2 + "\\OutputsDir\\"
            if not os.path.exists(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\CP\\"):
                   os.makedirs(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\CP\\")
            if not os.path.exists(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\H\\"):
                   os.makedirs(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\H\\")
            if not os.path.exists(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\S\\"):
                   os.makedirs(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\S\\")
            for d1 in (Temperatures2Plot1):
                CPs1.append(round(funcCP(d1, float(a1), float(a2), float(a3), float(a4), float(a5)),5))
                Hs1.append(round(funcH(d1, float(a1), float(a2), float(a3), float(a4), float(a5), float(a6)),5))
                Ss1.append(round(funcS(d1, float(a1), float(a2), float(a3), float(a4), float(a5), float(a7)),5))
            for d12 in (Temperatures2Plot12):
                CPs11.append(round(funcCP(d12, float(a11), float(a21), float(a31), float(a41), float(a51)),5))
                Hs11.append(round(funcH(d12, float(a11), float(a21), float(a31), float(a41), float(a51), float(a61)),5))
                Ss11.append(round(funcS(d12, float(a11), float(a21), float(a31), float(a41), float(a51), float(a71)),5))
            for d2 in (Temperatures2Plot2):
                CPs2.append(round(funcCP(d2, float(b1), float(b2), float(b3), float(b4), float(b5)),5))
                Hs2.append(round(funcH(d2, float(b1), float(b2), float(b3), float(b4), float(b5), float(b6)),5))
                Ss2.append(round(funcS(d2, float(b1), float(b2), float(b3), float(b4), float(b5), float(b7)),5))
            for d22 in (Temperatures2Plot22):
                CPs22.append(round(funcCP(d22, float(b11), float(b21), float(b31), float(b41), float(b51)),5))
                Hs22.append(round(funcH(d22, float(b11), float(b21), float(b31), float(b41), float(b51), float(b61)),5))
                Ss22.append(round(funcS(d22, float(b11), float(b21), float(b31), float(b41), float(b51), float(b71)),5))
            # printing on screen the Breaking point difference
            #IdxNumLT     = Temperatures2Plot1.index(1000)
            IdxNumLT      = Temperatures2Plot1.index(BPRowL)
            IdxNumLT2     = Temperatures2Plot12.index(BPRowL2)
            #IdxNumHT     = Temperatures2Plot2.index(1000)        
            IdxNumHT      = Temperatures2Plot2.index(BPRowL)        
            IdxNumHT2     = Temperatures2Plot22.index(BPRowL2)        
            #print(CPs1, CPs2)
            #print(CPs11, CPs22)
            print("\t" + "--"*10)        
            print("\t" + "--"*46)
            print("\t 1st file:")
            print("\t CP at "+str(BPRowLZ)+" 1st set NASA:\t%.3f" % (round(CPs1[IdxNumLT],3)))
            print("\t CP at "+str(BPRowLZ)+" 2nd set NASA:\t%.3f" % (round(CPs2[IdxNumHT],3)))
            AVECPs = abs((CPs2[IdxNumHT] - CPs1[IdxNumLT])/CPs1[IdxNumLT])*100
            print("\t AbsError CPs: %.3f " % (round(AVECPs,3)) + " %")
            print("\t 2nd file:")
            print("\t CP at "+str(BPRowLZ2)+" 1st set NASA:\t%.3f" % (round(CPs11[IdxNumLT2],3)))
            print("\t CP at "+str(BPRowLZ2)+" 2nd set NASA:\t%.3f" % (round(CPs22[IdxNumHT2],3)))
            AVECPs2 = abs((CPs22[IdxNumHT2] - CPs11[IdxNumLT2])/CPs11[IdxNumLT2])*100
            print("\t AbsError CPs: %.3f " % (round(AVECPs2,3)) + " %")
            print("\t" + "-"*50)
            print("\t 1st file:")
            print("\t H at "+str(BPRowLZ)+" 1st set NASA:\t%.3f" % (round(Hs1[IdxNumLT],3)))
            print("\t H at "+str(BPRowLZ)+" 2nd set NASA:\t%.3f" % (round(Hs2[IdxNumHT],3)))
            AVEHs = abs((Hs2[IdxNumHT] - Hs1[IdxNumLT])/Hs1[IdxNumLT])*100
            print("\t AbsError Hs: %.3f " % (round(AVEHs,3)) + " %")
            print("\t 2nd file:")
            print("\t H at "+str(BPRowLZ2)+" 1st set NASA:\t%.3f" % (round(Hs11[IdxNumLT2],3)))
            print("\t H at "+str(BPRowLZ2)+" 2nd set NASA:\t%.3f" % (round(Hs22[IdxNumHT2],3)))
            AVEHs2 = abs((Hs22[IdxNumHT2] - Hs11[IdxNumLT2])/Hs11[IdxNumLT2])*100
            print("\t AbsError Hs: %.3f " % (round(AVEHs2,3)) + " %")
            print("\t" + "-"*50)
            print("\t 1st file:")
            print("\t S at "+str(BPRowLZ)+" 1st set NASA:\t%.3f" % (round(Ss1[IdxNumLT],3)))
            print("\t S at "+str(BPRowLZ)+" 2nd set NASA:\t%.3f" % (round(Ss2[IdxNumHT],3)))
            AVESs = abs((Ss2[IdxNumHT] - Ss1[IdxNumLT])/Ss1[IdxNumLT])*100
            print("\t AbsError Ss: %.3f " % (round(AVESs,3)) + " %")
            print("\t 2nd file:")
            print("\t S at "+str(BPRowLZ2)+" 1st set NASA:\t%.3f" % (round(Ss11[IdxNumLT2],3)))
            print("\t S at "+str(BPRowLZ2)+" 2nd set NASA:\t%.3f" % (round(Ss22[IdxNumHT2],3)))
            AVESs2 = abs((Ss22[IdxNumHT2] - Ss11[IdxNumLT2])/Ss11[IdxNumLT2])*100
            print("\t AbsError Ss: %.3f " % (round(AVESs2,3)) + " %")
            print("\t")
            print("\t" + "--"*10)        
            for y in range(2):
                Title = Names2Plot[u]
                title    = "Specie name: " + Title + "\n" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1])
                plt.title(title, fontsize=14, weight= 'bold')
                plt.rcParams["font.weight"] = "bold"
                plt.tick_params(axis= 'both', direction='out', length=4, width=2, labelsize=10)
                #
                plt.rcParams["axes.labelweight"] = "bold"
                plt.rcParams["axes.linewidth"] = "3"
                #
                # 1st file
                plt.plot(Temperatures2Plot1,CPs1,'-', color="black",  label=" 1st LT",lw=3)
                # 2nd file
                plt.plot(Temperatures2Plot12,CPs11, '-', color="red", label=" 2nd LT",lw=3)
                # 1st file
                plt.plot(Temperatures2Plot2,CPs2, '--', color="black", label=" 1st HT",lw=3)
                # 2nd file
                plt.plot(Temperatures2Plot22,CPs22, '--', color="red", label=" 2nd HT",lw=3)
                IndexNumber1000K  = Temperatures2Plot1.index(BPRowL)
                IndexNumber1000K2 = Temperatures2Plot12.index(BPRowL2)
                #print("CP 2nd file")
                #plt.scatter(1000, CPs1[IndexNumber1000K], label=' Breaking point', s=200, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.scatter(BPRowL, CPs1[IndexNumber1000K], label=" 1st BP at "+str(BPRowLZ)+" K", s=150, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.scatter(BPRowL2, CPs11[IndexNumber1000K2], label=" 2nd BP at "+str(BPRowLZ2)+" K", s=150, facecolors='none', edgecolors='red', marker = "s", lw=2)
                plt.ylabel('Heat capacity (cal/mol*K)', fontsize=14, weight= 'bold'); 
                plt.xlabel('Temperature / K', fontsize=14, weight= 'bold'); 
                plt.legend();
                try:
                    Plot = plt.savefig(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\CP\\" + str(Title) + "_CP.png")
                    plt.close()
                    if y == 0:
                       print("\t " + str(Title) + "_CP.png saved.")
                    else:
                       pass
                except:
                    print("\t Something went wrong, error saving Cp plot...")
                #print(Hs1, Hs2)
                Title = Names2Plot[u]
                title    = "Specie name: " + Title + "\n" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1])
                plt.title(title, fontsize=14, weight= 'bold')
                # 1st file
                plt.plot(Temperatures2Plot1,Hs1, '-', color="black",label=" 1st LT",lw=3)
                # 2nd file
                plt.plot(Temperatures2Plot12,Hs11, '-', color="red",label=" 2nd LT",lw=3)
                # 1st file
                plt.plot(Temperatures2Plot2,Hs2, '--', color="black",label=" 1st HT",lw=3)
                # 2nd file                
                plt.plot(Temperatures2Plot22,Hs22, '--', color="red",label=" 2nd HT",lw=3)
                IndexNumber1000K  = Temperatures2Plot1.index(BPRowL)
                IndexNumber1000K2 = Temperatures2Plot12.index(BPRowL2)
                plt.scatter(BPRowL, Hs1[IndexNumber1000K], label=" 1st BP at "+str(BPRowLZ)+" K", s=150, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.scatter(BPRowL2, Hs11[IndexNumber1000K2], label=" 2nd BP at "+str(BPRowLZ2)+" K", s=150, facecolors='none', edgecolors='red', marker = "s", lw=2)
                plt.ylabel('Enthalpy (kcal/mol)', fontsize=14, weight= 'bold'); 
                plt.xlabel('Temperature / K', fontsize=14, weight= 'bold'); 
                plt.legend();
                #plt.close()
                try:
                    Plot = plt.savefig(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\H\\" + str(Title) + "_H.png")
                    plt.close()
                    if y == 0:
                       print("\t " + str(Title) + "_H.png  saved.")
                    else:
                       pass
                except:
                    print("\t Something went wrong, error saving H plot...")
                ##print(Ss1, Ss2)
                Title = Names2Plot[u]
                title    = "Specie name: " + Title + "\n" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1])
                plt.title(title, fontsize=14, weight= 'bold')
                # 1st file
                plt.plot(Temperatures2Plot1,Ss1, '-', color="black",label=" 1st LT",lw=3)
                # 2nd file
                plt.plot(Temperatures2Plot12,Ss11, '-', color="red",label=" 2nd LT",lw=3)
                # 1st file
                plt.plot(Temperatures2Plot2,Ss2, '--', color="black",label=" 1st HT",lw=3)
                # 2nd file
                plt.plot(Temperatures2Plot22,Ss22, '--', color="red",label=" 2nd HT",lw=3)
                IndexNumber1000K  = Temperatures2Plot1.index(BPRowL)
                IndexNumber1000K2 = Temperatures2Plot12.index(BPRowL2)
                plt.scatter(BPRowL, Ss1[IndexNumber1000K], label=" 1st BP at "+str(BPRowLZ)+" K", s=150, facecolors='none', edgecolors='black', marker = "s", lw=2)
                plt.scatter(BPRowL2, Ss11[IndexNumber1000K2], label=" 2nd BP at "+str(BPRowLZ2)+" K", s=150, facecolors='none', edgecolors='red', marker = "s", lw=2)
                plt.ylabel('Entropy (cal/mol*K)', fontsize=14, weight= 'bold'); 
                plt.xlabel('Temperature / K', fontsize=14, weight= 'bold'); 
                plt.legend();
                #plt.close()
                try:
                    Plot = plt.savefig(OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\S\\" + str(Title) + "_S.png")
                    plt.close()
                    if y == 0:
                       print("\t " + str(Title) + "_S.png  saved.")
                    else:
                       pass
                except:
                    print("\t Something went wrong, error saving S plot...")
            PathPlots = OutputPlotPath + "\\Plots\\multi\\" + str(dia) + "\\" + str(ListOfFiles[0]) + "_VS_" + str(ListOfFiles[1]) + "\\"
            GraphsChecker = glob.glob(PathPlots + "\\*\\*.png")
            if  GraphsChecker == 0:
                print("\t Not graphs were found on directory folder...")
                print("\t " + GraphsChecker)
            else:
                pass
            
        return self        


