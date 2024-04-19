=0= READ ME =0=

#--------------------------------------------------------------
#--------------------------------------------------------------
Some GAVs for 
BD
HC
HCO
have been updated by Manik Gosh with new values and labeled as
"MG optimisation". 11/08/2021.
# ------------------------------------------------------------
# ------------------------------------------------------------

This directory includes files containing groups (GAVs) and their corresponding values.
***********************************************************************************************************************************************************
List of GAVs files included in this folder(directory):
1.-	BD.grp
2.-	CDOT.grp
3.-	CLC.grp
4.-	CYCH.grp	
5.-	HC.grp
6.-	HCN.grp
7.-	HCO.grp
8.-	INT.grp
***********************************************************************************************************************************************************
Please, if you considered we have missed any another GAVs file available, let us know.
Sergio Martinez:
email us 
s.martinez3@nuigalway.ie
sergioesmartinez@gmail.com
or hit me up at 
+353 89 965 5912
***********************************************************************************************************************************************************

An example of these files and their format is shown next:
-----------------------------------------------------------------------------------------------------------------------------------------------------------
BOND_DISSOCIATION                                                                   
138               Hf   S   Cp:300   400    500    600    800   1000   1500          
P,             101.15   2.51  -1.14  -1.64  -2.14  -2.60  -3.32  -3.89  -4.69     SMB optimisation 
IC4H9,         101.60   1.64  -0.23  -0.90  -1.54  -2.09  -3.01  -3.60  -4.59     SMB optimisation 
C2H5,          100.77   2.27  -0.21  -0.77  -1.33  -1.89  -2.80  -3.48  -4.56     SMB optimisation
NEO-C5,        101.50   3.03  -0.59  -1.32  -2.05  -2.65  -3.50  -4.06  -4.87     S&Cp:04-94
S,              98.07   4.26  -1.78  -2.98  -3.49  -3.81  -4.40  -4.78  -5.40     SMB optimisation AUG 2012 Optimised to Exp and calculations
------------------------------------------------------------------------------------------------------------------------------------------------------------

The software read line by line the files with extension "*.grp", there are 8 GAVs files already included along with this software
and can be found at "\Therm21V1.0\GroupsDir" folder, as I already mentioned. DON'T FORGET THE "," AFTER THE GAV's NAME!!!, comma is used as a separator.
As well, if you have gaps for Cp values you will have to fill them with "0.00", DON'T LET EMPTY BLANKS, SOFTWARE NEEDS A FLOAT TO READ!!!
I would recomend to use this code with option 5 to fit a polynomial to your data and refill the gaps, please read the isntructions of "how to use it"
in the Main README.txt file provided with this software.

if you want to add old GAVs, or any kind of anotation regarding the old or new GAVs please add a "#" symbol at the beginning of each line, as follows:
-------------------------------------------------------------------------------------------------------------------------------------------------------------
# old values
# 
# CY/C4O,          5.63  24.84  -4.90  -4.52  -3.43  -2.52  -1.68  -1.01   0.00
# CY/CO2,         26.11  26.4   -0.47  -1.51  -2.16  -2.51  -2.52  -1.96  -1.78    Boz/92
# CY/C2O2,        25.43  24.42  -1.7   -2.54  -3.17  -3.07  -2.09  -1.2   -0.78    Boz/92
# CY/C2O2/E       29.25  25.72  -1.39  -1.21  -1.97  -1.85  -1.40  -0.97  -0.71    jwb cy/c4o - (cyc4 - cyc4e)
# CY/C3O2,         5.35  30.62  -2.93  -3.73  -4.33  -4.35  -3.44  -2.2   -1.59    Boz/92
--------------------------------------------------------------------------------------------------------------------------------------------------------------

Otherwise the software may understand you want to count double that GAV.



