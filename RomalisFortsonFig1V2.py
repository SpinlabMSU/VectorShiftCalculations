#===============================================================================================
# Name: RomalisFortsonFig1V2
# Author: Gordon Arrowsmith-Kron
# Purpose: To create figure 1 from Romalis-Fortson, this time with Hg as well
# Date of Creation: 03/22/2023
#===============================================================================================

import numpy as np
import matplotlib.pyplot as plt

#Values for energy levels

omega1Cs = 11178.2686 #cm^-1
omega3Cs = 11732.3079 #cm^-1

omega1Hg = 54068.781
omega3Hg = 39412.300


#Conversion factor for wave number to angular frequency
FrequencyFactor = 2*3.1415*(2.998)*(10**10) # Hz/cm^-1

omega1CsHz = omega1Cs*FrequencyFactor
omega3CsHz = omega3Cs*FrequencyFactor

omega1HgHz = omega1Hg*FrequencyFactor
omega3HgHz = omega3Hg*FrequencyFactor


#oscillator strengths 

f1Cs = .35
f3Cs = .72

f1Hg = 1.2
f3Hg = .025

#Gaussian conversion factors

ChargeFactor =  2.998 * (10**9) #qG/qI
EFieldFactor = (1/2.998)*(10**-4) #EG/EI
PowerFactor = 10**7 #erg/s / W

#Constants in cgs units

electronMass = 9.1093837 * (10**-28) #in grams
electronCharge = 1.602 * (10**-19) * ChargeFactor # in gaussian CGS units, Fr
BoltzmannConstant = 1.380659 * (10**-16) #cgs units, erg/K

TrapDepth = 100*(10**-6) #Kelvin
TrapEnergy = -TrapDepth*BoltzmannConstant #Note the minus sign- the energy is BELOW zero, due to it being a trap.

spotSize = 100*(10**-4) #Diameter of the beam, in cm

c = 2.998*(10**10) #Speed of light, cm/s


#sub function for determining electric field strength
def sumTerm(f, omegaE, omega):
	return(f/(omegaE**2 - omega**2))

#find the electric field strength for Cs
def ESquaredCs(U, omega):
	Sum = sumTerm(f1Cs, omega1CsHz, omega) + sumTerm(f3Cs, omega3CsHz, omega)
	return(-4*electronMass*U/(electronCharge**2 * Sum))

#find the electric field strength for Hg
def ESquaredHg(U, omega):
	Sum = sumTerm(f1Hg, omega1HgHz, omega) + sumTerm(f3Hg, omega3HgHz, omega)
	return(-4*electronMass*U/(electronCharge**2 * Sum))
	
#Make the frequency arrays
omegaArrayCs = []
omegaArrayHg = []
omegaRatioArray = []
i = 0
while (i < 1000):
	omegaArrayCs.append(i*omega1CsHz/1000)
	omegaArrayHg.append(i*omega3HgHz/1000)
	omegaRatioArray.append(i/1000)
	i = i + 1

#Make ESquared Arrays
ESquaredArrayCs = []
ESquaredArrayHg = []
i = 0
while (i < 1000):
	ESquaredArrayCs.append(ESquaredCs(TrapEnergy, omegaArrayCs[i]))
	ESquaredArrayHg.append(ESquaredHg(TrapEnergy, omegaArrayHg[i]))
	i = i + 1
	
#Now, convert to power necessary in beam:

PowerArrayCs = []
PowerArraySICs = []
PowerArrayHg = []
PowerArraySIHg = []
i = 0
while (i < 1000):
	Power = c * ESquaredArrayCs[i]*(spotSize**2)/16
	PowerArrayCs.append(Power) #In CGS Units
	PowerArraySICs.append(Power/PowerFactor) # In Watts
	Power = c*ESquaredArrayHg[i]*(spotSize**2)/16
	PowerArrayHg.append(Power)
	PowerArraySIHg.append(Power/PowerFactor)
	i = i + 1

	
	
plt.plot(omegaRatioArray, ESquaredArrayCs)
plt.show()
	
plt.plot(omegaRatioArray, PowerArraySICs)
plt.title("Cs Required ODT Laser Power")
plt.ylim((0, 20))
plt.xlabel(r"$\omega/\omega_1$")
plt.ylabel("Power for Cs (W)")
plt.show()
	
	
plt.plot(omegaRatioArray, PowerArraySIHg)
plt.title("Hg Required ODT Laser Power")
plt.ylim((0, 350))
plt.xlabel(r"$\omega/\omega_1$")
plt.ylabel("Power for Hg (W)")
plt.show()
	






