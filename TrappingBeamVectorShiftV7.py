#===========================================================
# Name: TrappingBeamVectorShiftV7.py
# Author: Gordon Arrowsmith-kron
# Date of Creation: 03/22/2023
#===========================================================
import sympy.physics.wigner as wig
import matplotlib.pyplot as plt
import numpy


#Get the relevant definitions

#term of summation for MPrime
def MSummationCR(J, JPrime, I, M, MPrime, F, FPrime, EpsilonL, EpsilonR):
	tempexp = 4*FPrime + 4*F - MPrime - M + JPrime + J + 2*I + 2
	epsilonSquared = 0
	if (M - MPrime) == 1:
		epsilonSquared = -1*(EpsilonL**2)
	if (M - MPrime == -1):
		epsilonSquared = -1*(EpsilonR**2)
	return ((-1)**tempexp)*epsilonSquared*(wig.wigner_3j(FPrime, 1, F, -MPrime, MPrime - M, M)**2)
	
def MSummationR(J, JPrime, I, M, MPrime, F, FPrime, EpsilonL, EpsilonR):
	tempexp = 4*FPrime + 4*F - MPrime - M + JPrime + J + 2*I + 2
	epsilonSquared = 0
	if (M - MPrime) == 1:
		epsilonSquared = -1*(EpsilonR**2)
	if (M - MPrime == -1):
		epsilonSquared = -1*(EpsilonL**2)
	return ((-1)**tempexp)*epsilonSquared*(wig.wigner_3j(FPrime, 1, F, -MPrime, MPrime - M, M)**2)


#note: Need to put a factor of -e^2 E_0^2/(4m)	
def CR_evaluate(J, JPrime, I, Osc, E, F, M, FPrime, Omega, EpsilonL, EpsilonR):
	tempexp = J - JPrime
	#Create list of MPrime dependent on what FPrime is
	MPrimeArray = []
	i = -FPrime
	while i <= FPrime:
		MPrimeArray.append(i)
		i = i + 1
	#Sum over the various values of MPrime
	Sum = 0
	i = 0
	while i < len(MPrimeArray):
		Sum = Sum + MSummationCR(J, JPrime, I, M, MPrimeArray[i], F, FPrime, EpsilonL, EpsilonR)
		i = i + 1
	return (2*F+1)*(2*FPrime+1)*((-1)**tempexp)*3*(2*J + 1)*Osc*Sum*(wig.wigner_6j(JPrime, FPrime, I, F, J, 1)**2)/(2*E*(E - Omega))

def R_evaluate(J, JPrime, I, Osc, E, F, M, FPrime, Omega, EpsilonL, EpsilonR):
	tempexp = J - JPrime
	#Create list of MPrime dependent on what FPrime is
	MPrimeArray = []
	i = -FPrime
	while i <= FPrime:
		MPrimeArray.append(i)
		i = i + 1
	#Sum over the various values of MPrime
	Sum = 0
	i = 0
	while i < len(MPrimeArray):
		Sum = Sum + MSummationR(J, JPrime, I, M, MPrimeArray[i], F, FPrime, EpsilonL, EpsilonR)
		i = i + 1
	return (2*F+1)*(2*FPrime+1)*((-1)**tempexp)*3*(2*J + 1)*Osc*Sum*(wig.wigner_6j(JPrime, FPrime, I, F, J, 1)**2)/(2*E*(E + Omega))

#Put it all together

def DeltaE(J, JPrime, I, Osc, E, F, M, FPrime, Omega, EpsilonL, EpsilonR):
	return R_evaluate(J, JPrime, I, Osc, E, F, M, FPrime, Omega, EpsilonL, EpsilonR) + CR_evaluate(J, JPrime, I, Osc, E, F, M, FPrime, Omega, EpsilonL, EpsilonR)
	
		
#Conversion factors, constants for cgs gaussian units

hbar = 1.054*(10**-27)
c = 299792458*100 # cm per second

FrequencyFactor = 2*3.1415*(2.998)*(10**10) # Hz/cm^-1
ChargeFactor = 2.998 * (10**9) #qG/qI
PowerFactor = (10**-7) #1erg per second/ Watts
	
electronCharge = 1.602 * (10**-19) * ChargeFactor # in gaussian CGS units, Fr
electronMass = 9.1093837 * (10**-28) #in grams

BoltzmannConstant = 1.380659 * (10**-16) #cgs units, erg/K

TrapDepth = 100*(10**-6) #Kelvin
TrapEnergy = -TrapDepth*BoltzmannConstant #Note the minus sign- the energy is BELOW zero, due to it being a trap.
	
	
#Define things

F = 4
M = {-4, -3, -2, -1, 0, 1, 2, 3, 4}

#First transition values
FPrime1 = 3
FPrime2 = 4
FPrimeArray1 = {3, 4}
Osc1 = .35

J = .5
JPrime1 = .5

I = 3.5

k1 = 11178.2686
omega1 = k1*FrequencyFactor


EpsilonL = 1

EpsilonR = numpy.sqrt(1 - EpsilonL**2)

#Second transition values
FPrime3 = 2
FPrime4 = 3
FPrime5 = 4
FPrime6 = 5

FPrimeArray2 = {2, 3, 4, 5}

Osc2 = .72

JPrime2 = 1.5

k2 = 11732.3079
omega2 = k2*FrequencyFactor

#Create Omega Array
omegaArray = []
i = 0
while i < 100:
	omega = omega1 * i /100.0
	omegaArray.append(omega)
	i = i + 1
	
#Make the constant that the complicated stuff is multiplied by. Going to use the profile of the ODT holding beam
P = 30 #Watts
Pcgs = P * PowerFactor
w_0 = 130*(10**-4) # cm

ESquared = Pcgs*16/(c*(w_0**2))

Constant = (electronCharge**2)*ESquared/(4*electronMass)

print(ESquared)
print(Constant)	
#Loop over all transitions
for m in M:
	deltaEArray = []
	sum = 0
	i = 0
	while(i < 100):
		sum = 0
		for FPrime1 in FPrimeArray1:
			sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F, m, FPrime1, omegaArray[i], EpsilonL, EpsilonR)
		for FPrime2 in FPrimeArray2:
			sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F, m, FPrime2, omegaArray[i], EpsilonL, EpsilonR)
		deltaEArray.append(Constant*sum/hbar)
		i = i + 1
	plt.plot(omegaArray, deltaEArray, label = m)
	
plt.legend()
plt.title("Delta Omega, different F, M")
plt.yscale("log")
plt.show()

fracArray = []
i = 0
while(i < 100):
	fracArray.append(omegaArray[i]/omega1)
	i = i + 1

#Think I found an error with all my previous codes; going to look at difference again
deltaE0Array = []
i = 0
while(i < 100):
	sum = 0
	for FPrime1 in FPrimeArray1:
		sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F, 0, FPrime1, omegaArray[i], EpsilonL, EpsilonR)
	for FPrime2 in FPrimeArray2:
		sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F, 0, FPrime2, omegaArray[i], EpsilonL, EpsilonR)
	deltaE0Array.append(Constant*sum/hbar)
	i = i + 1

deltaE1Array = []
i = 0
while(i < 100):
	sum = 0
	for FPrime1 in FPrimeArray1:
		sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F, 1, FPrime1, omegaArray[i], EpsilonL, EpsilonR)
	for FPrime2 in FPrimeArray2:
		sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F, 1, FPrime2, omegaArray[i], EpsilonL, EpsilonR)
	deltaE1Array.append(Constant*sum/hbar)
	i = i + 1

deltaE2Array = []
i = 0
while(i < 100):
	sum = 0
	for FPrime1 in FPrimeArray1:
		sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F, 2, FPrime1, omegaArray[i], EpsilonL, EpsilonR)
	for FPrime2 in FPrimeArray2:
		sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F, 2, FPrime2, omegaArray[i], EpsilonL, EpsilonR)
	deltaE2Array.append(Constant*sum/hbar)
	i = i + 1

deltaE3Array = []
i = 0
while(i < 100):
	sum = 0
	for FPrime1 in FPrimeArray1:
		sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F, 3, FPrime1, omegaArray[i], EpsilonL, EpsilonR)
	for FPrime2 in FPrimeArray2:
		sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F, 3, FPrime2, omegaArray[i], EpsilonL, EpsilonR)
	deltaE3Array.append(Constant*sum/hbar)
	i = i + 1

deltaE4Array = []
i = 0
while(i < 100):
	sum = 0
	for FPrime1 in FPrimeArray1:
		sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F, 4, FPrime1, omegaArray[i], EpsilonL, EpsilonR)
	for FPrime2 in FPrimeArray2:
		sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F, 4, FPrime2, omegaArray[i], EpsilonL, EpsilonR)
	deltaE4Array.append(Constant*sum/hbar)
	i = i + 1	
	
	
ratioArray = []
diffArray1 = []
diffArray2 = []
diffArray3 = []
diffArray4 = []
i = 0
while(i < 100): 
	diffArray1.append(deltaE1Array[i] - deltaE0Array[i])
	diffArray2.append(deltaE2Array[i] - deltaE0Array[i])
	diffArray3.append(deltaE3Array[i] - deltaE0Array[i])
	diffArray4.append(deltaE4Array[i] - deltaE0Array[i])
	ratioArray.append((deltaE2Array[i] - deltaE0Array[i])/(deltaE1Array[i] - deltaE0Array[i]))
	i = i + 1
	
plt.plot(fracArray[1:], diffArray4[1:], label =  r"$M_F$ = 4")
plt.plot(fracArray[1:], diffArray3[1:], label =  r"$M_F$ = 3")
plt.plot(fracArray[1:], diffArray2[1:], label =  r"$M_F$ = 2")
plt.plot(fracArray[1:], diffArray1[1:], label =  r"$M_F$ = 1")

plt.legend()
plt.yscale("log")
plt.title(r"Difference in Vector Shift of Various $M_F$ sublevels and $M_F$ = 0, $\epsilon_L = 1$")
plt.xlabel(r"$\omega / \omega_1$")
plt.ylabel(r"Difference in Vector Shift (Hz)")
plt.show()



plt.plot(omegaArray, ratioArray)
plt.show()





#Make a list of values from 0 to 1
PolarizationArray = []
i = 0
while(i <= 100):
	PolarizationArray.append(i/100)
	i = i + 1

	
HoldingOmega = 2*3.1415*c/(1550*(10**(-7)))
#Want to make a plot with the vector shift as a function of the circular polarization
#loop over all starting F
FArray = {4}
for F0 in FArray:
	#Create an Array of M for a given F
	m = -F0
	MArray = []
	while(m <= F0):
		MArray.append(m)
		m = m + 1
	for m0 in MArray:
		deltaEPolarizationArray = []
		i = 0
		while(i <= 100):
			sum = 0
			EpsilonL = PolarizationArray[i]
			EpsilonR =  numpy.sqrt(1 - (EpsilonL**2))
			for FPrime1 in FPrimeArray1:
				sum = sum + DeltaE(J, JPrime1, I, Osc1, omega1, F0, m0, FPrime1, HoldingOmega, EpsilonL, EpsilonR)
			for FPrime2 in FPrimeArray2:
				sum = sum + DeltaE(J, JPrime2, I, Osc2, omega2, F0, m0, FPrime2, HoldingOmega, EpsilonL, EpsilonR)
			deltaEPolarizationArray.append(Constant*sum/hbar)
			i = i + 1
		plt.plot(PolarizationArray, deltaEPolarizationArray, label = r"$M_F$ = " + str(m0))

plt.legend()
plt.xlabel(r"$\epsilon_L$")
plt.ylabel(r"Vector Shift $\omega_0$ (Hz)")
plt.title("Total Vector Shift as Function of Circular Polarization, F = 4, $\lambda  = 1550 nm$")
plt.show()

#Make array for omega/omega1
fracArray = []
i = 0
while(i < 100):
	fracArray.append(omegaArray[i]/omega1)
	i = i + 1



plt.plot(fracArray, deltaE0Array)
plt.xlabel(r"$\omega / \omega_1$")
plt.ylabel(r"vector Shift $\omega_0$ (Hz)")
plt.title(r"Total Vector Shift as Function of ODT Trapping Light, F = 4, $M_F$ = 0, $\omega_1$ = 210.588 THz, $\epsilon_L = 1$")
plt.show()
