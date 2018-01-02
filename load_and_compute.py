import sys
import numpy as np
import time
from scipy import interpolate
from scipy import integrate

#********************************************************************

# Load the polarizability data ( alpha_xx and alpha_zz)

alpha_xx=np.loadtxt("matrix_xxf_trimmed.txt")
alpha_zz=np.loadtxt("matrix_zzf_trimmed.txt")
omega=np.loadtxt("freq_trimmed.txt")

distance=np.loadtxt("distance.txt")

# compute the isotropy(mean polarizability) and anisotropy

isotropy=np.absolute(2*(np.array(alpha_xx))+np.array(alpha_zz))/3
anisotropy=np.absolute(np.array(alpha_zz)-np.array(alpha_xx))


# optional save : isotropy(mean polarizability) and anisotropy as txt files
np.savetxt("isotropy.txt",isotropy)
np.savetxt("gamma.txt",anisotropy)


# check size of the arrays
print("Dimensions of isotropy matrix :",isotropy.shape)
print("Dimensions of anisotropy matrix :",anisotropy.shape)


# check size of wavelength arrays (original) and final
print("Dimensions of omega 1D array (Hartree) :",len(omega))


# convert omega (original) to nanometers
omega_nm=(1e7/(omega*219474.6313702000))

print("Range of wavelength available :",omega_nm[0],omega_nm[len(omega)-1])




# comment out this sys.exit() to proceed with the execution.
#sys.exit()

#********************************************************************

# Load the 1D array of final wavelength (in nanometers), file name is 'omega_final.txt'
omega_final=np.loadtxt("omega_final.txt")

#-------------WAVEFUNCTIONS---------------
# load the required wavefunctions which are available in the 'wavefunctions' subdirectory
# wavefunctions provided are normalized.

H2v0j0=np.loadtxt("./wavefunctions/H2v0J0_norm.txt")
H2v0j1=np.loadtxt("./wavefunctions/H2v0J1_norm.txt")
H2v0j2=np.loadtxt("./wavefunctions/H2v0J2_norm.txt")
H2v0j3=np.loadtxt("./wavefunctions/H2v0J3_norm.txt")
H2v0j4=np.loadtxt("./wavefunctions/H2v0J4_norm.txt")
H2v0j5=np.loadtxt("./wavefunctions/H2v0J5_norm.txt")
H2v0j6=np.loadtxt("./wavefunctions/H2v0J6_norm.txt")
H2v0j7=np.loadtxt("./wavefunctions/H2v0J7_norm.txt")
H2v0j8=np.loadtxt("./wavefunctions/H2v0J8_norm.txt")
H2v0j9=np.loadtxt("./wavefunctions/H2v0J9_norm.txt")
H2v0j10=np.loadtxt("./wavefunctions/H2v0J10_norm.txt")


D2v0j0=np.loadtxt("./wavefunctions/D2v0J0_norm.txt")
D2v0j1=np.loadtxt("./wavefunctions/D2v0J1_norm.txt")
D2v0j2=np.loadtxt("./wavefunctions/D2v0J2_norm.txt")
D2v0j3=np.loadtxt("./wavefunctions/D2v0J3_norm.txt")
D2v0j4=np.loadtxt("./wavefunctions/D2v0J4_norm.txt")
D2v0j5=np.loadtxt("./wavefunctions/D2v0J5_norm.txt")
D2v0j6=np.loadtxt("./wavefunctions/D2v0J6_norm.txt")
D2v0j7=np.loadtxt("./wavefunctions/D2v0J7_norm.txt")
D2v0j8=np.loadtxt("./wavefunctions/D2v0J8_norm.txt")
D2v0j9=np.loadtxt("./wavefunctions/D2v0J9_norm.txt")
D2v0j10=np.loadtxt("./wavefunctions/D2v0J10_norm.txt")



rwave=np.loadtxt("./wavefunctions/r_wave.txt")
#----------------------------------------

#----------------------------------------
# function to interpolate a 2D wave along an axis (along columns)
def interpolate2Dw(input2D,omega,omega_final, OutputFileName):
    inputSize=input2D.shape
    print(inputSize)
    num_omega_final=(omega_final.shape)
    SizeOmegaFinal=(num_omega_final[0])
    outputArray=np.zeros((SizeOmegaFinal,inputSize[1]))
    for i in range(0,inputSize[1]):
        col=input2D[:,i]
        splc=interpolate.splrep(omega,col,s=0)
        interp=interpolate.splev(omega_final,splc,der=0)
        outputArray[:,i]=interp

    np.savetxt(OutputFileName,outputArray)

#----------------------------------------


# Use the interpolate2Dw function to interpolate the isotropy or anisotropy
# the interpolated output is saved as txt file.
interpolate2Dw( isotropy,omega_nm,omega_final,"isotropy_final.txt")
interpolate2Dw( anisotropy,omega_nm,omega_final,"anisotropy_final.txt")

isotropy_F=np.loadtxt("isotropy_final.txt")
isotropy_FT=isotropy_F.T        # transpose

anisotropy_F=np.loadtxt("anisotropy_final.txt")
anisotropy_FT=anisotropy_F.T    # transpose

# extract out polarizability along internuclear distance for certain wavelength
index=0 # index selects out certain column corresponding to certain wavelength.

omega0=omega_final[index]               # index governs the wavelength.
print("Selected wavelength:", omega0)

parameter1=isotropy_FT[:,index]
parameter2=anisotropy_FT[:,index]
#----------------------------------------
#----------------------------------------
#----------------------------------------
# function to the matrix element for certain rovibrational state
def MEcompute(psi1,psi2,psi_r, parameter, parameter_r ):

        # step 1: gen cubic spline coefs.
        splc=interpolate.splrep(parameter_r,parameter,s=0)

        # function defining the integrand which uses the spline coef array to give interpolated values
        def integrand(xpoint):
            return interpolate.splev(xpoint,splc,der=0)

        # generate interpolated parameter for same xaxis as psi
        parameter_interp=interpolate.splev(psi_r,splc,der=0)

        # compute the pointwise products
        p1=np.multiply(psi1,psi2)
        p2=np.multiply(p1,psi_r)
        p3=np.multiply(p2,psi_r)
        product=np.multiply(p3,parameter_interp)

        # step 1: gen cubic spline coefs
        splc=interpolate.splrep(psi_r,product,s=0)

        # compute the integral using adaptive Quadrature
        result=integrate.quadrature(integrand,0.2,4.48,tol=1.0e-9,maxiter=500)
        print("<psi1|parameter|psi2> = ",result)

#----------------------------------------
#******************************************************************************************


def integrand2(xpoint, splc):
    spline_array=splc
    result=interpolate.splev(xpoint,spline_array,der=0)
    return result


#----------------------------------------
# function to the matrix element for certain rovibrational state
def MEcompute2(psi1,psi2,psi_r, parameter, parameter_r ):

        # step 1: gen cubic spline coefs.
        splc=interpolate.splrep(parameter_r,parameter,s=0)

        # generate interpolated parameter for same xaxis as psi
        parameter_interp=interpolate.splev(psi_r,splc,der=0)

        # compute the pointwise products
        p1=np.multiply(psi1,psi2)
        p2=np.multiply(p1,psi_r)
        p3=np.multiply(p2,psi_r)
        product=np.multiply(p3,parameter_interp)

        # step 1: gen cubic spline coefs
        splc=interpolate.splrep(psi_r,product,s=0)

        # compute the integral using adaptive Quadrature
        result = integrate.quadrature(integrand2, 0.2, 4.48,tol=1.0e-9, maxiter=500,args=(splc,))
        print("<psi1|parameter|psi2> = ",result)

#----------------------------------------


#----------------------------------------
# computing the matrix element

# isotropy, 182.260 nm ; for H2v0j0
MEcompute(H2v0j0,H2v0j0,rwave,parameter1,distance)

# anisotropy, 182.260 nm ; for H2v0j0
MEcompute2(H2v0j0,H2v0j0,rwave,parameter1,distance)
#----------------------------------------
