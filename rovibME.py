# Purpose : Load the matrix containing polarizability and wavefunctions, interpolate the
# polarizability if needed, and compute the respective matrix elements.
# Author : Ankit

# Load necessary modules
import sys
import numpy as np
import time
from scipy import interpolate
from scipy import integrate

#********************************************************************
# Dimension check on polarizability data :

# Load the polarizability data ( alpha_xx and alpha_zz)
alpha_xx=np.loadtxt("matrix_xxf.txt")
alpha_zz=np.loadtxt("matrix_zzf.txt")
omega=np.loadtxt("freq.txt")
distance=np.loadtxt("distance.txt")

# check size of the arrays
print("Dimensions of isotropy matrix :",alpha_xx.shape)
print("Dimensions of anisotropy matrix :",alpha_zz.shape)
#print(alpha_xx.shape[0], alpha_xx.shape[1])
if not(alpha_xx.shape ==  alpha_zz.shape  or len(omega)==alpha_xx.shape[0] ):
    print("Dimension check on polarizability data matrices or wavelength file failed.")
    quit()
else :
    omega_nm=(1e7/(omega*219474.6313702000))
    omega_A=(omega_nm*10)
    print("Available wavelength range: ",round(omega[0],4),"-",round(omega[-1],4),"Hartree; ",round(omega_nm[0],4),"-",round(omega_nm[-1],4)," nm; ",round(omega_A[0],4),"-",round(omega_A[-1],4)," Angstrom.")
    print("Available wavelength range: {0} - {1} Hartree; {2} - {3} nm; {4} - {5} Angstrom".format(round(omega[0],4),round(omega[-1],4),round(omega_nm[0],4),round(omega_nm[-1],4),round(omega_A[0],4),round(omega_A[-1],4)   ))
    print("...ready.")
#********************************************************************

# function to the matrix element for certain rovibrational state
def compute(mol, vl, Jl, vr, Jr, wavelength, wavelength_unit, operator):
    '''#  parameters:
    # mol  =    molecule (for H2 enter "H2", for D2 enter "D2", for HD enter "HD")
    # vl   =    vibrational state for the bra, vl = [0,4]
    # Jl   =    rotational state for the bra,  Jl = [0,10]
    # vr   =    vibrational state for the ket, vr = [0,4]
    # Jr   =    rotational state for the ket,  Jr = [0,10]
    # wavelength =  wavelength ( can be Hartree, nanometers or Angstrom)
    # wavelength_unit = specify unit using the specifier
                                ( for  Hartree           use "H" or "h"  )
                                ( for  nanometers        use "n" or "nm"  )
                                ( for  Angstrom          use "a" or "A"  )

    # operator   = property namely alpha_xx, alpha_zz, mean polarizability
                                   (isotropy)[\bar{alpha}], anisotropy[\gamma]
                                   Specify operator using the specifier.
                                 ( for  alpha_xx          use "x"     or  "xx"  )
                                 ( for  alpha_zz          use "z"     or  "zz"  )
                                 ( for  isotropy          use "iso"   or  "mp" or "mean" )
                                 ( for  anisotropy        use "aniso" or  "g"  or "diff" )
                                 ( for  all the above     use "all"   or  "All"or "ALL"  )

    This function runs on both Python 2.7x and 3.x
    '''

    # set a dictionary for output array
    d={'output':[0]}
    #----------------------------------------------------------------
    # interpolation function defined here and is used later.
    def interpolate2D_common(input2D,originalx,finalx):
        inputSize=input2D.shape
        tempx=np.zeros((len(originalx),1))
        tempx[:,0]=originalx
        col=np.zeros((len(originalx),1))
        outputArray=np.zeros((1,inputSize[1]))
        for i in range(0,inputSize[1]):
            col[:,0]=input2D[:,i]
            spl=interpolate.splrep(tempx,col,k=3,s=0)
            interp=interpolate.splev(finalx,spl,der=0)
            outputArray[:,i]=interp
        d['output']=outputArray.T
    #----------------------------------------------------------------

    # Load the polarizability data ( alpha_xx and alpha_zz)
    alpha_xx=np.loadtxt("matrix_xxf.txt")
    alpha_zz=np.loadtxt("matrix_zzf.txt")
    omega=np.loadtxt("freq.txt")
    dist=np.loadtxt("distance.txt")
    distance= np.asarray(dist)
    omega_nm=(1e7/(omega*219474.6313702000)) # convert the original freq to nm

    # compute the isotropy(mean polarizability) and anisotropy (gamma)
    isotropy=np.absolute(2*(np.array(alpha_xx))+np.array(alpha_zz))/3
    anisotropy=np.absolute(np.array(alpha_zz)-np.array(alpha_xx))

    # step 1: load the required wavefunctions ------------------------
    Wfn1="./wavefunctions/{0}v{1}J{2}_norm.txt".format(mol,vl,Jl)
    Wfn2="./wavefunctions/{0}v{1}J{2}_norm.txt".format(mol,vr,Jr)
    r_wave="./wavefunctions/r_wave.txt"
    #print(Wfn1,Wfn2)
    if vl < 0 or vr < 0 or vl > 4 or vr > 4 :
        print("v value out of range. vl and vr = [0,4].")
        quit()

    if Jl < 0  or Jr < 0 or Jl > 10 or Jr > 10 : 
        print("J value out of range. Jl and Jr =[0,10].")
        quit()

    if not (mol == "H2"  or mol == "HD" or mol == "D2" ): 
        print("Incorrect molecule chosen. For H2 enter H2, for D2 enter D2, for HD enter HD. Use quotes." )
        quit()

    # Proceed to load wavefunctions.
    psi1=np.loadtxt(Wfn1)
    psi2=np.loadtxt(Wfn2)
    rwave=np.loadtxt(r_wave)
    #print(len(psi1),len(psi2),len(rwave))
    #----------------------------------------------------------------

    wv=float(wavelength) # entered wavelength is a float.

    if (wavelength_unit == "h" or wavelength_unit == "H"):
        omegaFinal=(1e7/(wv*219474.6313702000))
    elif (wavelength_unit == "n" or wavelength_unit == "nm"):
        omegaFinal=wv
    elif (wavelength_unit == "a" or wavelength_unit == "A"):
        omegaFinal=(wv/10)
    elif not (wavelength_unit == "h" or wavelength_unit == "H" or wavelength_unit == "n" or wavelength_unit == "nm"
              or wavelength_unit == "a" or wavelength_unit == "A" ):
        print("Default unit of nm will be used.")
        omegaFinal=wv

    print("Wavelength in nanometer : {0}".format(round(omegaFinal,6)))

    #print(omega_nm[0],omega_nm[-1])
    if omegaFinal < omega_nm[0] or omegaFinal > omega_nm[-1]:
        sys.exit("Requested wavelength is out of range.")
        #print("Requested wavelength is out of range.")
        #quit()
    #--------------------------------------------------------------
    n=0
    if (operator == "x" or operator == "xx" ):
        param=alpha_xx
        name=["alpha_xx"]
        n=1
    elif (operator == "z" or operator == "zz" ):
        param=alpha_zz
        name=["alpha_zz"]
        n=1
    elif (operator == "mean" or operator == "mp" or operator == "iso" ):
        param=isotropy
        name=["isotropy"]
        n=1
    elif (operator == "diff" or operator == "g" or operator == "aniso" ):
        param=anisotropy
        name=["anisotropy"]
        n=1
    elif (operator == "all" or operator == "All" or operator == "ALL") :
        list = [alpha_xx,alpha_zz,isotropy,anisotropy]
        name=["alpha_xx","alpha_zz","isotropy","anisotropy"]
        n=4
    else :
        print("Operator not correctly specified.")
        quit()

    for i in range(n):    # evaluation of  interpolation and integral
        #print(i)

        if not(n == 1):
            param=list[i]
            #print(param.shape)
            interpolate2D_common(param,omega_nm,omegaFinal)
        else :
            interpolate2D_common(param,omega_nm,omegaFinal)

        parameter=np.zeros((len(distance),1))
        temp=d['output'] 
        parameter[:,0]=temp[:,0] 

        # step 1: gen cubic spline coefs.
        splc=interpolate.splrep(distance,parameter,s=0)

        # step 2: generate interpolated parameter for same xaxis as psi
        parameter_interp=interpolate.splev(rwave,splc,der=0)

        # compute the pointwise products
        p1=np.multiply(psi1,psi2)
        p2=np.multiply(p1,rwave)
        p3=np.multiply(p2,rwave)
        product=np.multiply(p3,parameter_interp)

        # function defining the integrand which uses the spline coef array to give interpolated values
        def integrand(xpoint):
            result=interpolate.splev(xpoint,spline_array,der=0)
            return result

        # step 3: gen cubic spline coefs
        spline_array=interpolate.splrep(rwave,product,s=0)

        # compute the integral using adaptive Quadrature
        result=integrate.quadrature(integrand,0.2,4.48,tol=1.0e-9,maxiter=1000)
        print("{0} < v={1}J={2} | {3} | v={4}J={5} >  =  {6} a.u.". format(mol,vl,Jl,name[i],vr,Jr,abs(round(result[0],6)) ) )    

    

