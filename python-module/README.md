Python-module
----------------
Python module `rovibME` for the computation of matrix elements.


Requirements
----------------
Python 2.x or Python 3.x with Numpy and Scipy modules

Usage
----------------
1. After cloning the respository and moving in the `python-module` directory, add the current folder to path allowing to import the module in the current folder.
    > import sys
    
    > sys.path.append("..")
     
2. Import the `rovibME` which should be in your current folder.
    > import rovibME
3. If all requirements are met the following output should be produced.
    ```
    Give  rovibME.compute  command with parameters:
        rovibME.compute(molecule, bra_v, bra_J, ket_v, ket_J, lambda, unit of lambda, operator)
         for example:  rovibME.compute("H2",0,2,0,4,488,"n","mp")
                       rovibME.compute("D2",1,0,1,0,"static","n","all")

                molecule = for H2 enter "H2", for D2 enter "D2", for HD enter "HD"
                bra_v    = vibrational state, v=[0,4]
                bra_J    = rotataional state, J=[0,10]
                ket_v    = vibrational state, v=[0,4]
                ket_J    = rotataional state, J=[0,10]
                lambda   = wavelength in Hartree, nm or Angstrom, for static specify "s" or "static" here
                unit of lambda =  for  Hartree           use "H" or "h"
                                  for  nanometers        use "n" or "nm"
                                  for  Angstrom          use "a" or "A"
                                  if static property is asked then this parameter can be any of the three
                Available wavelength range: 0.25 - 0.0345 Hartree;
                                            182.2534 - 1320.6769 nm;
                                            1822.5341 - 13206.7688 Angstrom
                operator        = alpha(perpendicular) = alpha_xx given by "xx" or "x"
                                  alpha(parallel) = alpha_xx given by "zz" or "z"
                                  isotropy or mean polarizability given by "iso" or "mp" or "mean"
                                  anisotropy or polarizability difference or gamma given by "aniso" or "g"  or "diff"
                                  for all the above use "all"   or  "All"or "ALL"
    ...ready.

    ```
4. Use the following command to do computation of the matrix element.
    > rovibME.compute(mol, vl, Jl, vr, Jr, wavelength, wavelength_unit, operator):
    
    where the parameters are described below: 
      
    - mol  =    molecule specification (for H<sub>2</sub> enter "H2", for D<sub>2</sub> enter "D2", for HD enter "HD")
    - vl   =    vibrational state for the bra, vl = [0,4]
    - Jl   =    rotational state for the bra,  Jl = [0,10]
    - vr   =    vibrational state for the ket, vr = [0,4]
    - Jr   =    rotational state for the ket,  Jr = [0,10]
    - wavelength =  wavelength within the specified range ( 0.25 - 0.0345 Hartree;  182.2534 - 1320.6768  nm;  1822.5341 - 13206.7688  Angstrom ). Specify unit accordingly in the next parameter.
    - wavelength_unit = specify unit using the specifier, ( for  Hartree use "H" or "h" , for  nanometers use "n" or "nm" , for  Angstrom use "a" or "A"  )
    - operator   = property namely alpha_xx, alpha_zz, mean polarizability (isotropy) and anisotropy. Specify operator using the specifier. ( For  alpha_xx  use "x"     or  "xx" , for  alpha_zz  use "z"     or  "zz" , for  isotropy  use "iso"   or  "mp" or "mean" , for  anisotropy use "aniso" or  "g"  or "diff" and for  all the above 4 properties  use "all"   or  "ALL" .
    

**Examples**
---

A few matrix elements and their corresponding commands are shown below,

- ![f1] for H<sub>2</sub> 
 
```rovibME.compute("H2",0,0,0,0,488,"n","mp")``` 


- ![f2] for D<sub>2</sub>

```rovibME.compute("D2",2,1,1,1,0.15,"H","g")``` 


- ![f3] for HD

```rovibME.compute("HD",2,1,1,1,3550,"A","zz")``` 
 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0}\rangle
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\gamma|\psi_{v=1,J=1}\rangle
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{\parallel}|\psi_{v=1,J=1}\rangle
