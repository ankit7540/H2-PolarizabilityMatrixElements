# H2-PolarizabilityMatrixElements
Set of data and python programs for interpolation (of wavelength dependent polarizability) and computation of the matrix elements(for rovibrational states) within the ground electronic state. The program evaluates the following integral.

![integral image][img0]


This repository contains :
 - Internuclear distance dependent polarizability for molecualr hydrogen (Omega in the above integral).
 - Rovibrational wavefunctions for H2, HD and D2 for v=0--4 and J=0--10.
 - A python module which can be used to compute the wavelength dependent matrix elements. Wavelength range available is 182.5 to 1320.6 nm.
 
**Requirements**

Local installation of : 
  - python2 or python3
  - numpy
  - scipy
 
====== 

**To use the program (in a UNIX environment or MS-DOS environment) follow the steps.**
1. Download the zip file and unzip it to a folder.
2. Move to the unzipped folder on the console.
3. Initialize python by `python3`
4. In the python console, import the `sys` module and add the current folder to path allowing to import the module in the current folder.
    > import sys
    
    > sys.path.append("..")
     
5. Import the `ME-module` which should be in your current folder.
    > import rovibME
6. If all requirements are met the following output should be produced.
    ```
    Dimensions of isotropy matrix : (44, 176)
    Dimensions of anisotropy matrix : (44, 176)
    Available wavelength range:  0.25 - 0.0345 Hartree;  182.253410111 - 1320.67688486  nm;  1822.53410111 - 13206.7688486  Angstrom.
    ...ready.
    ```
7. Use the following command to do computation of the matrix element.
    > rovibME.MEcompute(mol, vl, Jl, vr, Jr, wavelength, wavelength_unit, operator):
    
    where the parameters are described below: 
      
    - mol  =    molecule specification (for H2 enter "H2", for D2 enter "D2", for HD enter "HD")
    - vl   =    vibrational state for the bra, vl = [0,4]
    - Jl   =    rotational state for the bra,  Jl = [0,10]
    - vr   =    vibrational state for the ket, vr = [0,4]
    - Jr   =    rotational state for the ket,  Jr = [0,10]
    - wavelength =  wavelength ( can be Hartree, nanometers or Angstrom)
    - wavelength_unit = specify unit using the specifier, ( for  Hartree use "H" or "h" , for  nanometers use "n" or "nm" , for  Angstrom use "a" or "A"  )
    - operator   = property namely alpha_xx, alpha_zz, mean polarizability (isotropy) and anisotropy. Specify operator using the specifier. ( For  alpha_xx  use "x"     or  "xx" , for  alpha_zz  use "z"     or  "zz" , for  isotropy  use "iso"   or  "mp" or "mean" , for  anisotropy use "aniso" or  "g"  or "diff" and for  all the above  use "all"   or  "ALL" .

**Examples**

A few matrix elements and their corresponding commands are shown below,

- ![f1] for H2 
 
```rovibME.MEcompute("H2",0,0,0,0,488,"n","mp")``` 


- ![f2] for D2

```rovibME.MEcompute("D2",2,1,1,1,0.15,"H","g")``` 


- ![f3] for HD

```rovibME.MEcompute("HD",2,1,1,1,3550,"A","zz")``` 
 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0}\rangle
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\gamma|\psi_{v=1,J=1}\rangle
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{zz}|\psi_{v=1,J=1}\rangle

[img0]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/01-05-2018_82.png "Logo Title Text 2"
