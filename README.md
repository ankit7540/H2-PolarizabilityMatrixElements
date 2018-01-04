# H2-PolarizabilityMatrixElements
Set of data and python programs for interpolation (of wavelength dependent polarizability) and computation of the matrix elements(for rovibrational states) within the ground electronic state.

This repository contains :
 - Internuclear distance dependent polarizability for molecualr hydrogen
 - Rovibrational wavefunctions for H2, HD and D2 for v=0--4 and J=0--10.
 - A python module which can be used to compute the wavelength dependent matrix elements. Wavelength range available is 182.5 to 1320.6 nm.
 
 **Requirements**
  - python3
  - numpy
  - scipy
 
 **To use the program (in a UNIX environment or MS-DOS environment) follow the steps.**
 1. Download the zip file and unzip it to a folder.
 2. Move to the unzipped folder on the console.
 3. Initialize python by `python3`
 4. Import the `sys` module and add the current folder to path allowing to import the module in the current folder.
    > import sys
    
    > sys.path.append("..")
     
 5. Import the presnet module.
    > import ME-module
 6. If all requirements are met the following output should be produced.
    ```
    Dimensions of isotropy matrix : (44, 176)
    Dimensions of anisotropy matrix : (44, 176)
    Available wavelength range:  0.25 - 0.0345 Hartree;  182.253410111 - 1320.67688486  nm;  1822.53410111 - 13206.7688486  Angstrom.
    ...ready.
    ```
 7. Use the following command to do computation of the matrix element.
    > ME-module.MEcompute(mol, vl, Jl, vr, Jr, wavelength, wavelength_unit, operator):
      where the parameters are described below: 
      
    - mol  =    molecule specification (for H2 enter "H2", for D2 enter "D2", for HD enter "HD")
    - vl   =    vibrational state for the bra, vl = [0,4]
    - Jl   =    rotational state for the bra,  Jl = [0,10]
    - vr   =    vibrational state for the ket, vr = [0,4]
    - Jr   =    rotational state for the bra,  Jr = [0,10]
    - wavelength =  wavelength ( can be Hartree, nanometers or Angstrom)
    - wavelength_unit = specify unit using the specifier
                                ( for  Hartree           use "H" or "h"  )
                                ( for  nanometers        use "n" or "nm"  )
                                ( for  Angstrom          use "a" or "A"  )

    - operator   = property namely alpha_xx, alpha_zz, mean polarizability
                                   (isotropy)[\bar{alpha}], anisotropy[\gamma]
                                   Specify operator using the specifier.
                                 ( for  alpha_xx          use "x"     or  "xx"  )
                                 ( for  alpha_zz          use "z"     or  "zz"  )
                                 ( for  isotropy          use "iso"   or  "mp" or "mean" )
                                 ( for  anisotropy        use "aniso" or  "g"  or "diff" )
                                 ( for  all the above     use "all"   or  "ALL" )

 *Example*
 For H2,the following matrix element would need,
 \langle \psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0} \rangle
 ![f1]
 
 
 ```ME-module.MEcompute("H2",0,0,0,0,720.26,"n","diff")``` 
 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle \psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0} \rangle
