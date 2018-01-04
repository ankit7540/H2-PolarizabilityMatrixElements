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
    > Dimensions of isotropy matrix : (44, 176)
    > Dimensions of anisotropy matrix : (44, 176)
    > Available wavelength range:  0.25 - 0.0345 Hartree;  182.253410111 - 1320.67688486  nm;  1822.53410111 - 13206.7688486  Angstrom.
...ready.

