# H<sub>2</sub>-PolarizabilityMatrixElements
Set of data and python programs for interpolation (of wavelength dependent polarizability) and computation of the matrix elements (for rovibrational states) within the ground electronic state. The program evaluates the following integral.

![integral image][img0]

This repository contains :
 - Internuclear distance dependent polarizability  ,
 
Property | Symbol
------------ | -------------
Polarizability perpendicular to the internuclear axis | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" width="30" height="15" />

Polarizability parallel to the internuclear axis | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png" width="23" height="20" />

Mean polarizability | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" width="195" height="28" />

Polarizability anisotropy | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" width="155" height="28" />



c
<img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" width="30" height="15" />
Content in the first column 
 
 <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" width="30" height="15" />
 
  <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png" width="23" height="20" />
  
   <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" width="195" height="28" />
   
    <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" width="155" height="28" />
 
 (![alpha_{perp}][img1])  (![alpha_{parallel}][img2])(![alpha_{mp}][img3])(![gamma][img4])   ( <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha_{\perp}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha_{\perp}" title="\alpha_{\perp}" /></a> , <a href="https://www.codecogs.com/eqnedit.php?latex=\alpha_{\parallel}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\alpha_{\parallel}" title="\alpha_{\parallel}" /></a>, 
 <a href="https://www.codecogs.com/eqnedit.php?latex=\bar{\alpha}=(2\alpha_{\perp}&plus;\alpha_{\parallel})/3" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\bar{\alpha}=(2\alpha_{\perp}&plus;\alpha_{\parallel})/3" title="\bar{\alpha}=(2\alpha_{\perp}+\alpha_{\parallel})/3" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=\gamma=\alpha_{\parallel}-\alpha_{\perp}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma=\alpha_{\parallel}-\alpha_{\perp}" title="\gamma=\alpha_{\parallel}-\alpha_{\perp}" /></a>) for molecular hydrogen (Omega in the above integral).
 - Rovibrational wavefunctions for H<sub>2</sub>, HD and D<sub>2</sub> for v=0--4 and J=0--10.
 - A python module which can be used to compute the wavelength dependent matrix elements. Wavelength range available is 182.5 to 1320.6 nm.
 
**Requirements**
---

Local installation of : 
  - python2 or python3 with numpy and scipy modules
 
 

**Using the program**
---
1. Download the zip file and unzip it to a folder, or simply clone the repository if using `git`by executing

``` git clone https://github.com/ankit7540/H2-PolarizabilityMatrixElements.git ```

2. Move to the unzipped folder on the console/terminal.
3. In the python console, import the `sys` module and add the current folder to path allowing to import the module in the current folder.
    > import sys
    
    > sys.path.append("..")
     
4. Import the `rovibME` which should be in your current folder.
    > import rovibME
5. If all requirements are met the following output should be produced.
    ```
    Dimensions of isotropy matrix : (44, 176)
    Dimensions of anisotropy matrix : (44, 176)
    Available wavelength range: 0.25 - 0.0345 Hartree; 182.2534 - 1320.6769 nm; 1822.5341 - 13206.7688 Angstrom
    ...ready.
    ```
6. Use the following command to do computation of the matrix element.
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
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{zz}|\psi_{v=1,J=1}\rangle

[img0]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/01-05-2018_82.png "Logo Title Text 2"
[img1]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png "Logo alpha_{perp}"
[img2]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png "Logo alpha_{paralell}"
[img3]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png "Logo alpha_{mp}"
[img4]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png "Logo alpha_{aniso}"
