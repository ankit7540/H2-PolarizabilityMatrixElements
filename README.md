# H<sub>2</sub>-PolarizabilityMatrixElements
Set of data and python programs for interpolation (of wavelength dependent polarizability) and computation of the matrix elements (for rovibrational states) within the ground electronic state. The program evaluates the following integral.

![integral image][img0]

This repository contains :
 - Internuclear distance dependent polarizability,
 
Property | Definition
------------ | -------------
Polarizability perpendicular to the internuclear axis | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png" width="30" height="15" />
Polarizability parallel to the internuclear axis | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png" width="23" height="20" />
Mean polarizability | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png" width="195" height="28" />
Polarizability anisotropy | <img src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" data-canonical-src="https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png" width="155" height="28" />

The above properties are available as Omega in the above integral for H<sub>2</sub> HD and D<sub>2</sub>.
 - Rovibrational wavefunctions for H<sub>2</sub>, HD and D<sub>2</sub> for v=0--4 and J=0--10.
 - A python module which can be used to compute the wavelength dependent matrix elements. Wavelength range available is 182.25 to 1320.6 nm.
 
**Available programs**
--- 
The programs for computation of matrix element are written in FORTRAN and Python. These are independent programs which do the same job.

In the case of FORTRAN, two different programs exist, *(i)* `rovibME_dynamic.f` for wavelength dependent matrix elements and *(ii)* `rovibME_static.f` for static ones.

In the case of Python, one program `rovibME.py` deals with both static and dynamic matrix elements.

**Usage**
---
Refer to the `README.md` in the FORTRAN-program folder and the Python-module folder respectively according to your usage.


**Comments of numerical accuracy**
---
The integral calculation is accurate usually to ~5e-7. The net numerical uncertainity in the computed matrix element however is  1e-4 which includes the uncertainities introduced by the accuracy of the wavefunctions, polarizability, spline procedures and physical constants. 

**Credits**
---
Cubic spline interpolation procedure used in FORTRAN and python codes has been adapted from Numerical Recipes in FORTRAN, William H. Press, Saul A. Teukolsky, William T. Vetterling, Brian P. Flannery, Michael Metcalf, Cambridge University Press; 2 edition (September 25, 1992).

Adaptive Gausssian Quadrature implemented in SciPy has been used.

FORTRAN code by Prof. Henryk A Witek (NCTU, Taiwan).

Python code by Ankit Raj.



[img0]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/01-05-2018_82.png "Logo Title Text 2"
[img1]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_perp.png "Logo alpha_{perp}"
[img2]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_parallel.png "Logo alpha_{paralell}"
[img3]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/alpha_mp.png "Logo alpha_{mp}"
[img4]: https://github.com/ankit7540/H2-PolarizabilityMatrixElements/blob/master/image/gamma.png "Logo alpha_{aniso}"
