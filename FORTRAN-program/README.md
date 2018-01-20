`rovibME_dynamic` and `rovibME_static` FORTRAN programs
-----------------------------
Programs for the computation of rovibrational matrix elements for for H<sub>2</sub>, HD and D<sub>2</sub> within v=0--4 and J=0--10.

Compiling the code
-----------------------------
Program has been tested with `gfortran`

Usage
-----------------------------
1. Under the FORTRAN-program directory compile the code using `gfortran`, for example,

```
gfortram -o ComputeDynamic rovibME_dynamic.f
```

to generate the executable `computeME`

2. Next, execute `./ComputeDynamic` which should produce the following message and ask for further input.


```
 Give input parameters: molecule name, v_bra, J_bra, v_ket, J_ket, laser wavelength lambda, wavelength unit, Omega
 (example: H2 1 0 2 3 532 nm alpha_mean)

 info: molecule name should be H2, D2, or HD
       v should be in the interval [0..4]
       J should be in the interval [0..10]
       units for lambda are: Hartree, nm, or A
       lambda should be in the interval [0.0345..0.25] Hartree
       lambda should be in the interval [182.26..1320.0] nm
       lambda should be in the interval [1822.6..13200] A
       possible values of Omega are: alpha_par, alpha_perp, alpha_mean, gamma

 now give your input:
```

3. Type in the parameters and press enter. For example,
 `H2 1 0 2 3 532 nm alpha_mean`

For static polarizabilities, compile as shown above and use the executable. The following message should be produced asking for user input. 

```
 COMPUTE STATIC POLARIZABILITY MATRIX ELEMENTS
 Give input parameters: molecule name, v_bra, J_bra, v_ket, J_ket, Omega
 (example: H2 1 0 2 3 alpha_mean)

 info: molecule name should be H2, D2, or HD
       v should be in the interval [0..4]
       J should be in the interval [0..10]
       possible values of Omega are: alpha_par, alpha_perp, alpha_mean, gamma

 now give your input:
 ```
 

**Examples**
---

A few matrix elements and their corresponding commands are shown below,

- ![f1] for H<sub>2</sub> at 488 nm
 
```H2 0 0 0 0 488 nm alpha_mean``` 


- ![f2] for D<sub>2</sub> at 0.15 Hartree

```D2 2 1 1 1 0.15 Hartree gamma ``` 


- ![f3] for HD at 3550 Angstrom

```HD 2 1 1 1 3550 A alpha_par ``` 


- Static mean polarizability for D2 in the ground vibrational and rotational state
``` D2 0 0 0 0 alpha_mean```

 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0}\rangle
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\gamma|\psi_{v=1,J=1}\rangle
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{\parallel}|\psi_{v=1,J=1}\rangle
