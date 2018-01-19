`rovibME` FORTRAN program
-----------------------------
A FORTRAN program for the computation of rovibrational matrix elements for for H<sub>2</sub>, HD and D<sub>2</sub> within v=0--4 and J=0--10.

Compiling the code
-----------------------------
Presently the program is tested with `gfortran`.

Usage
-----------------------------
1. Under the FORTRAN-program directory compile the code using `gfortran`, for example,

```
gfortram -o computeME rovibME
```

to generate the executable `computeME`

2. Next, execute `./computeME` which should produce the following message and input details.


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


**Examples**
---

A few matrix elements and their corresponding commands are shown below,

- ![f1] for H<sub>2</sub> 
 
```H2 0 0 0 0 488 nm alpha_mean``` 


- ![f2] for D<sub>2</sub>

```D2 2 1 1 1 0.15 Hartree  gamma )``` 


- ![f3] for HD

```HD 2 1 1 1 3550 A zz alpha_par )``` 
 
[f1]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=0,J=0}|\bar{\alpha}|\psi_{v=0,J=0}\rangle
[f2]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\gamma|\psi_{v=1,J=1}\rangle
[f3]: http://chart.apis.google.com/chart?cht=tx&chl=\langle\psi_{v=2,J=1}|\alpha_{\parallel}|\psi_{v=1,J=1}\rangle
