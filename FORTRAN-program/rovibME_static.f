	program rovibme

	implicit none
	integer maxviblev,maxrotlev,numdist,numfreq
	parameter(maxviblev=4,maxrotlev=15,numdist=176,numfreq=44)
	integer ngrid
	parameter(ngrid=1075)
	integer i,j,k,l,brav,braJ,ketv,ketJ
	integer neval,ier
	real*8 wf(3,maxviblev+1,maxrotlev+1,ngrid)
	real*8 gridval(ngrid),grand(ngrid),grand2(ngrid)
	real*8 rmin,rmax,r,lambda,der1,der2
	real*8 intstart,intstop,epsabs,epsrel,result,abserr
	parameter(rmin=0.2d0,rmax=4.48d0)
	real*8 dist(numdist),freq(numfreq)
        real*8 static_xx(numdist),static_zz(numdist)
        real*8 alpha_mean(numdist),gamma(numdist)
C	real*8 alphapar(numfreq,numdist),alphaperp(numfreq,numdist)
	real*8 alpha(numdist),alpha0(numfreq),alpha2(numfreq)
	real*8 omval(numdist),omval2(numdist)
	character*2 elname(3),molname
	character*16 filename,omega,unitname 
	data elname /'H2','D2','HD'/


C	read static_xx data
C	############################################################
	open(20,file="data/static_xx.txt")
	do i=1,numdist
	  read(20,*)static_xx(i)
	end do
	close(20)


C       read static_zz data
C       ############################################################
        open(20,file="data/static_zz.txt")
        do i=1,numdist
          read(20,*)static_zz(i)
        end do
        close(20)


C	read distances at which polarizability was computed
C	###################################################
	open(20,file="data/distance.txt")
	do i=1,numdist
	  read(20,*)dist(i)
	end do
	close(20)


C	read the distance vector at which the rovibrational wfs are tabulated
C	#####################################################################
	open(20,file="wavefunctions/r_wave.txt")
	do i=1,ngrid
	  read(20,*)gridval(i)
	end do
	close(20)

C	read rovibrational wave functions
C	#################################
	do i=1,3
	  write(filename(1:3),'(A2,A1)')elname(i),'v'
	  do j=0,maxviblev
   	    write(filename(4:5),'(I1,A1)')j,'J'
	    do k=0,maxrotlev
 	      if (k .lt. 10) then
   	        write(filename(6:15),'(I1,A)')k,'_norm.txt'
	        filename=filename(1:15)
	      else
   	        write(filename(6:16),'(I2,A)')k,'_norm.txt'
	        filename=filename(1:16)
	      end if
	      open(20,file="wavefunctions/"//filename)
	        do l=1,ngrid
		  read(20,*)wf(i,j+1,k+1,l)
	        end do
	      close(20)
	    end do
	  end do
	end do

C	ask user for the input parameters
C	#################################
	write(*,*)'COMPUTE STATIC POLARIZABILITY MATRIX ELEMENTS'
	write(*,*)'Give input parameters: molecule name, v_bra, J_bra,',
     *            ' v_ket, J_ket, Omega'
	write(*,*)'(example: H2 1 0 2 3 alpha_mean)'
	write(*,*)'   '
	write(*,*)'info: molecule name should be H2, D2, or HD'
	write(*,*)'      v should be in the interval [0..4]'
	write(*,*)'      J should be in the interval [0..15]'
	write(*,*)'      possible values of Omega are: ',
     *		  'alpha_par, alpha_perp, alpha_mean, gamma'
	write(*,*)'   '
	write(*,*)'now give your input:'
	write(*,*)'   '
	read(*,*)molname,brav,braJ,ketv,ketJ,omega

C	sanity check of the input parameters
C	####################################
	if (trim(molname) .ne. 'H2' .and. trim(molname) .ne. 'D2' 
     *      .and. trim(molname) .ne. 'HD') then
	    print*,'error exit: wrong molecule name'
	    stop
        end if
	if (trim(omega) .ne. 'alpha_par' 
     *      .and. trim(omega) .ne. 'alpha_perp' 
     *      .and. trim(omega) .ne. 'alpha_mean' 
     *      .and. trim(omega) .ne. 'gamma') then
	    print*,'error exit: wrong omega operator name'
	    stop
        end if
	if (braJ .lt. 0 .or. braJ .gt. maxrotlev) then
	    print*,'error exit: wrong value of J_bra'
	    stop
        end if
	if (ketJ .lt. 0 .or. ketJ .gt. maxrotlev) then
	    print*,'error exit: wrong value of J_ket'
	    stop
        end if
	if (brav .lt. 0 .or. brav .gt. maxrotlev) then
	    print*,'error exit: wrong value of v_bra'
	    stop
        end if
	if (ketv .lt. 0 .or. ketv .gt. maxrotlev) then
	    print*,'error exit: wrong value of v_ket'
	    stop
        end if
	

C	select the correct quantity for interpolation
C	#############################################
	if (trim(omega) .eq. 'alpha_par') then
	do j=1,numdist
	    alpha(j)=static_zz(j)
	 end do
	else if (trim(omega) .eq. 'alpha_perp') then 
	do j=1,numdist
	    alpha(j)=static_xx(j)
	end do
	else if (trim(omega) .eq. 'alpha_mean') then 
	do j=1,numdist
	    alpha(j)=(2*static_xx(j)+static_zz(j))/3.0d0
        end do
	else if (trim(omega) .eq. 'gamma') then
	do j=1,numdist
	    alpha(j)=static_zz(j)-static_xx(j)
	end do
	end if

C	find the value of omega at lambda for each distance
C	###################################################
C	do j=1,numdist
C	  do i=1,numfreq
C	    alpha0(i)=alpha(i,j)
C	  end do
C	  call deriv1a(freq,alpha0,der1)
C	  call deriv1b(freq(numfreq-3),alpha0(numfreq-3),der2)
C	  call spline(freq,alpha0,numfreq,der1,der2,alpha2)
C	  call splint(freq,alpha0,alpha2,numfreq,lambda,omval(j))
C	end do

C	find the spline parameters for omega at lambda
C	##############################################
	call deriv1a(dist,alpha,der1)
	call deriv1b(dist(numdist-3),alpha(numdist-3),der2)
	call spline(dist,alpha,numdist,der1,der2,omval2)

C	compute the value of the integrand on the wf grid
C	#################################################
	if (trim(molname) .eq. 'H2') then
	  i=1
	else if (trim(molname) .eq. 'D2') then 
	  i=2
	else if (trim(molname) .eq. 'HD') then
	  i=3
	end if
	do l=1,ngrid
	  call splint(dist,alpha,omval2,numdist,gridval(l),grand(l))
	  grand(l)=grand(l)*wf(i,brav+1,braJ+1,l)
	  grand(l)=grand(l)*wf(i,ketv+1,ketJ+1,l)
	  grand(l)=grand(l)*gridval(l)**2
	end do

C	find the spline parameters for the integrand
C	############################################
	call deriv1a(gridval,grand,der1)
	call deriv1b(gridval(ngrid-3),grand(ngrid-3),der2)
	call spline(gridval,grand,ngrid,der1,der2,grand2)

C	compute the value of the integral
C	#################################
	intstart=gridval(1)
	intstop=4.48d0
	epsabs=1.0d-05
	epsrel=1.0d-05
        call dqng(gridval,grand,grand2,ngrid,intstart,
     *            intstop,epsabs,epsrel,result,abserr,neval,ier)

	write(*,*) " "
	write(*,*) " "
	write(*,*) "      value      error     np_in_quadrature"
	write(*,*) "      #####      #####     ################" 
	write(*,'(f14.8,d10.2,i13)') result,abserr,neval


	end 

C	********************************************************
C	********************************************************


	SUBROUTINE spline(x,y,n,yp1,ypn,y2)
C	fit cubic spline taken from NUMERICAL RECIPIES IN FORTRAN 77
	implicit none
	INTEGER n,NMAX
	REAL*8 yp1,ypn,x(n),y(n),y2(n)
	PARAMETER (NMAX=1100)
	INTEGER i,k
	REAL*8 p,qn,sig,un,u(NMAX)
	if (yp1.gt..99e30) then 
	  y2(1)=0.0d0
	  u(1)=0.0d0
	else
	  y2(1)=-0.5d0
	  u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
	endif
	do i=2,n-1 
	  sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	  p=sig*y2(i-1)+2.0d0
	  y2(i)=(sig-1.0d0)/p
	  u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *          /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
	end do
	if (ypn.gt..99e30) then
	  qn=0.0d0
	  un=0.0d0
	else
	  qn=0.5d0
	  un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
	endif
	y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
	do k=n-1,1,-1
	  y2(k)=y2(k)*y2(k+1)+u(k)
	end do
	return
	END

        SUBROUTINE splint(xa,ya,y2a,n,x,y)
C	compute value of cubic spline taken from NUMERICAL RECIPIES IN FORTRAN 77
	implicit none
        INTEGER n
        REAL*8 x,y,xa(n),y2a(n),ya(n)
        INTEGER k,khi,klo
        REAL*8 a,b,h
        klo=1 
        khi=n
1       if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if (xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
          goto 1
        endif  
        h=xa(khi)-xa(klo)
        if (h.eq.0.d0) stop 'bad xa input in splint' 
          a=(xa(khi)-x)/h
          b=(x-xa(klo))/h
          y=a*ya(klo)+b*ya(khi)+
     *       ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0d0
        return
        END


	subroutine deriv1a(x,y,der)
C	compute the value of first derivative at x[1]
C	for a cubic fit to (x[1],y[1]),..,(x[4],y[4])
	implicit none
	real*8 x(4),y(4),der
	der=(x(2)*x(4)+x(3)*x(4)+x(2)*x(3)-2*x(1)*x(2)-2*x(1)*x(4)
     *     -2*x(1)*x(3)+3*x(1)**2)/(-x(4)+x(1))/(-x(2)+x(1))/(-x(3)
     *     +x(1))*y(1)-(-x(3)+x(1))*(-x(4)+x(1))/(x(2)-x(4))/(-x(2)
     *     +x(1))/(x(2)-x(3))*y(2)+(-x(2)+x(1))*(-x(4)+x(1))/(x(3)
     *     -x(4))/(x(2)-x(3))/(-x(3)+x(1))*y(3)-(-x(2)+x(1))*(-x(3)
     *     +x(1))/(x(3)-x(4))/(x(2)-x(4))/(-x(4)+x(1))*y(4)
	return
	end 


	subroutine deriv1b(x,y,der)
C	compute the value of first derivative at x[4]
C	for a cubic fit to (x[1],y[1]),..,(x[4],y[4])
	implicit none
	real*8 x(4),y(4),der
	der=(x(2)-x(4))*(x(3)-x(4))/(-x(4)+x(1))/(-x(2)+x(1))/(-x(3)
     *      +x(1))*y(1)-(-x(4)+x(1))*(x(3)-x(4))/(x(2)-x(4))/(-x(2)
     *      +x(1))/(x(2)-x(3))*y(2)+(-x(4)+x(1))*(x(2)-x(4))/(x(3)
     *      -x(4))/(x(2)-x(3))/(-x(3)+x(1))*y(3)-(-2*x(1)*x(4)
     *      -2*x(2)*x(4)-2*x(3)*x(4)+x(1)*x(2)+x(1)*x(3)+x(2)*x(3)
     *      +3*x(4)**2)/(x(3)-x(4))/(x(2)-x(4))/(-x(4)+x(1))*y(4)
	return
	end 


      subroutine dqng(gridval,grand,grand2,ngrid,a,b,epsabs,
     *                epsrel,result,abserr,neval,ier)
C	taken from https://github.com/scipy/scipy/blob/master/scipy/integrate/quadpack/dqng.f
c***begin prologue  dqng
c***date written   800101   (yymmdd)
c***revision date  810101   (yymmdd)
c***category no.  h2a1a1
c***keywords  automatic integrator, smooth integrand,
c             non-adaptive, gauss-kronrod(patterson)
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl math & progr. div. - k.u.leuven
c           kahaner,david,nbs - modified (2/82)
c***purpose  the routine calculates an approximation result to a
c            given definite integral i = integral of f over (a,b),
c            hopefully satisfying following claim for accuracy
c            abs(i-result).le.max(epsabs,epsrel*abs(i)).
c***description
c
c non-adaptive integration
c standard fortran subroutine
c double precision version
c
c           f      - double precision
c                    function subprogram defining the integrand function
c                    f(x). the actual name for f needs to be declared
c                    e x t e r n a l in the driver program.
c
c           a      - double precision
c                    lower limit of integration
c
c           b      - double precision
c                    upper limit of integration
c
c           epsabs - double precision
c                    absolute accuracy requested
c           epsrel - double precision
c                    relative accuracy requested
c                    if  epsabs.le.0
c                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
c                    the routine will end with ier = 6.
c
c         on return
c           result - double precision
c                    approximation to the integral i
c                    result is obtained by applying the 21-point
c                    gauss-kronrod rule (res21) obtained by optimal
c                    addition of abscissae to the 10-point gauss rule
c                    (res10), or by applying the 43-point rule (res43)
c                    obtained by optimal addition of abscissae to the
c                    21-point gauss-kronrod rule, or by applying the
c                    87-point rule (res87) obtained by optimal addition
c                    of abscissae to the 43-point rule.
c
c           abserr - double precision
c                    estimate of the modulus of the absolute error,
c                    which should equal or exceed abs(i-result)
c
c           neval  - integer
c                    number of integrand evaluations
c
c           ier    - ier = 0 normal and reliable termination of the
c                            routine. it is assumed that the requested
c                            accuracy has been achieved.
c                    ier.gt.0 abnormal termination of the routine. it is
c                            assumed that the requested accuracy has
c                            not been achieved.
c           error messages
c                    ier = 1 the maximum number of steps has been
c                            executed. the integral is probably too
c                            difficult to be calculated by dqng.
c                        = 6 the input is invalid, because
c                            epsabs.le.0 and
c                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
c                            result, abserr and neval are set to zero.
c
c***references  (none)
c***routines called  d1mach,xerror
c***end prologue  dqng
c
      double precision a,absc,abserr,b,centr,dabs,dhlgth,dmax1,dmin1,
     *  d1mach,epmach,epsabs,epsrel,f,fcentr,fval,fval1,fval2,fv1,fv2,
     *  fv3,fv4,hlgth,result,res10,res21,res43,res87,resabs,resasc,
     *  reskh,savfun,uflow,w10,w21a,w21b,w43a,w43b,w87a,w87b,x1,x2,x3,x4
      integer ier,ipx,k,l,neval
	integer ngrid
	real*8 gridval(ngrid),grand(ngrid),grand2(ngrid),myval
	real*8 fval1a,fval2a
C      external f
c
      dimension fv1(5),fv2(5),fv3(5),fv4(5),x1(5),x2(5),x3(11),x4(22),
     *  w10(5),w21a(5),w21b(6),w43a(10),w43b(12),w87a(21),w87b(23),
     *  savfun(21)
c
c           the following data statements contain the
c           abscissae and weights of the integration rules used.
c
c           x1      abscissae common to the 10-, 21-, 43- and 87-
c                   point rule
c           x2      abscissae common to the 21-, 43- and 87-point rule
c           x3      abscissae common to the 43- and 87-point rule
c           x4      abscissae of the 87-point rule
c           w10     weights of the 10-point formula
c           w21a    weights of the 21-point formula for abscissae x1
c           w21b    weights of the 21-point formula for abscissae x2
c           w43a    weights of the 43-point formula for abscissae x1, x3
c           w43b    weights of the 43-point formula for abscissae x3
c           w87a    weights of the 87-point formula for abscissae x1,
c                   x2, x3
c           w87b    weights of the 87-point formula for abscissae x4
c
c
c gauss-kronrod-patterson quadrature coefficients for use in
c quadpack routine qng.  these coefficients were calculated with
c 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
c
      data x1    (  1) / 0.9739065285 1717172007 7964012084 452 d0 /
      data x1    (  2) / 0.8650633666 8898451073 2096688423 493 d0 /
      data x1    (  3) / 0.6794095682 9902440623 4327365114 874 d0 /
      data x1    (  4) / 0.4333953941 2924719079 9265943165 784 d0 /
      data x1    (  5) / 0.1488743389 8163121088 4826001129 720 d0 /
      data w10   (  1) / 0.0666713443 0868813759 3568809893 332 d0 /
      data w10   (  2) / 0.1494513491 5058059314 5776339657 697 d0 /
      data w10   (  3) / 0.2190863625 1598204399 5534934228 163 d0 /
      data w10   (  4) / 0.2692667193 0999635509 1226921569 469 d0 /
      data w10   (  5) / 0.2955242247 1475287017 3892994651 338 d0 /
c
      data x2    (  1) / 0.9956571630 2580808073 5527280689 003 d0 /
      data x2    (  2) / 0.9301574913 5570822600 1207180059 508 d0 /
      data x2    (  3) / 0.7808177265 8641689706 3717578345 042 d0 /
      data x2    (  4) / 0.5627571346 6860468333 9000099272 694 d0 /
      data x2    (  5) / 0.2943928627 0146019813 1126603103 866 d0 /
      data w21a  (  1) / 0.0325581623 0796472747 8818972459 390 d0 /
      data w21a  (  2) / 0.0750396748 1091995276 7043140916 190 d0 /
      data w21a  (  3) / 0.1093871588 0229764189 9210590325 805 d0 /
      data w21a  (  4) / 0.1347092173 1147332592 8054001771 707 d0 /
      data w21a  (  5) / 0.1477391049 0133849137 4841515972 068 d0 /
      data w21b  (  1) / 0.0116946388 6737187427 8064396062 192 d0 /
      data w21b  (  2) / 0.0547558965 7435199603 1381300244 580 d0 /
      data w21b  (  3) / 0.0931254545 8369760553 5065465083 366 d0 /
      data w21b  (  4) / 0.1234919762 6206585107 7958109831 074 d0 /
      data w21b  (  5) / 0.1427759385 7706008079 7094273138 717 d0 /
      data w21b  (  6) / 0.1494455540 0291690566 4936468389 821 d0 /
c
      data x3    (  1) / 0.9993333609 0193208139 4099323919 911 d0 /
      data x3    (  2) / 0.9874334029 0808886979 5961478381 209 d0 /
      data x3    (  3) / 0.9548079348 1426629925 7919200290 473 d0 /
      data x3    (  4) / 0.9001486957 4832829362 5099494069 092 d0 /
      data x3    (  5) / 0.8251983149 8311415084 7066732588 520 d0 /
      data x3    (  6) / 0.7321483889 8930498261 2354848755 461 d0 /
      data x3    (  7) / 0.6228479705 3772523864 1159120344 323 d0 /
      data x3    (  8) / 0.4994795740 7105649995 2214885499 755 d0 /
      data x3    (  9) / 0.3649016613 4658076804 3989548502 644 d0 /
      data x3    ( 10) / 0.2222549197 7660129649 8260928066 212 d0 /
      data x3    ( 11) / 0.0746506174 6138332204 3914435796 506 d0 /
      data w43a  (  1) / 0.0162967342 8966656492 4281974617 663 d0 /
      data w43a  (  2) / 0.0375228761 2086950146 1613795898 115 d0 /
      data w43a  (  3) / 0.0546949020 5825544214 7212685465 005 d0 /
      data w43a  (  4) / 0.0673554146 0947808607 5553166302 174 d0 /
      data w43a  (  5) / 0.0738701996 3239395343 2140695251 367 d0 /
      data w43a  (  6) / 0.0057685560 5976979618 4184327908 655 d0 /
      data w43a  (  7) / 0.0273718905 9324884208 1276069289 151 d0 /
      data w43a  (  8) / 0.0465608269 1042883074 3339154433 824 d0 /
      data w43a  (  9) / 0.0617449952 0144256449 6240336030 883 d0 /
      data w43a  ( 10) / 0.0713872672 6869339776 8559114425 516 d0 /
      data w43b  (  1) / 0.0018444776 4021241410 0389106552 965 d0 /
      data w43b  (  2) / 0.0107986895 8589165174 0465406741 293 d0 /
      data w43b  (  3) / 0.0218953638 6779542810 2523123075 149 d0 /
      data w43b  (  4) / 0.0325974639 7534568944 3882222526 137 d0 /
      data w43b  (  5) / 0.0421631379 3519181184 7627924327 955 d0 /
      data w43b  (  6) / 0.0507419396 0018457778 0189020092 084 d0 /
      data w43b  (  7) / 0.0583793955 4261924837 5475369330 206 d0 /
      data w43b  (  8) / 0.0647464049 5144588554 4689259517 511 d0 /
      data w43b  (  9) / 0.0695661979 1235648452 8633315038 405 d0 /
      data w43b  ( 10) / 0.0728244414 7183320815 0939535192 842 d0 /
      data w43b  ( 11) / 0.0745077510 1417511827 3571813842 889 d0 /
      data w43b  ( 12) / 0.0747221475 1740300559 4425168280 423 d0 /
c
      data x4    (  1) / 0.9999029772 6272923449 0529830591 582 d0 /
      data x4    (  2) / 0.9979898959 8667874542 7496322365 960 d0 /
      data x4    (  3) / 0.9921754978 6068722280 8523352251 425 d0 /
      data x4    (  4) / 0.9813581635 7271277357 1916941623 894 d0 /
      data x4    (  5) / 0.9650576238 5838461912 8284110607 926 d0 /
      data x4    (  6) / 0.9431676131 3367059681 6416634507 426 d0 /
      data x4    (  7) / 0.9158064146 8550720959 1826430720 050 d0 /
      data x4    (  8) / 0.8832216577 7131650137 2117548744 163 d0 /
      data x4    (  9) / 0.8457107484 6241566660 5902011504 855 d0 /
      data x4    ( 10) / 0.8035576580 3523098278 8739474980 964 d0 /
      data x4    ( 11) / 0.7570057306 8549555832 8942793432 020 d0 /
      data x4    ( 12) / 0.7062732097 8732181982 4094274740 840 d0 /
      data x4    ( 13) / 0.6515894665 0117792253 4422205016 736 d0 /
      data x4    ( 14) / 0.5932233740 5796108887 5273770349 144 d0 /
      data x4    ( 15) / 0.5314936059 7083193228 5268948562 671 d0 /
      data x4    ( 16) / 0.4667636230 4202284487 1966781659 270 d0 /
      data x4    ( 17) / 0.3994248478 5921880473 2101665817 923 d0 /
      data x4    ( 18) / 0.3298748771 0618828826 5053371824 597 d0 /
      data x4    ( 19) / 0.2585035592 0216155180 2280975429 025 d0 /
      data x4    ( 20) / 0.1856953965 6834665201 5917141167 606 d0 /
      data x4    ( 21) / 0.1118422131 7990746817 2398359241 362 d0 /
      data x4    ( 22) / 0.0373521233 9461987081 4998165437 704 d0 /
      data w87a  (  1) / 0.0081483773 8414917290 0002878448 190 d0 /
      data w87a  (  2) / 0.0187614382 0156282224 3935059003 794 d0 /
      data w87a  (  3) / 0.0273474510 5005228616 1582829741 283 d0 /
      data w87a  (  4) / 0.0336777073 1163793004 6581056957 588 d0 /
      data w87a  (  5) / 0.0369350998 2042790761 4589586742 499 d0 /
      data w87a  (  6) / 0.0028848724 3021153050 1334156248 695 d0 /
      data w87a  (  7) / 0.0136859460 2271270188 8950035273 128 d0 /
      data w87a  (  8) / 0.0232804135 0288831112 3409291030 404 d0 /
      data w87a  (  9) / 0.0308724976 1171335867 5466394126 442 d0 /
      data w87a  ( 10) / 0.0356936336 3941877071 9351355457 044 d0 /
      data w87a  ( 11) / 0.0009152833 4520224136 0843392549 948 d0 /
      data w87a  ( 12) / 0.0053992802 1930047136 7738743391 053 d0 /
      data w87a  ( 13) / 0.0109476796 0111893113 4327826856 808 d0 /
      data w87a  ( 14) / 0.0162987316 9678733526 2665703223 280 d0 /
      data w87a  ( 15) / 0.0210815688 8920383511 2433060188 190 d0 /
      data w87a  ( 16) / 0.0253709697 6925382724 3467999831 710 d0 /
      data w87a  ( 17) / 0.0291896977 5647575250 1446154084 920 d0 /
      data w87a  ( 18) / 0.0323732024 6720278968 5788194889 595 d0 /
      data w87a  ( 19) / 0.0347830989 5036514275 0781997949 596 d0 /
      data w87a  ( 20) / 0.0364122207 3135178756 2801163687 577 d0 /
      data w87a  ( 21) / 0.0372538755 0304770853 9592001191 226 d0 /
      data w87b  (  1) / 0.0002741455 6376207235 0016527092 881 d0 /
      data w87b  (  2) / 0.0018071241 5505794294 8341311753 254 d0 /
      data w87b  (  3) / 0.0040968692 8275916486 4458070683 480 d0 /
      data w87b  (  4) / 0.0067582900 5184737869 9816577897 424 d0 /
      data w87b  (  5) / 0.0095499576 7220164653 6053581325 377 d0 /
      data w87b  (  6) / 0.0123294476 5224485369 4626639963 780 d0 /
      data w87b  (  7) / 0.0150104473 4638895237 6697286041 943 d0 /
      data w87b  (  8) / 0.0175489679 8624319109 9665352925 900 d0 /
      data w87b  (  9) / 0.0199380377 8644088820 2278192730 714 d0 /
      data w87b  ( 10) / 0.0221949359 6101228679 6332102959 499 d0 /
      data w87b  ( 11) / 0.0243391471 2600080547 0360647041 454 d0 /
      data w87b  ( 12) / 0.0263745054 1483920724 1503786552 615 d0 /
      data w87b  ( 13) / 0.0282869107 8877120065 9968002987 960 d0 /
      data w87b  ( 14) / 0.0300525811 2809269532 2521110347 341 d0 /
      data w87b  ( 15) / 0.0316467513 7143992940 4586051078 883 d0 /
      data w87b  ( 16) / 0.0330504134 1997850329 0785944862 689 d0 /
      data w87b  ( 17) / 0.0342550997 0422606178 7082821046 821 d0 /
      data w87b  ( 18) / 0.0352624126 6015668103 3782717998 428 d0 /
      data w87b  ( 19) / 0.0360769896 2288870118 5500318003 895 d0 /
      data w87b  ( 20) / 0.0366986044 9845609449 8018047441 094 d0 /
      data w87b  ( 21) / 0.0371205492 6983257611 4119958413 599 d0 /
      data w87b  ( 22) / 0.0373342287 5193504032 1235449094 698 d0 /
      data w87b  ( 23) / 0.0373610737 6267902341 0321241766 599 d0 /
c
c           list of major variables
c           -----------------------
c
c           centr  - mid point of the integration interval
c           hlgth  - half-length of the integration interval
c           fcentr - function value at mid point
c           absc   - abscissa
c           fval   - function value
c           savfun - array of function values which have already been
c                    computed
c           res10  - 10-point gauss result
c           res21  - 21-point kronrod result
c           res43  - 43-point result
c           res87  - 87-point result
c           resabs - approximation to the integral of abs(f)
c           resasc - approximation to the integral of abs(f-i/(b-a))
c
c           machine dependent constants
c           ---------------------------
c
c           epmach is the largest relative spacing.
c           uflow is the smallest positive magnitude.
c
c***first executable statement  dqng
      epmach = d1mach(4)
      uflow = d1mach(1)
c
c           test on validity of parameters
c           ------------------------------
c
      result = 0.0d+00
      abserr = 0.0d+00
      neval = 0
      ier = 6
      if(epsabs.le.0.0d+00.and.epsrel.lt.dmax1(0.5d+02*epmach,0.5d-28))
     *  go to 80
      hlgth = 0.5d+00*(b-a)
      dhlgth = dabs(hlgth)
      centr = 0.5d+00*(b+a)
	call splint(gridval,grand,grand2,ngrid,centr,fcentr)
C      fcentr = f(centr)
      neval = 21
      ier = 1
c
c          compute the integral using the 10- and 21-point formula.
c
      do 70 l = 1,3
      go to (5,25,45),l
    5 res10 = 0.0d+00
      res21 = w21b(6)*fcentr
      resabs = w21b(6)*dabs(fcentr)
      do 10 k=1,5
        absc = hlgth*x1(k)
 	call splint(gridval,grand,grand2,ngrid,centr+absc,fval1)
	call splint(gridval,grand,grand2,ngrid,centr-absc,fval2)
C       fval1 = f(centr+absc)
C       fval2 = f(centr-absc)
        fval = fval1+fval2
        res10 = res10+w10(k)*fval
        res21 = res21+w21a(k)*fval
        resabs = resabs+w21a(k)*(dabs(fval1)+dabs(fval2))
        savfun(k) = fval
        fv1(k) = fval1
        fv2(k) = fval2
   10 continue
      ipx = 5
      do 15 k=1,5
        ipx = ipx+1
        absc = hlgth*x2(k)
 	call splint(gridval,grand,grand2,ngrid,centr+absc,fval1)
	call splint(gridval,grand,grand2,ngrid,centr-absc,fval2)
C        fval1 = f(centr+absc)
C        fval2 = f(centr-absc)
        fval = fval1+fval2
        res21 = res21+w21b(k)*fval
        resabs = resabs+w21b(k)*(dabs(fval1)+dabs(fval2))
        savfun(ipx) = fval
        fv3(k) = fval1
        fv4(k) = fval2
   15 continue
c
c          test for convergence.
c
      result = res21*hlgth
      resabs = resabs*dhlgth
      reskh = 0.5d+00*res21
      resasc = w21b(6)*dabs(fcentr-reskh)
      do 20 k = 1,5
        resasc = resasc+w21a(k)*(dabs(fv1(k)-reskh)+dabs(fv2(k)-reskh))
     *                  +w21b(k)*(dabs(fv3(k)-reskh)+dabs(fv4(k)-reskh))
   20 continue
      abserr = dabs((res21-res10)*hlgth)
      resasc = resasc*dhlgth
      go to 65
c
c          compute the integral using the 43-point formula.
c
   25 res43 = w43b(12)*fcentr
      neval = 43
      do 30 k=1,10
        res43 = res43+savfun(k)*w43a(k)
   30 continue
      do 40 k=1,11
        ipx = ipx+1
        absc = hlgth*x3(k)
  	call splint(gridval,grand,grand2,ngrid,centr+absc,fval1a)
	call splint(gridval,grand,grand2,ngrid,centr-absc,fval2a)
	fval=fval1a+fval2a
C       fval = f(absc+centr)+f(centr-absc)
        res43 = res43+fval*w43b(k)
        savfun(ipx) = fval
   40 continue
c
c          test for convergence.
c
      result = res43*hlgth
      abserr = dabs((res43-res21)*hlgth)
      go to 65
c
c          compute the integral using the 87-point formula.
c
   45 res87 = w87b(23)*fcentr
      neval = 87
      do 50 k=1,21
        res87 = res87+savfun(k)*w87a(k)
   50 continue
      do 60 k=1,22
        absc = hlgth*x4(k)
  	call splint(gridval,grand,grand2,ngrid,centr+absc,fval1a)
	call splint(gridval,grand,grand2,ngrid,centr-absc,fval2a)
        res87 = res87+w87b(k)*(fval1a+fval2a)
C        res87 = res87+w87b(k)*(f(absc+centr)+f(centr-absc))
   60 continue
      result = res87*hlgth
      abserr = dabs((res87-res43)*hlgth)
   65 if(resasc.ne.0.0d+00.and.abserr.ne.0.0d+00)
     *  abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
      if (resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1
     *  ((epmach*0.5d+02)*resabs,abserr)
      if (abserr.le.dmax1(epsabs,epsrel*dabs(result))) ier = 0
c ***jump out of do-loop
      if (ier.eq.0) go to 999
   70 continue
   80 call xerror('abnormal return from dqng ',26,ier,0)
  999 return
      end


      DOUBLE PRECISION FUNCTION D1MACH(I)
C	taken from https://github.com/scipy/scipy/blob/2526df72e5d4ca8bad6e2f4b3cbdfbc33e805865/scipy/integrate/mach/d1mach.f

      INTEGER I
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C  D1MACH( 5) = LOG10(B)
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
C  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
C  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
C  MANY MACHINES YET.
C  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
C  ON THE NEXT LINE
      DATA SC/0/
C  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
C  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
C          mail netlib@research.bell-labs.com
C          send old1mach from blas
C  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS.
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532
     *       .AND. SMALL(1) .EQ. -448790528) THEN
*           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935
     *       .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943
     *       .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684
     *       .AND. SMALL(2) .EQ. -448790528) THEN
*           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074
     *       .AND. SMALL(2) .EQ. 58688) THEN
*           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
*                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
*    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
* /* Standard C source for D1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*double d1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return DBL_MIN;
*	  case 2: return DBL_MAX;
*	  case 3: return DBL_EPSILON/FLT_RADIX;
*	  case 4: return DBL_EPSILON;
*	  case 5: return log10(FLT_RADIX);
*	  }
*	fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
*	exit(1); return 0; /* some compilers demand return values */
*}
      END
      SUBROUTINE I1MCRY(A, A1, B, C, D)
**** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END


      SUBROUTINE XERROR(MESS,NMESS,L1,L2)
C	taken from https://github.com/scipy/scipy/blob/2526df72e5d4ca8bad6e2f4b3cbdfbc33e805865/scipy/integrate/mach/xerror.f
C
C     THIS IS A DUMMY XERROR ROUTINE TO PRINT ERROR MESSAGES WITH NMESS
C     CHARACTERS. L1 AND L2 ARE DUMMY PARAMETERS TO MAKE THIS CALL
C     COMPATIBLE WITH THE SLATEC XERROR ROUTINE. THIS IS A FORTRAN 77
C     ROUTINE.
C
      CHARACTER*(*) MESS
      NN=NMESS/70
      NR=NMESS-70*NN
      IF(NR.NE.0) NN=NN+1
      K=1
      PRINT 900
  900 FORMAT(/)
      DO 10 I=1,NN
        KMIN=MIN0(K+69,NMESS)
        PRINT *, MESS(K:KMIN)
        K=K+70
   10 CONTINUE
      PRINT 900
      RETURN
      END
