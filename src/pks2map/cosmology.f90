module cosmology

  implicit none

  real omegam,omegab,omegal,h,ns,sigma8,w,h0,fb ! h0 is in km/s/Mpc
  real rho_0, rhocrit_0, y0, tau0, tau_dm0, mue

  double precision pi
  parameter (pi = 3.14159265359)

  double precision sigmat,c,me,mp,msun, Mpc, kpc ! all in cgs
  double precision Mpc2km, H100, ckms, rho100    ! H100 is 100 km/s/Mpc in 1/s
  parameter (sigmat = 6.65246d-25  ,&
                  c = 2.99792458d10,&
                 me = 9.109383d-28 ,&
                mue = 1.136        ,&
                 mp = 1.6726e-24   ,&
               msun = 1.98847d33   ,&
             Mpc2km = 3.08567758d19,&
                Mpc = 3.08567758d24,&
                kpc = 3.08567758d21,&
               H100 = 1e2/Mpc2km   ,&
             rho100 = 1.88d-29   ,&
               ckms = c/1e5)

  double precision ypermass ! this is the quantity (in units of 1/Msun)
                            !   P_Delta * (P_e/P_th) * sigmat / (me*c^2) * R_Delta / M_Delta / h^2 / Delta
                            !     = (P_e/P_th) * G * sigmat / (me*c^2) * rhocrit / 2 / h^2
                            !     = 3 * sigmat * (P_e/P_th) * (100 km/sec/Mpc)^2 / (16*pi*me*c^2) 
                            ! where 
                            !   P_Delta = G * M_Delta * Delta * rhocrit * fb / (2*R_Delta)
                            !   rhocrit = 3 * H0^2 / (8*pi*G) = 3 * h^2 * (100 km/sec/Mpc) / (8*pi*G)
                            !   P_th = 1.932 * P_e 
                            ! therefore
                            !   y = ypermass * (fb*Delta*h^2) * M_Delta * int [(ds/Rvir) * (P_th/P_Delta)]
                            !   y = sigmat / (me*c^2) * int [ ds * P_e ] 
  parameter (ypermass = 3*sigmat*H100**2/16./pi/me/c**2*msun/1.932)

  double precision tau0persigma ! this is the dimensionless quantity 
                                ! sigmat/mu/mp*rho100^(2/3)
                                ! * (3/4/pi*msun)**(1./3)
  parameter (tau0persigma = sigmat/mue/mp * (3./4/pi*msun)**(1./3) * & 
       rho100**(2./3))
  
  double precision sigmatovermump ! this is sigmat/mu/mp
                                  ! in Mpc^2 / Msun
  parameter (sigmatovermump = sigmat/mue/mp * (msun / Mpc**2))


  double precision ylpermsun ! this is the quantity
                             ! 3*sigmat*ke*(100 km/sec/Mpc)^2/4/(me*c^2)
                             ! in units of 1/Msun
  parameter (ylpermsun = 3*sigmat*H100**2/4./me/c**2*msun/1.932)

  integer nzofrtable,nrofztable
  real zofrtable_min,zofrtable_max,dzofrtable
  real rofztable_min,rofztable_max,drofztable
  real, allocatable :: rofztable(:),zofrtable(:)

  logical testcase 

contains

  !=======================================================================

  subroutine set_rofztable
    
    implicit none
    
    real r,z,dz
    integer i,j
    
    rofztable_min = 0
    rofztable_max = 6.
    nrofztable = 10000
    drofztable = (rofztable_max-rofztable_min)/(nrofztable-1)

    dz = drofztable

    allocate(rofztable(nrofztable))

    rofztable(1)=0.

    do i=2,nrofztable
       z=(i-1)*dz+rofztable_min
       rofztable(i)=rofztable(i-1)+drdz(z)*dz
    enddo

    return

  end subroutine set_rofztable

  !=======================================================================  

  subroutine set_zofrtable
    
    implicit none
    
    real z,dr
    integer i,j
    
    zofrtable_min = 0.
    zofrtable_max = 1.3e4
    nzofrtable = 10000
    dzofrtable = (zofrtable_max-zofrtable_min)/(nzofrtable-1)

    dr=dzofrtable

    allocate(zofrtable(nzofrtable))

    zofrtable(1)=0.

    do i=2,nzofrtable
       z=zofrtable(i-1)
       zofrtable(i)=zofrtable(i-1)+dzdr(z)*dr
    enddo

    return

  end subroutine set_zofrtable

  !=======================================================================

  real function drdz(z)

    implicit none

    real z

    drdz = 1 / dzdr(z)

    return

  end function drdz

  !=======================================================================

  real function dzdr(z)

    implicit none

    real z

    dzdr = h*100/ckms*sqrt(omegam*(1.+z)**3+1.-omegam) ! dzdr is in 1/Mpc

    return
    
  end function dzdr

  !=======================================================================

  real function zofr(r)

    real r
    real fint
    integer itab

    itab = ((r-zofrtable_min)/dzofrtable)+1

    if(itab<1) then
       zofr = zofrtable(1)
       return
    elseif(itab>=nzofrtable) then
       zofr = zofrtable(nzofrtable)
       return
    endif

    fint=(r-real(itab-1)*dzofrtable-zofrtable_min)/dzofrtable

    zofr=(1-fint)*zofrtable(itab)+fint*zofrtable(itab+1)

  end function zofr

  !=======================================================================

  real function rofz(z)

    real z
    real fint
    integer itab

    itab = ((z-rofztable_min)/drofztable)+1

    if(itab<1) then
       rofz = rofztable(1)
       return
    elseif(itab>=nrofztable) then
       rofz = rofztable(nrofztable)
       return
    endif

    fint=(z-real(itab-1)*drofztable-rofztable_min)/drofztable

    rofz=(1-fint)*rofztable(itab)+fint*rofztable(itab+1)

  end function rofz

  !=======================================================================

  real function rofzslow(z)
    
    implicit none
    
    real r,z,dz,zcur
    integer i,j
    
    dz = 0.01
    r = 0.0
    do while(zcur<z)
       r    = r    + drdz(zcur)*dz
       zcur = zcur + dz
    enddo

    rofzslow = r

    return 

  end function rofzslow

  !=======================================================================  

  !=======================================================================

  real function rhocrit(z)

    ! This is the *comoving* critical density as a function of redshift

    real z

    rhocrit = rhocrit_0 * (omegam*(1+z)**3+omegal) / (1+z)**3

    return

  end function rhocrit

  !=======================================================================

  real function deltavirc(z)

    ! This is the post-spherical-collapse overdensity for an isolated 
    ! uniform sphere in virial equilibrium in units of the *critical*
    ! density; see Bryan & Norman (1997)

    real z,x,omegaz

    omegaz = omegam*(1+z)**3/(omegam*(1+z)**3+1-omegam)
    x = omegaz - 1

    deltavirc = 18*pi**2 + 82 * x - 39 * x**2

    return

  end function deltavirc

  !=======================================================================

  real function W_Kappa(z,chi,chist)

    implicit none

    ! Lensing Kernel

    real z, chi, chist

    W_Kappa = 3./2 * omegam * (h * 100/ckms)**2 * (1+z) * chi*(1-chi/chist)

    return

  end function W_Kappa

end module cosmology

