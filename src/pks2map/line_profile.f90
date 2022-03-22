module line_profile

  use cosmology
  use flatskyvars
  implicit none
  real maxMhtab, minMhtab, maxZttab, minZttab 
  real, allocatable  :: SFRtab(:,:)
  real, allocatable  :: Ztitab(:), Zttab(:), Mhtab(:), Mhitab(:)
  integer, parameter :: lenZttab = 137
  integer, parameter :: lenMhtab = 122
  real, parameter    :: pi_local = 4*atan(1.0)
  
  !CO line emission parameters
  real nu_rest, dnu, nuobs_i, nuobs_f !rest freq, bandwidth, range
  real nuobsi, nuobsj  !nu of halo, nu of individual slice (matches nested do j= and do i= loops)
  integer num_nuobs
  
contains

  subroutine SFRtabload(filename)
    !Zt is (1+z) from data
    !Mh is halo mass from data
    implicit none
  
    integer i,j,index
    character *128 filename
   

    open(unit=1,file=filename)
    
    !allocate data length
    allocate(SFRtab(lenMhtab,lenZttab))
    allocate(Ztitab(lenMHtab*lenZttab))
    allocate(Mhitab(lenMHtab*lenZttab))
    allocate(Mhtab(lenMHtab))
    allocate(Zttab(lenZttab))
    !read table
    
    index=0
    do i=1,lenMhtab
       do j=1,lenZttab
          index=index+1
          read(1,*) Ztitab(index), Mhitab(index), SFRtab(i,j)
       enddo
    !Z=Z-1
    enddo
    
    close(1)
    !sorting Zt and Mh vectors to not repeat
    Zttab=Ztitab(1:lenZttab)
    
    do i=1,lenMhtab
      Mhtab(i)=Mhitab(i*lenZttab)
    enddo
    
    maxMhtab=maxval(Mhtab)
    minMhtab=minval(Mhtab)
    maxZttab=maxval(Zttab)
    minZttab=minval(Zttab)

!    write(*,*) "Min/Max Halo Mass for allowed SFR is,", minMhtab, Maxmhtab
!    write(*,*) "Min/Max Redshift for allowed SFR is,", minZttab-1, MaxZttab-1

    return

  end subroutine SFRtabload

!===========================================================================
 
  real function SFR(mhaloi,zi)
    implicit none
    real mhaloi,zi
    real dZt, dMh, z1, z2, mh1, mh2, A
    integer indz, indmhalo

    !----------------------------------------------------------------
    ! The Star formation rate in units of log10(Msun/yr)
    ! from Peter Behroozi et al. fit
    ! data file from http://www.peterbehroozi.com/data.html
    ! Bilinear Interpolated
    ! lenMH = 122
    ! lenZt = 137
    !-------------------------------------------------------    
    !find interval (Mh is already equally spaced in log, but mhalo is not)
    dMh = (maxMhtab - minMhtab) / (lenMhtab - 1)
    dZt = (log(maxZttab)-log(minZttab)) / (lenZttab - 1)

    if (mhaloi>maxMhtab) mhaloi=maxMhtab 
    !--------------------------------------------------------------
    !  Bilinear Interpolation
    !  sfr(z,mhalo) = 1 / [(z2-z1)(mh2-mh1)] * { sfr(z1,mh1) * (z2-z)(mh2-mhalo) +
    !  sfr(z2,mh1) * (z-z1)(mh2-mhalo) + sfr(z2,m1) * (z2-z)(mhalo-mh1) +
    !  sfr(z2,mh2) * (z-z1)(mhalo-mh2) }
    !
    !  so given z, mhalo, we want to find the index above and below it for both 
    !--------------------------------------------------------------------
    
    !finding the index number below the value 

    indz = int( ( log(zi+1) - log(minZttab) ) / dZt ) + 1
    indmhalo = int( ( mhaloi - minMhtab ) / dMh) + 1
    
    ! to find the values of z1,z2,m1,m2 
    z1 = Zttab(indz)
    z2 = Zttab(indz+1)
    mh1 = mhtab(indmhalo)
    mh2 = mhtab(indmhalo+1)
    
    A=1.0 / ( (z2-z1)*(mh2-mh1) ) 
   
    sfr =A * ( SFRtab(indmhalo, indz)*(z2-(zi+1))*(mh2-mhaloi) + &
         SFRtab(indmhalo, indz+1)*((zi+1)-z1)*(mh2-mhaloi) + &
         SFRtab(indmhalo+1, indz)*(z2-(zi+1))*(mhaloi-mh1) + &
        SFRtab(indmhalo+1, indz+1)*((zi+1)-z1)*(mhaloi-mh1) )
         
    return
  end function SFR

!============================================================
 

  real function L_CO_orig(sfr,z)
    implicit none
    real angle, sfr, z
    real chi, rh
    
    ! ----------------------------------------------------------------
    ! The CO Luminosity for a given halo is given L_CO in units of Lsun
    ! L_CO = 3.2e4*(sfr)**(3.0/5)  
    ! ----------------------------------------------------------------

    ! Luminosity per unit area=L_tot*(1/(4*pi*Rvir**2))
    ! Flux is F_CO=L_CO/(4*pi*d_L)**2
    ! Where d_L=d_p*(1+z)
    L_CO_orig = 3.2e4*(sfr)**(3.0/5)

   return

  end function L_CO_orig

  real function L_CO(sfr,z)
    implicit none
    real sfr, z
    real chi
    real, parameter :: dmf = 1.0, beta_co=-1.74, alpha_co=1.37, sig_sfr=0.3, sig_lco=0.3
    ! ----------------------------------------------------------------
    ! The CO Luminosity for a given halo is given L_CO in units of Lsun
    ! L_CO = 4.9e-5*(nu_co,rest/115.27GHz)^3 * exp( (log(10^10/dmf*sfr)-beta)/alpha)
    ! units of L_sun
    ! ----------------------------------------------------------------
    

    L_CO = 4.9e-5*(nu_rest/115.27)**3 * 10**( (log10(1.0e10/dmf*sfr)-beta_co)/alpha_co)

   return

  end function L_CO

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! Neutral Hydrogen functions
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  real function MHtoMHI(mass,redshift)
    implicit none
    real mass, redshift
    real MminHI, MmaxHI
    real, parameter :: f3 = 0.000391003842023
    ! ----------------------------------------------------------------
    ! M_HI(M,z) = f3 * M/(1+M/M_max(z)) , if M>M_min. else M_HI(z)=0
    ! M_max and M_min are found from circular velocity arguments (see 0908.3796v3)
    ! M_max -> v_circ = 200km/s , M_min -> v_circ = 30km/s
    ! where M = 10^10 M_sun * (v_circ/60km/s)^3 * ((1+z)/4)^(-3/2) 
    ! ----------------------------------------------------------------

    MminHI = 1.e10 * (30./60)**3 * ((1+redshift)/4)**(-3./2)
    MmaxHI = 1.e10 * (200./60)**3 * ((1+redshift)/4)**(-3./2)
    if (mass<MminHI) then 
       MHtoMHI = 0.
    else
       MHtoMHI = f3 * mass/(1+mass/MmaxHI)
    endif

   return

 end function MHtoMHI

  real function L_HI(HImass)
    implicit none
    real HImass
    real, parameter :: convfac = 6.215e-9

    ! ----------------------------------------------------------------
    ! L_HI = 3*A_10*h*nu_0/4/m_p * M_HI
    ! A10 = 2.85e-15 s^-1 , nu_0 = 1420.405751786 MHz  
    ! so L_HI = 6.215e-9 L_sun * M_HI[M_sun]
    ! ----------------------------------------------------------------

    L_HI = convfac * HImass

   return

 end function L_HI




!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! General line functions
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  real function I_line(Lline,z,chi,dnu)
    implicit none
    real Lline, z
    real chi, dnu
    ! ----------------------------------------------------------------
    ! The Intensity
    ! I_line = L_line/4/pi/D_L^2/dnu
    ! Where D_L=D_p*(1+z)
    ! units of L_sun/Mpc^2/Hz
    ! ----------------------------------------------------------------

    ! Luminosity per unit area per steridian = L_tot*(1/(4*pi*Rvir**2))/sr
     I_line = Lline/4/pi_local/chi**2/(1+z)**2/dnu/Ompix

   return

 end function I_line
 
  real function T_line(iline,nuobs)
    implicit none
    real iline, nuobs
    real, parameter :: convfac = 2.63083e-6
    ! ----------------------------------------------------------------
    ! The line Temperature in Rayleigh-Jeans limit
    ! T_line = c^2*I_line/2/kb/nuobs^2 
    ! units of [L_sun/Mpc^2/GHz] * [(km/s)^2 / (J/K) / (GHz) ^2] * 1/sr
    !   = [ 3.48e26 W/Mpc^2/GHz ] * [ 6.50966e21 s^2/K/kg ] 
    !   = 2.63083e-6 K
    ! ----------------------------------------------------------------

    T_line = 1./2*convfac*iline/nuobs**2

   return

 end function T_line


end module line_profile


