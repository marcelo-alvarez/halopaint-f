module maptable

  use cosmology
  use diffuse_profiles
  use integrate_profiles
  use halomodels

  implicit none

  ! Table bounds and interval
  real   chit,   chimint,   chimaxt
  real    mht,    mhmint,    mhmaxt
  real     rt,     rmint,     rmaxt

  ! Table sizes and intervals
  integer nrt, nchit, nmht, nprofiles
  real    drt, dchit, dmht

  parameter (nprofiles=3)

  ! Table
  real, allocatable :: table(:,:,:,:),table_norm(:,:,:)

contains

  subroutine allocatetable

    allocate(table_norm(    nmht,nchit,nprofiles))
    allocate(table     (nrt,nmht,nchit,nprofiles))

    return

  end subroutine allocatetable

  !=======================================================================

  subroutine makemaptable(fileout)

    implicit none 

    character *128 fileout

    ! Integrated gas
    real sgas
    ! Integrated DM
    real sdm

    ! Other variables
    integer i,j,k,l,npos
    real    m,z,theta
    real    pfac,mfac,zfac,dp_rvir0,total,dtheta,thetap,da,mup,mu,rtp

    ! Set profile function, currently only gnfw available
    dp_function => gnfw
       
    ! Create distance-redshift tables
    call set_rofztable
    call set_zofrtable
  
    dp_delta = delta_vir(dp_modelp)
 
    y0 = fb * dp_delta * h**2 * ypermass  ! in 1/Msun

    ! Setup table bounds

    ! Radius
    nrt   = 100
    rmint = 1e-2
    rmaxt = 20
    
    drt = (log(rmaxt)-log(rmint))/(nrt-1)

    ! Mass
    nmht   = 500
    mhmint = 1e12
    mhmaxt = 1e16
    
    dmht = (log(mhmaxt)-log(mhmint))/(nmht-1)

    ! Distance
    nchit   = 500
    chimint = 1
    chimaxt = 8e3

    ! Ensure that observer is never inside integration region
    chimint = max(1.1*rmaxt,chimint)

    dchit = (log(chimaxt)-log(chimint))/(nchit-1)

    call allocatetable

    ! Calculate tabulated values and output
    
    open(unit=1,file=fileout,form='binary')
    write(1) nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt
    
    do i=3,1,-1

       dp_rmax = 1                   ! only go to Rvir by default
       pfac = 1
       if(i==1) then
          dp_rmax  =  4             ! for tSZ go to 4Rvir          
          dp_model =  dp_modelp
          pfac     =  y0
          profp    => pres_profp
       elseif(i==2) then
          dp_model =  dp_modeld
          profp    => rhog_profp
       elseif(i==3) then
          dp_model =  dp_modeld
          profp    => rhod_profp
       endif

       dp_delta = delta_vir(dp_model)

       do j=1,nchit
          write(*,*) i,j,nchit,nmht
          chit    = exp( log(chimint) + (j-1) * dchit )
          z       = zofr(chit)          
          dp_chi  = rofz(z)
          dp_chi  = chit

          zfac       = 1
          if(i==1) then
             zfac     = (omegam*(1+z)**3+omegal)
          elseif(dp_modeld<4) then
             zfac     = (1+omegal/omegam/(1+z)**3)
          endif
          dp_rvir0 = (3./4/pi/dp_delta/rhocrit(z))**(1./3.)

          do k=1,nmht
             m  = exp( log(mhmint) + (k-1) * dmht )

             dp_rvir = dp_rvir0 * m**(1./3.)
             if ((dp_modeld==6) .and. (i==3)) then
                xc    = gnfw_pofmz(m * (1.e14/(2.e12/h)) ,z,profp(dp_model)%    xc) ! is now xc = c200
                f0    = nfw_pofmz(xc, profp(dp_model)%    f0)
                xc    = 1./xc  !xc is now really xc
             else
                xc    = gnfw_pofmz(m,z,profp(dp_model)%    xc)
                f0    = gnfw_pofmz(m,z,profp(dp_model)%    f0)
             endif

             alpha = gnfw_pofmz(m,z,profp(dp_model)% alpha)
             beta  = gnfw_pofmz(m,z,profp(dp_model)%  beta)
             gamma = gnfw_pofmz(m,z,profp(dp_model)% gamma)

             mfac     = m
             if(i>1) then
                mfac  = 1
                beta  = (beta-gamma)/alpha
                gamma = -gamma
             endif
             
             total = 0 
             thetap = 0
             rtp = rmint
             npos = 0
             do l=1,nrt
                rt = exp( log(rmint) + (l-1) * drt )
                if(rt/chit<1e-2 .or. rt/chit>0.9) then
                   theta = rt/chit
                else
                   theta = asin(rt/chit)
                endif
                da = sin(theta) * theta * drt
                table(l,k,j,i) = pfac * mfac * zfac * integrate_profile(theta)
                total = total + da * table(l,k,j,i)
             enddo
             table_norm(k,j,i) = 2 * pi * total

          enddo
       enddo
    enddo

    write(1) table
    write(1) table_norm
    
    close(1)
    
  end subroutine makemaptable
  
  !=======================================================================

  subroutine loadmaptable(filein)

    implicit none

    character *128 filein

    open(unit=1,file=filein,form='binary')

    read(1) nchit,nmht,nrt,chimint,chimaxt,mhmint,mhmaxt,rmint,rmaxt

    call allocatetable

    read(1) table
    read(1) table_norm

    close(1)

    dchit = (log(chimaxt) - log(chimint)) / (nchit - 1)
    dmht  = ( log(mhmaxt) -  log(mhmint)) / (nmht - 1)
    drt   = (  log(rmaxt) -   log(rmint)) / (nrt - 1)

    return

  end subroutine loadmaptable

  !=======================================================================

  real function interpolate_table(theta,z,mh)

    implicit none

    real theta,z,chi,r,mh
    real fr,fc,fm
    integer ir,ichi,imh

    chi = rofz(z)
    r   = chi * sin(theta)

    interpolate_table = 0
    if(r>rmaxt .or. chi>chimaxt) return

    if( chi < chimint ) chi = chimint + 1e-5
    if(  mh < mhmint  )  mh =  mhmint + 1e-5
    if(   r < rmint   )   r =   rmint + 1e-5

    if(  mh > mhmaxt  )  mh =  mhmaxt - 1e-5

    ichi = int( ( log(chi) - log(chimint) ) / dchit ) + 1
    imh  = int( (  log(mh) -  log(mhmint) ) / dmht  ) + 1
    ir   = int( (   log(r) -   log(rmint) ) / drt   ) + 1
    
    fc = ( log(chi) - ( log(chimint) + (ichi - 1) * dchit ) ) / dchit
    fm = ( log(mh)  - (  log(mhmint) + ( imh - 1) *  dmht ) ) / dmht
    fr = ( log(r)   - (   log(rmint) + (  ir - 1) *   drt ) ) / drt

    if(fc<0) fc=0
    if(fc>1) fc=1
    if(fm<0) fm=0
    if(fm>1) fm=1
    if(fr<0) fr=0
    if(fr>1) fr=1
    interpolate_table = &
         table(ir  ,imh  ,ichi  ,profile) * (1-fr) * (1-fm) * (1-fc) + &
         table(ir+1,imh  ,ichi  ,profile) * (  fr) * (1-fm) * (1-fc) + &
         table(ir  ,imh+1,ichi  ,profile) * (1-fr) * (  fm) * (1-fc) + &
         table(ir+1,imh+1,ichi  ,profile) * (  fr) * (  fm) * (1-fc) + &
         table(ir  ,imh  ,ichi+1,profile) * (1-fr) * (1-fm) * (  fc) + &
         table(ir+1,imh  ,ichi+1,profile) * (  fr) * (1-fm) * (  fc) + &
         table(ir  ,imh+1,ichi+1,profile) * (1-fr) * (  fm) * (  fc) + &
         table(ir+1,imh+1,ichi+1,profile) * (  fr) * (  fm) * (  fc)
    
    if(interpolate_table<0) interpolate_table = 0

    !write(*,*) interpolate_table,ir,imh,ichi,fc,fm,fr,chi,mh,r,z,theta

  end function interpolate_table

  !=======================================================================

  real function interpolate_table_norm(z,mh)

    implicit none

    real theta,z,chi,r,mh
    real fr,fc,fm
    integer ir,ichi,imh

    chi = rofz(z)

    interpolate_table_norm = 0
    if(chi>chimaxt) return

    if( chi < chimint ) chi = chimint + 1e-5
    if(  mh < mhmint  )  mh =  mhmint + 1e-5

    if(  mh > mhmaxt  )  mh =  mhmaxt - 1e-5

    ichi = int( ( log(chi) - log(chimint) ) / dchit ) + 1
    imh  = int( (  log(mh) -  log(mhmint) ) / dmht  ) + 1

    fc = ( log(chi) - ( log(chimint) + (ichi - 1) * dchit ) ) / dchit
    fm = ( log(mh)  - (  log(mhmint) + ( imh - 1) *  dmht ) ) / dmht
    
    if(fc<0) fc=0
    if(fc>1) fc=1
    if(fm<0) fm=0
    if(fm>1) fm=1
    interpolate_table_norm = &
         table_norm(imh  ,ichi  ,profile) * (1-fm) * (1-fc) + &
         table_norm(imh+1,ichi  ,profile) * (  fm) * (1-fc) + &
         table_norm(imh  ,ichi+1,profile) * (1-fm) * (  fc) + &
         table_norm(imh+1,ichi+1,profile) * (  fm) * (  fc)
                
  end function interpolate_table_norm

end module maptable
