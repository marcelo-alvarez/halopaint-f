module haloproject

  use cosmology

  use mpivars
  use flatskyvars
  use healpixvars
  use maptable

  use diffuse_profiles
  use line_profile
  use halomodels

  use pix_tools
  use healpix_types
  use head_fits

  implicit none

  real, parameter :: rmax = 4.0
  double precision theta_halo,phi_halo,dtheta_halo,dtheta_max
  real angle,angle0,dangle

  ! For integration
  real tot_integrated,tot_pixelated
  integer na

  ! Halo angular properties
  real maxtheta,maxthetamass,maxthetachih  
  real maxthetal,maxthetamassl,maxthetachihl  

  ! Map moments
  double precision fsmom1,fsmom2,fsmom3
  double precision crmom1,crmom2,crmom3
  double precision hpmom1,hpmom2,hpmom3
  real curval

  real interpval, sfri
  ! Weird stuff
  real localpi

contains

!=======================================================================

  subroutine projecthalo_healpix(x,y,z,r,chi,m,redshift,v,mf)

    implicit none

    real x,y,z,r,chi,m,redshift,v,mf,interpval0
    integer i,j,profile

    localpi = 4*atan(1.0)

    dtheta_halo = asin(r/chi)
    dtheta_max  = dp_rmax * dtheta_halo

    ! Get unit vector    
    hpvec1(1) = x / chi 
    hpvec1(2) = y / chi
    hpvec1(3) = z / chi   

    ! Get integrated profile 
    tot_integrated = interpolate_table_norm(redshift,m)

    ! Profile paste
    call query_disc(hpnside,hpvec1,dtheta_max,hplist,nlist)
    if(nlist==0) then
       call vec2pix_ring(hpnside,hpvec1,j)
       hpmapl(j,1)=hpmapl(j,1) + tot_integrated / hpdomega * mf * v
       return
    endif
       
    tot_pixelated=0
    do i=0,nlist-1

       j=hplist(i)
       call pix2vec_ring(hpnside,j,hpvec2)
       call angdist(hpvec1,hpvec2,hpangle)
       
       angle = hpangle

       if(j>hpnpix-1.or.j<0) write(*,*) 'error',i,j,hpnpix,nlist
       
       if(angle<dtheta_max) then
          interpval   = interpolate_table(angle ,redshift,m)
          curval      =  mf * v * interpval 
          meantaul    = meantaul + interpval
          tot_pixelated = tot_pixelated + interpval
       endif

    enddo
    tot_pixelated = tot_pixelated * hpdomega

    if(tot_pixelated>0) then
       do i=0,nlist-1

          j=hplist(i)
          call pix2vec_ring(hpnside,j,hpvec2)
          call angdist(hpvec1,hpvec2,hpangle)
          
          angle = hpangle

          if(j>hpnpix-1.or.j<0) write(*,*) 'error',i,j,hpnpix,nlist
       
          if(angle<dtheta_max) then
             interpval   = interpolate_table(angle ,redshift,m)
             curval      =  mf * v * interpval 
             meantaul    = meantaul + interpval
             hpmapl(j,1) = hpmapl(j,1) + curval * tot_integrated / tot_pixelated
          endif

       enddo
    else
       j=hplist(0)
       hpmapl(j,1)=hpmapl(j,1) + tot_integrated / hpdomega * mf * v
    endif

    return

  end subroutine projecthalo_healpix

!=======================================================================

  subroutine projecthalo_flatsky(x,y,z,r,chi,m,redshift,v,mf)

    implicit none

    real x,y,z,r,chi,m,redshift,v,mf

    integer i,j,k,jmin,jmax,kmin,kmax


    localpi = 4*atan(1.0)

    dtheta_halo = dp_rmax*asin(r/chi)
    if (profile==itco) dtheta_halo = 0.

    ! Get unit vector
    x = x / chi
    y = y / chi
    z = z / chi

    if(z<0.0) return

    theta_halo = asin(x) 
    phi_halo   = asin(y)

    thetamin = theta_halo - 1.1 * dtheta_halo
    thetamax = theta_halo + 1.1 * dtheta_halo
    
    phimin   = phi_halo - 1.1 * dtheta_halo
    phimax   = phi_halo + 1.1 * dtheta_halo
    
    if(thetamax < -fov/2 .or. thetamin > fov/2) return
    if(  phimax < -fov/2 .or.   phimin > fov/2) return

    ! For Pasted Luminosities
    if ((profile==itco) .or. (profile==ithi)) then
       if(theta_halo < -fov/2 .or. theta_halo > fov/2) return
       if(  phi_halo < -fov/2 .or.   phi_halo > fov/2) return

       j = int((theta_halo + fov/2)/dpix) + 1
       k = int((phi_halo + fov/2)/dpix) + 1

       if (profile==itco) then
          sfri   = 10**SFR(log10(m),redshift)
          curval = L_CO(sfri,redshift)

       elseif (profile==ithi) then
          curval = MHtoMHI(m,redshift)
          curval = L_HI(curval)
       endif
       curval = I_line(curval,redshift,chi,dnu)
       curval = T_line(curval,nuobsj)
       fsmapl(j+npix*(k-1)) = fsmapl(j+npix*(k-1)) + curval
    else
    
    ! For Pasted Profiles
    if(thetamin <= -fov/2) thetamin = -fov/2 + 1e-10
    if(thetamax >=  fov/2) thetamax =  fov/2 - 1e-10
    
    if(  phimin <= -fov/2) phimin = -fov/2 + 1e-10
    if(  phimax >=  fov/2) phimax =  fov/2 - 1e-10
    
    jmin = max(1,int((thetamin + fov/2)/dpix)+1)
    jmax = min(int((thetamax + fov/2)/dpix)+1,npix)
    
    kmin = max(1,int((phimin + fov/2)/dpix)+1)
    kmax = min(int((phimax + fov/2)/dpix)+1,npix)
    
    if(dtheta_halo > maxtheta) then
       maxtheta = dtheta_halo
       maxthetachih = chi
       maxthetamass = m
    endif
    
    do k=kmin,kmax
       phip   = -fov/2 + (k-0.5)*dpix
       do j=jmin,jmax           
          thetap   = -fov/2 + (j-0.5)*dpix
          
          dphi   = phip   - phi_halo
          dtheta = thetap - theta_halo
          
          angle = sqrt(dphi**2+dtheta**2)
          
          if(angle<dtheta_halo) then
             curval =  mf * v *interpolate_table(angle,redshift,m)              
             fsmapl(j+npix*(k-1)) = fsmapl(j+npix*(k-1)) + curval
          endif
          
       enddo
    enddo
    endif
    return 
    
  end subroutine projecthalo_flatsky

!=======================================================================

  subroutine projecthalo_car(x,y,z,r,chi,m,redshift,v,mf)

    implicit none

    real x,y,z,r,chi,m,redshift,v,mf
    real lon,lat,dlon,dlat,lonmin,lonmax,latmin,latmax,xpix,ypix,zpix,rcyl,&
         lonpix,latpix
    real(kind=dp) thetacr,phicr,domega
    integer i,j,k,imin,imax,jmin,jmax


    localpi = 4*atan(1.0)

    domega = 4*localpi / npixc / npixc / 2.0

    dtheta_halo = dp_rmax*asin(r/chi)

    ! Get unit vector
    hpvec1(1) = x / chi
    hpvec1(2) = y / chi
    hpvec1(3) = z / chi

    ! Get integrated profile 
    tot_integrated = interpolate_table_norm(redshift,m)

    call vec2ang(hpvec1,thetacr,phicr)

    lon = phicr   - localpi
    lat = localpi / 2 - thetacr

    dlat = dtheta_halo

    dlon = min(localpi, dtheta_halo * cos(lat) / &
         min(cos(max(lat-dlat,-localpi/2)),cos(min(lat+dlat,localpi/2))))

    latmax = lat + dlat
    latmin = lat - dlat

    lonmin = lon - dlon
    lonmax = lon + dlon

    jmin = int((latmin + localpi/2) / dpixc) + 1
    jmax = int((latmax + localpi/2) / dpixc) + 1

    imin = int((lonmin + localpi  ) / dpixc) + 1
    imax = int((lonmax + localpi  ) / dpixc) + 1

    if(jmin<1)     jmin = 1
    if(jmax>npixc) jmax = npixc

    if(imin<1)       imin = 1
    if(imax>2*npixc) imax = 2*npixc

    tot_pixelated=0
    do i=imin,imax
       lonpix = -localpi + (i-0.5) * dpixc
       phicr  =  localpi + lonpix
       
       do j=jmin,jmax
          latpix  = -localpi/2 + (j-0.5) * dpixc
          thetacr = localpi/2 - latpix

          hpvec2(3) = cos(thetacr)
          rcyl      = sin(thetacr)

          hpvec2(1) = rcyl * cos(phicr)
          hpvec2(2) = rcyl * sin(phicr)       
          
          call angdist(hpvec1,hpvec2,hpangle)
          
          angle = hpangle
          
          if(angle<dtheta_halo) then
             curval =  mf * v *interpolate_table(angle,redshift,m)              
             tot_pixelated = tot_pixelated + curval 
          endif
          
       enddo
    enddo
    tot_pixelated = tot_pixelated * domega 

    do i=imin,imax
       lonpix = -localpi + (i-0.5) * dpixc
       phicr  =  localpi + lonpix
       
       do j=jmin,jmax
          latpix  = -localpi/2 + (j-0.5) * dpixc
          thetacr = localpi/2 - latpix

          hpvec2(3) = cos(thetacr)
          rcyl      = sin(thetacr)

          hpvec2(1) = rcyl * cos(phicr)
          hpvec2(2) = rcyl * sin(phicr)       
          
          call angdist(hpvec1,hpvec2,hpangle)
          
          angle = hpangle
          
          if(angle<dtheta_halo) then
             curval =  mf * v *interpolate_table(angle,redshift,m)              
             crmapl(npixc*(i-1)+j) = crmapl(npixc*(i-1)+j) + curval * tot_integrated / tot_pixelated
          endif
          
       enddo
    enddo

    if(tot_pixelated == 0) then
       i = int((lon + localpi  ) / dpixc) + 1
       j = int((lat + localpi/2) / dpixc) + 1
       crmapl(npixc*(i-1)+j) = crmapl(npixc*(i-1)+j) + tot_integrated / domega * mf * v
    endif

    return 
    
  end subroutine projecthalo_car

end module haloproject
