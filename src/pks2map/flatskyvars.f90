module flatskyvars
  
  !-----------------------------------------------------------------------
  ! FLATSKY VARIABLES
  !-----------------------------------------------------------------------
  
  ! Map arrays
  real, allocatable :: fsmap(:), fsmapl(:)
  real, allocatable :: crmap(:), crmapl(:)

  ! Map size and resolutions
  integer npix,npixc
  real fov,dpix,dpixc,Ompix,fovlon,fovlat
  real theta,phi,theta0,phi0,theta1
  real thetamin,thetamax,phimin,phimax
  real thetap,phip,dtheta,dphi

end module flatskyvars
