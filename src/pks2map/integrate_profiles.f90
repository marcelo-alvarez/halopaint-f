module integrate_profiles

  use cosmology
  use diffuse_profiles
  use halomodels

  implicit none

  real, parameter :: dsR=5e-2

contains

  real function integrate_profile(theta)

    ! integrate profile along line of sight theta away from halo center
    ! s = chi / rvir is the dimensionless integration variable

    implicit none

    integer i,nint
    real theta
    real b,s,x,s0
    real y
    double precision ds0,y0
    
    integrate_profile = 0

    b = dp_chi * sin(theta) / dp_rvir 

    if(b > dp_rmax) return

    ! this is for a spherical cut at rmax
    s0 = sqrt(dp_rmax**2-b**2) 
    nint = int(s0 / (dSR*dp_rvir))

    if (nint==0) return

    ds0 = s0/nint

    y = 0
    do i=1,nint
       s = (i-0.5) * ds0
       x = sqrt(s**2+b**2)
       y = y + dp_function(x)
    enddo

    y = y * ds0

    ! factor of two is because we only integrate half of it (int_0^rmax)
    integrate_profile = 2*y 

    return 

  end function integrate_profile

end module integrate_profiles
