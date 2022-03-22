module diffuse_profiles

  use cosmology

  implicit none

  integer :: dp_model, dp_modelp, dp_modeld
  integer :: nprof_pressure 
  integer :: nprof_density  

  real dp_rmax, dp_delta, dp_chi, dp_rvir

  real, dimension(6), parameter :: delta_vir   = (/  200,  500,  200,  200,  200,  200/)

  type :: dp_mzparams
     real :: amplitude, mexponent, zexponent
  end type dp_mzparams

  type :: dp_profparams
     type(dp_mzparams) :: f0, xc, alpha, beta, gamma
  end type dp_profparams

  type :: dp_gnfwparams
     real :: f0, xc, alpha, beta, gamma
  end type dp_gnfwparams

  type(dp_profparams), pointer :: pres_profp(:)
  type(dp_profparams), pointer :: rhog_profp(:), rhod_profp(:)

  type(dp_profparams), pointer :: profp(:) 

  abstract interface
     function gnfw_function(x)
       real :: func_ptr
       real, intent (in) :: x
     end function gnfw_function
  end interface

  procedure(gnfw_function), pointer :: dp_function

  real f0,xc,alpha,beta,gamma
  real amplitude, mexponent, zexponent
    
contains

  subroutine set_dpparams

    nprof_pressure = 3
    nprof_density  = 6

    allocate(pres_profp(nprof_pressure))
    allocate(rhog_profp(nprof_density))
    allocate(rhod_profp(nprof_density))

    ! diffuse pressure profile parameters
    !                              (     BBPS1,   BBPS2,    BBPS3)
    pres_profp%    f0% amplitude = (/ 18.100  ,  7.490 , 20.7   /)
    pres_profp%    f0% mexponent = (/  0.154  ,  0.226 , -0.074 /)
    pres_profp%    f0% zexponent = (/ -0.758  , -0.957 , -0.743 /)
  
    pres_profp%    xc% amplitude = (/  0.497  ,  0.710 ,  0.438 /)
    pres_profp%    xc% mexponent = (/ -0.00865, -0.0833,  0.011 /)
    pres_profp%    xc% zexponent = (/  0.731  ,  0.853 ,  1.01  /)
  
    pres_profp% alpha% amplitude = (/  1.0    ,  1.0   ,  1.0   /)
    pres_profp% alpha% mexponent = (/  0.0    ,  0.0   ,  0.0   /)
    pres_profp% alpha% zexponent = (/  0.0    ,  0.0   ,  0.0   /)

    pres_profp%  beta% amplitude = (/  4.35   ,  4.19  ,  3.82  /)
    pres_profp%  beta% mexponent = (/  0.0393 ,  0.048 ,  0.0375/)
    pres_profp%  beta% zexponent = (/  0.415  ,  0.615 ,  0.535 /)

    pres_profp% gamma% amplitude = (/ -0.3    , -0.3   , -0.3   /)
    pres_profp% gamma% mexponent = (/  0.0    ,  0.0   ,  0.0   /)
    pres_profp% gamma% zexponent = (/  0.0    ,  0.0   ,  0.0   /)

    ! diffuse gas density profile parameters
    !                              (     BBPS1,   BBPS2,   BBPS3,     SIS,   SUS,  NFW)
    rhog_profp%    f0% amplitude = (/  4.0e3  ,  1.9e4 ,  1.5e-4,       0.,    0.,    0./)
    rhog_profp%    f0% mexponent = (/  0.29   ,  0.09  ,  0.14  ,       0.,    0.,    0./)
    rhog_profp%    f0% zexponent = (/ -0.66   , -0.95  , -1.32  ,       0.,    0.,    0./)
  
    rhog_profp%    xc% amplitude = (/  0.5    ,  0.5   ,  0.5   ,       1.,    1.,    1./)
    rhog_profp%    xc% mexponent = (/  0.0    ,  0.0   ,  0.0   ,       0.,    0.,    0./)
    rhog_profp%    xc% zexponent = (/  0.0    ,  0.0   ,  0.0   ,       0.,    0.,    0./)
  
    rhog_profp% alpha% amplitude = (/  0.88   ,  0.70  ,  0.68  ,       1.,    1.,    1./)
    rhog_profp% alpha% mexponent = (/ -0.03   , -0.017 , -0.02  ,       0.,    0.,    0./)
    rhog_profp% alpha% zexponent = (/  0.19   ,  0.27  ,  0.29  ,       0.,    0.,    0./)

    rhog_profp%  beta% amplitude = (/  3.83   ,  4.43  ,  6.4   ,       2.,    0.,    3./)
    rhog_profp%  beta% mexponent = (/  0.04   ,  0.005 ,  0.028 ,       0.,    0.,    0./)
    rhog_profp%  beta% zexponent = (/ -0.025  ,  0.037 , -0.25  ,       0.,    0.,    0./)

    rhog_profp% gamma% amplitude = (/  0.2    ,  0.2   ,  0.2   ,       2.,    0.,    1./)
    rhog_profp% gamma% mexponent = (/  0.0    ,  0.0   ,  0.0   ,       0.,    0.,    0./)
    rhog_profp% gamma% zexponent = (/  0.0    ,  0.0   ,  0.0   ,       0.,    0.,    0./)

    ! diffuse dark matter density profile parameters
    !                              (     BBPS1,   BBPS2,   BBPS3,     SIS,   SUS   NFW-Duffy)
    rhod_profp%    f0% amplitude = (/ 4.1e3   ,  2.8e3 ,  1.5e-4, 200./3., 200.,    200./3 /)
    rhod_profp%    f0% mexponent = (/  -0.23  ,  0.09  ,  0.14  ,       0.,    0.,    0.   /)
    rhod_profp%    f0% zexponent = (/ -1.09   , -0.81  , -1.32  ,       0.,    0.,    0.   /)
  
    rhod_profp%    xc% amplitude = (/  0.5    ,  0.5   ,  0.5   ,       1.,    1.,  10.14  /)
    rhod_profp%    xc% mexponent = (/  0.0    ,  0.0   ,  0.0   ,       0.,    0., -0.081  /)
    rhod_profp%    xc% zexponent = (/  0.0    ,  0.0   ,  0.0   ,       0.,    0., -1.01   /)
  
    rhod_profp% alpha% amplitude = (/  0.77   ,  0.9   ,  0.68  ,       1.,    1.,    1.   /)
    rhod_profp% alpha% mexponent = (/  0.12   ,  0.07  , -0.02  ,       0.,    0.,    0.   /)
    rhod_profp% alpha% zexponent = (/  0.37   ,  0.34  ,  0.29  ,       0.,    0.,    0.   /)

    rhod_profp%  beta% amplitude = (/  3.75   ,  3.76  ,  6.4   ,       2.,    0.,    3.   /)
    rhod_profp%  beta% mexponent = (/  0.007  ,  0.004 ,  0.028 ,       0.,    0.,    0.   /)
    rhod_profp%  beta% zexponent = (/ -0.002  ,  0.024 , -0.25  ,       0.,    0.,    0.   /)

    rhod_profp% gamma% amplitude = (/  1.0    ,   1.0  ,  1.0   ,       2.,    0.,    1.   /)
    rhod_profp% gamma% mexponent = (/  0.0    ,   0.0  ,  0.0   ,       0.,    0.,    0.   /)
    rhod_profp% gamma% zexponent = (/  0.0    ,   0.0  ,  0.0   ,       0.,    0.,    0.   /)

  end subroutine set_dpparams

  !=======================================================================

  real function gnfw_pofmz(m,z,mzparams)

    implicit none

    real m,z

    type(dp_mzparams) :: mzparams

    amplitude = mzparams% amplitude 
    mexponent = mzparams% mexponent
    zexponent = mzparams% zexponent

    gnfw_pofmz = amplitude * (m/1e14)**mexponent * (1+z)**zexponent

    return

  end function gnfw_pofmz

  !=======================================================================

  real function nfw_pofmz(c200,mzparams)
    ! NFW paramters = Duffy et al. arXiV:2008 0804.2486
    implicit none

    real c200
    type(dp_mzparams) :: mzparams

    amplitude = mzparams% amplitude 
    mexponent = mzparams% mexponent
    zexponent = mzparams% zexponent
    
    nfw_pofmz = amplitude * c200**3 * (-c200/(1.+c200) - log(1./c200) + log(1./c200+1.))**(-1.)

    return

  end function nfw_pofmz

  !=======================================================================

  real function gnfw(x)

    implicit none

    real, intent(in) :: x

    gnfw = f0 * (x / xc)**gamma &
         * (1 + (x / xc)**alpha )**(-beta)

    return

  end function gnfw

end module diffuse_profiles
