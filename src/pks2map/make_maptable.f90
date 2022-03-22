program make_maptable

  use maptable
  use cosmology
  use textlib

  implicit none

  ! Output filename
  character *128 fileout

  ! Usage check
  if(iargc()/=3) then
     write(*,11) 
11   format('Usage: make_maptable <fileout> <pressure model> <density model>')
     stop
  endif

  ! Get commandline
  call getarg(1,fileout)
  dp_modelp = i4arg(2,1)
  dp_modeld = i4arg(3,1)

  ! Set background cosmology
  omegam = 0.31
  omegab = 0.049
  omegal = 1 - omegam
  h      = 0.68
  sigma8 = 0.82
  ns     = 0.965
  w      = -1
  fb     = omegab/omegam

  rho_0     = 2.775e11*omegam*h**2
  rhocrit_0 = 2.775e11*h**2  

  ! Set diffuse profile parameters
  call set_dpparams

  ! Make table 
  call makemaptable(fileout)

  stop
  
end program make_maptable
