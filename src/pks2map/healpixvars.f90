module healpixvars
  
  !-----------------------------------------------------------------------
  ! HEALPIX VARIABLES
  !-----------------------------------------------------------------------
  
  logical hpnest
  integer hpnside,hpnpix,hppix,nlist
  parameter (nest=.true.)
  double precision hpvec1(3),hpvec2(3),hptheta,hpphi,hpangle,hpdomega
  real, allocatable :: hpmap(:,:),hpmapl(:,:)
  integer, allocatable :: hplist(:)
  logical, allocatable :: hpincut(:)
  character(len=80), dimension(1:10) :: hpheader  

  double precision meantaul, meantau

end module healpixvars
