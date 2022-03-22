module pksc  

  use cosmology
  use mpivars 
  use random 
  implicit none

  real delta_pksc
  parameter(delta_pksc=200)

  ! Halo arrays
  real, allocatable :: posxyz(:,:),velxyz(:,:),rth(:),mass(:),vrad(:)
  integer nhalo, nhalol 
  real RTHmaxtot
  integer(kind=8) offset_num_bytes

  ! For random number generation                                                                               
  integer,dimension( IRandNumSize ) :: seed
  real(kind=8) :: randnum

contains

  subroutine loadpksc(filename)

  !=======================================================================
  !LOAD IN HALOS
  !=======================================================================

    implicit none
    
    integer i,j,idum,m

    character *128 filename

    real rho, mhtest, zhtest, rthtest

    if(testcase) then

       nhalo = 1

       allocate(posxyz(3,nhalo))
       allocate(velxyz(3,nhalo))
       allocate(rth(nhalo))
       allocate(mass(nhalo))
       allocate(vrad(nhalo))

       mhtest  = 1e14
       zhtest  = 3.2
       rho     = 2.775e11 * 0.31 * 0.68**2
       rthtest = (3*mhtest/4./3.14159/rho)**(1./3.)

       posxyz(3,1) = 0
       posxyz(2,1) = 0
       posxyz(1,1) = rofz(zhtest)
       velxyz(1,1) = 0.
       velxyz(2,1) = 0.
       velxyz(3,1) = 0.
       rth(1)      = rthtest
       vrad(1)     = 100. / 299792.458

       if(myid==0) write(*,*) 'TEST: mhtest, zhtest, rthtest, chitest = ',&
            mhtest,zhtest,rthtest,rofz(zhtest)

       return
       
    endif

    open(unit=1,file=filename,form='binary')
    read(1) nhalo, RTHmaxtot
    close(1)

    if(myid==0) write(*,*) 'nhalotot = ',nhalo    

    !Read in only the appropriate nhalo/ntasks, not whole catalogue
    nhalol = int((nhalo-1)/ntasks)+1          ! even number per processor
    nhalo  = min(nhalol, nhalo - nhalol*myid) ! last processor may not have 
                                              ! same number of halos

    if(myid==0) write(*,*) 'nhalo local = ',nhalo    

    allocate(posxyz(3,nhalo))
    allocate(velxyz(3,nhalo))
    allocate(rth(nhalo))
    allocate(mass(nhalo))
    allocate(vrad(nhalo))

    offset_num_bytes = int(10,8)*int(nhalol,8)*int(myid,8)*int(4,8)+int(12+1,8)
    open(unit=1,file=filename,access='stream',form='binary')
    read(1,pos=offset_num_bytes) ((posxyz(j,i),j=1,3),&
             (velxyz(j,i),j=1,3),&
             rth(i),idum,idum,idum,i=1,nhalo)
    close(1)

    vrad = (posxyz(1,:)*velxyz(1,:) + posxyz(2,:)*velxyz(2,:)+ posxyz(3,:)*velxyz(3,:))
    vrad = vrad/(posxyz(1,:)**2 + posxyz(2,:)**2 + posxyz(3,:)**2)**(1./2)
    vrad = vrad / 299792.458 ! in units of c 
    
    return

  end subroutine loadpksc


  subroutine scramble_halos
  !=======================================================================
  !Scramble halos in theta and phi
  !=======================================================================
    implicit none

    ! For random number generation   
    !integer,dimension( IRandNumSize ) :: seed
    integer i

    !real(kind=8) :: randnum
    real xh,yh,zh,rh,chih
    double precision muran, phiran

!    seed = 13579 * (myid+1)
!    call ranf(seed,randnum)
    do i=1,nhalo
        xh   = posxyz(1,i)
        yh   = posxyz(2,i)
        zh   = posxyz(3,i)

        chih = sqrt(xh**2+yh**2+zh**2)
!        write(*,*) '1',chih

        call ranf(seed,randnum)
        muran  = 2 * randnum - 1
        call ranf(seed,randnum)
        phiran = 2 * 3.14159 * randnum

        rh = sqrt(1-muran**2) * chih

!        xh = rh    * cos(phiran)
!        yh = rh    * sin(phiran)
        zh = muran * chih

        posxyz(1,i) = xh
        posxyz(2,i) = yh
        posxyz(3,i) = zh

        chih = sqrt(xh**2+yh**2+zh**2)
!        write(*,*) '2',chih

     enddo
     write(*,*) 'done scrambling'

  end subroutine scramble_halos

  subroutine center_halos
  !=======================================================================
  !Center halos on the most massive object in the map (x,0,0)
  !=======================================================================

    implicit none

    ! For random number generation   
    integer i, imaxval(1)

    real xh,yh,zh,rh,chih
    real xmmax, xmmaxl, ymmax, ymmaxl, zmmax, zmmaxl, rmmax
    real RTHmax, RTHmaxl, rmaxtheta, rmaxphi

    !find largest remaining halo across all processors
    RTHmaxl = maxval(rth(1:nhalo))
    imaxval = maxloc(rth(1:nhalo))
    !write(*,*) "DEBUG maxloc, maxval RTH", imaxval, rth(imaxval(1))
    call mpi_allreduce(RTHmaxl,RTHmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    !Distribute position of largest halo to all processors
    xmmaxl = -1e5
    ymmaxl = -1e5
    zmmaxl = -1e5
    if(RTHmaxl == RTHmax) then
       xmmaxl = posxyz(1,imaxval(1))
       ymmaxl = posxyz(2,imaxval(1))
       zmmaxl = posxyz(3,imaxval(1))

       rmmax = sqrt(xmmaxl**2+ymmaxl**2+zmmaxl**2)
       write(*,19) RTHmax,rmmax, xmmaxl, ymmaxl, zmmaxl
19     format(/,'Before rotation largest halo is at',/,&
            'RTH, distance:    ',2(1pe9.3,1x),/,&
            'x, y, z:          ',3(1e10.3,1x),/)
       
    endif
    call mpi_allreduce(xmmaxl,xmmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    call mpi_allreduce(ymmaxl,ymmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)
    call mpi_allreduce(zmmaxl,zmmax,1,mpi_real,mpi_max,mpi_comm_world,ierr)


    !rotate so largest halo is in first quadrant to make angles easier   
    if(xmmax <0.0)   posxyz(1,:) = -posxyz(1,:)
    if(ymmax <0.0)   posxyz(2,:) = -posxyz(2,:)
    if(zmmax <0.0)   posxyz(3,:) = -posxyz(3,:)
    xmmax = abs(xmmax)
    ymmax = abs(ymmax)
    zmmax = abs(zmmax)

    ! angle to rotate by in z-x plane
    rmaxtheta = asin(zmmax/sqrt(xmmax**2+zmmax**2))  
    xmmaxl    = xmmax/abs(xmmax)*sqrt(xmmax**2+zmmax**2)
    ! angle to rotate by in z-y plane
    rmaxphi   = asin(ymmax/sqrt(ymmax**2+xmmaxl**2))  

    ! Rotate by angles defined by largest halo                          
    do i=1,nhalo
       xh =   posxyz(1,i)*cos(rmaxtheta) + posxyz(3,i)*sin(rmaxtheta)
       zh =  -posxyz(1,i)*sin(rmaxtheta) + posxyz(3,i)*cos(rmaxtheta)

       yh =  -xh*sin(rmaxphi) + posxyz(2,i)*cos(rmaxphi)
       xh =   xh*cos(rmaxphi) + posxyz(2,i)*sin(rmaxphi)

       posxyz(1,i) = zh 
       posxyz(2,i) = yh
       posxyz(3,i) = xh

    enddo

    ! Double check rotation worked and find new distance to largest halo
    if(RTHmaxl == RTHmax) then
       xmmax = posxyz(3,imaxval(1))
       ymmax = posxyz(2,imaxval(1))
       zmmax = posxyz(1,imaxval(1))
       rmmax = sqrt(xmmax**2+ymmax**2+zmmax**2)
       write(*,20) RTHmax,rmmax, xmmax, ymmax, zmmax
20     format(/,'After rotation largest halo is at',/,&
            'RTH, distance:    ',2(1pe9.3,1x),/,&
            'x, y, z:          ',3(1e10.3,1x),/)
    endif

  end subroutine center_halos

  logical function select_halo(m50,m80,m,i)
    
    implicit none
    
    ! Masses
    real  m50,  m80,  m, lm, dfdlnm
    real lm50, lm80, lm0, lm100, f
    integer i

    lm   = log(m)    
    lm50 = log(m50)
    lm80 = log(m80)
    dfdlnm = 0.3 / (lm80 - lm50)
    
    lm100 = lm80 + 0.2 / dfdlnm
    lm0   = lm50 - 0.5 / dfdlnm
    
    select_halo = .false.
    if(lm < lm0) then
       return
    elseif(lm > lm100) then
       select_halo = .true.
       return
    else
        call ranf(seed,randnum)
       f = (lm - lm0) * dfdlnm
       if(f>randnum) select_halo = .true.
!       write(*,*) "seed, f, randnum = ",seed, f, randnum
    endif
    !write(*,*) "lm50,lm80,lm100,lm0,dfdlnm = ",lm50,lm80,lm100,lm0,dfdlnm

    return
    
  end function select_halo

end module pksc

