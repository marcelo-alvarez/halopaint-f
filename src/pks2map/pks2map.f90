program pks2map

  ! This code reads in a peak patch catalog and outputs a map.
  ! For Pasted Profiles the contribution from each halo to a given pixel is 
  ! stored in a 3D table as a function of halo mass, redshift, and l.o.s. angle

  use cosmology
  use line_profile
  use pksc
  use haloproject
  use halomodels

  use mpivars
  use fitsvars
  use healpixvars
  use flatskyvars
  use diffuse_profiles

  use fitstools
  use healpix_types
  use head_fits

  use textlib
  
  use random

  implicit none

  ! Filenames
  character *128 filein, fileout, fileout_bn,fileout_fs,fileout_cr,fileout_hp,&
       fileout_dp,tablefile,freqout
  character *4 outcode
  integer ilast, ilastnu

  ! Halo properties
  real xh,yh,zh,mh,rh,redshifth,chih,vh,mfh
  real xmmax,xmmaxl, ymmax,ymmaxl, zmmax,zmmaxl, rmmax

  ! dPdy
  integer, parameter :: ndpdy = 10000
  real,    parameter :: dpdy_max = 1e-4,dpdy_dy=dpdy_max/ndpdy
  real dpdy(ndpdy)

  ! Other variables
  integer i,j,k,jmin,jmax,kmin,kmax,m,nskip,scramble,virflag,testflag

  real mmin,chihview,zview,cut_low,cut_high,zmin,zmax,mmax,m80,m50,m20,m500c_m200m,m500c_m200m_fit,chicmb,omegamz
  real mypi
  real RTHmaxl,RTHmax,rmaxtheta,rmaxphi
  real halo_mean, halo_meanl, halo_meant, halo_meantl, halo_cur, halo_curt

  integer nhalotot, nhalomap 
  integer model, centre, kappa, PSZcut

  ! profile
  character *3 mapcode

  rootid=0
  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myid,ierr)
  call mpi_comm_size(mpi_comm_world,ntasks,ierr)

  ! Usage check
  if(iargc()<3) then
     if(myid==0) write(*,11) 
11   format('Usage: pks2map <filein> <fileout> <tablefile> <profilecode> ',&
          '[<zmin> <zmax> <mmin> <nside> <model> <scramble> <npix> <fov> <PSZcut> <npixcar> <virflag> ',&
          '<testflag>] map choices: tsz, ksz, kap, tau; virflag = 0 --> input mass is m200c, virflag = 1 --> input mass is mvir',&
          ' virflag = 2 --> input mass is m200m')
     call mpi_finalize(ierr)
     stop
  endif

  ! Get commandline
  call      getarg(1,filein)
  call      getarg(2,fileout)
  call      getarg(3,tablefile)
  call      getarg(4,mapcode)
  zmin     = r4arg(5,0.0)
  zmax     = r4arg(6,0.0)
  mmin     = r4arg(7,1.4e13)
  hpnside  = i4arg(8,4096)
  model    = i4arg(9,1)
  scramble = i4arg(10,0)
  npix     = i4arg(11,1024)
  fov      = r4arg(12,10.0)
  PSZcut   = i4arg(13,0)
  npixc    = i4arg(14,1024)
  virflag  = i4arg(15,0)
  testflag = i4arg(16,0)

  if (mapcode == 'kap') then
     profile = 3
     kappa   = 1
  elseif (mapcode == 'tsz') then
     profile = 1
  elseif (mapcode == 'ksz' .or. mapcode == 'tau') then
     profile = 2 
  endif

  num_nuobs = 1

  testcase = .false.
  if(testflag>0) testcase = .true.

  mypi = 4.*atan(1.0)

  if((scramble==1).or.(PSZcut==1)) then
     !Set up random number seeds
     seed = 13579 * (myid+1)
     call ranf(seed,randnum)
  endif
  ! Set field of view and resolution
  fovlon  = 2 * mypi
  fovlat  = mypi 
  fov     = fov / 360 * 2 * mypi ! radians
  dpix    = fov / npix
  dpixc   = mypi / npixc
  Ompix   = (dpix)**2  ! pixel size to convert to brightness temp

  if(myid==0) write(*,12) npix,npix,dpix,fov/2/mypi*360.
12 format(/,/,&
          'Resolution:         ',i0,' x ',i0,' pixels (',1pe8.2,' radians)',/,&
          'Field of view:      ',1pe8.2,' degrees')


  ! Healpix map resolution
  hpnpix   = nside2npix(hpnside)
  hpdomega = 4*mypi / hpnpix

  if(myid==0) write(*,13) hpnside,hpnpix
13 format('Healpix Nside:      ',i0,/,&
          'Number of pixels:   ',i0)

  ! Allocate maps  
  if(myid==0) allocate(fsmap(npix*npix))
  allocate(fsmapl(npix*npix))

  if(myid==0) allocate(crmap(npixc*npixc*2))
  allocate(crmapl(npixc*npixc*2))

  if(myid==0) allocate(hpmap(0:hpnpix-1,1:1))
  allocate(hpmapl(0:hpnpix-1,1:1))
  allocate(hplist(0:hpnpix-1))

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

  ! Create distance-redshift tables
  call set_rofztable
  call set_zofrtable

  chicmb = rofzslow(1100.)
  if(myid==0) write(*,21) chicmb/1e3
21 format('chi_cmb =          ',f5.2,' Gpc')

  cut_low  = rofz(zmin)
  cut_high = rofz(zmax)

  if(myid==0) write(*,14) cut_low/1e3,cut_high/1e3
14 format('cut_low, cut_high:  ',f4.2,1x,f4.2,' Gpc')

  ! Load profile table
  if (mapcode == 'tco')then
    call SFRtabload(tablefile)   ! load SFR Behroozi table
  else
    call loadmaptable(tablefile) ! load table created from make_maptable
    if (mapcode == 'kap') then
        ! add (gas profile) + (DM PROFILE) together in table level to make kappa_tot
        table(:,:,:,3) = fb*table(:,:,:,2) + (1-fb)*table(:,:,:,3)
    endif
  endif

  ! Load halos
  call loadpksc(filein)

  ! Scramble
  if(scramble>0) call scramble_halos

  ! Cut all the halos not within range
  m=0
  do i=1,nhalo 
     xh  = posxyz(1,i)
     yh  = posxyz(2,i) 
     zh  = posxyz(3,i)           

     chih = sqrt(xh**2+yh**2+zh**2)
     if(chih>cut_high.or.chih<cut_low) cycle

     m = m+1
     posxyz(3,m) = posxyz(3,i)
     posxyz(2,m) = posxyz(2,i)
     posxyz(1,m) = posxyz(1,i)
     rth(m)      = rth(i)
  enddo
  nhalo = m

  mass = 4. / 3. * mypi * rho_0 * rth ** 3  

  write(*,*) myid,minval(mass(1:m)),maxval(mass(1:m))

  ! theta0 and phi0 are the center of fov
  ! here we set the center to lie along z-axis
  theta0 = 0
  phi0   = 0
  theta1 = mypi/2 - theta0
  !open(unit=3,file='PSZ2cut_halos_M_z.bin',access='stream')
  ! Loop over halos and project each into map
  do j=1,num_nuobs
  crmapl=0.
  fsmapl=0.
  hpmapl=0.
  meantaul=0.
  maxtheta=0.
  m=0
  nuobsj = nuobs_i - (j-1./2)*dnu

  halo_meanl = 0.0
  halo_meantl = 0.0
  nskip=nhalo/50
  do i=1,nhalo

     xh  = posxyz(1,i)
     yh  = posxyz(2,i) 
     zh  = posxyz(3,i)           

     vh  = vrad(i)

     chih = sqrt(xh**2+yh**2+zh**2)
     redshifth=zofr(chih)
     xh = xh - chihview
     chih = sqrt(xh**2+yh**2+zh**2)

     mh = mass(i)

     !Mass Conversions and Halo Selection
     if(virflag==1) then
        ! convert mass to one with mean density = delta_vir(model) * rhocrit assuming SIS 
        ! deltavirc is the bryan and norman virial density in units of rhocrit     
        ! input catalog is assumed to be masses with rhohalo = deltavirc * rhocrit
        mh = mh * sqrt(deltavirc(redshifth)/delta_vir(model))
     elseif(virflag==2) then
        omegamz = omegam*(1+redshifth)**3/(omegam*(1+redshifth)**3+1-omegam)
        mh = mh * omegamz**0.35 ! roughly scaling for NFW with c ~ 5-10
     endif

     ! mass cut
     if(mass(i)<mmin) cycle
     m = m+1

     !comoving virial radius
     rh  = (3*mh/4/mypi/delta_vir(model)/rhocrit_0/ &
           (omegam*(1+redshifth)**3+1-omegam))**(1./3.)*(1+redshifth)

     ! default multiplication factors 
     mfh = 1.0 
     dp_rmax = 1.0 ! 2.0 original factor

     ! kappa
     if (mapcode == 'kap') then
        mfh = rh*W_Kappa(redshifth,chih,chicmb)

        halo_cur    = mh / chih**2 / 4 / mypi / rho_0
        halo_curt   = rh * interpolate_table_norm(redshifth,mh) / 4 / mypi

        halo_meanl  = halo_meanl  + halo_cur  * W_Kappa(redshifth,chih,chicmb)!14.2e3)
        halo_meantl = halo_meantl + halo_curt * W_Kappa(redshifth,chih,chicmb)!14.2e3)

     endif

     ! ksz and tau
     if (mapcode == 'ksz' .or. mapcode == 'tau') &
          mfh = fb*omegam*h**(4./3) * tau0persigma * (1+redshifth)**2 &
                             *(mh/delta_vir(model))**(1./3)
     if (mapcode == 'ksz') then
        vh = -vh
     else
        vh = 1.0
     endif
     if (mapcode == 'tau') then
        halo_meanl = halo_meanl &
             + mh * fb * sigmatovermump / chih**2 / 4 / mypi &
             * (1+redshifth)**2
        halo_meantl = halo_meantl + mfh * interpolate_table_norm(redshifth,mh) / 4 / mypi
        write(*,*) halo_meanl,halo_meantl
     endif

     if(mapcode=='tsz') dp_rmax = 4.0

     ! project halo profile
     call projecthalo_healpix(xh,yh,zh,rh,chih,mh,redshifth,vh,mfh)
     call projecthalo_flatsky(zh,yh,xh,rh,chih,mh,redshifth,vh,mfh)
     call projecthalo_car(    zh,yh,xh,rh,chih,mh,redshifth,vh,mfh)

  enddo
  !close(3)
  nhalol        = m
  maxthetal     = maxtheta
  maxthetamassl = maxthetamass
  maxthetachihl = maxthetachih
  call mpi_reduce(nhalol,nhalomap,1,mpi_int,mpi_sum,rootid,mpi_comm_world,ierr)
  call mpi_reduce(meantaul,meantau,1,mpi_double_precision,mpi_sum,rootid,&
       mpi_comm_world,ierr)

  call mpi_reduce(halo_meanl, halo_mean, 1,mpi_real,mpi_sum,rootid,&
       mpi_comm_world,ierr)
  call mpi_reduce(halo_meantl,halo_meant,1,mpi_real,mpi_sum,rootid,&
       mpi_comm_world,ierr)

  ! Gather maps from processors
  call mpi_reduce(fsmapl, fsmap,   npix**2,mpi_real,mpi_sum,rootid,mpi_comm_world,ierr)
  call mpi_reduce(crmapl, crmap,2*npixc**2,mpi_real,mpi_sum,rootid,mpi_comm_world,ierr)
  call mpi_reduce(hpmapl, hpmap,    hpnpix,mpi_real,mpi_sum,rootid,mpi_comm_world,ierr)

  ! Output maps from processor 0
  if(myid==0) then

  write(*,15) nhalomap
15 format('Number of halos in map: ',2x,i0)

  fsmom1 = sum(fsmap**1) / npix / npix
  fsmom2 = sum(fsmap**2) / npix / npix
  fsmom3 = sum(fsmap**3) / npix / npix

  crmom1 = sum(crmap**1) / npixc / npixc / 2
  crmom2 = sum(crmap**2) / npixc / npixc / 2
  crmom3 = sum(crmap**3) / npixc / npixc / 2

  hpmom1 = sum(hpmap**1) / hpnpix
  hpmom2 = sum(hpmap**2) / hpnpix
  hpmom3 = sum(hpmap**3) / hpnpix
  
  write(*,17) fsmom1,fsmom2,fsmom3,crmom1,crmom2,crmom3,hpmom1,hpmom2,hpmom3
17 format('Flatsky moments:        ',3(1pe11.3,1x),/,&
          ' CARsky moments:        ',3(1pe11.3,1x),/,&
          'Healpix moments:        ',3(1pe11.3,1x))

  write(*,18) sqrt(sum((fsmap-fsmom1)**2)/npix**2)   ,fsmom1,&
              sqrt(sum((crmap-crmom1)**2)/npixc**2/2),crmom1,&
              sqrt(sum((hpmap-hpmom1)**2)/hpnpix)    ,hpmom1
18 format('Flatsky RMS:            ',1pe11.3,/,&
          'Flatsky Mean:           ',1pe11.3,/,&
          ' CARsky RMS:            ',1pe11.3,/,&
          ' CARsky Mean:           ',1pe11.3,/,&
          'Healpix RMS:            ',1pe11.3,/,&
          'Healpix Mean:           ',1pe11.3,/)

  if(mapcode == 'kap'.or.mapcode == 'tau') write(*,19) halo_mean,halo_meant,&
       hpmom1/halo_mean,hpmom1/halo_meant
19 format('sum over halo ratio     ',1pe11.3,1x,1pe11.3,1x,1pe11.3,1x,1pe11.3/)

  write(*,20) minval(fsmap),maxval(fsmap)
20 format('min and max (healpix):  ',1pe11.3,1x,1pe11.3/)

  write(*,40) minval(crmap),maxval(crmap)
40 format('min and max (CAR):      ',1pe11.3,1x,1pe11.3/)
     
  ! Get dp/dy
  dpdy=0
  do i=0,hpnpix-1
     m = int(hpmap(i,1)/dpdy_dy)+1
     if(m>0.and.m<=ndpdy) dpdy(m)=dpdy(m)+1
  enddo
  dpdy = dpdy / real(hpnpix)

  ! Check mean
  fsmom1=0.0
  do i=1,ndpdy
     fsmom1=fsmom1+(i-0.5)*dpdy_dy*dpdy(i)
  enddo

  do i=2,ndpdy
     dpdy(i)=dpdy(i-1)+dpdy(i)
  enddo  

  ! Create filenames 
  ilast   = indlnb(fileout)
  write (freqout, "(I0.3)") j
  if ((mapcode == 'tco') .or. (mapcode == 'thi')) then
     fileout_fs=fileout(1:ilast)//'_fs.map'//trim(freqout)
     fileout_hp=fileout(1:ilast)//'_hp.fits'//trim(freqout)
     fileout_dp=fileout(1:ilast)//'_dp.bin'//trim(freqout)
  else
     fileout_fs=fileout(1:ilast)//'_fs.map'
     fileout_cr=fileout(1:ilast)//'_car.map'
     fileout_hp=fileout(1:ilast)//'_hp.fits'
     fileout_dp=fileout(1:ilast)//'_dp.bin'
  endif

  ! P(<y) file
  open(unit=1,file=fileout_dp,form='binary')
  write(1) ndpdy,dpdy_max
  write(1) dpdy
  close(1)

  ! Flat sky binary file
  open(unit=1,file=fileout_fs,form='binary')
  write(1) npix,npix,fov,fov
  write(1) fsmap
  close(1)
     
  ! CAR binary file
  open(unit=1,file=fileout_cr,form='binary')
  write(1) 2*npixc,npixc,fovlon,fovlat
  write(1) crmap
  close(1)
     
  ! Healpix FITS file
  hpheader(:)=''
  call add_card(hpheader,'NSIDE',hpnside,'the nside of the map')
  call add_card(hpheader,'ORDERING','RING','the nside of the map')
  call output_map(hpmap,hpheader,fileout_hp)

  endif
  enddo !end loop over freq for mapcode == itco and mapcode == ithi
  
  if(myid==0) write(*,*)

  call mpi_finalize(ierr)

  stop

81 format('000',i1)
82 format( '00',i2)
83 format(  '0',i3)
84 format(      i4)

end program pks2map
