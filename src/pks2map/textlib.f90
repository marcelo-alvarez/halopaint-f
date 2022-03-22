MODULE textlib

  implicit none

  real*8,  parameter :: r8_never = -1.0d-60
  integer, parameter :: i4_never = -2000000000
  
  !-----------------------------------------------------------
  !     Copyright Matthew W. Choptuik, 1988-1997
  !
  !     User routines:
  !
  !        i4arg
  !        r4arg
  !        getu
  !        indlnb
  !
  !     Support routines
  !
  !        onedot
  !        s2i4
  !        s2r8
  !-----------------------------------------------------------
  
  !-----------------------------------------------------------
  !     Converts argument 'argno' to integer*4 and returns 
  !     value or 'defval' if parsing error or if argument is 
  !     single dot ('.')
  !-----------------------------------------------------------

CONTAINS 
      
  integer function i4arg(argno,defval)
    
    implicit       none
    
    integer        argno
    integer        defval
    
    character*32   argi
    
    if( argno .le. command_argument_count() ) then
       call getarg(argno,argi)
       if( onedot(argi) ) then
          i4arg = defval
       else
          i4arg = s2i4(argi)
          if( i4arg .eq. i4_never ) then
             i4arg = defval
          end if
       end if
    else
       i4arg = defval
    end if
    
    return
    
  end function i4arg
  
  !-----------------------------------------------------------
  !     Converts argument 'argno' to real*8 and returns 
  !     value or 'defval' if parsing error or if argument is 
  !     single dot ('.')
  !-----------------------------------------------------------

  real function r4arg(argno,defval)

    implicit       none
    
    integer        argno
    real         defval
    
    character*32   argi
    
    if( argno .le. command_argument_count() ) then
       call getarg(argno,argi)
       if( onedot(argi) ) then
          r4arg = defval
       else 
          r4arg = s2r8(argi)
          if( r4arg .eq. r8_never ) then
             r4arg = defval
          end if
       end if
    else
       r4arg = defval
    end if
    
    return
    
  end function r4arg
  
  !-----------------------------------------------------------
  
  double precision function r8arg(argno,defval)
    
    implicit       none
    
    integer        argno
    real*8         defval
    
    character*32   argi   
    
    if( argno .le. command_argument_count() ) then
       call getarg(argno,argi)
       if( onedot(argi) ) then
          r8arg = defval
       else
          r8arg = s2r8(argi)
          if( r8arg .eq. r8_never ) then
             r8arg = defval
          end if
       end if
    else
       r8arg = defval
    end if
    
    return
    
  end function r8arg
  
  !-----------------------------------------------------------
  !     Returns index of last non-blank character in S.
  !-----------------------------------------------------------

  integer function indlnb(s)
 
    character*(*)    s
    
    do indlnb = len(s) , 1 , -1
       if( s(indlnb:indlnb) .ne. ' ' ) RETURN
    end do
    indlnb = 0
    
    return
    
  end function indlnb
  
  !-----------------------------------------------------------
  !     Returns first unit number .ge. umin not attached to 
  !     a file.
  !-----------------------------------------------------------

  integer function getu()
    
    implicit      none
    
    integer, parameter :: umin = 10
    integer, parameter :: umax = 99
    
    integer       u
    logical       opened
    
    getu = -1
    do u = umin , umax
       inquire(unit=u,opened=opened)
       if( .not. opened ) then
          getu = u
          return
       end if
    end do
    write(0,*) 'getu: Panic--no available unit number'
    stop
    
  end function getu
  
  !-----------------------------------------------------------
  !     Utility routine for parameter handling.
  !-----------------------------------------------------------

  logical function onedot(s)
    
    character*(*)     s
    
    onedot = s(1:1) .eq. '.'  .and.  indlnb(s) .eq. 1
    
    return
    
  end function onedot
  
  !-----------------------------------------------------------
  !     Converts string to integer.
  !-----------------------------------------------------------

  integer function s2i4(s)
    
    implicit      none
    
    character*(*) s     
    
    character*32  buffer
    integer       lens
    
    integer, parameter :: default = 0
    
    lens = indlnb(s)
    if( lens .gt. 99 ) then
       write(*,*) '>>> s2i4:: String too long for conversion.'
       s2i4 = default
    else 
       if( lens .le. 9 ) then
          write(buffer,100) lens
       else 
          write(buffer,101) lens
       end if
100    format('(I',i1,')')
101    format('(I',i2,')')
       read(s,fmt=buffer,end=900,err=900) s2i4
    end if
    
    return
    
900 s2i4 = i4_never
    return
    
  end function s2i4
  
  !-----------------------------------------------------------
  !     Converts string to real*8 value.
  !-----------------------------------------------------------

  double precision function s2r8(s)

    implicit      none
    
    character*(*) s     
    
    character*32  buffer
    integer       lens
    
    double precision, parameter :: default = 0.0d0
    
    lens = indlnb(s)
    if( lens .gt. 99 ) then
       write(*,*) '>>> s2r8:: String too long for conversion.'
       s2r8 = default
    else 
       if( lens .le. 9 ) then
          write(buffer,100) lens
       else 
          write(buffer,101) lens
       end if
100    format('(G',i1,'.0)')
101    format('(G',i2,'.0)')
       read(s,fmt=buffer,end=900,err=900) s2r8
    end if
    
    return
    
900 s2r8 = r8_never 
    return
    
  end function s2r8
  
END MODULE textlib
