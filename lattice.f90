
program lat
  real(kind=dp), parameter    :: eps5  = 1.0e-5_dp
  real(kind=dp) :: real_lattice(3,3)
  real(kind=dp) :: recip_lattice(3,3)
  real(kind=dp) :: cell_volume
  real(kind=dp) :: real_metric(3,3)
  real(kind=dp) :: recip_metric(3,3)
  integer       :: stdout 
  character(len=50) :: seedname
  integer :: chk_unit
  character(len=9) :: pos,stat
  
  sendname='lat'
  chk_unit=301
  stdout=302
  pos='append'

  inquire(file=trim(seedname)//'.wout',exist=wout_found)
  if (wout_found) then
     stat='old'
  else
     stat='replace'
  endif
  open(unit=stdout,file=trim(seedname)//'.wout',status=trim(stat),position=trim(pos))

  open(unit=chk_unit,file=trim(seedname)//'.in',form='unformatted')
  read(chk_unit) ((real_lattice(i,j),i=1,3),j=1,3)  ! Real lattice
  write(stdout,'(a)') "Real lattice: read."
    
  write(chk_unit, '(9G25.17)') ((real_lattice(i,j),i=1,3),j=1,3)        ! Real lattice

  cell_volume = real_lattice(1,1)*(real_lattice(2,2)*real_lattice(3,3)-real_lattice(3,2)*real_lattice(2,3)) +&
                real_lattice(1,2)*(real_lattice(2,3)*real_lattice(3,1)-real_lattice(3,3)*real_lattice(2,1)) +& 
                real_lattice(1,3)*(real_lattice(2,1)*real_lattice(3,2)-real_lattice(3,1)*real_lattice(2,2))
        
  call utility_recip_lattice(real_lattice,recip_lattice,cell_volume)
  call utility_metric(real_lattice,recip_lattice,real_metric,recip_metric)
  close(stdout)
  close(chk_unit)

  contains

  !===================================================================
  subroutine utility_recip_lattice (real_lat,recip_lat,volume)  !
    !==================================================================!
    !                                                                  !
    !  Calculates the reciprical lattice vectors and the cell volume   !
    !                                                                  !
    !===================================================================

    use w90_constants,  only : dp,twopi,eps5
    use w90_io,         only : io_error

    implicit none
    real(kind=dp), intent(in)  :: real_lat (3, 3)
    real(kind=dp), intent(out) :: recip_lat (3, 3)  
    real(kind=dp), intent(out) :: volume

    recip_lat(1,1)=real_lat(2,2)*real_lat(3,3)-real_lat(3,2)*real_lat(2,3)
    recip_lat(1,2)=real_lat(2,3)*real_lat(3,1)-real_lat(3,3)*real_lat(2,1)
    recip_lat(1,3)=real_lat(2,1)*real_lat(3,2)-real_lat(3,1)*real_lat(2,2)
    recip_lat(2,1)=real_lat(3,2)*real_lat(1,3)-real_lat(1,2)*real_lat(3,3)
    recip_lat(2,2)=real_lat(3,3)*real_lat(1,1)-real_lat(1,3)*real_lat(3,1)
    recip_lat(2,3)=real_lat(3,1)*real_lat(1,2)-real_lat(1,1)*real_lat(3,2)
    recip_lat(3,1)=real_lat(1,2)*real_lat(2,3)-real_lat(2,2)*real_lat(1,3)
    recip_lat(3,2)=real_lat(1,3)*real_lat(2,1)-real_lat(2,3)*real_lat(1,1)
    recip_lat(3,3)=real_lat(1,1)*real_lat(2,2)-real_lat(2,1)*real_lat(1,2)

    volume=real_lat(1,1)*recip_lat(1,1) + &
         real_lat(1,2)*recip_lat(1,2) + &
         real_lat(1,3)*recip_lat(1,3)  


    if( abs(volume) < eps5 ) then
       call io_error(' Found almost zero Volume in utility_recip_lattice')
    end if

    recip_lat=twopi*recip_lat/volume
    volume=abs(volume)

    return

  end subroutine utility_recip_lattice
  
  !===================================================================
  subroutine utility_metric(real_lat,recip_lat, &
       real_metric,recip_metric)
    !==================================================================!
    !                                                                  !
    !  Calculate the real and reciprical space metrics                 !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out) :: real_metric(3,3)
    real(kind=dp), intent(out) :: recip_metric(3,3)

    integer :: i,j,l

    real_metric=0.0_dp ; recip_metric=0.0_dp

    do j=1,3
       do i=1,j
          do l=1,3
             real_metric(i,j)=real_metric(i,j)+real_lat(i,l)*real_lat(j,l)
             recip_metric(i,j)=recip_metric(i,j)+recip_lat(i,l)*recip_lat(j,l)
          enddo
          if(i.lt.j) then
             real_metric(j,i)=real_metric(i,j)
             recip_metric(j,i)=recip_metric(i,j)
          endif
       enddo
    enddo

  end subroutine utility_metric
  
  
  
    !===================================================================
  subroutine utility_frac_to_cart(frac,cart,real_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from fractional to Cartesian coordinates                !
    !                                                                  !
    !===================================================================  
    implicit none

    real(kind=dp), intent(in)  :: real_lat(3,3)
    real(kind=dp), intent(in)  :: frac(3)
    real(kind=dp), intent(out) :: cart(3)

    integer :: i

    do i=1,3
       cart(i)=real_lat(1,i)*frac(1) + real_lat(2,i)*frac(2) + real_lat(3,i)*frac(3) 
    end do

    return

  end subroutine utility_frac_to_cart


  !===================================================================
  subroutine utility_cart_to_frac(cart,frac,recip_lat)
    !==================================================================!
    !                                                                  !
    !  Convert from Cartesian to fractional coordinates                !
    !                                                                  !
    !===================================================================  
    use w90_constants, only : twopi
    implicit none

    real(kind=dp), intent(in)  :: recip_lat(3,3)
    real(kind=dp), intent(out)  :: frac(3)
    real(kind=dp), intent(in)  :: cart(3)

    integer :: i

    do i=1,3
       frac(i)=recip_lat(i,1)*cart(1) + recip_lat(i,2)*cart(2) + recip_lat(i,3)*cart(3) 
    end do

    frac=frac/twopi


    return

  end subroutine utility_cart_to_frac
  
  
  
    subroutine utility_translate_home(vec,real_lat,recip_lat)
    !========================================================!
    !                                                        !
    !        Translate a vector to the home unit cell        !
    !                                                        !
    !========================================================!

    implicit none

    real(kind=dp), intent(inout) :: vec(3)
    real(kind=dp), intent(in)    :: real_lat(3,3)
    real(kind=dp), intent(in)    :: recip_lat(3,3)

    ! <<<local variables>>>
    integer       :: ind
    real(kind=dp) :: r_home(3),r_frac(3)
    real(kind=dp) :: shift

    r_home=0.0_dp;r_frac=0.0_dp

    ! Cartesian --> fractional
    call utility_cart_to_frac(vec,r_frac,recip_lat)
    ! Rationalise to interval [0,1]
    do ind=1,3
       if (r_frac(ind).lt.0.0_dp) then
          shift=real(ceiling(abs(r_frac(ind))),kind=dp)
          r_frac(ind)=r_frac(ind)+shift
       endif
       if (r_frac(ind).gt.1.0_dp) then
          shift=-real(int(r_frac(ind)),kind=dp)
          r_frac(ind)=r_frac(ind)+shift
       endif
    enddo
    ! Fractional --> Cartesian
    call utility_frac_to_cart(r_frac,r_home,real_lat)

    vec = r_home

    return
  end subroutine utility_translate_home
  


  end program lat
  
