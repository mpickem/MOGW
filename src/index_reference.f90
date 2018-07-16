module Mindex
  use Mglobal
  implicit none
  private
  public :: index_, inverse_l, inverse_r

  contains

  subroutine index_(x,y)  ! x - left compound index --- y - right compound index
  ! this can also be done via a function written like below - however it is faster to do one loop since every index is used anyways
    integer, intent(out) :: x(ndim,ndim),y(ndim,ndim)  ! x = 1 ... ndim**2 = y ---- Polarization arrays have to be changed
    integer              :: cnt,i,j

    cnt=1

    do i=1,ndim
    do j=1,ndim
      x(i,j) = cnt ! 11, 12, 13, ..., 21, 22, 23, ...
      y(j,i) = cnt ! 11, 21, 31, ..., 12, 22, 32, ...
      cnt = cnt+1
    enddo
    enddo

  end subroutine index_

  integer function inverse_l(i) ! i -> L (from LL_) --- j -> L_ (from L_L)
    integer,intent(in) :: i
    inverse_l=(i-1)/ndim+1 ! integer division
    return
  end function inverse_l

  integer function inverse_r(i) ! i -> L_ (from LL_) --- j -> L (from L_L)
    integer,intent(in) :: i
    inverse_r=mod(i-1,ndim)+1
    return
  end function inverse_r

end module Mindex
