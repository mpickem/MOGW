module lapack_module
  implicit none
  private
  public inverse_matrix
  integer, public :: lierr
  
!
!   Calculation of the inverse matrix
!
  interface inverse_matrix
    module procedure inverse_matrix_d,inverse_matrix_z
  end interface inverse_matrix
  
  contains
  
  subroutine inverse_matrix_z( a )
    implicit none
     complex(kind=8), intent (inout) :: a(:,:)
    integer                     :: dim,info
    integer,        allocatable :: ipiv(:)
    complex(kind=8), allocatable :: lwork(:)
  
      lierr=0
    dim = size( a,1 )
!      write(*,*)'dim',dim
    allocate( ipiv(dim),lwork(dim) )

    call zgetrf( dim,dim,a,dim,ipiv,info )
    if( info /= 0 )then
      write(66,"(//,'INVERSE_MATRIX: error in ZGETRF, info =',1x,i3)") info
      lierr=1
!      stop
    end if
    call zgetri( dim,a,dim,ipiv,lwork,dim,info )
    if( info /= 0 )then
      write(66,"(//,'INVERSE_MATRIX: error in ZGETRI, info =',1x,i3)") info
      lierr=1
!      stop
    end if
  
    deallocate( ipiv,lwork )
  
  end subroutine inverse_matrix_z
  
  subroutine inverse_matrix_d( a )
    implicit none
    double precision, intent (inout) :: a(:,:)
    integer                          :: dim,info
    integer,          allocatable    :: ipiv(:)
    double precision, allocatable    :: lwork(:)
  
    dim = size( a,1 )
    allocate( ipiv(dim),lwork(dim) )

    call dgetrf( dim,dim,a,dim,ipiv,info )
    if( info /= 0 )then
      write(66,"(//,'INVERSE_MATRIX: error in DGETRF, info =',1x,i3)") info
      stop
    end if
    call dgetri( dim,a,dim,ipiv,lwork,dim,info )
    if( info /= 0 )then
      write(66,"(//,'INVERSE_MATRIX: error in DGETRI, info =',1x,i3)") info
      stop
    end if
  
    deallocate( ipiv,lwork )
  
  end subroutine inverse_matrix_d
  
end module lapack_module
