module Mglobal
  public
  integer, parameter :: dp=kind(0.d0) ! double precision - use real(dp) :: instead of double precision
  integer, parameter :: sp=kind(0.0)  ! single precision

  real(dp)           :: pi
  complex(dp)        :: ci
  parameter(ci=(0.d0,1.d0)) ! complex unit i
  parameter(pi=3.1415926535897932384626433832795028841971693993751D0)

  integer,save            :: ikstart, ikend

  real(dp),save           :: beta
  integer,save            :: nw
  integer,save            :: mu
  integer,save            :: occ
  integer,save            :: ndim
  integer,save            :: nkp

  character(len=100),save   :: file_hmlt
  character(len=100),save   :: file_hmlt_kpq
  character(len=100),save   :: file_hmlt_mq

! variables (Greensfunctions, interactions, self-energy, ...)
  complex(dp),save, allocatable :: Giw(:,:,:,:),Gconv(:,:,:)
  complex(dp),save, allocatable :: P(:,:,:,:)
  complex(dp),save, allocatable :: V(:,:,:,:),Vend(:,:,:)
  complex(dp),save, allocatable :: W(:,:,:,:)
  complex(dp),save, allocatable :: SE(:,:,:,:)

  real(dp),save, allocatable    :: U(:,:,:,:)

  contains

  subroutine clear_file_arrays()
    implicit none
    file_hmlt = ''
    file_hmlt_kpq = ''
    file_hmlt_mq = ''
  end subroutine clear_file_arrays

end module Mglobal
