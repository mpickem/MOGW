module Mglobal
  integer, parameter :: dp=kind(0.d0) ! double precision - use real(dp) :: instead of double precision
  integer, parameter :: sp=kind(0.0)  ! single precision
  real(dp)           :: pi
  complex(dp)        :: ci
  parameter(ci=(0.d0,1.d0)) ! complex unit i
  parameter(pi=3.1415926535897932384626433832795028841971693993751D0)

  integer            :: ikstart, ikend

  real(dp)           :: beta
  integer            :: nw
  integer            :: mu
  integer            :: occ
  integer            :: ndim
  integer            :: nkp

  character(len=200)   :: file_hmlt
  character(len=200)   :: file_hmlt_kpq
  character(len=200)   :: file_hmlt_mq

! variables (Greensfunctions, interactions, self-energy, ...)
  complex(dp), allocatable     :: Giw(:,:,:,:),Gconv(:,:,:)
  complex(dp), allocatable     :: P(:,:,:,:)
  complex(dp), allocatable     :: V(:,:,:,:),Vend(:,:,:)
  complex(dp), allocatable     :: W(:,:,:,:)
  complex(dp), allocatable     :: SE(:,:,:,:)

  real(dp), allocatable        :: U(:,:,:,:)

end module Mglobal
