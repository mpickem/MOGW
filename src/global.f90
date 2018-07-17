module Mglobal
  use, intrinsic :: iso_c_binding, only: sp=>C_FLOAT, dp=>C_DOUBLE

  real(dp)    :: pi
  complex(dp) :: ci
  parameter(ci=(0.d0,1.d0)) ! complex unit i
  parameter(pi=3.1415926535897932384626433832795028841971693993751D0)

! variables (Greensfunctions, interactions, self-energy, ...)
  ! double complex, allocatable :: Giw(:,:,:,:),Gconv(:,:,:)
  ! double complex, allocatable :: P(:,:,:,:)
  ! double complex, allocatable :: V(:,:,:,:),Vend(:,:,:)
  ! double complex, allocatable :: W(:,:,:,:)
  ! double complex, allocatable :: SE(:,:,:,:)

  ! real(dp), allocatable    :: U(:,:,:,:)

end module Mglobal
