module aux
  integer, parameter :: dp=kind(0.d0) ! double precision - use real(dp) :: instead of double precision
  integer, parameter :: sp=kind(0.0)  ! single precision
  real(dp)           :: pi
  complex(dp)        :: ci
  parameter(ci=(0.d0,1.d0)) ! complex unit i
  parameter(pi=3.1415926535897932384626433832795028841971693993751D0)
  character(len=80)  :: outfolder, hamfolder
  character(len=80)  :: filename_umatrix, filename_vq, filename_dmft
  real(dp)           :: beta
  integer            :: nw
end module aux
