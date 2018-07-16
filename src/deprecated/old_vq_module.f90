module vq_module
  use aux
  use hamiltonian_module ! provides nkp,ndim
  use hdf5_module
  implicit none


contains

subroutine read_u(u_tmp, filename_umatrix)
  implicit none
  real(dp)          :: u_tmp(ndim,ndim,ndim,ndim), u_tilde_tmp(ndim,ndim,ndim,ndim)
  real(dp)          :: u_value
  ! complex(dp), intent(out) :: u(ndim**2, ndim**2), u_tilde(ndim**2, ndim**2)
  ! complex(dp), intent(out) :: u_tilde(ndim**2, ndim**2)
  integer           :: n,i,j,k,l,i1,i2
  character(len=80) :: filename_umatrix

  open(21,file=trim(filename_umatrix),status='old')
  read(21,*)
  do n=1,ndim**4
     read(21,*) i, j, k, l, u_value
     u_tmp(i,j,k,l) = u_value
     u_tilde_tmp(i,j,l,k) = u_value ! utilde for ph transverse channel
  enddo
  close(21)

end subroutine read_u
!======================================================================================================

end module vq_module



subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)
  implicit none
  integer,intent(in)  :: Nbands
  integer,intent(in)  :: b1, b2, b3, b4
  integer,intent(out) :: ind

  ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4

end subroutine component2index_band



! converting an index into a band pattern
subroutine index2component_band(Nbands, ind, b1, b2, b3, b4)
  implicit none
  integer,intent(in)  :: Nbands,ind
  integer,intent(out) :: b1, b2, b3, b4
  integer             :: tmp1, tmp2, tmp3, ind_tmp
  integer             :: g1, g2, g3, g4

  ! the proposed back conversion assumes the indices are
  ! given form 0 to max-1
  ind_tmp = ind - 1
  tmp1 = Nbands**3
  tmp2 = Nbands**2
  tmp3 = Nbands

  b1 = ind_tmp/tmp1 + 1
  b2 = (ind_tmp-tmp1*(b1-1))/tmp2 + 1
  b3 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1))/tmp3 + 1
  b4 = (ind_tmp-tmp1*(b1-1)-tmp2*(b2-1)-tmp3*(b3-1)) + 1

end subroutine index2component_band
