module vq_module
  use aux
  use hamiltonian_module ! provides nkp,ndim
  use hdf5_module
  implicit none


contains

!===============================================================================================================================
subroutine read_vq(iq, vq, filename_vq)
  implicit none
  integer, intent(in) :: iq
  ! complex(kind=8), intent(out) :: v(ndim**2,ndim**2)
  complex(kind=8) :: vq(ndim,ndim,ndim,ndim)
  integer(hid_t) :: vq_file_id, grp_id, iq_id, iq_space_id
  integer :: err, ind, i, j, k, l, i1, i2
  integer(hid_t) :: nmembers, imembers, itype
  character(len=20) :: name_buffer
  character(len=80) :: filename_vq
  integer(hsize_t), dimension(2) :: iq_dims, iq_maxdims
  integer(hsize_t), dimension(1) :: vq_dims
  double precision :: vq_tmp_r(nkp), vq_tmp_i(nkp)

  call h5fopen_f(trim(filename_vq), h5f_acc_rdonly_f, vq_file_id, err)

  call h5dopen_f(vq_file_id, ".axes/Q-points", iq_id, err)
  call h5dget_space_f(iq_id, iq_space_id, err)
  call h5sget_simple_extent_dims_f(iq_space_id, iq_dims, iq_maxdims, err)
  if(iq_dims(2) .ne. nkp) then
     write(*,*) 'Inconsistent number of q-points in V^q!', iq_dims(2),'/',nkp
     stop
  endif
 
  vq = 0.d0

  call h5gn_members_f(vq_file_id, "/", nmembers, err)
  do imembers = 1,nmembers - 1
     call h5gget_obj_info_idx_f(vq_file_id, "/", imembers, name_buffer, itype, err)
     
     read(name_buffer,'(I5.5)') ind
     call index2component_band(ndim**2,ind,i,j,k,l)
     
     call h5dopen_f(vq_file_id, name_buffer, grp_id, err)
     call h5dread_f(grp_id, type_r_id, vq_tmp_r, vq_dims, err)
     call h5dread_f(grp_id, type_i_id, vq_tmp_i, vq_dims, err)

     vq(i,j,k,l) = vq_tmp_r(iq)+ci*vq_tmp_i(iq)
      
     call h5dclose_f(grp_id, err)
  enddo
  call h5fclose_f(vq_file_id, err)

  ! v = 0.d0
  ! i2 = 0
  ! do l=1,ndim
  !    do j=1,ndim
  !       i2 = i2+1
  !       i1 = 0
  !       do i=1,ndim
  !          do k=1,ndim
  !             i1 = i1+1
  !             v(i1,i2) = vq(i,j,k,l)
  !          enddo
  !       enddo
  !    enddo
  ! enddo

end subroutine read_vq
!========================================================================================================


!========================================================================================================
subroutine read_u(u_tmp, filename_umatrix)
  implicit none
  real(kind=8) :: u_tmp(ndim,ndim,ndim,ndim), u_tilde_tmp(ndim,ndim,ndim,ndim)
  real(kind=8) :: u_value
  ! complex(kind=8), intent(out) :: u(ndim**2, ndim**2), u_tilde(ndim**2, ndim**2)
  ! complex(kind=8), intent(out) :: u_tilde(ndim**2, ndim**2)
  integer :: n,i,j,k,l,i1,i2
  character(len=80) :: filename_umatrix

  open(21,file=trim(filename_umatrix),status='old')
  read(21,*)
  do n=1,ndim**4
     read(21,*) i, j, k, l, u_value
     u_tmp(i,j,k,l) = u_value
     u_tilde_tmp(i,j,l,k) = u_value ! utilde for ph transverse channel
  enddo
  close(21)

  ! !go into compound index:
  ! u = 0.d0
  ! u_tilde = 0.d0
  ! i2 = 0
  ! do l=1,ndim
  !    do j=1,ndim
  !       i2 = i2+1
  !       i1 = 0
  !       do i=1,ndim
  !          do k=1,ndim
  !             i1 = i1+1
  !             u(i1,i2) = u_tmp(i,j,k,l)
  !             u_tilde(i1,i2) = u_tilde_tmp(i,j,k,l)
  !          enddo
  !       enddo
  !    enddo
  ! enddo

end subroutine read_u
!======================================================================================================


!subroutine read_v_r(v_r,r_data)
!  implicit none
!  real(kind=8) v_r(:,:,:)!(ndim**2,ndim**2,nr)
!  real(kind=8) r_data(3,nr),v_r_real(ndim**2,ndim**2)
!  integer :: nr_file,ir,i,j,nd

!  open(unit=2,file=filename_vr)
!  read(2,*) nr_file,nd,a,b,c
!  if (nr_file .ne. nr) then
!    write(*,*) 'V(r) file says there are',nr_file,'r points. '
!    write(*,*) 'Please adapt config file.'
!    stop
!  end if
!  if (nd .ne. ndim) then
!    write(*,*) 'V(r) file says there are',nd,'orbitals. '
!    write(*,*) 'Please adapt config file.'
!    stop
!  end if

!  do ir=1,nr
!    read(2,*) (r_data(i,ir),i=1,3)
!! TODO: correctly read multi-band components and go to compound index.
!    do i=1,nd**2
!       read(2,*) (v_r(i,j,ir),j=1,nd**2)
!    enddo
!  enddo

!  close(2)

!end subroutine read_v_r

!subroutine get_vq(v,q,v_r,r_data)
!  implicit none
!  complex(kind=8),intent(out) :: v(ndim**2,ndim**2)
!  real(kind=8),intent(in) :: q(3)
!  real(kind=8) :: v_r(:,:,:),r_data(:,:)
!  integer :: i

!  v=cmplx(0.d0,0.d0,kind=8)
!  if (nr.eq.0) then
!    v=u
!  end if

!  do i=1,nr
!    v = v + v_r(:,:,i)*exp(2.d0*pi*ci*(r_data(1,i)*q(1)+r_data(2,i)*q(2)+r_data(3,i)*q(3)))
!  end do

!end subroutine get_vq

end module vq_module



subroutine component2index_band(Nbands, ind, b1, b2, b3, b4)
  implicit none
  integer,intent(in) :: Nbands
  integer,intent(in) :: b1, b2, b3, b4
  integer,intent(out) :: ind

  ind =  Nbands**3*(b1-1) + Nbands**2*(b2-1) + Nbands*(b3-1) + b4

end subroutine component2index_band



! converting an index into a band pattern
subroutine index2component_band(Nbands, ind, b1, b2, b3, b4)
  implicit none
  integer,intent(in) :: Nbands,ind
  integer,intent(out) :: b1, b2, b3, b4
  integer :: tmp1,tmp2,tmp3,ind_tmp
  integer :: g1,g2,g3,g4

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

