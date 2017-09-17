module hdf5_module
  use aux
  use hdf5
  use hamiltonian_module ! for nkp, ndim
  implicit none
  integer                        :: iwfmax, iwbmax
  integer(hid_t)                 :: plist_id
  integer(hid_t)                 :: dspace_id
  integer(hid_t)                 :: compound_id, type_r_id, type_i_id
  integer(size_t)                :: compound_size, type_sized
  integer(hsize_t), dimension(2) :: dims
  real(dp), allocatable          :: tmp_r_1(:,:), tmp_i_1(:,:), tmp_err_1(:,:)
  integer                        :: err

 contains

!=====================================================================================
   subroutine create_complex_datatype

     integer(size_t), parameter :: zero = 0

     ! Set dataset transfer property to preserve partially initialized fields during write/read to/from dataset with compound datatype (necessary?)
     call h5pcreate_f(h5p_dataset_xfer_f, plist_id, err)
     call h5pset_preserve_f(plist_id, .true., err)

     ! create compound datatype for complex arrays:
     call h5tget_size_f(h5t_native_double, type_sized, err)
     compound_size = 2*type_sized
     call h5tcreate_f(h5t_compound_f, compound_size, compound_id, err)
     call h5tinsert_f(compound_id, "r", zero, h5t_native_double, err)
     call h5tinsert_f(compound_id, "i", type_sized, h5t_native_double, err)

     !complex type to write real and imaginary individually:
     call h5tcreate_f(h5t_compound_f, type_sized, type_r_id, err)
     call h5tinsert_f(type_r_id, "r", zero, h5t_native_double, err)
     call h5tcreate_f(h5t_compound_f, type_sized, type_i_id, err)
     call h5tinsert_f(type_i_id, "i", zero, h5t_native_double, err)

   end subroutine create_complex_datatype
!=====================================================================================

!=====================================================================================
   subroutine read_axes(file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
     real(dp), allocatable :: iwb_array(:),iwf_array(:)
     integer(hsize_t), dimension(1), intent(out) :: dim_iwb, dim_iwf
     integer(hsize_t), dimension(1) :: dim_iwf_max, dim_iwb_max
     integer :: err
     integer(hid_t) :: axes_id,iwb_id,iwf_id
     integer(hid_t), intent(out) :: dspace_iwb_id, dspace_iwf_id
     integer(hid_t) :: file_id
 
     ! read fermionic Matsubara frequencies iwf:
     call h5dopen_f(file_id, ".axes/iwf-g4", iwf_id, err)
     call h5dget_space_f(iwf_id, dspace_iwf_id, err)
     call h5sget_simple_extent_dims_f(dspace_iwf_id, dim_iwf, dim_iwf_max, err)
     iwfmax = dim_iwf(1)/2
     allocate(iwf_array(-iwfmax:iwfmax-1))
     call h5dread_f(iwf_id, h5t_native_double, iwf_array, dim_iwf, err)
     call h5dclose_f(iwf_id, err)

     ! read bosonic Matsubara frequencies iwf:
     call h5dopen_f(file_id, ".axes/iwb-g4", iwb_id, err)
     call h5dget_space_f(iwb_id, dspace_iwb_id, err)
     call h5sget_simple_extent_dims_f(dspace_iwb_id, dim_iwb, dim_iwb_max, err)
     iwbmax = dim_iwb(1)/2
     allocate(iwb_array(-iwbmax:iwbmax))
     call h5dread_f(iwb_id, h5t_native_double, iwb_array, dim_iwb, err)
     call h5dclose_f(iwb_id, err)
   
   end subroutine read_axes


   subroutine write_axes(file_id, iwb_array, iwf_array, dspace_iwb_id, dspace_iwf_id, dim_iwb, dim_iwf)
     integer(hid_t)                :: file_id
     real(dp), intent(in)          :: iwb_array(-iwbmax:iwbmax), iwf_array(-iwfmax:iwfmax-1)
     integer(hsize_t),dimension(1) :: dim_iwb, dim_iwf
     integer                       :: err
     integer(hid_t)                :: axes_id, iwb_id, iwf_id
     integer(hid_t), intent(in)    :: dspace_iwb_id, dspace_iwf_id

     !write Matsubara frequency axes:
     call h5gcreate_f(file_id, ".axes", axes_id, err)
     call h5dcreate_f(axes_id, "iwb-g4", h5t_native_double, dspace_iwb_id, iwb_id, err)
     call h5dcreate_f(axes_id, "iwf-g4", h5t_native_double, dspace_iwf_id, iwf_id, err)

     call h5dwrite_f(iwb_id, h5t_native_double, iwb_array, dim_iwb, err) 
     call h5dwrite_f(iwf_id, h5t_native_double, iwf_array, dim_iwf, err)
     
     call h5dclose_f(iwb_id, err)
     call h5dclose_f(iwf_id, err)
     call h5gclose_f(axes_id, err)

   end subroutine write_axes
!============================================================================================

!===========================================================================================
   subroutine create_channels(file_id)
     implicit none

     integer           :: iwb,err
     integer(hid_t)    :: grp_dens_id,grp_magn_id,iw_magn_id,iw_dens_id
     character(len=20) :: name_buffer
     integer(hid_t)    :: file_id
     
     !create dens and magn groups:
     call h5gcreate_f(file_id, "dens", grp_dens_id, err)
     call h5gcreate_f(file_id, "magn", grp_magn_id, err)

     do iwb=0,2*iwbmax

        write(name_buffer, '(I5.5)') iwb
        call h5gcreate_f(grp_dens_id, name_buffer, iw_magn_id, err)
        call h5gcreate_f(grp_magn_id, name_buffer, iw_dens_id, err)
        call h5gclose_f(iw_magn_id, err)
        call h5gclose_f(iw_dens_id, err)

     enddo

     call h5gclose_f(grp_dens_id, err)
     call h5gclose_f(grp_magn_id, err)

     return
   end subroutine create_channels
!=====================================================================================

!=====================================================================================
    subroutine create_component(file_id, ichannel, iwb, ind_orb)
      implicit none

      integer, intent(in) :: ichannel !1=magn, 2=dens
      integer, intent(in) :: iwb, ind_orb
      character(len=20)   :: grpname
      integer(hid_t)      :: grp_id, dset_id, dset_err_id
      integer(hid_t)      :: file_id

      if (ichannel==1) then
         write(grpname, '("magn/",I5.5,"/",i5.5)') iwb, ind_orb
      else
         write(grpname, '("dens/",I5.5,"/",i5.5)') iwb, ind_orb
      endif

      call h5gcreate_f(file_id, grpname, grp_id, err)
      call h5dcreate_f(grp_id, "value", compound_id, dspace_id, dset_id, err)
      call h5dcreate_f(grp_id, "error", h5t_native_double, dspace_id, dset_err_id, err)

      tmp_r_1 = 0.d0
      tmp_i_1 = 0.d0
      tmp_err_1 = 0.d0

      call h5dwrite_f(dset_id, type_r_id, tmp_r_1, dims, err)
      call h5dwrite_f(dset_id, type_i_id, tmp_i_1, dims, err)
      call h5dwrite_f(dset_err_id, h5t_native_double, tmp_err_1, dims, err)

      call h5dclose_f(dset_err_id, err)
      call h5dclose_f(dset_id, err)
      call h5gclose_f(grp_id, err)

    end subroutine create_component
!====================================================================================

!====================================================================================
   subroutine add_to_component(file_id, ichannel, iwb, ind_orb, g4iw_r, g4iw_i, g4err)
     implicit none
     integer, intent(in) :: ichannel, iwb, ind_orb
     real(dp),intent(in) :: g4iw_r(2*iwfmax, 2*iwfmax, 2*iwbmax+1)
     real(dp),intent(in) :: g4iw_i(2*iwfmax, 2*iwfmax, 2*iwbmax+1)
     real(dp),intent(in) :: g4err(2*iwfmax, 2*iwfmax, 2*iwbmax+1)
     character(len=20)   :: grpname
     integer(hid_t)      :: grp_id, dset_id, dset_err_id
     integer(hid_t)      :: file_id

     if (ichannel==1) then
        write(grpname, '(A5,(I5.5),A1,(I5.5))') "magn/", iwb, "/", ind_orb
     else
        write(grpname, '(A5,(I5.5),A1,(I5.5))') "dens/", iwb, "/", ind_orb
     endif
        
     call h5gopen_f(file_id, grpname, grp_id, err)
     call h5dopen_f(grp_id, "value", dset_id, err)
     call h5dopen_f(grp_id, "error", dset_err_id, err)
     
     call h5dread_f(dset_id, type_r_id, tmp_r_1, dims, err)
     call h5dread_f(dset_id, type_i_id, tmp_i_1, dims, err)
     call h5dread_f(dset_err_id, h5t_native_double, tmp_err_1, dims, err)

     tmp_r_1(:,:) = g4iw_r(:,:,iwb+1)+tmp_r_1(:,:)
     tmp_i_1(:,:) = g4iw_i(:,:,iwb+1)+tmp_i_1(:,:)
     tmp_err_1(:,:) = sqrt(g4err(:,:,iwb+1)**2+tmp_err_1(:,:)**2)

     call h5dwrite_f(dset_id, type_r_id, tmp_r_1, dims, err)
     call h5dwrite_f(dset_id, type_i_id, tmp_i_1, dims, err)
     call h5dwrite_f(dset_err_id, h5t_native_double, tmp_err_1, dims, err)

     call h5dclose_f(dset_err_id, err)
     call h5dclose_f(dset_id, err)
     call h5gclose_f(grp_id, err)

 end subroutine add_to_component
!===========================================================================================

end module hdf5_module
