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

   ! this is from the ADGA code
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

end module hdf5_module
