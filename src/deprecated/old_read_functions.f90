module read_functions
  use Mglobal
  use Mindex
  use Mham, only: ndim,nkp

  implicit none
  private
  public :: read_V, read_DMFT_SE

  contains

subroutine read_V
! input / output
  real(dp)                 :: u_tmp(ndim,ndim,ndim,ndim)
  complex(dp)              :: vq(ndim,ndim,ndim,ndim)

! auxiliaries
  integer                  :: VL(ndim,ndim), VR(ndim,ndim),iv
  integer                  :: i,j,k,l,iq,ina,inb, error
  double precision         :: tmp(ndim,nw),tmp2(ndim)
  complex(dp)              :: sumq

! initialization
  V=0.d0
  Vend=0.d0
  tmp=0.d0
  tmp2=0.d0
  sumq=0.d0


  call index_(VL,VR)

! real V ... from input files
  if (flagVfile .eqv. .true.) then

  call h5open_f(error)

  call create_complex_datatype

    ! if (myid .eq. master) then
    ! only master reads and broadcasts to everyone

    ! read non-local V(q)

    if (.false.) then

      if (myid .eq. master) write(*,*) 'Reading V(q) from ', trim(filename_vq)

      do iq=1,nkp
        call read_vq(iq, vq, filename_vq)
        do l=1,ndim
        do k=1,ndim
        do j=1,ndim
        do i=1,ndim
          V(VL(i,k),VR(j,l),iq,1) = V(VL(i,k),VR(j,l),iq,1) + vq(i,j,k,l)
        enddo
        enddo
        enddo
        enddo
      enddo

      ! check wether V(q) is purely non local
      do l=1,ndim
      do k=1,ndim
      do j=1,ndim
      do i=1,ndim
        sumq = 0.d0
        do iq=1,nkp
          sumq = sumq + V(VL(i,j),VR(k,l),iq,1)
        enddo
        sumq=sumq/nkp ! normalize
        ! write(*,*) i, j, k, l, sumq
        ! make it purley non local if thats not the case
        V(VL(i,j),VR(k,l),:,1) = V(VL(i,j),VR(k,l),:,1) - sumq
      enddo
      enddo
      enddo
      enddo
    endif ! no v(q)

    ! read local U
      if (myid .eq. master) write(*,*) 'Reading U from ', trim(filename_umatrix)

      call read_u(u_tmp, filename_umatrix)
      do l=1,ndim
      do k=1,ndim
      do j=1,ndim
      do i=1,ndim
        ! ATTENTION HERE
        ! ADGA UMATRIX FORMAT IS USED
        ! LMML - spinflip -- LMLM - U' -- LLMM - double hopping
        ! HERE
        ! LMML - spinflip -- LMLM - double hopping -- LLMM - U'
        V(VL(i,k),VR(j,l),:,1) = V(VL(i,k),VR(j,l),:,1) + u_tmp(i,j,k,l) ! static hubbard term
      enddo
      enddo
      enddo
      enddo


    ! endif

  ! broadcast from master to everyone else
    ! do l=1,ndim
    ! do k=1,ndim
    ! do j=1,ndim
    ! do i=1,ndim
    !   call &
    !   mpi_bcast(V(VL(i,j),VR(k,l),:,1),nkp,mpi_double_complex,master,mpi_comm_world)
    ! enddo
    ! enddo
    ! enddo
    ! enddo
    ! deallocate(mpi_cwork)


    ! no frequency dependency in the ADGA case
    ! everyone for himself
    do i=2,nw
      V(:,:,:,i) = V(:,:,:,1)
    enddo
    Vend(:,:,:) = V(:,:,:,1)


    if (myid .eq. master) then
    open(unit=10,file=trim(outfolder)//'/viw_test.dat')
    do iq=1,nkp
      write(10,*) real(V(1,1,iq,1)), aimag(V(1,1,iq,1))
    enddo
    close(10)
    endif

! #ifdef MPI
! ! parallelization - communication
!   allocate(mpi_cwork(nkp))
!     do ina=1,ndim**2
!       do inb=1,ndim**2
!         call MPI_ALLGATHERV(V(ina,inb,ikstart:ikend,1),ncount, MPI_DOUBLE_COMPLEX ,mpi_cwork(:),rcounts(1:nproc), displs(1:nproc),MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD ,mpierr)
!         V(ina,inb,:,i)=mpi_cwork(:)
!       enddo
!     enddo
!   deallocate(mpi_cwork)
! #endif



    ! call input_V(tmp,nw,"input/Uiv-LLLL.dat",3)
    ! ! ndim ... 1= 1111, 2= 2222, 3= 3333
    ! do iv=1,nw
    ! do i=1,ndim
    !   V(VL(i,i),VR(i,i),:,iv) = tmp(i,1)  !tmp(i,iv)
    ! enddo
    ! enddo

    ! tmp=0.d0

    ! call input_V(tmp,nw,"input/Uiv-LLMM.dat",3)
    ! ! ndim ... 1= 1122, 2= 1133, 3= 2233
    ! do iv=1,nw
    !   V(VL(1,1),VR(2,2),:,iv) = tmp(1,iv)!tmp(1,1)  !  !!!!!!
    !   V(VL(1,1),VR(3,3),:,iv) = tmp(2,iv)!tmp(2,1)  !
    !   V(VL(2,2),VR(3,3),:,iv) = tmp(3,iv)!tmp(3,1)  !
    ! enddo

    ! V(VL(2,2),VR(1,1),:,:) = V(VL(1,1),VR(2,2),:,:)
    ! V(VL(3,3),VR(1,1),:,:) = V(VL(1,1),VR(3,3),:,:)
    ! V(VL(3,3),VR(2,2),:,:) = V(VL(2,2),VR(3,3),:,:)

    ! tmp=0.d0

    ! call input_V(tmp,nw,"input/Uiv-LMML.dat",3)
    ! ! ndim ... 1= 1221, 2= 1331, 3= 2332
    ! do iv=1,nw
    !   V(VL(1,2),VR(2,1),:,iv) = tmp(1,iv)!tmp(1,1)  !
    !   V(VL(1,3),VR(3,1),:,iv) = tmp(1,iv)!tmp(2,1)  !
    !   V(VL(2,3),VR(3,2),:,iv) = tmp(1,iv)!tmp(3,1)  !
    ! enddo

    ! V(VL(2,1),VR(1,2),:,:) = V(VL(1,2),VR(2,1),:,:)
    ! V(VL(1,2),VR(1,2),:,:) = V(VL(1,2),VR(2,1),:,:)
    ! V(VL(2,1),VR(2,1),:,:) = V(VL(1,2),VR(2,1),:,:)

    ! V(VL(3,1),VR(1,3),:,:) = V(VL(1,3),VR(3,1),:,:)
    ! V(VL(1,3),VR(1,3),:,:) = V(VL(1,3),VR(3,1),:,:)
    ! V(VL(3,1),VR(3,1),:,:) = V(VL(1,3),VR(3,1),:,:)

    ! V(VL(3,2),VR(2,3),:,:) = V(VL(2,3),VR(3,2),:,:)
    ! V(VL(2,3),VR(2,3),:,:) = V(VL(2,3),VR(3,2),:,:)
    ! V(VL(3,2),VR(3,2),:,:) = V(VL(2,3),VR(3,2),:,:)

    ! !################################# nw-> oo
    ! tmp2=0.d0
    ! call input_Vend(tmp2,"input/V-LLLL.dat",0)
    ! do i=1,ndim
    !   Vend(VL(i,i),VR(i,i),:) =  tmp2(i)!V(VL(i,i),VR(i,i),:,1) !
    ! enddo

    ! tmp2=0.d0
    ! call input_Vend(tmp2,"input/V-LLMM.dat",0)
    ! Vend(VL(1,1),VR(2,2),:) = tmp2(1)!V(VL(1,1),VR(2,2),:,1)!
    ! Vend(VL(1,1),VR(3,3),:) = tmp2(2)!V(VL(1,1),VR(3,3),:,1)!
    ! Vend(VL(2,2),VR(3,3),:) = tmp2(3)!V(VL(2,2),VR(3,3),:,1)!
    ! Vend(VL(2,2),VR(1,1),:) = Vend(VL(1,1),VR(2,2),:)
    ! Vend(VL(3,3),VR(1,1),:) = Vend(VL(1,1),VR(3,3),:)
    ! Vend(VL(3,3),VR(2,2),:) = Vend(VL(2,2),VR(3,3),:)

    ! tmp2=0.d0
    ! call input_Vend(tmp2,"input/V-LMML.dat",0)
    ! Vend(VL(1,2),VR(2,1),:) = tmp2(1)!V(VL(1,2),VR(2,1),:,1)!
    ! Vend(VL(1,3),VR(3,1),:) = tmp2(2)!V(VL(1,3),VR(3,1),:,1)!
    ! Vend(VL(2,3),VR(3,2),:) = tmp2(3)!V(VL(2,3),VR(3,2),:,1)!
    ! Vend(VL(2,1),VR(1,2),:) = Vend(VL(1,2),VR(2,1),:)
    ! Vend(VL(1,2),VR(1,2),:) = Vend(VL(1,2),VR(2,1),:)
    ! Vend(VL(2,1),VR(2,1),:) = Vend(VL(1,2),VR(2,1),:)
    ! Vend(VL(3,1),VR(1,3),:) = Vend(VL(1,3),VR(3,1),:)
    ! Vend(VL(1,3),VR(1,3),:) = Vend(VL(1,3),VR(3,1),:)
    ! Vend(VL(3,1),VR(3,1),:) = Vend(VL(1,3),VR(3,1),:)
    ! Vend(VL(3,2),VR(2,3),:) = Vend(VL(2,3),VR(3,2),:)
    ! Vend(VL(2,3),VR(2,3),:) = Vend(VL(2,3),VR(3,2),:)
    ! Vend(VL(3,2),VR(3,2),:) = Vend(VL(2,3),VR(3,2),:)

  else

  ! V for testing purposes ... constant U at NNNN for every iv and ikq
    do i=1,ndim
      V(VL(i,i),VR(i,i),:,:) = 1.5d0
      Vend(VL(i,i),VR(i,i),:) = 1.5d0
    enddo

  endif !flagVfile


  return
end subroutine read_V

end module
