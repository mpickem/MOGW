module mpi_org
  use aux
#ifdef MPI
  use mpi!_f08 ! latest MPI Standard -- use with latest mpi and intel compiler
#endif

integer                  :: nproc, myid, master, mpierr
integer                  :: ikstart, ikend, nkthis, ncount
integer, allocatable     :: displs(:),rcounts(:)
character(len=80)        :: chmyid
complex(dp), allocatable :: mpi_cwork(:), mpi_cwork3(:,:,:)
complex(dp), allocatable :: sndbuf(:)

contains

  subroutine mpi_close()
    implicit none
#ifdef MPI
    deallocate(rcounts,displs)
    call MPI_FINALIZE( mpierr )
#endif
  end subroutine mpi_close

  subroutine mpi_env(nk)
    implicit none
    integer nk
#ifdef MPI
    integer i
    call MPI_INIT(mpierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,mpierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,mpierr)
    allocate(displs(nproc),rcounts(nproc))
    master=0
    do i=1,nproc
       call select_krange(i-1,nk,ikstart,ikend)
       displs(i)=ikstart-1
       rcounts(i)=ikend-ikstart+1
    enddo
    call select_krange(myid,nk,ikstart,ikend)
#ifdef DEBUG
    write(chmyid,'(I5.5)') myid
    open(99,file='debug_'//trim(chmyid),status='unknown',position='append')
    write(99,'(100I12)')ikstart,ikend,nk,displs(myid+1),rcounts(myid+1),nproc
    close(99)
#endif

#else
    nproc=1
    master=0
    myid=0
    ikstart=1
    ikend=nk
#endif

  end subroutine mpi_env

  subroutine select_krange(myidin,nk,ikstart,ikend)
    implicit none
    integer myidin,id,ikstart,ikend,dq,npp,nk

    id=myidin+1 !usually starting at 0                                                                                                         
    dq=nk/nproc+1
    npp=0
    if(id.gt.nproc+nk-dq*nproc) then
       npp=nproc+nk-dq*nproc
       dq=dq-1
    endif
    ikstart=npp*(dq+1)+(id-npp-1)*dq+1
    ikend=ikstart+dq-1
    ncount=ikend-ikstart+1
    return
  end subroutine select_krange

end module mpi_org
