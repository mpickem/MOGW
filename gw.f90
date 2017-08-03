program gw
  use computation_functions
  use read_functions
  use index_reference
  use hamiltonian_module
  use lapack_module
  use mpi_org

  implicit none

! parameters of the calculation
  real(dp) :: mu
  integer :: L
  complex(kind=8) :: ntot

! flags
  logical :: flagN,flagSE,flagVfile

! variables (Greensfunctions, interactions, self-energy, ...)
  complex(kind=8), allocatable :: Giw(:,:,:,:),Gconv(:,:,:)
  complex(kind=8), allocatable :: ncur(:,:),n_min(:,:),n_max(:,:),n_mid(:,:)
  real(dp), allocatable :: Gtau(:,:,:,:)
  complex(kind=8), allocatable :: P(:,:,:,:)
  complex(kind=8), allocatable :: V(:,:,:,:),Vend(:,:,:)
  complex(kind=8), allocatable :: W(:,:,:,:)
  complex(kind=8), allocatable :: SE_old(:,:,:,:), SE_new(:,:,:,:)

! auxiliaries
  complex(kind=8) :: n,nmin,nmax,nmid,dif,maxdif
  real(dp) :: mu_min,mu_max,mu_mid,tmp,dummy
  real(dp) :: tstart, tend  ! mpi times
  integer :: i,j,ikp,iw,cyc
  character(len=100) :: cmd_arg


  ! read command line argument -> file name of config file
  if (iargc() .ne. 1) then
    write(*,*) 'The program has to be executed with exactly one argument. (Name of config file)'
    stop
  end if

  call getarg(1,cmd_arg)
  write(*,*) 'Reading config: ', trim(cmd_arg)

! parameters
  open(unit=10,file=trim(cmd_arg), &
    status='old', action='read', position='rewind')
      read(10,*)
      read(10,*)
      read(10,*) nw    ! number of Matsubara frequencies
      read(10,*) L     ! number of time discretization point in [0,beta) 
      read(10,*) beta    ! inverse temperature
      read(10,*) mu    ! chemical potential
      read(10,*) dummy   ! total number of particles
      ntot=dummy ! cant read complex number properly
      read(10,*)
      read(10,*) flagN   ! flag for mu-cycle
      read(10,*) flagSE    ! flag for SE-cycle
      read(10,*) flagVfile   ! flag for V-array from files
      read(10,*)
      read(10,*) outfolder
      read(10,*)
      read(10,*) filename_umatrix
      read(10,*) filename_vq
      read(10,*)
      read(10,*) filename_dmft
    close(unit=10)

  call system("mkdir -p "//trim(outfolder))
  ! no error if existing

  if (flagN .eqv. .true.) mu=0.d0 !symmetric bisection start around 0

  call read_hamiltonian
  call read_bzindices
  call mpi_env(nkp)

  if(myid.eq.master) tstart=mpi_wtime()
! calls mpi environment ... associates different cores with start and enpoints of ikp .. ikstart/ikend

! program introduction
  if(myid.eq.master) then
    write(*,*) '# ----------------------------------------------------------------------'
    write(*,*) '# parameters of GW calculation'
    write(*,*) '# ----------------------------------------------------------------------'
    write(*,*) '# number of matsubara frequencies      | nw:    ', nw
    write(*,*) '# inverse temperature                  | beta:  ', beta
    write(*,*) '# chemical potential                   | mu:    ', mu
    write(*,*) '#                                      |'
    write(*,*) '# fix mu mode                          | flagN: ', flagN
    if (flagN .eqv. .true.) then
    write(*,*) '# number of particles                  | n_tot: ', real(ntot)
    endif
    write(*,*) '# SE self-consistency                  | flagSE:', flagSE
    write(*,*) '# ----------------------------------------------------------------------'
    write(*,*) '# ----------------------------------------------------------------------'
    write(*,*) '# number of orbitals                   | ndim:    ', ndim
    write(*,*) '# number of k-points                   | nkp:     ', nkp
    write(*,*) '# ----------------------------------------------------------------------'
  endif


  allocate(Giw(ndim,ndim,nkp,2*nw))
! allocate(Gtau(ndim,ndim,nkp,nw))
  allocate(SE_old(ndim,ndim,nkp,2*nw),SE_new(ndim,ndim,nkp,2*nw))


! setting up variables for first runthrough

  ! for first runthrough -> GDMFT(ik)
  ! call read_DMFT_SE(SE_new,mu,filename_dmft)
  ! if(myid.eq.master) write(*,*) 'mu: ', mu

  ! for first runthrough -> G0(ik)
  SE_new = 0.d0
  cyc=0     ! counting cycles

  !######################################
  !############# SECYCLE ################
  !######################################

  SEcycle: do
    SE_old = SE_new ! Saving results - used to calculate new Giw
    deallocate(SE_new)  ! resource management
    cyc=cyc+1   ! increasing cycle number by 1

    if(myid.eq.master) write(*,*) 'Iteration number: ', cyc

    !######################################
    !############ BISECTION ###############
    !######################################

    allocate(ncur(ndim,ndim),n_min(ndim,ndim),n_mid(ndim,ndim),n_max(ndim,ndim))
    if(flagN .eqv. .true.) then

      if(myid.eq.master) write(*,*) 'starting bisection'

        mu_max=mu+20.d0         ! Bisection start and end point
        mu_min=mu-20.d0         ! symmetric around mu of previous cycle (first cycle - mu=0)

        call compute_Giw(mu_min,Giw,SE_old) ! calculate Giw with mu_min
        call compute_n(n_min,nmin,Giw)    ! calculate ncur_min with Giw

        call compute_Giw(mu_max,Giw,SE_old) ! calculate Giw with mu_max
        call compute_n(n_max,nmax,Giw)    ! calculate ncur_max with Giw


      mucycle: do

       mu_mid = (mu_max + mu_min)/2.d0      !calculate mu_mid

       call compute_Giw(mu_mid,Giw,SE_old)  ! calculate Giw with mu_max+mu_min / 2
       call compute_n(n_mid,nmid,Giw)   ! calculate ncur_mid with Giw

       ! compare signs and choose mu_mid as new mu_min or mu_max
       if( (real(ntot-nmin) .gt. 0.d0) .and. (real(ntot-nmax) .lt. 0.d0) ) then ! ncur is increasing with mu
          if (abs(real(ntot-nmid)) .le. 1d-9 ) then
            mu = mu_mid
            ncur = n_mid    ! array
            n = nmid    ! trace
            if(myid.eq.master) write(*,*) "Bisection yielded mu = ", mu
            if(myid.eq.master) write(*,*) "Bisection yielded n = ", real(n)
            exit
          else if (real(ntot-nmid) .gt. 0.d0)  then
            mu_min = mu_mid
            n_min = n_mid
            nmin = nmid
          else
            mu_max = mu_mid
            n_max = n_mid
            nmax = nmid
          endif
       endif

      enddo mucycle
    endif   !flagN

    !######################################
    !######### END-BISECTION ##############
    !######################################

  ! if flagN = false ... Calculate with given mu
     if (flagN .eqv. .false.) then
       call compute_Giw(mu,Giw,SE_old)   ! calculate Giw with given mu
       call compute_n(ncur,n,Giw)    ! calculate ncur with Giw
       if(myid.eq.master) write(*,*) "Calculation yielded n = ", real(n)
     endif


     deallocate(ncur,n_min,n_max,n_mid)

  ! call Giw2Gtau(nw,L,beta,Giw,Gtau)
  ! Computation of Greenfunction (tau)

    if (cyc .eq. 1) then ! only at the first run
    ! with resource management
      allocate(Gconv(ndim,nkp,2*nw))
      call compute_Gconv(mu,Gconv)    ! Convergence Term for P with updated mu
      allocate(P(ndim**2,ndim**2,nkp,nw))
      call compute_P(mu,Giw,Gconv,P)  ! Computation of Polarization (iv,iq) with G*G
      deallocate(Gconv)
      allocate(V(ndim**2,ndim**2,nkp,nw),Vend(ndim**2,ndim**2,nkp))
      call read_V(V,Vend,flagVfile)   ! reading input files to creating V matrix
      allocate(W(ndim**2,ndim**2,nkp,nw))
      call compute_W(P,V,W)     ! Computation of screened interaction
      deallocate(P,V)
      allocate(SE_new(ndim,ndim,nkp,2*nw))    
      call compute_SE(Giw,W,Vend,SE_new)  ! Computation of Self-Energy
      deallocate(W,Vend)
    endif
   

  ! Difference between old and new Self-Energy - self consistency
    tmp=0.d0
    if(flagSE .eqv. .true.) then
      maxdif=0.d0
      do iw=1,nw
      do ikp=1,nkp
      do i=1,ndim
      do j=1,ndim
        dif=SE_old(i,j,ikp,iw) - SE_new(i,j,ikp,iw)
        if( (real(dif)**2.d0 + aimag(dif)**2.d0) .gt. (real(maxdif)**2.d0 + aimag(maxdif)**2.d0) ) maxdif=dif
      enddo
      enddo
      enddo
      enddo
     
      tmp = real(maxdif)**2.d0 + aimag(maxdif)**2.d0
      if(myid.eq.master) then
        write(*,*) '# ----------------------------------------------------------------------'
        write(*,*) "# max difference between this and previous cycle:", tmp
        write(*,*) '# ----------------------------------------------------------------------'
      endif
    endif

    ! for one shot calculations with updated G
    if (cyc .ge. 2) then
      if ( (tmp .lt. 1.d-8 ) .or. (flagSE .eqv. .false.)) then
      exit
      endif
    endif

  enddo SEcycle

  !######################################
  !########### END-SECycle ##############
  !######################################


  if(myid.eq.master) then
    write(*,*) '# number of cycles:', cyc
    tend=mpi_wtime()
    write(*,*) 'Time needed [s]: ',tend-tstart
    write(*,*) 'Closing mpi environment'
  endif

  call mpi_close()

  deallocate(Giw,SE_old)
end program gw
