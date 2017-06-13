module index_reference
  use aux
  use hamiltonian_module, only : ndim
  implicit none
  private
  public :: index_, inverse_l, inverse_r

  contains

subroutine index_(x,y)  ! x - left compound index --- y - right compound index
! this can also be done via a function written like below - however it is faster to do one loop since every index is used anyways
  integer, intent(out) :: x(ndim,ndim),y(ndim,ndim)  ! x = 1 ... ndim**2 = y ---- Polarization arrays have to be changed
  integer :: cnt,i,j

  cnt=1

  do i=1,ndim
  do j=1,ndim
    x(i,j) = cnt ! 11, 12, 13, ..., 21, 22, 23, ...
    y(j,i) = cnt ! 11, 21, 31, ..., 12, 22, 32, ...
    cnt = cnt+1
  enddo
  enddo

end subroutine index_

integer function inverse_l(i) ! i -> L (from LL_) --- j -> L_ (from L_L)
  integer,intent(in) :: i
  inverse_l=(i-1)/ndim+1 ! integer division
  return
end function inverse_l

integer function inverse_r(i) ! i -> L_ (from LL_) --- j -> L (from L_L)
  integer,intent(in) :: i
  inverse_r=mod(i-1,ndim)+1
  return
end function inverse_r

end module index_reference

!##################################################

module computation_functions
  use aux
  use index_reference
  use hamiltonian_module ! contains routines for reading the Hamiltonian
  use lapack_module      ! contains an easy interface for computing the inverse of a matrix
  use mpi_org    ! parallelization
  use hdf5_module
  use vq_module
  implicit none
  character(len=80) :: filen

  private :: filen
  public :: compute_Giw, compute_Gconv, compute_P, compute_W, compute_SE, compute_n, compute_V

contains

subroutine compute_Giw(mu,Giw,SE) 
! input/output
  real(dp), intent(in) :: mu
  complex(kind=8), intent(out) :: Giw(ndim,ndim,nkp,2*nw)
  complex(kind=8), intent(inout) :: SE(ndim,ndim,nkp,2*nw)  !inout because of Continuation

! auxiliaries
  integer :: ikp,iw,i,j,ina,inb
  complex(kind=8), allocatable :: cmat(:,:)

! initialization
  Giw = 0.d0
! continuation of Se from nw+1 to 2*nw - required to calculate P
  do iw=nw+1,2*nw
    SE(:,:,:,iw) = SE(:,:,:,nw)*real(2*(nw-1)+1,kind=8)/real(2*(iw-1)+1,kind=8)
  enddo

! compute G(k,iw) = [ iw + mu - H(k) - SE(k,iw)]^{-1}
  allocate(cmat(ndim,ndim))
  do iw=1,2*nw ! Matsubara frequencies
    do ikp=ikstart,ikend !1,nkp ! k-points
      cmat=0.d0
      do i=1,ndim ! orbitals
         cmat(i,i) = ci * real(2*(iw-1)+1,kind=8)*pi/beta + mu
      enddo
      ! iw + mu - H(k) - SE :
      cmat(:,:)=cmat(:,:)-h(:,:,ikp)-SE(:,:,ikp,iw)
      ! [iw + mu - H(k) - Sigma(w,k) ]^{-1} :
      call inverse_matrix( cmat ) ! inverts complex square matrix
      Giw(:,:,ikp,iw)=cmat
    enddo ! ikp
  enddo !iw
  deallocate(cmat)

! parallelization - communication
#ifdef MPI
  allocate(mpi_cwork(nkp))
  do i=1,2*nw
    do ina=1,ndim
      do inb=1,ndim
        call MPI_ALLGATHERV(Giw(ina,inb,ikstart:ikend,i),ncount, MPI_DOUBLE_COMPLEX ,mpi_cwork(:),rcounts(1:nproc), displs(1:nproc),MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD ,mpierr)
        ! call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
        Giw(ina,inb,:,i)=mpi_cwork(:)
      enddo
    enddo
  enddo
  deallocate(mpi_cwork)
#endif

! output - should go into io.f90
  if(myid.eq.master) then
    filen=trim(outfolder)//"/Gw.dat"
    call local_output_Matsub(Giw,1,filen) ! this routine is in io.f90
    filen=trim(outfolder)//"/Gw_diag.dat"
    call local_output_Matsub_diagonal(Giw,1,filen) ! this routine is in io.f90

! output - should go into io.f90
    
    filen=trim(outfolder)//"/Gkw.dat"
    open(unit=10,file=filen)
      write(10,*) "## ikp, iw, iomega, RE[G_ii(ikp,iw)], IM(G_ii(ikp,iw)]"
      write(10,*) "##"
      do ikp=1,nkp
      do iw=1, min(250,nw)
        write(10,'(I7,I7,F20.8)',advance='no') ikp,iw,pi/beta*real(2*(iw-1)+1,kind=8)
        do i=1,ndim
          write(10,'(E23.8)',advance='no') real(Giw(i,i,ikp,iw))
          write(10,'(E23.8)',advance='no') aimag(Giw(i,i,ikp,iw))
        enddo
        write(10,*) ! line break
      enddo
      enddo
    close(10)

    filen=trim(outfolder)//"/Hk.dat"
    open(unit=11,file=filen)
      write(11,*) "## ikp, RE[H_ii(ikp)], IM(H_ii(ikp)]"
      write(11,*) "##"
      do ikp=1,nkp
        write(11,'(I7)',advance='no') ikp
        do i=1,ndim
          write(11,'(E23.8)',advance='no') real(h(i,i,ikp))
          write(11,'(E23.8)',advance='no') aimag(h(i,i,ikp))
        enddo
        write(11,*) ! line break
      enddo
    close(11)

  endif ! master

  return
end subroutine compute_Giw

subroutine compute_Gconv(mu,Gconv)
! input/output
  real(dp), intent(in) :: mu
  complex(kind=8), intent(out) :: Gconv(ndim,nkp,2*nw) 

! auxiliaries
  integer :: ikp,iw,i,ina
  complex(kind=8), allocatable :: cmat(:)

! compute Gtilde_ii(k,iw) = [ iw + mu - real(H(k))]^{-1}
  allocate(cmat(ndim))
  do iw=1,2*nw ! Matsubara frequencies
    do ikp=ikstart,ikend ! 1,nkp ! k-points
      cmat=0.d0
        do i=1,ndim ! orbitals
          !iw+mu-real(h(k))
           cmat(i) = ci * real(2*(iw-1)+1,kind=8)*pi/beta + mu - real(h(i,i,ikp))
        enddo
      cmat = 1.d0 / cmat
      Gconv(:,ikp,iw) = cmat
    enddo ! ikp
  enddo !iw
  deallocate(cmat)

#ifdef MPI
! parallelization - communication
  allocate(mpi_cwork(nkp))
  do i=1,2*nw
    do ina=1,ndim
      call MPI_ALLGATHERV(Gconv(ina,ikstart:ikend,i),ncount, MPI_DOUBLE_COMPLEX ,mpi_cwork(:),rcounts(1:nproc), displs(1:nproc),MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD ,mpierr)
      ! call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
      Gconv(ina,:,i)=mpi_cwork(:)
    enddo
  enddo
  deallocate(mpi_cwork)
#endif

  return
end subroutine compute_Gconv

!########################################

subroutine compute_P(mu,Giw,Gconv,P)
!input/output
  real(dp), intent(in) :: mu
  complex(kind=8), intent(in) :: Giw(ndim,ndim,nkp,2*nw)
  complex(kind=8), intent(in) :: Gconv(ndim,nkp,2*nw)
  complex(kind=8), intent(out) :: P(ndim**2,ndim**2,nkp,nw)

!auxiliaries
  integer :: iv,ikq,iw,ikp,i,j,k,l,ina,inb
  complex(kind=8) :: G1,G2,ctmp
  integer :: IL(ndim,ndim),IR(ndim,ndim)

!initialization
  P = 0.d0
  G1 = 0.d0
  G2 = 0.d0
  ctmp = 0.d0

!compound indizes
  call index_(IL,IR)

  if(myid.eq.master) write(*,*) 'Computing P(q,iv)'

!multiorbital calculation
!compute P_ilkj(ikq,iv) = 2 / beta * G_ij(ikp,iw)* G_kl(ikp+ikq,iw+iv)
  do iv=1,nw
  do iw=1,(iv-1)
  do ikq=ikstart,ikend ! 1,nkp
  do ikp=1,nkp
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim
  do i=1,ndim
    G1=Giw(i,j,ikp,iw)
    G2=Giw(k,l,ikpq(ikp,ikq),iw+iv-1)
    ctmp=G1*G2
    ctmp = ctmp + conjg(Giw(j,i,ikp,iw)) * Giw(k,l,ikpq(ikp,ikq),iv-iw)   
    P(IL(i,l),IR(k,j),ikq,iv) = P(IL(i,l),IR(k,j),ikq,iv) + ctmp * wtkp(ikp) 
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  do iv=1,nw
  do iw=iv,nw
  do ikq=ikstart,ikend ! 1,nkp
  do ikp=1,nkp
  do j=1,ndim
  do k=1,ndim
  do l=1,ndim
  do i=1,ndim
    G1=Giw(i,j,ikp,iw)
    G1=Giw(i,j,ikp,iw)
    G2=Giw(k,l,ikpq(ikp,ikq),iw+iv-1)
    ctmp=G1*G2
    ctmp = ctmp + conjg(Giw(j,i,ikp,iw)) * conjg(Giw(l,k,ikpq(ikp,ikq),iw-iv+1))
    P(IL(i,l),IR(k,j),ikq,iv) = P(IL(i,l),IR(k,j),ikq,iv) + ctmp * wtkp(ikp)
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo


!#####################################
!for convergence reasons:
!P_ijji(ikq,iv) = P_nqqnn - 2/beta * Gtilde_ii(ikp,iw) * Gtilde_jj(ikp+ikq,iw+iv)
!##################################### 


  do iv=1,nw
  do iw=1,(iv-1)
  do ikq=ikstart,ikend ! 1,nkp
  do ikp=1,nkp
  do i=1,ndim
  do j=1,ndim
    G1=Gconv(i,ikp,iw)
    G2=Gconv(j,ikpq(ikp,ikq),iw+iv-1)
    ctmp=G1*G2
    ctmp = ctmp + conjg(G1) * Gconv(j,ikpq(ikp,ikq),iv-iw)
    P(IL(i,j),IR(j,i),ikq,iv) = P(IL(i,j),IR(j,i),ikq,iv) - ctmp * wtkp(ikp)
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  do iv=1,nw
  do iw=iv,nw
  do ikq=ikstart,ikend ! 1,nkp
  do ikp=1,nkp
  do i=1,ndim
  do j=1,ndim
    G1=Gconv(i,ikp,iw)
    G2=Gconv(j,ikpq(ikp,ikq),iw+iv-1)
    ctmp=G1*G2
    ctmp = ctmp + conjg(G1) * conjg(Gconv(j,ikpq(ikp,ikq),iw-iv+1))
    P(IL(i,j),IR(j,i),ikq,iv) = P(IL(i,j),IR(j,i),ikq,iv) - ctmp * wtkp(ikp)
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  P = P / beta  ! common prefactor, no factor of 2 for spin summation, since sum_ikp wtkp(ikp) = 2

!#######################################
!adding the previous term in form of an analytical solution
!#######################################

  do ikq=ikstart,ikend ! 1,nkp
  do ikp=1,nkp
  do i=1,ndim
  do j=1,ndim
    ! special case for iv=1 (no shift) ----- df/dw at w = e_q(k)-mu
    if( abs(real(h(i,i,ikpq(ikp,ikq))-h(j,j,ikp))) .lt. 1.d-8) then
    P(IL(i,j),IR(j,i),ikq,1) = P(IL(i,j),IR(j,i),ikq,1) - beta*exp(beta*(real(h(j,j,ikp))-mu)) * wtkp(ikp) &
              / (exp(beta*(real(h(j,j,ikp))-mu))+1.d0)**2.d0
    else    ! same formula as below with iv=1
    P(IL(i,j),IR(j,i),ikq,1) = P(IL(i,j),IR(j,i),ikq,1) + ((1.d0/(exp((real(h(i,i,ikpq(ikp,ikq)))-mu)*beta)+1.d0)) &
        - (1.d0/(exp((real(h(j,j,ikp))-mu)*beta)+1.d0))) &
        * (1.d0/real(h(i,i,ikpq(ikp,ikq))-h(j,j,ikp))) * wtkp(ikp)
    endif
  enddo
  enddo
  enddo
  enddo

  do iv=2,nw  ! + f(e_n(k+q)-mu)-f(e_q(k)-mu) * 1/(-iv + en(k+q) - eq(k)) ---- iv= i*2pi*n /beta
  do ikp=1,nkp
  do ikq=ikstart,ikend ! 1,nkp
  do i=1,ndim
  do j=1,ndim
    P(IL(i,j),IR(j,i),ikq,iv) = P(IL(i,j),IR(j,i),ikq,iv) + ((1.d0/(exp((real(h(i,i,ikpq(ikp,ikq)))-mu)*beta)+1.d0)) &
        - (1.d0/(exp((real(h(j,j,ikp))-mu)*beta)+1.d0))) &
        * (1.d0/(real(h(i,i,ikpq(ikp,ikq))-h(j,j,ikp))-ci*2.d0*(iv-1)*pi/beta)) * wtkp(ikp)
  enddo !iv
  enddo !j
  enddo !i
  enddo !ikp
  enddo !ikq

#ifdef MPI
! parallelization - communication
  allocate(mpi_cwork(nkp))
  do i=1,nw
    do ina=1,ndim**2
      do inb=1,ndim**2
        call MPI_ALLGATHERV(P(ina,inb,ikstart:ikend,i),ncount, MPI_DOUBLE_COMPLEX ,mpi_cwork(:),rcounts(1:nproc), displs(1:nproc),MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD ,mpierr)
        ! call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
        P(ina,inb,:,i)=mpi_cwork(:)
      enddo
    enddo
  enddo
  deallocate(mpi_cwork)
#endif

  if(myid.eq.master) then
    filen=trim(outfolder)//"/Pv.dat"
    write(*,*) 'output of Ploc(iv) to file Pv.dat'
    call local_output_Compound(P,0,filen) ! 0 ... bosonic
  endif
! bosonic frequencies because iv=1 is equivalent to no shift etc.

  return
end subroutine compute_P

!##########################################

subroutine Giw2Gtau(L,Giw,Gtau)
! input/output
  integer, intent(in) :: L
  complex(kind=8), intent(in) :: Giw(ndim,ndim,nkp,nw)
  real(dp), intent(out) :: Gtau(ndim,ndim,nkp,nw)

! auxiliaries
  integer :: ikp,ina,inb
  complex(kind=8), allocatable :: mat(:,:)

  Gtau = 0.d0 ! initialization

  write(*,*) 'FT G(iw) --> G(tau)'
! the following should be a subroutine: "call Giw2Gtau()"
  do ikp=1,nkp ! loop over k-points
    do ina=1,ndim
      ! diagonal elements of G:
      call invfourierhp(beta,L,nw,Giw(ina,ina,ikp,:),Gtau(ina,ina,ikp,:),0.d0) ! high precision version... 
      ! the FT routines are in four.f90
      do inb=ina+1,ndim
         ! off-diagonal elements of G:
         call invfourier(beta,L,nw,1,1,Giw(ina,inb,ikp,:),Gtau(ina,inb,ikp,:))
         call invfourier(beta,L,nw,1,1,Giw(inb,ina,ikp,:),Gtau(inb,ina,ikp,:))
         ! The second call (that for (inb,ina) ) can be avoided by exploiting the symmetry of G:
         ! Giw(ina,inb)=Giw^*(inb,ina)  --> Gtau(ina,inb,tau)=Gtau(inb,ina, -tau ) = -Gtau(inb,ina,beta-tau) ... I think ... to be checked:
         !Gtau(inb,ina,ikp,1) = Gtau(ina,inb,ikp,1) ! because there is not jump in offdiag G...
         !do iL=2,L
         !   Gtau(inb,ina,ikp,iL)= - Gtau(ina,inb,ikp,L+2-iL) 
         !enddo
       enddo
    enddo
  enddo ! ikp
! output of local part of Gloc to file
  write(*,*) 'output of G0loc(tau) to file Gt.dat'
! We have information about Gtau for tau in [0,beta)
! In order to have G(tau=beta) we can exploit properties of G:
! - diagonal elements jump by 1 (antiperiodic!)
! - offdiagonal elements are continous (no jump)
  allocate(mat(ndim,ndim))
  mat=0.d0
  do ina=1,ndim
     mat(ina,ina)=1.d0
  enddo
  call local_output_tau(Gtau,ndim,L,1,mat,trim(outfolder)//"//Gt.dat")
  deallocate(mat)

  return
end subroutine Giw2Gtau

!##########################################

subroutine compute_V(V,Vend,flagVfile)
! input / output
  logical, intent(in) :: flagVfile
  complex(kind=8), intent(out) :: V(ndim**2,ndim**2,nkp,nw),Vend(ndim**2,ndim**2,nkp)
  real(dp) :: u_tmp(ndim,ndim,ndim,ndim)
  complex(kind=8) :: vq(ndim,ndim,ndim,ndim)

! auxiliaries
  integer :: VL(ndim,ndim), VR(ndim,ndim),iv
  integer :: i,j,k,l,iq,ina,inb, error
  double precision :: tmp(ndim,nw),tmp2(ndim)
  complex(dp) :: sumq

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
    
    if (myid .eq. master) then
    ! only master reads and broadcasts to everyone

    ! read non-local V(q)
      do iq=1,nkp
        call read_vq(iq, vq, filename_vq)
        do l=1,ndim
        do k=1,ndim
        do j=1,ndim
        do i=1,ndim
          V(VL(i,j),VR(k,l),iq,1) = V(VL(i,j),VR(k,l),iq,1) + vq(i,j,k,l)
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

    ! read local U
      call read_u(u_tmp, filename_umatrix)
      do l=1,ndim
      do k=1,ndim
      do j=1,ndim
      do i=1,ndim
        V(VL(i,j),VR(k,l),:,1) = V(VL(i,j),VR(k,l),:,1) + u_tmp(i,j,k,l) ! static hubbard term
      enddo
      enddo
      enddo
      enddo


    ! broadcast from master to everyone else
      allocate(mpi_cwork(nkp))
      do l=1,ndim
      do k=1,ndim
      do j=1,ndim
      do i=1,ndim
        mpi_cwork(:) = V(VL(i,j),VR(k,l),:,1)
        call &
        mpi_bcast(mpi_cwork,nkp,mpi_double_complex,master,mpi_comm_world)
      enddo
      enddo
      enddo
      enddo
      deallocate(mpi_cwork)

    endif

    ! no frequency dependency in the ADGA case
    ! everyone for himself
    do i=2,nw
      V(:,:,:,i) = V(:,:,:,1)
    enddo
    Vend(:,:,:) = V(:,:,:,1)

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
end subroutine compute_V

!##########################################

subroutine compute_W(P,V,W)
! input/output
  complex(kind=8),intent(in) :: P(ndim**2,ndim**2,nkp,nw),V(ndim**2,ndim**2,nkp,nw)
  complex(kind=8),intent(out) :: W(ndim**2,ndim**2,nkp,nw)

! auxiliaries
  integer :: iv,ikq,i,ina,inb
  complex(kind=8), allocatable :: dmat(:,:)

 ! initialization
  W = 0.d0

  allocate(dmat(ndim**2,ndim**2))
  if(myid.eq.master) write(*,*) 'Computing W(q,iv)'


  do iv=1,nw
  do ikq=ikstart,ikend ! 1,nkp 
    dmat=0.d0
      do i=1,ndim**2 ! orbitals
         dmat(i,i) = 1.d0 ! 1 - UP
      enddo
    !multi-orbital
    dmat(:,:) = dmat(:,:) - matmul(V(:,:,ikq,iv),P(:,:,ikq,iv))
    call inverse_matrix( dmat )
    W(:,:,ikq,iv) = matmul(dmat(:,:),V(:,:,ikq,iv)) ! in this order
  enddo !ikq
  enddo !iv

  deallocate(dmat)

#ifdef MPI
! parallelization - communication
  allocate(mpi_cwork(nkp))
  do i=1,nw
    do ina=1,ndim**2
      do inb=1,ndim**2
        call MPI_ALLGATHERV(W(ina,inb,ikstart:ikend,i),ncount, MPI_DOUBLE_COMPLEX ,mpi_cwork(:),rcounts(1:nproc), displs(1:nproc),MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD ,mpierr)
        ! call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
        W(ina,inb,:,i)=mpi_cwork(:)
      enddo
    enddo
  enddo
  deallocate(mpi_cwork)
#endif

  if(myid.eq.master) then
    filen=trim(outfolder)//"/Wv.dat"
    write(*,*) 'output of Wloc(iv) to file Wv.dat'
    call local_output_Compound(W,0,filen)
  endif

  return
end subroutine compute_W

!##########################################

subroutine compute_SE(Giw,W,Vend,SE)
! input/output
  complex(kind=8), intent(in) :: Giw(ndim,ndim,nkp,2*nw)
  complex(kind=8), intent(in) :: Vend(ndim**2,ndim**2,nkp)
  complex(kind=8), intent(in) :: W(ndim**2,ndim**2,nkp,nw)
  complex(kind=8), intent(out) :: SE(ndim,ndim,nkp,2*nw)

! auxiliaries
  integer :: iw,ikp,iv,ikq,i,j,k,l,ina,inb
  complex(kind=8) :: tmpW(ndim,ndim,ndim,ndim,nkp,nw)
  complex(kind=8) :: tmpVend(ndim,ndim,ndim,ndim,nkp)
  complex(kind=8) :: FT(ndim,ndim,nkp,2*nw),ST(ndim,ndim,nkp,2*nw)
! initialization

  SE = 0.d0
  FT = 0.d0
  ST = 0.d0

  if(myid.eq.master) write(*,*) 'Computing SE(k,iw)'

! index -return ---   W_ij -> W_abcd // V_ij -> V_abcd

! this is only done once - doesnt matter which order these loops are.
  do ikq=1,nkp
  do i=1,ndim**2 
  do j=1,ndim**2
    tmpVend(inverse_l(i),inverse_r(i),inverse_r(j),inverse_l(j),ikq) = Vend(i,j,imq(ikq))
    do iv=1,nw
      tmpW(inverse_l(i),inverse_r(i),inverse_r(j),inverse_l(j),ikq,iv) = W(i,j,imq(ikq),iv)
    enddo
  enddo
  enddo
  enddo

!##########################################
!######### first term ####################
!##########################################

!! \Sigma_kl{k,iwn}
  do iw=1,nw
  do ikq=1,nkp
  do ikp=ikstart,ikend
  do k=1,ndim
  do l=1,ndim 
  do i=1,ndim
  do j=1,ndim
      FT(k,l,ikp,iw) = FT(k,l,ikp,iw) + (tmpW(k,i,j,l,ikq,1)-tmpVend(k,i,j,l,ikq)) * wtkp(ikq) &
                * Giw(i,j,ikpq(ikp,ikq),iw)  ! iv = 1
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  do iw=1,nw
  do iv=2,iw
  do ikq=1,nkp
  do ikp=ikstart,ikend
  do k=1,ndim
  do l=1,ndim 
  do i=1,ndim
  do j=1,ndim
        FT(k,l,ikp,iw) = FT(k,l,ikp,iw) + (tmpW(k,i,j,l,ikq,iv)-tmpVend(k,i,j,l,ikq)) * wtkp(ikq) &
                                        * ( Giw(i,j,ikpq(ikp,ikq),iw+iv-1) + Giw(i,j,ikpq(ikp,ikq),iw-iv+1) )
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  do iw=1,nw
  do iv=iw+1,nw
  do ikq=1,nkp
  do ikp=ikstart,ikend
  do k=1,ndim
  do l=1,ndim 
  do i=1,ndim
  do j=1,ndim
        FT(k,l,ikp,iw) = FT(k,l,ikp,iw) + (tmpW(k,i,j,l,ikq,iv)-tmpVend(k,i,j,l,ikq)) * wtkp(ikq) &
                * ( conjg(Giw(j,i,ikpq(ikp,ikq),iv-iw)) + Giw(i,j,ikpq(ikp,ikq),iw+iv-1) )
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  FT = FT/(-2.d0*beta)

!##########################################
!######### second term ####################
!##########################################

  do iw=1,nw
  do ikp=ikstart,ikend ! 1,nkp
  do ikq=1,nkp
  do l=1,ndim
  do k=1,ndim
  do i=1,ndim
  do j=1,ndim
    ST(k,l,ikp,1) = ST(k,l,ikp,1) + tmpVend(k,i,j,l,ikq) & 
        * (Giw(i,j,ikpq(ikp,ikq),iw) + conjg(Giw(j,i,ikpq(ikp,ikq),iw))) * wtkp(ikq)
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  do iw=1,2*nw
    ST(:,:,:,iw) = ST(:,:,:,1)
  enddo

  ST = ST/(2.d0*beta)   ! sum over ikq = 2

  do iw=1,2*nw
  do ikq=1,nkp
  do ikp=ikstart,ikend
  do l=1,ndim 
  do i=1,ndim
  do k=1,ndim
      ST(k,l,ikp,iw) = ST(k,l,ikp,iw) + tmpVend(k,i,i,l,ikq) * wtkp(ikq) / 4.d0
  enddo
  enddo
  enddo
  enddo
  enddo
  enddo

  SE = FT + ST

#ifdef MPI
 allocate(mpi_cwork(nkp))
  do i=1,nw
    do ina=1,ndim
      do inb=1,ndim
        call MPI_ALLGATHERV(SE(ina,inb,ikstart:ikend,i),ncount, MPI_DOUBLE_COMPLEX ,mpi_cwork(:),rcounts(1:nproc), displs(1:nproc),MPI_DOUBLE_COMPLEX , MPI_COMM_WORLD ,mpierr)
        ! call MPI_BARRIER( MPI_COMM_WORLD, mpierr )
        SE(ina,inb,:,i)=mpi_cwork(:)
      enddo
    enddo
  enddo
  deallocate(mpi_cwork)
#endif

  !this also works if compiled without MPI flag (-DMPI)
  if(myid.eq.master) then
    write(*,*) 'output of Sloc(iw) to file Sw.dat'
    filen=trim(outfolder)//"/Sw.dat"
    call local_output_Matsub(SE,1,filen) ! this routine is in io.f90
    filen=trim(outfolder)//"/Sw_diag.dat"
    call local_output_Matsub_diagonal(SE,1,filen) ! this routine is in io.f90
    filen=trim(outfolder)//"/Skw.dat"
    open(unit=10,file=filen)
      write(10,*) "## ikp, iw, iomega, RE[S_ii(ikp,iw)], IM(S_ii(ikp,iw)]"
      write(10,*) "##"
      do ikp=1,nkp
      do iw=1,min(250,nw)
        write(10,'(I7,I7,F20.8)',advance='no') ikp,iw,pi/beta*real(2*(iw-1)+1,kind=8)
        do i=1,ndim
          write(10,'(E23.8)',advance='no') real(SE(i,i,ikp,iw))
          write(10,'(E23.8)',advance='no') aimag(SE(i,i,ikp,iw))
        enddo
        write(10,*) ! line break
      enddo
      enddo
    close(10)
  endif ! master

  return
end subroutine compute_SE

!##########################################

subroutine compute_n(ncur,trace,Giw)
! input/output
  complex(kind=8), intent(in) :: Giw(ndim,ndim,nkp,2*nw)
  complex(kind=8), intent(out) :: ncur(ndim,ndim)
  complex(kind=8), intent(out) :: trace

! auxiliaries
  integer :: iw,ikp,i,j

! initialization
  ncur=0.d0

!  write(*,*) 'Computing ncur'
  do iw=1,nw
  do ikp=1,nkp
    do i=1,ndim
    do j=1,ndim
      ncur(i,j) = ncur(i,j) + (Giw(i,j,ikp,iw) + conjg(Giw(j,i,ikp,iw))) * wtkp(ikp)
    enddo
    enddo   
  enddo
  enddo

  ncur = ncur/beta  ! common prefactor, no factor of 2 for spin summation, since sum_ikp wtkp(ikp) = 2

  do i=1,ndim
    ncur(i,i)= ncur(i,i) + 1.d0
  enddo

  trace=0.d0
  do i=1,ndim
    trace=trace+ncur(i,i)
  enddo

!  write(*,*) 'number of bands', ndim
!  write(*,'(A,E15.7,E15.7)') ' number of current particles', trace

  return
end subroutine compute_n

end module computation_functions



program gw
! before running this code, you need to:
! - prepare a HMLT Hamiltonian file, for example with PRO/HMLT/3dCUBIC
! - prepare the HMLT.index.* files with PRO/HMLT/INDEX which contain tables for the indices of ik+iq and -iq

  use aux                 ! contains i, pi 
  use computation_functions   ! contains all computation functions
  use index_reference     ! contains index functions for multiorbitals  
  use hamiltonian_module  ! contains routines for reading the Hamiltonian
  use lapack_module       ! contains an easy interface for computing the inverse of a matrix
  use mpi_org     ! contains parallelization routines

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

! parameters
  open(unit=10,file="input/parameters.in", &
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
      if (flagVfile .eqv. .true.) then
        read(10,*) filename_umatrix
        read(10,*) filename_vq
      endif
    close(unit=10)

  call system("mkdir -p "//trim(outfolder))
  ! no error if existing

if (flagN .eqv. .true.) mu=0.d0 !symmetric bisection start around 0

!  nw=500 ! number of Matsubara frequencies
!  L=64    ! number of time discretization point in [0,beta) 
!  U=2.5d0 ! constant repulsion for testing purposes
!  beta=10.d0 ! inverse temperature
!  mu=0.d0 ! chemical potential
!  ntot=2.d0 ! total number of particles

! flags for enabling mu / SE cycle / V array from files
!  flagN=.true.
!  if (flagN .eqv. .true.) mu=0.d0  ! symmetric start of bisection around 0
!  flagSE=.false.
!  flagVfile=.false.

! read Hamiltonian from file "HMLT"
!  write(*,*) 'reading Hamiltonian'
  call read_hamiltonian
! this gives access to complex Hamiltonian: ndim, nkp, h(ndim,ndim,nkp)
! ndim : number of orbitals
! nkp  : number of k-points

! read index-files "HMLT.index.*"
!  write(*,*) 'reading k-point index tables'
  call read_bzindices
! this gives access to the integer tables: ikpq(nkp,nkp),imq(nkp)
! For example: given the two k-points k and q indexed by ik and iq, ikpq(ik,iq) is the index of k+q
! imq(iq) is the index for -q 

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
  SE_new = 0.d0   ! for first runthrough -> G0(ik)
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
  call compute_V(V,Vend,flagVfile)   ! reading input files to creating V matrix
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
    write(*,*) 'Time needed: ',tend-tstart
    write(*,*) 'Closing mpi environment'
  endif

  call mpi_close()

deallocate(Giw,SE_old)
end program gw
