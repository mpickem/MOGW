module computation_functions
  use aux
  use index_reference
  use hamiltonian_module
  use lapack_module
  use mpi_org
  implicit none
  character(len=80) :: filen
  private :: filen
  public :: compute_Giw, compute_Gconv, compute_P, compute_W, compute_n

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
      write(10,*) "## ikp, kx, ky, kz, iw, iomega, RE[S_ij(ikp,iw)], IM(S_ij(ikp,iw)]"
      write(10,*) "##"
      do ikp=1,nkp
      do iw=1,min(250,nw)
        write(10,'(I7,3F20.8,I7,F20.8)',advance='no') ikp, bk(1,ikp), bk(2,ikp), bk(3,ikp) ,iw,pi/beta*real(2*(iw-1)+1,kind=8)
        write(10,'(20E23.8)',advance='no') ((real(SE(i,j,ikp,iw)), aimag(SE(i,j,ikp,iw)), j=1,ndim),i=1,ndim)
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
