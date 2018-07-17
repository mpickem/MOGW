module Mcomp
  use Mglobal
  use Mindex
  use Mhamil
  use Mlapack

contains

subroutine compute_Giw(mu_loc)
! input/output
  real(dp), intent(in)       :: mu_loc

! auxiliaries
  integer                    :: ikp,iw,i,j,ina,inb
  complex(dp), allocatable   :: cmat(:,:)

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
         cmat(i,i) = ci * real(2*(iw-1)+1,kind=8)*pi/beta + mu_loc
      enddo
      ! iw + mu - H(k) - SE :
      cmat(:,:)=cmat(:,:)-h(:,:,ikp)-SE(:,:,ikp,iw)
      ! [iw + mu - H(k) - Sigma(w,k) ]^{-1} :
      call inverse_matrix( cmat ) ! inverts complex square matrix
      Giw(:,:,ikp,iw)=cmat
    enddo ! ikp
  enddo !iw
  deallocate(cmat)

  return
end subroutine compute_Giw

subroutine compute_Gconv
! auxiliaries
  integer                  :: ikp,iw,i,ina
  complex(dp), allocatable :: cmat(:)

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
      Gconv(:,ikp,iw) = cmat(:)
    enddo ! ikp
  enddo !iw
  deallocate(cmat)

  return
end subroutine compute_Gconv

!########################################

subroutine compute_P
!auxiliaries
  integer                  :: iv,ikq,iw,ikp,i,j,k,l,ina,inb
  complex(dp)              :: G1,G2,ctmp
  integer                  :: IL(ndim,ndim),IR(ndim,ndim)

!initialization
  P = 0.d0
  G1 = 0.d0
  G2 = 0.d0
  ctmp = 0.d0

!compound indizes
  call index_(IL,IR)

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
    ! write(*,*) 'P special case needed'
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
  enddo
  enddo
  enddo
  enddo
  enddo

  return
end subroutine compute_P

subroutine compute_W
! auxiliaries
  integer                  :: iv,ikq,i,ina,inb
  complex(dp), allocatable :: dmat(:,:)

 ! initialization
  W = 0.d0

  allocate(dmat(ndim**2,ndim**2))

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

  return
end subroutine compute_W

!##########################################

subroutine compute_SE
! auxiliaries
  integer                  :: iw,ikp,iv,ikq,i,j,k,l,ina,inb
  complex(dp)              :: tmpW(ndim,ndim,ndim,ndim,nkp,nw)
  complex(dp)              :: tmpVend(ndim,ndim,ndim,ndim,nkp)
  complex(dp)              :: FT(ndim,ndim,nkp,2*nw),ST(ndim,ndim,nkp,2*nw)
! initialization

  SE = 0.d0
  FT = 0.d0
  ST = 0.d0

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

  return
end subroutine compute_SE

!##########################################

subroutine compute_n(ncur, trace, ndim_loc)
! input/output
  integer, intent(in) :: ndim_loc
  complex(dp), intent(inout) :: ncur(ndim_loc,ndim_loc)
  complex(dp), intent(inout) :: trace

! auxiliaries
  integer                      :: iw,ikp,i,j
! initialization
  ncur=0.d0

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

  return
end subroutine compute_n

end module Mcomp
