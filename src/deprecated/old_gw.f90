program gw

  if (cyc .eq. 1) then ! configuration for one shot calculations
  ! with resource management
    allocate(Gconv(ndim,nkp,2*nw))
    call compute_Gconv(mu,Gconv)    ! Convergence Term for P with updated mu
    allocate(P(ndim**2,ndim**2,nkp,nw))
    call compute_P(mu,Giw,Gconv,P)  ! Computation of Polarization (iv,iq) with G*G
    deallocate(Gconv)
    allocate(V(ndim**2,ndim**2,nkp,nw),Vend(ndim**2,ndim**2,nkp))
    call read_V(V,Vend,flagVfile)   ! reading input files to create V matrix
    allocate(W(ndim**2,ndim**2,nkp,nw))
    call compute_W(P,V,W)     ! Computation of screened interaction
    deallocate(P,V)
    allocate(SE_new(ndim,ndim,nkp,2*nw))
    call compute_SE(Giw,W,Vend,SE_new)  ! Computation of Self-Energy
    deallocate(W,Vend)
  endif

  deallocate(Giw,SE_old)
end program gw
