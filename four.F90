! FOURIER transforms

subroutine invfourier(beta,L,nw,ifermi,ioff,cindata,routdata)
! ifermi=1 : fermionic case
!            G(iw_(-n))=G(iw_(n-1))*
! ifermi=0 : bosonic case
!            W(iv_(-n))=W(iv_n)*  BEWARE the irregular part of W should be substracted beforehand!!                                                                                       
  implicit none
  integer i,j,L,iwmax,nw,ifermi,ioff !offdiag element?                                                                                                                               
  double precision xpi,beta,tau,om,dummy
  DOUBLE PRECISION routdata(L)
  COMPLEX(kind=8) cindata(0:nw-1),cdummy
  iwmax=nw-1
  xpi=acos(-1.D0)
  do i=1,L
     routdata(i)=0.D0
     tau=(real(i-1,kind=8))*beta/real(L,kind=8)
     do j=0,Iwmax
        if (ifermi.eq.1) then
           om=mod((2.d0*real(j,kind=8)+real(ifermi,kind=8))*xpi/Beta*tau,2*xpi)
        else
           om=(2.d0*real(j,kind=8))*xpi/Beta*tau
        endif
        cdummy=cmplx(0.D0,1.D0)*om
        dummy=real(cindata(j)*exp(-cdummy),kind=8)
        routdata(i)=routdata(i)+2.D0/beta*dummy
     enddo
  enddo
  
  if (ifermi.eq.0) then
     routdata(:)=routdata(:)-1.d0/beta*real(cindata(0),kind=8)
  endif
  !     special treatment for tau=0:
  if (ioff.ne.1) then
     routdata(1)=-1.D0/2.D0*real(ifermi,kind=8)+routdata(1)
  endif

  RETURN
end subroutine invfourier


subroutine invfourierhp(beta,L,nw,cindata,routdata,shift) !high precision  
  !shift would be -(mu-realSigma)   !FERMI ONLY ...  
  implicit none
  integer i,j,L,iwmax,nw,ioff
  double precision xpi,beta,tau,om,dummy,dtau
  DOUBLE PRECISION routdata(L),shift
  COMPLEX(kind=8) ci,cindata(0:nw-1),cdummy,cintmp(0:nw-1)
  DOUBLE PRECISION zero,one,two

  iwmax=nw-1
  one=1.d0
  zero=0.d0
  two=2.d0
  xpi=acos(-One)
  ci=(0.d0,1.d0)

  dtau = Beta / real(L,kind=8)

  cintmp=cindata

  do i=0,iwmax
     cindata(i)=cindata(i)   -1.d0/(ci*real(2*i+1,kind=8)*xpi/beta) -shift/(ci*real(2*i+1,kind=8)*xpi/beta)**2.d0 !first, second order correction                                    
  enddo

  do i=1,L
     routdata(i)=Zero
     tau=real(i-1)*beta/real(L)
     do j=0,Iwmax
        om=mod((2*j+One)*xpi/Beta*tau,2*xpi)
        cdummy=cmplx(Zero,One)*om
        dummy=cindata(j)*exp(-cdummy)
        routdata(i)=routdata(i)+Two/beta*dummy
     enddo
  enddo
  
  do i=1,L   !!! it is -1/2 because of convention ... G(tau)<0    
     routdata(i)=routdata(i)-1.d0/2.d0 +shift*( ((+dtau*real(i-1)))/2.d0 -beta/4.d0 )
  enddo

  cindata=cintmp
  RETURN
end subroutine invfourierhp
