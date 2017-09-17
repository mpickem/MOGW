
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Some subroutines for input/output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE local_output_Matsub(G,ifermi,ch)
  use aux
  use hamiltonian_module, only : nkp, ndim, wtkp
  IMPLICIT NONE
  character(len=80) :: ch
  integer           :: ina,inb
  integer           :: i,ifermi  ! ifermi=1 ... fermionic ... ifermi=0 ... bosonic
  complex(DP)       :: G(ndim,ndim,nkp,2*nw)
  complex(DP)       :: GL(ndim,ndim,2*nw)

  GL=0.d0
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo
  GL=GL/2.d0 !SPIN                                                                                                                                                                    
  open(unit=10,file=trim(ch),status='unknown')
  write(10,*) 'i, iomega, RE / IM [Variable(i,j)] in order: 11, 12, ... , 21, 22, ...'
  if (ifermi.eq.1) then
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i,real(2*(i-1)+1,dp)*pi/beta,((GL(ina,inb,i),inb=1,ndim),ina=1,ndim)
     enddo
  else
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1),dp)*pi/beta,((GL(ina,inb,i),inb=1,ndim),ina=1,ndim)
     enddo
  endif
  close(10)
  RETURN
END SUBROUTINE local_output_Matsub

SUBROUTINE local_output_Matsub_diagonal(G,ifermi,ch)
  use aux
  use hamiltonian_module, only : nkp, ndim, wtkp
  IMPLICIT NONE
  character(len=80) :: ch
  integer           :: ina
  integer           :: i,ifermi
  complex(DP)       :: G(ndim,ndim,nkp,2*nw)
  complex(DP)       :: GL(ndim,ndim,2*nw)

  GL=0.d0
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo 
GL=GL/2.d0! SPIN                                                                                                                                                                    
  open(unit=10,file=trim(ch),status='unknown')
  write(10,*) 'i, iomega, RE / IM [Variable(i,j)] in order: 11, 22, ... '
  if (ifermi.eq.1) then
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1)+1,dp)*pi/beta,(GL(ina,ina,i),ina=1,ndim)
     enddo
  else
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1),dp)*pi/beta,(GL(ina,ina,i),ina=1,ndim)
     enddo
  endif
  close(10)
  RETURN
END SUBROUTINE local_output_Matsub_diagonal

SUBROUTINE local_output_Compound(G,ifermi,ch)
  use aux
  use hamiltonian_module, only : nkp, ndim, wtkp
  IMPLICIT NONE
  character(len=80) :: ch
  integer           :: ina,inb
  integer           :: i,ifermi
  complex(DP)       :: G(ndim**2,ndim**2,nkp,nw)
  complex(DP)       :: GL(ndim**2,ndim**2,nw)

  GL=0.d0
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo
  GL=GL/2.d0 !SPIN                                                                                                                                                                    
  open(unit=10,file=trim(ch),status='unknown')
  write(10,*) 'i, iomega, RE / IM [Variable(i,j)] in order: 11, 12, ... , 21, 22, ...'
  if (ifermi.eq.1) then
     do i=1,nw
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1)+1,dp)*pi/beta, &
        ((GL(ina,inb,i),inb=1,ndim**2),ina=1,ndim**2)
     enddo
  else
     do i=1,nw
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1),dp)*pi/beta, &
        ((GL(ina,inb,i),inb=1,ndim**2),ina=1,ndim**2)
!            ((((GL(ina,inb,inc,ind,i),ind=1,ndim),inc=1,ndim), inb=1,ndim),ina=1,ndim)
     enddo
  endif
  close(10)
  RETURN
END SUBROUTINE local_output_Compound


SUBROUTINE input_V(tmp,ch,skip)
  use aux
  use hamiltonian_module, only: ndim
  implicit none
  character(len=80) :: ch
  integer           :: i,j,skip
  real(dp)          :: f
  real(dp)          :: tmp(ndim,nw)

  open(15,file=(ch),form='formatted',status='old',    &
             action='read',position='rewind' )
  
  do i=1,skip
    read(15,*)
  enddo
  
  do i=1,nw
    read(15,*) f,(tmp(j,i),j=1,ndim)
  enddo

  close(15)
  
END SUBROUTINE

SUBROUTINE input_Vend(tmp,ch,skip)
  use aux
  use hamiltonian_module, only: ndim
  implicit none
  character(len=80) :: ch
  integer           :: i,j,skip
  real(dp)          :: tmp(ndim)

  open(16,file=(ch),form='formatted',status='old',    &
             action='read',position='rewind' )
  
  do i=1,skip
    read(16,*)
  enddo
  
  read(16,*) (tmp(j),j=1,ndim)

  close(16)
  
END SUBROUTINE



! SUBROUTINE local_output_tau(G,nlm,L,beta,ifermi,jump,ch)
!   use aux
!   use hamiltonian_module
!   implicit none
!   character(len=80) :: ch
!   integer           :: nlm,ina,inb
!   integer           :: L,i,ifermi
!   real(dp)          :: beta,jump(nlm,nlm)
!   real(dp)          :: G(nlm,nlm,nkp,L),GL(nlm,nlm,L)
  
!   GL=0.d0
  
!   do i=1,nkp
!      GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
!   enddo !i                                                                                                                                                                            
!   GL=GL/2.d0  !SPIN                                                                                                                                                                   
  
!   open(unit=10,file=(ch),status='unknown')
!   do i=1,L
!      write(10,'(3000F12.6)')real(i-1,dp)*beta/real(L,kind=8),((GL(ina,inb,i),inb=1,nlm),ina=1,nlm)
!   enddo
!   if (ifermi.eq.1) then
!      write(10,'(3000F12.6)')beta,(((-jump(ina,inb)-GL(ina,inb,1)),inb=1,nlm),ina=1,nlm)
!   else
!      write(10,'(3000F12.6)')beta,((GL(ina,inb,1),inb=1,nlm),ina=1,nlm)
!   endif
!   close(10)
  
!   RETURN
! END SUBROUTINE local_output_tau
