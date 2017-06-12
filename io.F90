
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Some subroutines for input/output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



SUBROUTINE local_output_Matsub(G,ifermi,ch)
  use aux
  use hamiltonian_module, only : nkp, ndim, wtkp
  IMPLICIT NONE
  CHARACTER(LEN=80) ch
  integer ina,inb
  INTEGER i,ifermi  ! ifermi=1 ... fermionic ... ifermi=0 ... bosonic 
  COMPLEX(kind=8)G(ndim,ndim,nkp,2*nw)
  COMPLEX(kind=8)GL(ndim,ndim,2*nw)

  GL=0.d0
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo
  GL=GL/2.d0 !SPIN                                                                                                                                                                    
  open(unit=10,file=trim(ch),status='unknown')
  write(10,*) 'i, iomega, RE / IM [Variable(i,j)] in order: 11, 21, ... , 12, 22, ...'
  if (ifermi.eq.1) then
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i,real(2*(i-1)+1,kind=8)*pi/beta,((GL(ina,inb,i),inb=1,ndim),ina=1,ndim)
     enddo
  else
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1),kind=8)*pi/beta,((GL(ina,inb,i),inb=1,ndim),ina=1,ndim)
     enddo
  endif
  close(10)
  RETURN
END SUBROUTINE local_output_Matsub

SUBROUTINE local_output_Matsub_diagonal(G,ifermi,ch)
  use aux
  use hamiltonian_module, only : nkp, ndim, wtkp
  IMPLICIT NONE
  CHARACTER(LEN=80) ch
  integer ina
  INTEGER i,ifermi
  COMPLEX(kind=8)G(ndim,ndim,nkp,2*nw)
  COMPLEX(kind=8)GL(ndim,ndim,2*nw)

  GL=0.d0
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo 
GL=GL/2.d0! SPIN                                                                                                                                                                    
  open(unit=10,file=trim(ch),status='unknown')
  write(10,*) 'i, iomega, RE / IM [Variable(i,j)] in order: 11, 22, ... '
  if (ifermi.eq.1) then
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1)+1,kind=8)*pi/beta,(GL(ina,ina,i),ina=1,ndim)
     enddo
  else
     do i=1,min(250,nw)
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1),kind=8)*pi/beta,(GL(ina,ina,i),ina=1,ndim)
     enddo
  endif
  close(10)
  RETURN
END SUBROUTINE local_output_Matsub_diagonal

SUBROUTINE local_output_Compound(G,ifermi,ch)
  use aux
  use hamiltonian_module, only : nkp, ndim, wtkp
  IMPLICIT NONE
  CHARACTER(LEN=80)ch
  integer ina,inb
  INTEGER i,ifermi
  COMPLEX(kind=8)G(ndim**2,ndim**2,nkp,nw)
  COMPLEX(kind=8)GL(ndim**2,ndim**2,nw)

  GL=0.d0
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo
  GL=GL/2.d0 !SPIN                                                                                                                                                                    
  open(unit=10,file=trim(ch),status='unknown')
  write(10,*) 'i, iomega, RE / IM [Variable(i,j)] in order: 11, 21, ... , 12, 22, ...'
  if (ifermi.eq.1) then
     do i=1,nw
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1)+1,kind=8)*pi/beta, &
        ((GL(ina,inb,i),inb=1,ndim**2),ina=1,ndim**2)
     enddo
  else
     do i=1,nw
        write(10,'(I7,F20.8,1000E23.8)') i, real(2*(i-1),kind=8)*pi/beta, &
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
  IMPLICIT NONE
  character(len=80) :: ch
  integer :: i,j,skip
  double precision :: f
  
  double precision :: tmp(ndim,nw)

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
  IMPLICIT NONE
  character(len=80) :: ch
  integer :: i,j,skip
  double precision :: tmp(ndim)

  open(16,file=(ch),form='formatted',status='old',    &
             action='read',position='rewind' )
  
  do i=1,skip
    read(16,*)
  enddo
  
  read(16,*) (tmp(j),j=1,ndim)

  close(16)
  
END SUBROUTINE



SUBROUTINE local_output_tau(G,nlm,L,beta,ifermi,jump,ch)
  use hamiltonian_module
  IMPLICIT NONE
  CHARACTER(LEN=80)ch
  integer nlm,ina,inb
  INTEGER L,i,ifermi
  
  DOUBLE PRECISION beta,jump(nlm,nlm)
  DOUBLE PRECISION G(nlm,nlm,nkp,L),GL(nlm,nlm,L)
  
  GL=0.d0
  
  do i=1,nkp
     GL(:,:,:)=GL(:,:,:)+wtkp(i)*G(:,:,i,:)
  enddo !i                                                                                                                                                                            
  GL=GL/2.d0  !SPIN                                                                                                                                                                   
  
  open(unit=10,file=(ch),status='unknown')
  do i=1,L
     write(10,'(3000F12.6)')real(i-1,kind=8)*beta/real(L,kind=8),((GL(ina,inb,i),inb=1,nlm),ina=1,nlm)
  enddo
  if (ifermi.eq.1) then
     write(10,'(3000F12.6)')beta,(((-jump(ina,inb)-GL(ina,inb,1)),inb=1,nlm),ina=1,nlm)
  else
     write(10,'(3000F12.6)')beta,((GL(ina,inb,1),inb=1,nlm),ina=1,nlm)
  endif
  close(10)
  
  RETURN
END SUBROUTINE local_output_tau

