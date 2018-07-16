module Mhamil
  use Mglobal
  implicit none

  integer                       :: nsham,ntet
  double precision              :: efermi
  double precision, allocatable :: wtkp(:)
  complex(dp), allocatable      :: h(:,:,:)
  double precision, allocatable :: bk(:,:)
  integer, allocatable          :: ikpq(:,:),imq(:)

contains
  
  subroutine read_hamiltonian
!   fills wtkp, h, bk
!
!   Reads hamiltonian form disk file 
!
    implicit none
    
    integer                       :: i, j, ikp, ios, ntet
    integer, allocatable          :: itt(:,:)
    double precision, parameter   :: Ry2eV = 1.36058d+1
    double precision              :: dum
    double precision, allocatable :: hr(:,:,:), hi(:,:,:)

       open( 77,file=file_hmlt,form='formatted',status='old',iostat=ios,    &
             action='read',position='rewind' )
       if( ios /= 0 )then
         write(6,*)
         write(6,*)' Cannot open file "HMLT"'
         stop
       end if
    
       read(77,*) nkp,ntet
       read(77,*) nsham,ndim

       allocate( wtkp(nkp),bk(3,nkp) )

!   Reading weights and inequivalent k-points from "HMLT"

       read(77,*) efermi
       read(77,*) (wtkp(i),i=1,nkp)                             
       read(77,*) ((bk(i,j),i=1,3),j=1,nkp)                

!   Reading inequivalent tetrahedra from file "HMLT" 

       if( ntet /= 0 )then
         allocate( itt(5,ntet) )
         read(77,*) ((itt(i,j),i=1,5),j=1,ntet)
       end if
           
!   Reading Hamiltonian from file "HMLT"
       allocate( h(ndim,ndim,nkp),hr(ndim,ndim,nkp),hi(ndim,ndim,nkp) )
   
       do ikp = 1,nkp
         read(77,*)((hr(i,j,ikp),j=i,ndim),i=1,ndim)
         read(77,*)((hi(i,j,ikp),j=i,ndim),i=1,ndim)
       end do
   
       h = cmplx( hr,hi )

       if ( ntet /= 0 )  then
          deallocate( itt ) 
          deallocate( hr,hi,bk )
       endif
    
       h = Ry2eV * h
       

       forall( i = 1:ndim ) h(i,i,:) = h(i,i,:) - efermi * Ry2eV

       dum=0.d0
       do i=1,nkp
          dum=dum+wtkp(i)
       enddo
       wtkp=wtkp*2.d0/dum


       close(77)
!!    end if
!
!   Hermitian conjugated
!
    do i = 1,ndim-1
      do j = i+1,ndim
        h(j,i,:) =  conjg(h(i,j,:))
      end do  
    end do 
    
  end subroutine read_hamiltonian
  

  subroutine read_hamiltonian_n(ch)
!
!   Reads hamiltonian form disk file 
!
    implicit none
    
    integer                       :: i, j, ikp, ios, ntet
    integer, allocatable          :: itt(:,:)
    double precision, parameter   :: Ry2eV = 1.36058d+1
    double precision              :: efermi
    double precision, allocatable :: hr(:,:,:), hi(:,:,:)
    character(len=*), intent(in)  :: ch

       open( 77,file=ch,form='formatted',status='old',iostat=ios,    &
             action='read',position='rewind' )
       if( ios /= 0 )then
         write(6,*)
         write(6,*)' Cannot open file '//ch
         stop
       end if
    
       read(77,*) nkp,ntet
       read(77,*) nsham,ndim

       allocate( wtkp(nkp),bk(3,nkp) )

!   Reading weights and inequivalent k-points from "HMLT"

       read(77,*) efermi
       read(77,*) (wtkp(i),i=1,nkp)                             
       read(77,*) ((bk(i,j),i=1,3),j=1,nkp)                

!   Reading inequivalent tetrahedra from file "HMLT" 

       if( ntet /= 0 )then
         allocate( itt(5,ntet) )
         read(77,*) ((itt(i,j),i=1,5),j=1,ntet)
       end if
           
!   Reading Hamiltonian from file "HMLT" 
       allocate( h(ndim,ndim,nkp),hr(ndim,ndim,nkp),hi(ndim,ndim,nkp) )
   
       do ikp = 1,nkp
         read(77,*)((hr(i,j,ikp),j=i,ndim),i=1,ndim)
         read(77,*)((hi(i,j,ikp),j=i,ndim),i=1,ndim)
       end do
   
       h = cmplx( hr,hi )

       if ( ntet /= 0 )  then
          deallocate( itt ) 
          deallocate( hr,hi,bk )
       endif
    
       h = Ry2eV * h
       

       forall( i = 1:ndim ) h(i,i,:) = h(i,i,:) - efermi * Ry2eV

       close(77)
!
!   Hermitian conjugated
!
    do i = 1,ndim-1
      do j = i+1,ndim
        h(j,i,:) =  conjg(h(i,j,:))
      end do  
    end do 
    
  end subroutine read_hamiltonian_n
  
 subroutine deallocate_ham
      implicit none
      deallocate(h,bk,wtkp)
end subroutine deallocate_ham



SUBROUTINE read_bzindices
      ! fills imq, ikpq
      implicit none
      integer :: ikp,jkp

      allocate(imq(nkp),ikpq(nkp,nkp))

      ikpq=0
      open(11,file=file_hmlt_kpq,status='unknown',form='unformatted')
      do ikp=1,nkp
         read(11)(ikpq(ikp,jkp),jkp=ikp,nkp)
  do jkp=1,ikp-1
    ikpq(ikp,jkp)=ikpq(jkp,ikp)   ! ikpq(j,i) = ikp(i,j)
  enddo
      enddo
      close(11)

      open(11,file=file_hmlt_mq,status='unknown',form='unformatted')
      read(11)(imq(ikp),ikp=1,nkp)
      close(11)

return
end SUBROUTINE read_bzindices


end module Mhamil
