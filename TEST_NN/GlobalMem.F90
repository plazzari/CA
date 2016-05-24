module GlobalMem

type domain
INTEGER                               :: L        ! Linear dimensions
#include <define_array.h>
INTEGER                               :: S        ! Size of the vectorized domain
INTEGER, allocatable, dimension(:)    :: v
CHARACTER (LEN=10)                    :: lattice
INTEGER                               :: NN       ! number of near nieghbours
INTEGER, allocatable, dimension(:,:)  :: vNN      ! matrix of vectorized nearneighbours
end type domain

type(domain)                       :: CA_dom1, prova(10,10)

contains

subroutine allocate_dom(CA_dom)
implicit none
type(domain)                       :: CA_dom
! local variable
integer   :: i,j
integer   :: L,S,D



D = CA_dom%D


L = CA_dom%L
#include <allocate_array.h>


SELECT CASE (TRIM(CA_dom%lattice))

   CASE('SC')
 
       write(*,*) 'SC'

       CA_dom%S=CA_dom%L**CA_dom%D

       allocate(CA_dom%v(CA_dom%S))

       FORALL (i=1:CA_dom%S) CA_dom%v(i) = i

       CA_dom%m=UNPACK(CA_dom%v,CA_dom%m==CA_dom%m,CA_dom%m)! D-dimensional matrix of indexes

   CASE('FCC')


       write(*,*) 'FCC'

       CA_dom%n = 0
       write(*,*) shape(CA_dom%n)

       do i=1,D
          CA_dom%m = 0
          do j=1,floor(REAL(CA_dom%L)/2)*2,2
             CA_dom%m=EOSHIFT(CA_dom%m,SHIFT= 1,BOUNDARY=1,DIM=i)
             CA_dom%m=CSHIFT(CA_dom%m,SHIFT= 1,DIM=i)
          enddo
          CA_dom%n = CA_dom%m + CA_dom%n 
       enddo

       CA_dom%n = mod(CA_dom%n,2)

       CA_dom%S=SUM(CA_dom%n)

       allocate(CA_dom%v(CA_dom%S))

       FORALL (i=1:CA_dom%S) CA_dom%v(i) = i

       CA_dom%m = UNPACK(CA_dom%v,CA_dom%n == 1,CA_dom%n) ! D-dimensional matrix of indexes

       write(*,*) 'm',CA_dom%m

   CASE DEFAULT

       STOP     'Lattice not understood, select SC or FCC'

END SELECT

end subroutine
subroutine deallocate_dom(CA_dom)
implicit none
type(domain)   :: CA_dom
deallocate(CA_dom%m)
deallocate(CA_dom%n)
deallocate(CA_dom%v)
deallocate(CA_dom%vNN)

end subroutine


end module
