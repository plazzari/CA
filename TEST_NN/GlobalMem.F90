module GlobalMem
INTEGER                            :: L ! Linear dimensions
#include <define_array.h>
INTEGER                            :: S ! Size of the vectorized domain
INTEGER, allocatable, dimension(:) :: v

CHARACTER (LEN=10)                 :: lattice

contains

subroutine global_allocate()

! local variable
integer   :: i,j

L=4

lattice='FCC'

#include <allocate_array.h>


SELECT CASE (TRIM(lattice))

   CASE('SC')
 
       write(*,*) 'SC'

       S=L**D

       allocate(v(S))

       FORALL (i=1:S) v(i) = i

       m=UNPACK(v,m==m,m) ! D-dimensional matrix of indexes

   CASE('FCC')


       write(*,*) 'FCC'

       n = 0

       do i=1,D
          m = 0
          do j=1,(L/2)*2,2
             m=EOSHIFT(m,SHIFT= 1,BOUNDARY=1,DIM=i)
             m=CSHIFT(m,SHIFT= 1,DIM=i)
          enddo
          n = m + n 
       enddo

       n = mod(n,2)

       S=SUM(n)

       allocate(v(S))

       FORALL (i=1:S) v(i) = i

       m = UNPACK(v,n == 1,n) ! D-dimensional matrix of indexes

       write(*,*) 'm',m

   CASE DEFAULT

       STOP     'Lattice not understood, select SC or FCC'

END SELECT

end subroutine

end module
