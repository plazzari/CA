module GlobalMem

type domain
INTEGER                               :: L        ! Linear dimensions
#include <define_array.h>
INTEGER                               :: S        ! Size of the vectorized domain
CHARACTER (LEN=10)                    :: lattice
LOGICAL                               :: shiftby1 ! Shift by 1 the lattice FCC 
INTEGER                               :: NN       ! number of near nieghbours
INTEGER, allocatable, dimension(:,:)  :: v1DtoND  ! map of 1D to ND
INTEGER, allocatable, dimension(:,:)  :: vNN      ! matrix of vectorized nearneighbours
end type domain

type(domain)                          :: CA_dom1
REAL   , parameter                    :: PI=acos(-1.)
REAL(4), allocatable, dimension(:,:)  :: examap
LOGICAL, allocatable, dimension(:)    :: v_aux

contains

subroutine allocate_dom(CA_dom)

type(domain)                       :: CA_dom
! local variable
integer   :: i,j
integer   :: L,D
INTEGER, allocatable, dimension(:)    :: v        ! state_vector


D = CA_dom%D
L = CA_dom%L

#include <allocate_array.h>


SELECT CASE (TRIM(CA_dom%lattice))

   CASE('SC','HC')
 
       write(*,*) 'SC'

       CA_dom%S=CA_dom%L**CA_dom%D

       allocate(v(CA_dom%S))

       FORALL (i=1:CA_dom%S) v(i) = i

       CA_dom%m=UNPACK(v,CA_dom%m==CA_dom%m,CA_dom%m)! D-dimensional matrix of indexes

   CASE('FCC')


       write(*,*) 'FCC'

       CA_dom%n = 0
       write(*,*) shape(CA_dom%n)

       do i=1,D
          CA_dom%m = 0
          do j=1,(CA_dom%L/2)*2,2
             CA_dom%m=EOSHIFT(CA_dom%m,SHIFT= 1,BOUNDARY=1,DIM=i)
             CA_dom%m=CSHIFT(CA_dom%m,SHIFT= 1,DIM=i)
          enddo
          CA_dom%n = CA_dom%m + CA_dom%n 
       enddo

       CA_dom%n = mod(CA_dom%n,2)

       CA_dom%S=SUM(CA_dom%n)

       allocate(v(CA_dom%S))

       FORALL (i=1:CA_dom%S) v(i) = i

       CA_dom%m = UNPACK(v,CA_dom%n == 1,CA_dom%n) ! D-dimensional matrix of indexes

       if (CA_dom%shiftby1)  CA_dom%m = -CA_dom%m + 1 ! Zero cells became 1 and viceversa


   CASE DEFAULT

       STOP     'Lattice not understood, select SC or FCC'

END SELECT

! construct map between 1D world and ND world

allocate(CA_dom%v1DtoND(CA_dom%S,D))

do i=1,D
   CA_dom%n = 0
     do j=1,CA_dom%L
        CA_dom%n=EOSHIFT(CA_dom%n,SHIFT= 1,BOUNDARY=j,DIM=i)
     enddo
   CA_dom%v1DtoND(:,i) =PACK(CA_dom%n,CA_dom%m>0)
enddo


deallocate(v)

end subroutine

subroutine dump_lattice(CA_dom)
implicit none
!local variables
type(domain)   :: CA_dom
integer        :: i,j
character(len=1024)  :: filename
character(len=1024)  :: cube_side

write (cube_side, "(I3.3)") CA_dom%L
filename =  TRIM(CA_dom%lattice)//TRIM(cube_side)//'.xyz'

OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
    write(unit=333,FMT=*) CA_dom%S
    write(unit=333,FMT=*)


SELECT CASE (CA_dom%lattice)

     CASE('HC')
         allocate(examap(CA_dom%S,2))
         allocate(v_aux(CA_dom%S))
         examap = 0.
         v_aux  = .FALSE.
         call create_exa(CA_dom,1,0.,0.)
         do i=1,CA_dom%S
             write(unit=333,FMT=*) 'Au', examap(i,1),examap(i,2), CA_dom%v1DtoND(i,3:CA_dom%D)
         enddo
     CASE DEFAULT
         do i=1,CA_dom%S
             write(unit=333,FMT=*) 'Au', CA_dom%v1DtoND(i,:)
         enddo

END SELECT

CLOSE(UNIT=333)

end subroutine

subroutine deallocate_dom(CA_dom)
implicit none
type(domain)   :: CA_dom

deallocate(CA_dom%m)
deallocate(CA_dom%n)
deallocate(CA_dom%v1DtoND)
deallocate(CA_dom%vNN)

end subroutine

recursive subroutine create_exa(CA_dom,point0,x0,y0)
implicit none
type(domain)   :: CA_dom
integer    :: i,point0,point
real(4)    :: x0,y0
real(4)    :: x,y

if (point0 .eq. 0) then
   return
elseif (v_aux(point0)) then
   return
else
   v_aux(point0)    = .TRUE.
   examap(point0,1) = x0
   examap(point0,2) = y0
   do i =1,6
        point = CA_dom%vNN(point0,i)
        x     = x0 + cos(-2*PI/6.*real(i-1) + 5./6. * PI )
        y     = y0 + sin(-2*PI/6.*real(i-1) + 5./6. * PI )
       call create_exa(CA_dom,point,x,y)
   enddo
   do i =7,CA_dom%NN
        point = CA_dom%vNN(point0,i)
        x     = x0 
        y     = y0 
       call create_exa(CA_dom,point,x,y)
   enddo
endif

end subroutine
end module
