module NearNeighbour
USE GlobalMem

implicit none

contains
!--------von Neumann----------!
!     ________________________________
!     |     |     |     |     |     |
!     |     |     |     |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |     |  1  |     |     |      -1,1 --> Near Neighbour
!     |_____|_____|_____|_____|_____|__       O --> Site considered
!     |     |     |     |     |     |
!     |     | -1  |  O  |  1  |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |     | -1  |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |     |     |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
subroutine compute_von_Neumann(CA_dom)
   implicit none
type(domain)                       :: CA_dom
!Local variables
   INTEGER                   :: i,j
   INTEGER                   :: NNvN

   CA_dom%NN=2*CA_dom%D      ! von Neumann


   allocate(CA_dom%vNN(CA_dom%S,CA_dom%NN))
    
   j=1
   DO i=1,CA_dom%D
     CA_dom%n=EOSHIFT(CA_dom%m,SHIFT= 1,BOUNDARY=0,DIM=i) ;CA_dom%vNN(:,j) =PACK(CA_dom%n,CA_dom%n==CA_dom%n) ; j=j+1
     CA_dom%n=EOSHIFT(CA_dom%m,SHIFT=-1,BOUNDARY=0,DIM=i) ;CA_dom%vNN(:,j) =PACK(CA_dom%n,CA_dom%n==CA_dom%n) ; j=j+1
   ENDDO

end subroutine
!---------Moore---------------!
!     ________________________________
!     |     |     |     |     |     |       (-1,-1)  --> 1 Near Neighbour
!     |     |     |     |     |     |       (-1, 0)  --> 2 Near Neighbour
!     |_____|_____|_____|_____|_____|__     (-1, 1)  --> 3 Near Neighbour
!     |     |     |     |     |     |       ( 0,-1)  --> 4 Near Neighbour
!     |     |  6  |  7  |  8  |     |       ( 0, 0)  --> O Site considered
!     |_____|_____|_____|_____|_____|__     ( 0, 1)  --> 5 Near Neighbour
!     |     |     |     |     |     |       ( 1,-1)  --> 6 Near Neighbour
!     |     |  4  |  O  |  5  |     |       ( 1, 0)  --> 7 Near Neighbour
!     |_____|_____|_____|_____|_____|__     ( 1, 1)  --> 8 Near Neighbour
!     |     |     |     |     |     |
!     |     |  1  |  2  |  3  |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |     |     |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
subroutine compute_Moore(CA_dom)
   implicit none
type(domain)                       :: CA_dom
!Local variables
   INTEGER                          :: i,j,a,res,NNMo
   INTEGER,DIMENSION(CA_dom%D)      :: move                      

   CA_dom%NN=3**CA_dom%D-1   ! Moore 

   allocate(CA_dom%vNN(CA_dom%S,CA_dom%NN))

   a = 0
   DO i=1,CA_dom%NN
     res=CodeBase(REAL(i-1+a,4),3,move,CA_dom%D)
     if (res .EQ. 0) then
       a = 1
       res=CodeBase(REAL(i-1+a,4),3,move,CA_dom%D)
     endif
     write(*,*) 'move-->', move
     CA_dom%n=CA_dom%m
     Do j=1,CA_dom%D
       CA_dom%n=EOSHIFT(CA_dom%n,SHIFT= move(j),BOUNDARY=0,DIM=j)
     ENDDO
     CA_dom%vNN(:,i) =PACK(CA_dom%n,CA_dom%n==CA_dom%n)
   ENDDO

end subroutine
!---------Honeycomb----------!
!     ________________________________
!     |     |     |     |     |     |
!     |     |     |     |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |     |  2  |  3  |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |       1     --> Site considered
!     |     |  7  |  1  |  4  |     |       [2,7] --> Near Neighbour
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |  6  |  5  |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
!     |     |     |     |     |     |
!     |_____|_____|_____|_____|_____|__
!     |     |     |     |     |     |
! Scheme from Libbrecht K. G. Physically Derived Rules for Simulating
! Faceted Crystal Growth using Cellular Automata pg.16
! http://arxiv.org/pdf/0807.2616v1.pdf
subroutine compute_Honeycomb(CA_dom)
   implicit none
type(domain)                       :: CA_dom
!Local variables
   INTEGER                   :: i,j,a,res
   INTEGER,DIMENSION(6,2)    :: moveHC

   if (CA_dom%D .LT. 2) STOP 'Honeycomb requires dimension > 1'

   CA_dom%NN=6 + 2*(CA_dom%D-2)   ! Honeycomb

   moveHC(1,1)= 0; moveHC(1,2)= 1 ! Position --> 2
   moveHC(2,1)= 1; moveHC(2,2)= 1 ! Position --> 3
   moveHC(3,1)= 1; moveHC(3,2)= 0 ! Position --> 4
   moveHC(4,1)= 0; moveHC(4,2)=-1 ! Position --> 5
   moveHC(5,1)=-1; moveHC(5,2)=-1 ! Position --> 6
   moveHC(6,1)=-1; moveHC(6,2)= 0 ! Position --> 7

   allocate(CA_dom%vNN(CA_dom%S,CA_dom%NN))
      
   FORALL (i=1:CA_dom%S) CA_dom%v(i) = i

   CA_dom%m=UNPACK(CA_dom%v,CA_dom%m==CA_dom%m,CA_dom%m) ! D-dimensional matrix of indexes
   CA_dom%n=CA_dom%m
   DO i=1,6
     CA_dom%n=EOSHIFT(CA_dom%n,SHIFT= moveHC(i,1),BOUNDARY=0,DIM=1)
     CA_dom%n=EOSHIFT(CA_dom%n,SHIFT= moveHC(i,2),BOUNDARY=0,DIM=2)
     CA_dom%vNN(:,i) =PACK(CA_dom%n,CA_dom%n==CA_dom%n)
     CA_dom%n=CA_dom%m
   ENDDO

   j=7
   DO i=3,CA_dom%D
     CA_dom%n=EOSHIFT(CA_dom%m,SHIFT= 1,BOUNDARY=0,DIM=i) ;CA_dom%vNN(:,j) =PACK(CA_dom%n,CA_dom%n==CA_dom%n) ; j=j+1
     CA_dom%n=EOSHIFT(CA_dom%m,SHIFT=-1,BOUNDARY=0,DIM=i) ;CA_dom%vNN(:,j) =PACK(CA_dom%n,CA_dom%n==CA_dom%n) ; j=j+1
   ENDDO

end subroutine
!---------FCC-----------------!
subroutine compute_FCC(CA_dom)
   implicit none
   type(domain)                       :: CA_dom
!Local variables
   INTEGER                             :: i,j,a,res
   INTEGER,allocatable, DIMENSION(:,:) :: moveFCC

   CA_dom%NN = CA_dom%D*(2*(CA_dom%D-1)) ! FCC Nearneighbours

   allocate(moveFCC(CA_dom%NN,CA_dom%D))
   allocate(CA_dom%vNN(CA_dom%S,CA_dom%NN))

   call CodeBaseFCC(moveFCC,CA_dom%D,CA_dom%NN)

   DO i=1, CA_dom%NN
     CA_dom%n=CA_dom%m
     Do j=1,CA_dom%D
       CA_dom%n=EOSHIFT(CA_dom%n,SHIFT= moveFCC(i,j),BOUNDARY=0,DIM=j)
     ENDDO
     CA_dom%vNN(:,i) =PACK(CA_dom%n,CA_dom%m>0)
   ENDDO

   deallocate(moveFCC)
end subroutine
!-----------------------------!
subroutine CodeBaseFCC(moveFCC,D,NNFCC)
implicit none
INTEGER                    :: D,NNFCC
INTEGER                    :: i,j,res
INTEGER,DIMENSION(D)       :: move
INTEGER,DIMENSION(NNFCC,D) :: moveFCC
   j=0
   DO i=0,3**D-1

      res=CodeBase(REAL(i,4),3,move,D)

!     write(*,*) 'moveA', i, move
      if (CheckMoveFCC(move,D)) then
         j=j+1
         moveFCC(j,:)=move
!        write(*,*) 'moveB',i,  move
      endif

   ENDDO

if (j .NE. NNFCC) STOP 'Problem in computing FCC near neighbours'

return
end subroutine

logical Function CheckMoveFCC(move,D)
implicit none
INTEGER                   :: i,a,b,D
INTEGER,DIMENSION(D)      :: move
a = 0
b = 0
CheckMoveFCC = .FALSE.

do i=1,D
     if ( move(i) .EQ. -1) a = a +1
     if ( move(i) .EQ.  1) b = b +1
enddo

if (( a .EQ. 1) .AND. ( b .EQ. 1)  ) CheckMoveFCC = .TRUE.
if (( a .EQ. 2) .AND. ( b .EQ. 0)  ) CheckMoveFCC = .TRUE.
if (( a .EQ. 0) .AND. ( b .EQ. 2)  ) CheckMoveFCC = .TRUE.

return
end

!Convert a number from base 10 to base b. The function returns
!FALSE if b not in [2..36] or if string x contains invalid
!characters in base 10 or if number x is too big
integer Function CodeBase(x,b,out_vect,D)
implicit none
real*4 x
integer b, n, D
integer, dimension(D) :: out_vect
character(30):: y
  CodeBase=0
  if (b<2.or.b>36) then
    print *,' Base must be between 2 and 36 !'
    return
  end if
  y=''
  do while (x>0.)
    n=INT(x-b*INT(x/b))
    if (n<10) then 
          y=ACHAR(IACHAR('0')+n)//y
    else 
          y=ACHAR(IACHAR('A')+n-10)//y
    end if
    x=INT(x/b)           
  end do
  do n=1,D-LEN(TRIM(y))
          y=ACHAR(IACHAR('0'))//y
  enddo
  do n=1,LEN(TRIM(y))
          read(y(n:n),'(i10)') out_vect(n)
  enddo
  out_vect = out_vect -1 
  CodeBase=1
  if (DOT_PRODUCT(out_vect,out_vect) .EQ. 0) then
    CodeBase=0
  endif
  return
end
end module
