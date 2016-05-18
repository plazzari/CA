module NearNeighbour
USE GlobalMem

implicit none

INTEGER :: NNvN,NNMo,NNHC 
INTEGER, allocatable, dimension(:,:)  :: vNN   ! matrix of vectorized near neigh...


contains
!--------Von Neumann----------!
subroutine compute_Von_neumann()
!Local variables

   implicit none

   INTEGER                   :: i,j

   NNvN=2*D      ! von Neumann


   allocate(vNN(S,NNvN))
    
   FORALL (i=1:S) v(i) = i

   m=UNPACK(v,m==m,m) ! D-dimensional matrix of indexes

! Start computing V.N. Near Neigh D-dimensional case
   j=1
   DO i=1,D
     n=EOSHIFT(m,SHIFT= 1,BOUNDARY=0,DIM=i) ;vNN(:,j) =PACK(n,n==n) ; j=j+1
     n=EOSHIFT(m,SHIFT=-1,BOUNDARY=0,DIM=i) ;vNN(:,j) =PACK(n,n==n) ; j=j+1
   ENDDO
! End 

end subroutine
!---------Moore---------------!
subroutine compute_Moore()

   implicit none
!Local variables
   INTEGER                   :: i,j,a,res,NNMo
   INTEGER,DIMENSION(D)      :: move                      

   NNMo=3**D-1   ! Moore 

   allocate(vNN(S,NNMo))
      
   FORALL (i=1:S) v(i) = i

   m=UNPACK(v,m==m,m) ! D-dimensional matrix of indexes
   
   a = 0
   DO i=1,NNMo
     res=CodeBase(REAL(i-1+a,4),3,move)
     if (res .EQ. 0) then
       a = 1
       res=CodeBase(REAL(i-1+a,4),3,move)
     endif
     write(*,*) 'move-->', move
     n=m
     Do j=1,D
       n=EOSHIFT(n,SHIFT= move(j),BOUNDARY=0,DIM=j)
     ENDDO
     vNN(:,i) =PACK(n,n==n)
   ENDDO

end subroutine
!---------Honey_Comb----------!
subroutine compute_Honey_Comb()

   implicit none
!Local variables
   INTEGER                   :: i,j,a,res
   INTEGER,DIMENSION(6,2)    :: moveHC

   if (D .LT. 2) then
       STOP
   endif

! Scheme from Libbrecht K. G. Physically Derived Rules for Simulating
! Faceted Crystal Growth using Cellular Automata pg.16 
! http://arxiv.org/pdf/0807.2616v1.pdf

   NNHC=6 + 2*(D-2)   ! Honey_Comb
   moveHC(1,1)= 0; moveHC(1,2)= 1 ! Position --> 2
   moveHC(2,1)= 1; moveHC(2,2)= 1 ! Position --> 3
   moveHC(3,1)= 1; moveHC(3,2)= 0 ! Position --> 4
   moveHC(4,1)= 0; moveHC(4,2)=-1 ! Position --> 5
   moveHC(5,1)=-1; moveHC(5,2)=-1 ! Position --> 6
   moveHC(6,1)=-1; moveHC(6,2)= 0 ! Position --> 7

   allocate(vNN(S,NNHC))
      
   FORALL (i=1:S) v(i) = i

   m=UNPACK(v,m==m,m) ! D-dimensional matrix of indexes
   n=m
   DO i=1,6
     n=EOSHIFT(n,SHIFT= moveHC(i,1),BOUNDARY=0,DIM=1)
     n=EOSHIFT(n,SHIFT= moveHC(i,2),BOUNDARY=0,DIM=2)
     vNN(:,i) =PACK(n,n==n)
     n=m
   ENDDO

   j=7
   DO i=3,D
     n=EOSHIFT(m,SHIFT= 1,BOUNDARY=0,DIM=i) ;vNN(:,j) =PACK(n,n==n) ; j=j+1
     n=EOSHIFT(m,SHIFT=-1,BOUNDARY=0,DIM=i) ;vNN(:,j) =PACK(n,n==n) ; j=j+1
   ENDDO

end subroutine

!-----------------------------!

!Convert a number from base 10 to base b. The function returns
!FALSE if b not in [2..36] or if string x contains invalid
!characters in base 10 or if number x is too big
integer Function CodeBase(x,b,out_vect)
implicit none
real*4 x
integer b, n
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
