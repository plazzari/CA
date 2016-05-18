PROGRAM test_pack_1
USE NearNeighbour

implicit none

!Local variables
INTEGER                   :: i,j

call global_allocate()

  call compute_von_Neumann()
! call compute_Moore()
! call compute_Honeycomb()

   DO i=1,S
    write(*,*) 'point=',i, 'NNidx',vNN(i,:)
   ENDDO


END PROGRAM
