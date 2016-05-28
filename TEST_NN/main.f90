PROGRAM test_pack_1
USE GlobalMem
USE NearNeighbour

implicit none

!Local variables
INTEGER                   :: i

!TEST FCC--------------------

CA_dom1%L=10
CA_dom1%lattice='FCC'

call allocate_dom(CA_dom1)

call compute_FCC(CA_dom1)

  DO i=1,CA_dom1%S
    write(*,*) 'point=',i, 'NNidx',CA_dom1%vNN(i,:)
  ENDDO

call  dump_lattice(CA_dom1)

call deallocate_dom(CA_dom1)

!TEST vN--------------------
write(*,*) 'von Neumann'

CA_dom1%L=13
CA_dom1%lattice='SC'

call allocate_dom(CA_dom1)

call compute_von_Neumann(CA_dom1)

  DO i=1,CA_dom1%S
    write(*,*) 'point=',i, 'NNidx',CA_dom1%vNN(i,:)
  ENDDO

call  dump_lattice(CA_dom1)

call deallocate_dom(CA_dom1)

!TEST Moore--------------------
write(*,*) 'Moore'

CA_dom1%L=40
CA_dom1%lattice='SC'

call allocate_dom(CA_dom1)

call compute_Moore(CA_dom1)

  DO i=1,CA_dom1%S
    write(*,*) 'point=',i, 'NNidx',CA_dom1%vNN(i,:)
  ENDDO

call  dump_lattice(CA_dom1)

call deallocate_dom(CA_dom1)

!TEST Honeycomb--------------------
write(*,*) 'Honeycomb'

CA_dom1%L=7
CA_dom1%lattice='SC'

call allocate_dom(CA_dom1)

call compute_Honeycomb(CA_dom1)

  DO i=1,CA_dom1%S
    write(*,*) 'point=',i, 'NNidx',CA_dom1%vNN(i,:)
  ENDDO

call  dump_lattice(CA_dom1)

call deallocate_dom(CA_dom1)


END PROGRAM
