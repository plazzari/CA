module GlobalMem
INTEGER, parameter        :: L=5        ! Linear dimensions
#include <define_array.h>
INTEGER, parameter        :: S=L**D ! Size of the vectorized domain
INTEGER, dimension(S)     :: v

contains

subroutine global_allocate()
#include <allocate_array.h>

end subroutine

end module
