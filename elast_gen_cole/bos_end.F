	SUBROUTINE BOS_END(recname)
	character*8 recname
	integer ierr

#include "bcs.inc"     

        call fwbos(iw,33,'0',ierr)

c       prints information about data banks 
c        call bosta
        call fparm('CLOSE recname')
        
        RETURN
        END
