       program project_band
      character (len = 40) :: char1
      character (len = 2) :: cob
      character (len = 3) :: temp3
      real, allocatable :: band(:,:,:),proj(:,:,:)
      real :: fermi 
      integer :: nbnd,nkp,nat,ibnd,ik,iorb,iat
      integer :: tep1,temp2,temp4,temp5
      nbnd = 780 
      nkp = 121 
      nat = 155
      fermi =  10.4295 
      allocate(band(nbnd,nkp,7),proj(nbnd,nkp,5))

       proj = 0.0d0
       band = 0.0d0
      open(24,file='bands.dat.gnu')
        do ibnd= 1, nbnd
          do ik = 1, nkp
            read(24,*)  band(ibnd,ik,1:2)
          enddo
        enddo
      close(24)
        do ibnd= 1, nbnd
          do ik = 1, nkp
            band(ibnd,ik,2) = band(ibnd,ik,2) -fermi
          enddo
        enddo


       open(30,file='bands.projwfc_up.new' )
       do iat = 1,nat
!   1    1  Ta1 5S     1    0    1
        read(30,*) tmp1,temp2, temp3 , cob, temp4,temp5, iorb

        if(cob.ne.'5D') then
         do ibnd =1,(nbnd*nkp)
           read(30,*)
         enddo
        else
            write(6,*)temp3,cob,iorb

        do ik = 1, nkp
         do ibnd= 1, nbnd
            read(30,*) temp1, temp2, proj(ibnd,ik,iorb)
          enddo
         enddo
         band(1:nbnd,1:nkp,iorb+2) =band(1:nbnd,1:nkp,iorb+2)+proj(1:nbnd,1:nkp,iorb)
       endif
       enddo

       close(30)

       open(46,file='proj.band.txt')
         do ibnd= 1, nbnd
          do ik = 1, nkp
           write(46,119) band(ibnd,ik,1:7)
          enddo 
         enddo 
      close(46)
 119  format(9(F10.5))
!

      stop
      end
      
