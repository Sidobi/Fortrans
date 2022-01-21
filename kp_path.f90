      PROGRAM average
          real :: one
         real,allocatable :: kpath(:,:) 
         integer nk, nkpath,ik


          one = 1.00

         open(13,file='kp.txt') 
          read(13,*) nk,nkpath 
          allocate( kpath(3,nk))
           do ik = 1, nk
            read(13,*) kpath(1:3,ik)     
           enddo
         close(13)

 
         open(13,file='kp_path.txt')
          write(13,*) (nk-1)*(nkpath+1)-(nk-2)
         do ik = 1,nk-1
         do ipath = 0,nkpath
           if((ik.ne.1).and.(ipath==0)) go to 40
           write(13,1588) kpath(1:3,ik)+(kpath(1:3,ik+1)-kpath(1:3,ik))*dble(ipath)/dble(nkpath),one
40    continue
         enddo
         enddo

         close(13)


1588 format ( 4(3x, f17.8)  )
        stop

         end program
