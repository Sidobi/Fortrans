      PROGRAM average
         real  current(3)
         real current_temp(4)
         integer ix,iy,iz,nx,ny,nz,totn

         write(6,*)' numkp- x,y,z '
         read(5,*) nx,ny,nz
         totn = nx*ny*nz

         open(13,file='kp.txt')
         write(13,*)totn
         do ix = 1,nx
         do iy = 1,ny
         do iz = 1,nz
          write(13,1588) 1.0d0*dble(ix-1)/dble(nx), 1.0d0*dble(iy-1)/dble(ny) , 1.0d0*dble(iz-1)/dble(nz) ,1.d0/dble(totn)
         enddo
         enddo
         enddo

         close(13)


1588 format ( 4(3x, f17.8)  )
        stop

         end program
