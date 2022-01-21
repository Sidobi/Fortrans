      real, allocatable :: eigup(:,:),eigdn(:,:),feigup(:,:)
      real, allocatable :: feigdn(:,:),xrk(:,:)
      real :: fermi,spin,dk,temp
      open(13,file='EIGENVAL')
      open(24,file='band.txt')
      open(25,file='band_up.txt')
      open(26,file='band_dw.txt')
 118  format(900('      ', F10.5))
 119  format(16x,3(F10.5))
!
!   please remember the number of writting field in the format should be larger than the number of bands
!
      write(6,*) "spin unpolar(type '0') or spin polar(type '1')"
      read(5,*) spin
      write(6,*) "fermi level"
      read(5,*) fermi

       If(spin.eq.0) then
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)Nelectron,nk,nband
      allocate (feigup(nk,nband),eigup(nk,nband))
      allocate (xrk(nk,3))

      do i =1,nk
         read(13,*)
!        read(13,*)
         read(13,*)xrk(i,1),xrk(i,2),xrk(i,3),temp
        do j = 1,nband 
            read(13,*)itmp,eigup(i,j) 
         enddo
      enddo
      
      feigup(:,:) = eigup(:,:) - fermi

      do k = 1,nband
        dk =0.0d0
         do i = 1,nk
          write(24,118)dk,feigup(i,k)
        if(i.ne.1) dk = dk + sqrt((xrk(i,1)-xrk(i-1,1))**2+(xrk(i,2)-xrk(i-1,2))**2+(xrk(i,3)-xrk(i-1,3))**2)
         enddo
         write(24,*)
      enddo
       else

      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)
      read(13,*)Nelectron,nk,nband
      allocate (feigup(nk,nband),feigdn(nk,nband),eigup(nk,nband))
      allocate (eigdn(nk,nband),xrk(nk,3))

      do i =1,nk
         read(13,*)
         read(13,*)xrk(i,1),xrk(i,2),xrk(i,3),temp
!        read(13,*)
         do j = 1,nband 
            read(13,*)itmp,eigup(i,j),eigdn(i,j)
         enddo
      enddo
      
      feigup(:,:) = eigup(:,:) - fermi
      feigdn(:,:) = eigdn(:,:) - fermi

        do k = 1,nband
         dk = 0.d0
         do i = 1,nk
          write(24,118)dk,feigup(i,k),feigdn(i,k)
          write(25,118)dk,feigup(i,k)
          write(26,118)dk,feigdn(i,k)
        if(i.ne.nk) dk = dk + sqrt((xrk(i+1,1)-xrk(i,1))**2+(xrk(i+1,2)-xrk(i,2))**2+(xrk(i+1,3)-xrk(i,3))**2)
         enddo
         write(24,*)
         write(25,*)
         write(26,*)
      enddo

      endif

      stop
      end
      
