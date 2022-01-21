      PROGRAM green_function
        implicit none
         integer nkx,nky,nk, ik,jk,nstep,istep,jstep,kstep,lstep,kk
         integer nsigma,isigma,nomega,iomega,ntime,itime
         complex,allocatable :: green0(:,:,:),green(:,:,:,:),D0(:,:)
         complex,allocatable :: self(:,:) 
         real,allocatable :: band(:),spec(:,:,:),time(:)
         real temp,elph,Vpot,mu,freq,tpi,dt,domega,ph
         real temp_eig,temp_fer ,temp_bos,prove,add_bos
         real temp_t1,temp_t2 ,temp_t3

         temp = 0.013
         elph = 0.14
         Vpot = 0.25
         mu = -0.255 
         nkx= 10
         nky= 10
         nk = nkx*nky
         nstep = 300
         freq = 0.1  
         nsigma = 2       
         prove = 5.0
         dt = 1.0
         ph = 0
         nomega =200
         domega =0.010
         ntime = 1
         tpi = 3.141592 *2.0d0
         allocate(green0(nstep,nstep,nk),green(nstep,nstep,nk,nsigma),time(ntime))
         allocate(D0(nstep,nstep),self(nstep,nstep),band(nk),spec(nk,nomega,ntime))
 
          green0 = cmplx(0.0d0,0.0d0)
          green = cmplx(0.0d0,0.0d0)
          D0 = cmplx(0.0d0,0.0d0)
          self =cmplx(0.0d0,0.0d0)  
          spec= 0.0d0
          band = 0.0d0
          time(1) = 0.5d0

!!!!!!!! evaluate the band
         open(14,file='band')
          do jk = 0,nky-1
          do ik = 0,nkx-1
            kk = jk*nkx+ik+1
             band(kk) = -2.0d0 * Vpot*cos(tpi*dble(ik)/dble(nkx)) &
                                     *cos(tpi*dble(jk)/dble(nky)) -mu
            write(14,*) kk,band(kk) 
          enddo
          enddo 
         close(14)


!!!!!!!!  evaluate the non-interacting electron  green function
         open(14,file='green0')
          do ik = 1, nk 
           do istep =1, nstep 
             do jstep =1, nstep 
               temp_eig = band(ik)*dble(istep-jstep)*dt 
               temp_fer = 1.0d0/(exp(band(ik)/temp)+1)
               if(istep.ge.jstep)then 
                green0(istep,jstep,ik) = CMPLX(0.0d0,1.0d0)*CMPLX(cos(temp_eig),-sin(temp_eig)) &
                    *(temp_fer - 1.0d0)
               else
                green0(istep,jstep,ik) = CMPLX(0.0d0,1.0d0)*CMPLX(cos(temp_eig),-sin(temp_eig)) &
                    *(temp_fer)
               endif

               write(14,*)istep,jstep,ik,green0(istep,jstep,ik)
             enddo
           enddo
          enddo
          close(14)
!!!!!! evaluate the non-interacting phonon green function
         open(14,file='phonon')
               write(14,*) 1.0d0/(exp(freq/temp)-1)
           do istep =1, nstep 
             do jstep =1, nstep 
               temp_eig = freq*dble(istep-jstep)*dt 
               temp_bos = 1.0d0/(exp(freq/temp)-1)
               add_bos = 0.5d0*(dble(istep-jstep)*dt-cos(2.0d0*freq*dble(istep)*dt) &
                                +cos(2.0d0*freq*dble(jstep)*dt))

               if(istep.ge.jstep)then 
                D0(istep,jstep) =-CMPLX(0.0d0,1.0d0)*CMPLX(cos(temp_eig),sin(temp_eig))*(temp_bos +add_bos) &
                                    -CMPLX(0.0d0,1.0d0)*CMPLX(cos(temp_eig),-sin(temp_eig))*(temp_bos+add_bos-1) 
               else
                D0(istep,jstep) =-CMPLX(0.0d0,1.0d0)*CMPLX(cos(temp_eig),sin(temp_eig))*(temp_bos + 1+add_bos) &
                                    -CMPLX(0.0d0,1.0d0)*CMPLX(cos(temp_eig),-sin(temp_eig))*(temp_bos+add_bos) 
               endif 
               write(14,*) istep,jstep, D0(istep,jstep)
             enddo
           enddo
          close(14)

!!!!!!  evaluate the self energy
         do  ik = 1 , nk
           do istep =1, nstep 
             do jstep =1, nstep 
               self(istep,jstep) = self(istep,jstep) + CMPLX(0.0d0,1.0d0) *elph**2 &
                                   *green0(istep,jstep,ik)*D0(istep,jstep) 
             enddo
           enddo
         enddo
         open(14,file='self')
           do istep =1, nstep 
             do jstep =1, nstep 
              write(14,*) istep,jstep, self(istep,jstep)
             enddo
           enddo
          close(14)
 
!!!!!!  evaluate the interacting electron green function
           do ik = 1, nk 
            do istep =1, nstep 
             do jstep =1, nstep 
              do kstep =1, nstep 
                do lstep =1, nstep 
                if((lstep.ge.jstep).and.(kstep.ge.lstep).and.(istep.ge.kstep))then
                green(istep,jstep,ik,1) = green(istep,jstep,ik,1) + &
                green0(istep,kstep,ik)*self(kstep,lstep)*green0(lstep,jstep,ik)
                endif
                enddo
              enddo
             enddo
            enddo
           enddo
           green(:,:,:,1) = green(:,:,:,1)+green0(:,:,:)
        write(6,*) 'first green'

         do isigma = 2, nsigma
           do ik = 1, nk 
            do istep =1, nstep 
             do jstep =1, nstep 
              do kstep =1, nstep 
                do lstep =1, nstep 
                if((lstep.ge.jstep).and.(kstep.ge.lstep).and.(istep.ge.kstep))then
                green(istep,jstep,ik,isigma) = green(istep,jstep,ik,isigma) +&
                green0(istep,kstep,ik)*self(kstep,lstep)*green(lstep,jstep,ik,isigma-1)
                endif
                enddo
              enddo
             enddo
            enddo
            enddo
           green(:,:,:,isigma) = green(:,:,:,isigma)+green0(:,:,:)
        write(6,*) isigma,'th green'
         enddo

!!!!!!  evaluate the tr-arpes spectra

           do itime = 1, ntime
           do iomega = 1, nomega
           do ik = 1, nk 
             do jstep =1, nstep 
               do istep =jstep, nstep 
               temp_t1 = dble(dble(iomega)-dble(nomega)/2.0d0)*domega*dble(istep-jstep)*dt 
               temp_t2 = (dble(istep)*dt-time(itime))**2/2.0d0/prove**2 
               temp_t3 = (dble(jstep)*dt-time(itime))**2/2.0d0/prove**2 
                spec(ik,iomega,itime) = spec(ik,iomega,itime) &
                + aimag(green(istep,jstep,ik,nsigma)*cmplx(cos(temp_t1),sin(temp_t1)))&
                 *exp(-temp_t2)*exp(-temp_t3)

             enddo
            enddo
           enddo
           enddo
           enddo
        write(6,*) 'done spec'

            spec = spec/tpi/prove**2


         open(14,file='spec')
           do ik = 1, nk 
            do iomega =1,nomega
              write(14,1588) ik,dble(iomega-dble(nomega)/2.0d0)*domega, spec(ik,iomega,1)
            enddo
              write(14,*) 
           enddo
         
         close(14) 
 
1588 format ( 1(I5,2x) 4(3x, E21.9)  )
        stop 
    
         end program
