      PROGRAM phonon_ext 
         character (LEN =4) typ_name(100) 
         logical lphonon
         integer ntype,ext_cell(3),typ_num(100),nat,iat
         integer itype,idime,jdime,kdime
         real,allocatable :: posi(:,:),vec(:,:) ,kin(:,:)
         real  at(3,3),phase_k(3),tmp
         real current_temp(4),dis(3)
         complex phase

         write(6,*) 'write the number of atom type'
         read(5,*) ntype

         open(13,file='POSCAR')
          read(13,*)
          read(13,*)
          read(13,*)at(1:3,1)
          read(13,*)at(1:3,2)
          read(13,*)at(1:3,3)
          read(13,*) typ_name(1:ntype)
          read(13,*) typ_num(1:ntype)
          read(13,*)
          nat = 0
          do itype = 1,ntype
            nat = nat + typ_num(itype) 
          enddo 
          allocate(posi(3,nat),vec(3,nat),kin(3,nat))
         posi = 0.0d0
          do iat = 1,nat
           read(13,*) posi(1:3,iat)
          enddo 
         close(13)


         write(6,*) 'what is extend cell?'
         read(5,*) ext_cell(1:3) 

 
        ndeme = ext_cell(1)*ext_cell(2)*ext_cell(3)

         write(6,*) 'generate phonon q?(.true. or .false.)'
         read(5,*) lphonon 

         if(lphonon) then 
         open(13,file='phonon.txt')
          do iat = 1,nat
           read(13,*) vec(1:3,iat),kin(1:3,iat)
          enddo 
         close(13)
         endif

   
         open(14,file='POSCAR_new')

          write(14,*)'asdf'
          write(14,*)'1.0000'
          write(14,*)at(1:3,1)*dble(ext_cell(1))
          write(14,*)at(1:3,2)*dble(ext_cell(2))
          write(14,*)at(1:3,3)*dble(ext_cell(3))
          write(14,*) typ_name(1:ntype)
          write(14,*) typ_num(1:ntype)*ext_cell(1)*ext_cell(2) *ext_cell(3)
          write(14,*) 'Cartesian'
          do iat = 1, nat
           do idime = 0, ext_cell(1)-1
           do jdime = 0, ext_cell(2)-1
           do kdime = 0, ext_cell(3)-1
             write(14,*)posi(1:3,iat) +at(1:3,1)*dble(idime)+at(1:3,2)*dble(jdime)+at(1:3,3)*dble(kdime)
           enddo
           enddo
           enddo
          enddo

         close(14)


       if(lphonon)then

           phase_k(1) = 6.283185307 /dble(ext_cell(1)) 
           phase_k(2) = 6.283185307 /dble(ext_cell(2)) 
           phase_k(3) = 6.283185307 /dble(ext_cell(3)) 


         open(14,file='POSCAR_dis')
          write(14,*)'asdf'
          write(14,*)'1.0000'
          write(14,*)at(1:3,1)*dble(ext_cell(1))
          write(14,*)at(1:3,2)*dble(ext_cell(2))
          write(14,*)at(1:3,3)*dble(ext_cell(3))
          write(14,*) typ_name(1:ntype)
          write(14,*) typ_num(1:ntype)*ext_cell(1)*ext_cell(2) *ext_cell(3)
          write(14,*) 'Cartesian'


           do iat = 1, nat
           do idime = 0, ext_cell(1)-1
           do jdime = 0, ext_cell(2)-1
           do kdime = 0, ext_cell(3)-1

              tmp = phase_k(1)*dble(idime)+phase_k(2)*dble(jdime)+ phase_k(3)*dble(kdime)+3.141592*0.5d0
             phase = DCMPLX(cos(tmp),sin(tmp))
               dis(1:3) = dble(phase) * vec(1:3,iat) 
             write(14,*)posi(1:3,iat) +at(1:3,1)*dble(idime)+at(1:3,2)*dble(jdime)+at(1:3,3)*dble(kdime)+dis(1:3)
           enddo
           enddo
           enddo
           enddo
             write(14,*)
             write(14,*)
           do iat = 1, nat
           do idime = 0, ext_cell(1)-1
           do jdime = 0, ext_cell(2)-1
           do kdime = 0, ext_cell(3)-1

              tmp = phase_k(1)*dble(idime)+phase_k(2)*dble(jdime)+ phase_k(3)*dble(kdime)+3.141592*0.5d0
             phase = DCMPLX(cos(tmp),sin(tmp))
               dis(1:3) = dble(phase) * vec(1:3,iat) 
             write(14,*)dis(1:3)
           enddo
           enddo
           enddo
           enddo
             write(14,*)
             write(14,*)

           do iat = 1, nat
           do idime = 0, ext_cell(1)-1
           do jdime = 0, ext_cell(2)-1
           do kdime = 0, ext_cell(3)-1

              tmp = phase_k(1)*dble(idime)+phase_k(2)*dble(jdime)+ phase_k(3)*dble(kdime)
             phase = DCMPLX(cos(tmp),sin(tmp))
               dis(1:3) = dble(phase) * kin(1:3,iat) 
             write(14,*)dis(1:3)
           enddo
           enddo
           enddo
           enddo
         close(14)

        endif







          
       
1588 format ( 4(3x, E13.6)  )
        stop 
    
         end program
