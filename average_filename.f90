      PROGRAM average 
         character (LEN =20) filename 
         character (LEN =20) outname 
         real  current(3)
         real current_temp(4)
         integer ntime,nband,i,step,istep
         current(1:3) =0.0d0 

         write(6,*) 'write the filename'
         read(5,*) filename


       outname = 'avg_' // TRIM( filename ) 
 

         open(13,file=filename)

        step = 0
        do
          read(13,*,IOSTAT=io)
          if(io == 0 ) then
            step = step + 1
          elseif(io < 0) then
            exit
          endif
        enddo
       rewind(13)

      current = 0.0d0
   
         open(14,file=outname)

         do  istep=1,step        

          read(13,*) current_temp(1:4) 
          current(1:3) = current(1:3)*real(istep-1)/real(istep) + current_temp(2:4)/real(istep) 
          write(14,1588) current_temp(1) , current(1:3)

         enddo

         close(13)
         close(14)
          
       
1588 format ( 4(3x, E13.6)  )
        stop 
    
         end program
