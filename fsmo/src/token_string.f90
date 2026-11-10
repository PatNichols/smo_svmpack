   integer function tokenize_string(str, delims, tokens) result(nt)
      implicit none
      character(len=*), intent(in) :: str
      character(len=*), intent(in) :: delims
      character(len=:), dimension(:), allocatable :: tokens
      integer :: slen, pos, new, maxlen

      if (allocated(tokens)) then
         deallocate (tokens)
      end if
      maxlen = 0
      slen = len_trim(str)
      pos = 1
      nt = 0
      do
         if ( pos > slen) then
           exit
         end if
         if ( str(pos:pos) == ' ') then
            pos = pos + 1
            cycle
         end if
         new = scan(str(pos:slen), delims)
         if (new == 1) then
            pos = pos + 1
            cycle
         end if
         nt = nt + 1
         if (new == 0) then
            maxlen = max(maxlen,(slen-pos+1))
            exit
         end if
         maxlen = max(maxlen, (new-pos+1))
         pos = pos + new
         if (pos > slen) exit
      end do
      allocate (character(len=maxlen) :: tokens(nt))
      pos = 1
      nt = 0
      do
         if ( pos > slen) then
            return
         end if
         if ( str(pos:pos) == ' ') then
            pos = pos + 1
            cycle
         end if
         new = scan(str(pos:slen), delims)
         if (new == 1) then
            pos = pos + 1
            cycle
         end if
         nt = nt + 1
         if (new == 0) then
            tokens(nt) = trim(str(pos:slen))
            write(error_unit,*)'final token = ',tokens(nt),' from ',trim(str(pos:slen))
            write(error_unit,*)'# of tokens = ',nt,' pos = ',pos,' slen = ',slen
            return
         else
            tokens(nt) = trim(str(pos:(pos + new - 1)))
            write(error_unit,*)'next token = ',tokens(nt),' from ',trim(str(pos:(pos+new-1)))
            write(error_unit,*)'# of tokens = ',nt,' pos = ',pos,' slen = ',slen,' new ',new
            pos = pos + new
         end if
         if (pos > slen) exit
      end do
   end function tokenize_string   
