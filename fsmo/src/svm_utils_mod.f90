module svm_utils_mod
   use, intrinsic :: iso_fortran_env
contains
   subroutine svm_anal(ntp,nfp,ntn,nfn)
      implicit none
      integer,intent(in) :: nfp,ntp,ntn,nfn
      integer :: err,ntot
      real(kind=real64) :: acc,tpr,tnr,ppv,npv,fscore

      err = nfp + nfn
      ntot = real(err+ntp+ntn,kind=real64)
      acc = 1.0 - (real(err,kind=real64)/real(ntot,kind=real64))
      tpr = real(ntp,8)/real(ntp+nfn,kind=real64)
      tnr = real(ntn,8)/real(ntn+nfp,kind=real64)
      ppv = real(ntp)/real(ntp+nfp,kind=real64)
      npv = real(ntn)/real(ntn+nfn,kind=real64)
      fscore = real(ntp+ntp,kind=real64)/real(ntp+ntp+nfp+nfn,kind=real64)
      write(output_unit,*)'# of predictions = ',ntot
      write(output_unit,*)'# of errors      = ',err
      write(output_unit,*)'accuracy         = ',acc
      write(output_unit,*)'sensitivity      = ',tpr
      write(output_unit,*)'specificity      = ',tnr
      write(output_unit,*)'precision        = ',ppv
      write(output_unit,*)'balanced acc     = ',(tpr+tnr)*0.5d0
      write(output_unit,*)'f1 score         = ',fscore
   end subroutine

   subroutine readline(inp,line,ier)
   implicit none
   integer,intent(in) :: inp
   character(len=:),allocatable,intent(out) :: line
   integer,intent(out)                      :: ier
  
   integer,parameter                     :: buflen=1024
   character(len=buflen)                 :: buffer
   integer                               :: last
   integer                               :: isize

   line=''
   ier=0

   ! read characters from line and append to result
   INFINITE: do
      ! read next buffer (an improvement might be to use stream I/O
      ! for files other than stdin so system line limit is not
      ! limiting)
      read(inp,iostat=ier,fmt='(a)',advance='no',size=isize) buffer
      ! append what was read to result
      if(isize.gt.0)line=line//buffer(:isize)
      ! if hit EOR reading is complete unless backslash ends the line
      if(is_iostat_eor(ier))then
         last=len(line)
         if(last.ne.0)then
            ! if line ends in backslash it is assumed a continued line
            ! if(line(last:last).eq.'\')then
            !   No escape character in Fortran, so only one backslash here.
            !   Ref: https://fortranwiki.org/fortran/files/character_handling_in_Fortran.html
            if(line(last:last).eq.'\')then
               ! remove backslash
               line=line(:last-1)
               ! continue on and read next line and append to result
               cycle INFINITE
            endif
         endif
         ! hitting end of record is not an error for this routine
         ier=0
         ! end of reading line
         exit INFINITE
     ! end of file or error
     elseif(ier.ne.0)then
        exit INFINITE
     endif
   enddo INFINITE

   line=trim(line)
   end subroutine readline

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
         new = scan(str(pos:slen), delims)
         if (new == 1) then
            pos = pos + 1
            cycle
         end if
         nt = nt + 1
         if (new == 0) then
            tokens(nt) = trim(str(pos:slen))
!            write(output_unit,*)'final token = ',trim(tokens(nt)),' from ',trim(str(pos:slen))
!            write(output_unit,*)'# of tokens = ',nt,' pos = ',pos,' slen = ',slen
            return
         else
            tokens(nt) = trim(str(pos:(pos + new - 2)))
!            write(output_unit,*)nt,' token = ',trim(tokens(nt)),' from ',trim(str(pos:(pos+new-2)))
!            write(output_unit,*)'# of tokens = ',nt,' pos = ',pos,' slen = ',slen
            pos = pos + new
         end if
         if (pos > slen) exit
      end do
   end function tokenize_string
end module
