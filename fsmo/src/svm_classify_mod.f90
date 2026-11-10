module svm_classify_mod
use,intrinsic :: iso_fortran_env
use svm_data_mod
use svm_options_mod
use svm_model_mod
use svm_utils_mod
use svm_timer_mod

contains
   
   subroutine svm_validate(opts)
      implicit none
      class(svm_options_t) :: opts
      type(svm_data_t) :: data
      type(svm_model_t) :: model
      type(svm_timer_t) :: timer
      integer :: ntp,nfp,ntn,nfn,i,io
      integer,parameter :: outp = 58 
      real(kind=8) :: fx
      
      ntp = 0
      nfp = 0
      nfn = 0
      ntn = 0
      open(unit=outp,file=opts%output_file,action='write',status='new',iostat=io)
      call timer%start()
      call data%read(opts%data_file)
      call model%init(opts%model_file)
      call timer%stop()
      write(output_unit,*)'init time = ',timer%acc,' secs'
      call timer%clear()
      call timer%start()
      do i=1,data%nvecs
         fx = model%score(data%vecs(1:data%nfeat,i))
         if ( fx >= 0.0d0) then
            if ( data%y(i) >= 0.d0) then
               ntp = ntp + 1
            else
               nfp = nfp + 1
            end if
         else
            if ( data%y(i) >= 0.d0) then
               nfn = nfn + 1
            else
               ntn = ntn + 1
            end if
         end if
         write(outp)data%y(i)
         write(outp)fx
      end do
      call timer%stop()
      write(output_unit,*)'validation time = ',timer%acc,' secs'
      close(outp)
      call svm_anal(ntp,nfp,ntn,nfn)   
   end subroutine

end module

