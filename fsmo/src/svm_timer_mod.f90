module svm_timer_mod
use,intrinsic :: iso_fortran_env
   type svm_timer_t
      integer(kind=int64) :: rate
      integer(kind=int64) :: st
      integer(kind=int64) :: fn
      real(kind=real64) :: acc
   contains
      procedure :: start => svm_timer_start
      procedure :: stop => svm_timer_stop
      procedure :: clear => svm_timer_clear
   end type
contains
   subroutine svm_timer_start(this)
     implicit none
     class(svm_timer_t) :: this
     call system_clock(count_rate=this%rate)
     call system_clock(count=this%st)
   end subroutine

   subroutine svm_timer_stop(this)
     implicit none
     class(svm_timer_t) :: this
     real(kind=real64) :: atime
     call system_clock(count=this%fn)
     atime = (this%fn - this%st)
     this%acc = this%acc + (atime/this%rate)
   end subroutine

   subroutine svm_timer_clear(this)
     implicit none
     class(svm_timer_t) :: this
     this%acc = 0.d0
   end subroutine
end module