module smo_solver_mod
use,intrinsic :: iso_fortran_env
use svm_kernel_matrix_mod
use svm_options_mod
use svm_data_mod
use svm_utils_mod
use svm_timer_mod

   type smo_solver_t
      type(svm_timer_t) :: gap_timer, kmat_timer, fstep_timer, tstep_timer, grad_timer
      type(svm_kernel_matrix_t) :: kmat
      real(kind=8), dimension(:),allocatable :: alf
      real(kind=8), dimension(:),allocatable :: grad
      real(kind=8), dimension(:),pointer :: y
      integer(kind=4), dimension(:), allocatable :: stat
      real(kind=8) :: bias, asum, csum, gap, fun, fold, del_gap
      real(kind=8) :: eps, cost
      integer :: nvecs, nits, it, max_its
   contains
      final smo_solver_final
      procedure :: train => smo_solver_train
      procedure :: take_step => smo_solver_take_step
      procedure :: find_step => smo_solver_find_step
      procedure :: calc_gap => smo_solver_calc_gap
      procedure :: output_model_file => smo_solver_output_model_file
   end type
contains
   subroutine smo_solver_final(this)
      implicit none
      type(smo_solver_t) :: this
      deallocate(this%stat)
      deallocate(this%grad)
      deallocate(this%alf)
      this%y => null() 
   end subroutine

   integer function smo_solver_train(this, options, data) result(r)
      implicit none
      class(smo_solver_t) :: this
      class(svm_options_t),intent(in) :: options
      class(svm_data_t),intent(in),target :: data
      type(svm_timer_t) :: train_timer
      integer :: it,i, jt
      real(kind=8) :: fun,fold 

      call this%kmat%init(options, data)
      this%nvecs = data%nvecs 
      allocate (this%alf(this%nvecs))
      allocate (this%grad(this%nvecs))
      allocate (this%stat(this%nvecs))
      this%max_its = options%max_its
      if ( this%max_its == 0) then
         this%max_its = 10 
      end if      
      this%y => data%y(1:this%nvecs)
      call this%fstep_timer%clear()
      call this%gap_timer%clear()
      call this%tstep_timer%clear()
      call this%kmat_timer%clear() 
      call this%grad_timer%clear()
      call train_timer%clear()
      do i = 1, this%nvecs
         this%alf(i) = 0.d0
         this%grad(i) = data%y(i)
         this%stat(i) = -1
      end do
      this%cost = options%cost
      this%eps = options%eps
      fold = 0.0d0
      r = 0
      
      write(output_unit,*)'SMO Solver'
      write(output_unit,*)' max iteration = ',this%max_its
      write(output_unit,*)' # of vectors  = ',this%nvecs
      write(output_unit,*)' # of cache    = ',this%kmat%cache_size
      write(output_unit,*)' cost          = ',this%cost
      write(output_unit,*)' eps           = ',this%eps
      write(output_unit,*)' '
      call train_timer%start()
      do it = 1, this%max_its
         do jt=1,this%nvecs
            r = this%find_step()
            if ( r /= 0) then
               write(error_unit,*)' converged on jt = ',jt 
               exit	
            end if
         end do
         fun = this%calc_gap()
         write (output_unit, *) 'iteration = ', it
         write (output_unit, *) 'obj fun   = ', fun
         write (output_unit, *) 'deltaF    = ', dabs(fun - fold)
         write (output_unit, *) 'gap       = ', this%gap
         write (output_unit, *) 'deltaGap  = ', this%del_gap
         write (output_unit, *) 'bias      = ', this%bias
         fold = fun
         if (r /= 0) then
            exit
         end if
         if (dabs(this%del_gap) < this%eps) then
            r = 1
            exit
         end if
      end do
      call train_timer%stop()
      write (output_unit, *) "train     time = ", train_timer%acc
      write (output_unit, *) "grad      time = ", this%grad_timer%acc
      write (output_unit, *) "grad      time = ", this%kmat_timer%acc
      write (output_unit, *) "find_step time = ", this%fstep_timer%acc
      write (output_unit, *) "take_step time = ", this%tstep_timer%acc
      write (output_unit, *) "gap       time = ", this%gap_timer%acc
      if (r == 0) then
         write (output_unit, *) ' solver did not converge!'
         stop -1
      else
         write (output_unit, *) ' solver converged in ',it,' iterations'
      end if
      r = 0
      call this%output_model_file(options)
   end function

   subroutine smo_solver_output_model_file(this, options)
      implicit none
      class(smo_solver_t) :: this
      class(svm_options_t),intent(in) :: options
      integer :: nbsv,nsv,i,j,ntp,nfp,ntn,nfn, io
      real(kind=8) :: dtmp,fx
      integer,parameter :: outp = 48
      integer,parameter :: outp2 = 50
      integer,parameter :: inp = 49      
      nbsv = 0
      nsv = 0
      do i = 1, this%nvecs
         if (this%stat(i) >= 0) then
            nsv = nsv + 1
            if (this%stat(i) > 0) then
               nbsv = nbsv + 1
            end if
         end if
      end do
      open (unit=outp, file=options%model_file, action='write', access='STREAM')
      write (outp) nsv
      write (outp) this%kmat%nfeat
      write (outp) options%kernel_type
      write (outp) options%kernel_power
      write (outp) options%kernel_c1
      write (outp) options%kernel_c2
      write (outp) this%bias
      write (outp) options%kernel_scale
      do i = 1, this%nvecs
         if (this%stat(i) >= 0) then
            dtmp = this%y(i)*this%alf(i)*this%kmat%scale(i)
            write (outp) dtmp
         end if
      end do
      do i = 1, this%nvecs
         if (this%stat(i) >= 0) then
            do j = 1, this%kmat%nfeat
               write (outp) this%kmat%vecs(j, i)
            end do
         end if
      end do
      close (outp)
      open (unit=outp2, file=options%output_file, action='write', status='new', iostat=io, access='STREAM')
      do i = 1, this%nvecs
         fx = (this%grad(i) + this%bias)*this%y(i)
         if (fx > 0) then
            if (this%y(i) >= 0.0d0) then
               ntp = ntp + 1
            else
               nfp = nfp + 1
            end if
         else
            if (this%y(i) >= 0.0d0) then
               nfn = nfn + 1
            else
               ntn = ntn + 1
            end if
         end if
         write (outp2) this%y(i)
         write (outp2) fx
      end do
      close (outp2)
      call svm_anal(ntp, nfp, ntn, nfn)
   end subroutine

   integer function smo_solver_find_step(this) result(r)
      implicit none
      class(smo_solver_t) :: this
      integer :: imax, imin, i
      real(kind=real64) :: gmax, gmin, gx, ys

      imax = -1
      imin = -1
      gmax = -1.e30
      gmin =  1.e30
      call this%fstep_timer%start()
      do i = 1, this%nvecs
         ys = this%y(i)*this%stat(i)
         gx = this%grad(i)
         if (ys > 0.0d0) then
            if (gx < gmin) then
               gmin = gx
               imin = i
            end if
         else
            if (gx > gmax) then
               gmax = gx
               imax = i
            end if
         end if
      end do
      call this%fstep_timer%stop()
            write(error_unit,*)'imax =',imax,' imin =',imin
            write(error_unit,*)'gmax =',gmax,' gmin =',gmin
      if ( (imax /= -1) .and. (imin /= -1) .and. (imax /= imin)) then
         if (dabs(gmax - gmin) > this%eps) then
            write(error_unit,*)' take step'
            r = this%take_step(imax, imin)
         else
            write(error_unit,*)' gradients are equal'
            r = -1
         end if
      else
         write(error_unit,*)' no feasible step found'
         r = -1
      end if
      write(error_unit,*)' return ',r
   end function

   integer function smo_solver_take_step(this, imax, imin) result(r)
      implicit none
      class(smo_solver_t) :: this
      integer, intent(in) :: imax, imin
      real(kind=8) :: a1, a2, ai, aj, qc, ds, L, H, gam, da1, da2
      real(kind=8), dimension(:), pointer :: qmax
      real(kind=8), dimension(:), pointer :: qmin
      real(kind=8), parameter :: tau = 2.d-14
      integer :: i, st1, st2

      call this%kmat_timer%start()
      call this%kmat%get_row(imax,qmax)
      call this%kmat%get_row(imin,qmin)
      call this%kmat_timer%stop()
      call this%tstep_timer%start()
      ds = this%y(imax)*this%y(imin)
      a1 = this%alf(imax)
      a2 = this%alf(imin)
      ai = a1
      aj = a2
      gam = a2 + ds*a1
      if (ds >= 0.0) then
         L = dmax1(0.0d0, (gam - this%cost))
         H = dmin1(this%cost, gam)
      else
         L = dmax1(0.0d0, gam)
         H = dmin1(this%cost, (gam + this%cost))
      end if
      qc = dmax1( (qmax(imax) + qmin(imin) - qmax(imin) - qmax(imin)), tau)
      a2 = a2 + (((this%grad(imax) - this%grad(imin))*this%y(imin))/qc)
      a2 = dmin1(dmax1(a2, L), H)
      a1 = a1 + ds*(aj - a2)
      if (a1 > tau) then
         if (a1 < (this%cost - tau)) then
            st1 = 0
         else
            st1 = 1
            a1 = this%cost
            a2 = gam - ds*this%cost
         end if
      else
         a1 = 0.0d0
         a2 = gam
         st1 = -1
      end if
      if (a2 > tau) then
         if (a2 < (this%cost - tau)) then
            st2 = 0
         else
            st2 = 1
         end if
      else
         st2 = -1
      end if
      call this%tstep_timer%stop()
      r = 0
      if (dabs(a2 - aj) > tau) then
         this%stat(imax) = st1
         this%stat(imin) = st2
         this%alf(imax) = a1
         this%alf(imin) = a2
         call this%grad_timer%start()
         da1 = this%y(imax)*(a1 - ai)
         da2 = this%y(imin)*(a2 - aj)
         do i = 1, this%nvecs
            this%grad(i) = this%grad(i) - (da1*qmax(i) + da2*qmin(i))
         end do
         call this%grad_timer%stop()
      else
         r = -1
      end if
   end function

   real(kind=8) function smo_solver_calc_gap(this) result(fsum)
      class(smo_solver_t) :: this
      real(kind=8) :: asum, csum, bsum
      integer(kind=8) :: nsum, i

      call this%gap_timer%start()
      asum = 0.0d0
      fsum = 0.0d0
      bsum = 0.0d0
      nsum = 0
      csum = 0.0d0
      do i = 1, this%nvecs
         asum = asum + this%alf(i)
         fsum = fsum + this%alf(i)*this%grad(i)*this%y(i)
      end do
      fsum = (fsum + asum)*0.5d0
      do i = 1, this%nvecs
         if (this%stat(i) == 0) then
            bsum = bsum + this%grad(i)
            nsum = nsum + 1
         end if
      end do
      bsum = -bsum
      if (nsum /= 0) then
         bsum = bsum/nsum
      end if
      do i = 1, this%nvecs
         csum = csum + dmax1(0.0d0, this%y(i)*(this%grad(i) + bsum))
      end do
      csum = csum*this%cost
      gap = (csum + asum - fsum - fsum)/(1.0 + asum + csum - fsum)
      this%del_gap = gap - this%gap
      this%gap = gap
      this%fun = fsum
      this%bias = bsum
      call this%gap_timer%stop()
   end function
end module
