module svm_options_mod
use,intrinsic :: iso_fortran_env
use program_options_mod

   type svm_options_t
      real(kind=8) :: eps, cost
      real(kind=8) :: kernel_c1, kernel_c2
      integer :: kernel_type, kernel_scale, kernel_power, kernel_cache_rows
      integer :: nthreads, task, max_its
      character(len=64) :: output_file, model_file, data_file
   contains
      procedure :: init => svm_options_init
      procedure :: write => svm_options_write
   end type
contains
   subroutine svm_options_write(this,outp)
      class(svm_options_t),intent(in) :: this
      integer,intent(in) :: outp
      write(outp,*)'# of threads = ',this%nthreads
      write(outp,*)'data         = ',this%data_file
      write(outp,*)'model        = ',this%model_file
      write(outp,*)'out          = ',this%output_file
      write(outp,*)'task         = ',this%task
      if ( this%task == 0) then 
         write(outp,*)'kernel parameters'
         write(outp,*)'type = ',this%kernel_type
         write(outp,*)'power= ',this%kernel_power
         write(outp,*)'scale= ',this%kernel_scale
         write(outp,*)'c1   = ',this%kernel_c1
         write(outp,*)'c2   = ',this%kernel_c2
         write(outp,*)'cache= ',this%kernel_cache_rows
         write(outp,*)'solver parameters'
         write(outp,*)'eps     = ',this%eps
         write(outp,*)'cost    = ',this%cost
         write(outp,*)'max its = ',this%max_its
         write(outp,*)'eps = ',this%eps
      else
         write(outp,*)'validation task'
      end if    
   end subroutine
   
   subroutine svm_options_init(this)
      class(svm_options_t),intent(out) :: this
      type(prog_opts_t) :: popts
      character(len=254) :: vstr,task_str,fstr
      
      call popts%init()
      call popts%add_options_from_file('svm_opts.in')
      call popts%write(output_unit)
      call popts%parse_command_line()
      if (popts%has_value('config') .eqv. .true.) then
         call popts%get_value('config', fstr)
         call popts%parse_file(fstr)
      end if
      call popts%write(output_unit)
      call popts%get_value('task', task_str)
      call popts%get_value('model_file', this%model_file)
      call popts%get_value('data_file', this%data_file)
      call popts%get_value('output_file', this%output_file)
      call popts%get_value('nthreads', vstr)
      read (vstr, *) this%nthreads
      if (trim(task_str) == 'classify') then
         this%task = 1
      else
         this%task = 0
         call popts%get_value('max_its', vstr)
         read (vstr, *) this%max_its
         call popts%get_value('kernel_type', vstr)
         read (vstr, *) this%kernel_type
         call popts%get_value('kernel_power', vstr)
         read (vstr, *) this%kernel_power
         call popts%get_value('kernel_scale', vstr)
         read (vstr, *) this%kernel_scale
         call popts%get_value('kernel_c1', vstr)
         read (vstr, *) this%kernel_c1
         call popts%get_value('kernel_c2', vstr)
         read (vstr, *) this%kernel_c2
         call popts%get_value('kernel_cache_rows', vstr)
         read (vstr, *) this%kernel_cache_rows
         call popts%get_value('eps', vstr)
         read (vstr, *) this%eps
         call popts%get_value('cost', vstr)
         read (vstr, *) this%cost
         call this%write(output_unit)
      end if
   end subroutine
end module
