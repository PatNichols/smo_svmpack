module program_options_mod
   use, intrinsic :: iso_fortran_env
   use :: svm_utils_mod
   integer,parameter :: LINE_SZ = 80
   
   type opt_t
      character(len=LINE_SZ) :: key
      character(len=LINE_SZ) :: des
      character(len=LINE_SZ) :: val
      integer :: req, stat
   contains
      procedure :: init => opt_init
      procedure :: set_value => opt_set_value
      procedure :: write => opt_write
      procedure :: valid => opt_valid
      procedure :: matchs => opt_matchs 
      procedure,private,pass(this) :: opt_copy
      generic,public :: assignment(=) => opt_copy 
   end type

   type prog_opts_t
      type(opt_t), dimension(:), pointer :: opts => NULL()
      integer :: allow_extra_options, allow_extra_args
      integer :: ncap, nopts
   contains
      final prog_opts_final
      procedure :: init => prog_opts_init
      procedure :: parse_command_line => prog_opts_parse_command_line
      procedure :: parse_file => prog_opts_parse_file
      procedure :: add_option => prog_opts_add_option
      procedure :: add_options_from_file => prog_opts_add_options_from_file
      procedure :: get_value => prog_opts_get_value
      procedure :: set_value => prog_opts_set_value
      procedure :: has_value => prog_opts_has_value
      procedure :: has_option => prog_opts_has_option
      procedure :: write => prog_opts_write
   end type

contains
   subroutine opt_init(this,key,des,val,req)
      implicit none
      class(opt_t),intent(out) :: this
      character(len=*),intent(in) :: key
      character(len=*),intent(in) :: des
      character(len=*),intent(in) :: val
      integer,optional,intent(in) :: req
      integer :: sz
      
!      if (allocated(this%val)) deallocate(this%val)
!      if (allocated(this%des)) deallocate(this%des)
!      if (allocated(this%key)) deallocate(this%key)
      sz = len_trim(key) 
      if ( sz /= 0) then
!         allocate(character(len=len_trim(key)) :: this%key)
         this%key = trim(key)
      else
         write(error_unit,*)'all option must have a valid keyword'
         stop -1
      end if
      
      sz = len_trim(des) 
      if ( sz /= 0) then
!         allocate(character(len=len_trim(des)) :: this%des)
         this%des = trim(des)
      else
         this%des = repeat(' ',LINE_SZ)
      end if
      
      sz = len_trim(val) 
      if ( sz /= 0) then
!         allocate(character(len=len_trim(val)) :: this%val)
	 this%val = trim(val)
         this%stat = 0
      else
         this%val = repeat(' ',LINE_SZ) 
         this%stat = -1
      end if
      if ( present(req) ) then
          this%req = req
      else
          this%req = 0
      end if
   end subroutine
   
   subroutine opt_copy(this,o2)
      implicit none
      class(opt_t),intent(out) :: this
      class(opt_t),intent(in) :: o2
      call this%init(o2%key,o2%des,o2%val,o2%req)      
   end subroutine

   pure logical function opt_matchs(this,key) result(r)
      implicit none
      class(opt_t),intent(in) :: this
      character(len=*),intent(in) :: key
      r = ( trim(this%key) == trim(key) )
   end function

   subroutine opt_set_value(this, val)
      implicit none
      class(opt_t) :: this
      character(len=*), intent(in) :: val
      integer :: sz
      
      if (this%stat /= 1) then
         this%stat = 1
 !        if (allocated(this%val)) then
 !           deallocate(this%val)
 !           allocate(character(len=len_trim(val)) :: this%val)
 !        end if
         this%val = trim(val)
      end if
   end subroutine

   subroutine opt_write(this, outp)
      implicit none
      class(opt_t), intent(in) :: this
      integer, intent(in) :: outp

      if ( len_trim(this%key) /= 0) then
         write (outp, '(a,a)', advance='no') '-', trim(this%key)
         if (len_trim(this%des) /= 0) then
            write (outp, '(x,a)', advance='no') trim(this%des)
         end if
         if (this%req > 0) then
            write (outp, '(a)', advance='no') ' required'
         end if
         write (outp, *)
         if ( this%stat /= -1) then
             write (outp, '(3x,a)', advance='no') trim(this%val)
             if (this%stat > 0) then
                write (outp, '(a)', advance='no') ' set by user'
             else
                write (outp, '(a)', advance='no') ' default'
             end if
         else
             write (outp, '(a)', advance='no') ' no value set'
         end if
         write (outp, *)
      else
         write(error_unit,*)'attempt to write unitialized option'
      end if
   end subroutine

   subroutine prog_opts_write(this, outp)
      implicit none
      class(prog_opts_t), intent(in) :: this
      integer, intent(in) :: outp
      integer :: i

      write (outp, *) 'there are ',this%nopts,' options'
      write (outp, *) 'options are :'
      do i = 1, this%nopts
         call this%opts(i)%write(outp)
      end do
      write (outp, *)
   end subroutine

   logical function opt_valid(this)
     implicit none
     class(opt_t) :: this
    
     opt_valid = .true.
     if(this%stat == -1) then
       opt_valid = .false.
     else 
       if(len_trim(this%val) .eq. 0) then
         opt_valid = .false.
       end if
     end if
   end function

   subroutine prog_opts_init(this)
      implicit none
      class(prog_opts_t) :: this
!      class(opt_t), dimension(:), allocatable, target :: opts
      
      this%ncap = 16
      this%nopts = 0
      allocate(this%opts(this%ncap))
      this%allow_extra_options = 0
      this%allow_extra_args = 0
   end subroutine

   subroutine prog_opts_final(this)
      implicit none
      type(prog_opts_t) :: this 
      deallocate (this%opts)
      this%opts => null()
      this%ncap = 0
      this%nopts = 0
   end subroutine

   subroutine prog_opts_add_option(this, key, des, val, req)
      implicit none
      class(prog_opts_t) :: this
      character(len=*), intent(in) :: key, des, val
      integer, intent(in), optional :: req
      type(opt_t), dimension(:), allocatable, target :: new_opts
      integer :: i

      write(error_unit,*)'# of opts = ',this%nopts,' cap = ',this%ncap
      if (this%ncap == this%nopts) then
         write(error_unit,*)'reallocating'
         this%ncap = this%ncap + this%ncap
         allocate (new_opts(this%ncap))
         do i = 1, this%nopts
            new_opts(i) = this%opts(i)
         end do
         deallocate (this%opts)
         this%opts => new_opts
      end if
      this%nopts = this%nopts + 1
      call this%opts(this%nopts)%init(key, des, val, req)
! check      
      call this%opts(this%nopts)%write(error_unit)
   end subroutine

   logical function prog_opts_has_value(this, key) result(fnd)
      implicit none
      class(prog_opts_t) :: this
      character(len=*), intent(in) :: key
      integer :: i

      fnd = .false.
      do i = 1, this%nopts
         if ( this%opts(i)%matchs(key) ) then
            fnd = this%opts(i)%valid()
         end if
      end do
   end function

   logical function prog_opts_has_option(this, key) result(fnd)
      implicit none
      class(prog_opts_t) :: this
      character(len=*), intent(in) :: key
      integer :: i
      fnd = .false.
      do i = 1, this%nopts
         if ( this%opts(i)%matchs(key) ) then
            fnd = .true.
            exit
         end if
      end do
   end function

   subroutine prog_opts_parse_command_line(this)
      implicit none
      class(prog_opts_t) :: this
      integer :: i, pos, narg
      character(len=255) :: arg, arg2, key, val
      logical :: next_line = .false.
 
      narg = command_argument_count()
      write(error_unit,*)'# of args = ',narg 
      do i= 1,narg
         arg = repeat(' ',LEN(arg))
         call get_command_argument(i, arg)
         write(error_unit,*) i,' ',trim(arg)
      end do
      do i = 1, narg
         if (next_line) then
            write(error_unit,*)' pcl ',arg,' ',arg2    
            call this%set_value(arg, arg2)
            next_line = .false.
            arg2(1:) = repeat(' ',LEN(arg2))
            cycle
         end if
         arg = repeat(' ',LEN(arg))
         call get_command_argument(i, arg)
         if (scan(arg, '-') == 1) then
            arg = trim(arg(2:))  
            if ( scan(arg,'-') == 1) then
               arg = trim(arg(2:))
            end if
            write(error_unit,*)' arg ', i, ' ' , arg
            pos = scan(arg, '=')
            if (pos == 0) then
               if (i /= narg) then
                  call get_command_argument(i + 1, arg2)
                  if (scan(arg2, '-') /= 1) then
                     next_line = .true.
                  else
                     arg2 = ' '
                     write(error_unit,*)' pcl ',arg,' def true' 
                     call this%set_value(arg, 'true')
                  end if
               else
                  write(error_unit,*)' pcl ',arg,' def true'
                  call this%set_value(arg, 'true')
               end if
            else
               key = ' '
               val = ' '
               key = trim(arg(1:(pos - 1)))
               val = trim(arg((pos + 1):len_trim(arg)))
               write(error_unit,*)' pcl ',key,' ',val
               call this%set_value(key,val)
            end if
         else
            write (error_unit, *) 'parse error command line'
            stop
         end if
      end do
      write(error_unit,*)' parsed command line'
   end subroutine
   
   subroutine prog_opts_set_value(this, key, val)
      implicit none
      class(prog_opts_t) :: this
      character(len=*), intent(in) :: key
      character(len=*), intent(in) :: val
      integer :: i
      logical :: fnd = .false.

      do i = 1, this%nopts
         if ( this%opts(i)%matchs(key) ) then
            write(error_unit,*)'po set_value ',trim(key),' ',trim(val)
            call this%opts(i)%set_value(val)
            fnd = .true.
            exit
         end if
      end do
      if (fnd .eqv. .false.) then
         write (error_unit, *) 'set value no such option as ', trim(key)
         if (this%allow_extra_options == 0) then
            stop -1
         end if
      end if
   end subroutine

   subroutine prog_opts_get_value(this, key, v)
      implicit none
      class(prog_opts_t), intent(in) :: this
      character(len=*), intent(in) :: key
      character(len=*), intent(out) :: v
      integer :: i
      logical :: fnd = .false.

      v = ''
      do i = 1, this%nopts
         if ( this%opts(i)%matchs(key) ) then
            if ( this%opts(i)%valid() ) then
                write(error_unit,*)' value of ',trim(key),' is ',this%opts(i)%val
                v = trim(this%opts(i)%val)
                fnd = .true.
                exit
            else
                write(error_unit,*)' the option ',trim(key),' has not valid value'
                stop -1 
            end if
         end if
      end do
      if (fnd .eqv. .false.) then
         write (error_unit, *) 'get value no such option as ', trim(key)
         stop
      end if
   end subroutine

   subroutine prog_opts_parse_file(this, inp_name)
      implicit none
      class(prog_opts_t) :: this
      character(len=*), intent(in) :: inp_name
      character(len=LINE_SZ) :: line_str
      character(len=:), dimension(:), allocatable :: tokens
      character(len=LINE_SZ) :: des
      integer :: req
      integer :: p, q, nt, r
      integer, parameter :: inp = 98
      logical :: file_exists
      character(len=128) :: acc
	
      inquire(file=inp_name,exist=file_exists)
      if ( file_exists .eqv. .false.) then
         write(error_unit,*)inp_name,' does not exist'
         stop -1
      end if
      open (unit=inp, file=inp_name, action='read', status='old')
      do
         line_str = repeat(' ',LINE_SZ)
         read(inp,'(a)',iostat=r)line_str
         if (r == iostat_end) then
            exit
         end if
         if ( r /= 0) then
            write(error_unit,*)'reading error in config file'
            stop -1
         end if
         if ( len_trim(line_str) == 0) then
            write(error_unit,*)'found blank line in config file'
            exit
         end if
         p = scan(line_str,'=')
         if ( p /= 0) then
             line_str(p:p) = ' '
         end if
         nt = tokenize_string(line_str,' ',tokens)
         select case (nt)
           case (0)
              write(error_unit,*)'blank line'
              stop
           case (1) 
              call this%set_value(trim(tokens(1)),'true')
           case (2)
	      call this%set_value(trim(tokens(1)),trim(tokens(2)))
           case default
              write(error_unit,*)'too many tokens on line in options file'
              stop            
         end select 
      end do
      close (inp)
   end subroutine

   subroutine prog_opts_add_options_from_file(this,inp_name)
      implicit none
      class(prog_opts_t) :: this
      character(len=*), intent(in) :: inp_name
      character(len=:),dimension(:),allocatable :: tokens
      character(len=*),parameter :: default = 'false'
      character(len=LINE_SZ) :: des
      integer,parameter :: inp = 34
      integer :: r,k,nt
      logical :: file_exists
      character(len=128) :: acc
      
      inquire(file=inp_name,exist=file_exists)
      if ( file_exists .eqv. .false.) then
         write(error_unit,*)inp_name,' does not exist'
         stop -1
      end if
      open (unit=inp, file='opts.in', action='read', status='old')
      do
         des = ' '
         r = process_input(inp, des, tokens, nt)
         if (r == iostat_end) then
            exit
         end if
         write (error_unit, *)' in prog_opts_add_options_from_file'
         write (error_unit, *) 'des = (', trim(des), ')'
         do k = 1, nt
            write (error_unit, *) k, ' (', tokens(k), ')'
         end do
         select case (nt)
           case (0)
              write(error_unit,*)'blank line'
              stop
           case (1) 
              call this%add_option(trim(tokens(1)),des,'',0)
           case (2)
              if ( tokens(2) == 'required') then
	         call this%add_option(trim(tokens(1)),des,'',1)
              else
	         call this%add_option(trim(tokens(1)),des,tokens(2),0)
              end if  
           case (3)
              if ( tokens(3) == 'required') then
	         call this%add_option(trim(tokens(1)),des,tokens(2),1)
              else
	         write(error_unit,*)'parse error in options file'
	         stop
              end if  
           case default
              write(error_unit,*)'too many tokens on line in options file'
              stop            
         end select
         deallocate(tokens) 
      end do
      write(error_unit,*)'added options from file' 
      close(inp)
   end subroutine

   integer function process_input(inp, des, tokens, nt) result(io)
      implicit none
      integer, intent(in) :: inp
      character(len=LINE_SZ), intent(inout) :: des
      character(len=:), dimension(:), allocatable, intent(out) :: tokens
      integer, intent(out) :: nt
      character(len=LINE_SZ) :: line_str
      character(len=LINE_SZ) :: str
      character(len=LINE_SZ) :: fstr
      character(len=LINE_SZ) :: ostr
      integer :: p, q, k, sz

      des = repeat(' ',LINE_SZ)
      line_str = repeat(' ',LINE_SZ)
      fstr = line_str
      ostr = line_str
      str = line_str
      read (inp, '(a)', iostat=io) line_str
      write (error_unit, *) 'read ', len_trim(line_str), ' chars'
      if (io == iostat_end) then
         write (error_unit, *) 'eof reached'
         return
      end if
      if (io > 0) then
         write (error_unit, *) 'reading error'
         stop -1
      end if
      str = trim(line_str)
      sz = len_trim(line_str)
      if (sz == 0) then
         io = iostat_end
         write (error_unit, *) 'blank line detected'
         return
      end if
      write (error_unit, *) 'input string = (', trim(line_str), ')'
      p = scan(str, '"')
      if (p /= 0) then
         q = scan(str(p + 1:sz), '"')
         write (error_unit, *) 'p = ', p, ' q= ', q, ' p+q = ', (p + q), ' sz = ', sz
         if (q == 0) then
            write (error_unit, *) 'unterminated string'
            stop - 1
         end if
         fstr = trim(str(1:(p - 1)))
         des = trim(str((p + 1):(p + q - 1)))
         if ((p + q + 1) > sz) then
            ostr = trim(fstr)
         else
            ostr = trim(fstr)//' '//trim(str((p + q + 1):sz))
         end if
      else
         ostr = trim(str)
      end if
      write (error_unit, *) 'ostr = (', trim(ostr), ')'
      nt = tokenize_string(ostr, ' ', tokens)
      if ( len_trim(des) /= 0) then
         write (error_unit, *) 'description = (', trim(des), ')'
      else
         write (error_unit, *) 'description = (none)'
      end if
      do k = 1, nt
         write (error_unit, *) k, ' (', trim(tokens(k)), ')'
      end do
      write (error_unit, *) 'next input'
   end function
end module
