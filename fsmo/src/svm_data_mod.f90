module svm_data_mod
use,intrinsic :: iso_fortran_env
use svm_utils_mod

  type svm_data_t
    real(real64),dimension(:,:),allocatable :: vecs
    real(real64),dimension(:),allocatable :: y
    integer(int32) :: nvecs,nfeat
  contains 
    final svm_data_final
    procedure :: read => svm_data_read
    procedure :: write => svm_data_write
    procedure :: translate => svm_data_translate
    procedure :: libsvm_read => svm_data_libsvm_read
    procedure :: libsvm_write => svm_data_libsvm_write
    procedure :: tdo_read => svm_data_tdo_read
    procedure :: tdo_write => svm_data_tdo_write
  end type
contains
   subroutine svm_data_final(this)
     implicit none
     type(svm_data_t),intent(inout) :: this
     if (allocated(this%vecs) .eqv. .true.) deallocate(this%vecs)
     if (allocated(this%y) .eqv. .true.) deallocate(this%y)
     this%nvecs = 0
     this%nfeat = 0
   end subroutine

   subroutine svm_data_tdo_read(this,inp_name)
     implicit none
     class(svm_data_t) :: this
     character(len=*),intent(in) :: inp_name
     integer,parameter :: inp = 43    
     integer :: i,j 

     open(unit=inp,file=inp_name,status='old',action='read',ACCESS='STREAM')     
     if (allocated(this%vecs) .eqv. .true.) deallocate(this%vecs)
     if (allocated(this%y) .eqv. .true.) deallocate(this%y)
     read(inp)this%nvecs
     read(inp)this%nfeat
     allocate(this%y(this%nvecs))
     allocate(this%vecs(this%nfeat,this%nvecs))
     read(inp)this%y
     do i=1,this%nvecs
        do j=1,this%nfeat
          read(inp)this%vecs(j,i)
        end do
     end do
     close(inp)
   end subroutine
   
   subroutine svm_data_tdo_write(this,outp_name)
     implicit none
     class(svm_data_t),intent(in) :: this
     character(len=*),intent(in) :: outp_name
     integer,parameter :: outp = 44 
     integer :: i,j 
     
     open(unit=outp,file=outp_name,action='write',access='STREAM')
     write(outp)this%nvecs
     write(outp)this%nfeat
     write(outp)this%y
     do i=1,this%nvecs
        do j=1,this%nfeat
          write(outp)this%vecs(j,i)
        end do
     end do
     close(outp)
   end subroutine
  
   subroutine svm_data_libsvm_read(this,inp_name)
     implicit none
     class(svm_data_t) :: this
     character(len=*),intent(in) :: inp_name
     integer,parameter :: inp = 43
     character(len=65536) :: linestr
     character(len=:),dimension(:),allocatable :: tokens
     character(len=*),parameter :: delims =' :\n\t'
     integer :: i,j,nfmax,ix,ntokens,io,nf
               
     open(unit=inp,file=inp_name,status='old',action='read')
     this%nvecs = 0
     this%nfeat = 0
     nfmax = 0
     do while (.true.)
        read(inp,'(a)',iostat=io)linestr
        if (io == iostat_end) then
           exit
        end if
        if (io /= 0) then
           write(error_unit,*)'error reading libsvm file'
           stop -1
        end if    
        if ( len_trim(linestr) == 0) then
           exit 
        end if      
        ntokens = tokenize_string(linestr,delims,tokens)
        if ( (ntokens/2) == 0 .and. ntokens/=1) then
            write(error_unit,*)'file format error in libsvm file!\n'
            stop
        end if
        this%nvecs = this%nvecs + 1
        nfmax = max(nfmax,ntokens) 
        deallocate(tokens)
     end do
     rewind(inp)
     nfmax = (nfmax-1)/2   
     if ( nfmax == 0) then
        write(error_unit,*)'No vectors found in libsvm file!\n'
        stop
     end if 
     this%nfeat = nfmax
     if (allocated(this%vecs) .eqv. .true.) deallocate(this%vecs)
     if (allocated(this%y) .eqv. .true.) deallocate(this%y)
     allocate(this%y(this%nvecs))
     allocate(this%vecs(nfmax,this%nvecs))
     write(error_unit,*)'found ',this%nvecs,' vectors of length ',this%nfeat
     do i=1,this%nvecs
        do j=1,this%nfeat
           this%vecs(j,i) = 0.0d0
        end do
     end do
     do i=1,this%nvecs
        read(inp,'(a)',iostat=io)linestr
        if (io /= 0) then
           exit
        end if          
        ntokens = tokenize_string(linestr,delims,tokens)
        nf = (ntokens-1)/2
!        write(error_unit,*)'vec = ',i,' ',tokens(1)
        read(tokens(1),*)this%y(i)
!        write(error_unit,*)'label = ',this%y(i)
        do j=2,ntokens,2
!           write(error_unit,*)i,' ',j/2,' ',trim(tokens(j)),' ',trim(tokens(j+1))
           read(tokens(j),*)ix
           read(tokens(j+1),*)this%vecs(ix,i)
        end do               
     end do
     close(inp)
   end subroutine
   
   subroutine svm_data_libsvm_write(this,outp_name)
     class(svm_data_t),intent(in) :: this
     character(len=*),intent(in) :: outp_name
     integer,parameter :: outp = 44  
     open(unit=outp,file=outp_name,action='write')

     write(error_unit,*)'svm write ',this%nvecs,' ',this%nfeat
     do i=1,this%nvecs
        write(outp,'(d20.10)',advance='no')this%y(i)
 !       write(error_unit,'(d20.10)',advance='no')this%y(i)
        do j=1,this%nfeat
           if ( dabs(this%vecs(j,i)) > 1.d-12) then
              write(outp,'(1x,i12,a1,d20.10)',advance='no')j,':',this%vecs(j,i)
 !             write(error_unit,'(1x,i12,a1,d20.10)',advance='no')j,':',this%vecs(j,i)
           end if
        end do
        write(outp,*)' '
        write(error_unit,*)' '
     end do
     close(outp)
   end subroutine	   

   subroutine svm_data_read(this,outp_name)
     implicit none
     class(svm_data_t),intent(in) :: this
     character(len=*),intent(in) :: outp_name
     logical :: find_suffix = .false.
     
     find_suffix = is_tdo_name(outp_name)
     if ( find_suffix) then
         write(error_unit,*)' reading tdo file'
         call this%tdo_read(outp_name)
     else
         write(error_unit,*)' reading libsvm file'
         call this%libsvm_read(outp_name)
     end if     
   end subroutine
   
   subroutine svm_data_write(this,inp_name)
     implicit none
     class(svm_data_t),intent(in) :: this
     character(len=*),intent(in) :: inp_name
     logical :: find_suffix = .false.

     find_suffix = is_tdo_name(inp_name)
     if ( find_suffix) then
         call this%tdo_write(inp_name)
     else
         call this%libsvm_write(inp_name)
     end if     
   end subroutine

   subroutine svm_data_translate(this,inp_name)
     implicit none
     class(svm_data_t) :: this
     character(len=*),intent(in) :: inp_name
     logical :: find_suffix = .false.
     character(len=4) :: suffix
     character(len=255) :: outp_name
     integer :: slen

     slen = len(inp_name)
     if ( slen > 4) then
        suffix = outp_name(slen-3:slen)
        if ( suffix == '.tdo') then
             find_suffix = .true.
        end if
     end if
     if ( find_suffix) then
	 outp_name = outp_name(1:slen-4)	     
         call this%tdo_read(inp_name)
         call this%libsvm_write(outp_name)
     else
         outp_name = inp_name // '.tdo'
         call this%libsvm_read(inp_name)
         call this%tdo_write(outp_name)
     end if        
   end subroutine
 
   logical function is_tdo_name(str) result(f)
     character(len=*),intent(in) :: str
     character(len=4) :: suffix
     integer :: slen
     f = .false.
     slen = len_trim(str)
     if ( slen > 3) then
        suffix = trim(str((slen-3):slen))
        if ( suffix == '.tdo') then
           f = .true.
        end if
     end if
   end function
end module
