module svm_model_mod
use,intrinsic :: iso_fortran_env

  type :: svm_model_t 
    real(kind=real64),dimension(:),allocatable :: yalf
    real(kind=real64),dimension(:,:),allocatable :: vecs
    real(kind=real64) :: bias,c1,c2
    integer(kind=4) :: kernel_scale,kernel_type,kernel_power
    integer(kind=int32) :: nvecs,nfeat
  contains
    final svm_model_final
    procedure :: init => svm_model_init
    procedure :: score => svm_model_score
    procedure :: scale => svm_model_scale   
  end type
  
contains
  subroutine svm_model_final(this)
    implicit none
    type(svm_model_t) :: this  
    deallocate(this%vecs)
    deallocate(this%yalf)
    this%nvecs = 0
    this%nfeat = 0
  end subroutine

  subroutine svm_model_init(this,inp_name)
    implicit none
    class(svm_model_t) :: this
    character(len=*) :: inp_name
    integer :: io
    integer,parameter :: inp = 45
    
    open(unit=inp,file=inp_name,action='read',status='old',iostat=io)
    read(inp)this%nvecs
    read(inp)this%nfeat
    read(inp)this%kernel_type
    read(inp)this%kernel_power
    read(inp)this%kernel_scale
    read(inp)this%c1
    read(inp)this%c2
    read(inp)this%bias
    allocate(this%yalf(this%nvecs))
    allocate(this%vecs(this%nfeat,this%nvecs))
    read(inp)this%yalf
    read(inp)this%vecs
    close(inp)  
  end subroutine


  real(kind=real64) function svm_model_scale(this,v) result(scal)
    implicit none
    class(svm_model_t),intent(in) :: this
    real(kind=real64),intent(in) :: v(*)
    integer :: i     
    
    scal = 1.0d0  
    if ( this%kernel_scale /= 0) then
       scal = 0.0d0
       do i=1,this%nvecs
          scal = scal + v(i) * v(i)
       end do
       scal = this%c1 * scal + this%c2
       select case (this%kernel_type)
         case (0)  
              scal = scal
         case (1)  
              scal = scal ** this%kernel_power   
         case (2)  
              scal = 1.0d0
         case (3)  
              scal = datanh(scal)
         case default
             write(error_unit,*)'bad kernel type' 
             stop -1
       end select
       scal = 1.0d0/dsqrt(scal)
    end if
  end function
    
  real(kind=real64) function svm_model_score(this,v) result(scor)
    implicit none
    class(svm_model_t),intent(in) :: this
    real(kind=real64),intent(in) :: v(*)
    real(kind=real64) :: scal = 1.0d0
    real(kind=real64) :: sum,t
    integer :: i,j       
  
    if ( this%kernel_type /= 2) then
       scal = this%scale(v)
       do i=1,this%nvecs
          sum = 0.0d0
          do j=1,this%nfeat
             sum = sum + this%vecs(j,i) * v(j)
          end do
          sum = this%c1 * sum + this%c2
          select case (this%kernel_type)
            case (0)  
            case (1)  
              sum = sum ** this%kernel_power   
            case (2)  
              sum = 1.0d0
            case (3)  
              sum = datanh(sum)
            case default 
              write(error_unit,*)'bad kernel type'
              stop -1
          end select
          sum = sum * this%yalf(i) * scal
          scor = scor + sum 
       end do
       scor = scor + this%bias
    else
       do i=1,this%nvecs
          sum = 0.0d0
          do j=1,this%nfeat
             t = this%vecs(j,i) - v(j)
             sum = sum + t * t
          end do
          sum = this%c1 * sum
          sum = dexp(-sum) 
          sum = sum * this%yalf(i)
          scor = scor + sum 
       end do
       scor = scor + this%bias       
    end if 
  end function   

end module