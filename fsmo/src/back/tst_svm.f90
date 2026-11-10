program tst_svm_read
use,intrinsic :: iso_fortran_env
use svm_options_mod
use smo_solver_mod
use svm_data_mod

type(svm_options_t) :: opts
type(svm_data_t) :: data,data2
type(smo_solver_t) :: solver
integer :: r

    call opts%init()
    call data%read(opts%data_file)
    call data%tdo_write('data_data.tdo') 
    r =  solver%train(opts,data)
    call compare_data(data)
    call compare_data2(data)
contains
    subroutine compare_data(data)
       implicit none
       type(svm_data_t),intent(in) :: data
       type(svm_data_t) :: data2,data3
       integer :: i,j,cnt
       real(kind=real64) :: diff
              
       write(output_unit,*)' compare1 '
       call data2%read('data_data.tdo')
       call data3%read('data_kmat.tdo')
       write(output_unit,*)'# of vectors  = ',data%nvecs,' ',data2%nvecs
       write(output_unit,*)'# of features = ',data%nfeat,' ',data2%nfeat
       if ( data2%nvecs /= data3%nvecs .or. data2%nvecs /= data%nvecs) then
            write(error_unit,*)'# of vectors does not agree'
            write(error_unit,*) data%nvecs,' ',data2%nvecs,' ',data3%nvecs
       end if       
       if ( data2%nfeat /= data3%nfeat .or. data2%nfeat /= data%nfeat) then
            write(error_unit,*)'# of vectors does not agree'
            write(error_unit,*) data%nfeat,' ',data2%nfeat,' ',data3%nfeat
       end if       
       do i=1,data%nvecs
            diff = dabs(data%y(i) - data2%y(i))
            if (diff > 1.e-8) then
               write(error_unit,*)'label2 ',i,' ',data%y(i),' ',data2%y(i)
            end if
       end do
       do i=1,data%nvecs
            diff = dabs(data%y(i) - data3%y(i))
            if (diff > 1.e-8) then
               write(error_unit,*)'label3 ',i,' ',data%y(i),' ',data3%y(i)
            end if
       end do
       do i=1,data%nvecs
          do j=1,data%nfeat
            diff = dabs(data%vecs(j,i) - data2%vecs(j,i))
            if (diff > 1.e-8) then
               write(error_unit,*)'vecs2 ',i,' ',j,' ',data%vecs(j,i),' ',data2%vecs(j,i)
            end if
          end do
       end do
       do i=1,data%nvecs
          do j=1,data%nfeat
            diff = dabs(data%vecs(j,i) - data3%vecs(j,i))
            if (diff > 1.e-8) then
               write(error_unit,*)'vecs2 ',i,' ',j,' ',data%vecs(j,i),' ',data3%vecs(j,i)
            end if
          end do
       end do        	

       cnt = 0
       do i=1,data%nvecs
          do j=1,data%nfeat
             if ( dabs(data%vecs(j,i)) > 1.d-8) cnt = cnt + 1
          end do
       end do
       write(output_unit,*)'origianl # of non zero features = ',cnt
       cnt = 0
       do i=1,data%nvecs
          do j=1,data%nfeat
             if ( dabs(data2%vecs(j,i)) > 1.d-8) cnt = cnt + 1
          end do
       end do
       write(output_unit,*)'data2 # of non zero features = ',cnt       
    end subroutine

    subroutine compare_data2(data)
       implicit none
       type(svm_data_t),intent(in) :: data
       type(svm_data_t) :: data2
       integer :: i,j
       real(kind=real64) :: diff
              
       write(output_unit,*)' compare2 '
       write(output_unit,*)' # of vectors = ',data%nvecs
       call data%write('data_data.lsvm')
       call data2%read('data_data.lsvm')
       write(output_unit,*)'# of vectors  = ',data%nvecs,' ',data2%nvecs
       write(output_unit,*)'# of features = ',data%nfeat,' ',data2%nfeat
       if ( data2%nvecs /= data%nvecs ) then
            write(error_unit,*)'# of vectors does not agree'
            write(error_unit,*) data%nvecs,' ',data2%nvecs
       end if       
       if ( data%nfeat /= data%nfeat) then
            write(error_unit,*)'# of vectors does not agree'
            write(error_unit,*) data%nfeat,' ',data2%nfeat
       end if       
       do i=1,data%nvecs
            diff = dabs(data%y(i) - data2%y(i))
            if (diff > 1.e-8) then
               write(error_unit,*)'label2 ',i,' ',data%y(i),' ',data2%y(i)
            end if
       end do
       do i=1,data%nvecs
          do j=1,data%nfeat
            diff = dabs(data%vecs(j,i) - data2%vecs(j,i))
            if (diff > 1.e-8) then
               write(error_unit,*)'vecs2 ',i,' ',j,' ',data%vecs(j,i),' ',data2%vecs(j,i)
            end if
          end do
       end do
    end subroutine
end program