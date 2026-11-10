program svmtask
use svm_options_mod
use smo_solver_mod
use svm_classify_mod
use svm_data_mod

type(svm_options_t) :: opts
type(svm_data_t) :: data
type(smo_solver_t) :: solver
integer :: r

    call opts%init()
    if (opts%task == 0) then
       call data%read(opts%data_file)
       r =  solver%train(opts,data)
    else
       call svm_validate(opts)
    end if
end program