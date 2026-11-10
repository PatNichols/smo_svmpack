
FFLAGS='-g'

all:
	${F90} -c ${FFLAGS} svm_utils_mod.f90
	${F90} -c ${FFLAGS} svm_timer_mod.f90
	${F90} -c ${FFLAGS} prog_opts_mod.f90
	${F90} -c ${FFLAGS} svm_data_mod.f90
	${F90} -c ${FFLAGS} svm_options_mod.f90
	${F90} -c ${FFLAGS} svm_kernel_matrix_mod.f90
	${F90} -c ${FFLAGS} smo_solver_mod.f90
	${F90} -c ${FFLAGS} svm_model_mod.f90
	${F90} -c ${FFLAGS} svm_classify_mod.f90
	${F90} -c ${FFLAGS} svmtask.f90
	${F90} -c ${FFLAGS} tst_svm.f90
	${F90} ${FFLAGS} -o svmtask svmtask.o svm_classify_mod.o svm_model_mod.o smo_solver_mod.o svm_kernel_matrix_mod.o svm_options_mod.o svm_data_mod.o prog_opts_mod.o svm_timer_mod.o svm_utils_mod.o
	${F90} ${FFLAGS} -o tst_svm tst_svm.o svm_classify_mod.o svm_model_mod.o smo_solver_mod.o svm_kernel_matrix_mod.o svm_options_mod.o svm_data_mod.o prog_opts_mod.o svm_timer_mod.o svm_utils_mod.o
        
clean:
	rm -f *.o *.mod svmtask svm_test
		