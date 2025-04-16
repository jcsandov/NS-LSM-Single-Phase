# This bash file compiles all the necessary routines to generate output files
gfortran -c f2kcli.f90
gfortran -c str_int.f90

gfortran -o mktseries get_tseries.f90 f2kcli.o str_int.o
gfortran -o mk_solu.restart make_solu.restart.f90 f2kcli.o str_int.o
gfortran -o mk_phi.restart make_phi.restart.f90 f2kcli.o str_int.o
gfortran -o mk_output_tecfiles postproc_output.f90 tecio64.a f2kcli.o str_int.o -lstdc++