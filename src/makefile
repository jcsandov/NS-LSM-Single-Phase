# ================== CONFIGURATION ==================

CONFIG ?= fondecyt

# Compiler
ifeq ($(CONFIG),fondecyt)
    F90COMP = mpif90
    FFLAGS = -O3 -fbounds-check -g -fallow-argument-mismatch -ffpe-trap=invalid,zero,overflow
    LIBS = -L/usr/local/lib -llapack -lblas
endif

ifeq ($(CONFIG),leftraru)
    F90COMP = mpif90
    FFLAGS = -O0 -g
    OPENBLAS_DIR = $(shell spack location -i 'openblas%intel@18.0.5')
    LIBS = -L$(OPENBLAS_DIR)/lib -lopenblas -lpthread -lm
endif

# ================== OBJECTS ==================

OBJS = precision.o global.o global_param.o global_app.o global_mpi.o\
	global_osg.o checksum.o wf_mpi.o global_lsm.o global_obstacle.o global_debug.o\
	phi_inlet.o\
	NewtonRaphsonSolver.o\
	TetrahedronMethods.o\
	DataTypes.o\
	AdvectionMethods.o\
	InterpolationMethods.o\
	pascal_tdma.o\
	main.o \
	init.o \
	mg_nlevel.o \
	mg_driver.o \
	mg_inject.o \
	mg_metrics.o \
	solver_daf.o \
	bcond_fm.o \
	les_dynamic_smagorinsky.o \
	calc_phi_gradient.o \
	calc_phi_gradient_ENO2.o \
	calc_h.o \
	calc_h_gradient_ENO2.o \
	ghost_fluid_extrapolation.o \
	des_eddy.o \
	sgtsv.o \
	xerbla.o \
	rsmooth.o \
	init_LSM.o \
	levelsetmethod.o\
	init_obstacle.o\
	bcond_obstacle.o\
	phi_obstacle.o\
	conservation.o

# ================== MODULE GROUPS ==================

PACKMODULE  = global.mod global_param.mod global_mpi.mod
PACKMODULE2 = global_param.mod global_mpi.mod global_app.mod
PACKMODULE3 = global.mod global_param.mod global_app.mod global_mpi.mod global_osg.mod
PACKMODULE4 = ad_RKTVD3.F90 calc_RH_AD_WENO3.F90 calc_RH_AD_WENO3_test2.F90 calc_RH_AD.F90 bcond_lsm.F90 bcond_lsm_obstacle.F90 \
			  reinitialization_benchmark.F90 reinitialisation_test2.F90 calc_sign_ini.F90 calc_sign_ini_2.F90 calc_RH_RN.F90 calc_RH_RN_2.F90 \
			  calc_lambda_RN.F90 narrowband.F90 calc_NormGrad_d.F90 calc_sign_ini_obstacle.F90 calc_RH_RN_obstacle.F90

PACKMODULE5 = near_free_surface_q_update.F90 VelocityGradientTensor.F90 PhiGradientVector.F90 \
			  pressure_extrapolation.F90 p_correction_post_advection.F90

PACKMODULE6 = geometric_reinitialisation.F90 GetTetrahedraList.F90 GetIsosurfaceAndDistances.F90 \
			  PiecewiseConstantFunction.F90 OrthogonalProjection.F90 NeighboursPhiCorrection2.F90 \
			  WriteTriangulation.F90 TotalWaterVolumeGeom.F90 TotalWaterVolumeSussman.F90 TotalWaterVolume.F90 \
			  WriteTotalWaterVolume.F90

# ================== TARGET ==================

channel : $(OBJS)
	$(F90COMP) $(FFLAGS) -o channel $(OBJS) $(LIBS)

# ================== COMPILE RULES ==================

precision.o precision.mod : precision.F90
	$(F90COMP) $(FFLAGS) -c precision.F90 

global.o global.mod : global.F90 precision.mod 
	$(F90COMP) $(FFLAGS) -c global.F90

global_app.o global_app.mod : global_app.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_app.F90 

global_mpi.o global_mpi.mod : global_mpi.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_mpi.F90 

global_param.o global_param.mod : global_param.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_param.F90

global_osg.o global_osg.mod : global_osg.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_osg.F90

checksum.o checksum.mod : checksum.F90 precision.mod global.mod global_mpi.mod
	$(F90COMP) $(FFLAGS) -c checksum.F90 

wf_mpi.o wf_mpi.mod : wf_mpi.F90 global.mod global_param.mod global_mpi.mod
	$(F90COMP) $(FFLAGS) -c wf_mpi.F90

global_lsm.o global_lsm.mod : global_lsm.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_lsm.F90 

global_obstacle.o global_obstacle.mod : global_obstacle.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_obstacle.F90 

global_debug.o global_debug.mod : global_debug.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c global_debug.F90 

phi_inlet.o phi_inlet.mod: phi_inlet.F90 $(PACKMODULE) global_lsm.mod precision.mod
	$(F90COMP) $(FFLAGS) -c phi_inlet.F90

TetrahedronMethods.o TetrahedronMethods.mod : TetrahedronMethods.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c TetrahedronMethods.F90 

DataTypes.o DataTypes.mod : DataTypes.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c DataTypes.F90 

NewtonRaphsonSolver.o NewtonRaphsonSolver.mod : NewtonRaphsonSolver.F90 precision.mod DataTypes.mod TetrahedronMethods.mod
	$(F90COMP) $(FFLAGS) -c NewtonRaphsonSolver.F90 

AdvectionMethods.o AdvectionMethods.mod : AdvectionMethods.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c AdvectionMethods.F90 

InterpolationMethods.o InterpolationMethods.mod : InterpolationMethods.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c InterpolationMethods.F90 

pascal_tdma.o pascal_tdma.mod : pascal_tdma.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c pascal_tdma.F90 

SolitaryWave.o SolitaryWave.mod : SolitaryWave.F90 precision.mod
	$(F90COMP) $(FFLAGS) -c SolitaryWave.F90 

main.o : main.F90 $(PACKMODULE) global_lsm.mod global_debug.mod
	$(F90COMP) $(FFLAGS) -c main.F90

init.o : init.F90 $(PACKMODULE3) checksum.mod
	$(F90COMP) $(FFLAGS) -c init.F90

mg_driver.o : mg_driver.F90 $(PACKMODULE) global_osg.mod global_app.mod wf_mpi.mod global_lsm.mod global_debug.mod
	$(F90COMP) $(FFLAGS) -c mg_driver.F90

mg_nlevel.o : mg_nlevel.F90 global.mod global_param.mod global_mpi.mod
	$(F90COMP) $(FFLAGS) -c mg_nlevel.F90

mg_inject.o : mg_inject.F90 global.mod global_mpi.mod
	$(F90COMP) $(FFLAGS) -c mg_inject.F90

mg_metrics.o : mg_metrics.F90 $(PACKMODULE) global_app.mod
	$(F90COMP) $(FFLAGS) -c mg_metrics.F90

solver_daf.o : solver_daf.F90 $(PACKMODULE2) $(PACKMODULE5) checksum.mod global_debug.mod global_lsm.mod global_obstacle.mod
	$(F90COMP) $(FFLAGS) -c solver_daf.F90

bcond_fm.o : bcond_fm.F90 global_param.mod global_app.mod global_mpi.mod
	$(F90COMP) $(FFLAGS) -c bcond_fm.F90

les_dynamic_smagorinsky.o : les_dynamic_smagorinsky.F90 global_app.mod global_mpi.mod
	$(F90COMP) $(FFLAGS) -c les_dynamic_smagorinsky.F90

ghost_fluid_extrapolation.o : ghost_fluid_extrapolation.F90 $(PACKMODULE5) global_app.mod global_mpi.mod global_lsm.mod
	$(F90COMP) $(FFLAGS) -c ghost_fluid_extrapolation.F90 

calc_phi_gradient.o : calc_phi_gradient.F90 global_app.mod global_mpi.mod 
	$(F90COMP) $(FFLAGS) -c calc_phi_gradient.F90 

calc_phi_gradient_ENO2.o : calc_phi_gradient_ENO2.F90 global_app.mod global_mpi.mod 
	$(F90COMP) $(FFLAGS) -c calc_phi_gradient_ENO2.F90 

calc_h.o : calc_h.F90 global_app.mod global_mpi.mod 
	$(F90COMP) $(FFLAGS) -c calc_h.F90 

calc_h_gradient_ENO2.o : calc_h_gradient_ENO2.F90 global_app.mod global_mpi.mod 
	$(F90COMP) $(FFLAGS) -c calc_h_gradient_ENO2.F90 

des_eddy.o : des_eddy.F90 $(PACKMODULE2) checksum.mod
	$(F90COMP) $(FFLAGS) -c des_eddy.F90

rsmooth.o : rsmooth.F90 $(PACKMODULE2) checksum.mod
	$(F90COMP) $(FFLAGS) -c rsmooth.F90

sgtsv.o : sgtsv.f
	$(F90COMP) $(FFLAGS) -c sgtsv.f

xerbla.o : xerbla.f
	$(F90COMP) $(FFLAGS) -c xerbla.f

init_LSM.o : init_LSM.F90 $(PACKMODULE3) checksum.mod global_lsm.mod global_debug.mod
	$(F90COMP) $(FFLAGS) -c init_LSM.F90

levelsetmethod.o : levelsetmethod.F90 $(PACKMODULE4) $(PACKMODULE6) $(PACKMODULE2) checksum.mod global_debug.mod global_lsm.mod
	$(F90COMP) $(FFLAGS) -c levelsetmethod.F90

init_obstacle.o : init_obstacle.F90 $(PACKMODULE3) checksum.mod global_obstacle.mod global_debug.mod
	$(F90COMP) $(FFLAGS) -c init_obstacle.F90

bcond_obstacle.o : bcond_obstacle.F90 global_app.mod global_mpi.mod global_obstacle.mod
	$(F90COMP) $(FFLAGS) -c bcond_obstacle.F90

phi_obstacle.o : phi_obstacle.F90 global_app.mod global_mpi.mod global_obstacle.mod
	$(F90COMP) $(FFLAGS) -c phi_obstacle.F90

conservation.o : conservation.F90  $(PACKMODULE2) global_lsm.F90 
	$(F90COMP) $(FFLAGS) -c conservation.F90

# ================== CLEAN ==================

.PHONY : clean

clean : 
	rm -f *.mod
	rm -f *.o
	rm -f channel
	rm -f $(OBJS)

# ================== CONFIG SHORTCUT TARGETS ==================

fondecyt:
	$(MAKE) CONFIG=fondecyt

# This target runs a bash script to load necessary modules, configure the Spack environment, 
# load OpenBLAS, set library paths, and then run 'make' with the 'leftraru' configuration.
# The '@bash -c' command executes the following script as a new bash shell. The '-c' flag allows 
# us to pass the script as a string argument, executing it directly within the shell.
leftraru:
	@bash -c ' \
		echo "Loading necessary modules..." && \
		source /etc/profile.d/modules.sh && \
		ml purge && \
		ml iimpi/2018.04 && \
		echo "Loading Spack environment..." && \
		. ~/spack/share/spack/setup-env.sh && \
		echo "Ensuring Spack recognizes Intel compiler..." && \
		spack compiler find && \
		echo "Loading OpenBLAS..." && \
		spack load openblas%intel@18.0.5 && \
		echo "Setting LD_LIBRARY_PATH..." && \
		export LD_LIBRARY_PATH=$$(spack location -i openblas%intel@18.0.5)/lib:$$LD_LIBRARY_PATH && \
		echo "Running make with CONFIG=leftraru..." && \
		$(MAKE) CONFIG=leftraru \
	'

	