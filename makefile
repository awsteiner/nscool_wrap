# LIBS is the list of libraries
# LFC is the local fortran compiler
# LCXX is the local C++ compiler
# LCFLAGS are the local C++ compiler flags
# LFFLAGS are the local fortran compiler flags

#----------------------------------------------------------------------
# Organize data files by date

MONTH = $(shell date +'%m')
DAY = $(shell date +'%d')
YEAR = $(shell date +'%y')

# ----------------------------------------------------------------
# Various user-specific settings
# ----------------------------------------------------------------

ifeq ($(MACHINE),bridges)

# For bridges
LIB_DIRS = -L/sw/cs400_centos7.3\
_acfsoftware/gsl/2.3/centos7.3_intel17.2.174/lib \
-L/nics/d/home/asteine1/install/o2scl-0.922/lib
LIBS = -fPIC -L/pylon5/ph561dp/asteiner/install/o2scl-0.922/lib \
	-L/pylon5/ph561dp/asteiner/install/gsl-2.4/lib -lgfortran \
	-lo2scl_eos -lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 \
	-lgsl -lgslcblas -lm -lreadline
LFC = mpif90
LCXX = mpiicpc
LCFLAGS = -std=gnu++11 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
	-O3 -I/pylon5/ph561dp/asteiner/install/o2scl-0.922/include \
	-I/pylon5/ph561dp/asteiner/install/gsl-2.4/include -DO2SCL_OPENMP \
	-fopenmp -DO2SCL_MPI -DO2SCL_HDF_SVAR -DO2SCL_OLDER_COMPILER
LFFLAGS = -O3
LUSTRE_COPY = cp -r bhso qlmxb *.o2 $(SCRATCH); \
	mkdir -p $(SCRATCH)/script; cp script/* $(SCRATCH)/script
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

else
ifeq ($(MACHINE),acf)

ifeq ($(USER),sbeloin)

# For Spencer's login on beacon
LIBS = -fPIC -L/nics/d/home/asteine1/install/o2scl-0.922/lib \
	$(HDF5PARALLEL_LIB) $(GSL_LIB) -lo2scl_eos \
	-lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 \
	-lgsl -lgslcblas -lm -lreadline
LFC = mpif90
LCXX = mpiicpc
LCFLAGS = -std=gnu++11 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
	-O3 -I/nics/d/home/asteine1/install/o2scl-0.922/include \
	$(GSL_INC) $(HDF5PARALLEL_INC) -DO2SCL_OPENMP -fopenmp \
	-DO2SCL_HDF_SVAR -DO2SCL_OLDER_COMPILER
LFFLAGS = -O3
LUSTRE_COPY = cp -r bhso qlmxb *.o2 $(SCRATCHDIR); mkdir -p \
	$(SCRATCHDIR)/script; cp script/* $(SCRATCHDIR)/script
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

else

# For Andrew's login on beacon
LIBS = -fPIC -L$(O2SCL_LIB) $(HDF5PARALLEL_LIB) $(GSL_LIB) -lo2scl_eos \
	-lo2scl_part -lo2scl_hdf -lo2scl -lhdf5 \
	-lgsl -lgslcblas -lm -lreadline
LFC = mpif90
LCXX = mpiicpc
LCFLAGS = -std=gnu++11 -DGSL_RANGE_CHECK=0 -DO2SCL_NO_RANGE_CHECK \
	-O3 -I$(O2SCL_INC) $(GSL_INC) $(HDF5PARALLEL_INC) \
	-DO2SCL_OPENMP -fopenmp \
	-DO2SCL_HDF_SVAR -DO2SCL_OLDER_COMPILER
LFFLAGS = -O3
LUSTRE_COPY = cp -r bhso qlmxb *.o2 $(SCRATCHDIR); mkdir -p \
	$(SCRATCHDIR)/script; cp script/* $(SCRATCHDIR)/script
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

endif

else
ifeq ($(USER),hanjun)

# For Sophia's laptop
LIBS = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl
LFC = gfortran
LCXX = g++ -L/usr/local/gfortran/lib -DNO_MPI
LCFLAGS = -O3 -std=c++11
LFFLAGS = -O3 
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

else
ifeq ($(USER),spencerbeloin)

# For Spencer's laptop
LIBS = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl \
	-L/usr/local/lib -L/usr/local/gfortran/lib
LFC = gfortran
LCXX = mpic++ 
LCFLAGS = -O3 -std=c++11 -I/usr/local/include
LFFLAGS = -O3 
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

else
ifeq ($(HOSTNAME),isospin.roam.utk.edu)

ifeq ($(USER),sbeloin)

# On isospin for Spencer
LIBS = -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
	-lo2scl_hdf -lo2scl_eos \
	-lo2scl_part -lo2scl -lhdf5 -lgsl -lreadline -lgfortran
LFC = mpif90
LCXX = mpic++ 
LCFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/serial/include \
	-I/usr/include/eigen3 -Wno-deprecated-declarations \
	-O3 -std=c++11
LFFLAGS = -O3
STACK_CMD = ulimit -s 65536
OUT_DIR = data

else

# On isospin for Andrew
LIBS = -L$(HDF5_LIB) -lo2scl_hdf -lo2scl_eos \
	-lo2scl_part -lo2scl -lhdf5 -lgsl -lreadline -lgfortran
LFC = mpif90
LCXX = mpic++
LCFLAGS = -I$(HDF5_INC) \
	-I/usr/include/eigen3 -Wno-deprecated-declarations \
	-O0 -std=c++11 -DO2SCL_MPI -DO2SCL_OPENMP -fopenmp 
LFFLAGS = -O3 -fopenmp 
STACK_CMD = ulimit -s 65536 
OUT_DIR = $(HOME)/data/$(YEAR)/$(MONTH)/$(DAY)/

endif

else
ifeq ($(MACHINE),mimosa)

# On antares
LFC = gfortran -DNO_MPI
LCXX = mpic++
LIBS = -L$(O2SCL_LIB) -L$(GSL_LIB) -L$(HDF5_LIB) -L/usr/local/gfortran/lib \
	-lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl
LCFLAGS = -O3 -std=c++11 -I$(O2SCL_INC) -I$(EIGEN_INC) -I$(GSL_INC) \
	-Wno-deprecated-declarations -I$(HDF5_INC)
LFFLAGS = -O3
OUT_DIR = data

else
ifeq ($(MACHINE),hedgehog)

# On Andrew's laptop
LFC = mpif90
LCXX = mpic++
LIBS = -L$(O2SCL_LIB) -L$(GSL_LIB) -L$(HDF5_LIB) \
	-lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl \
	-lreadline -lgfortran 
LCFLAGS = -O3 -std=c++11 -I$(O2SCL_INC) -I$(EIGEN_INC) -I$(GSL_INC) \
	-Wall -Wno-deprecated-declarations -I$(HDF5_INC) \
	-Wno-ignored-attributes -Wno-unused \
	-DO2SCL_MPI -DO2SCL_OPENMP -fopenmp 
LFFLAGS = -O3 
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

else

# Default settings
LFC = $(FC)
LCXX = $(CXX)
LIBS = -lo2scl_hdf -lo2scl_eos -lo2scl_part -lo2scl -lhdf5 -lgsl
LCFLAGS = -O3
LFFLAGS = -O3 
STACK_CMD = echo "No stack command necessary"
OUT_DIR = data

endif
endif
endif
endif
endif
endif
endif

# ----------------------------------------------------------------
# Main targets
# ----------------------------------------------------------------

OBJS  = precool.o conductivity.o conductivity_core.o \
           opacity.o neutrino.o neutrino_core.o neutrino_crust.o \
           spec_heat.o density.o tc.o tc_Ioffe.o \
           Tools.o conductivity_crust.o 

OBJS_PG = precool_pg.o conductivity_pg.o conductivity_core_pg.o \
           opacity_pg.o neutrino_pg.o neutrino_core_pg.o neutrino_crust_pg.o \
           spec_heat_pg.o density_pg.o tc_pg.o tc_Ioffe_pg.o \
           Tools_pg.o conductivity_crust_pg.o 

help:
	@echo "test"

empty:

doc: empty
	cd sphinx/static; cat bib_header.txt > ../bib.rst
	cd sphinx/static; btmanip -parse nscool.bib -rst ../bib_temp.rst
	cd sphinx; cat bib_temp.rst >> bib.rst; rm -f bib_temp.rst
	cd doc; doxygen doxyfile
	cd sphinx; make html

test: $(OBJS) NSCool.o test.o nscool_wrap.h
	$(LCXX) $(LCFLAGS) -o test test.o NSCool.o \
		$(OBJS) $(LIBS)

test.o: test.cpp
	$(LCXX) $(LCFLAGS) -o test.o -c test.cpp 

test_pg: $(OBJS_PG) NSCool_pg.o test_pg.o nscool_wrap.h
	$(LCXX) -g -pg -o test_pg test_pg.o NSCool_pg.o $(OBJS_PG) \
		$(LIBS) -lgfortran

test_pg.o: test.cpp
	$(LCXX) $(LCFLAGS) -o test_pg.o -g -pg -c test.cpp 

aws_x1:
	bhso -reprocess ~/data/18/06/20/nc3_combined_0_out \
		~/data/18/06/20/nc3_reprocess_0_out
	acol -read ~/data/18/06/20/nc3_reprocess_0_out \
		-select-rows "log_wgt>(-830)" -internal \
		~/data/18/06/20/nc3_reprocess_0_out 

ot1h:
	mpirun -np 1 bhso \
		-threads 1 -initial-point-best last \
		-set file_update_time 1800 \
		-set include_cooling 0 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix newmcmc -mcmc

let:
	$(STACK_CMD) && bhso -threads 1 -load-estimate last_0_out \
		-initial-point-best last \
		-set file_update_time 1800 \
		-set include_cooling 1 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix newmcmc -mcmc > let.out 2> let.err &

nc1:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 1 \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix newmcmc -mcmc > nc1.scr 2> nc1.err &

nc2:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 2 -load-estimate data/last_nc_0_out \
		-initial-point-best data/last_nc \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/nc2 -mcmc

nc2_old:
	$(STACK_CMD) && mpirun -np 1 bhso_old \
		-threads 1 -load-estimate data/last_nc_0_out \
		-initial-point-best data/last_nc \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/nc2 -mcmc

#> data/nc2.scr 2> data/nc2.err &

nc3:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 6 -load-estimate data/last_nc_0_out \
		-initial-point-best data/last_nc \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/nc3 -mcmc > data/nc3.scr 2> data/nc3.err &

nc4:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 1 -read-prev-results data/last \
		-initial-point-best data/last \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/nc4 -mcmc
#> data/nc4.scr 2> data/nc4.err &

nc5:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 1 -read-prev-results data/last \
		-initial-point-best data/last \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/nc5 -mcmc \
		> $(OUT_DIR)/nc5.scr 2> $(OUT_DIR)/nc5.err &

nc5b:
	$(STACK_CMD) && mpirun -np 4 bhso \
		-threads 1 -read-prev-results data/last \
		-initial-point-best data/last \
		-set file_update_time 300 \
		-set verbose 2 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/nc5 -mcmc \
		> $(OUT_DIR)/nc5.scr 2> $(OUT_DIR)/nc5.err &

mg_2:
	bhso -make-gaussian2 data/temp mg_2

nc6:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results data/last \
		-initial-point-best data/last \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/nc6 -mcmc

ce1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results data/last \
		-set ptype mixed \
		-initial-point-best data/best \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/ce1 -mcmc

ce4:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results data/last \
		-set ptype mixed \
		-initial-point-best data/best \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/ce4 -mcmc \
		> $(OUT_DIR)/ce4.scr 2> $(OUT_DIR)/ce4.err &

mixed1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results data/last_0_out \
		-set refine_estimates 1 \
		-set ptype mixed \
		-initial-point-best data/best_0_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/mixed1 -mcmc

mixed1b:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set ptype mixed \
		-initial-point-best data/best \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mixed1b -mcmc \
		> $(OUT_DIR)/mixed1b.scr 2> $(OUT_DIR)/mixed1b.err &

mixed1c:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set ptype mixed \
		-initial-point-best data/best \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mixed1c -mcmc \
		> $(OUT_DIR)/mixed1c.scr 2> $(OUT_DIR)/mixed1c.err &

mixed1d:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set ptype mixed \
		-initial-point-best data/best \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mixed1d -mcmc \
		> $(OUT_DIR)/mixed1d.scr 2> $(OUT_DIR)/mixed1d.err &

cct:
	-$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 0 \
		-initial-point-best data/last_out \
		-set step_fac 3.0e1 \
		-set file_update_time 300 \
		-set verbose 2 -set mcmc_verbose 1 \
		-mcmc > temp_old
	-$(STACK_CMD) && bhso -threads 1 \
		-set compute_estimates 0 \
		-initial-point-best data/last_out \
		-set step_fac 3.0e1 \
		-set file_update_time 300 \
		-set verbose 2 -set mcmc_verbose 1 \
		-set new_cc 1 -mcmc > temp_new

testy:
	-$(STACK_CMD) && bhso -threads 1 \
		-initial-point-best data/best_one_out \
		-set verbose 2 -set new_cc 1 -mcmc > data/testx_new.txt

testx:
	-$(STACK_CMD) && bhso -threads 1 \
		-initial-point-best data/best_one_out \
		-set verbose 2 -mcmc > data/testx_old.txt
	-$(STACK_CMD) && bhso -threads 1 \
		-initial-point-best data/best_one_out \
		-set verbose 2 -set new_cc 1 -mcmc > data/testx_new.txt
	tail -n 40 data/testx_old.txt
	tail -n 40 data/testx_new.txt

point1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-initial-point-last smallrad.o2 \
		-set verbose 2 -point > data/point1.scr 2> data/point1.err &

step1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 0 \
		-initial-point-best data/best_one_out \
		-set step_fac 3.0e1 \
		-set file_update_time 300 \
		-set verbose 2 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/step1 -mcmc

gauss1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set ptype gauss -initial-point-best \
		"/home/awsteiner/data/18/newmcmc/all_limitK_out" \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix gauss1 -mcmc

gauss2:
	$(STACK_CMD) && mpirun -np 2 bhso -threads 2 \
		-set ptype gauss -initial-point-best \
		"/home/awsteiner/data/18/09/28/temp" \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/gauss2 -mcmc \
		> $(OUT_DIR)/gauss2.scr 2> $(OUT_DIR)/gauss2.err &

gauss4:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set covar_dec_factor 10.0 \
		-set include_vela 1 \
		-set ptype gauss \
		-initial-point-best \
		"/home/awsteiner/data/18/newmcmc2/w_vela_out" \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/gauss4 -mcmc \
		> $(OUT_DIR)/gauss4.scr 2> $(OUT_DIR)/gauss4.err &

gauss4em:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results \
		/home/awsteiner/data/18/newmcmc2/w_vela_out \
		-set covar_dec_factor 5.0 \
		-set ptype gauss -set include_vela 1 \
		-initial-point-best \
		/home/awsteiner/data/18/newmcmc2/w_vela_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/gauss4em -mcmc \
		> $(OUT_DIR)/gauss4em.scr 2> $(OUT_DIR)/gauss4em.err &

gauss2em:
	$(STACK_CMD) && mpirun -np 2 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results \
		/home/awsteiner/data/18/newmcmc2/wo_vela_out \
		-set covar_dec_factor 5.0 \
		-set ptype gauss -set include_vela 0 \
		-initial-point-best \
		/home/awsteiner/data/18/newmcmc2/wo_vela_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/gauss2em -mcmc \
		> $(OUT_DIR)/gauss2em.scr 2> $(OUT_DIR)/gauss2em.err &

#TEMPFILE = $(HOME)/wo_vela_out_nofail
TEMPFILE = $(HOME)/data/18/newmcmc2/w_vela_out_nofail

gauss1em:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results $(TEMPFILE) \
		-set covar_dec_factor 100.0 \
		-set ptype gauss -set include_vela 0 \
		-initial-point-best $(TEMPFILE) \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/gauss1em -mcmc \
		> gauss1em.scr 2> gauss1em.err &

gauss4em2:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results $(TEMPFILE) \
		-set covar_dec_factor 100.0 \
		-set ptype gauss -set include_vela 1 \
		-initial-point-best $(TEMPFILE) \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/gauss4em2 -mcmc \
		> $(OUT_DIR)/gauss4em2.scr 2> $(OUT_DIR)/gauss4em2.err &

gauss3:
	$(STACK_CMD) && mpirun -np 3 bhso -threads 1 \
		-set ptype gauss -initial-point-best "data/last_out" \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/gauss3 -mcmc \
		> $(OUT_DIR)/gauss3.scr 2> $(OUT_DIR)/gauss3.err &

step4:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results data/last_out \
		-initial-point-best data/last_out \
		-set step_fac 3.0e1 \
		-set file_update_time 1000 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/step4 -mcmc \
		> $(OUT_DIR)/step4.scr 2> $(OUT_DIR)/step4.err &

step3:
	$(STACK_CMD) && mpirun -np 3 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results data/all_out \
		-initial-point-best data/all_out \
		-set ptype gauss \
		-set step_fac 3.0e1 \
		-set file_update_time 1000 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/step3 -mcmc \
		> $(OUT_DIR)/step3.scr 2> $(OUT_DIR)/step3.err &

pdma1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set compute_estimates 0 \
		-set refine_estimates 0 \
		-read-prev-results data/last_out \
		-set ptype pdma \
		-initial-point-best data/last_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set new_cc 1 \
		-set prefix $(OUT_DIR)/pdma1 -mcmc

#	bhso -make-gaussian \
#/home/awsteiner/data/18/10/15/w_vela_out gaussian

mixed4:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set ptype mixed -set include_vela 0 \
		-set covar_dec_factor 30.0 \
		-initial-point-best \
		/home/awsteiner/data/18/newmcmc2/wo_vela_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mixed4 -mcmc \
		> $(OUT_DIR)/mixed4.scr 2> $(OUT_DIR)/mixed4.err &

mixed4em_gaussian:
	bhso -make-gaussian \
		/home/awsteiner/data/18/newmcmc2/wo_vela_out gaussian

mixed4em:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set compute_estimates 1 \
		-read-prev-results \
		/home/awsteiner/data/18/newmcmc2/wo_vela_out \
		-set covar_dec_factor 30.0 \
		-set ptype mixed -set include_vela 0 \
		-initial-point-best \
		/home/awsteiner/data/18/newmcmc2/wo_vela_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mixed4em -mcmc \
		> $(OUT_DIR)/mixed4em.scr 2> $(OUT_DIR)/mixed4em.err &

pdma2:
	$(STACK_CMD) && mpirun -np 2 bhso -threads 1 \
		-set compute_estimates 0 \
		-set refine_estimates 0 \
		-read-prev-results data/last_out \
		-set ptype pdma \
		-initial-point-best data/last_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set new_cc 1 \
		-set prefix $(OUT_DIR)/pdma2 -mcmc \
		> $(OUT_DIR)/pdma2.scr 2> $(OUT_DIR)/pdma2.err &

pdma4:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set compute_estimates 0 \
		-set refine_estimates 0 \
		-read-prev-results data/last_out \
		-set ptype pdma \
		-initial-point-best data/last_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set new_cc 1 \
		-set prefix $(OUT_DIR)/pdma4 -mcmc \
		> $(OUT_DIR)/pdma4.scr 2> $(OUT_DIR)/pdma4.err &

mixed2:
	$(STACK_CMD) && mpirun -np 2 bhso -threads 1 \
		-set ptype mixed -set include_vela 0 \
		-initial-point-best \
		/home/awsteiner/data/18/10/15/wo_vela_out \
		-set file_update_time 300 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mixed2 -mcmc \
		> $(OUT_DIR)/mixed2.scr 2> $(OUT_DIR)/mixed2.err &

mv1:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 1 -read-prev-results data/exact_0_out \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix data/mv1 -max-var 

mv3:
	$(STACK_CMD) && mpirun -np 3 bhso \
		-threads 1 -read-prev-results data/exact_0_out \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mv3 -max-var \
		> $(OUT_DIR)/mv3.scr 2> $(OUT_DIR)/mv3.err &

mv2:
	$(STACK_CMD) && mpirun -np 2 bhso \
		-threads 1 -read-prev-results data/exact_0_out \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mv2 -max-var \
		> $(OUT_DIR)/mv2.scr 2> $(OUT_DIR)/mv2.err &

mvo:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 4 -read-prev-results data/last \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mvo -max-var \
		> $(OUT_DIR)/mvo.scr 2> $(OUT_DIR)/mvo.err &

mvm:
	$(STACK_CMD) && mpirun -np 4 bhso \
		-threads 1 -read-prev-results data/last \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mvm -max-var \
		> $(OUT_DIR)/mvm.scr 2> $(OUT_DIR)/mvm.err &

sun1:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 1 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/sun1 -sunday \
		> $(OUT_DIR)/sun1.scr 2> $(OUT_DIR)/sun1.err &

sunt:
	$(STACK_CMD) && mpirun -np 1 bhso -threads 2 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/sun1 -sunday 
#		> $(OUT_DIR)/sun1.scr 2> $(OUT_DIR)/sun1.err &

sun4:
	$(STACK_CMD) && mpirun -np 4 bhso -threads 1 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/sun4 -sunday \
		> $(OUT_DIR)/sun4.scr 2> $(OUT_DIR)/sun4.err &

mm1:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 1 -initial-point-best data/best \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mmin -mmin2 \
		> $(OUT_DIR)/mm1.scr 2> $(OUT_DIR)/mm1.err &

mm1b:
	$(STACK_CMD) && mpirun -np 1 bhso \
		-threads 1 -initial-point-best data/best \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mm1b -mmin2 \
		> $(OUT_DIR)/mm1b.scr 2> $(OUT_DIR)/mm1b.err &

mm4:
	$(STACK_CMD) && mpirun -np 4 bhso \
		-threads 1 -initial-point-best data/best_0_out \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mm4 -mmin \
		> $(OUT_DIR)/mm4.scr 2> $(OUT_DIR)/mm4.err &

mm2:
	$(STACK_CMD) && mpirun -np 2 bhso \
		-threads 1 -initial-point-best data/last_out \
		-set new_cc 1 \
		-set verbose 1 -set mcmc_verbose 1 \
		-set prefix $(OUT_DIR)/mm2 -mmin \
		> $(OUT_DIR)/mm2.scr 2> $(OUT_DIR)/mm2.err &

paramsurv:
	ulimit -s 65536 && ./newmcmc -threads 1 \
		-set include_cooling 1 -param-survey > ps.out 2> ps.err &

psaws:
	ulimit -s 65536 && ./bhso -threads 1 -set verbose 1 \
		-set include_cooling 1 -param-survey > ps.out 2> ps.err & 

mat_test:
	o2graph -set logy 1 -read mg_2 covar -diag -plot1 color=black \
		-read mg_3 covar -diag -plot1 color=blue \
		-save var_plot.pdf -show

mass:
	mpirun -np 1 bhso \
		-initial-point-best last \
		-set verbose 1 -set mcmc_verbose 1 \
		-mass-survey > mass_survey.scr 2> mass_survey.err &

maxvar:
	rm -f newmcmc.scr newmcmc.err newmcmc_?_scr newmcmc_?_out
	mpirun -np 1 bhso \
		-initial-point-best last \
		-set verbose 1 -set mcmc_verbose 1 \
		-max-var

bhso: $(OBJS) NSCool.o bhso.o nscool_wrap.h 
	$(LCXX) $(LCFLAGS) -o bhso bhso.o NSCool.o \
		$(OBJS) $(LIBS)

bhso.o: bhso.cpp 
	$(LCXX) $(LCFLAGS) -o bhso.o -c bhso.cpp 

clean:
	-rm -f *.o test test_pg bhso

# ----------------------------------------------------------------
# ACF targets
# ----------------------------------------------------------------

lustre_copy:
	$(LUSTRE_COPY)

# ----------------------------------------------------------------
# Object files
# ----------------------------------------------------------------

NSCool.o: NSCool.f nscool_wrap.h
	$(LFC) $(LFFLAGS) -c NSCool.f
precool.o: precool.f
	$(LFC) $(LFFLAGS) -c precool.f
conductivity.o: conductivity.f
	$(LFC) $(LFFLAGS) -c conductivity.f
conductivity_core.o: conductivity_core.f
	$(LFC) $(LFFLAGS) -c conductivity_core.f
conductivity_crust.o: conductivity_crust.f
	$(LFC) $(LFFLAGS) -c conductivity_crust.f
opacity.o: opacity.f
	$(LFC) $(LFFLAGS) -c opacity.f
neutrino.o: neutrino.f
	$(LFC) $(LFFLAGS) -c neutrino.f
neutrino_core.o: neutrino_core.f
	$(LFC) $(LFFLAGS) -c neutrino_core.f
neutrino_crust.o: neutrino_crust.f
	$(LFC) $(LFFLAGS) -c neutrino_crust.f
spec_heat.o: spec_heat.f
	$(LFC) $(LFFLAGS) -c spec_heat.f
density.o: density.f
	$(LFC) $(LFFLAGS) -c density.f
tc.o: tc.f
	$(LFC) $(LFFLAGS) -c tc.f
tc_Ioffe.o: tc_Ioffe.f
	$(LFC) $(LFFLAGS) -c tc_Ioffe.f
Tools.o: Tools.f
	$(LFC) $(LFFLAGS) -c Tools.f


# ----------------------------------------------------------------
# Object files
# ----------------------------------------------------------------

NSCool_pg.o: NSCool.f nscool_wrap.h
	$(LFC) $(LFFLAGS) -g -pg -c NSCool.f -o NSCool_pg.o
precool_pg.o: precool.f
	$(LFC) $(LFFLAGS) -g -pg -c precool.f -o precool_pg.o
conductivity_pg.o: conductivity.f
	$(LFC) $(LFFLAGS) -g -pg -c conductivity.f -o conductivity_pg.o
conductivity_core_pg.o: conductivity_core.f
	$(LFC) $(LFFLAGS) -g -pg -c conductivity_core.f \
		-o conductivity_core_pg.o
conductivity_crust_pg.o: conductivity_crust.f
	$(LFC) $(LFFLAGS) -g -pg -c conductivity_crust.f \
		-o conductivity_crust_pg.o
opacity_pg.o: opacity.f
	$(LFC) $(LFFLAGS) -g -pg -c opacity.f -o opacity_pg.o
neutrino_pg.o: neutrino.f
	$(LFC) $(LFFLAGS) -g -pg -c neutrino.f -o neutrino_pg.o
neutrino_core_pg.o: neutrino_core.f
	$(LFC) $(LFFLAGS) -g -pg -c neutrino_core.f -o neutrino_core_pg.o
neutrino_crust_pg.o: neutrino_crust.f
	$(LFC) $(LFFLAGS) -g -pg -c neutrino_crust.f -o neutrino_crust_pg.o
spec_heat_pg.o: spec_heat.f
	$(LFC) $(LFFLAGS) -g -pg -c spec_heat.f -o spec_heat_pg.o
density_pg.o: density.f
	$(LFC) $(LFFLAGS) -g -pg -c density.f -o density_pg.o
tc_pg.o: tc.f
	$(LFC) $(LFFLAGS) -g -pg -c tc.f -o tc_pg.o
tc_Ioffe_pg.o: tc_Ioffe.f
	$(LFC) $(LFFLAGS) -g -pg -c tc_Ioffe.f -o tc_Ioffe_pg.o
Tools_pg.o: Tools.f
	$(LFC) $(LFFLAGS) -g -pg -c Tools.f -o Tools_pg.o
