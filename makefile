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

sxrt2: $(OBJS) NSCool.o sxrt2.o nscool_wrap.h 
	$(LCXX) $(LCFLAGS) -o sxrt2 sxrt2.o NSCool.o \
		$(OBJS) $(LIBS)

sxrt2.o: sxrt2.cpp 
	$(LCXX) $(LCFLAGS) -o sxrt2.o -c sxrt2.cpp 

clean:
	-rm -f *.o test sxrt2

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
