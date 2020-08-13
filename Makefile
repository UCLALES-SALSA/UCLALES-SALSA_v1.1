###############################################################
#
# Location of code (in $ROOT) and location where model is to be built $BIN
#
ROOT      :=$(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
BIN       = $(ROOT)/bin
ARCH      := $(shell uname)
#
# Generic Variables
#
SRC         =$(ROOT)/src
SRC_UTIL    =$(SRC)/src_util
SRC_LES     =$(SRC)/src_LES
SRC_SALSA   =$(SRC)/src_salsa
SRC_RAD     =$(SRC)/src_rad
SRC_EMIS    =$(SRC)/src_emission
SRC_SHARED  =$(SRC)/src_shared
SRC_IO      =$(SRC)/src_IO
VPATH = $(SRC_LES):$(SRC_SALSA):$(SRC_UTIL):$(SRC_RAD):   \
        $(SRC_EMIS):$(SRC_SHARED):$(SRC_IO):$(SRC)
ECHO    = /bin/echo
RM      = /bin/rm -f
NETCDROOT   = $(NETCDF_FORTRAN_INSTALL_ROOT)
NETCDF_LIB  = -L$(NETCDFROOT)lib -lnetcdff
NETCDF_INCLUDE  = -I$(NETCDFROOT)/include
ARCHIVE = ar rs
RANLIB =:
SEQFFLAGS = -I$(SRC) $(NETCDF_INCLUDE)
MPIFFLAGS = -I$(SRC) $(NETCDF_INCLUDE)
MPI = /usr
MPILIB =
MPIINC =
NCDF = /usr
NCDFLIB =
NCDFINC =
LIBS = '$(NETCDF_LIB)'
F90 = ftn
MPIF90 = mpif90
#FFLAGS = -O2 -s real64 -Onoomp -h flex_mp=strict # -R bcps   # CRAY
FFLAGS = -O2 -march=native -real-size 64 -fp-model precise -convert big_endian -fpe0 # Puhti Intel                 # INTEL
#FFLAGS = -O2 -fdefault-real-8 -std=f2008 -fbounds-check # GCA
F77FLAGS = $(FFLAGS)
LES_OUT_MPI=$(BIN)/les.mpi
LES_OUT_SEQ=$(BIN)/les.seq
default: mpi
all:  mpi seq
seq: $(LES_OUT_SEQ)
mpi: $(LES_OUT_MPI)
$(LES_OUT_SEQ): 
	cd $(SRC); $(MAKE) LES_ARC=seq                      \
	FFLAGS='$(FFLAGS) $(SEQFFLAGS)' F90=$(F90)          \
	OUT=$(LES_OUT_SEQ) LIBS=$(LIBS) SRCUTIL=$(SRC_UTIL) \
	SRCLES=$(SRC_LES) SRCSALSA=$(SRC_SALSA)             \
	SRCRAD=$(SRC_RAD) SRCEMIS=$(SRC_EMIS)               \
	SRCSHARED=$(SRC_SHARED) SRCIO=$(SRC_IO)
$(LES_OUT_MPI):
	cd $(SRC); $(MAKE) LES_ARC=mpi                      \
	FFLAGS='$(FFLAGS) $(MPIFFLAGS)' F90=$(MPIF90)       \
	OUT=$(LES_OUT_MPI) LIBS=$(LIBS) SRCUTIL=$(SRC_UTIL) \
	SRCLES=$(SRC_LES) SRCSALSA=$(SRC_SALSA)             \
	SRCRAD=$(SRC_RAD) SRCEMIS=$(SRC_EMIS)               \
	SRCSHARED=$(SRC_SHARED) SRCIO=$(SRC_IO)
.PHONY: $(LES_OUT_SEQ) 
.PHONY: $(LES_OUT_MPI)
#
# cleaning
# --------------------
#
clean: cleanmpi cleanseq 
	$(RM) $(SRC)/*mod $(SRC)/*.o
cleanmpi:
	$(ECHO) "cleaning mpi model"
	$(RM) core $(LES_OUT_MPI) $(SRC)/mpi/*mod $(LES_ARC_MPI)
cleanseq:
	$(ECHO) "clean sequential model"
	$(RM) core $(LES_OUT_SEQ) $(SRC)/seq/*mod $(LES_ARC_SEQ)
FORCE: 
.PRECIOUS: $(LIBS)

