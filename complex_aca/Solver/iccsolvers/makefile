
ifndef compiler
compiler=intel
endif

ifndef interface
interface=em64t
endif

ifndef threading
threading=parallel
threading=sequential
endif

ifndef _IA
_IA=em64t
endif

ifndef EXT
EXT=so
endif

ICCROOT = /opt/intel/Compiler/11.1/064
MKLROOT = $(ICCROOT)/mkl
MKL_PATH = $(MKLROOT)/lib/$(_IA)

override CC = icc
CCOPTS = -w -I$(MKLROOT)/include
override FC = ifort
FCOPTS = -w -I$(MKLROOT)/include -fpp

ifeq ($(interface),ilp64)
IFACE_LIB=$(MKL_PATH)/libmkl_intel_ilp64.$(EXT)
FOPTS = -i8
SOLVER_INTERFACE=_ilp64
COPTS= -DMKL_ILP64
else
IFACE_LIB=$(MKL_PATH)/libmkl_intel_lp64.$(EXT)
FOPTS =
SOLVER_INTERFACE=_lp64
COPTS= -g
endif

ifeq ($(_IA),32)
   ifeq ($(compiler),intel)
       SPEC_OPT=-xK
#This option is required by Intel(R) 11.0 compiler to produce workable binaries for Pentium(R) III.
#If you don't need it, you can remove this option.
   endif
IFACE_LIB=$(MKL_PATH)/libmkl_intel.$(EXT)
FOPTS =
SOLVER_INTERFACE=
COPTS=
endif

ifeq ($(threading),sequential)
THREADING_LIB=$(MKL_PATH)/libmkl_sequential.$(EXT)
OMP_LIB =
SOLVER_THREADING=_sequential
else
THREADING_LIB=$(MKL_PATH)/libmkl_intel_thread.$(EXT)
OMP_LIB = -L$(MKL_PATH) -liomp5
OMP_LIB = -L$(ICCROOT)/lib/intel64 -liomp5
SOLVER_THREADING=
endif

CORE_LIB=$(MKL_PATH)/libmkl_core.$(EXT)

ifeq ($(EXT),so)
MKL_LAPACKLIB=$(MKL_PATH)/libmkl_lapack.$(EXT)
else
MKL_LAPACKLIB=
endif

SOLVER_LIB=$(MKL_PATH)/libmkl_solver$(SOLVER_INTERFACE)$(SOLVER_THREADING).a

MKL_LIBS= $(SOLVER_LIB) $(IFACE_LIB) -Wl,--start-group $(THREADING_LIB) $(MKL_LAPACKLIB) $(CORE_LIB) -Wl,--end-group $(OMP_LIB)


.c.o:
	$(CC) $(SPEC_OPT) $(CCOPTS) $(COPTS) -c $<  -o $@

test_iccgmres: test_iccgmres.o iccsolvers.o
	$(CC) \
		$(SPEC_OPT) \
		$(CCOPTS) \
		$(COPTS) \
		test_iccgmres.o iccsolvers.o \
		-L$(MKL_PATH) \
		$(MKL_LIBS) \
		-lpthread -lm -o test_iccgmres
	#export LD_LIBRARY_PATH=$(MKL_PATH):$(ICCROOT)/lib/intel64:$(LD_LIBRARY_PATH); test_iccgmres 

test_zgmres1: test_zgmres1.o zgmres1.o
	$(CC) \
		$(SPEC_OPT) \
		$(CCOPTS) \
		$(COPTS) \
		test_zgmres1.o zgmres1.o \
		-L$(MKL_PATH) \
		$(MKL_LIBS) \
		-lpthread -lm -o test_iccgmres
	#export LD_LIBRARY_PATH=$(MKL_PATH):$(ICCROOT)/lib/intel64:$(LD_LIBRARY_PATH); test_iccgmres 
clean:
	rm -rf test_iccgmres.o iccsolvers.o test_iccgmres


#-------------------------------------------------------------------------------
