##### User configurable options #####

MAKE          = make --no-print-directory

MPI_CC        = /usr/bin/mpicc
MPI_CLINKER   = $(MPI_CC)
#
DEBUG_LIBS    =
MPI_INC       =
MPI_LIBS      =
MPE_CFLAGS    =  -DMPI_LINUX -DUSE_STDARG -DHAVE_PROTOTYPES
MPE_FFLAGS    =
MPE_INCDIR    = /usr/lib/openmpi/include/mpe
MPE_LIBDIR    = /usr/lib/openmpi/lib/mpe
PROF_LIBS     = -L$(MPE_LIBDIR) -lmpe $(MPI_LIBS)
LOG_LIBS      = -L$(MPE_LIBDIR) -llmpe -lmpe $(MPI_LIBS)
TRACE_LIBS    = -L$(MPE_LIBDIR) -ltmpe $(MPI_LIBS)
ANIM_LIBS     = -L$(MPE_LIBDIR) -lampe -lmpe -L/usr/X11R6/lib -lX11 $(MPI_LIBS)
XGRAPH_LIBS   = -L$(MPE_LIBDIR) -lmpe -L/usr/X11R6/lib -lX11 $(MPI_LIBS)
F2CMPI_LIBS   = -lfmpich
FLIB_PATH     = -L$(MPE_LIBDIR) $(F2CMPI_LIBS) -llmpe -lmpe \
                $(MPI_LIBS)
#
INCLUDE_DIR   = -I$(MPE_INCDIR) $(MPI_INC)
LDFLAGS       =
NO_X          = 0
NOF77         = 0
RM            = /bin/rm

### End User configurable options ###

VPATH=.:$(srcdir)


CFLAGS        = $(MPE_CFLAGS) $(INCLUDE_DIR)
FFLAGS        = $(MPE_FFLAGS) $(INCLUDE_DIR)
CC_LDFLAGS    = $(LDFLAGS) $(MPE_CFLAGS)
FC_LDFLAGS    = $(LDFLAGS) $(MPE_FFLAGS)
EXECS         = BucketSort


.c.o:
	$(MPI_CC) $(CFLAGS) -c $<

default: $(EXECS)

BucketSort: BucketSort.o
	$(MPI_CLINKER) $(CC_LDFLAGS)  -o $@ $? $(XGRAPH_LIBS) $(DEBUG_LIBS) -lm

