#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  This is the generic part of the Makefile. Use this
#  template within each project.
#
#  General definitions
#
#  Paths
#
bindir = ./bin
obsdir = ./obj
moduledir = ./mod
srcdir = ./src
AR = ar
F95 = mpif90


ifeq ($(notdir $(F95)),g95)
        FFLAGS = -O3 -Wunused -fmod=$(moduledir)
else
	FFLAGS =   -J$(moduledir) -fimplicit-none -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check #-Wall
#	FFLAGS =   -O3 -J$(moduledir) -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check -Wunused
#        ifeq ($(notdir $(MPIF90)),gfortran)
#        FFLAGS = -O3 -J$(moduledir) -fbacktrace -fbounds-check -Wunused  -Wline-truncation -fimplicit-none 
#		FFLAGS =   -O3 -J$(moduledir) -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check
#		FFLAGS =   -O3 -J$(moduledir) -fimplicit-none -funroll-all-loops -ffixed-line-length-none -ffree-line-length-none -fbacktrace -fbounds-check
 #       else#  ($(notdir $(MPIF90)),ifort)
	#	FFLAGS =  -module $(moduledir) -heap-arrays -funroll-loops
     #   endif
endif

obj_prog = 	program_dg2d.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	triangulation.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	timeloopMod.o \
	matrixMod.o \
	waveMod.o \
	dtMod.o \
	stfMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o \
	externalModelMod.o \
	mpiMod.o\
	fileunitMod.o\
	fileParameterMod.o\
	logoMod.o\
	pmlMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o

obj_mesher = program_mesher.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
        vtkMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	triangulation.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	matrixMod.o \
	dtMod.o \
	stfMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o\
	externalModelMod.o \
	linearSystemMod.o \
	fileunitMod.o \
	logoMod.o\
	fileParameterMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o \
	progressbarMod.o

obj_solver = program_solver.o \
	mpiincludeMod.o \
	constantsMod.o \
	parameterMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
        vtkMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	triangulation.o \
	plotMod.o \
	liftMod.o \
	derMod.o \
	normalsMod.o \
	timeloopMod.o \
	matrixMod.o \
	waveMod.o \
	dtMod.o \
	stfMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o \
	externalModelMod.o \
	mpiMod.o \
	linearSystemMod.o \
	fileunitMod.o \
	fileParameterMod.o\
	logoMod.o\
	pmlMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o \
	progressbarMod.o

obj_movie = program_movie.o \
	constantsMod.o \
	parameterMod.o \
	collectMovieMod.o \
        vtkMod.o \
	gllMod.o \
	warpfactorMod.o \
	nodesMod.o \
	triTrafoMod.o \
	jacobiMod.o \
	simplexMod.o \
	vandermondeMod.o \
	dmatricesMod.o \
	geometricFactorsMod.o \
	rosettaGammaMod.o \
	meshMod.o \
	triangulation.o \
	plotMod.o \
	normalsMod.o \
	matrixMod.o \
	dtMod.o \
	liftMod.o \
	sourceReceiverMod.o \
	triangleNeighborMod.o \
	externalModelMod.o \
	mpiMod.o  \
	linearSystemMod.o \
	fileunitMod.o \
	logoMod.o\
	fileParameterMod.o \
	errorMessage.o \
	realloc.o \
	materialsMod.o \
	progressbarMod.o



#-------------------------------------------------------
#  Direcory search
#
vpath %.o $(obsdir)
#--------------------------------------------------------
#  additional directories to be searched for module or include dependencies
#  default is search in ./ only
#
DEPDIRS = 
#-------------------------------------------------------
#  Implicit rule to compile .o files from .f90 files.
#  Because of vpath, targets and dependencies need not be
#  in the current directory.
#
%.o: $(srcdir)/%.f90
	$(F95) -c $(FFLAGS) $< -o $(obsdir)/$@
#--------------------------------------------------------------
#  Object string for linking:
#  Adds object dir as prefix and removes directory part
#  of $^ (all dependencies)
#
obstring = $(addprefix $(obsdir)/,$(notdir $^))
obstringtest = $(addprefix $(obsdir)/,$(notdir $^))
#
#   End of generic part of Makefile
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#  Library paths
#
#pgplot = -lpgplot -L/usr/lib -lpng -lz -lX11
#pgplot =  -L/usr/X11R6/lib -lX11 -Wl,-framework -Wl,Foundation -L/sw/lib -lpng -lz -L/usr/local/lib -laquaterm  -L/sw/lib/pgplot -lpgplot
la = -llapack -lblas
#lib = /usr/lib/libmetis.a /usr/lib/libscotch.a /usr/lib/libscotcherr.a
lib = metis-4.0.3/libmetis.a 
#---------------------------------------------------------
.PHONY: all

#
#  create dependencies on modules 
#  make.incdep is a Makefile because it is included. Such files are first updated
#  before anything else and make starts from anew including all updated makefiles
#
make.incdep:
	./tools/scripts/makeDepFromUseInclude.py $(srcdir) $(DEPDIRS) > $@
-include make.incdep
#
#       Targets
#
all: required mesher solver movie

allclean :
	rm -f $(bindir)/mesher $(bindir)/movie $(moduledir)/*.mod $(obsdir)/*.o out/* make.incdep

clean :
	rm -f $(bindir)/solver $(bindir)/mesher $(bindir)/movie $(moduledir)/*.mod $(obsdir)/*.o make.incdep


dg2d: $(obj_prog)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la) $(lib)

mesher: $(obj_mesher)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

solver: $(obj_solver)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

movie: $(obj_movie)
	$(F95) $(FFLAGS) -o $(bindir)/$@ $(obstring) $(pgplot) $(la)  $(lib)

bin:
	mkdir -p $(bindir)

mod:
	mkdir -p $(moduledir)

obj:
	mkdir -p $(obsdir)

required: mod obj bin
