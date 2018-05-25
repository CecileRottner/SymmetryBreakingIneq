SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

VPATH += src

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /home/rottner/Programmes/Cplex12.8/cplex
CONCERTDIR    = /home/rottner/Programmes/Cplex12.8/concert

include cplex_dir_128.mk

# ---------------------------------------------------------------------
# Compiler selection 
# ---------------------------------------------------------------------

CCC = g++ -O0

# ---------------------------------------------------------------------
# Compiler options 
# ---------------------------------------------------------------------

CCOPT = -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXJARDIR   = $(CPLEXDIR)/lib/cplex.jar
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNDIRS  = -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR)
CLNDIRS   = -L$(CPLEXLIBDIR)
CCLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -ldl
CLNFLAGS  = -lcplex -lm -lpthread -ldl


all:
	make all_cpp

execute: all
	make execute_cpp

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

EXDIR         = $(CPLEXDIR)/examples
EXINC         = $(EXDIR)/include
EXDATA        = $(EXDIR)/data
EXSRCCPP      = $(EXDIR)/src/cpp

CFLAGS  = $(COPT)  -I$(CPLEXINCDIR)
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
JCFLAGS = $(JOPT)


#------------------------------------------------------------
#  make all      : to compile the examples. 
#  make execute  : to compile and execute the examples.
#------------------------------------------------------------

CPP_EX = mf

all_cpp: $(CPP_EX)


execute_cpp: $(CPP_EX)
	 ./mf





# ------------------------------------------------------------
#
# Obj
#

OBJ_LIST := InstanceUCP.o \
	    ModeleUCP.o \
	    ModeleFlot.o \
            Process.o \
	    main.o \

OBJ:= $(addprefix build/, $(OBJ_LIST))

# ------------------------------------------------------------
.PHONY: clean mrproper

clean:
	@$(remove) $(OBJ)

mrproper: clean
	@$(remove) $(BIN)

# ------------------------------------------------------------
# link

mf: $(OBJ)
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o $@ $^ $(CCLNFLAGS)


# ------------------------------------------------------------
# compil

build/%.o: src/%.cc
	$(CCC) -c $(CCFLAGS) $^ -o $@


# List of dependencies
build/InstanceUCP.o   : InstanceUCP.h
build/ModeleFlot.o : InstanceUCP.h ModeleFlot.h
build/ModeleUCP.o : InstanceUCP.h ModeleUCP.h
build/Process.o : Process.h InstanceUCP.h ModeleUCP.h
build/main.o      : InstanceUCP.h ModeleUCP.h Process.h 

# Local Variables:
# mode: makefile
# End:
