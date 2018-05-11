SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

#------------------------------------------------------------
#
# When you adapt this makefile to compile your CPLEX programs
# please copy this makefile and set CPLEXDIR and CONCERTDIR to
# the directories where CPLEX and CONCERT are installed.
#
#------------------------------------------------------------

CPLEXDIR      = /home/rottner/Programmes/Cplex12.8/cplex
CONCERTDIR    = /home/rottner/Programmes/Cplex12.8/concert

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

clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm -rf $(C_EX) $(CX_EX) $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
#
# The examples
#


mf: main.o
	$(CCC) $(CCFLAGS) $(CCLNDIRS) -o mf main.o $(CCLNFLAGS)
main.o: src/main.cpp
	$(CCC) -c $(CCFLAGS) src/main.cpp -o main.o


# Local Variables:
# mode: makefile
# End:
