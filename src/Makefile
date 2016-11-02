EXECDIR      = ../run
OBJDIR       = ../obj
EFFECTIVE    = $(EXECDIR)/eff
INTERACTING  = $(EXECDIR)/inter
ANNEAL       = $(EXECDIR)/anneal
SINGLE       = $(EXECDIR)/singlenocell
WIDOM        = $(EXECDIR)/widom
UMBRELLA     = $(EXECDIR)/umbrella
UMBRELLABIAS = $(EXECDIR)/umbrellabias
INC          = $(MKLROOT)/include
LIB          = $(MKLROOT)/lib/intel64 
CC           = gcc
CLINKER      = gcc

#CFLAGS       = -std=c99 -Wall  -O3 -w1 -mp1 -g  #-static
CFLAGS       = -std=c99 -Wall  -O3  -g  #-static
MYFLAGS      = -DLJ=0 #-DLUT_FENE #-DLUT_MORSE 
LIBS         = -lm
DEPEND       = makedepend
## profiling: compile and link with prof-gen, run, c & l with prof-use
#PROFDIR      = ../prof
#PROFILE     = -prof-use -prof-dir$(PROFDIR) -ipo
#PROFILE     = -prof-genx -prof-dir$(PROFDIR)
## parallel report and auto
#PARALLEL = -parallel -par-report2

# INTEL MATH LIBRARY

MKLPATH=$(MKLROOT)/lib/ia32
MKLINCLUDE=$(MKLROOT)/include

# default src files
SRC = readInput.c system.c initialize.c dendrimers.c \
 ran250.c potential.c moves.c energy.c  mc.c cellLists.c
OBJS = $(patsubst %.c,$(OBJDIR)/%.o,$(SRC))
# single source files
SINGLE_SRC = singledend.c $(SRC)
SINGLE_OBJS =$(patsubst %.c,$(OBJDIR)/%.o,$(SINGLE_SRC))
# effective interaction fixed CMass
EFFECTIVE_SRC  = MCdendeff.c $(SRC)
EFFECTIVE_OBJS =$(patsubst %.c,$(OBJDIR)/%.o,$(EFFECTIVE_SRC))
# interacting dendrimers
INTERACTING_SRC  = MCinterdend.c $(SRC)
INTERACTING_OBJS =$(patsubst %.c,$(OBJDIR)/%.o,$(INTERACTING_SRC))
# simulated annealing objects
ANNEAL_SRC  = MCanneal.c $(SRC)
ANNEAL_OBJS =$(patsubst %.c,$(OBJDIR)/%.o,$(ANNEAL_SRC))
# widom insertion particle method
WIDOM_SRC  = widom.c $(SRC)
WIDOM_OBJS =$(patsubst %.c,$(OBJDIR)/%.o,$(WIDOM_SRC))
# umbrella sampling
UMBRELLA_SRC  = umbrella.c $(SRC)
UMBRELLA_OBJS =$(patsubst %.c,$(OBJDIR)/%.o,$(UMBRELLA_SRC))

default: umbrella

all:
widom:  CFLAGS+=  
widom:  MYFLAGS+= -DUSE_PBC=0 -DUSE_CELL_LIST=0 -DCELL_SEC_NEI=0 -DFIXED_CMASS=0
widom:  $(WIDOM_OBJS)
	$(CLINKER) $(MYFLAGS) $(CFLAGS)-o $(WIDOM) $(WIDOM_OBJS) $(LIBS)
umbrella:  CFLAGS+= 
umbrella:  MYFLAGS+=  -DUSE_PBC=0 -DUSE_CELL_LIST=1 -DCELL_SEC_NEI=0 -DFIXED_CMASS=0
umbrella:  $(UMBRELLA_OBJS)
	$(CLINKER) $(MYFLAGS) $(CFLAGS)-o $(UMBRELLA) $(UMBRELLA_OBJS) $(LIBS)
umbrellabias:  MYFLAGS+=  -DUSE_PBC=0 -DBIAS=1 -DUSE_CELL_LIST=1 -DCELL_SEC_NEI=0 -DFIXED_CMASS=0
umbrellabias:  $(UMBRELLA_OBJS)
	$(CLINKER) $(MYFLAGS) $(CFLAGS)-o $(UMBRELLABIAS) $(UMBRELLA_OBJS) $(LIBS)
single: MYFLAGS+= -DUSE_PBC=0 -DUSE_CELL_LIST=1 -DCELL_SEC_NEI=0 -DFIXED_CMASS=0
single: $(SINGLE_OBJS)
	$(CLINKER) $(MYFLAGS) $(CFLAGS)-o $(SINGLE) $(SINGLE_OBJS) $(LIBS)

eff:    MYFLAGS+= -DUSE_CELL_LIST=0 -DUSE_PBC=0 -DFIXED_CMASS=1
eff:	$(EFFECTIVE_OBJS)
	$(CLINKER) $(PROFILE) $(MYFLAGS) $(CFLAGS) -o $(EFFECTIVE) $(EFFECTIVE_OBJS) $(LIBS) $(PROFILE)

inter:  MYFLAGS+=-DUSE_PBC=1 -DUSE_CELL_LIST=1 -DFIXED_CMASS=0 -DCELL_SEC_NEI=0 
inter:  $(INTERACTING_OBJS) 
		$(CLINKER) $(MYFLAGS) $(CFLAGS) -o $(INTERACTING) $(INTERACTING_OBJS) $(LIBS)
anneal:  MYFLAGS+=-DUSE_PBC=1  -DUSE_CELL_LIST=1 -DFIXED_CMASS=0 -DCELL_SEC_NEI=0 
anneal:  $(ANNEAL_OBJS) 
		 $(CLINKER) $(MYFLAGS) $(CFLAGS) -o $(ANNEAL) $(ANNEAL_OBJS) $(LIBS)
# clean objects only
clean: 
	rm -f $(OBJDIR)/*.o
veryclean:
	rm -f *.o *~ $(SINGLE) $(EFFECTIVE) $(INTERACTING) $(WIDOM) $(ANNEAL)
	rm  -f *.o *~ *.dat *.xyz *.vmd $(SINGLE) $(INTERACTING) $(WIDOM)
	rm  -f $(OBJDIR)/*.o

$(OBJDIR)/%.o : %.c system.h ran250.h
	@echo $<
	@$(CC) $(PROFILE) $(MYFLAGS) $(CFLAGS) -c $(PARALLEL) $< -o $@

#MCinterdend.o   : system.h ran250.h preprocessor.h
#MCdendeff.o     : system.h ran250.h preprocessor.h
#singledend.o    : system.h ran250.h preprocessor.h
#dendrimers.o    : system.h ran250.h preprocessor.h
#system.o        : system.h preprocessor.h
#initialize.o    : system.h preprocessor.h ran250.h
#readInput.o     : system.h preprocessor.h
#mc.o            : system.h ran250.h preprocessor.h
#potential.o     : system.h preprocessor.h
#ran250.o        : ran250.h
#potential.o     : system.h preprocessor.h
#moves.o         : system.h ran250.h preprocessor.h
#energy.o        : system.h preprocessor.h
#cellLists.o     : system.h preprocessor.h
#preprocessor.h  :

