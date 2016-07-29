#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
CDIR = $(addprefix ../../,$(COMP))
LIBDIR = $(addsuffix /build/lib,$(CDIR))
MODDIR = $(addsuffix /build/mod,$(CDIR))
OBJDIR = $(MODDIR)
BDIR = $(addsuffix /bin,$(CDIR))
BINDIR = $(strip $(BDIR))
$(shell [ ! -d $(MODDIR) ] && mkdir -p $(MODDIR))
$(shell [ ! -d $(LIBDIR) ] && mkdir -p $(LIBDIR))
$(shell [ ! -d $(BINDIR) ] && mkdir -p $(BINDIR))
#---------------------------------------------------------------------
#-
LIBRARY = $(LIBDIR)/libslm.a
#---------------------------------------------------------------------
RM      = rm -f
STRIP   = strip
SIZE    = size
#---------------------------------------------------------------------
#-
VPATH = ../atmland:../echam:../base:../drivers:../shared:../assimilation:../radiation:../soil:../vegetation:$(OBJDIR)
MODS = $(shell find ../atmland/ -name '*.f90' -print)
MODS += $(shell find ../echam/ -name '*.f90' -print)
MODS += $(shell find ../base/ -name '*.f90' -print)
MODS += $(shell find ../drivers/ -name 'mo*.f90' -print)
MODS += $(shell find ../shared/ -name '*.f90' -print)
MODS += $(shell find ../assimilation/ -name '*.f90' -print)
MODS += $(shell find ../radiation/ -name '*.f90' -print)
MODS += $(shell find ../soil/ -name '*.f90' -print)
MODS += $(shell find ../vegetation/ -name '*.f90' -print)

OBJFILES = $(notdir $(MODS))
OBJSMODS = $(OBJFILES:.f90=.o)
.SUFFIXES: .o .f90

all: $(BINDIR)/land.x
$(BINDIR)/land.x: $(LIBRARY) land_model_driver.f90
	$(F_C) $(F_O) land_model_driver.f90
	$(F_L) -o $(BINDIR)/land.x land_model_driver.o $(LIBRARY)
	@$(RM) land_model_driver.o
	@echo "land model is ok"


lib: $(LIBRARY)
	@echo library is updated

$(LIBRARY): $(OBJSMODS)

clean:
	$(RM) $(LIBRARY) $(MODDIR)/*.mod $(OBJDIR)/*.o

help:
	@echo "commands available with make:"
	@echo "help:              this help" 
	@echo "all (or nothing):  stand alone land model (complete)" 
	@echo "clean:             clean up library files" 

%.o: %.f90
	$(F_C) $(F_O) $<
	@$(A_C) $(LIBRARY) $*.o
	@cp $*.o $(OBJDIR)
	@$(RM) $*.o
#
# All dependencies
#
include ../echam/Make.inc.deps
include ../base/Make.inc.deps
include ../atmland/Make.inc.deps
include ../shared/Make.inc.deps
include ../assimilation/Make.inc.deps
include ../radiation/Make.inc.deps
include ../soil/Make.inc.deps
include ../vegetation/Make.inc.deps
include ../drivers/Make.inc.deps
