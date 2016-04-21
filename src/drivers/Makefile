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
VPATH = ../atmland:../base:../drivers:../land:../shared:../vegetation:$(OBJDIR)
MODS = $(shell find ../atmland/ -name '*.f90' -print)
MODS += $(shell find ../base/ -name '*.f90' -print)
MODS += $(shell find ../drivers/ -name 'mo*.f90' -print)
MODS += $(shell find ../land/ -name '*.f90' -print)
MODS += $(shell find ../shared/ -name '*.f90' -print)
MODS += $(shell find ../vegetation/ -name '*.f90' -print)

OBJFILES = $(notdir $(MODS))
OBJSMODS = $(OBJFILES:.f90=.o)
.SUFFIXES: .o .f90

all: $(BINDIR)/land.x
$(BINDIR)/land.x: $(LIBRARY) land_model_driver_test.f90
	$(F_C) $(F_O) land_model_driver_test.f90
	$(F_L) -o $(BINDIR)/land.x land_model_driver.o $(LIBRARY)
	@$(RM) land_model_driver.o
	@echo "land model is ok"

test_canopy: $(BINDIR)/test_canopy.x
$(BINDIR)/test_canopy.x: $(LIBRARY) canopy_test_model_driver.f90
	$(F_C) $(F_O) canopy_test_model_driver.f90
	$(F_L) -o $(BINDIR)/test_canopy.x canopy_test_model_driver.o $(LIBRARY)
	@$(RM) canopy_test_model_driver.o
	@echo "test canopy model is ok"

canopy: $(BINDIR)/canopy.x
$(BINDIR)/canopy.x:  $(LIBRARY) canopy_model_driver.f90 
	$(F_C) $(F_O) canopy_model_driver.f90
	$(F_L) -o $(BINDIR)/canopy.x canopy_model_driver.o $(LIBRARY)
	@$(RM) canopy_model_driver.o
	@echo "canopy model is ok"

plant: $(BINDIR)/plant.x
$(BINDIR)/plant.x:  $(LIBRARY) plant_model_driver.f90 
	$(F_C) $(F_O) plant_model_driver.f90
	$(F_L) -o $(BINDIR)/plant.x plant_model_driver.o $(LIBRARY)
	@$(RM) plant_model_driver.o
	@echo "plant model is ok"

lib: $(LIBRARY)
	@echo library is updated

$(LIBRARY): $(OBJSMODS)

clean:
	$(RM) $(LIBRARY) $(MODDIR)/*.mod $(OBJDIR)/*.o

help:
	@echo "commands available with make:"
	@echo "help:              this help" 
	@echo "all (or nothing):  stand alone land model (complete)" 
	@echo "canopy:            stand alone canopy model"
	@echo "test_canopy:       stand alone canopy model for testing"  
	@echo "plant:             stand alone plant model" 
	@echo "clean:             clean up library files" 

%.o: %.f90
	$(F_C) $(F_O) $<
	@$(A_C) $(LIBRARY) $*.o
	@cp $*.o $(OBJDIR)
	@$(RM) $*.o
#
# All dependencies
#
include ../shared/Make.inc.deps
include ../base/Make.inc.deps
include ../atmland/Make.inc.deps
include ../land/Make.inc.deps
include ../vegetation/Make.inc.deps
include ../drivers/Make.inc.deps
