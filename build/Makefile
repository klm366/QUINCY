#
# main folder is trunk/build
# with Makefile location: trunk/build/Makefile
# compiler specific folders will be created during build process, such as trunk/x86_64-intel
#
#---------------------------------------------------------------------
# include compiler specific settings using the symlink 'Make.inc.local' to any file Make.inc.local.*
# default: cd trunk/build && ln -s Make.inc.local Make.inc.local.ifort_BGC
include ./Make.inc.local
#---------------------------------------------------------------------
# COMP = x86_64-intel or x86_64-gfortran
# CDIR = ../x86_64-intel  = trunk/x86_64-intel
CDIR = $(addprefix ../,$(COMP))
# LIBDIR trunk/x86_64-intel/build/lib
LIBDIR = $(addsuffix /build/lib,$(CDIR))
# MODDIR trunk/x86_64-intel/build/mod
MODDIR = $(addsuffix /build/mod,$(CDIR))
# OBJDIR = trunk/x86_64-intel/build/mod
OBJDIR = $(MODDIR)
# BDIR = trunk/x86_64-intel/bin
BDIR = $(addsuffix /bin,$(CDIR))
# $(strip string)  Removes leading and trailing whitespace from string and replaces each internal sequence of one or more whitespace characters with a single space.
BINDIR = $(strip $(BDIR))
# if not present, create following folder plus parent folders
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
AR      = ar -r
#---------------------------------------------------------------------
# The VPATH is a list of directories to be searched for missing source files.
# The list can be separated by spaces or colons - the best approach is probably to use spaces as this works across all platforms.
VPATH = ../src/adapters ../src/atmland ../src/echam ../src/base ../src/drivers ../src/shared ../src/assimilation ../src/radiation ../src/soil ../src/soil_biogeochemistry ../src/phenology ../src/vegetation $(OBJDIR)
# 
MODS = $(shell find ../src/adapters/ -name '*.f90' -print)
MODS += $(shell find ../src/atmland/ -name '*.f90' -print)
MODS += $(shell find ../src/echam/ -name '*.f90' -print)
MODS += $(shell find ../src/base/ -name '*.f90' -print)
MODS += $(shell find ../src/drivers/ -name 'mo*.f90' -print)
MODS += $(shell find ../src/shared/ -name '*.f90' -print)
MODS += $(shell find ../src/assimilation/ -name '*.f90' -print)
MODS += $(shell find ../src/radiation/ -name '*.f90' -print)
MODS += $(shell find ../src/soil/ -name '*.f90' -print)
MODS += $(shell find ../src/soil_biogeochemistry/ -name '*.f90' -print)
MODS += $(shell find ../src/phenology/ -name '*.f90' -print)
MODS += $(shell find ../src/vegetation/ -name '*.f90' -print)

# 'notdir' takes away the path from the file name leaving just the file name (so /x/y/foo.a becomes foo.a)
OBJFILES = $(notdir $(MODS))
OBJSMODS = $(OBJFILES:.f90=.o)
.SUFFIXES: .o .f90


# main target of make; removes all trunk/build/*.o files after build process has finished successfully
all: $(BINDIR)/land.x
	@rm -f ./*.o

# rules for creating the land.x binary
$(BINDIR)/land.x: $(LIBRARY) land_model_driver.f90
	$(F_C) $(F_O) ../src/drivers/land_model_driver.f90
	$(F_L) -o $(BINDIR)/land.x land_model_driver.o $(LIBRARY)
	@$(RM) land_model_driver.o
	@echo -e "\n land.x binary was built successfully"

# create libslm.a from .o files
$(LIBRARY): $(OBJSMODS)
	@cp $(OBJDIR)/*.o ./
	$(AR) $@ $(OBJSMODS)

# removes all .o .mod files and the libslm.a
clean:
	$(RM) $(LIBRARY) $(MODDIR)/*.mod $(OBJDIR)/*.o ./*.o
	@echo -e "\n all *.o and *.mod files and the libslm.a were removed"

help:
	@echo "## Commands available with make ##"
	@echo " "
	@echo " help:              this help" 
	@echo " all (or nothing):  stand alone model (land.x binary including all SLM models)" 
	@echo " clean:             clean up library files" 

# create .o files from .f90 files and copy to OBJDIR which is in the VPATH so that only changed/new .f90 files have to be compiled
%.o: %.f90
	$(F_C) $(F_O) $<
	@cp $*.o $(OBJDIR)

#
# All dependencies
#
include ../src/echam/Make.inc.deps
include ../src/base/Make.inc.deps
include ../src/atmland/Make.inc.deps
include ../src/shared/Make.inc.deps
include ../src/assimilation/Make.inc.deps
include ../src/radiation/Make.inc.deps
include ../src/soil/Make.inc.deps
include ../src/soil_biogeochemistry/Make.inc.deps
include ../src/phenology/Make.inc.deps
include ../src/vegetation/Make.inc.deps
include ../src/drivers/Make.inc.deps
include ../src/adapters/Make.inc.deps
