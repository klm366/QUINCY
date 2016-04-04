#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
SHARED_LIB = $(LIBDIR)/lib_shared.a
SXSHARED_LIB = $(SHARED_LIB)
BASE_LIB =$(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
LAND_LIB = $(LIBDIR)/lib_land.a
SXLAND_LIB = $(LAND_LIB)
VEG_LIB = $(LIBDIR)/lib_vegetation.a
SXVEG_LIB = $(VEG_LIB)
ATM_LIB = $(LIBDIR)/lib_atmland.a
SXATM_LIB = $(ATM_LIB)
MODEL_LIB = $(LIBDIR)/lib_drivers.a
SXMODEL_LIB = $(MODEL_LIB)

RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_land_model_interface.f90 

OBJSMODS1 = $(MODS1:.f90=.o)

all: 
	@$(M_K) all_libs 
	@$(M_K) m_all
	@echo main is OK

m_all: $(BINDIR)/land.x
$(BINDIR)/land.x: $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) $(SXSHARED_LIB) \
  $(MODEL_LIB)(mo_land_model_interface.o) $(MODEL_LIB)(land_model_driver.o)
	$(A_G) $(MODEL_LIB) land_model_driver.o
	$(F_L) -o $(BINDIR)/land.x land_model_driver.o \
            $(SXMODEL_LIB) $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) $(SXSHARED_LIB)
	$(RM) land_model_driver.o

canopy: 
	@$(M_K) all_libs 
	@$(M_K) m_canopy
	@echo canopy is OK

m_canopy: $(BINDIR)/canopy.x
$(BINDIR)/canopy.x: $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) \
  $(MODEL_LIB)(canopy_model_driver.o)
	$(A_G) $(MODEL_LIB) canopy_model_driver.o
	$(F_L) -o $(BINDIR)/canopy.x canopy_model_driver.o \
            $(SXMODEL_LIB) $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) 
	$(RM) canopy_model_driver.o 

transient_canopy: 
	@$(M_K) all_libs 
	@$(M_K) m_transient_canopy
	@echo transient canopy is OK

m_transient_canopy: $(BINDIR)/transient_canopy.x
$(BINDIR)/transient_canopy.x: $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) $(SXSHARED_LIB) \
  $(MODEL_LIB)(transient_canopy_model_driver.o)
	$(A_G) $(MODEL_LIB) transient_canopy_model_driver.o
	$(F_L) -o $(BINDIR)/transient_canopy.x transient_canopy_model_driver.o \
            $(SXMODEL_LIB) $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) $(SXSHARED_LIB) 
	$(RM) transient_canopy_model_driver.o 

plant: 
	@$(M_K) all_libs 
	@$(M_K) m_plant
	@echo plant is OK

m_plant: $(BINDIR)/plant.x
$(BINDIR)/plant.x: $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) $(SXSHARED_LIB) \
  $(MODEL_LIB)(plant_model_driver.o)
	$(A_G) $(MODEL_LIB) plant_model_driver.o
	$(F_L) -o $(BINDIR)/plant.x $(SXDRIVER_LIB) \
            $(SXMODEL_LIB) $(SXVEG_LIB) $(SXLAND_LIB) $(SXATM_LIB) $(SXBASE_LIB) $(SXSHARED_LIB) 
	$(RM) plant_model_driver.o

all_libs:
	@$(M_K) libshared
	@$(M_K) libbase
	@$(M_K) libatmland
	@$(M_K) libland
	@$(M_K) libvegetation

libatmland:
	(cd ../atmland; $(M_K) -f Makefile atmland)

libland:
	(cd ../land; $(M_K) -f Makefile land)

libvegetation:
	(cd ../vegetation; $(M_K) -f Makefile vegetation)

libbase:
	(cd ../base; $(M_K) -f Makefile)

libshared:
	(cd ../shared; $(M_K) -f Makefile shared)

clean:
	(cd ../shared; $(M_K) -f Makefile clean)
	(cd ../base; $(M_K) -f Makefile clean)
	(cd ../atmland; $(M_K) -f Makefile clean)
	(cd ../land; $(M_K) -f Makefile clean)
	(cd ../vegetation; $(M_K) -f Makefile clean)
	$(RM) $(MODEL_LIB)

help:
	@echo "commands available with make:"
	@echo "help:              this help" 
	@echo "all (or nothing):  stand alone land model (complete)" 
	@echo "canopy:            stand alone canopy model" 
	@echo "plant:             stand alone plant model without canopy" 
	@echo "clean:             clean up library files" 

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	@$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

$(MODEL_LIB)(mo_land_model_interface.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(VEG_LIB)(mo_veg_parameters.o) 
