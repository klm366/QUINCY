#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
BASE_LIB = $(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
LAND_LIB = $(LIBDIR)/lib_land.a
SXLAND_LIB = $(LAND_LIB)
MODEL_LIB = $(LIBDIR)/lib_vegetation.a
SXMODEL_LIB = $(MODEL_LIB)

#-
RM      = rm -f
STRIP   = strip
SIZE    = size

#-
MODS1 = mo_veg_memory.f90 \
        mo_veg_parameters.f90 \
        mo_veg_interface.f90 \
        mo_veg_respiration.f90 \
        mo_veg_photosynthesis.f90 \
        mo_veg_turnover.f90 \
        mo_veg_growth.f90 \
        mo_veg_radiation.f90 \
        mo_veg_output.f90

OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libland
	@$(M_K) libbase
	@$(M_K) m_all
	@echo  vegetation is OK

vegetation:
	@$(M_K) m_all
	@echo  vegetation is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libland:
	(cd ../land; $(M_K) -f Makefile land)

libbase:
	(cd ../base; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_veg_memory.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) 

$(MODEL_LIB)(mo_veg_interface.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(MODEL_LIB)(mo_veg_memory.o) \
  $(MODEL_LIB)(mo_veg_growth.o) \
  $(MODEL_LIB)(mo_veg_respiration.o) \
  $(MODEL_LIB)(mo_veg_turnover.o)

$(MODEL_LIB)(mo_veg_growth.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(MODEL_LIB)(mo_veg_memory.o)

$(MODEL_LIB)(mo_veg_respiration.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o) 

$(MODEL_LIB)(mo_veg_turnover.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(MODEL_LIB)(mo_veg_memory.o)

$(MODEL_LIB)(mo_veg_parameters.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_constants.o) 

$(MODEL_LIB)(mo_veg_radiation.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o)

$(MODEL_LIB)(mo_veg_photosynthesis.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o) \
  $(MODEL_LIB)(mo_veg_respiration.o) \
  $(MODEL_LIB)(mo_veg_radiation.o)

$(MODEL_LIB)(mo_veg_output.o): \
  $(BASE_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(BASE_LIB)(mo_base_constants.o) \
  $(BASE_LIB)(mo_base_configuration.o)

