#---------------------------------------------------------------------
include ../../util/Make.inc.local
#---------------------------------------------------------------------
SHARED_LIB = $(LIBDIR)/lib_shared.a
SXSHARED_LIB = $(SHARED_LIB)
BASE_LIB = $(LIBDIR)/lib_base.a
SXBASE_LIB = $(BASE_LIB)
ATM_LIB = $(LIBDIR)/lib_atmland.a
SXATM_LIB = $(ATM_LIB)
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
        mo_veg_soil.f90 \
        mo_veg_output.f90 \
        mo_veg_util.f90 \
        mo_veg_config_class.f90


OBJSMODS1 = $(MODS1:.f90=.o)
#-
.PRECIOUS : $(MODEL_LIB)

#-
all:
	@$(M_K) libland
	@$(M_K) libatm
	@$(M_K) libbase
	@$(M_K) libshared
	@$(M_K) m_all
	@echo  vegetation is OK

vegetation:
	@$(M_K) m_all
	@echo  vegetation is OK

m_all: $(MODEL_LIB)($(OBJSMODS1))

libland:
	(cd ../land; $(M_K) -f Makefile land)

libatm:
	(cd ../atmland; $(M_K) -f Makefile atmland)

libshared:
	(cd ../shared; $(M_K) -f Makefile shared)

libbase:
	(cd ../base; $(M_K) -f Makefile)

$(MODEL_LIB)(%.o): %.f90
	$(F_C) $(F_O) $*.f90
	$(A_C) $(MODEL_LIB) $*.o
	@$(RM) $*.o

clean:
	$(RM) $(MODEL_LIB)

$(MODEL_LIB)(mo_veg_memory.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) 

$(MODEL_LIB)(mo_veg_interface.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(MODEL_LIB)(mo_veg_memory.o) \
  $(MODEL_LIB)(mo_veg_growth.o) \
  $(MODEL_LIB)(mo_veg_photosynthesis.o) \
  $(MODEL_LIB)(mo_veg_respiration.o) \
  $(MODEL_LIB)(mo_veg_turnover.o)

$(MODEL_LIB)(mo_veg_growth.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(ATM_LIB)(mo_atmland_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o) \
  $(MODEL_LIB)(mo_veg_memory.o)

$(MODEL_LIB)(mo_veg_respiration.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(MODEL_LIB)(mo_veg_parameters.o) 

$(MODEL_LIB)(mo_veg_soil.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(MODEL_LIB)(mo_veg_parameters.o) 

$(MODEL_LIB)(mo_veg_turnover.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(MODEL_LIB)(mo_veg_memory.o)

$(MODEL_LIB)(mo_veg_parameters.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_constants.o) 

$(MODEL_LIB)(mo_veg_radiation.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(ATM_LIB)(mo_atmland_memory.o) \
  $(ATM_LIB)(mo_atmland_forcing.o) \
  $(ATM_LIB)(mo_atmland_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o)

$(MODEL_LIB)(mo_veg_photosynthesis.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(ATM_LIB)(mo_atmland_memory.o) \
  $(ATM_LIB)(mo_atmland_forcing.o) \
  $(ATM_LIB)(mo_atmland_constants.o) \
  $(ATM_LIB)(mo_atmland_util.o) \
  $(LAND_LIB)(mo_land_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o) \
  $(MODEL_LIB)(mo_veg_respiration.o) \
  $(MODEL_LIB)(mo_veg_radiation.o)

$(MODEL_LIB)(mo_veg_output.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(BASE_LIB)(mo_base_constants.o) \
  $(BASE_LIB)(mo_base_config.o)

$(MODEL_LIB)(mo_veg_util.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_memory.o) \
  $(BASE_LIB)(mo_base_util_functions.o) \
  $(BASE_LIB)(mo_base_constants.o) \
  $(MODEL_LIB)(mo_veg_parameters.o)

$(MODEL_LIB)(mo_veg_config_class.o): \
  $(SHARED_LIB)(mo_kind.o) \
  $(BASE_LIB)(mo_base_config.o) \
  $(SHARED_LIB)(mo_namelist.o)

